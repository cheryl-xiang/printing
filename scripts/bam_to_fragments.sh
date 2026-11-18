#!/usr/bin/env bash
# Convert bulk ATAC-seq BAM files into one scPrinter fragments file.
#
# Usage:
#   bash scripts/bam_to_fragments.sh OUT.frags.tsv.gz  label1:file1.bam  label2:file2.bam ...
# Example:
#   bash scripts/bam_to_fragments.sh /data/GR_dex/GR_all.frags.tsv.gz  no_dex:/data/bam/nodex.bam  dex:/data/bam/dex.bam
#
# - Each label becomes the "barcode" (column 4) for every fragment from that BAM.
# - Paired-end BAMs -> one line per fragment (chrom, start, end, barcode, 1).
# - Single-end BAMs -> one line per read with strand (chrom, start, end, barcode, 1, strand);
#   scPrinter treats 6-column files as single-end.
# - Filters: mapped, primary, MAPQ >= MIN_MAPQ (default 30; set MIN_MAPQ=0 if the BAM was already
#   MAPQ-filtered, e.g. ENCODE pipeline *.nodup.* BAMs), not QC-fail, duplicates removed, main chromosomes only.
# - Output is sorted by barcode, then position (use sorted_by_barcode=True in import_fragments).
set -euo pipefail

OUT="$1"; shift
THREADS="${THREADS:-${N_JOBS:-8}}"
MIN_MAPQ="${MIN_MAPQ:-30}"
TMP="$(dirname "$OUT")/tmp_frags"      # per-sample results kept here so finished samples are reused if a later step fails
mkdir -p "$TMP"

# check every BAM exists before doing any work
for spec in "$@"; do
  bam="${spec#*:}"
  if [ ! -s "$bam" ]; then
    echo "!! BAM not found: $bam"; echo "   files in $(dirname "$bam"):"; ls -lh "$(dirname "$bam")"; exit 1
  fi
done
MAIN_CHROM='^chr([0-9]+|X|Y)$'

for spec in "$@"; do
  label="${spec%%:*}"; bam="${spec#*:}"
  echo "== $label: $bam"
  if [ -f "$TMP/$label.done" ] && [ -s "$TMP/$label.frags.tsv" ]; then
    echo "   already converted, reusing $(wc -l < "$TMP/$label.frags.tsv") fragments"
    continue
  fi

  # paired-end if any read has the "paired" flag (0x1); samtools does the flag test (no awk and() needed)
  is_paired=$( (samtools view -f 1 "$bam" 2>/dev/null || true) | head -n 1 | wc -l )
  if [ "$is_paired" -gt 0 ]; then
    echo "   paired-end"
    # properly paired, both mapped, primary, not QC-fail, not duplicate, not supplementary (-F 3852)
    samtools view -@ "$THREADS" -b -f 2 -F 3852 -q "$MIN_MAPQ" "$bam" \
      | samtools sort -n -@ "$THREADS" -o "$TMP/$label.namesorted.bam" -
    bedtools bamtobed -bedpe -i "$TMP/$label.namesorted.bam" 2>/dev/null \
      | awk -v lab="$label" -v re="$MAIN_CHROM" 'BEGIN{OFS="\t"}
          $1==$4 && $1 ~ re {
            s = ($2 < $5 ? $2 : $5); e = ($3 > $6 ? $3 : $6);
            print $1, s, e, lab, 1
          }' \
      | LC_ALL=C sort -k1,1 -k2,2n -k3,3n -u -S 4G --parallel="$THREADS" -T "$TMP" \
      > "$TMP/$label.frags.tsv"          # -u on chrom/start/end removes PCR duplicates
    rm -f "$TMP/$label.namesorted.bam"
  else
    echo "   single-end"
    samtools view -@ "$THREADS" -b -F 3844 -q "$MIN_MAPQ" "$bam" \
      | bedtools bamtobed -i - \
      | awk -v lab="$label" -v re="$MAIN_CHROM" 'BEGIN{OFS="\t"} $1 ~ re {print $1,$2,$3,lab,1,$6}' \
      | LC_ALL=C sort -k1,1 -k2,2n -k3,3n -k6,6 -u -S 4G --parallel="$THREADS" -T "$TMP" \
      > "$TMP/$label.frags.tsv"
  fi
  echo "   $(wc -l < "$TMP/$label.frags.tsv") fragments"
  touch "$TMP/$label.done"
done

echo "== combining -> $OUT"
labels=(); for spec in "$@"; do labels+=("$TMP/${spec%%:*}.frags.tsv"); done
cat "${labels[@]}" \
  | LC_ALL=C sort -k4,4 -k1,1 -k2,2n -S 4G --parallel="$THREADS" -T "$TMP" \
  | bgzip -@ "$THREADS" > "$OUT.partial"
mv "$OUT.partial" "$OUT"                # only gets its final name once complete
rm -rf "$TMP"
echo "== done"; { zcat "$OUT" || true; } | head -3      # preview; "|| true" stops head closing the pipe from counting as an error
