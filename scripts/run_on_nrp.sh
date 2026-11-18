#!/usr/bin/env bash
# Runs inside the NRP job container (started by nrp/gr-footprint-job.yaml).
#  1. builds the conda env once and caches it on the PVC
#  2. converts the BAMs to a fragments file (first run only)
#  3. executes the notebook with papermill
#  4. saves the executed notebook + plots to /data/GR_dex/results/<run-id>/ for you to download
set -euo pipefail

REPO_DIR="${REPO_DIR:-/data/code}"          # notebook, environment.yml and scripts copied here
NOTEBOOK="${NOTEBOOK:-dex_footprinting.ipynb}"
DATA_DIR="${DATA_DIR:-/data/GR_dex}"          # working folder on the PVC (wd in the notebook)
export SCPRINTER_DATA="/data/scprinter_cache" # genome + Tn5 bias downloaded once, reused
RUN_ID="$(date +%Y%m%d-%H%M%S)"
export N_JOBS="${N_JOBS:-8}"                   # match the job's CPU request (nproc can report the whole node)
export GR_WD="$DATA_DIR"                       # the notebook reads wd from this

cd "$REPO_DIR"
mkdir -p "$DATA_DIR" "$SCPRINTER_DATA"

# ---- 1. environment
# Python imports thousands of small files, which is very slow on network storage (the PVC).
# So the env always lives on the pod's local disk (/tmp):
#   first run : build it locally, then save it to the PVC as ONE archive (fast to write/read)
#   later runs: unpack that archive to local disk (a minute or two)
ENV_TAR="/data/envs/scprinter.tar.gz"
LOCAL_ENV="/tmp/scprinter_env"
mkdir -p /data/envs

build_env() {
  echo "== building conda env on local disk (~15-25 min)"
  rm -rf "$LOCAL_ENV"
  mamba env create -p "$LOCAL_ENV" -f environment.yml
  PATH="$LOCAL_ENV/bin:$PATH" "$LOCAL_ENV/bin/pip" install --no-deps "git+https://github.com/buenrostrolab/scPrinter.git"
  mamba clean -afy                          # drop the download cache to save local disk
  echo "== saving env to $ENV_TAR for future runs"
  mamba install -y -n base -c conda-forge conda-pack
  conda-pack -p "$LOCAL_ENV" -o /tmp/scprinter_env.tar.gz --format tar.gz --compress-level 1 \
             --ignore-missing-files --ignore-editable-packages
  cp /tmp/scprinter_env.tar.gz "$ENV_TAR.partial"
  sync "$ENV_TAR.partial"                   # make sure every byte is on the storage before naming it final
  if gzip -t "$ENV_TAR.partial"; then       # verify the copy is complete
    mv "$ENV_TAR.partial" "$ENV_TAR"
    echo "== env saved and verified"
  else
    echo "!! saved env failed verification; future runs will rebuild"
    rm -f "$ENV_TAR.partial"
  fi
  rm -f /tmp/scprinter_env.tar.gz
}

if [ -s "$ENV_TAR" ]; then
  echo "== unpacking saved env to local disk"
  rm -rf "$LOCAL_ENV"; mkdir -p "$LOCAL_ENV"
  if tar -xzf "$ENV_TAR" -C "$LOCAL_ENV" && "$LOCAL_ENV/bin/conda-unpack"; then
    echo "== env unpacked"
  else
    echo "!! saved env is damaged; deleting it and rebuilding"
    rm -f "$ENV_TAR"
    build_env
  fi
else
  build_env
fi
export PATH="$LOCAL_ENV/bin:$PATH"
echo "== env ready: $(which python)"
python -c "import numpy, pyBigWig, scprinter; print('numpy', numpy.__version__, '| pyBigWig numpy:', pyBigWig.numpy, flush=True)"

# ---- 2. BAM -> fragments (skipped if the fragments file already exists)
FRAGS="$DATA_DIR/GR_all.frags.tsv.gz"
if [ ! -s "$FRAGS" ]; then
  echo "== converting BAMs to fragments"
  bash scripts/bam_to_fragments.sh "$FRAGS" \
       "no_dex:/data/bam/${NODEX_BAM:-nodex.bam}" \
       "dex:/data/bam/${DEX_BAM:-dex.bam}"
fi

# ---- 3. run the notebook
OUT_DIR="$DATA_DIR/results/nrp_${RUN_ID}"
mkdir -p "$OUT_DIR"
echo "== running $NOTEBOOK"
set +e
papermill "$NOTEBOOK" "$OUT_DIR/${NOTEBOOK%.ipynb}_executed.ipynb" \
  -k python3 \
  -p wd "$DATA_DIR" \
  -p n_jobs "$N_JOBS" \
  --log-output
STATUS=$?
set -e

# ---- 4. collect small outputs for download (the big data files stay where they are)
cp -r "$DATA_DIR/plots" "$OUT_DIR/" 2>/dev/null || true
find "$DATA_DIR" -maxdepth 1 \( -name "*.csv" -o -name "*.tsv" \) -size -20M -exec cp {} "$OUT_DIR/" \; 2>/dev/null || true
echo "== results saved in $OUT_DIR (notebook exit status $STATUS)"
ls -lh "$OUT_DIR"

exit $STATUS