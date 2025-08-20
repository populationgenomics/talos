# Script for gathering large inputs required by Talos
# This runs a download from multiple different sources
# Output file names created/expected by this script match the initial configuration file content in `talos.config` and `annotation.config`

TMX=$(command -v tmux)
POLL_INTERVAL=1
SESSION_NAME="download_manager"
WINDOW_NAME="Downloads"
TMX_WINDOW_ID=""

if [ ! -z "$TMX" ] && [ -z "$TMUX" ]; then
  # tmux installed, but not in a tmux session. restart in tmux.
  tmux new-session -d -s "$SESSION_NAME" sh "$0" "$@"
  tmux attach-session -t "$SESSION_NAME"
  exit 0
fi

cleanup() {
  echo "[INFO] Caught exit signal. Cleaning up..."
  CURRENT_PGID=$(ps -o pgid= -p $$ | tr -d '[:space:]')
  pkill -SIGTERM -g "$CURRENT_PGID" -f curl

  if [ -z $TMX ]; then
	wait
  else
    # Kill the tmux window we created for the downloads
    if $TMX list-windows | grep -q "$WINDOW_NAME"; then
      $TMX kill-window -t "$TMX_WINDOW_ID"
    fi
  fi
  echo "[INFO] Cleanup complete."
  exit 0
}

start_download() {
  local url="$1"
  local output="$2"
  local banner="$3"

  if [ -z "$output" ]; then
    output=$(basename "$url")
  fi
  if [ -z "$banner" ]; then
    banner="Downloading [$output]"
  fi
  local cmd="echo \"$banner\"; echo; curl -C - -# -L --fail -o \"$output\" \"$url\";"
  if [ -z $TMX ]; then
     ( eval "$cmd" ) &
  else
    echo "[JOB START] Starting download for: $url to $output."
    if [ -z "$TMX_WINDOW_ID" ]; then
      TMX_WINDOW_ID=$($TMX new-window -P -n "$WINDOW_NAME" "$cmd")
    else
      $TMX split-window -v -t "$TMX_WINDOW_ID" "$cmd"
      $TMX select-layout -t "$TMX_WINDOW_ID" even-vertical
    fi
  fi
}

await() {
  if [ -z $TMX ]; then
    wait
  else
  last_pane_count=-1
    while true; do
      pane_count=$($TMX list-panes -t "$TMX_WINDOW_ID" | wc -l)
      if [ $last_pane_count != $pane_count ]; then
        $TMX select-layout -t "$TMX_WINDOW_ID" even-vertical
        last_pane_count=$pane_count
      fi
      if [ $pane_count == 0 ]; then
        break;
      fi
      sleep "$POLL_INTERVAL"
    done
  fi
}

trap cleanup SIGINT SIGTERM

# Echtvar-encoded gnomAD 4.1 population frequencies - this is a big one (~6GB) so it's started early and backgrounded
ECHTVAR_FILE="gnomad_4.1_region_merged_GRCh38_whole_genome"
start_download https://zenodo.org/records/15222100/files/gnomad_4.1_region_merged_GRCh38_whole_genome?download=1 "${ECHTVAR_FILE}" "\
Downloading Echtvar from https://zenodo.org/records/15222100"

# Monarch phenotype DB - another large download. 17GB decompressed
COMPRESSED_PHENIO="phenio.db.gz"
DECOMPRESSED_PHENIO="phenio.db"
# if it doesn't exist, download it
if [ ! -f "${DECOMPRESSED_PHENIO}" ]; then
  echo "[INFO] phenio.db not found, downloading..."
  start_download https://data.monarchinitiative.org/monarch-kg/latest/phenio.db.gz
else
  echo "[INFO] phenio.db already exists, skipping download."
fi

# GRCh38 reference genome
GRCh38="ref.fa.gz"
GRCh38_decompressed="ref.fa"
# if it doesn't exist, download it
if [ ! -f "${GRCh38_decompressed}" ]; then
  echo "[INFO] compressed ref.fa not found, attempting download..."
    GRCh38_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    start_download ${GRCh38_URL} ${GRCh38} "Downloading compressed GRCh38 reference genome from UCSC"
else
  echo "[INFO] ${GRCh38_decompressed} already exists, skipping download."
fi

# MANE gene data
start_download https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz

# Ensembl GFF3 data
start_download https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh38.113.gff3.gz

# Jax lab file for phenotype matching
start_download https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2025-05-06/hp.obo

start_download https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2025-05-06/genes_to_phenotype.txt

# latest clinvarbitration data
CLINVAR="clinvarbitration.tar.gz"
start_download https://zenodo.org/records/16792026/files/clinvarbitration_Aug_2025_clinvar_decisions.release.tar.gz?download=1 "${CLINVAR}"

# AlphaMissense raw data
AM="AlphaMissense_hg38.tsv.gz"
start_download "https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz?download=1" "${AM}"

await

#if compressed exists, but decompressed doesn't, gunzip it
if [ ! -f "${GRCh38_decompressed}" ] && [ -f "${GRCh38}" ]; then
    gunzip ${GRCh38}
fi

# same for the phenio files
if [ ! -f "${DECOMPRESSED_PHENIO}" ] && [ -f "${COMPRESSED_PHENIO}" ]; then
    gunzip ${COMPRESSED_PHENIO}
fi
