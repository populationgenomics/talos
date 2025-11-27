# Script for gathering large inputs required by Talos
# This runs a download from multiple different sources
# Output file names created/expected by this script match the initial configuration file content in `talos.config` and `annotation.config`

TMX=$(command -v tmux)
POLL_INTERVAL=1
SESSION_NAME="download_manager"
WINDOW_NAME="Downloads"
TMX_WINDOW_ID=""
declare -a DOWNLOAD_TARGETS=()

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
  # Use a temporary partial file so incomplete downloads never appear at the final path.
  # We keep resume ability by always resuming (-C -) against the .part file and only mv to final name on success.
  local final="$output"
  local part_file="${final}.part"

  # Record expected final file for later summary (avoid duplicates)
  DOWNLOAD_TARGETS+=("$final")

  # If the final file already exists, skip (assume complete). Could add checksum logic here if desired.
  if [ -f "$final" ]; then
    echo "[INFO] $final already present, skipping download."
    return 0
  fi

  # Banner shown once per (re)attempt.
  # curl exit codes: 18 = partial file; we keep part_file for later resume.
  # On success (exit 0) we atomically mv into place.
  # NOTE: All internal $ variables are escaped (\$) so they are evaluated when the command runs, not now.
  local script_name="$(basename "$0")"
  local cmd="echo \"$banner\"; echo; curl -C - -# -L --fail -o \"$part_file\" \"$url\" && mv -f \"$part_file\" \"$final\" || { rc=\$?; if [ \"\$rc\" -ne 0 ]; then echo \"[WARN] $final Download failed for (exit \\${rc}). Restart ${script_name} to gracefully resume download.\"; fi; }"
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
start_download https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz

# Jax lab file for phenotype matching
start_download https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2025-05-06/hp.obo

start_download https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2025-05-06/genes_to_phenotype.txt

# latest clinvarbitration data
CLINVAR="clinvarbitration.tar.gz"
start_download https://zenodo.org/records/17577130/files/ClinvArbitration_Nov_2025_clinvar_decisions.release.tar.gz?download=1 "${CLINVAR}"

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

# Final status summary
summary_fail=0
script_name_summary="$(basename "$0")"
for target in "${DOWNLOAD_TARGETS[@]}"; do
  if [ ! -f "$target" ]; then
    echo "[MISSING] $target"
    summary_fail=$((summary_fail+1))
  fi
done

if [ $summary_fail -eq 0 ]; then
  echo "[SUCCESS] All ${#DOWNLOAD_TARGETS[@]} downloads completed successfully."
else
  echo "[SUMMARY] $summary_fail of ${#DOWNLOAD_TARGETS[@]} downloads missing. Restart ${script_name_summary} to gracefully resume."
fi
