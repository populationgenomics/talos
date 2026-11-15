#!/usr/bin/env bash

# Script for gathering large inputs required by Talos
# This runs a download from multiple different sources
# Output file names created/expected by this script match the initial configuration file content in `talos.config` and `annotation.config`
set -o pipefail

TMX=$(command -v tmux)
POLL_INTERVAL=1
SESSION_NAME="download_manager"
WINDOW_NAME="Downloads"
TMX_WINDOW_ID=""
declare -a DOWNLOAD_TARGETS=()

# All log lines (banners, INFO/WARN/JOB START/SUMMARY) are tee'd here so they
# survive and stay readable even if the tmux session that ran them is
# detached, killed, or never attached to. Appended, not truncated, so a
# restart after a crash keeps the download history.
LOGFILE="${LOGFILE:-$(pwd)/gather_files.log}"
touch "$LOGFILE" 2>/dev/null

log() {
  echo "$@" | tee -a "$LOGFILE"
}

if [ ! -z "$TMX" ] && [ -z "$TMUX" ]; then
  # tmux installed, but not in a tmux session. restart in tmux.
  # Re-exec with bash, not sh: this script uses arrays, `local` and `declare`.
  log "[INFO] Running downloads in tmux session '$SESSION_NAME'. Log persisted to $LOGFILE"
  log "[INFO] To follow progress from outside tmux: tail -f $LOGFILE"
  LOGFILE="$LOGFILE" tmux new-session -d -s "$SESSION_NAME" bash "$0" "$@"
  tmux attach-session -t "$SESSION_NAME"
  # The tmux session owns the real output; replay the log so the invoking
  # terminal still holds a record after the session detaches or is killed.
  echo
  echo "[INFO] tmux session ended. Log retained at $LOGFILE"
  exit 0
fi

cleanup() {
  log "[INFO] Caught exit signal. Cleaning up..."
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
  log "[INFO] Cleanup complete."
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
    log "[INFO] $final already present, skipping download."
    return 0
  fi

  # Banner shown once per (re)attempt.
  # curl exit codes: 18 = partial file; we keep part_file for later resume.
  # On success (exit 0) we atomically mv into place.
  # NOTE: All internal $ variables are escaped (\$) so they are evaluated when the command runs, not now.
  # Status lines are appended to the logfile as well as printed to the pane, because the pane
  # is destroyed on completion/cleanup and its scrollback goes with it. The curl progress bar
  # (-# writes to stderr) is deliberately NOT logged - it is carriage-return spam in a file.
  local script_name="$(basename "$0")"
  local cmd="echo \"$banner\" | tee -a \"$LOGFILE\"; echo; curl -C - -# -L --fail -o \"$part_file\" \"$url\" && mv -f \"$part_file\" \"$final\" && echo \"[DONE] $final download complete.\" | tee -a \"$LOGFILE\" || { rc=\$?; if [ \"\$rc\" -ne 0 ]; then echo \"[WARN] $final Download failed for (exit \\${rc}). Restart ${script_name} to gracefully resume download.\" | tee -a \"$LOGFILE\"; fi; }"
  if [ -z $TMX ]; then
     ( eval "$cmd" ) &
  else
    log "[JOB START] Starting download for: $url to $output."
    if [ -z "$TMX_WINDOW_ID" ]; then
      TMX_WINDOW_ID=$($TMX new-window -P -n "$WINDOW_NAME" "$cmd")
    else
      # tmux refuses a split once panes get too small ("no space for new pane").
      # A failed split must not lose the download - fall back to a plain
      # background job whose output goes to the logfile only.
      if $TMX split-window -v -t "$TMX_WINDOW_ID" "$cmd" 2>/dev/null; then
        $TMX select-layout -t "$TMX_WINDOW_ID" tiled
      else
        log "[INFO] No tmux pane space left for $output; running in background (see $LOGFILE)."
        ( eval "$cmd" ) >/dev/null 2>&1 &
      fi
    fi
  fi
}

await() {
  if [ -z $TMX ]; then
    wait
  else
    # No downloads were started in tmux — nothing to wait on.
    if [ -z "$TMX_WINDOW_ID" ]; then
      wait
      return 0
    fi
  last_pane_count=-1
    while true; do
      pane_count=$($TMX list-panes -t "$TMX_WINDOW_ID" 2>/dev/null | wc -l)
      if [ $last_pane_count != $pane_count ]; then
        $TMX select-layout -t "$TMX_WINDOW_ID" tiled
        last_pane_count=$pane_count
      fi
      if [ $pane_count == 0 ]; then
        break;
      fi
      sleep "$POLL_INTERVAL"
    done
    # Downloads that couldn't get a pane ran as plain background jobs.
    wait
  fi
}

trap cleanup SIGINT SIGTERM

# Echtvar-encoded gnomAD 4.1 population frequencies - this is a big one (~4GB)
ECHTVAR_FILE="gnomad_4.1_22297079_GRCh38_whole_genome"
start_download https://zenodo.org/records/22297079/files/gnomad_4.1_whole_genome.zip?download=1 "${ECHTVAR_FILE}" "\
Downloading Echtvar from https://zenodo.org/records/22297079"

# Monarch phenotype DB - another large download. 17GB decompressed
COMPRESSED_PHENIO="phenio.db.gz"
DECOMPRESSED_PHENIO="phenio.db"
if [ ! -f "${DECOMPRESSED_PHENIO}" ]; then
  log "[INFO] phenio.db not found, downloading..."
  start_download https://data.monarchinitiative.org/monarch-kg/latest/phenio.db.gz
else
  log "[INFO] phenio.db already exists, skipping download."
fi

# GRCh38 reference genome
GRCh38="ref.fa.gz"
GRCh38_decompressed="ref.fa"
# if it doesn't exist, download it
if [ ! -f "${GRCh38_decompressed}" ]; then
  log "[INFO] compressed ref.fa not found, attempting download..."
    GRCh38_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    start_download ${GRCh38_URL} ${GRCh38} "Downloading compressed GRCh38 reference genome from UCSC"
else
  log "[INFO] ${GRCh38_decompressed} already exists, skipping download."
fi

# MANE gene data
start_download https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.summary.txt.gz

# Ensembl GFF3 data
GFF3_ORIGINAL="Homo_sapiens.GRCh38.116.chr.gff3.gz"
GFF3_chrM_RENAMED="Homo_sapiens.GRCh38.116.MTtoM.gff3.gz"
start_download https://ftp.ensembl.org/pub/release-116/vertebrates/gff3/homo_sapiens/Homo_sapiens.GRCh38.116.chr.gff3.gz

# Jax lab file for phenotype matching
start_download https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-02-16/hp.obo

start_download https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-02-16/genes_to_phenotype.txt

# mitochondrial annotations
start_download https://www.fr.mitomap.org/downloads/mitotip_scores.txt

start_download "https://mitimpact.mcb2lab.org/cdn/nAPOGEE_v1.0.0.txt.zip"

start_download "https://mitimpact.mcb2lab.org/cdn/MitImpact_db_3.1.3.txt.zip"

# AlphaMissense raw data
AM="AlphaMissense_hg38.tsv.gz"
start_download "https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz?download=1" "${AM}"

# Structural Variant annotation resources
# SVAFotate population frequency BED - native GRCh38 gnomAD v4.1, merged with CCDG, TOPMed and 1000G.
# Contigs are Ensembl-style (1, not chr1) and must stay that way: SVAFotate strips the `chr` prefix from the
# query VCF but not from this file, so adding one here silently zeroes every frequency. See docs.
SVA="SVAFotate_SV_popAFs.GRCh38.v4.1.bed.gz"
SKINNY_SVAF="SVAFotate_reduced_gnomAD.bed.gz"
start_download "https://zenodo.org/records/11642574/files/SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz?download=1" "${SVA}"

# Non-coding elements BED, used by GATK SVAnnotate to add PREDICTED_NONCODING_* annotations.
# 4 columns and no header, despite SVAnnotate's help describing 6 columns with a header - it reads it fine.
# Contigs are UCSC-style, matching ref.fa. Left uncompressed, SVAnnotate has no gzip codec for it.
start_download "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/noncoding.sort.hg38.bed"

# MANE GTF, used by GATK SVAnnotate for gene consequences.
# SVAnnotate requires one transcript per gene, which MANE Select satisfies by construction.
# This is the same MANE release as the summary file above, so the gene symbols agree with RunSmallFilteringSv
start_download "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz"

THIS_MONTH=$(date '+%Y-%m')

submission_summary="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/submission_summary.txt.gz"
variant_summary="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

start_download $submission_summary "submissions_${THIS_MONTH}.txt.gz"
start_download $variant_summary "variants_${THIS_MONTH}.txt.gz"

await

# if the reduced SVAFotate input doesn't exist, generate it
# this is a reduction of the canonical SVAFotate reference BED file, removing:
# * any non-gnomAD rows
# * any ancestry-specific columns
# Loading the whole file exhausts memory on most VMs, so reduction in advance to only the fields we care about is ideal.
if [ ! -f "${SKINNY_SVAF}" ] && [ -f "${SVA}" ]; then
    log "Minimising SVAFotate reference BED"
    gunzip -c ${SVA} \
      | awk -F'\t' 'NR == 1 || $6 == "gnomAD"' \
      | cut -f1-21,176 \
      | gzip -c > ${SKINNY_SVAF}
fi

#if compressed exists, but decompressed doesn't, gunzip it
if [ ! -f "${GRCh38_decompressed}" ] && [ -f "${GRCh38}" ]; then
    log "Decompressing reference Fasta"
    gunzip ${GRCh38}
fi

# same for the phenio files
if [ ! -f "${DECOMPRESSED_PHENIO}" ] && [ -f "${COMPRESSED_PHENIO}" ]; then
    log "Decompressing Phenio.db"
    gunzip ${COMPRESSED_PHENIO}
fi

# rinse the Ensembl GFF3 file through a chrMT -> chrM rename
if [ ! -f "${GFF3_chrM_RENAMED}" ] && [ -f "${GFF3_ORIGINAL}" ]; then
    log "Renaming MT to M in ensembl GFF3 file"
    gunzip -c "${GFF3_ORIGINAL}" | sed 's/^MT/M/' | gzip > "${GFF3_chrM_RENAMED}"
fi

# Final status summary
summary_fail=0
script_name_summary="$(basename "$0")"
for target in "${DOWNLOAD_TARGETS[@]}"; do
  if [ ! -f "$target" ]; then
    log "[MISSING] $target"
    summary_fail=$((summary_fail+1))
  fi
done

if [ $summary_fail -eq 0 ]; then
  log "[SUCCESS] All ${#DOWNLOAD_TARGETS[@]} downloads completed successfully."
else
  log "[SUMMARY] $summary_fail of ${#DOWNLOAD_TARGETS[@]} downloads missing. Restart ${script_name_summary} to gracefully resume."
fi
