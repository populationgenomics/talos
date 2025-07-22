# Script for gathering large inputs required by Talos
# This runs a download from multiple different sources
# Output file names created/expected by this script match the initial configuration file content in `talos.config` and `annotation.config`

# Echtvar-encoded gnomAD 4.1 population frequencies - this is a big one (~6GB) so it's started early and backgrounded
ECHTVAR_FILE="gnomad_4.1_region_merged_GRCh38_whole_genome"
if [ ! -f ${ECHTVAR_FILE} ]; then
    echo "Downloading Echtvar from https://zenodo.org/records/15222100"
    curl https://zenodo.org/records/15222100/files/gnomad_4.1_region_merged_GRCh38_whole_genome?download=1 -o "${ECHTVAR_FILE}" &
else
    echo "${ECHTVAR_FILE} already exists"
fi

# Monarch phenotype DB - another large download
PHENIO_DB="phenio.db.gz"
if [ ! -f ${PHENIO_DB} ]; then
    curl https://data.monarchinitiative.org/monarch-kg/latest/phenio.db.gz -o "${PHENIO_DB}" &
else
    echo "${PHENIO_DB} already exists"
fi

# GRCh38 reference genome
GRCh38="ref.fa"
if [ ! -f ${GRCh38} ]; then
    echo "Downloading GRCh38 reference genome from the Broad Institute public data bucket"
    echo "n.b. this is a GRCh38 reference genome, and your data must also be aligned to this reference genome for Talos to work correctly"
    curl https://storage.cloud.google.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta -o "${GRCh38}"
else
    echo "${GRCh38} already exists"
fi

# MANE gene data
MANE="MANE.GRCh38.v1.4.summary.txt.gz"
if [ ! -f ${MANE} ]; then
    curl https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz -o "${MANE}"
else
    echo "${MANE} already exists"
fi

# Ensembl GFF3 data
GFF3="Homo_sapiens.GRCh38.113.gff3.gz"
if [ ! -f ${GFF3} ]; then
    curl https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh38.113.chr.gff3.gz -o "${GFF3}"
else
    echo "${GFF3} already exists"
fi

# Jax lab file for phenotype matching
OBO="hp.obo"
if [ ! -f ${OBO} ]; then
    curl https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2025-05-06/hp.obo -o "${OBO}"
else
    echo "${OBO} already exists"
fi

G2P="genes_to_phenotype.txt"
if [ ! -f ${G2P} ]; then
    curl https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2025-05-06/genes_to_phenotype.txt -o "${G2P}"
else
    echo "${G2P} already exists"
fi

# latest clinvarbitration data
CLINVAR="clinvarbitration.tar.gz"
if [ ! -f ${CLINVAR} ]; then
    curl https://zenodo.org/records/15896156/files/clinvarbitration_25-07_clinvar_decisions.release.tar.gz?download=1 -o "${CLINVAR}"
else
    echo "ClinvArbitration data already exists (${CLINVAR}), please check it's up to date (released monthly)"
fi

# AlphaMissense raw data
AM="AlphaMissense_hg38.tsv.gz"
if [ ! -f ${AM} ]; then
    curl https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz?download=1 -o "${AM}"
else
    echo "${AM} already exists"
fi
