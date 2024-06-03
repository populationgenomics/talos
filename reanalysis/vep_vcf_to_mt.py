import hail as hl

csq = [
    'af',
    'allele',
    'allele_num',
    'am_class',
    'am_pathogenicity',
    'amino_acids',
    'amr_af',
    'biotype',
    'cdna_position',
    'cds_position',
    'codons',
    'consequence',
    'ensp',
    'exon',
    'feature',
    'feature_type',
    'gene',
    'gene_pheno',
    'gnomade_af',
    'gnomadg_af',
    'hgnc_id',
    'hgvs_offset',
    'hgvsc',
    'hgvsp',
    'impact',
    'lof',
    'mane_plus_clinical',
    'mane_select',
    'max_af',
    'polyphen',
    'protein_position',
    'sift',
    'symbol',
    'symbol_source',
    'variant_class',
]
RELEVANT_FIELDS = [
    'af',
    'allele',
    'am_class',
    'am_pathogenicity',
    'biotype',
    'canonical',
    'cdna_position',
    'cdna_start',
    'cds_end',
    'cds_position',
    'cds_start',
    'codons',
    'consequence',
    'ensp',
    'exon',
    'feature',
    'feature_type',
    'gene',
    'gene_id',
    'gene_symbol',
    'gene_symbol_source',
    'gnomade_af',
    'gnomadg_af',
    'hgvsc',
    'hgvsp',
    'impact',
    'lof',
    'mane_select',
    'max_af',
    'protein_end',
    'protein_id',
    'protein_position',
    'protein_start',
    'symbol',
    'symbol_source',
    'transcript_id',
    'variant_class',
]
VCF = 'part3.vcf.gz'

hl.init()
hl.default_reference('GRCh38')

mt = hl.import_vcf(VCF, array_elements_required=False, force_bgz=True)

# get the '|'-delimited String of all header names
headers = hl.get_vcf_metadata(VCF)
csq_string = headers['info']['CSQ']['Description'].split('Format: ').lower()[-1].split('|')
# get the CSQ contents as a list of lists of strings, per variant
split_csqs = mt.info.CSQ.map(lambda csq_entry: csq_entry.split('\|'))
# take those lists of lists and build them into a list of structs based on the known fields
mt = mt.annotate_rows(
    info=mt.info.annotate(
        CSQ=split_csqs.map(
            lambda x: hl.struct(
                **{csq_string[n]: x[n] for n in range(len(csq_string)) if csq_string[n] in RELEVANT_FIELDS},
            ),
        ),
    ),
)
mt = mt.explode_rows(mt.info.CSQ)
mt.describe()
mt.rows().show(5)

# expl = mt.annotate_rows(csq=mt.info.CSQ.map(lambda csq_entry: csq_entry.split('\|')))
# expl.describe()
# expl = expl.explode_rows(expl.csq)
# expl.describe()
#
# expl = expl.annotate_rows(csq_struct=hl.struct(**{csq_string[n]: expl.csq[n] for n in range(len(csq_string))}))
# expl.rows().show()
