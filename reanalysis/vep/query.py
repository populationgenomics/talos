"""
Hail Query function to parse JSON VEP results.
"""

import hail as hl


def vep_json_to_ht(vep_results_paths, out_path):
    """
    Parse results from VEP with annotations formatted in JSON,
    and write into a Hail Table.
    """
    # Defining schema inside the function, so we can submit
    # the function to the Batch Job:
    schema = hl.dtype(
        'struct{assembly_name:str,allele_string:str,ancestral:str,'
        'colocated_variants:array<struct{aa_allele:str,aa_maf:float64,'
        'afr_allele:str,afr_maf:float64,allele_string:str,amr_allele:str,'
        'amr_maf:float64,clin_sig:array<str>,end:int32,eas_allele:str,'
        'eas_maf:float64,ea_allele:str,ea_maf:float64,eur_allele:str,'
        'eur_maf:float64,exac_adj_allele:str,exac_adj_maf:float64,'
        'exac_allele:str,exac_afr_allele:str,exac_afr_maf:float64,'
        'exac_amr_allele:str,exac_amr_maf:float64,exac_eas_allele:str,'
        'exac_eas_maf:float64,exac_fin_allele:str,exac_fin_maf:float64,'
        'exac_maf:float64,exac_nfe_allele:str,exac_nfe_maf:float64,'
        'exac_oth_allele:str,exac_oth_maf:float64,exac_sas_allele:str,'
        'exac_sas_maf:float64,id:str,minor_allele:str,minor_allele_freq:float64,'
        'phenotype_or_disease:int32,pubmed:array<int32>,sas_allele:str,'
        'sas_maf:float64,somatic:int32,start:int32,strand:int32}>,context:str,'
        'end:int32,id:str,input:str,intergenic_consequences:array<struct{'
        'allele_num:int32,consequence_terms:array<str>,impact:str,minimised:int32,'
        'variant_allele:str}>,most_severe_consequence:str,motif_feature_consequences'
        ':array<struct{allele_num:int32,consequence_terms:array<str>,high_inf_pos:str,'
        'impact:str,minimised:int32,motif_feature_id:str,motif_name:str,motif_pos:'
        'int32,motif_score_change:float64,strand:int32,variant_allele:str}>,'
        'regulatory_feature_consequences:array<struct{allele_num:int32,biotype:str,'
        'consequence_terms:array<str>,impact:str,minimised:int32,regulatory_feature_id:'
        'str,variant_allele:str}>,seq_region_name:str,start:int32,strand:int32,'
        'transcript_consequences:array<struct{allele_num:int32,amino_acids:str,'
        'appris:str,biotype:str,canonical:int32,mane_select:str,mane_plus_clinical:str,'
        'ccds:str,cdna_start:int32,cdna_end:int32,cds_end:int32,cds_start:int32,codons:'
        'str,consequence_terms:array<str>,distance:int32,domains:array<struct{db:str,'
        'name:str}>,exon:str,gene_id:str,gene_pheno:int32,gene_symbol:str,'
        'gene_symbol_source:str,hgnc_id:str,hgvsc:str,hgvsp:str,hgvs_offset:int32,'
        'impact:str,intron:str,lof:str,lof_flags:str,lof_filter:str,lof_info:str,'
        'minimised:int32,polyphen_prediction:str,polyphen_score:float64,protein_end:'
        'int32,protein_start:int32,protein_id:str,sift_prediction:str,sift_score:'
        'float64,strand:int32,swissprot:str,transcript_id:str,trembl:str,tsl:int32,'
        'uniparc:str,variant_allele:str}>,variant_class:str}'
    )
    ht = hl.import_table(vep_results_paths, no_header=True, types={'f0': schema})
    ht = ht.transmute(vep=ht.f0)
    # Can't use ht.vep.start for start because it can be modified by VEP (e.g. it
    # happens for indels). So instead parsing POS from the original VCF line stored
    # as ht.vep.input field.
    original_vcf_line = ht.vep.input
    start = hl.parse_int(original_vcf_line.split('\t')[1])
    chrom = ht.vep.seq_region_name
    ht = ht.annotate(locus=hl.locus(chrom, start))
    ht = ht.key_by(ht.locus)

    ht.write(str(out_path), overwrite=True)
