"""
Query function to parse JSON VEP results.
"""

import hail as hl


def vep_json_to_ht(vep_results_paths, out_path):
    """
    Parse results from VEP with annotations formatted in JSON,
    and write into a Hail Table.
    """
    # Defining schema inside the function, so we can submit
    # the function to the Batch Job:
    json_schema = (
        'Struct{assembly_name:String,allele_string:String,ancestral:String,'
        'colocated_variants:Array[Struct{aa_allele:String,aa_maf:Float64,'
        'afr_allele:String,afr_maf:Float64,allele_string:String,'
        'amr_allele:String,amr_maf:Float64,clin_sig:Array[String],'
        'end:Int32,eas_allele:String,eas_maf:Float64,ea_allele:String,'
        'ea_maf:Float64,eur_allele:String,eur_maf:Float64,'
        'exac_adj_allele:String,exac_adj_maf:Float64,exac_allele:String,'
        'exac_afr_allele:String,exac_afr_maf:Float64,'
        'exac_amr_allele:String,exac_amr_maf:Float64,'
        'exac_eas_allele:String,exac_eas_maf:Float64,'
        'exac_fin_allele:String,exac_fin_maf:Float64,exac_maf:Float64,'
        'exac_nfe_allele:String,exac_nfe_maf:Float64,'
        'exac_oth_allele:String,exac_oth_maf:Float64,'
        'exac_sas_allele:String,exac_sas_maf:Float64,id:String,'
        'minor_allele:String,minor_allele_freq:Float64,'
        'phenotype_or_disease:Int32,pubmed:Array[Int32],sas_allele:String,'
        'sas_maf:Float64,somatic:Int32,start:Int32,strand:Int32}],'
        'context:String,end:Int32,id:String,input:String,'
        'intergenic_consequences:Array[Struct{allele_num:Int32,'
        'consequence_terms:Array[String],impact:String,minimised:Int32,'
        'variant_allele:String}],most_severe_consequence:String,'
        'motif_feature_consequences:Array[Struct{allele_num:Int32,'
        'consequence_terms:Array[String],high_inf_pos:String,impact:String,'
        'minimised:Int32,motif_feature_id:String,motif_name:String,'
        'motif_pos:Int32,motif_score_change:Float64,strand:Int32,'
        'variant_allele:String}],regulatory_feature_consequences:Array['
        'Struct{allele_num:Int32,biotype:String,consequence_terms:Array['
        'String],impact:String,minimised:Int32,'
        'regulatory_feature_id:String,variant_allele:String}],'
        'seq_region_name:String,start:Int32,strand:Int32,'
        'transcript_consequences:Array[Struct{allele_num:Int32,'
        'amino_acids:String,appris:String,biotype:String,canonical:Int32,'
        'mane_select:String,mane_plus_clinical:String,ccds:String,'
        'cdna_start:Int32,cdna_end:Int32,cds_end:Int32,cds_start:Int32,'
        'codons:String,consequence_terms:Array[String],distance:Int32,'
        'domains:Array[Struct{db:String,name:String}],exon:String,'
        'gene_id:String,gene_pheno:Int32,gene_symbol:String,'
        'gene_symbol_source:String,hgnc_id:String,hgvsc:String,'
        'hgvsp:String,hgvs_offset:Int32,impact:String,intron:String,'
        'lof:String,lof_flags:String,lof_filter:String,lof_info:String,'
        'minimised:Int32,polyphen_prediction:String,polyphen_score:Float64,'
        'protein_end:Int32,protein_start:Int32,protein_id:String,'
        'sift_prediction:String,sift_score:Float64,strand:Int32,'
        'swissprot:String,transcript_id:String,trembl:String,tsl:Int32,'
        'uniparc:String,variant_allele:String}],variant_class:String}'
    )
    for k, v in [
        ('String', 'str'),
        ('Array', 'array'),
        ('Set', 'set'),
        ('Tuple', 'tuple'),
        ('Struct', 'struct'),
        ('[', '<'),
        (']', '>'),
        ('Int32', 'int32'),
        ('Int64', 'int64'),
        ('Float64', 'float64'),
        ('Float32', 'float32'),
    ]:
        json_schema = json_schema.replace(k, v)

    schema = hl.dtype(json_schema)
    ht = hl.import_table(vep_results_paths, no_header=True, types={'f0': schema})
    ht = ht.transmute(vep=ht.f0)
    # Can't use ht.vep.start for start because it can be modified by VEP (e.g. it
    # happens for indels). So instead parsing POS from the original VCF line stored
    # as ht.vep.input field.
    start = hl.parse_int(ht.vep.input.split('\t')[1])
    chrom = ht.vep.seq_region_name
    ht = ht.annotate(
        locus=hl.locus(chrom, start),
        alleles=[ht.vep.input.split('\t')[3], ht.vep.input.split('\t')[4]],
    )
    ht = ht.key_by(ht.locus, ht.alleles)
    ht.write(str(out_path), overwrite=True)
