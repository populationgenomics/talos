"""
Query function to parse JSON VEP results.
"""

import hail as hl


def vep_json_to_ht(vep_result_paths: list[str], out_path: str):
    """
    Parse results from VEP with annotations formatted in JSON,
    and write into a Hail Table.
    """
    # Defining schema inside the function, so we can submit
    # the function to the Batch Job:
    json_schema = hl.dtype(
        """struct{
        minimised:int32,
        assembly_name:str,
        allele_string:str,
        ancestral:str,
        colocated_variants:array<struct{
            allele_string:str,
            clin_sig:array<str>,
            clin_sig_allele:str,
            end:int32,
            id:str,
            minimised:int32,
            minor_allele:str,
            minor_allele_freq:float64,
            phenotype_or_disease:int32,
            pubmed:array<int32>,
            seq_region_name:str,
            somatic:int32,
            start:int32,
            strand:int32
        }>,
        context:str,
        end:int32,
        id:str,
        input:str,
        intergenic_consequences:array<struct{
            allele_num:int32,
            consequence_terms:array<str>,
            impact:str,minimised:int32,
            variant_allele:str
        }>,
        most_severe_consequence:str,
        motif_feature_consequences:array<struct{
            allele_num:int32,
            consequence_terms:array<str>,
            high_inf_pos:str,
            impact:str,
            minimised:int32,
            motif_feature_id:str,
            motif_name:str,
            motif_pos:int32,
            motif_score_change:float64,
            strand:int32,
            transcription_factors:array<str>,
            variant_allele:str
        }>,
        regulatory_feature_consequences:array<struct{
            allele_num:int32,
            biotype:str,
            consequence_terms:array<str>,
            impact:str,
            minimised:int32,
            regulatory_feature_id:str,
            variant_allele:str
        }>,
        seq_region_name:str,
        start:int32,
        strand:int32,
        transcript_consequences:array<struct{
            allele_num:int32,
            amino_acids:str,
            appris:str,
            biotype:str,
            canonical:int32,
            mane_select:str,
            mane_plus_clinical:str,
            ccds:str,
            cdna_start:int32,
            cdna_end:int32,
            cds_end:int32,
            cds_start:int32,
            codons:str,
            consequence_terms:array<str>,
            distance:int32,
            domains:array<struct{
                db:str,
                name:str
            }>,
            exon:str,
            gene_id:str,
            gene_pheno:int32,
            gene_symbol:str,
            gene_symbol_source:str,
            hgnc_id:str,
            hgvsc:str,
            hgvsp:str,
            hgvs_offset:int32,
            impact:str,
            intron:str,
            lof:str,
            lof_flags:str,
            lof_filter:str,
            lof_info:str,
            minimised:int32,
            mirna:array<str>,
            polyphen_prediction:str,
            polyphen_score:float64,
            protein_end:int32,
            protein_start:int32,
            protein_id:str,
            sift_prediction:str,
            sift_score:float64,
            strand:int32,
            swissprot:array<str>,
            transcript_id:str,
            trembl:array<str>,
            tsl:int32,
            uniparc:array<str>,
            uniprot_isoform:array<str>,
            variant_allele:str,
            am_class:str,
            am_pathogenicity:float64,
            source:str,
            flags:array<str>
        }>,
        variant_class:str
    }"""
    )
    ht = hl.import_table(
        paths=vep_result_paths, no_header=True, types={'f0': json_schema}
    )
    ht = ht.transmute(vep=ht.f0)

    # Can't use ht.vep.start for start because it can be modified by VEP (e.g. it
    # happens for indels). So instead parsing POS from the original VCF line stored
    # as ht.vep.input field.
    original_vcf_line = ht.vep.input
    start = hl.parse_int(original_vcf_line.split('\t')[1])
    chrom = ht.vep.seq_region_name
    ht = ht.annotate(
        locus=hl.locus(chrom, start),
        alleles=[ht.vep.input.split('\t')[3], ht.vep.input.split('\t')[4]],
    )
    ht = ht.key_by(ht.locus, ht.alleles)
    ht.write(str(out_path), overwrite=True)
