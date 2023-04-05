"""
creates a data model for generation of validation results
"""


# pylint: disable=invalid-name,too-many-instance-attributes


import json
import os.path
from dataclasses import dataclass, field

import hail as hl

from reanalysis.utils import CustomEncoder


schema = hl.dtype(
    'struct{'
    'locus:str,alleles:array<str>,AC:int32,AF:float64,AN:int32,'
    'gnomad_genomes:struct{AF:float64,AN:int32,AC:int32,Hom:int32,Hemi:int32},'
    'gnomad_exomes:struct{AF:float64,AN:int32,AC:int32,Hom:int32,Hemi:int32},'
    'splice_ai:struct{delta_score:float64,splice_consequence:str},'
    'cadd:struct{PHRED:float64},'
    'dbnsfp:struct{REVEL_score:str,Mutationtaster_pred:str},'
    'clinvar:struct{clinical_significance:str,gold_stars:int32,allele_id:str},'
    'vep:struct{'
    'transcript_consequences:array<struct{gene_symbol:str,gene_id:str,'
    'variant_allele:str,consequence_terms:array<str>,transcript_id:str,'
    'protein_id:str,gene_symbol_source:str,canonical:int32,cdna_start:int32,'
    'cds_end:int32,biotype:str,protein_start:int32,protein_end:int32,'
    'sift_score:float64,sift_prediction:str,polyphen_score:float64,'
    'mane_select:str,lof:str}>,variant_class:str},geneIds:set<str>'
    '}'
)


@dataclass
class BaseFields:
    """
    base row fields
    """

    locus: str
    alleles: list[str]
    AC: int = field(default=1)
    AF: float = field(default=0.001)
    AN: int = field(default=1)


@dataclass
class AFGeneric:
    """
    generic Allele Frequency data model
    """

    AF: float = field(default=0.001)
    AN: int = field(default=1)
    AC: int = field(default=1)
    Hom: int = field(default=0)


@dataclass
class AFData:
    """
    specific Allele Frequency data model
    """

    gnomad_genomes: AFGeneric = field(default_factory=AFGeneric)
    gnomad_exomes: AFGeneric = field(default_factory=AFGeneric)


@dataclass
class Splice:
    """
    Splice data model
    """

    delta_score: float = field(default=0.01)
    splice_consequence: str = field(default='none')


@dataclass
class CADD:
    """
    CADD data model
    """

    PHRED: float = field(default=0.01)


@dataclass
class DBnsfp:
    """
    DBnsfp data model
    """

    REVEL_score: str = field(default='0.0')
    Mutationtaster_pred: str = field(default='n')


@dataclass
class Clinvar:
    """
    Clinvar data model
    """

    clinical_significance: str = field(default_factory=str)
    gold_stars: int = field(default_factory=int)
    allele_id: str = field(default_factory=str)


@dataclass
class TXFields:
    """
    TX fields data model
    """

    gene_symbol: str
    gene_id: str
    variant_allele: str = field(default_factory=str)
    consequence_terms: list = field(default_factory=list)
    transcript_id: str = field(default_factory=str)
    protein_id: str = field(default_factory=str)
    gene_symbol_source: str = field(default_factory=str)
    canonical: int = field(default=1)
    cdna_start: int = field(default=1)
    cds_end: int = field(default=1)
    biotype: str = field(default_factory=str)
    protein_start: int = field(default=1)
    protein_end: int = field(default=1)
    sift_score: int = field(default=1.0)  # lowest possible score
    sift_prediction: str = field(default_factory=str)
    polyphen_score: float = field(default=0.01)
    mane_select: str = field(default_factory=str)
    lof: str = field(default_factory=str)


@dataclass
class VepVariant:
    """
    class object to sweep up all the data models
    """

    def __init__(
        self,
        tx: list[TXFields],
        base: BaseFields,
        af: AFData | None = None,
        cadd: CADD | None = None,
        dbnsfp: DBnsfp | None = None,
        clinvar: Clinvar | None = None,
        splice: Splice | None = None,
        var_class: str = 'SNV',
    ):
        """
        VEP variant data model
        """
        self.data = (
            {
                'vep': {'transcript_consequences': tx, 'variant_class': var_class},
                'geneIds': {tx.gene_id for tx in tx},
                'dbnsfp': dbnsfp or DBnsfp(),
                'cadd': cadd or CADD(),
                'clinvar': clinvar or Clinvar(),
                'splice_ai': splice or Splice(),
            }
            | base.__dict__
            | (af or AFData()).__dict__
        )

    def to_string(self) -> str:
        """
        convert this object to a string
        Returns:
            a single line JSON string
        """
        return json.dumps(self.data, cls=CustomEncoder) + '\n'


class SneakyTable:
    """
    class to take multiple individual variants
    and generate a Hail Matrix Table from them
    """

    def __init__(self, variants: list[VepVariant], tmp_path: str):
        """
        Args:
            variants (list[VepVariant]): list of VepVariant objects
            tmp_path (str): where to write the hail table
        """
        self.variants = variants
        self.tmp_path = tmp_path

    def to_hail(self) -> hl.Table:
        """
        write the data model to a hail table

        Returns:
            the table of the faux annotation data
        """

        hl.init(default_reference='GRCh38')

        # write this object
        json_temp = os.path.join(self.tmp_path, 'vep.json')
        with open(json_temp, 'w', encoding='utf-8') as f:
            for variant in self.variants:
                f.write(variant.to_string())

        # read it back in as a hail table
        # field must be f0 if no header
        ht = hl.import_table([json_temp], no_header=True, types={'f0': schema})

        # unwrap the annotation data
        ht = ht.transmute(**ht.f0)

        # transmute the locus and alleles, set as keys
        ht = ht.transmute(locus=hl.parse_locus(ht.locus), alleles=ht.alleles)
        ht = ht.key_by('locus', 'alleles')

        # checkpoint out to a local dir
        ht.write(os.path.join(self.tmp_path, 'vep.ht'), overwrite=True)

        # send it
        return ht
