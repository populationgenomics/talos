"""
creates a data model for generation of validation results
"""


# pylint: disable=invalid-name,too-many-instance-attributes


from dataclasses import dataclass, field
import hail as hl


# from reanalysis.utils import CustomEncoder


field_of_1 = field(default=1)
field_of_0 = field(default=0)
field_of_1_percent = field(default=0.01)
field_of_milli = field(default=0.001)
field_of_a_blank_string = field(default=' ')


@dataclass
class BaseFields:
    """
    base row fields
    """

    locus: str
    alleles: list[str]
    AC: int = field_of_1
    AF: float = field_of_milli
    AN: int = field_of_1


@dataclass
class AFGeneric:
    """
    generic Allele Frequency data model
    """

    AF: float = field_of_milli
    AN: int = field_of_1
    AC: int = field_of_1
    Hom: int = field_of_0


@dataclass
class AFData:
    """
    specific Allele Frequency data model
    """

    genomes: AFGeneric = field(default=AFGeneric)
    exomes: AFGeneric = field(default=AFGeneric)


@dataclass
class Splice:
    """
    Splice data model
    """

    delta_score: float = field_of_1_percent
    splice_consequence: str = field(default='none')


@dataclass
class CADD:
    """
    CADD data model
    """

    PHRED: float = field_of_1_percent


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
    allele_id: str = field(default=field_of_a_blank_string)


@dataclass
class TXFields:
    """
    TX fields data model
    """

    gene_symbol: str
    gene_id: str
    variant_allele: str = field_of_a_blank_string
    consequence_terms: list = field(default_factory=list)
    transcript_id: str = field_of_a_blank_string
    protein_id: str = field_of_a_blank_string
    gene_symbol_source: str = field_of_a_blank_string
    canonical: int = field_of_1
    cdna_start: int = field_of_1
    cds_end: int = field_of_1
    biotype: str = field_of_a_blank_string
    protein_start: int = field_of_1
    protein_end: int = field_of_1
    sift_score: int = field_of_1  # lowest possible score
    sift_prediction: str = field_of_a_blank_string
    polyphen_score: float = field_of_1_percent
    mane_select: str = field_of_a_blank_string
    lof: str = field_of_a_blank_string


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
    ):
        """
        VEP variant data model
        """
        self.tx = tx
        self.base = base
        self.af = af or AFData()
        self.cadd = cadd or CADD()
        self.dbnsfp = dbnsfp or DBnsfp()
        self.clinvar = clinvar or Clinvar()
        self.splice = splice or Splice()

    def to_hail(self, table_path: str):
        """
        write the data model to a hail table
        Args:
            table_path ():

        Returns:

        """

        # write this object


schema = hl.dtype(
    'struct{locus:str,alleles:array<str>,AC:int32,AF:float64,AN:int32,'
    'gnomad_genomes:struct{AF:float64,AN:int32,AC:int32,Hom:int32,Hemi:int32},'
    'gnomad_exomes:struct{AF:float64,AN:int32,AC:int32,Hom:int32,Hemi:int32},'
    'splice_ai:struct{delta_score:float64,splice_consequence:str},'
    'cadd:struct{PHRED:float64},'
    'dbnsfp:struct{REVEL_score:str,Mutationtaster_pred:str},'
    'clinvar:struct{clinical_significance:str,gold_stars:int32,allele_id:str},'
    'vep:struct{transcript_consequences:array<struct{gene_symbol:str,gene_id:str,'
    'variant_allele:str,consequence_terms:array<str>,transcript_id:str,'
    'protein_id:str,gene_symbol_source:str,canonical:int32,cdna_start:int32,'
    'cds_end:int32,biotype:str,protein_start:int32,protein_end:int32,'
    'sift_score:float64,sift_prediction:str,mane_select:str,lof:str}>,'
    'most_severe_consequence:str}}'
)
