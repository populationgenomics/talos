from sqlalchemy import update
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship
from sqlalchemy.orm.session import Session
from sqlalchemy.schema import ForeignKey, UniqueConstraint, Index
from sqlalchemy.types import Boolean, Float, Integer, PickleType, String

from talos.static_values import get_granular_date


class Base(DeclarativeBase):
    pass


class Family(Base):
    """
    a family ID, used as a foreign key to collect all family members, nothing major here, just simplifies queries
    """

    __tablename__ = 'family'
    id: Mapped[str] = mapped_column(String(50), primary_key=True, unique=True)


class Participant(Base):
    """
    Purpose: atomic description of each participant:
      - no relationships to other members
      - can be used to collect external IDs, affected status, HPO terms, etc.
      - we could add name, age, etc. but these aren't important for Talos

    The "trio" table (previously "probands", I don't like either name) has proband~mother~father relationships
    """

    __tablename__ = 'participant'

    # the sample ID we would expect in the VCF
    id: Mapped[str] = mapped_column(String(50), primary_key=True, unique=True, index=True)

    # these two attributes are probably also unique, but not primary keys as they can be null
    ext_id: Mapped[str | None] = mapped_column(String(50), nullable=True)
    # seqr ID is actually at the family level, but that's not important here, it still applies to the proband
    seqr_id: Mapped[str | None] = mapped_column(String(50), nullable=True)

    # each participant belongs to a family, non-nullable foreign key. This is mandatory from the Pedigree
    family: Mapped[str] = mapped_column(String(50), ForeignKey('family.id'), nullable=False)

    # boolean
    affected: Mapped[bool] = mapped_column(nullable=False, default=False)

    # embedded python set type, nullable. to be populated from the pedigree
    hpo_terms: Mapped[set[str] | None] = mapped_column(PickleType, nullable=True)


class Trio(Base):
    """
    Purpose: Affected members and immediate family context
      - Each row represents a nuclear trio used in inheritance filtering
      - This is centered around the affected participant, who may or may not have parents in the system
      - The proband is the affected participant, mother and father are optional
      - We can have multiple trios in the same family (i.e. affected sibs), each being a separate entity here
      - This would be more important in the report/presentation stages
          - retrieve all members and affected status which played a role in the filtering decision
          - present genotypes for the trio members

    TODO I don't like the name of this table - the family unit acted on by Talos is not always a trio, but is limited
    TODO to proband, mother, father.
    """

    __tablename__ = 'trio'

    # multiple foreign keys to the participant table, we need to explicitly define the foreign_keys in the relationship
    id: Mapped[str] = mapped_column(String(50), ForeignKey('participant.id'), primary_key=True, unique=True, index=True)
    proband = relationship('Participant', foreign_keys=[id])
    mother_id: Mapped[str | None] = mapped_column(String(50), ForeignKey('participant.id'), nullable=True)
    mother = relationship('Participant', foreign_keys=[mother_id])
    father_id: Mapped[str | None] = mapped_column(String(50), ForeignKey('participant.id'), nullable=True)
    father = relationship('Participant', foreign_keys=[father_id])

    # at the affected participant/proband level, we can track whether the case has been solved
    # this is discrete from the particpant level affected status, and the family (could have multiple affected members)
    # this can be accompanied by a solved_by foreign key to a decision - the variant plus notes
    solved: Mapped[bool | None] = mapped_column(Boolean, default=False)
    solved_by: Mapped[int | None] = mapped_column(Integer, ForeignKey('decision.id'), nullable=True)

    # a timestamp if the family existed, and was updated to contain additional members
    family_updated: Mapped[str | None] = mapped_column(String(50), nullable=True)


class Variant(Base):
    """
    Purpose: an entity to store the details of a given variant, participant/trio agnostic
      - the coordinates and alleles defining the variant
      - static annotations which do not change between participants, e.g., gnomAD, ClinVar, SpliceAI, etc.
      - relationships to transcript consequences, which may be multiple per variant

    we do not store per-event information here, e.g. categories, genotypes, that is in the ReportEvent table
    """

    __tablename__ = 'variant'

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True, nullable=False)
    chromosome: Mapped[str] = mapped_column(String(5), nullable=False)
    position: Mapped[int] = mapped_column(Integer, nullable=False)
    reference: Mapped[str] = mapped_column(String(500), nullable=False)
    alternate: Mapped[str] = mapped_column(String(500), nullable=False)

    # small, str, mito, cnv, sv, etc.
    var_type: Mapped[str | None] = mapped_column(String(20), default='small')

    # relatively fixed annotations. Not currently checking multiple discrete populations, though we could.
    gnomad_af: Mapped[float | None] = mapped_column(Float, nullable=True)
    gnomad_ac: Mapped[int | None] = mapped_column(Integer, nullable=True)
    gnomad_an: Mapped[int | None] = mapped_column(Integer, nullable=True)
    gnomad_homalt: Mapped[int | None] = mapped_column(Integer, nullable=True)

    # multiple foreign key connections to the TX consequences for this variant
    transcript_consequence: Mapped[list['TranscriptConsequence']] = relationship(back_populates='variant')

    # SpliceAI can be 4 separate scores, but currently we only have one (or none by default)
    # we're currently getting SpliceAI scores at the variant level, not the transcript level
    spliceai: Mapped[float | None] = mapped_column(Float, nullable=True)

    # clinvar related attributes - significance, stars, date of review status increase (all nullable)
    clinvar_clinsig: Mapped[str | None] = mapped_column(String(20), nullable=True)
    clinvar_stars: Mapped[int | None] = mapped_column(Integer, nullable=True)
    clinvar_increased: Mapped[str | None] = mapped_column(String(20), nullable=True)

    # specific to SVs - potentially multiple implicated genes
    predicted_lof: Mapped[set[str] | None] = mapped_column(PickleType, nullable=True)

    # ensure uniqueness on the combination [chromosome, position, reference, alternate], use as an index
    __table_args__ = (Index('var_index', 'chromosome', 'position', 'reference', 'alternate', unique=True),)

    def __repr__(self) -> str:
        return f'{self.chromosome}-{self.position}-{self.reference}-{self.alternate}'


class TranscriptConsequence(Base):
    """
    Each transcript consequence, many-to-one:variant relationship.
    We could add more fields if we want to track them.
    """

    __tablename__ = 'transcript_consequence'

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    variant_id: Mapped[int] = mapped_column(ForeignKey('variant.id'), index=True)
    variant: Mapped['Variant'] = relationship(back_populates='transcript_consequence')

    # broadly this matches the fields obtained from bcftools, VEP, and/or ANNOVAR
    gene_symbol: Mapped[str | None] = mapped_column(String(50), nullable=True)
    ensg: Mapped[str | None] = mapped_column(String(50), nullable=True)
    enst: Mapped[str | None] = mapped_column(String(50), nullable=True)
    ensp: Mapped[str | None] = mapped_column(String(50), nullable=True)
    dna_change: Mapped[str | None] = mapped_column(String(50), nullable=True)
    amino_acid_change: Mapped[str | None] = mapped_column(String(50), nullable=True)
    codon: Mapped[int | None] = mapped_column(Integer, nullable=True)
    exon: Mapped[int | None] = mapped_column(Integer, nullable=True)
    mane_id: Mapped[str | None] = mapped_column(String(20), nullable=True)
    mane_status: Mapped[str | None] = mapped_column(String(20), nullable=True)

    alphamissense_class: Mapped[str | None] = mapped_column(String(50), nullable=True)
    alphamissense_pathogenicity: Mapped[float | None] = mapped_column(Float, nullable=True)

    # not including alphagenome yet, it has too many scores and too much variability, requires more thought

    # the bcftools/VEP consequence term(s), comma-separated if multiple
    consequence: Mapped[str | None] = mapped_column(String(50), nullable=True)

    # contains the clinvar PM5 details, if applicable
    clinvar_pm5: Mapped[str | None] = mapped_column(String(100), nullable=True)


class ReportEvent(Base):
    """
    A reportable event; a variant, the family it was seen in, genotypes, and categories.
    We can revisit the same event later, if it is seen with different (additional) categories.
    """

    __tablename__ = 'report_event'

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    timestamp: Mapped[str] = mapped_column(String(50), nullable=False, default=get_granular_date())
    trio_id: Mapped[str] = mapped_column(String(50), ForeignKey('trio.id'), nullable=False)
    variant_id: Mapped[str] = mapped_column(String(50), ForeignKey('variant.id'), nullable=False)

    # optional, for comp-het pairs
    second_variant_id: Mapped[str | None] = mapped_column(String(50), ForeignKey('variant.id'), nullable=True)

    # the MOI model which was satisfied for this event, singular, as we don't merge variants any more
    moi_satisfied: Mapped[str] = mapped_column(String(20), nullable=False)

    # optional, any filters or AB test failures which should be noted for the variant
    flags: Mapped[set[str] | None] = mapped_column(PickleType, nullable=True)

    # a collection of all the panels applied to this trio containing this gene
    panels: Mapped[set[int] | None] = mapped_column(PickleType, nullable=True)

    # # splitting the nuclear family genotypes out into separate fields. This is a bit limiting, but simpler than a dict
    genotype_proband: Mapped[str] = mapped_column(String(10), nullable=False)
    genotype_mother: Mapped[str | None] = mapped_column(String(10), nullable=True)
    genotype_father: Mapped[str | None] = mapped_column(String(10), nullable=True)

    # todo store depth and AB info?

    # a dictionary of categories, and dates of first discovery, can be updated.
    # this format is something SQLalchemy can handle fine, but requires the field being replaced with each change.
    categories: Mapped[dict[str, str]] = mapped_column(PickleType, nullable=False)

    callset_af: Mapped[float | None] = mapped_column(Float, nullable=True)
    callset_ac: Mapped[int | None] = mapped_column(Integer, nullable=True)
    callset_an: Mapped[int | None] = mapped_column(Integer, nullable=True)

    # collect HPO terms and/or panels that matched this variant gene to the proband
    forced_match_panel: Mapped[set[str] | None] = mapped_column(PickleType, nullable=True)
    pheno_match_hpo: Mapped[set[str] | None] = mapped_column(PickleType, nullable=True)
    pheno_match_panels: Mapped[set[str] | None] = mapped_column(PickleType, nullable=True)

    # date of first phenomatch discovery, if any
    first_pheno_match: Mapped[str | None] = mapped_column(String(50), nullable=True)
    evidence_last_updated: Mapped[str | None] = mapped_column(String(50), nullable=True)

    __table_args__ = (Index('variant_moi_unique', 'variant_id', 'second_variant_id', 'moi_satisfied', unique=True),)

    def add_categories(self, category: str | list[str], session: Session) -> None:
        """
        Add one or more new categories with the current date if it doesn't already exist.
        Example of chucking this update logic into the model rather than the application code.
        Dubious about whether this would cause longer runtimes than a single mass update at the end.
        """
        # make it iterable
        categories = category if isinstance(category, list) else [category]
        updated = False
        for cat in categories:
            if cat not in self.categories:
                updated = True
                self.categories[cat] = get_granular_date()

        if not updated:
            return
        session.execute(
            update(ReportEvent),
            [
                {
                    'evidence_last_updated': max(self.categories.values()),
                    'categories': self.categories,
                    'id': self.id,
                },
            ],
        )
        session.commit()


class Decision(Base):
    """
    A record of a decision made on a report event, with a timestamp and optional notes.
    In case of an interactive application, the use-case here would be to track changes to decisions over time.
    A user could select the variant event, make a decision, and add notes.
    We could add a user ID if we want to track who made the decision, though that would require a user table and login.
    """

    __tablename__ = 'decision'

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)

    # the report event table connects the trio, variant, genotypes, categories, etc.
    report_event_id: Mapped[int] = mapped_column(Integer, ForeignKey('report_event.id'), nullable=False)

    timestamp: Mapped[str] = mapped_column(String(50), default=get_granular_date())

    # e.g., 'causative', 'phenotype expansion', 'artefact'. Could be an enum if we want to restrict values.
    decision: Mapped[str] = mapped_column(String(50), nullable=False)

    # extended opportunity for notes, nullable. This is analogous to the Seqr notes field.
    notes: Mapped[str | None] = mapped_column(String(500), nullable=True)


class Run(Base):
    """
    A record of Talos runs on a cohort, with parameters and a timestamp.
    We could add a JSON field for parameters if we want to track them.
    We could also tag each first-time variant appearance with the run ID.
    """

    __tablename__ = 'run'

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    timestamp: Mapped[str] = mapped_column(String(50), default=get_granular_date())
    version: Mapped[str] = mapped_column(String(50), nullable=False)
