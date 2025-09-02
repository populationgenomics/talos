"""
demonstration in two parts for the creation of a db schema, and the creation of the db itself
the accompanying file create_sqlite_schema.py contains the schema definition
"""

from talos.db.create_sqlite_schema import (
    Base,
    Family,
    Participant,
    Decision,
    Trio,
    Variant,
    TranscriptConsequence,
    ReportEvent,
    Run,
)

from sqlalchemy import create_engine, select, join, update
from sqlalchemy.orm.session import Session

from cloudpathlib.anypath import to_anypath

_SESSION: Session | None = None


def main():
    """
    if the db exists, just connect to it. If this is the first time creation, create the tables.
    Returns:

    """
    db_name = 'test.db'
    # record whether the db already exists, or is a first time creation
    db_exists = to_anypath(db_name).exists()

    engine = create_engine('sqlite+pysqlite:///test.db')
    session = Session(engine)

    if not db_exists:
        print('creating new db tables from imported schema')
        Base.metadata.create_all(engine)

    global _SESSION
    _SESSION = session


main()

# some examples of how this would work
f1 = Family(id='FAM1')
_SESSION.add(f1)
p1 = Participant(
    id='P1',
    family=f1.id,
    ext_id='EXT1',
    affected=True,
)
_SESSION.add(p1)

pro1 = Trio(
    id=p1.id,
    solved=False,
)
_SESSION.add(pro1)

var1 = Variant(
    chromosome='chr1',
    position=123456,
    reference='A',
    alternate='T',
)
_SESSION.add(var1)

txcsq1 = TranscriptConsequence(
    variant=var1,
    gene_symbol='GENE1',
    consequence='missense_variant',
)
_SESSION.add(txcsq1)

# autoincrementing IDs are assigned on commit, so this needs to be pushed before the foreign key relationships can exist
_SESSION.commit()

reportevent1 = ReportEvent(
    trio_id=pro1.id,
    variant_id=var1.id,
    moi_satisfied='Initial population',
    panels={1, 2, 3},
    genotype_proband='0/1',
    categories={'alphamissense': '2025-02-14'},
)

_SESSION.add(reportevent1)

_SESSION.commit()

# query for the family, scalar_one() raises if not exactly one result
family_result = _SESSION.execute(select(Family).where(Family.id == 'FAM1')).scalar_one()
print(f'Family ID: {family_result.id}')

# join probands to participants, and filter by family id
result = _SESSION.execute(select(Trio, Participant).join(Trio.proband).where(Participant.family == family_result.id))

# iterate over the probands in the family
for participant in result:
    print(f'Participant ID: {participant.Participant.id}, EXT ID: {participant.Participant.ext_id}')
    for var in _SESSION.execute(
        select(ReportEvent).where(ReportEvent.trio_id == participant.Participant.id),
    ):
        print(f'  Reported variant ID: {var.ReportEvent.variant_id}')
        print(f'  Categories: {var.ReportEvent.categories}')

        # query back down to the transcript consequences - would be populated at runtime, and pulled to generate HTML
        for tx in _SESSION.execute(
            select(TranscriptConsequence).where(TranscriptConsequence.variant_id == var.ReportEvent.variant_id),
        ):
            print(
                f'    Gene: {tx.TranscriptConsequence.gene_symbol}, consequence: {tx.TranscriptConsequence.consequence}'
            )

        print(var.ReportEvent.categories)

        # now maybe I want to update the report event with an additional category
        if 'new_category' not in var.ReportEvent.categories:
            var.ReportEvent.add_categories(['new_category', 'another_category'], _SESSION)

# just to make sure it updated
reportevent_check = _SESSION.execute(select(ReportEvent).where(ReportEvent.id == reportevent1.id)).scalar_one()
print(f'Updated categories: {reportevent_check.categories}')
print(f'Updated date: {reportevent_check.evidence_last_updated}')

# record a new decision
decision = Decision(
    report_event_id=reportevent_check.id,
    decision='pathogenic',
    notes='I decided that this is bad',
)
_SESSION.add(decision)
_SESSION.commit()

# select it back
decision_check = _SESSION.execute(select(Decision).where(Decision.id == decision.id)).scalar_one()
print(f'Decision recorded: {decision_check.decision}, notes: {decision_check.notes}')


# sqlalchemy doesn't have a built-in get_or_create, but it's easy enough to do with a select and an add if not found
# https://stackoverflow.com/questions/2546207/does-sqlalchemy-have-an-equivalent-of-djangos-get-or-create

# get variants and corresponding report events in a single query
result = _SESSION.execute(select(Variant, ReportEvent).join(ReportEvent, Variant.id == ReportEvent.variant_id))
for var, report in result:
    print(f'Variant {var.id} reported in event {report.id}')

_SESSION.close()
