from talos.HPOFlagging import find_genes_in_these_results
from talos.models import Coordinates, ParticipantMeta, ParticipantResults, ReportVariant, ResultData, SmallVariant

TEST_COORDS = Coordinates(chrom='1', pos=1, ref='A', alt='C')
SMALL_1 = SmallVariant(
    info={'gene_id': 'ENSG1'},
    coordinates=TEST_COORDS,
    transcript_consequences=[],
)
SV_1 = SmallVariant(
    info={'gene_id': 'ENSG2', 'lof': 'ENSG3,ENSG4'},
    coordinates=TEST_COORDS,
    transcript_consequences=[],
)


def test_find_genes_in_these_results():
    pm = ParticipantMeta(ext_id='male', family_id='family_1')
    test_results = ResultData(
        results={
            'male': ParticipantResults(metadata=pm, variants=[ReportVariant(sample='male', var_data=SMALL_1)]),
            'female': ParticipantResults(metadata=pm, variants=[ReportVariant(sample='male', var_data=SV_1)]),
        },
    )
    identified = find_genes_in_these_results(test_results)
    assert identified == {'ENSG1', 'ENSG2', 'ENSG3', 'ENSG4'}
