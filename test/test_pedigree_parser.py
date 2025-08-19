from unittest.mock import mock_open, patch

import pytest

from talos.pedigree_parser import PedigreeParser

PED_CONTENT = 'FAM1\tS1\t0\t0\t1\t2\tHP:0000118,HP:0001250\nFAM1\tS2\t0\t0\t2\t1\t\nFAM1\tS3\tS1\tS2\t1\t1\t\n'


@patch('talos.pedigree_parser.to_anypath')
def test_read_pedigree(mock_to_anypath):
    mock_to_anypath.return_value.open = mock_open(read_data=PED_CONTENT)
    parser = PedigreeParser('dummy_path')
    participants = parser.participants
    assert 'S1' in participants
    assert participants['S1'].affected == 2
    assert participants['S1'].is_affected
    assert participants['S1'].hpo_terms == {'HP:0000118', 'HP:0001250'}
    assert participants['S2'].sex == 2
    assert participants['S2'].is_female
    assert participants['S3'].father_id == 'S1'
    assert participants['S3'].mother_id == 'S2'


@patch('talos.pedigree_parser.to_anypath')
def test_read_pedigree_and_strip(mock_to_anypath):
    mock_to_anypath.return_value.open = mock_open(read_data=PED_CONTENT)
    parser = PedigreeParser('dummy_path')
    participants = parser.strip_pedigree_to_samples(['S1'])
    assert 'S1' in participants
    assert 'S2' not in participants
    assert {'S1'} == set(participants.keys())


@patch('talos.pedigree_parser.to_anypath')
def test_get_affected_members(mock_to_anypath):
    mock_to_anypath.return_value.open = mock_open(read_data=PED_CONTENT)
    parser = PedigreeParser('dummy_path')
    affected = parser.get_affected_members()
    assert 'S1' in affected
    assert affected['S1'].affected == 2
    assert affected['S1'].is_affected
    assert 'S2' not in affected


@patch('talos.pedigree_parser.to_anypath')
def test_as_singletons(mock_to_anypath):
    mock_to_anypath.return_value.open = mock_open(read_data=PED_CONTENT)
    parser = PedigreeParser('dummy_path')
    singletons = parser.as_singletons()
    for p in singletons.values():
        assert p.father_id == '0'
        assert p.mother_id == '0'


@patch('talos.pedigree_parser.to_anypath')
def test_validate_sex_valid(mock_to_anypath):
    mock_to_anypath.return_value.open = mock_open(read_data=PED_CONTENT)
    parser = PedigreeParser('dummy_path')
    assert parser.validate_sex('1', 'S1') == 1
    assert parser.validate_sex('2', 'S2') == 2
    assert parser.validate_sex('0', 'S3') == 0
    assert parser.validate_sex('male', 'S4') == 1
    assert parser.validate_sex('female', 'S5') == 2
    assert parser.validate_sex('unknown', 'S6') == 0
    assert len(parser.parsing_issues) == 0


@patch('talos.pedigree_parser.to_anypath')
def test_validate_sex_invalid(mock_to_anypath):
    mock_to_anypath.return_value.open = mock_open(read_data=PED_CONTENT)
    parser = PedigreeParser('dummy_path')
    assert parser.validate_sex('not_a_sex', 'S7') == 0
    assert 'Invalid Sex provided! Sample S7: not_a_sex' in parser.parsing_issues


@patch('talos.pedigree_parser.to_anypath')
def test_validate_affected_valid(mock_to_anypath):
    mock_to_anypath.return_value.open = mock_open(read_data=PED_CONTENT)
    parser = PedigreeParser('dummy_path')
    assert parser.validate_affected('0', 'S1') == 0
    assert parser.validate_affected('1', 'S2') == 1
    assert parser.validate_affected('2', 'S3') == 2
    assert parser.validate_affected('unaffected', 'S4') == 1
    assert len(parser.parsing_issues) == 0


@patch('talos.pedigree_parser.to_anypath')
def test_validate_affected_grudgingly_valid(mock_to_anypath, caplog):
    mock_to_anypath.return_value.open = mock_open(read_data=PED_CONTENT)
    parser = PedigreeParser('dummy_path')
    assert parser.validate_affected('affected', 'S5') == 2
    assert parser.validate_affected('true', 'S6') == 2
    assert len(parser.parsing_issues) == 0
    assert 'Grudgingly valid Affected status provided, please correct data! Sample S5: affected' in caplog.text
    assert 'Grudgingly valid Affected status provided, please correct data! Sample S6: true' in caplog.text


@patch('talos.pedigree_parser.to_anypath')
def test_validate_affected_invalid(mock_to_anypath):
    mock_to_anypath.return_value.open = mock_open(read_data=PED_CONTENT)
    parser = PedigreeParser('dummy_path')
    assert parser.validate_affected('not_status', 'S7') == 0
    assert 'Invalid Affected status provided! Sample S7: not_status' in parser.parsing_issues


@patch('talos.pedigree_parser.to_anypath')
def parsing_failure(mock_to_anypath):
    mock_to_anypath.return_value.open = mock_open(read_data='FAM1\tS1\t0\t0\t1\tIDK')
    with pytest.raises(ValueError, match='Errors found during pedigree parsing, see log for details'):
        PedigreeParser('dummy_path')


def test_read_write(tmp_path):
    real_ped = str(tmp_path / 'test.ped')
    with open(real_ped, 'w') as f:
        f.write('FAM1\tS1\t0\t0\t1\t2\tHP:0000118,HP:0001250\nFAM1\tS2\t0\t0\t2\t1')
    parser = PedigreeParser(real_ped)
    output_path = str(tmp_path / 'output.ped')
    parser.write_pedigree(output_path)
    with open(output_path) as f:
        content = f.read().strip()
    expected = 'FAM1\tS1\t0\t0\t1\t2\nFAM1\tS2\t0\t0\t2\t1'
    assert content == expected, f'Expected: {expected}, but got: {content}'
