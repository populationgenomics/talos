"""
test seqr interactions
"""


from reanalysis.seqr_interactions import filter_aip_to_new_flags


def test_filter_aip_to_new_flags_no_family_match():
    """
    check that the flag reduction works
    """

    aip_flags = {'family1': {'variant1': {'tag1', 'tag2'}}}
    seqr_flags = {'family2': {}}

    assert [
        x in filter_aip_to_new_flags(aip_flags, seqr_flags)
        for x in [
            {'family': 'family1', 'tag': 'tag1', 'variant': 'variant1'},
            {'family': 'family1', 'tag': 'tag2', 'variant': 'variant1'},
        ]
    ]


def test_filter_aip_to_new_flags_no_variant_match():
    """
    check that the flag reduction works
    """

    aip_flags = {'family1': {'variant1': {'tag1', 'tag2'}}}
    seqr_flags = {'family1': {'variant2': {}}}

    assert [
        x in filter_aip_to_new_flags(aip_flags, seqr_flags)
        for x in [
            {'family': 'family1', 'tag': 'tag1', 'variant': 'variant1'},
            {'family': 'family1', 'tag': 'tag2', 'variant': 'variant1'},
        ]
    ]


def test_filter_aip_to_new_flags_one_matches():
    """
    check that the flag reduction works
    """

    aip_flags = {'family1': {'variant1': {'tag1', 'tag2'}}}
    seqr_flags = {'family1': {'variant1': {'tag2'}}}

    assert filter_aip_to_new_flags(aip_flags, seqr_flags) == [
        {'family': 'family1', 'tag': 'tag1', 'variant': 'variant1'}
    ]


def test_filter_aip_to_new_flags_nothing_needed():
    """
    check that the flag reduction works
    """

    aip_flags = {'family1': {'variant1': {'tag1', 'tag2'}}}
    seqr_flags = {'family1': {'variant1': {'tag1', 'tag2'}}}

    assert not filter_aip_to_new_flags(aip_flags, seqr_flags)


def test_filter_aip_to_new_flags_1_var_each_fam():
    """
    check that the flag reduction works
    """

    aip_flags = {
        'family1': {'variant1': {'tag1', 'tag2'}},
        'family2': {'variant2': {'tag2'}},
    }
    seqr_flags = {
        'family1': {'variant1': {'tag2'}},
        'family2': {'variant2': {'tag1'}},
    }

    assert [
        x in filter_aip_to_new_flags(aip_flags, seqr_flags)
        for x in [
            {'family': 'family1', 'tag': 'tag1', 'variant': 'variant1'},
            {'family': 'family2', 'tag': 'tag2', 'variant': 'variant2'},
        ]
    ]
