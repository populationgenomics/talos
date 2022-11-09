"""
tests for clinvar manual summaries
"""


from copy import deepcopy
from datetime import datetime

from clinvar.conflict_huntr import (
    ACMG_THRESHOLD,
    Consequence,
    Submission,
    # county_county,
    check_stars,
    acmg_filter_submissions,
)


CURRENT_TIME = datetime.now()
BASIC_SUB = Submission(CURRENT_TIME, 'submitter', Consequence.UNKNOWN, 'review')


def test_check_stars_practice():
    """

    Returns:

    """
    sub1 = deepcopy(BASIC_SUB)
    sub1.review_status = 'practice guideline'
    subs = [BASIC_SUB, sub1]
    assert check_stars(subs) == 4


def test_check_stars_expert():
    """

    Returns:

    """
    sub1 = deepcopy(BASIC_SUB)
    sub1.review_status = 'reviewed by expert panel'
    subs = [BASIC_SUB, sub1]
    assert check_stars(subs) == 3


def test_check_stars_criteria_prov():
    """

    Returns:

    """
    sub1 = deepcopy(BASIC_SUB)
    sub1.review_status = 'something, something, criteria provided'
    subs = [BASIC_SUB, sub1]
    assert check_stars(subs) == 1


def test_check_stars_none():
    """

    Returns:

    """
    sub1 = deepcopy(BASIC_SUB)
    sub1.review_status = 'smithsonian'
    subs = [BASIC_SUB, sub1, BASIC_SUB]
    assert check_stars(subs) == 0


def tests_acmg_filter_neutral():
    """
    filter submissions against ACMG date threshold
    """
    subs = [BASIC_SUB, BASIC_SUB]
    assert acmg_filter_submissions(subs) == subs


def tests_acmg_filter_removes():
    """
    filter submissions against ACMG date threshold
    """
    sub1 = deepcopy(BASIC_SUB)
    sub1.date = datetime(year=1970, month=1, day=1)
    sub2 = deepcopy(BASIC_SUB)
    sub2.date = datetime(year=2000, month=1, day=1)
    subs = [BASIC_SUB, sub1, sub2]
    assert acmg_filter_submissions(subs) == [BASIC_SUB]


def tests_acmg_filter_gte():
    """
    filter submissions against ACMG date threshold
    """
    sub1 = deepcopy(BASIC_SUB)
    sub1.date = ACMG_THRESHOLD
    subs = [BASIC_SUB, sub1]
    assert acmg_filter_submissions(subs) == subs
