"""
This is a collection of methods to attempt the parsing of BCFTools CSQ
format protein consequence strings into HGVS p. notation

BCFTools outputs c.nnnnnn DNA positions, and long form single-letter AA protein changes
The DNA entries are not as simple to translate from chromosomal to transcript locations

The HGVS python package was supposed to be able to do all this, but it's rubbish
"""

import re

from loguru import logger

IUPAC_LOOKUP = {
    'A': 'Ala',
    'C': 'Cys',
    'D': 'Asp',
    'E': 'Glu',
    'F': 'Phe',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'K': 'Lys',
    'L': 'Leu',
    'M': 'Met',
    'N': 'Asn',
    'P': 'Pro',
    'Q': 'Gln',
    'R': 'Arg',
    'S': 'Ser',
    'T': 'Thr',
    'V': 'Val',
    'W': 'Trp',
    'Y': 'Tyr',
    '*': 'Ter',
}

FRAMESHIFT_RE_1 = re.compile(r'^(?P<codon1>\d+)(?P<ref>[A-Z.*]+)>(?P<codon2>\d+)(?P<alt>[A-Z.*]+)$')
FRAMESHIFT_RE_2 = re.compile(r'^(?P<codon1>\d+)(?P<ref>[\w.*]+)>(?P<codon2>\d+)(?P<alt>[\w.*]+)$')
MISSENSE_RE = re.compile(r'^(?P<codon1>\d+)(?P<ref>\w)>(?P<codon2>\d+)(?P<alt>\w)$')
STOP_GAIN_RE = re.compile(r'^(?P<codon1>\d+)(?P<ref>\w)>(?P<codon2>\d+)\*$')
INT_END_RE = re.compile(r'..(?P<term>\d+)$')
TYPES_RE = re.compile('frameshift|missense|stop_gained|synonymous')


def process_missense(alteration: str) -> str:
    """We expect a missense, so look for a missense."""
    missense_groups = MISSENSE_RE.search(alteration)
    if not missense_groups:
        logger.warning(f'No missense found in {alteration}')
        return alteration

    return f'p.{IUPAC_LOOKUP[missense_groups["ref"]]}{missense_groups["codon1"]}{IUPAC_LOOKUP[missense_groups["alt"]]}'


def process_frameshift(alteration: str) -> str:
    """Tougher one, arbitrary length, may be condensed..."""
    fs_groups = FRAMESHIFT_RE_1.search(alteration)
    if not fs_groups:
        fs_groups = FRAMESHIFT_RE_2.search(alteration)
        if not fs_groups:
            logger.warning(f'No frameshift found in {alteration}')
            return alteration

    codon = int(fs_groups['codon1'])
    ref = fs_groups['ref']
    alt = fs_groups['alt']

    # HGVS nomenclature should start with the first altered base
    codon_bump = False
    if ref[0] == alt[0]:
        codon += 1
        codon_bump = True

    # find a terminating change
    first_ref = IUPAC_LOOKUP[ref[1 if codon_bump else 0]]
    first_alt = IUPAC_LOOKUP[alt[1 if codon_bump else 0]]

    if alt.endswith('*'):
        fs_length = len(alt) - 1 if codon_bump else len(alt)
        return f'p.{first_ref}{codon}{first_alt}fs*{fs_length}'

    # e.g. 834PYVMVRERESFLAPSSGVQP..941>834PVRHGEGKRELPCPQLRCAA..859
    if match := re.search(INT_END_RE, alt):
        end = int(match.group('term'))
        difference = end - codon
        return f'p.{first_ref}{codon}{first_alt}fs*{difference}'

    return alteration


def process_stop_gained(alteration: str) -> str:
    """Example: 128W>128*"""
    sg_groups = STOP_GAIN_RE.search(alteration)
    if not sg_groups:
        logger.warning(f'No stop gained found in {alteration}')
        return alteration

    return f'p.{IUPAC_LOOKUP[sg_groups["ref"]]}{sg_groups["codon1"]}Ter'


def classify_change(alteration: str, consequence: str | None = None) -> str:
    """Detect some known variant types."""
    # trust correctly assigned consequences
    if alteration.startswith('p.'):
        logger.debug(f'Detected {alteration} as a p., no processing')
        return alteration

    if consequence:
        match consequence:
            case 'missense':
                return process_missense(alteration)
            case 'frameshift':
                return process_frameshift(alteration)
            case 'stop_gained':
                return process_stop_gained(alteration)
            case 'synonymous':
                return 'p.='

    return alteration