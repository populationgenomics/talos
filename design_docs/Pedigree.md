# Pedigree

For this analysis I'm using a lightweight Pedigree representation, extensible to all likely MOI tests.

Each line of a well-formed PED file, representing one participant, is parsed into 'PedEntry', and stored in a dictionary:

```python
from typing import Optional


class PedEntry:
    family: str
    sample_id: str
    father: Optional[str]
    mother: Optional[str]
    affected: bool
    is_female: bool

member_dict = {
    'sample_1': PedEntry(),
    'sample_2': PedEntry(),
}
```

The PedEntry is a placeholder for the per-participant data presented in the raw file. The next step is to build
relationships into the representation. This is done through creation of `Participant` objects. A Participant represents
a node in the family, the PED file details from that node, and the immediate relationships to other members in the PED.

```python
from dataclasses import dataclass
from typing import List, Optional, Set, Type


@dataclass
class Participant:
    """
    dataclass representing a person within a family
    Type['Participant'] has to be used in order to have
    a self-referential class...
    """

    details: PedEntry
    mother: Optional[Type['Participant']]
    father: Optional[Type['Participant']]
    children: List[Type['Participant']]
    affected_parents: Set[str]
    unaffected_parents: Set[str]
```

This structure allows a family to be traversed up/down, starting from any node. e.g. from a proband we can find the
parents easily, and from a parent we can reach all children.

One use case this doesn't allow for is incomplete information in the PED file, i.e. siblings with no parents. For that
we create a final structure, the `families` dictionary. This is a dictionary indexed on family ID, which links to a list
of all the contained members. For most MOI tests we are looking for consistency within a family group, rather than
checking family members with specific directionality. e.g. when considering the impact of `Variant1`, all members of
`Family1` should either have a variant and be affected by it, or be `WT` and not be affected. Such a check doesn't care
about the directionality of child-parent linkage, and so most tests can be implemented as a flat view of a family - the
only requirement is that assumptions hold independently for all members.

The presence of both a flat & linked family structure means that for basic checks, e.g. autosomal dominant, we can run
a flat test on all family members. In future if inheritance tests become more niche, e.g. maternal imprinted disease, we
can use the existing `Pedigree` object to run specific tests per node that takes sample sex and relationship
directionality into account, with no further changes to this representation.
