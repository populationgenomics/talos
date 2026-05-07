# Pedigree

In Talos, we needed a way for users to supply both family structure and per-sample phenotype data. We tried to use a combination of a standard Pedigree and [PhenoPackets](https://www.ga4gh.org/product/phenopackets/), which could contain all the required information. Unfortunately this was not a data format which any of our users were familiar with, so we have reversed that decision, and now combine both sets of information in a single Pedigree format file.

## Format

Broadly this file is based on the [PLINK `.fam` file format](https://www.cog-genomics.org/plink/1.9/formats#fam).

* There are 6 mandatory columns, and 1 optional column.
* The file is a TSV format (each column is separated by a tab character).
* The file should not contain any header lines.
* Each row represents a single participant in the pedigree, and contains their relationship to other participants in the immediate family, sex, and affection status (the disease relevant to the analysis of this family is present/absent)
* All participants should have a separate row in the file, e.g. a mother-father-child trio would have three rows, one for each individual. We tolerate parent IDs only defined on their child's row, but they will be removed, as Talos requires a specific affection status for each participant to make use of them.
* Generally a value of zero (`0`) is used to indicate where data is missing, i.e. parental IDs for samples without parents, or parents with unknown affection status. "-" or an empty string are also accepted, but the parser will issue a warning.
* The parser has a fixed range of possible values for each column, and will fail if invalid values are used, with a descriptive error message.

| Column      | Description                                                                                                                                                                                                                                                                                                               |
|-------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `Family ID` | A unique identifier for each family. This can be any string, but should be consistent across all samples in the family. `0` is not permitted.                                                                                                                                                                             |
| `Sample ID` | A unique identifier for each Sample/participant. This must match the Sample ID used in the VCF/Variant data. `0` is not permitted.                                                                                                                                                                                        |
| `Father ID` | The Sample ID of the father of this participant. If unknown, use `0`.                                                                                                                                                                                                                                                     |
| `Mother ID` | The Sample ID of the mother of this participant. If unknown, use `0`.                                                                                                                                                                                                                                                     |
| `Sex`       | The sex of this row's participant. `1` is Male, `2` is Female, `0` is Unknown. "male" and "female" are accepted as `1` and `2` respectively.                                                                                                                                                                              |
| `Affected`  | The affection status of this participant. `1` is unaffected, `2` is affected, `0` is unknown. "affected" and "unaffected" are accepted as `1` and `2` respectively.                                                                                                                                                       |
| `HPO`       | Optional. A comma- or semicolon-separated list of HPO terms (e.g. `HP:0000118,HP:0001250`). This is used to provide additional phenotypic information for the participant. This column can be completely absent if phenotypic data is not available; the parser will handle column 7 being present for only some samples. |

## Pedigree Module

Talos contains a standalone [pedigree parser](../src/talos/pedigree_parser.py) to satisfy this need. Just what the world needed... Here the decision was driven by needing specific validation (e.g. allowing each row to contain only 6 columns, but validating the HPO terms if provided, providing data in a format Hail can parse) and specific functionality (being able to re-write the pedigree as a singletons).

The parser can be used to read a pedigree file, and will return a `Pedigree` object containing parsed data from the file. The parser will also validate the file format, and raise an error if any issues are found.

This whole module is contained within the single `pedigree_parser.py` file, and has minimal dependencies, so it can be stolen and used in any Python project without needing to install Talos.

Example data:
```tsv
1	SAM1	SAM2	SAM3	1	2	HP:0002779,HP:0004322,HP:0009145
1	SAM2	0	0	1	1
1	SAM3	0	0	2	1
```

```python
from talos.pedigree_parser import PedigreeParser, Participant

pedigree_file = "path/to/pedigree.ped"

pedigree = PedigreeParser(pedigree_file)

# Each pedigree row is parsed into a Participant object, accessible via the `participants` attribute
for sample_id, participant in pedigree.participants.items():
    print(f"""
Sample ID: {sample_id}, Family ID: {participant.family_id}, \
Father ID: {participant.father_id}, Mother ID: {participant.mother_id}, \
Sex: {participant.sex}. Affected: {participant.affected}
    """)

    # You can also access the HPO terms for each participant
    if participant.hpo_terms:
        print(f"HPO Terms: {', '.join(participant.hpo_terms)}")

# You can also access participants grouped by family
for family, participant_list in pedigree.by_family.items():
    print(f"Family ID: {family}")
    for participant in participant_list:
        print(f"  Sample ID: {participant.sample_id}")

# Access to the Sample IDs of all participants in the pedigree
ids = pedigree.get_all_sample_ids()

# Or access to only the affected individuals and their sample IDs
affected_participants = pedigree.get_affected_members()
affected_ids = pedigree.get_affected_member_ids()

# The data can be written back out to a clean 6-column pedigree file
pedigree.write_pedigree("path/to/output.ped")

# there are a few methods to update/edit the parsed pedigree content. These return the participants data, but by default
# do not update the instance's internal data.

# If you want to reduce the pedigree content to only a subset of participants (by sample ID)
subset_data: dict[str, Participant] = pedigree.strip_pedigree_to_samples(['SAM1', 'SAM2', 'SAM3'])

# If you want the parsed pedigree content, but as singletons (remove relationships between samples)
singleton_data = pedigree.as_singletons()

# the write method can take a dictionary of participants and will write that instead of the class instance content
pedigree.write_pedigree("path/to/singletons.ped", participants=singleton_data)

# if you want to permanently update the original pedigree data to singleton/subset data, you can use set_participants
pedigree.set_participants(singleton_data)

# or to only write the subset of the data to the file - these are equivalent:
# a write with an on-the-fly subset
pedigree.write_pedigree("path/to/subset.ped", only_participants=['SAM1', 'SAM2', 'SAM3'])

# overwriting the instance data and then writing
pedigree.set_participants(pedigree.strip_pedigree_to_samples(['SAM1', 'SAM2', 'SAM3']))
pedigree.write_pedigree("path/to/subset.ped")
```
