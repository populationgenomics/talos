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

Talos uses a companion library [MendelBrot](https://github.com/MattWellie/mendelbrot) to parse Pedigrees. Just what the world needed, yet another pedigree parsing utility... Here the decision was driven by needing specific validation (e.g. allowing for 6 or 7 column data, and validating the HPO terms if provided) and specific functionality (being able to re-write the pedigree as a singletons, or as a strict 6-column for Hail compatibility).

The parser can be used to read a pedigree file, and will return a `Pedigree` object containing parsed data from the file. The parser will also validate the file format, and raise an error if any issues are found.
