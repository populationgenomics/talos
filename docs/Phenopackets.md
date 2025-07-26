# Phenopackets

In Talos we have chosen to use a standardised format for encoding phenotype terms for cohort members. The chosen file format is a [GA4GH Cohort](https://phenopacket-schema.readthedocs.io/en/latest/cohort.html), containing a [GA4GH Phenopacket](https://phenopacket-schema.readthedocs.io/en/latest/phenopacket.html#rstphenopacket) for each participant in the analysis. An example is provided [here](../nextflow/inputs/test_cohort_phenopackets.json).

This README documents the expected content if a Cohort/Phenopacket file is provided directly, and the script available to assist generating this file format from plain Pedigrees.

---

In this data structure:

- A single `Cohort` encapsulates all individuals in the analysis
- Each individual in the analysis is a `Phenopacket` in `Cohort.members`
- The `Member.id` for each individual must match the Sample ID in the VCF/variant data
- Each `Member` has an [Individual](https://phenopacket-schema.readthedocs.io/en/latest/individual.html#rstindividual) entity as `Member.subject`
- `Member.subject.id` is the Family ID, matching the Pedigree (The `Cohort` structure doesn't accommodate Family IDs, and a [Family object](https://phenopacket-schema.readthedocs.io/en/latest/family.html) doesn't accommodate multiple separate families)
- `Member.subject.alternate_ids` is a list containing a single String - if participants are known by two identifiers (e.g. anonymised internal ID and project external ID) this should be the external ID. If the sample/participant only has one ID this should be a repetition of `Member.id`.
- `Member.phenotypic_features` is a list of [PhenotypicFeature](https://phenopacket-schema.readthedocs.io/en/latest/phenotype.html#phenotypicfeature) objects. The minimum content for each element is `{"type": {"id": String}}`, where the String is an HPO term.

---

If a process for generating this file format is not available, the Talos repository contains a script to convert a customised [Pedigree file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format) into a Phenopacket JSON file. Ths script is [ConvertPedToPhenopackets.py](../src/talos/ConvertPedToPhenopackets.py), and a guide to its usage is here:

This script requires two CLI flags:

* `--input`: path to a `.ped` format file
* `--output`: root output to write both a standard Pedigree and a Phenopacket file
  * Phenopackets written to `{OUTPUT}_phenopackets.json`
  * Pedigree written to `{OUTPUT}_pedigree.ped`

The input Pedigree can be a standard 6-column pedigree file, in which case a Phenopacket JSON file is created which contains no phenotypic information.

To create a full-featured Phenopacket file, a Pedigree with extra columns can be provided. The standard columns:

1. Family ID
2. Individual ID
3. Father ID
4. Mother ID
5. Sex
6. Phenotype/Affected

Additional columns, where columns [8, 9, ...] can be specified as many times as required to include a wide range of HPO terms:

7. External ID
8. HPO terms

If additional columns are present, the first (External ID) is mandatory. This can be a repetition of the Individual ID if the participant is only known by a single ID, or a different ID (e.g. to connect the pseudonymised sample ID in the callset to a sample/individual ID). This will be stored in the result metadata as `sample.ext_id`.

Columns 8+ can be HPO terms, with an arbitrary number of terms per line. These HPO terms will be added to the Cohort file, and listed within the Phenopacket for the individual participant. There is no obligation to include HPO terms for participants, even if they are probands within their respective families.

This snippet demonstrates an extended pedigree which would be parsed as:

* 4 individuals across 2 families
* sample ID `proband`, family ID `1`, external ID `proband1` has 3 HPO terms
* `mother` and `father` in family `1` have no phenotypic terms, and are unaffected
* sample ID `proband_2`, family ID `2`, external ID `proband_2` has 1 HPO term

```text
1	proband	father	mother	1	2	proband1	HP:0002779	HP:0004322	HP:0009145
1	father	0	0	1	1	father1
1	mother	0	0	2	1	mother1
2	proband_2	0	0	2	2	proband_2    HP:0009873
```
