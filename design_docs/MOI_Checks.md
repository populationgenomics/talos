# MOI Checks

Following the category labelling stage, we have a heavily filtered VCF. Combining this information with the pedigree, we
iterate through each variant and check if the gene's PanelApp MOI fits for any samples. For each participant, we compile
a list of variants where the inheritance model matches the PanelApp MOI, e.g.

- Variants in Monoallelic genes must pass additional population frequency filters
- Biallelic variants must be Homozygous or supported by a 'second-hit'
- X-Hemizygous variants must be monoallelic in males, and biallelic in females

---

## Logic

This process uses two types of object - static reference objects, and the variant objects parsed directly from the VCF.
Variant objects are dynamic; each one read from the source contains all the relevant annotations to help interpret it.
Static references (e.g. PanelApp content) are used as an accessory, allowing us to interpret the annotations on the
variant to decide whether we include it in the final report.

Static Objects:

- Configuration File: used throughout the pipeline, this contains all runtime settings
  - When we are assessing each variant, we may choose to apply thresholds (e.f. MAF, frequency in a joint call). These
  may change with each run, so we take those parameters from a central configuration file
- PanelApp data: associates genes with evidenced inheritance patterns
  - For each 'green' (high evidence) gene contains inheritance pattern to be applied & panels where the gene is 'new'
  (if any).

Dynamic Objects:

- AbstractVariant: Each variant in turn is read from the source (VCF) and cast in a custom AbstractVariant Class
representation. This is for a couple of reasons:

  1. CyVCF2 was chosen as a parsing library due to speed relative to pyVCF, but each library provides a slightly different
  representation. This normalises all logic to a common class, and allows for different sources/parsers to generate a
  common format.
  2. CyVCF2 and PyVCF variant object formats contain structures which cannot be pickled by default. This leads to issues
  with introducing parallelisation into the code. A dictionary/set/list-based class can be pickled easily, so async/await
  structure can be added in at any level.
  3. Making an abstract object from simple types will lead to simpler unit testing.

- Compound Het lookup: This object is built per-gene as we iterate over all genes independently

We open and read through the VCF, and for each variant create an AbstractVariant representation. This holds the
following details:

1. A Coordinate object (CHR, POS, REF, ALT)
2. The INFO content as a flat `key:value` Dictionary
3. Set of all samples Heterozygous at this locus
4. Set of all samples Homozygous at this locus
5. Lists of all present 'categories', separated by type (boolean, sample, support) - these are checked using methods
   rather than storing static booleans
6. If physical phasing is evident for any sample calls, record the numerical phase set ID, and sample-specific GT
    1. Within a phase set ID, two variants are in-phase if their GT is also the same, i.e. for phase set #1, variants
   with the GT `0|1` and `1|0` are explicitly out of phase, so knowing the PS ID alone is not sufficient

The AbstractVariant has a few internal methods:

- `has_boolean_categories`: True if any binary category is assigned to this variant
- `has_sample_categories`: True if any sample category is assigned to this variant, takes a sample ID as an argument
- `has_support`: returns True if the supporting category is assigned
- `support_only`: returns True if ONLY the supporting category is assigned
- `category_non_support`: returns True if any non-support category is assigned
- `is_classified`: returns True if any category is assigned
- `category_values`: returns a list of strings for each assigned category, i.e.
  - if a variant is True for category 2 and 3, the return will be `['2', '3']`
  - if a variant is True for category 2 and 4, the return will be `['2', 'de_novo']`
- `sample_specific_category_check`: returns True if the specific sample is _de novo_,
or if the variant has a boolean category assigned - accepts a switch to also check for `Support`
- `get_sample_flags`: checks for any variant flags - currently this only implements one check - AB Ratio test

## Variant Gathering

The Variant gathering strategy observed in this application is:

 - parse the VCF header to obtain all chromosome names
 - for each contig in turn, extract all variants & create an Abstract representation of each
 - group variants by gene ID, forming a structure `{'gene_name': [Variant1, Variant2, ...], }`
   - each variant contains exactly one gene annotation, from an earlier split, so grouping is accurate
   - each gene can be processed separately if required, allowing logical parallelisation breaks
 - once all variants on a contig are extracted, iterate over variants grouped by gene

Justification:

1. Use of the AbstractVariant structure means each variant can be pickled, so each group can be processed in parallel
(Cyvcf2 and pyvcf objects can't be pickled directly)
2. Allows us to collect the panelapp details once per gene, instead of once per variant; minor efficiency gain
3. Gathering all variants relevant to an MOI test means that when generating the results, we can access all attributes
of the relevant variants for presentation, e.g. instead of just variant coordinates for the variant pair, we can show
the exact consequence(s) that led to category labelling

## Compound-Heterozygous checks

Compound-Het testing has been migrated from Hail to Python for stability reasons. As we process a list of all variants
associated with a gene, we first run all variant-pair permutations, and check for compound hets. During this check we
use a couple of additional rules:

- If the variant is on Y or MT, don't consider for compound-het
- As we check per-sample genotypes, don't consider Male hets on X
- For each variant pair, check that they aren't in-phase from the read backed data

Assemble a lookup dictionary, indexed primarily on SampleID, with a second layer of indexing which is the String repr.
of the variant coordinates. If we are assessing a variant, we would find out it's string representation, and if that
exists in the dictionary as a key, the terminal value is a list of all the Variant objects where a compound-het exists
for this sample. These objects can be retrieved directly from this dictionary, and pairs can be passed to the familial
inheritance checks.

```python
from reanalysis.utils import AbstractVariant
_comp_het = {
    "sample": {
        "coords_string_1": [AbstractVariant, AbstractVariant],
        "coords_string_2": [AbstractVariant, AbstractVariant]
    }
}
```

When evaluating biallelic MOI, we consider the presence of a compound-het within the family group. See the Family Checks
section below for implementation detail.

## MOI Tests

Flexible implementation for calculating MOI.

Includes one controlling class `MoiRunner`, and one functional class `BaseMoi` and its children. The `BaseMoi`
derivatives each define a single Mode of Inheritance e.g.

- DominantAutosomal - requires a variant to be exceptionally rare, and homozygotes to be absent from population
databases. All samples with a heterozygous variant call are passed as fitting with the MOI

- XRecessive - Male samples are passed so long as the AF is below a stringent threshold, Female samples must be
Homozygous or supported in a trans compound-het

The separation between the methods defining the filter algorithm, and the MoiRunner using one or more algorithms
allows multiple filters to be applied for a given MOI string e.g.

- "BOTH monoallelic and biallelic, autosomal or pseudoautosomal" from PanelApp is easy to interpret as
2 filters, monoallelic, and biallelic
- "X-LINKED: hemizygous mutation in males, biallelic mutations in females" for male samples we call a
DominantMOI filter, for females we call a Recessive filter

The usage paradigm is:

1. Create an instance of the `MoiRunner`, passing a target MOI string and some configuration parameters
2. During the setup, the MOI string is used to determine which filters to add into the filter list
3. The `MoiRunner` implements a `.run()` method, which takes a single variant and the lookup of all compound-hets within
this gene, and passes through each filter in turn
4. Where a variant passes all conditions within a filter class, a 'result' object is created
5. The result of `.run()` is a list of valid modes of inheritance (ReportedVariant), each including details of the
variant(s), and samples

### Family Checks

Instead of tree-traversal within families, MOI checks represent each family as a flat pool
of participants. The inheritance rules are applied unilaterally to every member of a family group, rather than in a
directional manner, e.g. when considering each candidate variant:

- for a monoallelic variant, inherited with complete penetrance, we enforce a unilateral rule - every person with the
variant must also be affected. If this rule is violated for any individual member, the variant does not pass the
inheritance test.
- for a biallelic variant pair, all participants with the variant pair must also be affected. If the variant considered
is Homozygous, we permit unaffected members to be heterozygous. Similarly, for compound-het inheritance, unaffected
participants may have either variant but not both.

Complete and Incomplete penetrance modes are available in this module. For the Complete penetrance
inheritance model:

- all affected individuals must have the same variant
- unaffected members cannot have the variant

The alternative is Incomplete penetrance:

- all affected members must have the candidate variant
- unaffected members may still have the same variant

A situation not permitted under either model is when some, but not all, affected members have a candidate variant.
That would be plausible under a situation of multiple underlying diseases within a single family, but as the application
doesn't currently permit providing family or member-specific disease data, this is assumed not to be the case.

### _De Novo_ checks

Within the hail category-labelling process, we implement a
[ported version of K. Samocha's _de novo_ filtering test](https://github.com/ksamocha/de_novo_scripts) in
[Hail](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.de_novo). This test uses a number of factors in
either confirming or refuting the _de novo_ call, above and beyond just checking genotype calls:

- Checks the Allelic Balance and read depth of the Proband
- Requires a strong Phred-scaled likelihood of a true call
- Assesses the probability that the call is truly _de novo_, instead of a missed call in a parent

Elements of this calculation include the naive probability of this variant being a naturally occurring variant (~1 in
30M), and the likelihood of this variant being a missed call based on the known population frequency. Once this process
runs the different factors are combined into a probability(_de novo_), & binned into HIGH, MEDIUM, and LOW confidences.
For this process, we retain only the variant calls which are HIGH confidence.

Where a _de novo_ variant is identified, the flag annotated is a comma-delimited list of the samples which have a passed
the various _de novo_ tests; where no variants are _de novo_ the annotated value is the String `missing`.

When the VCF is being processed, the `category_4` field is parsed by splitting on commas (if `missing`, this is an empty
list). Where there is at least a single sample in this list, the variant is parsed using a separate `De_Novo` MOI test,
which will generate a reported variant structure for the sample(s) which passed the Hail confirmation process.

*Note*: if multiple samples have variants at a locus, Cat. 4 is only confirmed where the sample IDs are within that Cat.
4 list, showing that they were subject to all tests relevant to _de novo_ status.

## Flags

When a reportable event is found, a JSON blob representing the variant and sample is created. If the relevant gene was
in one or more of the additional panels requested (see [additional panels](PanelApp_interaction.md#per-participant-panels)),
the names/IDs of those panels are appended to the list of any variant-specific panels.

This leaves us the flexibility to mark individual samples/families as having disease-relevant panels, which can then be
cross-referenced against these flags for visual emphasis in the final report.
