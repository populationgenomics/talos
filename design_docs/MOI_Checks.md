# MOI Checks

Following the labelling and VCF re-headering, we have a heavily filtered VCF file containing all the required
labels and annotations, accompanied by a JSON file describing all the Compound-Het variant pairs. Combining this
information with the pedigree, we can iterate through each of the variants, and check if the evidence supports the
PanelApp suggested Mode Of Inheritance for the corresponding gene.

For each participant, we compile a list of variants where the inheritance model matches the PanelApp MOI, e.g.

- Variants in Monoallelic genes must pass additional population frequency filters
- Biallelic variants must be Homozygous or supported by a 'second-hit'
- X-Hemizygous variants must be monoallelic in males, and biallelic in females

## Extensions post MVP

1. Enable partially penetrant disease/affection-status
2. Familial inheritance checks using a supplied pedigree file

---

## Logic

This process uses two types of object - static reference objects, and the variant objects parsed directly from the VCF.
Variant objects are dynamic; each one read from the source contains all the relevant annotations to help interpret it.
Static references (e.g. PanelApp content) are used as an accessory, allowing us to interpret the annotations on the
variant to decide whether we include it in the final report.

Static Objects:

- Comp-Het mapping: a lookup object of the strings `Var_1: [Var_2, Var_3, ...]` for all variant pairs
  - As each variant is assessed, we first check if the variant fits with the expected MOI alone. If
  the variant alone doesn't pass the relevant MOI test, we can then check if it was seen as a compound-het with a
  `second hit`. If this is the case then we can report the variant pair as a plausible variant combination.
- Configuration File: used throughout the pipeline, this contains all runtime settings
  - When we are assessing each variant, we may choose to apply thresholds (e.f. MAF, frequency in a joint call). These
  may change with each run, so we take those parameters from a central configuration file
- PanelApp data: associates genes with evidenced inheritance patterns
  - For each 'green' (high evidence) gene contains inheritance pattern to be applied & if the gene as 'new'
  (since a given date, or specific gene list).
- Gene Lookup: a dictionary of all variants in this gene, indexed on `chr-pos-ref-alt` representation
  - When we check the partner variants of a compound-het, we also need to check whether or not other family members
  also have this same compound het pair. To do that, we must be able to reach the representation of the variant. For
  this purpose we pass a collection of variants to the `MOI.run()` method, so that we can access the MOI for all vars
  being considered. We could also do this same test using the comp-het dictionary, but there are other reasons to group
  variants by gene (parsing variants as a gene group is a logical level for parallelisation, as all X-variant impacts
  will be covered by grouping at this level) so this logic feels more versatile. Note: parsing variants as a gene-group
  removes need for separately gathering a compound-het dictionary

Dynamic Objects: AbstractVariant

Each variant in turn is read from the source (VCF) and cast in a custom AbstractVariant Class representation. This is
for a couple of reasons:

1. CyVCF2 was chosen as a parsing library due to speed relative to pyVCF, but each library provides a slightly different
representation. This normalises all logic to a common class, and allows for different sources/parsers to generate a
common format.
2. CyVCF2 and PyVCF variant object formats contain structures which cannot be pickled by default. This leads to issues
with introducing parallelisation into the code. A dictionary/set/list-based class can be pickled easily, so async/await
structure can be added in at any level.
3. Making an abstract object from simple types will lead to simpler unit testing.

We open and read through the VCF, and for each variant create an AbstractVariant representation. This holds the
following details:

1. A Coordinate object (CHR, POS, REF, ALT)
2. The INFO content as a flat `key:value` Dictionary
3. Set of all samples Heterozygous at this locus
4. Set of all samples Homozygous at this locus
5. A Boolean for each category-label used in reanalysis
   1. `Category_4` (_de novo_) is slightly different; this value is a list of Strings, identifying all samples which were
   confirmed to be _de novo_ for this variant call. This doesn't preclude other non-reference samples at this site, but
   sample ID presence in this list means all criteria of the _de novo_ test in Hail were satisfied.

The AbstractVariant has a few internal methods

- `is_classified`: returns True if any category was assigned to this variant
- `category_1_2_3`: returns True if any of category 1, 2, or 3 were assigned
- `support_only`: returns True if only the supporting category was assigned
- `category_ints`: returns a list of strings for each assigned category, i.e.
  - if a variant is True for category 2 and 3, the return will be `[2, 3]`
  - if a variant is True for category 2 and 4, the return will be `[2, de_novo]`
- `sample_specific_category_check`: returns True if the specific sample is _de novo_,
or if the variant has category 1, 2, or 3 assigned

## Variant Gathering

The Variant gathering strategy observed in this application is:

 - parse the VCF header to obtain all chromosome names
 - for each contig in turn, extract all variants & create an Abstract representation of each
 - group variants by gene ID
   - each variant contains exactly one gene annotation, from an earlier MT split
   - each gene can be processed separately, allowing a comp-het workaround
 - once all variants on a contig are extracted, iterate over variants grouped by gene

Justification:

1. Grouping variants by gene is a fallback for the compound-het process, which is currently causing problems, allowing
all possibly comp-het variants to be loaded together and evaluated
2. Use of the AbstractVariant structure means each variant can be pickled, so each group can be processed in parallel
(Cyvcf2 and pyvcf objects can't be pickled directly)
3. Allows us to collect the panelapp details once per gene, instead of once per variant; minor efficiency gain
4. Gathering all variants relevant to an MOI test means that when generating the results, we can access all attributes
of the relevant variants for presentation, e.g. instead of just variant coordinates for the variant pair, we can show
the exact consequence(s) that led to category labelling

## Compound-Heterozygous checks

Currently the compound-het checks are done in Hail, resulting in a simple result format:

```json
{
    "sample": {
        "gene": {
            "chr-pos-ref-alt": ["chr-pos2-ref-alt", "chr-pos3-ref-alt", "..."]
        },
        "gene2": {
            "chr-pos-ref-alt": ["chr-pos2-ref-alt", "chr-pos3-ref-alt", "..."]
        }
    }
}
```

This is a highly condensed representation, and doesn't hold any annotations from the relevant variants. It simply
contains information to state that the named sample(s) have co-located variants within the same gene.

Compound-het checks are triggered once we find  consists of parsing the above format to
At this stage, compound-het checking simply takes a variant read from the VCF, and for the given sample, gene, and
position, checks if a paired variant was found. If so, the 2nd var coordinates are logged as 'supporting' the 'primary'
variant.

When evaluating a biallelic MOI, we are able to consider the presence of a compound-het within the family group. See the
Family Checks section below for implementation detail.

Post MVP it will make sense to check that compounded variants are not in phase (where possible to determine during
variant calling).

## MOI Tests

Flexible implementation for calculating MOI.

Includes one controlling class `MoiRunner`, and one functional class `BaseMoi` and its children. The `BaseMoi`
derivatives each define a single Mode of Inheritance e.g.

- DominantAutosomal - requires a variant to be exceptionally rare, and homozygotes to be absent from population
databases. All samples with a heterozygous variant call are passed as fitting with the MOI

- XRecessive - Male samples are passed so long as the AF is below a stringent threshold, Female samples must be
Homozygous or supported in a compound-het

The separation between the methods defining the filter algorithm, and the MoiRunner using one or more algorithms
allows multiple filters to be applied for a given MOI string e.g.

- "BOTH monoallelic and biallelic, autosomal or pseudoautosomal" from PanelApp is easy to interpret as
2 filters, monoallelic, and biallelic
- "X-LINKED: hemizygous mutation in males, biallelic mutations in females" for male samples we call a
DominantMOI filter, for females we call a Recessive filter

The usage paradigm is:

1. Create an instance of the `MoiRunner`, passing a target MOI string and some configuration parameters
2. During the setup, the MOI string is used to determine which filters to add into the filter list
3. The `MoiRunner` implements a `.run()` method, which takes a single variant and passes through each filter in turn
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
