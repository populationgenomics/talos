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

- Comp-Het mapping: a lookup object of `Var_1: Var_2` for all variant pairs
  - When each variant is assessed, we check if the variant fits with the expected MOI alone. If
  the variant alone doesn't pass the relevant MOI test, we can then check if it was seen as a compound-het with a
  `second hit`. If this is the case then we can report the variant pair as a plausible variant combination.
- Configuration File: used throughout the pipeline, this contains all runtime settings
  - When we are assessing each variant, we may choose to apply thresholds (e.f. MAF, frequency in a joint call). These
  may change with each run, so we take those parameters from a central configuration file
- PanelApp data: associates genes with evidenced inheritance patterns
  - For each 'green' (high evidence) gene contains inheritance pattern to be applied & if the gene as 'new'
(since a given date, or specific gene list).

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

The AbstractVariant has a few internal methods

 - `is_classified`: returns True if any of the category booleans are set to True
 - `category_4_only`: returns True if only the category_4 boolean is True
 - `category_ints`: returns an ordered list of integers for each True boolean
   - i.e. if a variant is True for category 2 and 4, the return will be `[2, 4]`

Based on rough calculations, each AbstractVariant consumes approximately 2.5KB, forming a smaller memory footprint than
the same variant call in the VCF. There are a number of design factors that limit the size of the variants:

 - any variants not on Green genes are discarded
 - any transcripts not MANE or protein_coding are discarded
 - any remaining consequences which are unlikely to be relevant are discarded

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

At this stage, compound-het checking simply takes a variant read from the VCF, and for the given sample, gene, and
position, checks if a paired variant was found. If so, the 2nd var coordinates are logged as 'supporting' the 'primary'
variant.

Post MVP it will make sense to load more logic here, e.g.

 - for a compound het in a proband, we must also check that the same pair doesn't occur in an unaffected parent (
currently possible, but not implemented)
 - for 2 compound het variants, it is important to check that they were not in phase (where possible to determine during
variant calling)

These checks will require loading the full variant-representation pair, so that genotypes can be checked in more detail.
To do this, collecting all the variants together and iterating over in gene-groups suitable provides granularity.

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
