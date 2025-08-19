# Adding New Categories

1. `Boolean`
    * Naming convention: categoryboolean[NAME]
    * Content: the value 0, or 1
    * The category is a binary flag, either the variant has the flag assigned or does not. These flags are based on the *variant* annotations, so the flag will apply equally to all samples with the variant.
2. `Samples`
    * Naming convention: categorysample[NAME]
    * Content: a "comma,delimited" list of sample IDs, or "missing" if none
    * This type indicates that the flag has been assigned to only the identified samples, rather than all samples with the variant call. An example of this is _de novo_, where the assignment of the flag is conditional on the MOI, so this won't apply to all samples with a variant call. When processing these variants, only variant calls for samples in this list are treated as being categorised.
3. `Details`
    * Naming convention: categorydetails[NAME]
    * Content: bespoke
    * Any flag starting with _categorydetails_ is processed in some way upon ingestion of the VCF. An example is the `pm5` category, where the flag is not a boolean, or per-sample, but includes compound data to be digested when the VCF is read [here](https://github.com/populationgenomics/talos/blob/bc9385a21fcd9a0292d37b9ecf0f2faccae24f9e/src/talos/utils.py#L383) and [here](https://github.com/populationgenomics/talos/blob/bc9385a21fcd9a0292d37b9ecf0f2faccae24f9e/src/talos/utils.py#L309). The intention is that once parsed, it is converted into a simple Boolean or Sample label, with any other relevant data stored in the info dict.

This framework is designed to make the addition of new categories super simple. The minimal changes required to create a new category are:

1. Add new Category name/number and description to the config file (
   e.g. [here](https://github.com/populationgenomics/talos/blob/main/src/talos/example_config.toml#L61))
2. Add a new category method in the [RunHailFiltering.py script](../reanalysis/RunHailFiltering.py), (
   e.g. [here](https://github.com/populationgenomics/talos/blob/main/src/talos/RunHailFiltering.py#L410-L443)).
   This method should stand independently, and contain all the logic to decide whether the label is applied or not. This encapsulation should also include the decision about whether a classification is Boolean (True/False once per variant, annotate with `0/1`), Sample (only relevant to a subset of Samples, annotate with a comma-delimited list of Sample IDs), or Support (a lesser level of significance). Name your category accordingly.
3. Annotate the data using your new method, i.e. a new call [here](https://github.com/populationgenomics/talos/blob/main/src/talos/RunHailFiltering.py#L919)
4. Add your new category to the filtering method [here](https://github.com/populationgenomics/talos/blob/main/src/talos/RunHailFiltering.py#L671)
5. If required (details category), add some new parsing logic to the `create_small_variant` ingestion method [here](https://github.com/populationgenomics/talos/blob/bc9385a21fcd9a0292d37b9ecf0f2faccae24f9e/src/talos/utils.py#L362)
