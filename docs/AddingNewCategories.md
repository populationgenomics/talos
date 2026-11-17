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
    * Any flag starting with _categorydetails_ is processed in some way upon ingestion of the VCF. An example is the `pm5` category, where the flag is not a boolean, or per-sample, but includes compound data to be digested when the VCF is read. The intention is that once parsed, it is converted into a simple Boolean or Sample label, with any other relevant data stored in the info dict.

This framework is designed to make the addition of new categories simple. The minimal changes required to create a new category are:

1. Add new Category name/number and preferred String representation in the models.py file
2. Add a new classification decision in the write_gene_rows method, in the `run_small_filtering.py` script. This must also include a decision about whether a classification is Boolean (True/False once per variant, annotate with `0/1`), Sample (only relevant to a subset of Samples, annotate with a comma-delimited list of Sample IDs), or Details - Name your category accordingly.
3. If required (details category), add some new parsing logic to the `create_small_variant` ingestion method in utils.py
