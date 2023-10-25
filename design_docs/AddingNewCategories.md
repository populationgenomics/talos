# New Categories

This framework is designed to make the addition of new categories super simple. The minimal changes
required to create a new category are:

1. Add new Category name/number and description to the config file (
   e.g. [here](https://github.com/populationgenomics/automated-interpretation-pipeline/blob/afcf1bfa2acc30803558fa2092fab4fd8b0a58a5/reanalysis/reanalysis_global.toml#L54))
2. If new fields are acted upon (e.g. a new annotation field), add those to the `CSQ` field in config to ensure the
   values are exported in the labelled VCF (
   e.g. [here](https://github.com/populationgenomics/automated-interpretation-pipeline/blob/afcf1bfa2acc30803558fa2092fab4fd8b0a58a5/reanalysis/reanalysis_global.toml#L36))
   Without this change, the values will not be pulled from the MT, and cannot be presented in the report
3. Add a new category method in the [hail_filter_and_label.py script](../reanalysis/hail_filter_and_label.py), (
   e.g. [here](https://github.com/populationgenomics/automated-interpretation-pipeline/blob/afcf1bfa2acc30803558fa2092fab4fd8b0a58a5/reanalysis/hail_filter_and_label.py#L622-L658).
   This method should stand independently, and contain all the logic to decide whether the label is applied or not.
   This encapsulation should also include the decision about whether a classification is Boolean (once per variant),
   Sample (only relevant to a subset of Samples), or Support (a lesser level of significance)
4. Add a new diagram describing the decision tree to the [images folder](images), and reference it in the
   [README](Hail_Filter_and_Label.md)
5. If you require new fields to be displayed in the HTML report, make the appropriate changes to the templates
6. If you need additional logic (e.g. when this category is assigned we should interpret the variant under a partial
   penetrance model), that's... more complicated. Get in touch with the team to discuss.
