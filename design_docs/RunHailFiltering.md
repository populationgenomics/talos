# RunHailFiltering.py

This stage of the workflow uses hail query to load, filter, and apply provisional categories to variants.

---

Takes an annotated VCF in Hail MatrixTable format, adds category classifications based on the variant annotations. These
diagrams describe the logic used to select variants into each category. Assignment of each category is independent, and
variants can have multiple categories assigned.

### Input

Annotated Joint-called MatrixTable. Annotations applied either by, or consistent with, the prior annotation step.

---

### Types of Category flags

To make downstream operations easier, the different categories are grouped into types. There are currently 3 types:

1. Boolean - The category is a binary flag, either the variant has the flag assigned or does not. These flags are based
   on the *variant* annotations, so the flag will apply equally to all samples with an alt call at that position
2. Samples - The category is a list of sample IDs or 'missing'. This type indicates that the flag has been assigned to
   only the identified samples, rather than all samples with the variant call. An example of this is _de novo_, where
   the assignment of the flag is conditional on the MOI, so this won't apply to all samples with a variant call. When
   processing these variants, only variant calls for samples in this list are treated as being categorised
3. Support - Any flag starting with _CategorySupport_ is treated equally, but inferior to all other Categories. This
   means that the Support flags are never enough to categorise a variant alone, but may support a separate categorised
   variant in a compound inheritance MOI. If a variant has a Support flag & non-support flags, it will be treated as an
   independent variant
4. Details - Any flag starting with _CategoryDetails_ is processed in some way upon ingestion of the VCF. The content of
   the category label can be anything, with the intention that when each variant is read it is converted into a Boolean
   or Sample label. The only current implementation of this is the PM5 label. The flag content in this case is a list of
   all clinvar Pathogenic missense alleles which affect the same residue as this current variant. Upon ingestion of the
   variant, the collected AlleleIDs are split, filtered to unique, and any exact matches to the current variant are
   removed. If there are remaining AlleleIDs in the list, an entry is made in the variant's info dictionary (for
   rendering
   downstream in the report) and a flag CategoryBooleanPM5 is assigned. This then acts as a regular boolean flag, with
   the pm5 data in the dictionary available for display if appropriate. This approach means that upon creation and init
   processing of the AbstractVariant, the CategoryDetails flag no longer exists, and no advanced logic is needed to
   process it within the MOI logic.
