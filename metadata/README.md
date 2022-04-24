This directory stores metadata used in analysis such as info about samples, color palettes, gene lists.

## Contents

* `souporcell_results/`: barcodes from potential maternal contamination 
* `old_annotations/`: collection of uniformed cell type annotations from previous publications and original authors for scRNA-seq datasets included in the study
* `marker_genes/`: tables of marker genes for each annotated cell type
* `TCR_genes/`: tables of TCR genes by chain and segment types, used for segment usage analysis
* `220621_FACs_gating_proportions_correct.csv`: table of CD45 FACS gate quantifications (used for adjusting sorting effects in differential abunance analysis) 
* `abTCR_metadata_23032022.csv`: table matching abTCR library names to gene expression library names
* `gdTCR_metadata_23032022.csv`: table matching gdTCR library names to gene expression library names
* `BCR_metadata_29032022.csv`: table matching BCR library names to gene expression library names
* `ATO_metadata_26102021.csv`: table of metadata on ATO scRNA-seq libraries
* `cc_genes.csv`: list of cell cycle genes used for feature selection
* `organ_colors.csv`: color palette for organs
* `anno_groups.json`: dictionary of assignment of cell type labels to broad lineage groups
* `broad_annotation_dict.json`: dictionary of assignment of cell type labels to broad lineage groups (as used in Fig1)
* `broad_annotation_colors.csv`: color palette for broad lineage groups (as used in Fig1)
* `HGNC_chemokine_receptors.csv`: list of chemokines and receptor genes from https://www.genenames.org/
* `panfetal_visium_metadata_v2.csv`: metadata for visium slides
* `sorted_B_sample_manifest.csv`: table of metadata on sorted B1 cells scRNA-seq libraries
* `progenitor_groups.json`: dictionary of assignment of cell type labels to lineage progenitor groups (used for widespread hematopoiesis analysis)


