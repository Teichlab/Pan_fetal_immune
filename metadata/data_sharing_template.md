# Mapping the developing human immune system across organs

## scRNA-seq datasets

Full dataset of single-cell RNA-seq profiles from 9 developmental tissues across gestation (4-17 pcw).

[Interactive viewer](PAN.A01.v01.raw_count.20210429.PFI.embedding.cellxgene.h5ad) | [download .h5ad file](PAN.A01.v01.raw_count.20210429.PFI.embedding.h5ad) | [download scVI model](scRNA_data/scVI_models/scvi_PFI_model) | [download .csv file of cell annotations](scRNA_data/PAN.A01.v01.entire_data_normalised_log.20210429.full_obs.annotated.clean.csv) | [download CellTypist model](scRNA_data/celltypist_model.Pan_Fetal_Human.pkl)

### Cell type lineages

- HSC/immune cells (all hematopoietic-derived cells) | [Interactive viewer](scRNA_data/PAN.A01.v01.raw_count.20210429.HSC_IMMUNE.embedding.cellxgene.h5ad) | [download .h5ad file](scRNA_data/PAN.A01.v01.raw_count.20210429.HSC_IMMUNE.embedding.h5ad) | [download scVI model](scRNA_data/scVI_models/scvi_HSC_IMMUNE_model)
- Stromal cells (all non-immune cells) | [Interactive viewer](scRNA_data/PAN.A01.v01.raw_count.20210429.STROMA.embedding.cellxgene.h5ad) | [download .h5ad file](scRNA_data/PAN.A01.v01.raw_count.20210429.STROMA.embedding.h5ad) | [download scVI model](scRNA_data/scVI_models/scvi_STROMA_model)
- Megakaryocyte/erythroid cells | [Interactive viewer](scRNA_data/PAN.A01.v01.raw_count.20210429.MEM_PROGENITORS.embedding.cellxgene.h5ad) | [download .h5ad file](scRNA_data/PAN.A01.v01.raw_count.20210429.MEM_PROGENITORS.embedding.h5ad) | [download scVI model](scRNA_data/scVI_models/scvi_MEM_PROGENITORS_model)
- HSC/progenitor cells | [Interactive viewer](scRNA_data/PAN.A01.v01.raw_count.20210429.HSC_PROGENITORS.embedding.cellxgene.h5ad) | [download .h5ad file](scRNA_data/PAN.A01.v01.raw_count.20210429.HSC_PROGENITORS.embedding.h5ad) | [download scVI model](scRNA_data/scVI_models/scvi_HSC_PROGENITORS_model)
- Myeloid cells | [Interactive viewer](scRNA_data/PAN.A01.v01.raw_count.20210429.MYELOID_V2.embedding.cellxgene.h5ad) | [download .h5ad file](scRNA_data/PAN.A01.v01.raw_count.20210429.MYELOID_V2.embedding.h5ad) | [download scVI model](scRNA_data/scVI_models/scvi_MYELOID_V2_model)
- Lymphoid cells | [Interactive viewer](scRNA_data/PAN.A01.v01.raw_count.20210429.LYMPHOID.embedding.cellxgene.h5ad) | [download .h5ad file](scRNA_data/PAN.A01.v01.raw_count.20210429.LYMPHOID.embedding.h5ad) | [download scVI model](scRNA_data/scVI_models/scvi_LYMPHOID_model)
- NK/T cells | [Interactive viewer](scRNA_data/PAN.A01.v01.raw_count.20210429.NKT.embedding.cellxgene.h5ad) | [download .h5ad file](scRNA_data/PAN.A01.v01.raw_count.20210429.NKT.embedding.h5ad) | [download scVI model](scRNA_data/scVI_models/scvi_NKT_model)

See our [code repository](https://github.com/Teichlab/Pan_fetal_immune/tree/master/src/utils/scArches_utils) for example use of scVI models for reference-based mapping of new single-cell RNA-seq experiments. 

### _In vitro_ T cell development

Single-cell RNA-seq dataset of _in vitro_ derived T cells, from Artificial Thymic Organoid protocol (10459 cells) 

[Interactive viewer](ATO_adata.cellxgene.h5ad) | [download .h5ad file](ATO_adata.h5ad)

## VDJ datasets

Datasets of combined single-cell RNA-seq profiles and adaptive immune receptor repertoire data.

- abTCR data | [download .h5ad file](scVDJ_data/PAN.A01.v01.raw_count.20210429.NKT.embedding.abTCR.h5ad)
- gdTCR data | [download .h5ad file](scVDJ_data/PAN.A01.v01.raw_count.20210429.NKT.embedding.gdTCR.h5ad)
- BCR data | [download .h5ad file](scVDJ_data/PAN.A01.v01.raw_count.20210429.LYMPHOID.embedding.BCR.h5ad)

## Spatial transcriptomics datasets (Visium 10X)

Datasets of spatial transcriptomics profiles for 3 developing tissues (18 pcw). 

- Fetal liver (3 slides) | [download .h5ad file](Visium_data/Visium10X_data_LI.h5ad) | [download high res images](Visium_data/LI_img/)
- Fetal thymus (3 slides) | [download .h5ad file](Visium_data/Visium10X_data_SP.h5ad) | [download high res images](Visium_data/TH_img/)
- Fetal spleen (4 slides) | [download .h5ad file](Visium_data/Visium10X_data_TH.h5ad) | [download high res images](Visium_data/SP_img/)

See [this notebook]() for example exploration and analysis of spatial gene expression data and stored outputs from spatial cell type mapping.

## Contact information
This data has been generated and processed as a collaboration between the [Teichmann Lab](http://www.teichlab.org/) and [Haniffa Lab](https://haniffalab.com/), based at the [Wellcome Sanger Insitute](https://www.sanger.ac.uk/) and [Newcastle University](https://www.ncl.ac.uk/) 

We welcome feedback and suggestions on how to improve this data portal to further facilitate re-use of our dataset. Please send any queries to Emma Dann `ed6<at>sanger.ac.uk` or Chenqu Suo `cs43<at>sanger.ac.uk`.
