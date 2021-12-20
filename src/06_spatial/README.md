This folder contains the main scripts for analysis of Visium slides and spatial cell type mapping from scRNA-seq data with [cell2location](https://cell2location.readthedocs.io/en/latest/)

- `c2l_mapping/` - contains the main scripts for cell type mapping with cell2location (more versions in https://github.com/vitkl/pan_fetal_immune_mapping)
- `spatial_c2l_output_EDA.ipynb` - downstream analysis and plotting notebook

## Spatial celltype mapping workflow

1. Subset scRNA-seq data to organs profiled with visium (for thymus we include TEC annotation from Park et al.):
```
python c2l_mapping/prep_c2l_reference.py --subset_organ SP
python c2l_mapping/prep_c2l_reference.py --subset_organ LI
python c2l_mapping/prep_c2l_reference.py --subset_organ TH --add_TECs True --keep_fetal_TECs True
```

2. Train reference model for single-cell data - see `c2l_mapping/cell2location_estimating_signatures_*.ipynb`

3. Train spatial mapping model - see ``c2l_mapping/cell2location_*.ipynb``


