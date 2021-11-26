This folder contains scripts for cell type mapping with cell2location. More in https://github.com/vitkl/pan_fetal_immune_mapping

### Single-organ reference

1. Subset scRNA-seq data to organs profiled with visium (for thymus we include TEC annotation from Park et al.):
```
python ./prep_c2l_reference.py --subset_organ SP
python ./prep_c2l_reference.py --subset_organ LI
python ./prep_c2l_reference.py --subset_organ TH --add_TECs True --keep_fetal_TECs True
```
2. Train reference model for single-cell data 

### Multi-organ reference

1. Subset scRNA-seq data to organs profiled with visium (for thymus we include TEC annotation from Park et al.):
```
python ./prep_c2l_reference_multiorgan.py --split_stroma True
```

