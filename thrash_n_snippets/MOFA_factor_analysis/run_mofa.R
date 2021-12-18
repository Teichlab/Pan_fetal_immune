### Run group mofa ###
library(MOFA2)
library(Matrix)
library(SingleCellExperiment)
library(scran)
library(glue)

split = "LYMPHOID"
indir <- glue("/nfs/team205/ed6/data/Fetal_immune/LMM_data/LMM_input_{split}_PBULK/")
mofa <- readRDS( glue('{indir}LYMPHOID_mofa_obj_organCorrected.RDS'))

## Prepare for training
data_opts <- get_default_data_options(mofa)
data_opts$use_float32 <- TRUE
data_opts$center_groups <- FALSE

model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 30
model_opts$

train_opts <- get_default_training_options(mofa)
train_opts$seed <- 2020
train_opts$convergence_mode <- "medium" # use "fast" for faster training
train_opts$stochastic <- FALSE

# mefisto_opts <- get_default_mefisto_options(mofa)
# mefisto_opts$warping <- FALSE
# mefisto_opts$sparseGP <- TRUE

mofa <- prepare_mofa(
  object = mofa,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
#   mefisto_options = mefisto_opts
) 

outfile <- glue('{indir}{split}_mofa_model_oneview_organCorrected.hdf5')
mofa_trained <- run_mofa(mofa, outfile = outfile)