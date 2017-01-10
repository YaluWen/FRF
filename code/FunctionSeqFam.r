library(MASS)
library(kinship2)
library(coxme)

# Simulation 1 #
source(paste(progdir,'GeneDrop.r',sep=''));
source(paste(progdir,'pedigree_gen_heter.r',sep=''));

source(paste(progdir,'SmallFun.r',sep=''));
source(paste(progdir,'SampleStrategiesF.r',sep=''));
source(paste(progdir,'DiscodeFinal.r',sep=''));
source(paste(progdir,'FGRF_pred.r',sep=''));
source(paste(progdir,'Simulation1.r',sep=''));

# Simulation 2 #
source(paste(progdir,'Simulation2.r',sep=''));
source(paste(progdir,'FGRF_pred_heter.r',sep=''));
source(paste(progdir,'DiscodeFinal.Heter.r',sep=''));
