#!/usr/bin/env Rscript
source ("tool.R")
source("dataManager.R")


################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
paff = args[2] #to take affinity or class
prout = args[3]
valcor = args[4]
maxquantile = args[5]
proptraintest = args[6]
logaff = args[7]

pdesc = "/home/aborrel/fluoroquinolones/results/desc/desc_compound.csv"
pdata = "/home/aborrel/fluoroquinolones/MIC_currated.csv"
prout = "/home/aborrel/fluoroquinolones/results/QSARS/"


valcor = 0.8
logaff = 1
maxquantile = 85
proptraintest = 0.15



##############################
# Process descriptors matrix #
##############################

dglobal = openData(pdesc, valcor, prout, c(1,2))
dglobal = dglobal[[1]]

rownames(dglobal) = dglobal[,1]
dglobal = dglobal[,-1]
dglobal = dglobal[,-1]# remove SMILES

print("==== Preprocessing ====")
print(paste("Data initial: dim = ", dim(dglobal)[1], dim(dglobal)[2], sep = " "))

##########
# filter #
##########

dglobal = delnohomogeniousdistribution(dglobal, maxquantile)
print(paste("Data after filtering: dim = ", dim(dglobal)[1], dim(dglobal)[2], sep = " "))

#######################
# order with affinity #
#######################
# Opening
daffinity = read.csv(pdata, sep = ",", header = TRUE)
rownames(daffinity) = daffinity[,1]
daffinity = daffinity[,-1]

# transform #
if(logaff == 1){
  daffinity = -log10(daffinity)
}

# merge with data descriptors and remove data remove from the manual curation
lID = intersect(rownames(daffinity), rownames(dglobal))
dglobal = dglobal[lID,]
daffinity = daffinity[lID,]

##################
# divide dataset #
##################
ltraintest = samplingDataFraction(dglobal, proptraintest)
dtrain = ltraintest[[1]]
dtest = ltraintest[[2]]


############################
# Add affinity by bacteria #
############################

for (bacteria in colnames(daffinity)){
  lcontrol = list()
  
  # training set
  Aff = daffinity[rownames(dtrain),bacteria]
  dtrainglobal = cbind(dtrain, Aff)
  lcontrol = append(lcontrol, dtrainglobal)
  write.csv(dtrainglobal, paste(prout, "train_", bacteria, ".csv", sep = ""))
  
  # test set
  Aff = daffinity[rownames(dtest),bacteria]
  dtestglobal = cbind(dtest, Aff)
  lcontrol = append(lcontrol, dtestglobal)
  write.csv(dtestglobal, paste(prout, "test_", bacteria, ".csv", sep = ""))
  
  lcontrol = list(dtrainglobal, dtestglobal)
  controlDatasets(lcontrol, paste(prout, "qualitySplit_", bacteria, sep = ""))
}







