#!/usr/bin/env Rscript
source ("tool.R")
source("MachinLearning.R")
source("performance.R")
source("dataManager.R")

library(chemmodlab)
library (rpart)
################
#     MAIN     #
################

#./QSARs.R /home/aborrel/fluoroquinolones/results/QSARS/train_Escherichia.coli.csv /home/aborrel/fluoroquinolones/results/QSARS/test_Escherichia.coli.csv /home/aborrel/fluoroquinolones/results/desc_analysis/0.8/Table_hclust_ward.D2_gap_stat.csv /home/aborrel/fluoroquinolones/results/QSARS//Escherichia.coli/ 1 0 >/home/aborrel/fluoroquinolones/results/QSARS//Escherichia.coli/perf.txt


args <- commandArgs(TRUE)
ptrain = args[1]
ptest = args[2]
pcluster = args[3]
prout = args[4]
nbCV = as.integer(args[5])

# to test
#ptrain = "/home/aborrel/fluoroquinolones/results/QSARS/train_Escherichia.coli.csv"
#ptest = "/home/aborrel/fluoroquinolones/results/QSARS/test_Escherichia.coli.csv"
#pcluster = "/home/aborrel/fluoroquinolones/results/desc_analysis/0.8/Table_hclust_ward.D2_gap_stat.csv"
#prout = "/home/aborrel/fluoroquinolones/results/QSARS//Escherichia.coli/"

# cross validation 10
#nbCV = 10


# model regression #
####################
modelPCRreg = 1
modelPLSreg = 1
modelSVMreg = 1
modelRFreg = 1
modelCartreg = 1
chemmodlabreg = 1


#########################
#    PRINT PARAMETERS   #
#########################

print("=====PARAMETERS=====")
print (paste("Train csv: ", ptrain, sep = ""))
print (paste("Test csv: ", ptest, sep = ""))
print (paste("Folder out: ", prout, sep = ""))
print(paste("Nb of CV: ", nbCV, sep = ""))
print("")

print("=====Machine learning=====")
print("---Regression model----")
print(paste("PCR: ", modelPCRreg, sep = ""))
print(paste("PLS: ", modelPLSreg, sep = ""))
print(paste("SVM: ", modelSVMreg, sep = ""))
print(paste("CART: ", modelCartreg, sep = ""))
print(paste("RF: ", modelRFreg, sep = ""))
print(paste("Chemmodlab: ", chemmodlabreg, sep = ""))
print("")



##############################
# Process descriptors matrix #
##############################

# training set
dtrain = read.csv(ptrain, header = TRUE)
rownames(dtrain) == dtrain[,1]
dtrain = dtrain[,-1]

# test set
dtest = read.csv(ptest, header = TRUE)
rownames(dtest) == dtest[,1]
dtest = dtest[,-1]

print("==== Dataset ====")
print(paste("Data train: dim = ", dim(dtrain)[1], dim(dtrain)[2], sep = " "))
print(paste("Data test: dim = ", dim(dtest)[1], dim(dtest)[2], sep = " "))
print("")

# sampling data for CV #
########################
lgroupCV = samplingDataNgroup(dtrain, nbCV)
controlDatasets(lgroupCV, paste(prout, "ChecksamplingCV", nbCV, sep = ""))


##### REGRESSION MODELS #########
#################################
print("**************************")
print("*****  REGRESSSION   *****")
print("**************************")

### PCR ####
############

if (modelPCRreg == 1){
  nbCp = PCRgridCV(lgroupCV, prout)
  print(nbCp)
  PCRCV(lgroupCV, nbCp, prout)
  PCRTrainTest(dtrain, dtest, nbCp)
}

### PLS  ####
#############

if (modelPLSreg == 1){
  nbcp = PLSCV(lgroupCV, prout)
  PLSTrainTest(dtrain, dtest, nbcp)
}

### SVM ###
###########

if(modelSVMreg == 1){
  vgamma = 2^(-1:1)
  vcost = 2^(2:8)
  SVMRegCV(lgroupCV, vgamma, vcost, prout)
  
}

######
# RF #
######

if (modelRFreg == 1){
  vntree = c(10,50,100,200,500, 1000)
  vmtry = c(1,2,3,4,5,10,15,20, 25, 30)
  
  parameters = RFGridRegCV(vntree, vmtry, lgroupCV, prout)
  RFregCV(lgroupCV, parameters[[1]], parameters[[2]], prout)
  
  RFreg(dtrain, dtest, parameters[[1]], parameters[[2]], prout)
  
  #### case of external set #
  ###########################
  #print ("======External SET======")
  #dtest = read.table(ptest, sep = "\t", header = TRUE)
  #Aff = rep(0, dim(dtest)[1])
  #rownames(dtest) = dtest[,1]
  #dtest = dtest[,-1]
  #dtest = dtest[,-1]
  #dtest = cbind(dtest, Aff)
  #  
  #lfoldtrain = samplingDataNgroup(dglobalAff, nbCV)
  #parameters = RFGridRegCV(vntree, vmtry, lfoldtrain, paste(prout, "trainReg", sep = ""))
  #RFreg(dtrain, dtest, parameters[[1]], parameters[[2]], prout)
  
}



############
#   CART   #
############

if(modelCartreg == 1){
  CARTRegCV(lgroupCV, prout)
  CARTreg(dtrain, dtest, prout)
}


#############
# CHEMMOLAB #
#############

if(chemmodlabreg == 1){
  dchem = cbind(dglobalAff[,dim(dglobalAff)[2]],dglobalAff[,-dim(dglobalAff)[2]] )
  colnames(dchem)[1] = "Aff"
  
  pdf(paste(prout, "Reg_chemmolab.pdf", sep = ""))
  fit = ModelTrain(dchem, ids = FALSE)
  CombineSplits(fit, metric = "R2")
  CombineSplits(fit, metric = "rho")
  dev.off()
}



