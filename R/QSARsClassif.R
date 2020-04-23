
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
nbCV = as.integer(args[7])

# to test
ptrain = "/home/aborrel/fluoroquinolones/results/QSARS/train_Escherichia.coli.csv"
ptest = "/home/aborrel/fluoroquinolones/results/QSARS/test_Escherichia.coli.csv"
pcluster = "/home/aborrel/fluoroquinolones/results/desc_analysis/0.8/Table_hclust_ward.D2_gap_stat.csv"
prout = "/home/aborrel/fluoroquinolones/results/QSARS//Escherichia.coli/"

# cross validation 10
nbCV = 10



# model classification #
########################
modelPCRclass = 0
modelPLSclass = 0
modelSVMclass = 0
modelRFclass = 0
modelLDA = 0
modelCartclass = 0
chemmodlabclass = 0


#########################
#    PRINT PARAMETERS   #
#########################

print("=====PARAMETERS=====")
print (paste("Train csv: ", ptrain, sep = ""))
print (paste("Test csv: ", ptest, sep = ""))
print (paste("Folder out: ", prout, sep = ""))
print (paste("Cutoff activity in case of binary model: ", cutoffAff, sep = ""))
print(paste("Nb of CV: ", nbCV, sep = ""))
print("")

print("=====Machine learning=====")
print("---Classification----")
print(paste("PCR: ", modelPCRclass, sep = ""))
print(paste("PLS: ", modelPLSclass, sep = ""))
print(paste("LDA: ", modelLDA, sep = ""))
print(paste("SVM: ", modelSVMclass, sep = ""))
print(paste("CART: ", modelCartclass, sep = ""))
print(paste("RF: ", modelRFclass, sep = ""))
print(paste("Chemmodlab: ", chemmodlabclass, sep = ""))
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


############# NEVER USED NEED TO DEVELOP ################
#########################################################



# DATA TRANSFORMATION -> 0 and 1 #
##################################

if (cutoffAff == "med"){
  cutoffAff = median(dglobalAff[,c("Aff")])
}else{
  cutoffAff = as.double(cutoffAff)
}
print (paste("Cutoff -> ", cutoffAff, sep = ""))
i1 = which(dglobalAff[,c("Aff")] > cutoffAff)
dglobalAff[,c("Aff")] = 0 
dglobalAff[i1,c("Aff")] = 1

lgroupCV = samplingDataNgroupClass(dglobalAff, nbCV, "Aff")

#ltraintest = sampligDataFractionCluster(dglobalAff, proptraintest, )
ltraintest = samplingDataFraction(dglobalAff, proptraintest) 



####################
## CLASSIFICATION ##
####################

print("**************************")
print("***  CLASSIFICATION  *****")
print("**************************")



#########
#  RF  #
########

if (modelRFclass == 1){
  
  vntree = c(10,50,100,200,500, 1000)
  vmtry = c(1,2,3,4,5,10,15,20, 25, 30)
  parameters = RFGridClassCV(vntree, vmtry, lgroupCV, prout)
  RFClassCV(lgroupCV, parameters[[1]], parameters[[2]], prout)
  
  dtrain = ltraintest[[1]]
  dtest = ltraintest[[2]]
  
  lfoldtrain = samplingDataNgroup(dtrain, nbCV)
  parameters = RFGridClassCV(vntree, vmtry, lfoldtrain, paste(prout, "trainClass", sep = ""))
  RFClass(dtrain, dtest, parameters[[1]], parameters[[2]], paste(prout, "traintestperf", sep = ""))
  
  #### case of external set #
  ###########################
  print ("======External SET======")
  dtest = read.table(ptest, sep = "\t", header = TRUE)
  Aff = rep(0, dim(dtest)[1])
  rownames(dtest) = dtest[,1]
  dtest = dtest[,-1]
  dtest = dtest[,-1]
  dtest = cbind(dtest, Aff)
  
  lfoldtrain = samplingDataNgroup(dglobalAff, nbCV)
  parameters = RFGridClassCV(vntree, vmtry, lfoldtrain, paste(prout, "trainReg", sep = ""))
  RFClass(dtrain, dtest, parameters[[1]], parameters[[2]], prout)
  
}

############
#   CART   #
############

if(modelCartclass == 1){
  CARTClassCV(lgroupCV, prout)
}

#######
# SVM #
#######

if(modelSVMclass == 1){
  vgamma = 2^(-1:1)
  vcost = 2^(2:4)
  SVMClassCV(lgroupCV, vgamma, vcost, prout)
  
}

#############
# CHEMMOLAB #
#############


if(chemmodlabclass == 1){
  dchem = cbind(dglobalAff[,dim(dglobalAff)[2]],dglobalAff[,-dim(dglobalAff)[2]] )
  colnames(dchem)[1] = "Aff"
  
  pdf(paste(prout, "class_chemmolab.pdf", sep = ""))
  fit = ModelTrain(data = dchem, ids = FALSE)
  CombineSplits(fit, metric = "error rate")
  CombineSplits(fit, metric = "sensitivity")
  CombineSplits(fit, metric = "specificity")
  dev.off()
}




#pdf(paste(prout, "FLVbin_optimization_chemmolab.pdf", sep = ""))

# put data


#medaff = median(dglobalAff[,1])
#print (medaff)
#i1 = which(dglobalAff[,1] >= medaff)
#print(i1)
#dglobalAff[,1] = 0 
#dglobalAff[i1,1] = 1


#print(dglobalAff[,1])
# remove outlier aff > 20
#dglobalAff = dglobalAff[-which(dglobalAff[,1] > 20), ]


#fit = ModelTrain(data = dglobalAff, ids = FALSE)

#params <- MakeModelDefaults(n = nrow(dglobalAff), p = ncol(dglobalAff) -1, nfolds = 10, classify = T)
#params$SVM$gamma <- .5
#params$SVM$cost <- 4


#fit = ModelTrain(data = dglobalAff, ids = FALSE, user.params = params)

#CombineSplits(fit, metric = "error rate")
#CombineSplits(fit, metric = "sensitivity")
#CombineSplits(fit, metric = "specificity")
#dev.off()
