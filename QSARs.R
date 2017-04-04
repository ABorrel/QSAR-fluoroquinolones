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

#args <- commandArgs(TRUE)
#pdesc = args[1]
#pdata = args[2] #to take affinity or class
#prout = args[3]

# regression model
#modelPLS = args[4]
#modelSVM = args[5]

# Machine learning parameters
#nbCV = args[6]
#valcor = 0.7

#pdesc = "/home/aborrel/fluoroquinolones/results/desc/desc_Staphylococcus-aureus.csv"
#pdata = "/home/aborrel/fluoroquinolones/results/compound_filtered_MIC_MIC_Staphylococcus-aureus.csv"
#prout = "/home/aborrel/fluoroquinolones/results/model/Staphylococcus-aureus/" 


pdesc = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/desc/desc.csv"
pdata = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/Carbamate.csv"
prout = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/desc/OD" 



# model regression #
####################
modelPCRreg = 0
modelPLSreg = 0
modelSVMreg = 0
modelRFreg = 0
modelCartreg = 1
chemmodlabreg = 0


# model classification #
########################
modelPCRclass = 0
modelPLSclass = 0
modelSVMclass = 0
modelRFclass = 0
modelLDA = 0
modelCartclass = 1
chemmodlabclass = 0


# Paramaters global #
#####################
nbCV = 10
valcor = 0
logaff = 0
maxquantile = 85
cutoffAff = "med" 
proptraintest = 0.20


#########################
#    PRINT PARAMETERS   #
#########################

print("=====PARAMETERS=====")
print (paste("Descriptors: ", pdesc, sep = ""))
print (paste("Predictive variable: ", pdata, sep = ""))
print (paste("Folder out: ", prout, sep = ""))
print (paste("Cutoff activity in case of binary model: ", cutoffAff, sep = ""))
print(paste("Nb of CV: ", nbCV, sep = ""))
print("")

print("=====Data filtering=====")
print(paste("Cor value for clustering: ", valcor, sep = ""))
print(paste("Max percentage data by decile: ", maxdata, sep = ""))
print("")

print("=====Machine learning=====")
print("---Regression model----")
print(paste("PCR: ", modelPCRreg, sep = ""))
print(paste("PLS: ", modelPLSreg, sep = ""))
print(paste("SVM: ", modelSVMreg, sep = ""))
print(paste("CART: ", modelCartreg, sep = ""))
print(paste("RF: ", modelRFreg, sep = ""))
print(paste("Chemmodlab: ", chemmodlabreg, sep = ""))

print("---Classification----")
print(paste("PCR: ", modelPCRclass, sep = ""))
print(paste("PLS: ", modelPLSclass, sep = ""))
print(paste("LDA: ", modelLDA, sep = ""))
print(paste("SVM: ", modelSVMclass, sep = ""))
print(paste("CART: ", modelCartclass, sep = ""))
print(paste("RF: ", modelRFclass, sep = ""))
print(paste("Chemmodlab: ", chemmodlabclass, sep = ""))



##############################
# Process descriptors matrix #
##############################
dglobal = openData(pdesc, valcor, prout, c(1,2))
dglobal = dglobal[[1]]

rownames(dglobal) = dglobal[,1]
dglobal = dglobal[,-1]
dglobal = dglobal[,-1]# remove SMILES

# filter
dglobal = delnohomogeniousdistribution(dglobal, maxdata)

#######################
# order with affinity #
#######################
# Opening
ddata = read.table(pdata, sep = "\t", header = TRUE)
print(dim(ddata))
#daffinity = ddata[,c("CMPD_CHEMBLID", "STANDARD_VALUE")] # !!!!!!!!!!!!!
daffinity = ddata[,c("ID", "Odavg")]
#daffinity = ddata[,c("ID", "FLDavg")]

#print(daffinity)
rownames(daffinity) = daffinity[,1]

dglobalAff = cbind(dglobal[rownames(daffinity),], activity=daffinity[,2])


# remove NA
dglobalAff = na.omit(dglobalAff)


##### IMPORTANT RENAME AFF  #####
#################################
colnames(dglobalAff) = c(colnames(dglobal), "Aff")

# control distribution affinity and log10
#########################################

png(paste(prout, "histAffglobal.png", sep = ""), width = 800, height = 400)
par(mfrow = c(1,2))
hist(dglobalAff[,c("Aff")], col = "grey", main = "Hist aff", xlab = "Aff", ylab = "Occurencies")

if(logaff == 1){
  dglobalAff[,c("Aff")] = log10(dglobalAff[,c("Aff")])
  print(dglobalAff[,c("Aff")])
  hist(dglobalAff[,c("Aff")], col = "grey", main = "Hist of log10 Aff", xlab = "log10 aff", ylab = "Occurencies")
}

dev.off()

# sample data with fraction #
#############################
ltraintest = samplingDataFraction(dglobalAff, proptraintest)
controlDatasets(ltraintest, paste(prout, "CheckTrainTest", sep = ""))

# sampling data for CV #
########################
lgroupCV = samplingDataNgroup(dglobalAff, nbCV)
controlDatasets(lgroupCV, paste(prout, "ChecksamplingCV", sep = ""))


##### REGRESSION MODELS #########
#################################


### PCR ####
############

if (modelPCRreg == 1){
  nbCp = PCRCV(lgroupCV, prout)
  PCRTrainTest(ltraintest[[1]], ltraintest[[2]], nbCp)
}

### PLS  ####
#############

if (modelPLSreg == 1){
  nbCp = PCRgridCV(lgroupCV, prout)
  PCRCV (lgroupCV, nbCp, prout)
  PLSTrainTest(ltraintest[[1]], ltraintest[[2]], nbCp)
}

### SVM ###
###########

if(modelSVMreg == 1){
  vgamma = 2^(-1:1)
  vcost = 2^(2:4)
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
  
  dtrain = ltraintest[[1]]
  dtest = ltraintest[[2]]
  lfoldtrain = samplingDataNgroup(dtrain, nbCV)
  parameters = RFGridRegCV(vntree, vmtry, lfoldtrain, paste(prout, "trainReg", sep = ""))
  RFreg(dtrain, dtest,parameters[[1]], parameters[[2]], prout)
}

#############
# CHEMMOLAB #
#############

if(chemmodlabreg == 1){
  dchem = cbind(dglobalAff[,dim(dglobalAff)[2]],dglobalAff[,-dim(dglobalAff)[2]] )
  colnames(dchem)[1] = "Aff"
  
  pdf(paste(prout, "Reg_chemmolab.pdf", sep = ""))
  fit = ModelTrain(data = dchem, ids = FALSE)
  CombineSplits(fit, metric = "R2")
  CombineSplits(fit, metric = "rho")
  CombineSplits(fit, metric = "auc")
  dev.off()
}


############
#   CART   #
############

if(modelCartreg == 1){
  print("dd")
  
}


####################
## CLASSIFICATION ##
####################

# DATA TRANSFORMATION -> 0 and 1 #
##################################

if (cutoffAff == "med"){
  cutoffAff = median(dglobalAff[,c("Aff")])
}else{
  cutoffAff = as.double(cutoffAff)
}
i1 = which(dglobalAff[,c("Aff")] > cutoffAff)
dglobalAff[,c("Aff")] = 0 
dglobalAff[i1,c("Aff")] = 1

lgroupCV = samplingDataNgroupClass(dglobalAff, nbCV, "Aff")
ltraintest = samplingDataFraction(dglobalAff, proptraintest) 


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
