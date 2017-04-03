#!/usr/bin/env Rscript
source ("tool.R")
source("MachinLearning.R")
source("performance.R")
source("dataManager.R")

library(chemmodlab)

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

pdesc = "/home/aborrel/fluoroquinolones/results/desc/desc_Staphylococcus-aureus.csv"
pdata = "/home/aborrel/fluoroquinolones/results/compound_filtered_MIC_MIC_Staphylococcus-aureus.csv"
prout = "/home/aborrel/fluoroquinolones/results/model/Staphylococcus-aureus/" 



pdesc = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/desc/desc.csv"
pdata = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/Carbamate.csv"
prout = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/desc/" 


modelPCR = 1
modelSVM = 1 
nbCV = 10
valcor = 0
logaff = 0
valcor = 0.9
maxdata = 85
cutoffAff = "med" 
proptraintest = 0.20
jeremyLib = 1

maxcompPLS = 20

#########################
#    PRINT PARAMETERS   #
#########################

print("=====PARAMETERS=====")
print (paste("Descriptors: ", pdesc, sep = ""))
print (paste("Predictive variable: ", pdata, sep = ""))
print (paste("Folder out: ", prout, sep = ""))
print (paste("Cutoff activity in case of binary model: ", cutoffAff, sep = ""))
print("")

print("=====Data filtering=====")
print(paste("Cor value for clustering: ", valcor, sep = ""))
print(paste("Max percentage data by decile: ", maxdata, sep = ""))
print("")

print("=====Machine learning=====")
print(paste("Nb of CV: ", nbCV, sep = ""))
print("---Regression model----")
print(paste("PCR: ", modelPCR, sep = ""))
print(paste("PLS: ", modelPCR, sep = ""))
print(paste("SVM: ", modelPCR, sep = ""))
print(paste("CART: ", modelPCR, sep = ""))
print(paste("Chemmodlab: ", jeremyLib, sep = ""))

print("---Classification----")




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
#daffinity = ddata[,c("ID", "Odavg")]
daffinity = ddata[,c("ID", "FLDavg")]

#print(daffinity)
rownames(daffinity) = daffinity[,1]

dglobalAff = cbind(dglobal[rownames(daffinity),], activity=daffinity[,2])


# remove NA
dglobalAff = na.omit(dglobalAff)

# rename col
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

##### IMPORTANT RENAME AFF  #####
#################################

# sample data with fraction
ltraintest = samplingDataFraction(dglobalAff, proptraintest)
controlDatasets(ltraintest, paste(prout, "CheckTrainTest", sep = ""))

# sampling data for CV
lgroupCV = samplingDataNgroup(dglobalAff, nbCV)
controlDatasets(lgroupCV, paste(prout, "ChecksamplingCV", sep = ""))


##### REGRESSION MODELS #########
#################################


### PCR ####
############

#Machine learning - PCR
#nbCp = PCRCV(lgroupCV, prout)
#PCRTrainTest(ltraintest[[1]], ltraintest[[2]], nbCp)


### PLS  ####
#############

#Machine learning - PLS
#nbCp = PLSCV(lgroupCV, prout)
#PLSTrainTest(ltraintest[[1]], ltraintest[[2]], nbCp)


######
# RF #
######

vntree = c(10,50,100,200,500, 1000)
vmtry = c(1,2,3,4,5,10,15,20, 25, 30)

# based on CV - 5 or 10  #
##########################
parameters = RFGridRegCV(vntree, vmtry, lgroupCV, prout)
RFregCV(lgroupCV, parameters[[1]], parameters[[2]], prout)

# based on external set #
#########################
#dtrain = ltraintest[[1]]
#dtest = ltraintest[[2]]
#lfoldtrain = samplingDataNgroup(dtrain, 5)
#parameters = RFGridRegCV(vntree, vmtry, lfoldtrain, paste(prout, "trainReg", sep = ""))

#RFreg(dtrain, dtest, ntree, mtry, prout)



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

vntree = c(10,50,100,200,500, 1000)
vmtry = c(1,2,3,4,5,10,15,20, 25, 30)

# based on CV - 5 or 10  #
##########################
parameters = RFGridClassCV(vntree, vmtry, lgroupCV, prout)
RFClassCV(lgroupCV, parameters[[1]], parameters[[2]], prout)

# based on external set #
#########################
#dtrain = ltraintest[[1]]
#dtest = ltraintest[[2]]

#lfoldtrain = samplingDataNgroup(dtrain, nbCV)
#parameters = RFGridClassCV(vntree, vmtry, lfoldtrain, paste(prout, "trainClass", sep = ""))
#RFClass(dtrain, dtest, parameters[[1]], parameters[[2]], paste(prout, "traintestperf", sep = ""))



##### JEREMY PERFORMANCE  #####
###############################
# USEFULL TO COMPARE SEVERAL METHODS

#pdf(paste(prout, "FLVbin_optimization_chemmolab.pdf", sep = ""))

# put data
#dglobalAff = cbind(dglobalAff[,dim(dglobalAff)[2]],dglobalAff[,-dim(dglobalAff)[2]] )
#colnames(dglobalAff)[1] = "Aff"

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

#CombineSplits(fit, metric = "error rate")
#CombineSplits(fit, metric = "sensitivity")
#CombineSplits(fit, metric = "specificity")
#dev.off()
