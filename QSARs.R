#!/usr/bin/env Rscript
source ("tool.R")
source("MachinLearning.R")

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
modelPCR = 1
modelSVM = 1 
nbCV = 10
valcor = 0
logaff = 1
valcor = 0.7
maxdata = 80
cutoffAff = 0 

maxcompPLS = 20

#########################
#    PRINT PARAMETERS   #
#########################

print("=====PARAMETERS=====")
print (paste("Descriptors: ", pdesc, sep = ""))
print (paste("Predictive variable: ", pdesc, sep = ""))
print (paste("Folder out: ", pdesc, sep = ""))
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
daffinity = ddata[,c("CMPD_CHEMBLID", "STANDARD_VALUE")]

if(logaff == 1){
  daffinity[,2] = log10(daffinity[,2])
}
rownames(daffinity) = daffinity[,1]


print(dim(dglobal))
dglobalAff = cbind(dglobal, daffinity[rownames(dglobal),2])

##### IMPORTANT RENAME AFF  #####
#################################

colnames(dglobalAff) = c(colnames(dglobal), "Aff")
#print (dglobalAff[,dim(dglobalAff)[2]])

print(dim(dglobalAff))


# control distribution affinity
###############################

png(paste(prout, "histAffglobal.png", sep = ""), width = 400, height = 400)
hist(daffinity[,2], col = "grey", main = "Hist of log10 MIC", xlab = "log10 MIC", ylab = "Occurencies")
dev.off()

# sample data with fraction
ltraintest = samplingDataFraction(dglobalAff, 0.33)
controlDatasets(ltraintest, paste(prout, "CheckTrainTest", sep = ""))


# sampling data for CV
lgroupCV = samplingDataNgroup(dglobalAff, 10)
controlDatasets(lgroupCV, paste(prout, "ChecksamplingCV", sep = ""))


#Machine learning - PCR
#nbCp = PCRCV(lgroupCV, prout)
#PCRTrainTest(ltraintest[[1]], ltraintest[[2]], nbCp)


#Machine learning - PLS
nbCp = PLSCV(lgroupCV, prout)
PLSTrainTest(ltraintest[[1]], ltraintest[[2]], nbCp)




