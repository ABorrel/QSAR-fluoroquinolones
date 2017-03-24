#!/usr/bin/env Rscript
source ("tool.R")
source("cardMatrix.R")
source("PCAplot.R")
source("dendocircular.R")



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pdata = args[2] #to take affinity
prout = args[3]
valcor = as.double(args[4])
plotPCA = args[5]
corMatrix = args[6]
histplot = args[7]
circularDendo = args[8]
logaff = args[9]

#pdesc = "/home/aborrel/fluoroquinolones/results/desc/desc_Staphylococcus-aureus.csv"
#pdata = "/home/aborrel/fluoroquinolones/results/compound_filtered_MIC_MIC_Staphylococcus-aureus.csv"
#prout = "/home/aborrel/fluoroquinolones/results/desc/Staphylococcus-aureus" 

#plotPCA = 1
#corMatrix = 0
#histplot = 1
#circularDendo = 1
#valcor = 0.70


# Process descriptors matrix #
##############################
dglobal = openData(pdesc, valcor, prout, c(1,2))
dglobal = dglobal[[1]]

rownames(dglobal) = dglobal[,1]
dglobal = dglobal[,-1]
dglobal = dglobal[,-1]

# order with affinity #
#######################
# Opening
ddata = read.table(pdata, sep = "\t", header = TRUE)
print(dim(ddata))
daffinity = ddata[,c("CMPD_CHEMBLID", "STANDARD_VALUE")]

#print(dim(daffinity))
#print(dim(dglobal))



if(logaff == 1){
  daffinity[,2] = log10(daffinity[,2])
}


rownames(daffinity) = daffinity[,1]


ord = NULL
for(i in seq(1,dim(dglobal)[1])){
  ipos = which(daffinity[,1] == rownames(dglobal)[i])
  ord = append (ord,ipos)
}

print (ord)
daffinity = daffinity[ord,]


orderaff = order(daffinity[,2],decreasing=T)
daffinity = daffinity[orderaff,]


dglobal = dglobal[orderaff,]
#print(rownames(dglobal))
#print(colnames(dglobal))

#print(dim(dglobal))
#print(dim(daffinity))

#print(rownames(dglobal))
#print(rownames(daffinity))

if (corMatrix == 1){
  cardMatrixCor(cor(dglobal), paste(prout, "matrixCor_", valcor, sep = ""), 6)
}

if(plotPCA == 1){
  PCAplot(dglobal, paste(prout, "global_", valcor, sep = ""))
}

if (histplot == 1){
  histDataOne(data1 = cbind(dglobal, daffinity[,2]), paste(prout, "histDesc_", valcor, ".pdf", sep = ""))
}

if (circularDendo == 1){
  dendogramCircle(dglobal, daffinity, paste(prout, "dendo_", valcor, ".png", sep = ""))
}
