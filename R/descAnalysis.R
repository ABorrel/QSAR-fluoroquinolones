#!/usr/bin/env Rscript
source ("tool.R")
source("cardMatrix.R")
source("PCAplot.R")
source("dendocircular.R")
source("clustering.R")
source("distributions.R")

################
#     MAIN     #
################

args <- commandArgs(TRUE)
pdesc = args[1]
pdata = args[2] #to take affinity
prout = args[3]
valcor = as.double(args[4])
maxQuantile = as.double(args[5])
logaff = as.integer(args[6])
plotPCA = as.integer(args[7])
corMatrix = as.integer(args[8])
histplot = as.integer(args[9])
circularDendo = as.integer(args[10])
optimal_clustering = as.integer(args[11])


#pdesc = "./../../results/DESC/desc_compound.csv"
#pdata = "./../../results/dataset/MIC-curated_mol.csv"
#prout = "./../../results/Clustering-0.85-85/"

#plotPCA = 0
#corMatrix = 0
#histplot = 0
#circularDendo = 1
#valcor = 0.85
#maxQuantile = 85
#logaff = 1
#optimal_clustering = 1


# Process descriptors matrix #
##############################
dglobal = openData(pdesc, valcor, prout, c(1,2))
dglobal = dglobal[[1]]

rownames(dglobal) = dglobal[,1]
dglobal = dglobal[,-1]
dglobal = dglobal[,-1]

# remove descriptor with a distribution on one quantile
dglobal = delnohomogeniousdistribution(dglobal, maxQuantile)

#######################
# order with affinity #
#######################
# Opening
daffinity = read.csv(pdata, sep = "\t", header = TRUE)
rownames(daffinity) = daffinity[,1]
daffinity = daffinity[,-1]
daffinity = daffinity[,-which(colnames(daffinity) == "SMILES")] # remove SMILES 
print(dim(daffinity))



if(logaff == 1){
  daffinity = -log10(daffinity)
  multipleHist(daffinity, prout)
}

# remove no selected affinity or bad quality from manual curation
print(dim(dglobal))
lID = intersect(rownames(daffinity), rownames(dglobal))

dglobal = dglobal[lID,]
daffinity = daffinity[lID,]

write.csv(dglobal, paste(prout, "globalTable.csv", sep = ""))

#ord = NULL
#for(i in seq(1,dim(dglobal)[1])){
#  ipos = which(daffinity[,1] == rownames(dglobal)[i])
#  ord = append (ord,ipos)
#}

#daffinity = daffinity[ord,]


#orderaff = order(daffinity[,2],decreasing=T)
#daffinity = daffinity[orderaff,]


#dglobal = dglobal[orderaff,]
#print(rownames(dglobal))
#print(colnames(dglobal))

#print(dim(dglobal))
#print(dim(daffinity))

#print(rownames(dglobal))
#print(rownames(daffinity))



if (corMatrix == 1){
  cardMatrixCor(cor(dglobal), paste(prout, "matrixCor_", valcor, sep = ""), 6)
  
  # matrix cor with pMIC
  dcor = NULL
  lbact = colnames(daffinity)
  ldesc = colnames(dglobal)
  for(bact in lbact){
    lcor = NULL
    for(desc in ldesc){
      print(bact)
      print (desc)
      
      corval = cor(daffinity[,bact], dglobal[,desc])
      lcor = append(lcor, corval)
    }
    dcor = cbind(dcor, lcor)
  }
  print (dcor)
  rownames(dcor) = ldesc
  colnames(dcor) = lbact
  
  write.csv(dcor, file = paste(prout, "corMIC-Desc.csv"))
}


if(plotPCA == 1){
  PCAplot(dglobal, daffinity, prout)
}

if (histplot == 1){
  histDataOne(data1 = dglobal, paste(prout, "histDesc_", valcor, ".pdf", sep = ""))
}

if (circularDendo == 1){
  MultipleDendogramCircle(dglobal, daffinity, prout)
}


if (optimal_clustering ==1 ){
  lmetclustering = c("hclust", "kmeans")
  lmetagregation = c("ward.D2", "ward.D", "complete", "single", "average")
  lmetoptimal = c("silhouette", "wss", "gap_stat")
  for(metclustering in lmetclustering){
    if (metclustering == "kmeans"){
      lmetagregation = c("ward.D2")
    }else{
      lmetagregation = c("ward.D2", "complete", "single", "average")
    }
    for (metagregation in lmetagregation){
      for(metoptimal in lmetoptimal){
        pclust = paste(prout, "/", metclustering, "-", metagregation, "-", metoptimal, "/", sep = "")
        dir.create(pclust)
        dclust = optimalCluters(dglobal, pclust, metclustering, metoptimal, metagregation)
        print(dim(dclust))
        dMCluster = NULL
        lclust = unique(dclust[,2])
        for(clust in lclust){
          daffclust = daffinity[dclust[which(dclust[,2] == clust),1],]
          
          M1 = mean(daffclust[,1])
          SD1 = sd(daffclust[,1])
          M2 = mean(daffclust[,2])
          SD2 = sd(daffclust[,2])
          M3 = mean(daffclust[,3])
          SD3 = sd(daffclust[,3])
          M4 = mean(daffclust[,4])
          SD4 = sd(daffclust[,4])



          dMCluster = rbind(dMCluster, c(clust, M1, SD1, M2, SD2, M3, SD3, M4, SD4))
        }
        lbact = colnames(daffinity)
        colnames(dMCluster) = c("cluster", paste("M_", lbact[1], sep = ""),  paste("SD_", lbact[1], sep = ""),  paste("M_", lbact[2], sep = ""),  paste("SD_", lbact[2], sep = ""), paste("M_", lbact[3], sep = ""),  paste("SD_", lbact[3], sep = ""),  paste("M_", lbact[4], sep = ""),  paste("SD_", lbact[4], sep = ""))
        write.csv(dMCluster, paste(pclust, "Mclusters_", metclustering, metoptimal, metagregation, ".csv", sep = ""))
      }
    }
  }
}

