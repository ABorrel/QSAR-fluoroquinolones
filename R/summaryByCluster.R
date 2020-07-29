#!/usr/bin/env Rscript



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pcluster = args[1]
pMIC = args[2] #to take affinity
prout = args[3]

  
pcluster = "./../../results/Clustering_selected/hclust-ward.D2-gap_stat/Table_hclust_ward.D2_gap_stat.csv"
pMIC = "./../../results/dataset/MIC-curated_mol.csv"
prout = "./../../results/SumCluster-0.85-85/"


dcluster = read.csv(pcluster, sep = ",", header = TRUE)
rownames(dcluster) = dcluster[,1]

dMIC = read.csv(pMIC, sep = "\t", header = TRUE)
rownames(dMIC) = dMIC[,1]
dMIC = dMIC[,-1]
dMIC = dMIC[,-1] # remove SMILES too
dMIC = -log10(dMIC)

lcluster = unique(dcluster[,2])
lcluster =  seq(1,length(lcluster))

lID = intersect(rownames(dMIC), rownames(dcluster))
dMIC = dMIC[lID,]
dcluster = dcluster[lID,]

dMcluster = NULL


lbacteria = colnames(dMIC)

for(cluster in lcluster){
  lichem = which(dcluster$cluster == cluster)
  dtemp = dMIC[lichem,]
  
  lw = NULL
  for (bacteria in lbacteria){
    M = mean(dtemp[,bacteria])
    SD = sd(dtemp[,bacteria])
    lw = append(lw, M)
    lw = append(lw, SD)
  }
  dMcluster = rbind(dMcluster, lw)
}

eff = table(dcluster$cluster)
dMcluster = cbind(dMcluster, eff)

rownames(dMcluster) = lcluster
colnames(dMcluster) = c("M-Pseudomonas.aeruginosa", "SD-Pseudomonas.aeruginosa", "M-Staphylococcus.aureus", "SD-Staphylococcus.aureus", "M-Escherichia.coli", "SD-Escherichia.coli", "M-Streptococcus.pneumoniae", "SD-Streptococcus.pneumoniae", "Effective")
write.csv(dMcluster, paste(prout, "SumTableM-SD.csv", sep = ""))


dprobcluster = NULL

for(cluster in lcluster){
  lichem = which(dcluster$cluster == cluster)
  dtemp = dMIC[lichem,]
  
  dprob = c(0,0,0,0)
  for(ichem in seq(1,dim(dtemp)[1])){
    print(ichem)
    maxMIC = max(dtemp[ichem,])
    ibact = which(maxMIC == dtemp[ichem,])
    dprob[ibact] = dprob[ibact] + 1
  }
  dprob = round(dprob/sum(dprob),2)*100
  dprobcluster = rbind(dprobcluster, dprob)
}


rownames(dprobcluster) = lcluster
eff = table(dcluster$cluster)
dprobcluster = cbind(dprobcluster, eff)
colnames(dprobcluster) = c(lbacteria, "Effective")
write.csv(dprobcluster, paste(prout, "SumProb.csv", sep = ""))







