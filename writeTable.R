#!/usr/bin/env Rscript



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pcluster = args[1]
pMIC = args[2] #to take affinity
prout = args[3]


pcluster = "/home/borrela2/fluoroquinolones/results/desc_analysis/0.8/Table_hclust_ward.D2_gap_stat.csv"
pMIC = "/home/borrela2/fluoroquinolones/results/MIC-full"
prout = "/home/borrela2/fluoroquinolones/results/"


dcluster = read.csv(pcluster, sep = ",", header = TRUE)
rownames(dcluster) = dcluster[,1]
dMIC = read.csv(pMIC, sep = "\t", header = TRUE)
rownames(dMIC) = dMIC[,1]
dMIC = dMIC[,-1]

lcluster = unique(dcluster[,2])


dMcluster = data.frame()
dprobcluster = NULL
imax = length(lcluster)
jmax = length(colnames(dMIC))

print(imax)
print(jmax)

for(i in seq(1,imax)){
  k = 1
  lID = dcluster[which(dcluster[,2] == i),1]
  
  for(j in seq(1, jmax)){
    M = mean(-log10(dMIC[lID,j]))
    SD = sd(-log10(dMIC[lID,j]))
    dMcluster[i, k] = M
    dMcluster[i, k+1] = SD
    k = k + 2
  }  

  dprob = c(0,0,0,0)
  for(chem in lID){
    minMIC = min(dMIC[chem,])
    ibact = which(minMIC == dMIC[chem,])
    dprob[ibact] = dprob[ibact] + 1
  }
  dprob = round(dprob/sum(dprob),2)*100
  dprobcluster = rbind(dprobcluster, dprob)
}



rownames(dMcluster) = lcluster
colnames(dMcluster) = c("M-Pseudomonas.aeruginosa", "SD-Pseudomonas.aeruginosa", "M-Staphylococcus.aureus", "SD-Staphylococcus.aureus", "M-Escherichia.coli", "SD-Escherichia.coli", "M-Streptococcus.pneumoniae", "SD-Streptococcus.pneumoniae")
write.csv(dMcluster, paste(prout, "SumTableM-SD", sep = ""))


rownames(dprobcluster) = lcluster
colnames(dprobcluster) = colnames(dMIC)
write.csv(dprobcluster, paste(prout, "SumProb", sep = ""))
print(dprobcluster)

