#!/usr/bin/env Rscript
source("tool.R")

################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_desc = args[1]
p_aff = args[2]
pr_out = args[3]

p_desc = "./../../results/DESC/desc_compound.csv"
p_aff = "./../../results/dataset/MIC-curated_mol.csv"
pr_out = "./../../results/COR_DESCvsPMIC/"



# open affinity file #
######################

daffinity = read.csv(p_aff, sep = "\t", header = TRUE)
rownames(daffinity) = daffinity[,1]
daffinity = daffinity[,-1]
daffinity = daffinity[,-which(colnames(daffinity) == "SMILES")] # remove SMILES 
daffinity = -log10(daffinity)

# open descriptors #
####################
dglobal = openData(p_desc, 0, pr_out, c(1,2))
dglobal = dglobal[[1]]

rownames(dglobal) = dglobal[,1]
dglobal = dglobal[,-1]# remove name
dglobal = dglobal[,-1]# remove SMILES

# order aff
daffinity = daffinity[intersect(rownames(daffinity), rownames(dglobal)),]

# compute cor #
###############
m_out = NULL
for(col_name in colnames(daffinity)){
  m_temp = cbind(dglobal, daffinity[,col_name])
  m_cor = cor(m_temp)
  l_cor = m_cor[, dim(m_cor)[2]]
  m_out = cbind(m_out, l_cor)
}
m_out = m_out[-dim(m_out)[1],]

colnames(m_out) = colnames(daffinity)
m_out = round(m_out, 2)

write.csv(m_out, file=paste(pr_out, "corDescVSpMIC.csv", sep = ""))

