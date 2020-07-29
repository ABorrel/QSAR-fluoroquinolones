#!/usr/bin/env Rscript
source("tool.R")
library(ggplot2)

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

# plot correlation #
####################

for(col_name in colnames(daffinity)){
  for(desc in colnames(dglobal)){
    dplot = cbind(daffinity[,col_name], dglobal[,desc])
    colnames(dplot) = c("pMIC", "Desc")
    dplot = as.data.frame(dplot)
    
    valr2 = cor(dplot[,1], dplot[,2])

    p = ggplot(dplot, aes(pMIC, Desc))+
      geom_point(size=1.5, colour="black", shape=21) + 
      geom_text(x=7.5, y=min (dplot$Desc) + ((max(dplot$Desc) - min (dplot$Desc))*0.90), label = paste("R2=",round(valr2,2), sep = ""), size = 6)+
      labs(x = "pMIC (mol/l)", y = desc)  +
      theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
      xlim (c(3, 8))
    ggsave(paste(pr_out, col_name, "_", desc, ".png", sep = ""), width = 6,height = 6, dpi = 300)    
  }
}




