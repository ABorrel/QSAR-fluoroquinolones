#!/usr/bin/env Rscript
library(ggplot2)

################
#     MAIN     #
################

args <- commandArgs(TRUE)
paff = args[1]
prout = args[2]


paff = "/home/aborrel/fluoroquinolones/MIC_currated.csv"
prout = "/home/aborrel/fluoroquinolones/results/CorrMICAnalysis/"


# open affinity file #
######################

daffinity = read.csv(paff, sep = ",", header = TRUE)
rownames(daffinity) = daffinity[,1]
daffinity = daffinity[,-1]
daffinity = -log10(daffinity)
print(dim(daffinity))

print(colnames(daffinity))

corTable = data.frame()
pvalTable = data.frame()

i = 1
while(i <= dim(daffinity)[2]){
  j = i 
  while(j <= dim(daffinity)[2]){
    cor.val = cor(daffinity[,i], daffinity[,j])
    p.val = cor.test (daffinity[,i], daffinity[,j])$p.value
    print (p.val)
    corTable[i,j] =  cor.val
    pvalTable[i,j] = p.val
    j = j + 1
  }
  i = i + 1
}

rownames(corTable) = colnames(daffinity)
colnames(corTable) = colnames(daffinity)
write.csv(corTable, paste(prout, "cortable.csv", sep = ""))

rownames(pvalTable) = colnames(daffinity)
colnames(pvalTable) = colnames(daffinity)
write.csv(pvalTable, paste(prout, "pvaltable.csv", sep = ""))


p = ggplot(daffinity, aes(Pseudomonas.aeruginosa, Staphylococcus.aureus))+
  geom_point(size=1.5, col="black", shape=21) + 
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
  labs(x = expression(paste("pMIC ", italic("P. aeruginosa"), sep = "")), y =expression( paste("pMIC ", italic("S. aureus"), sep = ""))) + 
  xlim (c(-2.5, 2.5)) +
  geom_segment(aes(x = -2.5, y = -2.5, xend = 2.5, yend = 2.5)) + 
  ylim (c(-2.5, 2.5)) 
#print(p)
ggsave(paste(prout, "Pseudomonas.aeruginosa_Staphylococcus.aureus.png",sep=""), width = 6,height = 6, dpi = 300)



p = ggplot(daffinity, aes(Pseudomonas.aeruginosa, Escherichia.coli))+
  geom_point(size=1.5, col="black", shape=21) +
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
  labs(x = expression(paste("pMIC ", italic("P. aeruginosa"), sep = "")), y = expression(paste("pMIC ", italic("E. coli"), sep = ""))) + 
  xlim (c(-2.5, 2.5))+
  geom_segment(aes(x = -2.5, y = -2.5, xend = 2.5, yend = 2.5))+ 
  ylim (c(-2.5, 2.5)) 
#print(p)
ggsave(paste(prout, "Pseudomonas.aeruginosa_Escherichia.coli.png",sep=""), width = 6,height = 6, dpi = 300)



p = ggplot(daffinity, aes(Pseudomonas.aeruginosa, Streptococcus.pneumoniae))+
  geom_point(size=1.5, col="black", shape=21) +
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
  labs(x = expression(paste("pMIC ", italic("P. aeruginosa"), sep = "")), y = expression(paste("pMIC ", italic("S. pneumoniae"), sep = ""))) + 
  xlim (c(-2.5, 2.5)) +
  geom_segment(aes(x = -2.5, y = -2.5, xend = 2.5, yend = 2.5)) + 
  ylim (c(-2.5, 2.5)) 
#print(p)
ggsave(paste(prout, "Pseudomonas.aeruginosa_Streptococcus.pneumoniae.png",sep=""), width = 6,height = 6, dpi = 300)


p = ggplot(daffinity, aes(Escherichia.coli, Staphylococcus.aureus))+
  geom_point(size=1.5, col="black", shape=21) +
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
  labs(x = expression(paste("pMIC ", italic("E. coli"), sep = "")), y = expression(paste("pMIC ", italic("S. aureus"), sep = ""))) + 
  xlim (c(-2.5, 2.5)) +
  geom_segment(aes(x = -2.5, y = -2.5, xend = 2.5, yend = 2.5)) + 
  ylim (c(-2.5, 2.5)) 
#print(p)
ggsave(paste(prout, "Escherichia.coli_Staphylococcus.aureus.png",sep=""), width = 6,height = 6, dpi = 300)


p = ggplot(daffinity, aes(Streptococcus.pneumoniae, Staphylococcus.aureus))+
  geom_point(size=1.5, col="black", shape=21) +
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
  labs(x = expression(paste("pMIC ", italic("S. pneumoniae"), sep="")), y = expression(paste("pMIC ", italic("S. aureus")))) + 
  xlim (c(-2.5, 2.5)) +
  geom_segment(aes(x = -2.5, y = -2.5, xend = 2.5, yend = 2.5)) + 
  ylim (c(-2.5, 2.5)) 
#print(p)
ggsave(paste(prout, "Streptococcus.pneumoniae_Staphylococcus.aureus.png",sep=""), width = 6,height = 6, dpi = 300)


p = ggplot(daffinity, aes(Escherichia.coli, Streptococcus.pneumoniae))+
  geom_point(size=1.5, col="black", shape=21) +
  theme(axis.text.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 25, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 25, hjust = 0.5, vjust =0.1))+
  labs(x = expression(paste("pMIC ", italic("E. coli"), sep = "")), y = expression(paste("pMIC ", italic("S. pneumoniae"), sep = ""))) + 
  xlim (c(-2.5, 2.5)) +
  geom_segment(aes(x = -2.5, y = -2.5, xend = 2.5, yend = 2.5)) + 
  ylim (c(-2.5, 2.5)) 
#print(p)
ggsave(paste(prout, "Escherichia.coli_Streptococcus.pneumoniae.png",sep=""), width = 6,height = 6, dpi = 300)








