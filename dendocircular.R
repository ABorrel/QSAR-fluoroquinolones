#!/usr/bin/env Rscript
library(ape)
library(phangorn)
library(ggtree)
library(factoextra)




MultipleDendogramCircle = function(ddes, daff, prout){
  
  #calibrate affinity for color
  minMatrix = min(daff)
  maxMatrix = max(daff)
  
  for(i in seq(1, dim(daff)[2])){
    daff[which(daff[,i] == max(daff[,i])), i] = maxMatrix
    daff[which(daff[,i] == min(daff[,i])), i] = minMatrix
  }
  
  daff = as.data.frame(daff)
  daff = cbind(rownames(daff), daff)
  
 
  matTrans1 <- scale(ddes)
  d <- dist(matTrans1, method = "euc")
  tupgma2 <- upgma(d, method="ward.D2")
  legend.title = "MIC"
  
  pfilout = paste(prout, "Pseudomonas-aeruginosa_dendo.png", sep = "")
  t4 <- ggtree(tupgma2, layout="circular", size=1)
  t4 <- t4 %<+% daff + geom_text(aes(color=Pseudomonas.aeruginosa, label=label, angle=angle, fontface="bold"), hjust=-0.15, size=1.3) +
    geom_tippoint(aes(color=Pseudomonas.aeruginosa), alpha=0.75, size=2)+
    scale_color_continuous(low='red', high='lightgreen') +
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_treescale(x = 20, y = 20, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  ggsave(pfilout, dpi=300, height = 6, width = 7)
  
    
  pfilout = paste(prout, "Escherichia-coli_dendo.png", sep = "")
  t4 <- ggtree(tupgma2, layout="circular", size=1)
  t4 <- t4 %<+% daff + geom_text(aes(color=Escherichia.coli, label=label, angle=angle, fontface="bold"), hjust=-0.15, size=1.3) +
    geom_tippoint(aes(color=Escherichia.coli), alpha=0.75, size=2)+
    scale_color_continuous(low='red', high='lightgreen') +
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_treescale(x = 20, y = 20, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  ggsave(pfilout, dpi=300, height = 6, width = 7)
  
  
  pfilout = paste(prout, "Streptococcus-pneumoniae_dendo.png", sep = "")
  t4 <- ggtree(tupgma2, layout="circular", size=1)
  t4 <- t4 %<+% daff + geom_text(aes(color=Streptococcus.pneumoniae, label=label, angle=angle, fontface="bold"), hjust=-0.15, size=1.3) +
    geom_tippoint(aes(color=Streptococcus.pneumoniae), alpha=0.75, size=2)+
    scale_color_continuous(low='red', high='lightgreen') +
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_treescale(x = 20, y = 20, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  ggsave(pfilout, dpi=300, height = 6, width = 7)
  
  
  pfilout = paste(prout, "Staphylococcus.aureus_dendo.png", sep = "")
  t4 <- ggtree(tupgma2, layout="circular", size=1)
  t4 <- t4 %<+% daff + geom_text(aes(color=Staphylococcus.aureus, label=label, angle=angle, fontface="bold"), hjust=-0.15, size=1.3) +
    geom_tippoint(aes(color=Staphylococcus.aureus), alpha=0.75, size=2)+
    scale_color_continuous(low='red', high='lightgreen') +
    theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm")) +
    geom_treescale(x = 20, y = 20, width = NULL, offset = NULL,
                   color = "white", linesize = 1E-100, fontsize = 1E-100)
  ggsave(pfilout, dpi=300, height = 6, width = 7)
}





dendogramCluster = function(ddes, daff, dcluster, prout){
  
  
  #calibrate affinity for color
  minMatrix = min(daff)
  maxMatrix = max(daff)
  
  for(i in seq(1, dim(daff)[2])){
    daff[which(daff[,i] == max(daff[,i])), i] = maxMatrix
    daff[which(daff[,i] == min(daff[,i])), i] = minMatrix
  }
  
  daff = as.data.frame(daff)
  daff = cbind(rownames(daff), daff)
  dcluster[,2] = as.character(dcluster[,2])
  
  
  matTrans1 <- scale(ddes)
  d <- dist(matTrans1, method = "euc")
  tupgma2 <- upgma(d, method="ward.D2")
  tupgma2 <- groupOTU(tupgma2, 35)

  pfilout = paste(prout, "dendo_cluster.png", sep = "")
    
  t1 <- ggtree(tupgma2, layout="circular", size=0.8)
  t1 <- t1 %<+% dcluster + geom_text(aes(color=cluster, label=cluster, angle=angle,  fontface="bold"), hjust=-1.5, size=1.2) +
    geom_tippoint(aes(color=cluster), alpha=0.75, size=1)+
    #theme(legend.position="right")+
    theme(plot.margin = unit(c(0,0,0,0), "cm"))
    #geom_treescale(x = 40, y = 40, width = NULL, offset = NULL,
    #               color = "white", linesize = 1E-100, fontsize = 1E-100)
  
  
  daff = daff[,-1]
  t2 = gheatmap(t1, daff, font.size = 0, offset = 3, width = 0.2, colnames_offset_x = 2, colnames_offset_y = -0.5, low = "red", high = "lightgreen") +
    #scale_color_continuous(low='red', high='lightgreen') +
    theme_tree()
  #print (t2)
  open_tree(t2, 15) %>% rotate_tree(15)
  ggsave(pfilout, dpi=300, height = 11, width = 11)
}




