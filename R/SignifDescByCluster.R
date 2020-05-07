#!/usr/bin/env Rscript
source("tool.R")

################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_desc = args[1]
p_cluster = args[2]
pr_out = args[3]

p_desc = "./../../results/DESC/desc_compound.csv"
p_cluster = "./../../results/Clustering_selected/hclust-ward.D2-gap_stat/Table_hclust_ward.D2_gap_stat.csv"
pr_out = "./../../results/SignifDesc_byCluster/"


# open descriptors #
####################
dglobal = openData(p_desc, 0, pr_out, c(1,2))
dglobal = dglobal[[1]]

rownames(dglobal) = dglobal[,1]
dglobal = dglobal[,-1]# remove name
dglobal = dglobal[,-1]# remove SMILES

# cluster
dcluster = read.csv(p_cluster, sep = ",")
rownames(dcluster) = dcluster[,1]
dcluster = dcluster[rownames(dglobal),]


l_cluster = seq(1, max(dcluster$cluster))
l_desc = colnames(dglobal)

m_out = NULL

for (desc in l_desc){
  l_out = NULL
  for(cluster in l_cluster){
    i_chem = which(dcluster$cluster == cluster)
    v_desc = dglobal[,desc]
    v_cluster = v_desc[i_chem]
    v_nocluster = v_desc[-i_chem]
    
    if(length(i_chem) <= 2){
      l_out = append(l_out, "-")
    }else{
      parametric = conditionTtest(v_cluster, v_nocluster)
      if(parametric == 1){
        pval = comparisonTest (v_cluster, v_nocluster, "parametric")
      }else{
        pval = comparisonTest (v_cluster, v_nocluster, "no-parametric")
      }
      signif = signifPvalue(pval)
      
      l_out = append(l_out, signif)
      #print("====")
      #print(length(v_desc))
      #print(length(v_cluster))
      #print(length(v_nocluster))
    }
  }
  m_out = rbind(m_out, l_out)
}


N_CHEM = NULL
for(cluster in l_cluster){
  i_chem = which(dcluster$cluster == cluster)
  N_CHEM = append(N_CHEM, length(i_chem))  
}

rownames(m_out) = l_desc
colnames(m_out) = l_cluster

m_out = rbind(m_out, N_CHEM)


write.csv(m_out, file=paste(pr_out, "signif_clusters.csv", sep = ""))
