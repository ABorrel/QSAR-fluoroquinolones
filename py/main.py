from os import listdir, path, remove
from shutil import copyfile
#from math import log10

import dataset # main functions to build the dataset
import CHEMBLTable
import compoundProcessing
import runExternalSoft

#import tableParse
#import liganddescriptors
import pathFolder
#import runExternalSoft
#import toolbox
#import bycluster
#import parseSDF


################
# define PATH
################

PR_ROOT = "./../../"
PR_DATA = pathFolder.createFolder(PR_ROOT + "data/")
PR_RESULT = pathFolder.createFolder(PR_ROOT + "results/")

P_DATASET = PR_DATA + "compound-bioactiveCHEMBL23.txt"

###########
# 1. Define clean dataset from chembl
pr_dataset = pathFolder.createFolder(PR_RESULT + "dataset/")
p_currated_dataset = dataset.CleanCHEMBLFile(P_DATASET, pr_dataset)

##############
# 2. Load dataset
c_dataset = CHEMBLTable.CHEMBLTable(P_DATASET)
c_dataset.parseCuratedDataset(pr_dataset) # put folder with cleaned dataset

##############
# 3. Process compound

# 3.1. Compute molecular descriptors
pr_desc = pathFolder.createFolder(PR_RESULT + "DESC/")
pdesc = compoundProcessing.MolecularDesc(c_dataset.tableorgafull, pr_desc)

# 3.2. Compute PNG
#compoundProcessing.get_PNGAndSMI(pdesc, PR_RESULT)

################
# 4. Analyse chemicals
corcoef = 0.85 # pairwise correlation by descriptor
maxQuantile = 85 # descriptor distribution in 85% on the same quantile

# 4.1. PCA
pr_PCA = pathFolder.createFolder("%sPCA-%s-%s/"%(PR_RESULT, corcoef, maxQuantile))
#runExternalSoft.DescAnalysis(pdesc=pdesc, paffinity=p_currated_dataset, prout=pr_PCA, valcor=corcoef, maxquantile=maxQuantile, logaff=1, PCA=1, corMatrix=0, hist=0, dendo=0, cluster=0) #used to find the different stable cluster

# 4.2. correlation descriptor
pr_COR_DESC = pathFolder.createFolder("%sCOR_DESC-%s-%s/"%(PR_RESULT, corcoef, maxQuantile))
#runExternalSoft.DescAnalysis(pdesc=pdesc, paffinity=p_currated_dataset, prout=pr_COR_DESC, valcor=corcoef, maxquantile=maxQuantile, logaff=1, PCA=0, corMatrix=1, hist=0, dendo=0, cluster=0) #used to find the different stable cluster

# 4.3. clustering
pr_cluster = pathFolder.createFolder("%sClustering-%s-%s/"%(PR_RESULT, corcoef, maxQuantile))
#runExternalSoft.DescAnalysis(pdesc=pdesc, paffinity=p_currated_dataset, prout=pr_cluster, valcor=corcoef, maxquantile=maxQuantile, logaff=1, PCA=0, corMatrix=0, hist=0, dendo=0, cluster=1)

# 4.4. select manually best clustering
pr_cluster_selected = PR_RESULT + "Clustering_selected/"

################
# 5. analysis pMIC and cluster

# 5.1. Correaltion pMIC by organism
pr_COR_MIC = pathFolder.createFolder("%sCOR_PMI/"%(PR_RESULT))
runExternalSoft.corAnalysis(p_currated_dataset, pr_COR_MIC)

# 5.2. Correlation descriptors by pMIC


# 5.3. most significatif desc by cluster 
p_cluster = ""

# 5.4. Build pdf for pub with chem, cluster and pMIC






###############
# 6. machine learning






ddd










#############
# for final #
#############
pCHEMBL = "/home/borrela2/fluoroquinolones/compound-bioactiveCHEMBL.txt"
pdescfinal = "/home/borrela2/fluoroquinolones/results/QSARS5/globalSet.csv"
pMICmol = "/home/borrela2/fluoroquinolones/MIC_currated_Mol.csv"
pMIC = "/home/borrela2/fluoroquinolones/MIC_currated.csv"
prSI = pathFolder.createFolder("/home/borrela2/fluoroquinolones/SI/")
pdescSMI = "/home/borrela2/fluoroquinolones/results/desc/desc_compound.csv"

formatSI(pdescSMI, pdescfinal, pMICmol, pMIC, prSI)
kkk

##########
#  MAIN  #
##########

#pCHEMBL = "/home/borrela2/fluoroquinolones/compound-bioactiveCHEMBL.txt"
#presult = pathFolder.createFolder("/home/borrela2/fluoroquinolones/results/")

pMIC_molar = CleanCHEMBLFile(pCHEMBL, presult)
ggg

########################
# convert g/l to Molar #
########################

pMIC_molar = "/home/borrela2/fluoroquinolones/MIC_currated_Mol.csv"
toolbox.convertUgLtoMol(paffinity_currated, dtab, pMIC_molar)

############ DESC  ###############
##################################
prdesc = pathFolder.createFolder(pathFolder.PR_DESC)
plog = prdesc + "log.txt"


# Descriptors computation and variable #
########################################
#pdesc = MolecularDesc(dtab.tableorgafull, prdesc, plog)
corcoef = 0.85
maxQuantile = 85


#######################
# Cross MIC analysis  #
#######################
#prcorAnalysis = pathFolder.createFolder(pathFolder.PR_RESULT + "CorrMICAnalysis/")
#runExternalSoft.corAnalysis(pMIC_molar, prcorAnalysis)


#######################
# Analysis descriptor #
# and clustering      #
#######################
pdesc = "/home/borrela2/fluoroquinolones/results/desc/desc_compound.csv"
pranalysis = pathFolder.createFolder(pathFolder.PR_RESULT + "desc_analysis/" + str(corcoef) + "/")

# desc and visualisation #
##########################
runExternalSoft.DescAnalysis(pdesc=pdesc, paffinity=paffinity_currated, prout=pranalysis, valcor=corcoef, maxquantile=maxQuantile, logaff=1, PCA=1, corMatrix=1, hist=1, dendo=1, cluster=1) #used to find the different stable cluster


# analysis by cluster and cross #
#################################
pcluster = "/home/borrela2/fluoroquinolones/results/desc_analysis/0.85/hclust-ward.D2-gap_stat/Table_hclust_ward.D2_gap_stat.csv"


# cross analysis
prclusteranalysis = pathFolder.createFolder(pathFolder.PR_RESULT + "CrossClusterAnalysis/" + str(pcluster.split("/")[-1][6:-4]) + "/")
#runExternalSoft.clusterAnalysis(pdesc, pMIC_molar, pcluster, prclusteranalysis, corcoef, maxQuantile)

# extract compound by cluster #
###############################
#bycluster.extractByCluster(pcluster, prclusteranalysis, pdesc)


# organise compounds #
######################
prrank = pathFolder.createFolder(pathFolder.PR_RESULT + "topChemical/")
rankCompounds(pcluster, prclusteranalysis, pMIC_molar, prrank)

gggg


#################
#  QSAR models  #cetone
#################
# run 5 times

valSplit = 0.15
for i in range(1,6):
    prQSAR = pathFolder.createFolder(pathFolder.PR_RESULT + "QSARS" + str(i) + "/")

    # prepare data
    pdesc = "/home/borrela2/fluoroquinolones/results/desc/desc_compound_withoutSMILES.csv"
    dfileTrainTest = runExternalSoft.prepareDataset(pdesc, pMIC_molar, prQSAR, corcoef=corcoef, maxQuantile=maxQuantile, valSplit=valSplit)

    # run specifically by bacteria
    for bacteria in dfileTrainTest.keys():
        print dfileTrainTest[bacteria]
        prQSARorga = pathFolder.createFolder(prQSAR + bacteria + "/")
        runExternalSoft.QSARsReg(dfileTrainTest[bacteria]["train"], dfileTrainTest[bacteria]["test"], pcluster, prQSARorga)






