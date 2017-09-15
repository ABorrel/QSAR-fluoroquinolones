from os import listdir, path, remove

import tableParse
import liganddescriptors
import pathFolder
import runExternalSoft
import bycluster


def CleanCHEMBLFile(pfilin, presult):

    # add short cut if filtered table exist !!!!!

    table = tableParse.CHEMBL(pfilin)
    table.parseCHEMBLFile()
    print len(table.table), "Init"

    table.getOnlyExactConstant()
    print len(table.table), "strict value"

    table.getOnlyMIC()
    print len(table.table), "MIC"
    table.writeTable(presult + "filtered_MIC.txt")

    table.MICbyOrganisms(presult + "MIC-byorga.txt")

    table.completeMatrixMICByorganism(presult + "MIC_", nborga=4)
    table.compareMICorga(presult)# maybe analyse

    return table



def MolecularDesc(dtable, prdesc, plog, clean = 0):

    logfile = open(plog, "w")

    # useless to compute descriptor for four bacteria, same compound
    pfilout = prdesc + "desc_compound.csv"

    if path.exists(pfilout) and path.getsize(pfilout) > 50:
        return pfilout
    if clean == 1:
        remove(pfilout)

    orga = dtable.keys()[0]# take the first orga

    for dcompound in dtable[orga]:
        desc = liganddescriptors.Descriptors(dcompound, writecheck=0, logfile=logfile)
        desc.get_descriptorOD1D()
        desc.get_descriptor2D()

        if desc.log == "ERROR":
            continue

        desc.writeTablesDesc(pfilout)
    return pfilout


def converToSDF(dtable, prsdf):


    for orga in dtable.keys():
        pfiloutorga = prsdf + orga.replace(" ", "-") + ".sdf"
        filoutorga = open(pfiloutorga, "w")

        for dcompond in dtable[orga]:
            pfilesmi = prsdf + dcompond["CMPD_CHEMBLID"] + ".smi"
            filesmi = open(pfilesmi, "w")
            filesmi.write(dcompond["CANONICAL_SMILES"])
            filesmi.close()

            psdf = runExternalSoft.babelConvertSMItoSDF(pfilesmi)
            fsdf = open(psdf, "r")
            sdf = fsdf.read()
            fsdf.close()

            sdfwrite = str(dcompond["CMPD_CHEMBLID"]) + sdf[0:-5] + "> <name>\n" + str(dcompond["CMPD_CHEMBLID"])  + "\n\n$$$$\n"
            filoutorga.write(sdfwrite)




##########
#  MAIN  #
##########

pCHEMBL = "/home/aborrel/fluoroquinolones/compound-bioactiveCHEMBL.txt"
presult = "/home/aborrel/fluoroquinolones/results/"
paffinity_currated = "/home/aborrel/fluoroquinolones/MIC_currated.csv"

dtab = CleanCHEMBLFile(pCHEMBL, presult)



############ DESC  ###############
##################################
prdesc = pathFolder.createFolder(pathFolder.PR_DESC)
plog = prdesc + "log.txt"


# Descriptors computation and variable #
########################################
pdesc = MolecularDesc(dtab.tableorgafull, prdesc, plog)
corcoef = 0.8
maxQuantile = 85


#######################
# Cross MIC analysis  #
#######################
prcorAnalysis = pathFolder.createFolder(pathFolder.PR_RESULT + "CorrMICAnalysis/")
#runExternalSoft.corAnalysis(paffinity_currated, prcorAnalysis)



#######################
# Analysis descriptor #
# and clustering      #
#######################

pranalysis = pathFolder.createFolder(pathFolder.PR_RESULT + "desc_analysis/" + str(corcoef) + "/")

# desc and visualisation #
##########################
#runExternalSoft.DescAnalysis(pdesc=pdesc, paffinity=paffinity_currated, prout=pranalysis, valcor=corcoef, maxquantile=maxQuantile, logaff=1, PCA=1, corMatrix=1, hist=1, dendo=1, cluster=1) #used to find the different stable cluster


# analysis by cluster and cross #
#################################
pcluster = "/home/aborrel/fluoroquinolones/results/desc_analysis/0.8/Table_hclust_ward.D2_gap_stat.csv"

# cross analysis
prclusteranalysis = pathFolder.createFolder(pathFolder.PR_RESULT + "CrossClusterAnalysis/" + str(pcluster.split("/")[-1][6:-4]) + "/")
#runExternalSoft.clusterAnalysis(pdesc, paffinity_currated, pcluster, prclusteranalysis, corcoef, maxQuantile)

# extract compound by cluster #
###############################
bycluster.extractByCluster(pcluster, prclusteranalysis, pdesc)

#################
#  QSAR models  #
#################
valSplit = 0.15
prQSAR = pathFolder.createFolder(pathFolder.PR_RESULT + "QSARS/")


# split dataset based on descriptors
dfileTrainTest = runExternalSoft.prepareDataset(pdesc, paffinity_currated, prQSAR, corcoef=corcoef, maxQuantile=maxQuantile, valSplit=valSplit)

# run specifically by bacteria
for bacteria in dfileTrainTest.keys():
    print dfileTrainTest[bacteria]
    prQSARorga = pathFolder.createFolder(prQSAR + bacteria + "/")
    runExternalSoft.QSARsReg(dfileTrainTest[bacteria]["train"], dfileTrainTest[bacteria]["test"], pcluster, prQSARorga)


