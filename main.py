from os import listdir, path, remove
from shutil import copyfile

import tableParse
import liganddescriptors
import pathFolder
import runExternalSoft
import toolbox
import bycluster
import parseSDF
from math import log10


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

    table.checkIdenticSMIonFullMatrix()
    print len(table.tableorgafull['Escherichia coli']), "Identic SMI check"

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




def rankCompounds(ptableCluster, prcluster, pMIC_molar, prrank):


    drank = {}
    dMIC = toolbox.loadMatrix(pMIC_molar)
    dcluster = toolbox.loadMatrix(ptableCluster, ",")
    lorga = dMIC[dMIC.keys()[0]].keys()
    del lorga[lorga.index("CMPD_CHEMBLID")]

    for orga in lorga:
        drank[orga] = []

        for chem in dMIC.keys():
            drank[orga].append(float(dMIC[chem][orga]))

    for orga in drank.keys():
        drank[orga] = list(sorted(drank[orga], reverse = False))

    for orga in drank.keys():
        prdata = pathFolder.createFolder(prrank + orga + "/")
        lchem = []
        r = 1
        for MIC in drank[orga]:
            for chem in dMIC.keys():
                if not chem in dcluster.keys():
                    continue
                if float(dMIC[chem][orga]) == float(MIC) and not chem in lchem:
                    print dcluster[chem]
                    print prcluster + "cluster" + str(dcluster[chem]["cluster"]) +  "/" + chem + ".jpeg"
                    copyfile(prcluster + "cluster" + str(dcluster[chem]["cluster"]) +  "/" + chem + ".jpeg", prdata + str(r) + "_" + chem + "_" + str(dcluster[chem]["cluster"]) + ".jpeg")
                    lchem.append(chem)
            r = r + 1







def formatSI(pSMI, pdescfinal, pMICmol, pMIC, prout):


    ddescSMI = toolbox.loadMatrix(pSMI, sep = "\t")

    cdesc = toolbox.loadMatrix(pdescfinal, sep=",")
    cMICmol = toolbox.loadMatrix(pMICmol, sep="\t")
    cMIC = toolbox.loadMatrix(pMIC, sep="\t")

    print cMIC
    print cMICmol["CHEMBL192226"]
    #ddd

    psdfAll = prout + "all.sdf"
    fsdfAll = open(psdfAll, "w")
    prsdf = pathFolder.createFolder(prout + "SDF/")
    for chemID in cdesc.keys():
        print chemID
        pSMI = prsdf + chemID + ".smi"
        fsmi = open(pSMI, "w")
        for chem in ddescSMI.keys():
            if chem == chemID:
                fsmi.write(ddescSMI[chem]["SMILES"])
                break
        fsmi.close()
        psdf = runExternalSoft.babelConvertSMItoSDF(pSMI, H=0)
        fsdf = open(psdf, "r")
        rsdf = fsdf.readlines()
        fsdf.close()


        fsdfAll.write("%s\n%s\n%s"
                      ">  <MIC (mg) Escherichia coli>\n%s\n\n>  <MIC (mg) Pseudomonas aeruginosa>\n%s\n\n"
                      ">  <MIC (mg) Staphylococcus aureus>\n%s\n\n>  <MIC (mg) Streptococcus pneumoniae>\n%s\n\n"
                      ">  <MIC (M) Escherichia coli>\n%s\n\n>  <MIC (M) Pseudomonas aeruginosa>\n%s\n\n"
                      ">  <MIC (M) Staphylococcus aureus>\n%s\n\n>  <MIC (M) Streptococcus pneumoniae>\n%s\n\n"
                      ">  <pMIC Escherichia coli>\n%.2f\n\n>  <pMIC Pseudomonas aeruginosa>\n%.2f\n\n"
                      ">  <pMIC Staphylococcus aureus>\n%.2f\n\n>  <pMIC Streptococcus pneumoniae>\n%.2f\n\n"
                      "$$$$\n"%(chemID,chemID,"".join(rsdf[2:-1]),
                                      cMIC[chemID]["Escherichia coli"], cMIC[chemID]["Pseudomonas aeruginosa"],
                                      cMIC[chemID]["Staphylococcus aureus"], cMIC[chemID]["Streptococcus pneumoniae"],
                                      cMICmol[chemID]["Escherichia coli"], cMICmol[chemID]["Pseudomonas aeruginosa"],
                                      cMICmol[chemID]["Staphylococcus aureus"], cMICmol[chemID]["Streptococcus pneumoniae"],
                                      log10(float(cMICmol[chemID]["Escherichia coli"])), log10(float(cMICmol[chemID]["Pseudomonas aeruginosa"])),
                                      log10(float(cMICmol[chemID]["Staphylococcus aureus"])), log10(float(cMICmol[chemID]["Streptococcus pneumoniae"]))))


    fsdfAll.close()


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
paffinity_currated = presult + "MIC-full"

#dtab = CleanCHEMBLFile(pCHEMBL, presult)
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






