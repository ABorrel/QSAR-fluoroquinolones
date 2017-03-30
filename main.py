from os import listdir

import tableParse
import liganddescriptors
import pathFolder
import runExternalSoft



def CleanCHEMBLFile(pfilin, pfilout):

    # add short cut if filtered table exist !!!!!

    table = tableParse.CHEMBL(pfilin)
    table.parseCHEMBLFile()
    print len(table.table), "Init"

    #table.selectConfidencecore(cutoff=9)
    #print len(table.table), "prot confidence"

    table.getOnlyExactConstant()
    print len(table.table), "strict value"

    table.getOnlyMIC()
    print len(table.table), "MIC"

    table.MICbyOrganisms(pfilout[:-4] + "_MIC-orga.txt")
    table.writeTable(pfilout)

    table.completeMatrixMICByorganism(pfilout[:-4] + "_MIC_", nborga=4)

    return table



def MolecularDesc(dtable, pfilout, plog):

    logfile = open(plog, "w")

    for orga in dtable.keys():
        pfiloutorga = pfilout + "desc_" + orga.replace(" ", "-") + ".csv"

        for dcompound in dtable[orga]:
            desc = liganddescriptors.Descriptors(dcompound, writecheck=0, logfile=logfile)
            desc.get_descriptorOD1D()
            desc.get_descriptor2D()

            if desc.log == "ERROR":
                continue
            desc.writeTablesDesc(pfiloutorga)


def AnalyseDesc(pdesc, pdata, prout, PCA="1", dendo="1", cormatrix="1", hist="1", corcoef=0.0, logaff="0"):

    runExternalSoft.DescAnalysis(pdesc, pdata, prout, corcoef, PCA, cormatrix, hist, dendo, logaff)

    return

##########
#  MAIN  #
##########

pCHEMBL = "/home/aborrel/fluoroquinolones/compound-bioactiveCHEMBL.txt"
pCHEMBLClean = "/home/aborrel/fluoroquinolones/compound_filtered_target.txt"

pCHEMBLCleanMIC = "/home/aborrel/fluoroquinolones/results/compound_filtered_MIC.txt"


#ltab = CleanCHEMBLFile(pCHEMBL, pCHEMBLClean)
#dtab = CleanCHEMBLFile(pCHEMBL, pCHEMBLCleanMIC)

pdesc = pathFolder.PR_DESC
plog = pathFolder.PR_DESC + "log.txt"

#MolecularDesc(dtab.tableorgafull, pdesc, plog)

lfiledesorga = listdir(pathFolder.PR_DESC)
for fileorga in lfiledesorga:
    pfiledesc = pathFolder.PR_DESC + fileorga
    print fileorga[-3:]
    if fileorga[-3:] == "csv":
        orga = fileorga.split("_")[1][:-4]
        AnalyseDesc(pfiledesc, pCHEMBLCleanMIC[:-4] + "_MIC_" + orga + ".csv", pathFolder.PR_DESC + orga, corcoef=0.7,
                    logaff=1)



#lcompound = [ltab[i]["CMPD_CHEMBLID"] for i in range(0,len(ltab))]
#lcompound = list(set(lcompound))
#print len(lcompound), "NB different"

#pdesc = pathFolder.analyses(psub="desc") + "tableDesc.csv"
#plog = pathFolder.analyses(psub="desc") + "log.txt"
#MolecularDesc(ltab, pdesc, plog)





