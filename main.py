from numpy.f2py.auxfuncs import l_and

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







def MolecularDesc(ltable, pfilout, plog):

    logfile = open(plog, "w")

    for compound in ltable:
        dcompound = liganddescriptors.Descriptors(compound, logfile)

        if dcompound.log == "ERROR":
            continue

        dcompound.get_descriptorOD1D()
        dcompound.get_descriptor2D()

        dcompound.writeTablesDesc(pfilout)


def Analyse(pdesc, pdata, PCA=1, dendo=1, corcoef=0.0):

    if PCA==1:
        runExternalSoft.PCAplot(pdesc, pdata, corcoef)

    return

##########
#  MAIN  #
##########

pCHEMBL = "/home/aborrel/fluoroquinolones/compound-bioactiveCHEMBL.txt"
pCHEMBLClean = "/home/aborrel/fluoroquinolones/compound_filtered_target.txt"

pCHEMBLCleanMIC = "/home/aborrel/fluoroquinolones/compound_filtered_MIC.txt"

#ltab = CleanCHEMBLFile(pCHEMBL, pCHEMBLClean)
ltab = CleanCHEMBLFile(pCHEMBL, pCHEMBLCleanMIC)

#lcompound = [ltab[i]["CMPD_CHEMBLID"] for i in range(0,len(ltab))]
#lcompound = list(set(lcompound))
#print len(lcompound), "NB different"

#pdesc = pathFolder.analyses(psub="desc") + "tableDesc.csv"
#plog = pathFolder.analyses(psub="desc") + "log.txt"
#MolecularDesc(ltab, pdesc, plog)

