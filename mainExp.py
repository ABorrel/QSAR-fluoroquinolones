import parseSDF
import liganddescriptors
import runExternalSoft
import tableParse

psdfin = "/home/aborrel/lig8-15/all.sdf"
#psdfin = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/all_carbamates.sdf"
#paff = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/Carbamate.csv"




def MolecularDesc(lmolecules, pfilout, plog):

    logfile = open(plog, "w")


    lsmile = []
    for compound in lmolecules:
        desc = liganddescriptors.Descriptors(compound, writecheck=0, logfile=logfile, kSMILES="SMILES", kID="ID")
        print desc.compound["SMILES"]
        if not desc.compound["SMILES"] in lsmile:
            desc.get_descriptorOD1D()
            desc.get_descriptor2D()
            lsmile.append(desc.compound["SMILES"])
        else:
            print "ddd"
            continue

        if desc.log == "ERROR":
            continue
        desc.writeTablesDesc(pfilout, kID="ID", kSMILES="SMILES")



#############
#   MAIN    #
#############

#pdesc = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/desc/desc.csv"
#prdesc = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/desc/"
#plog = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/desc/log.txt"


pdesc = "/home/aborrel/lig8-15/desc.csv"
prdesc = "/home/aborrel/lig8-15/desc/"
plog = "/home/aborrel/lig8-15/desc/log.txt"


dsdf = parseSDF.sdf(psdf=psdfin)
dsdf.parseSDF()
dsdf.removeDuplicate("ID")

MolecularDesc(dsdf.lc, pdesc, plog)





# for ZINC compound #
#####################
#ptableZinc = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/zinc_compo.txt"
#pdesctest = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/desc/zinc_desc.csv"

#dcompoundtest = tableParse.CHEMBL(ptableZinc)
#dcompoundtest.parseCHEMBLFile()

#print dcompoundtest.table

#MolecularDesc(dcompoundtest.table, pdesctest, plog)


#runExternalSoft.DescAnalysis(pdesc, paff, prdesc, valcor=0.9, PCA=1, corMatrix=1, hist=1, dendo=1, logaff=0)

