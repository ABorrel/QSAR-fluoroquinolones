import parseSDF
import liganddescriptors
import runExternalSoft


psdfin = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/all_carbamates.sdf"
paff = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/Carbamate.csv"





def MolecularDesc(lmolecules, pfilout, plog):

    logfile = open(plog, "w")


    for compound in lmolecules:
        desc = liganddescriptors.Descriptors(compound, writecheck=0, logfile=logfile, kSMILES="SMILES")
        desc.get_descriptorOD1D()
        desc.get_descriptor2D()

        if desc.log == "ERROR":
            continue
        desc.writeTablesDesc(pfilout, kID="Name", kSMILES="SMILES")



#############
#   MAIN    #
#############

pdesc = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/desc/desc.csv"
prdesc = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/desc/"
plog = "/home/aborrel/fluoroquinolones/fluoroquinole_melender/desc/log.txt"

#dsdf = parseSDF.sdf(psdf=psdfin)
#dsdf.parseSDF()

#print dsdf.lc[0].keys()

#MolecularDesc(dsdf.lc, pdesc, plog)



runExternalSoft.DescAnalysis(pdesc, paff, prdesc, valcor=0.9, PCA=1, corMatrix=1, hist=1, dendo=1, logaff=0)

