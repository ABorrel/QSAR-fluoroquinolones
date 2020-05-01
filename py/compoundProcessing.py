from os import path

import Descriptors
import toolbox
import pathFolder
import runExternalSoft

def MolecularDesc(dtable, pr_desc, clean = 0):

    # useless to compute descriptor for four bacteria, same compound
    pfilout = pr_desc + "desc_compound.csv"

    if path.exists(pfilout) and path.getsize(pfilout) > 50:
        return pfilout
    if clean == 1:
        remove(pfilout)

    # log file
    p_log = pr_desc + "log.txt"
    logfile = open(p_log, "w")

    l_chemID = dtable[dtable.keys()[0]].keys()

    for chemID in l_chemID:
        desc = Descriptors.Descriptors(dtable[dtable.keys()[0]][chemID], writecheck=0, logfile=logfile)
        desc.get_descriptorOD1D()
        desc.get_descriptor2D()

        if desc.log == "ERROR":
            continue

        desc.writeTablesDesc(pfilout)
    return pfilout


def get_PNGAndSMI(p_desc, pr_results):

    # load the ddesc to have SMILES cleanned 
    ddesc = toolbox.loadMatrix(p_desc)
    pr_smi = pathFolder.createFolder(pr_results + "SMI/")
    pr_png = pathFolder.createFolder(pr_results + "PNG/")

    for chemID in ddesc.keys():
        print chemID
        SMILES = ddesc[chemID]["SMILES"]
        p_fsmi = pr_smi + chemID + ".smi"
        fsmi = open(p_fsmi, "w")
        fsmi.write(SMILES)
        fsmi.close()

        p_fpng = pr_png + chemID + ".png"
        runExternalSoft.molconvert(p_fsmi, p_fpng)

    return




