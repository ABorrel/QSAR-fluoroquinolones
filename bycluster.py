import toolbox
import runExternalSoft
import pathFolder


def extractByCluster(pcluster, prout, pdesc):

    dcluster = toolbox.formatClusterTable(pcluster)
    lcpddesc = toolbox.tableTolist(pdesc)


    for cluster in dcluster.keys():
        prcluster = pathFolder.createFolder(prout + "cluster" + str(cluster) + "/")

        for compound in dcluster[cluster]:
            for cpddec in lcpddesc:
                if compound == cpddec["ID"]:
                    smile = cpddec["SMILES"]
                    psmi = prcluster + str(compound) + ".smi"
                    filesmi = open(psmi,"w")
                    filesmi.write(smile)
                    filesmi.close()
                    runExternalSoft.molconvert(psmi)
                    break
