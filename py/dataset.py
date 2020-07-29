from os import path

import CHEMBLTable
import pathFolder

def CleanCHEMBLFile(pfilin, pr_out, nborga=4):
    """
    Main function to prepare the dataset
    - out: file in csv with the dataset
    """
    pr_result = pathFolder.createFolder(pr_out)
    pfilout = pr_result + "MIC-curated_mol.csv"
    if path.exists(pfilout):
        return pfilout    

    # add short cut if filtered table exist !!!!!
    table = CHEMBLTable.CHEMBLTable(pfilin)
    table.parseCHEMBLFile()
    print len(table.table), "Init"

    table.getOnlyExactConstant()
    print len(table.table), "strict value"

    table.getOnlyMIC()
    print len(table.table), "MIC"
    table.writeTable(pr_result + "filtered_MIC_all.csv")
    table.MICbyOrganisms(pr_result + "MIC_byorga_all.csv")

    # Optimize with complete matrix
    print "Optimze complete matrix with", nborga, "organisms"
    table.completeMatrixMICByorganism(nborga)

    table.checkIdenticSMIonFullMatrix()
    print len(table.tableorgafull['Escherichia coli']), "Identic SMI checked"

    # write the dataset by organ
    table.writeByOrganism(pr_result)
    print "Write table by organism"
    
    # convert in M
    table.convertUgmLtoMol(pfilout)

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
