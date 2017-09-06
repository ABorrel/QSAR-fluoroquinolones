from copy import deepcopy
from numpy import mean, std


def mergeDict(l_dict):

    dout = {}

    lk = l_dict[0].keys()

    for k in lk:
        lval = [l_dict[i][k] for i in range(0,len(l_dict))]
        #identic
        if lval.count(lval[0]) == len(lval):
            dout[k] = lval[0]
        else:
            dout[k] = "----".join(lval)

    return dout



def convertUg(valin, unit):

    if unit == "mg l-1":
        return float(valin)*1000


def convertUnit(l_values, l_units):
    """Convert list of affinity un uM"""

    #print l_values
    #print l_units
    lout = []
    i = 0
    imax = len(l_values)
    while i < imax:
        if l_units[i] == "uM" or l_units[i] == "10'-6M" or l_units[i] == "um" or l_units[i] == "":
            lout.append(l_values[i])
            i += 1
        elif l_units[i] == "nM":
            val = float(l_values[i])
            val = val/1000
            lout.append(val)
            i += 1
        else:
            print l_units[i], "ssss"
            ffff

    #print lout

    return lout

def organiseTableMIC(din, minNumberOrganisms, pfilin):

    dwork = deepcopy(din)
    lM = []
    lorga = []
    for compoundID in dwork.keys():
        if len(dwork[compoundID].keys()) < minNumberOrganisms:
            del dwork[compoundID]
            continue
        else:
            lorga = [i for i in dwork[compoundID].keys()]
            del lorga[lorga.index("SMILES")]
            print lorga
            lMMIC = [mean(dwork[compoundID][orga]) for orga in lorga]
            print lMMIC
            dwork[compoundID]["M"] = mean(lMMIC)
            for orga in dwork[compoundID].keys():
                if not orga in lorga and not orga == "SMILES" and not orga == "M":
                    lorga.append(orga)
            lM.append(mean(lMMIC))

    LM = sorted(lM)
    #print lM

    filout = open(pfilin, "w")
    filout.write("ID\tSMILES\tM\t" + "\t".join(lorga) + "\n")
    for M in LM:
        for compound in dwork.keys():
            if dwork[compound]["M"] == M:
                filout.write(compound + "\t" + str(dwork[compound]["SMILES"]) + "\t" + str(dwork[compound]["M"]))
                for orga in lorga:
                    if orga in dwork[compound].keys():
                        filout.write("\t" + str(mean(dwork[compound][orga])))
                    else:
                        filout.write("\tNA")
                filout.write("\n")
                del dwork[compound]
    filout.close()


def sameMagnitude(lval):
    """Test magnitude"""

    #Mval = mean(lval)
    SDval = std(lval)

    magnitudeval = str("%e"%min(lval))[-1]
    magnitudeSD = str("%e"%SDval)[-1]


    if magnitudeval != magnitudeSD:
        return 0
    else:
        return 1










