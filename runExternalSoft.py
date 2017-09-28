from os import system, path, remove, listdir
from pymol.gui import cmd
from re import search
from time import sleep

LIGPREP = "/opt/schrodinger2016-4/ligprep"
PADEL = "/home/aborrel/softwares/padel/PaDEL-Descriptor.jar"


def runLigprep(psmilin, forcefield="OPLS3", stereoisoster=1):

    """Maybe fix more option"""

    if forcefield == "OPLS3":
        bff = "16"
    else:
        bff = "14"

    cmd = LIGPREP + " -ismi " + psmilin + " -osd " + psmilin[0:-4] + ".sdf" + " -bff " + str(bff) + " -epik -s " + str(stereoisoster) + " -WAIT -NJOBS 3"

    print cmd
    system(cmd)

    # control if file exist
    if not path.exists(psmilin[0:-4] + ".sdf"):
        return "Ligprep ERROR"
    else:
        try:remove("tem.log")
        except:pass
        try:remove("tem-dropped.smi")
        except: pass
        try:remove("tem-dropped-indices.txt")
        except: pass

    return psmilin[0:-4] + ".sdf"

def runPadel(prin=""):
    """Input include a folder of sdf file"""
    if prin == "":
        return "ERROR - Padel Input"
    else:
        cmd = "java -jar " + PADEL + " -maxruntime 100 -3d -dir " + str(prin) + " -file " + prin + "tem.desc"
        print cmd
        system(cmd)

    return prin + "tem.desc"


def babelConvertSDFtoSMILE(sdfread, clean_smi=0, rm_smi=1):

    tempsdf = open("tempsdf.sdf", "w")
    tempsdf.write(sdfread)
    tempsdf.close()

    psmile = "tempsmile.smi"

    cmd_convert = "babel tempsdf.sdf " + psmile + " 2>/dev/null"
    system(cmd_convert)

    try : filin = open (psmile, "r")
    except : return "0"
    l_Fline = filin.readlines()
    filin.close()
    try : smile = l_Fline[0].split("\t")[0]
    except : return "0"

    # rewrite path in filout
    if clean_smi == 1:
        filout = open (psmile, "w")
        filout.write (str (smile))
        filout.close ()

    if rm_smi == 1:
        system("rm " + psmile)


    return smile



def PCAplot(pfildesc, pfildata, corcoef, prout):

    cmdplotPCA = "./PCAplot.R " + str(pfildesc) + " " + str(pfildata) + " " + str(corcoef) + " " + str(prout)

    print cmdplotPCA
    system(cmdplotPCA)

    return


def DescAnalysis(pdesc, paffinity, prout, valcor, maxquantile, logaff, PCA, corMatrix, hist, dendo, cluster):

    cmdDescAnalysis = "./descAnalysis.R " + str(pdesc) + " " + str(paffinity) + " " + str(prout) + " " + str(valcor) + " " + \
        str(maxquantile) + " " + str(logaff) + " " + str(PCA) + " " + str(corMatrix) + " " + str(hist) + \
        " " + str(dendo) + " " + str(cluster)

    print cmdDescAnalysis
    system(cmdDescAnalysis)


def clusterAnalysis(pdesc, paffinity, pcluster, prout, valcor, maxquantile, logaff=1):

    cmdClusterAnalysis = "./clusterAnalysis.R " + str(pdesc) + " " + str(paffinity) + " " + str(pcluster) + " " + \
        str(prout) + " " + str(valcor) + " " + str(maxquantile) + " " + str(logaff)

    print cmdClusterAnalysis
    system(cmdClusterAnalysis)


def plotMICByCpd(pfilin):

    cmdplot = "./plotMIC.R " + str(pfilin)
    print cmdplot
    system(cmdplot)



def babelConvertSMItoSDF(pfilesmi):

    cmd_convert = "babel " + pfilesmi + " " + pfilesmi[0:-4] + ".sdf 2>/dev/null"
    system(cmd_convert)
    return pfilesmi[0:-4] + ".sdf"


def QSARsReg(ptrain, ptest, pcluster, prout, nbfold=10):

    cmd_QSAR = "./QSARsReg.R " + ptrain + " " + ptest + " " + pcluster + " " + prout + " " + str(nbfold) + " >" + prout + "perfRF.txt"
    print cmd_QSAR
    system(cmd_QSAR)

    return prout + "perf.txt"



def molconvert(pfilin):

    if path.exists(pfilin[:-3] + "jpeg"):
        return pfilin[:-3] + "jpeg"
    cmdconvert = "molconvert \"jpeg:w500,Q95,#ffffff\" " + pfilin + " -o " + pfilin[:-3] + "jpeg"
    system(cmdconvert)



def corAnalysis(paffinity_currated, prcorAnalysis):

    cmd = "./corMIC.R " + str(paffinity_currated) + " " + str(prcorAnalysis)
    print cmd
    system(cmd)


def prepareDataset(pdesc, paff, prout, corcoef, maxQuantile, valSplit, logaff=1):


    # extract train and test file
    dfile = {}
    lfile = listdir(prout)
    for filedir in lfile:
        if search("^train_", filedir):
            bacteria = filedir.split("_")[-1][0:-4]
            if not bacteria in dfile.keys():
                dfile[bacteria] = {}
            dfile[bacteria]["train"] = prout + filedir

        elif search("^test_", filedir):
            bacteria = filedir.split("_")[-1][0:-4]
            if not bacteria in dfile.keys():
                dfile[bacteria] = {}
            dfile[bacteria]["test"] = prout + filedir

    if dfile == {}:
        cmd = "/QSARsPrep.R " + str(pdesc) + " " + str(paff) + " " + prout + " " + str(corcoef) + " " + str(
            maxQuantile) + " " + str(valSplit) + " " + str(logaff)
        print cmd
        system(cmd)
        return prepareDataset(pdesc, paff, prout, corcoef, maxQuantile, valSplit)
    else:
        return dfile





