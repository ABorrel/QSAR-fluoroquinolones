from os import system, path, remove, listdir, getcwd, chdir
from re import search
from time import sleep

LIGPREP = "/opt/schrodinger2016-4/ligprep"
PADEL = "/home/aborrel/softwares/padel/PaDEL-Descriptor.jar"
PRSCRIPT = "./../R/"


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


def babelConvertSMItoSDF(pfilesmi, H=0):

    if H ==1:
        cmd_convert = "babel " + pfilesmi + " " + pfilesmi[0:-4] + ".sdf -h --gen2d 2>/dev/null"
    else:
        cmd_convert = "babel " + pfilesmi + " " + pfilesmi[0:-4] + ".sdf --gen2d 2>/dev/null"

    system(cmd_convert)
    return pfilesmi[0:-4] + ".sdf"



def molconvert(p_smi, p_png):

    if path.exists(p_png):
        return p_png
    cmdconvert = "molconvert \"jpeg:w500,Q95,#ffffff\" " + p_smi + " -o " + p_png
    print cmdconvert
    system(cmdconvert)


def pngtopdf(ppng):
    cmd = "convert -quality 100 -density 50 " + ppng + " " + ppng[:-3] + "pdf"
    print cmd
    system(cmd)
    return ppng[:-3] + "pdf"

def mergepdfs(lpdfs, pout):

    cmd = "pdfunite " + " ".join(lpdfs) + " " + pout
    print cmd
    system(cmd)


#################
# R scripts run #
#################

def runRscript(cmd, out=0):

    chdir("./../R/")
    print cmd
    if out == 0:
        #system(cmd)  ####### to run
        output = 0
    else:
        import subprocess
        #output = subprocess.check_output(cmd, shell=True)  ####### to run
    chdir("./../py/")
    return output


def corBypMICAnalysis(paffinity_currated, prcorAnalysis):
    cmd = "./corBetweenMIC.R " + str(paffinity_currated) + " " + str(prcorAnalysis)
    runRscript(cmd)


def plotMICByCpd(pfilin, pfilout):
    cmdplot = "./plotMIC.R " + str(pfilin) + " " + pfilout
    print cmdplot
    runRscript(cmdplot)


def clusterAnalysis(pdesc, paffinity, pcluster, prout, valcor, maxquantile, logaff=1):
    cmdClusterAnalysis = "./clusterAnalysis.R " + str(pdesc) + " " + str(paffinity) + " " + str(pcluster) + " " + \
        str(prout) + " " + str(valcor) + " " + str(maxquantile) + " " + str(logaff)

    runRscript(cmdClusterAnalysis)


def PCAplot(pfildesc, pfildata, corcoef, prout):
    cmdplotPCA = "./PCAplot.R " + str(pfildesc) + " " + str(pfildata) + " " + str(corcoef) + " " + str(prout)
    runRscript(cmdplotPCA)


def DescAnalysis(pdesc, paffinity, prout, valcor, maxquantile, logaff, PCA, corMatrix, hist, dendo, cluster):
    cmdDescAnalysis = "./descAnalysis.R " + str(pdesc) + " " + str(paffinity) + " " + str(prout) + " " + str(valcor) + " " + \
        str(maxquantile) + " " + str(logaff) + " " + str(PCA) + " " + str(corMatrix) + " " + str(hist) + \
        " " + str(dendo) + " " + str(cluster)

    runRscript(cmdDescAnalysis)


def corDescVSpMIC(p_desc, p_currated_dataset, pr_out): 
    cmd = "./corDESCvsPMIC.R %s %s %s"%(p_desc, p_currated_dataset, pr_out)
    runRscript(cmd)


def signifDescByCluster(p_desc, p_cluster, pr_out):
    cmd = "./SignifDescByCluster.R %s %s %s"%(p_desc, p_cluster, pr_out)
    runRscript(cmd)


def sumByCluster(p_cluster, p_dataset, pr_out):
    cmd = "./summaryByCluster.R %s %s %s"%(p_cluster, p_dataset, pr_out)
    runRscript(cmd)

##################
# QSAR modeling  #
##################

def runRQSARModeling(cmd):

    workdir = getcwd()
    chdir("/home/borrela2/development/QSARPR/source/")
    print(cmd)
    system(cmd)
    chdir(workdir)



def prepareDataset(pdesc, paff, prout, corcoef, maxQuantile, valSplit, typeAff="All", logaff=1, nbNA = 10):

    # extract train and test file
    dfile = {}
    lfile = listdir(prout)
    for filedir in lfile:
        if search("trainSet", filedir):
            bacteria = filedir.split("_")[0]
            if not bacteria in dfile.keys():
                dfile[bacteria] = {}
            dfile[bacteria]["train"] = prout + filedir

        elif search("testSet", filedir):
            bacteria = filedir.split("_")[0]
            if not bacteria in dfile.keys():
                dfile[bacteria] = {}
            dfile[bacteria]["test"] = prout + filedir

    if dfile == {}:
        cmd = "./QSARsPrep.R " + str(pdesc) + " " + str(paff) + " " + prout + " " + str(corcoef) + " " + str(
            maxQuantile) + " " + str(valSplit) + " " + str(logaff) + " " + str(typeAff) + " " + str(nbNA)
        runRQSARModeling(cmd)
        return prepareDataset(pdesc, paff, prout, corcoef, maxQuantile, valSplit, typeAff, logaff, nbNA)
    else:
        return dfile



def QSARsReg(ptrain, ptest, pcluster, prout, nbfold=10):

    cmd_QSAR = "./QSARsReg.R " + ptrain + " " + ptest + " " + pcluster + " " + prout + " " + str(nbfold) + " 1 >" + prout + "perf.txt"
    runRQSARModeling(cmd_QSAR)

    return prout + "perf.txt"