from os import listdir, remove, makedirs

PR_REF = "/home/aborrel/fluoroquinolones/"
PR_RESULT = "/home/aborrel/fluoroquinolones/results/"
PR_TEMP3D = "/home/aborrel/fluoroquinolones/results/temp3D/"
PR_COMPOUNDS = "/home/aborrel/fluoroquinolones/results/compounds/"
PR_ANALYSIS = "/home/aborrel/fluoroquinolones/results/analysis/"
PR_DESC = "/home/aborrel/fluoroquinolones/results/desc/"



def cleanFolder(prin=PR_TEMP3D):
    lfiles = listdir(prin)
    if len(lfiles) != 0:
        for filin in lfiles:
            # problem with folder
            remove(prin + filin)
    return prin

def analyses(psub):

    if psub == "":
        return PR_ANALYSIS
    else:
        try: makedirs(PR_ANALYSIS + psub + "/")
        except: pass

    return PR_ANALYSIS + psub + "/"

def createFolder(prin):

    try: makedirs(prin)
    except: pass

    return prin

