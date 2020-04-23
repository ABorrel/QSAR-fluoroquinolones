from os import listdir, remove, makedirs

PR_REF = "/home/borrela2/fluoroquinolones/"
PR_RESULT = "/home/borrela2/fluoroquinolones/results/"
PR_TEMP3D = "/home/borrela2/fluoroquinolones/results/temp3D/"
PR_COMPOUNDS = "/home/borrela2/fluoroquinolones/results/compounds/"
PR_ANALYSIS = "/home/borrela2/fluoroquinolones/results/analysis/"
PR_DESC = "/home/borrela2/fluoroquinolones/results/desc/"



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

