from copy import deepcopy
from numpy import mean, std

import runExternalSoft
import toolbox

class CHEMBL:
    def __init__(self, pfilin):
        self.pfilin = pfilin


    def parseCHEMBLFile(self):

        lout = []
        filin = open(self.pfilin, "r")
        llines = filin.readlines()
        filin.close()

        lhead = llines[0].split("\t")[0:-1]
        print lhead

        i = 1
        while i < len(llines):
            lelem = llines[i].split("\t")
            #print lelem, i
            #print len(lelem), len(lhead), "***"
            #print lelem
            dout = {}
            j = 0
            while j < len(lhead):
                dout[lhead[j]] = lelem[j]
                j += 1
            lout.append(dout)
            i += 1

        self.table = lout



    def getOnlyExactConstant(self, MIC=1):
        """Remove inhibition in percentage and not exact value"""

        if not "table" in dir(self):
            self.parseCHEMBLFile()


        i = 0
        imax = len(self.table)
        while i < imax:
            row = self.table[i]
            if row["RELATION"] != "=":
                del self.table[i]
                imax = imax - 1
                continue
            #elif row["PCHEMBL_VALUE"] == "":
            #    del self.table[i]
            #    imax = imax - 1
            #    continue
            else:
                i += 1

        print len(self.table)

    def getOnlyIC50(self):

        if not "table" in dir(self):
            self.parseCHEMBLFile()

        i = 0
        imax = len(self.table)
        while i < imax:
            row = self.table[i]
            if row["PUBLISHED_TYPE"] != "IC50":
                del self.table[i]
                imax = imax - 1
                continue
            else:
                i += 1

    def getOnlyMIC(self):
        """Take only MIC and MIC in uM"""
        if not "table" in dir(self):
            self.parseCHEMBLFile()

        i = 0
        imax = len(self.table)
        while i < imax:
            row = self.table[i]
            if row["STANDARD_TYPE"] != "MIC":
                del self.table[i]
                imax = imax - 1
                continue
            else:
                if row["PUBLISHED_UNITS"] != "ug ml-1" and row["PUBLISHED_UNITS"] != "mg l-1":
                    del self.table[i]
                    imax = imax - 1
                    continue
                i += 1



    def MergeIdenticCHEMBLIDforACtivity(self):
        """Control quality of constant available and compute mean"""

        if not "table" in dir(self):
            self.parseCHEMBLFile()

        d_CHEMBLID = {}
        for row in self.table:
            if not row["CMPD_CHEMBLID"] in d_CHEMBLID.keys():
                d_CHEMBLID[row["CMPD_CHEMBLID"]] = []
            d_CHEMBLID[row["CMPD_CHEMBLID"]].append(deepcopy(row))

        # case no copy
        if len(d_CHEMBLID.keys()) == len(self.table):
            return

        i = 0
        imax = len(d_CHEMBLID.keys())
        while i < imax:
            # case not problem
            #print d_CHEMBLID.keys()[i]
            if len(d_CHEMBLID[d_CHEMBLID.keys()[i]]) == 1:
                i += 1
                continue
            else:
                # Control the published value
                l_PUBLISHED_VALUE = [float(d_CHEMBLID[d_CHEMBLID.keys()[i]][k]["PUBLISHED_VALUE"]) for k in range(0, len(d_CHEMBLID[d_CHEMBLID.keys()[i]]))]
                l_PUBLISHED_UNITS = [d_CHEMBLID[d_CHEMBLID.keys()[i]][k]["PUBLISHED_UNITS"] for k in range(0, len(d_CHEMBLID[d_CHEMBLID.keys()[i]]))]

                if not l_PUBLISHED_UNITS.count(l_PUBLISHED_UNITS[0]) == len(l_PUBLISHED_UNITS):
                    l_PUBLISHED_VALUE = toolbox.convertUnit(l_PUBLISHED_VALUE, l_PUBLISHED_UNITS)
                    #print l_PUBLISHED_UNITS

                MPUBLISHED_VALUE = mean(l_PUBLISHED_VALUE)
                SDPUBLISHED_VALUE = std(l_PUBLISHED_VALUE)

                magnitudeval = len(str(int(min(l_PUBLISHED_VALUE))))
                magnitudeSD = len(str(int(SDPUBLISHED_VALUE)))

                #print magnitudeSD, magnitudeval, "l110", d_CHEMBLID.keys()[i]

                if magnitudeval != magnitudeSD:
                    del d_CHEMBLID[d_CHEMBLID.keys()[i]]
                    imax = imax - 1
                    continue

                else:

                    l_STANDARD_VALUE = [float(d_CHEMBLID[d_CHEMBLID.keys()[i]][k]["STANDARD_VALUE"]) for k in
                                        range(0, len(d_CHEMBLID[d_CHEMBLID.keys()[i]]))]
                    l_PCHEMBL_VALUE = [float(d_CHEMBLID[d_CHEMBLID.keys()[i]][k]["PCHEMBL_VALUE"]) for k in
                                       range(0, len(d_CHEMBLID[d_CHEMBLID.keys()[i]]))]
                    MSTANDARD_VALUE = mean(l_STANDARD_VALUE)


                    d_CHEMBLID[d_CHEMBLID.keys()[i]] = [toolbox.mergeDict(d_CHEMBLID[d_CHEMBLID.keys()[i]])]
                    d_CHEMBLID[d_CHEMBLID.keys()[i]][0]["STANDARD_VALUE"] = str(MSTANDARD_VALUE)
                    d_CHEMBLID[d_CHEMBLID.keys()[i]][0]["PCHEMBL_VALUE"] = str(mean(l_PCHEMBL_VALUE))
                    d_CHEMBLID[d_CHEMBLID.keys()[i]][0]["PUBLISHED_VALUE"] = str(MPUBLISHED_VALUE)
                    i += 1

        # reformate the table
        self.table = []
        for k in d_CHEMBLID.keys():
            self.table.append(d_CHEMBLID[k][0])



    def delIdenticCHEMBLIDByCuration(self):
        """Prefered expert curation or most recent publication"""

        if not "table" in dir(self):
            self.parseCHEMBLFile()

        d_CHEMBLID = {}
        for row in self.table:
            if not row["CMPD_CHEMBLID"] in d_CHEMBLID.keys():
                d_CHEMBLID[row["CMPD_CHEMBLID"]] = []
            d_CHEMBLID[row["CMPD_CHEMBLID"]].append(deepcopy(row))

        # case no copy
        if len(d_CHEMBLID.keys()) == len(self.table):
            return

        for CHEMBLID in d_CHEMBLID.keys():
            if len(d_CHEMBLID[CHEMBLID]) == 1:
                continue
            else:
                l_curation = [d_CHEMBLID[CHEMBLID][i]["CURATED_BY"] for i in range(0,len(d_CHEMBLID[CHEMBLID]))]
                l_years = [d_CHEMBLID[CHEMBLID][i]["YEAR"] for i in range(0,len(d_CHEMBLID[CHEMBLID]))]

                # Curation
                if "Expert" in l_curation:
                    i = 0
                    imax = len(l_curation)
                    while i < len(l_curation):
                        if d_CHEMBLID[CHEMBLID][i]["CURATED_BY"] != "Expert":
                            del self.table[self.table.index(d_CHEMBLID[CHEMBLID][i])]
                            del l_years[i]
                            del l_curation[i]
                            del d_CHEMBLID[CHEMBLID][d_CHEMBLID[CHEMBLID].index(d_CHEMBLID[CHEMBLID][i])]
                            imax = imax -1
                            continue
                        i += 1

                #years
                if len(l_years) == 1:
                    continue
                else:
                    recent = max(l_years)
                    i = 0
                    imax = len(l_years)
                    while i < imax:
                        if d_CHEMBLID[CHEMBLID][i]["YEAR"] != str(recent):
                            del self.table[self.table.index(d_CHEMBLID[CHEMBLID][i])]
                            del l_years[i]
                            del d_CHEMBLID[CHEMBLID][d_CHEMBLID[CHEMBLID].index(d_CHEMBLID[CHEMBLID][i])]
                            imax -= 1
                            continue
                        i += 1

                # case where several identic - take first
                if len(d_CHEMBLID[CHEMBLID]) != 1:
                    for i in range (1,len(d_CHEMBLID[CHEMBLID])) :
                        del self.table[self.table.index(d_CHEMBLID[CHEMBLID][i])]


                #print l_years
                recent = max(l_years)

                #print recent, "eeeee"


    def selectConfidencecore(self, cutoff = 0):

        if not "table" in dir(self):
            self.parseCHEMBLFile()

        i = 0
        imax = len(self.table)
        while i < imax:
            score = float(self.table[i]["CONFIDENCE_SCORE"])
            if score < cutoff:
                del self.table[i]
                imax = imax - 1
                continue
            else:
                i += 1


    def MICbyOrganisms(self, pfilout, verbose=0):

        if not "table" in dir(self):
            self.parseCHEMBLFile()

        dout = {}
        lcompound = [self.table[i]["CMPD_CHEMBLID"] for i in range(0, len(self.table))]
        lcompound = list(set(lcompound))

        for compounID in lcompound:
            if not compounID in dout.keys():
                dout[compounID] = {}


        dorga = {}
        for dcompound in self.table:
            organism = dcompound["ORGANISM"]
            if organism == "":
                continue
            else:
                if dcompound["PUBLISHED_VALUE"] == "" or dcompound["PUBLISHED_TYPE"] != "MIC":
                    continue
                if not organism in dorga.keys():
                    dorga[organism] = 0
                if not organism in dout[dcompound["CMPD_CHEMBLID"]].keys():
                    dout[dcompound["CMPD_CHEMBLID"]][organism] = []
                if not "SMILES" in dout[dcompound["CMPD_CHEMBLID"]]:
                    dout[dcompound["CMPD_CHEMBLID"]]["SMILES"] = dcompound["CANONICAL_SMILES"]
                #print dcompound["PUBLISHED_VALUE"]
                dout[dcompound["CMPD_CHEMBLID"]][organism].append(float(dcompound["PUBLISHED_VALUE"]))
                #print dcompound["STANDARD_UNITS"]

        filout = open(pfilout, "w")
        filout.write("ID\t" + "\t".join(dorga.keys()) + "\n")
        for compoundID in dout.keys():
            filout.write(compoundID)
            for orga in dorga.keys():
                if orga in dout[compoundID].keys():
                    if toolbox.sameMagnitude(dout[compoundID][orga]) == 1:
                        filout.write("\t" + str(mean(dout[compoundID][orga])))
                        dorga[orga] += 1
                    else:
                        filout.write("\tNA")
                        del dout[compoundID][orga]
                else:
                    filout.write("\tNA")
            filout.write("\n")
        filout.close()

        if verbose == 1:
            for orga in dorga.keys():
                print orga, "===", dorga[orga]

        #toolbox.organiseTableMIC(dout, 3, pfilout[:-4] + "best.csv")

        # format orga in first keys
        dformat = {}
        lorga = dorga.keys()
        for orga in lorga:
            if not orga in dformat.keys():
                dformat[orga] = []
            for dcompound in self.table:
                if dcompound["ORGANISM"] == orga:
                    if dcompound["CMPD_CHEMBLID"] in dout.keys():
                        dformat[orga].append(dcompound)
        
        self.tablebyorga = dformat


    def completeMatrixMICByorganism(self, prout, nborga=3):

        if not "tablebyorga" in dir(self):
            self.MICbyOrganisms()


        # keep only best complete matrix
        dorga = {}
        for orga in self.tablebyorga.keys():
            dorga[orga] = len(self.tablebyorga[orga])

        lnumcompound = [dorga[k] for k in dorga.keys()]
        lnumcompound = sorted(lnumcompound, reverse=True)
        #print lnumcompound

        lcompoundID = []
        for orga in self.tablebyorga.keys():
            if len(self.tablebyorga[orga]) < lnumcompound[nborga - 1]:
                del self.tablebyorga[orga]
            else:
                if lcompoundID == []:
                    lcompoundID = [self.tablebyorga[orga][i]["CMPD_CHEMBLID"] for i in range(0,len(self.tablebyorga[orga]))]
                else:
                    lcompoundID = list(
                        set(lcompoundID).intersection([self.tablebyorga[orga][i]["CMPD_CHEMBLID"] for i in range(0,len(self.tablebyorga[orga]))]))

        #print len(lcompoundID)
        #print lcompoundID


        dout = {}
        for orga in self.tablebyorga.keys():
            pfilout = prout + str(orga.replace(" ", "-")) + ".csv"
            filout = open(pfilout, "w")
            lheader = ["CMPD_CHEMBLID", "STANDARD_VALUE"]
            print lheader
            filout.write("\t".join(lheader) + "\n")

            if not orga in dout.keys():
                dout[orga] = []

            for compoundID in lcompoundID:
                ltemp = []
                for compound in self.tablebyorga[orga]:
                    if compound["CMPD_CHEMBLID"] == compoundID:
                        ltemp.append(compound)
                if len(ltemp) > 1:
                    compoundtemp = toolbox.mergeDict(ltemp)
                    #print compoundtemp
                    #print ltemp
                    #print len(ltemp)
                    lvalpublished = []
                    lvalstandard = []
                    for i in range(0,len(ltemp)):
                        lvalpublished.append(float(ltemp[i]["PUBLISHED_VALUE"]))
                        lvalstandard.append(float(ltemp[i]["STANDARD_VALUE"]))
                    compoundtemp["PUBLISHED_VALUE"] = mean(lvalpublished)
                    compoundtemp["STANDARD_VALUE"] = mean(lvalstandard)
                else:
                    compoundtemp = ltemp[0]
                dout[orga].append(compoundtemp)
                filout.write("\t".join(str(compoundtemp[k]) for k in lheader) + "\n")
            filout.close()

        self.tableorgafull = dout



    def compareMICorga(self, prout):

        if not "tableorgafull" in dir(self):
            print "ERROR INPUT - tableorgafull"

        lorga = self.tableorgafull.keys()
        pfilout = prout + "MIC-full.txt"
        filout = open(pfilout, "w")
        filout.write("CMPD_CHEMBLID\t" + "\t".join(lorga) + "\n")

        nbComp = len(self.tableorgafull[lorga[0]])
        i = 0
        while i < nbComp:
            nameCpd = self.tableorgafull[lorga[0]][i]["CMPD_CHEMBLID"]
            filout.write(nameCpd + "\t" + str(self.tableorgafull[lorga[0]][i]["STANDARD_VALUE"]))

            for orga in lorga[1:]:
                for cpd in self.tableorgafull[orga]:
                    if cpd["CMPD_CHEMBLID"] == nameCpd:
                        filout.write("\t" + str(cpd["STANDARD_VALUE"]))
                        break

            filout.write("\n")
            i += 1
        filout.close()
        runExternalSoft.plotMICByCpd(pfilout)



    def writeTable(self, pfilout):

        print len(self.table)
        print self.table[0]

        filout = open(pfilout, "w")
        lheader = self.table[0].keys()
        filout.write("\t".join(lheader) + "\n")
        for row in self.table:
            lw = [row[h] for h in lheader]
            filout.write("\t".join(lw) + "\n")
        filout.close()




