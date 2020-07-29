import toolbox
from math import log10


class activityCliff:
    def __init__(self, p_MIC, p_cluster, pr_out):
        self.p_MIC = p_MIC
        self.p_cluster = p_cluster
        self.pr_out = pr_out
    
    def loadData(self):
        self.d_MIC = toolbox.loadMatrix(self.p_MIC)
        self.d_cluster = toolbox.loadMatrix(self.p_cluster, sep = ",")

    def formatForSALI(self):
        # define list of bacteria
        l_bacteria = list(self.d_MIC[list(self.d_MIC.keys())[0]].keys())
        l_bacteria.remove("SMILES")
        l_bacteria.remove("CMPD_CHEMBLID")
        
        dout = {}
        # by cluster and chemical
        for chem in self.d_cluster.keys():
            print(self.d_cluster[chem])
            cluster = self.d_cluster[chem]["cluster"]
            if not cluster in list(dout.keys()):
                dout[cluster] = []
            dout[cluster].append(chem)
        
        p_filout = self.pr_out + "Activity-Cliff.txt"
        filout = open(p_filout, "w")
        filout.write("SMILES\tID\tcluster\t%s\n"%("\t".join("pMIC " + bact for bact in l_bacteria)))
        for cluster in dout.keys():
            for chem in dout[cluster]:
                lpMIC = []
                for bacteria in l_bacteria:
                    lpMIC.append(str(-log10(float(self.d_MIC[chem][bacteria]))))
                print("\t".join(lpMIC))
                filout.write("%s\t%s\t%s\t%s\n"%(self.d_MIC[chem]["SMILES"], self.d_MIC[chem]["CMPD_CHEMBLID"], cluster, "\t".join(lpMIC)))
        filout.close()


    def filterSALIFile(self, cutoff_deltaAct, pSALI):
        l_SALI = toolbox.tableTolist(pSALI)

        pfilout = pSALI[:-4] + "_filtered.txt"
        filout = open(pfilout, "w")
        filout.write("ID1\tID2\tDelta activity\tSALI\tCluster\n")
        dout = {}
        for d_SALI in l_SALI:
            cluster1 = int(self.d_cluster[d_SALI["ID 1"]]["cluster"])
            cluster2 = int(self.d_cluster[d_SALI["ID 2"]]["cluster"])

            if cluster1 == cluster2:
                if not cluster1 in dout.keys():
                    dout[cluster1] = []
                dout[cluster1].append(d_SALI)

        for clust in dout.keys():
            for d_SALI in dout[clust]:
                if float(d_SALI["Delta Activity"]) >=cutoff_deltaAct and float(d_SALI["SALI"]) != 0.0:
                    filout.write("%s\t%s\t%s\t%s\t%s\n"%(d_SALI["ID 1"], d_SALI["ID 2"], d_SALI["Delta Activity"], d_SALI["SALI"], clust))

        filout.close()        









