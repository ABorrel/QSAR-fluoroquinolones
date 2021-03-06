from re import search




class sdf:

    def __init__(self, psdf):
        self.psdf = psdf


    def parseSDF(self):

        lout = []

        filin = open(self.psdf, "r")
        handle_read = filin.read()
        print handle_read

        l_compound = handle_read.split("$$$$\n")[:-1]
        print l_compound
        ddd

        for compound in l_compound:
            dcompound = {}
            dcompound["sdf"] = compound + "$$$$"
            llines = compound.split("\n")
            i = 0
            nblines = len(llines)
            while i < nblines:
                if search("> <", llines[i]):
                    kin = llines[i].split("> <")[1]
                    kin = kin.split(">")[0]
                    valuek = llines[i + 1].strip()
                    #valuek = valuek.replace("\n", "")
                    dcompound[kin] = valuek
                elif search(">  <", llines[i]):
                    kin = llines[i].split(">  <")[1]
                    kin = kin.split(">")[0]
                    valuek = llines[i + 1].strip()
                    #valuek = valuek.replace("\n", "")
                    dcompound[kin] = valuek

                i += 1
            lout.append(dcompound)
        self.lc = lout
        print len(lout)


    def get_dockingscore(self):
        """Keep smaller score"""

        if not "lc" in dir(self):
            self.parseSDF()

        if "docking" in dir(self):
            return self.docking

        dscore = {}
        for compound in self.lc:
            print compound.keys()
            print compound["s_m_entry_name"]
            print compound["r_i_docking_score"]

            chemblID = compound["s_m_entry_name"].split(".")[0]
            print chemblID

            if not chemblID in dscore.keys():
                dscore[chemblID] = float(compound["r_i_docking_score"])
            else:
                if float(compound["r_i_docking_score"]) < dscore[chemblID]:
                    dscore[chemblID] = float(compound["r_i_docking_score"])

        self.docking = dscore

        return self.docking



    def keepOnlyBestPoses(self):

        return

    def splitPoses(self, prDockingPose):

        self.prposes = prDockingPose
        self.lposefiles = []

        if not "lc" in dir(self):
            self.parseSDF()

        for pose in self.lc:
            namepose = pose["s_m_entry_name"]
            pfilout = prDockingPose + namepose + ".sdf"
            filout = open(pfilout, "w")
            filout.write(pose["sdf"])
            filout.close()

            self.lposefiles.append(pfilout)


    def removeDuplicate(self, colID):

        lduplicate = []
        i = 0
        imax = len(self.lc)
        while i < imax:
            if self.lc[i][colID] in lduplicate:
                del self.lc[i]
                imax = imax-1
                continue
            else:
                lduplicate.append(self.lc[i][colID])
                i += 1




        return