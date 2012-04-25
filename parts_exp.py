import EM
import cData
import RepPoints
import SimulateUser
import random as R
import numpy as np
import Starts
import Cons
import utils
import sys
import pickle
import copy

numOuterTrials = 10
numInnerTrials = 10

# runs scenario 1 tests (just EM alone)
# * takes a cData object and a stub name for output files
# * returns a list of EMs from each of the outer trials
def scenario1(D, sFileStub):
    f = open(sFileStub + ".scen1.results", "w")

    EMResults = []
    for iOutertrial in range(numOuterTrials):
        f.write("outertrial: %d\n" % iOutertrial)
        f.write("likelihood,NMI\n")

        bestEM = []
        bestLikelihood = 0
        for iRestart in range(numInnerTrials):
            EMAlg = EM.cEM(D)
            EMAlg.bPPC = False
            EMAlg.EM(len(D.classlist))
            if iRestart == 0 or EMAlg.dEMLikelihood > bestLikelihood:
                bestLikelihood = EMAlg.dEMLikelihood
                bestEM = EMAlg

        EMResults.append(bestEM)

        f.write("%f,%f\n" % (bestLikelihood,
                             utils.evaluateEM_NMI(D, EMAlg) ) )
        f.flush()
    f.close()
    return EMResults



# runs scenario 2 tests (our cluster comparison method for
#  starting points)
# * takes a cData object and a stub name for output files
# * returns a dictionary from option to list of EMs from each
#   of the outer trials
class StarterOptions:
    # ConsWeights - can use different weights for the constraints
    # SplittingOnly - only feedback is taking items out of clusters
    # CrossClusterCons - feedback can say some point should be in
    #                    differenct cluster from where it is
    ConsWeights, SplittingOnly, CrossClusterCons = range(3)
    lOptions = [ #StarterOptions.ConsWeights,
                 SplittingOnly ] #,
                 #StarterOptions.CrossClusterCons ]
    lOptionNames = ["ConsWeights", "SplittingOnly", "CrossClusterCons"]
def scenario2(D, EMStarts, sFileStub):
    f = open(sFileStub + ".scen2.results", "w")

    dEMResults = {} # from option string to list of EMs
    f.write("queries,cons,likelihood,NMI\n")
    for option in StarterOptions.lOptions:
        EMOptionResults = []
        f.write(StarterOptions.lOptionNames[option] + "\n")
        for iOutertrial in range(numOuterTrials):
            f.write("outertrial: %d\n" % iOutertrial)
            f.flush()
            # run the initial points algorithm several
            # times for this set of options and this initial point
            #for iRestart in range(numInnerTrials):
            # get starting EM
            em = copy.deepcopy(EMStarts[iOutertrial])
            em.bPPC = True
            
            # setup and run goodInitial
            RP = RepPoints.RepPoints()
            emclusts = RP.createClusters(em)
            RP.repPoints(em, emclusts)
            
            starter = Starts.starts()
            newEM = starter.goodInitial(D, em, emclusts, RP, f)

            EMOptionResults.append(newEM)

        dEMResults[StarterOptions.lOptionNames[option]] = EMOptionResults
    f.close()
    return dEMResults


class TripConsOptions:
    # ConsWeights - can use different weights for the constraints
    # MidCons - make constraints to midpoints only
    # CenterChunkCons - make constraints to center 20%
    ConsWeights, MidCons, CenterChunkCons = range(3)
    lOptions = [ #ConsWeights,
                 #MidCons,
                 CenterChunkCons ]
    lOptionNames = ["ConsWeights", "MidCons", "CenterChunkCons"]
# use TripCons with just the EMStarts
# this is used for scenario3 and scenario4, thus the sNum indicator for logging
def TripConsTest(D, sNum, EMStarts, fp):
    fp.write("trips,queries,cons,likelihood,NMI\n")
    for option in TripConsOptions.lOptions:
        optname = TripConsOptions.lOptionNames[option]
        fp.write(optname + "\n")
        for iOutertrial in range(numOuterTrials):
            print "scenario ", sNum, " options ", optname, " outertrial ", iOutertrial
            fp.write("outertrial: %d\n" % iOutertrial)

            em = copy.deepcopy(EMStarts[iOutertrial])
            em.bPPC = True

            prevTrips = 0
            totalCons = 0
            nmiResult = utils.evaluateEM_NMI(D, em)
            fp.write("Initial nmi: %f\n" % nmiResult)
            consobj = Cons.cCons(D)
            for numTrips in range(1,len(D.data)/4,1):
                if option == TripConsOptions.CenterChunkCons:
                    consobj.constype = Cons.cCons.eConsType.TripCenterChunk
                elif option == TripConsOptions.MidCons:
                    consobj.constype = Cons.cCons.eConsType.TripMids

                print em.mLikelihood_il

                cons = consobj.tripCons(em.mGammas,numTrips-prevTrips)
                prevTrips = numTrips
                totalCons += len(cons)
                for i in cons:
                    em.mCij[i[0]][i[1]] = i[2]
                    em.mCij[i[1]][i[0]] = i[2] 
                em.EM(len(D.classlist))

                nmiresult = utils.evaluateEM_NMI(D, em)
                fp.write("%d,%d,%d,%f,%f\n" % (numTrips,
                                              numTrips*14,
                                              totalCons,
                                              em.dEMLikelihood,
                                              nmiresult) )
                fp.flush()

                if (nmiresult > 0.999 or len(D.data)==numTrips):
                    break

# scenario4 is just scenario3 with different initial points,
# but those points come in a dictionary, and we must redo
# the test for each key
def scenario4(D, dEMStarts, sFileStub):
    for k in dEMStarts:
        f = open(sFileStub + ".scen4." + k + ".results", "w")
        TripConsTest(D, 4, dEMStarts[k], f)
        f.close()


def errmsg():
    print "  Correct usage:"
    print "    %s result_nickname [input filename stub] [commands]" % sys.argv[0]
    print "  Where optional input filename stub (IFS) chooses data file:"
    print "     data/IFS.csv"
    print "  and [commands] is a string with pieces of the experiment"
    print "  to run, separated by _.  Commands available:"
    print "    s1 - scenario 1, including pickle"
    print "    p1 - load existing pickle for scenario 1"
    print "    s2, p2, s3, s4 as you'd expect"
    print "  for example, you would use p1_s2_s3_s4"
    print "  to load pickle from previous run of s1 then run"
    print "  the rest of the phases"

if __name__ == "__main__":
    sNickname = ""
    lDataFiles = ["winenorm3_pyre"]
    sCommands = "s1_s2_s3_s4"

    argc = len(sys.argv)
    if argc >= 2:
        sNickname = sys.argv[1]
    if argc >= 3: # includes filename stub of input
        lDataFiles = [ sys.argv[2] ]
    if argc == 4:
        sCommands = sys.argv[3]
    if argc < 2 or argc > 4:
        errmsg()
        exit(1)

    lCommands = sCommands.split("_")
    for dfname in lDataFiles:
        D = cData.cData("data/" + dfname + ".csv")

        # scenario 1
        if "s1" in lCommands:
            print "running scenario 1"
            EMStartsS1 = scenario1(D, "results/" + sNickname + "." + dfname)
            print "pickling results of scenario 1"
            pf = open("results/" + sNickname + "." + dfname + ".scen1.pickle", "w")
            pickle.dump(EMStartsS1, pf)
            pf.close()

        # load results from scenario 1 pickle
        if "p1" in lCommands:
            print "loading scenario 1 results from pickle"
            pf = open("results/" + sNickname + "." + dfname + ".scen1.pickle")
            EMStartsS1 = pickle.load(pf)
            pf.close()
        
        # scenario 2
        if "s2" in lCommands:
            print "running scenario 2"
            dEMStartsS2 = scenario2(D, EMStartsS1, "results/" + sNickname + "." + dfname)
            print "pickling results of scenario 2"
            pf = open("results/" + sNickname + "." + dfname + ".scen2.pickle", "w")
            pickle.dump(dEMStartsS2, pf)
            pf.close()

        # load results from scenario 2 pickle
        if "p2" in lCommands:
            pf = open("results/" + sNickname + "." + dfname + ".scen2.pickle")
            dEMStartsS2 = pickle.load(pf)
            pf.close()
            
        # scenario 3
        if "s3" in lCommands:
            print "running scenario 3"
            sFileStub =  "results/" + sNickname + "." + dfname
            f = open(sFileStub + ".scen3" + ".results", "w")
            TripConsTest(D, 3, EMStartsS1, f)
            
        # scenario 4
        if "s4" in lCommands:
            print "running scenario 4"
            scenario4(D, dEMStartsS2, "results/" + sNickname + "." + dfname)
        
