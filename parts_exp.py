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

numOuterTrials = 10

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
        for iRestart in range(10):
            EMAlg = EM.cEM(D)
            EMAlg.bPPC = False
            EMAlg.EM(len(D.classlist))
            if iRestart == 0 or EMAlg.dEMLikelihood > bestLikelihood:
                bestLikelihood = EMAlg.dEMLikelihood
                bestEM = EMAlg

        EMResults.append(bestEM)

        f.write("%f,%f\n" % (bestLikelihood,
                             utils.evaluateEM_NMI(D, EMAlg) ) )
    f.close()
    return EMResults



# runs scenario 1 tests (our cluster comparison method for
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
    lOptions = [ StarterOptions.ConsWeights,
                 StarterOptions.SplittingOnly,
                 StarterOptions.CrossClusterCons ]
    lOptionNames = ["ConsWeights", "SplittingOnly", "CrossClusterCons"]
def scenario2(D, sFileStub):
    f = open(sFileStub + ".scen2.results", "w")

    EMResults = {}
    for option in StarterOptions.lOptions:
        EMOptionResults = []
        for iOutertrial in range(numOuterTrials):
            f.write("outertrial: %d\n" % iOutertrial)
            f.write("likelihood,NMI\n")

            bestEM = []
            bestLikelihood = 0
            for iRestart in range(10):
                EMAlg = EM.cEM(D)
                EMAlg.bPPC = False
                EMAlg.EM(len(D.classlist))
                if iRestart == 0 or EMAlg.dEMLikelihood > bestLikelihood:
                    bestLikelihood = EMAlg.dEMLikelihood
                    bestEM = EMAlg

            EMResults.append(bestEM)

            f.write("%f,%f\n" % (bestLikelihood,
                                 utils.evaluateEM_NMI(D, EMAlg) ) )

        EMResults[StarterOptions.lOptionNames[option]] = EMOptionResults
    f.close()
    return EMResults



if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Please specify a nickname for the results"
        exit(1)
    sNickname = sys.argv[1]
    
    for dfname in ["winenorm3_pyre"]:
        D = cData.cData("data/" + dfname + ".txt")

        #EMStartsS1 = scenario1(D, "results/" + sNickname + "." + dfname)
        #pf = open("results/" + sNickname + "." + dfname + ".scen1.pickle", "w")
        #pickle.dump(EMStartsS1, pf)
        #pf.close()

        # scenario 2
        EMStartsS2 = scenario2(D, "results/" + sNickname + "." + dfname)
        pf = open("results/" + sNickname + "." + dfname + ".scen2.pickle", "w")
        pickle.dump(EMStartsS2, pf)
        pf.close()

        
        # load results from scenario 1 pickle
        #pf = open("results/" + sNickname + "." + dfname + ".scen1.pickle")
        #EMStarts = pickle.load(pf)
        #pf.close()
