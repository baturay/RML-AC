import EM
from NMI import nmi
import numpy as np
from utils import evaluateEM_NMI
from utils import maybeWrite

class starts:
    def __init__(self):
        bOptions = True
                 
    # get an cEM object that had the best log likelihood
    # of numRuns runs of EM given cData D and EM's #clusters = k
    def emLLRestarts(self,D,numRuns,k):
        centers = []
        maxEM = cEM(D)
        maxEM.bPPC = False
        maxEM.EM(k)
        maxLL = maxEM.dEMLikelihood
        for i in range(numRuns - 1):
            iteration = cEM(D)
            iteration.EM(k)
            ll = iteration.dEMLikelihood
            if maxLL < ll:
                maxLL = ll
                maxEM = iteration
        return maxEM

    # using the representative points metric, find good initial centers
    # given cData object D and emclusters as [ emcluster ],
    # and RepPts as RepPoints object (with options filled in)
    # and fp as file to write results each iteration
    def goodInitial (self,D,em,emclusters,RepPts,fp):
        consistent = 0
        # Consistent means all the midpoints are same with the center.
        constraints = []
        iters = 0
        while consistent != len(emclusters) and iters < 5:
            consistent = 0
            resetCenters = []
            for ind,cl in enumerate(emclusters):
                realpoints = [D.data[i.index] for i in cl.midpoints]
                realcenter = D.data[cl.center.index]
                # Gets the real classes so we can do the simulation.
                rightclass = filter(lambda x: x.cl==realcenter.cl,realpoints)
                rightclass.append(realcenter)
                # Filters the midpoints same with center and then adds center.
                wrongclass = filter(lambda x: x.cl!=realcenter.cl,realpoints)
                # All the leftovers...
                if len(wrongclass) == 0:
                    consistent += 1
                else:
                    resetCenters.append(ind)
                # Cross constraints between right and wrong classes.
                for i in rightclass:
                    for j in realpoints:
                        if j in wrongclass:
                            constraints.append([i.index,j.index,-1])
                        elif j!= i:
                            constraints.append([i.index,j.index,1])
                for i in constraints:
                    em.mCij[i[0]][i[1]] = i[2]
                    em.mCij[i[1]][i[0]] = i[2]  
            print consistent
            # If all classes are not right, restart.
            if consistent != len(emclusters):    
                em.resetSomeCenters(em.lInitialCenters,resetCenters)
            em.EM(len(emclusters))
            emclusters = RepPts.createClusters(em)
            RepPts.repPoints(em, emclusters)
            print "goodInitial iter nmi: ", evaluateEM_NMI(D,em)
            iters += 1

            # queries,cons,likelihood,NMI
            maybeWrite(fp,
                       "%d,%d,%f,%f\n" % (iters*len(emclusters)*6,
                                          len(constraints),
                                          em.dEMLikelihood,
                                          evaluateEM_NMI(D,em) ) )
            
        return em   

    # returns a set of initial centers based on a clustering of
    # the centers of several initial clusterings
    # * D is the data (cData) object
    # * k is the number of classes (0 -> use #classes from D)
    def JLStartingPoint(D, k):
        if k == 0:
            k = len(D.classes)

        M = cEM(D)
        M.bPPC = False
        llCenters = []
        # get 20 different centers from running random-restarts of EM
        for iRestart in range(20):
            print "running EM ", iRestart, "..."
            M.lInitialCenters = []
            M.EM(k)
            llCenters.append(M.lCenters)

        # **** horrible hack - assumes this file exists because
        # it is not trivial to add a constructor that takes
        # a different type of data so a filename is needed
        D2 = cData("data/winenorm3_pyre.txt")
        D2.data = []
        print llCenters
        i = 0
        for center in llCenters:
            for V in center:
                D2.addDatum([0] + list(V), i)  # add 0 to beginning as class
                i += 1
            
        M2 = cEM(D2)
        M2.EM(k)
        return M2.lCenters    
        

