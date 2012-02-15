import numpy as np
import scipy as sp
from numpy import double
from numpy import matrix
from numpy import array
import scipy.spatial as sp_s
#from basis import cData
import random

# this module contains a class that performs EM
# clustering.  The input data is in the 'machine' format
# given in data.py

class EM:
    def __init__(self, _mData):
        self.mData = _mData

        self.lInitialCenters = []  # [ centers as [values] ]
        
        # iteration information
        self.numSteps = 0
        self.lLastCenters = []  # [ centers as [values] ]
        self.lCenters = [] # [ centers as [values] ]

        self.bPPC = True # use PPC
        self.mCij = []

        self.mGammas = [] # matrix prob of instance i in cluster j

    # returns true if we have reached convergence criteria
    def convergence(self):
        if len(self.lLastCenters) == 0:
            return False
        
        threshold = 0.00005
            
        diffs = [ sp_s.distance.euclidean(array(self.lLastCenters[i]),array(self.lCenters[i]))
                  for i in range(len(self.lCenters)) ]
        if max(diffs) < threshold:
            print "CONVERGED!!"
            return True

    
    # returns an array that is gamma_i,l values,
    # ie the probability that instance i is in cluster l
    # where i and l are indices into the data and centers lists resp.
    # inputs are number of data points, the current centers,
    # the P_l and the Sigma_l
    def clusterMembership(self, nData, nDataDim, lCenters, lPl, lSig, bPPC):
        coef = 1/(double(2*np.pi)**(nDataDim/double(2)))
        print "lsig: " , lSig[0].shape
        sigCoefs = [ 1/np.sqrt(np.linalg.det(sig)) for sig in lSig ]

        print "lpl: " , lPl

        ppc_lambda = 2
        # compute single gamma value
        def g():
            # gamma value for standard EM
            A = lPl[l] * coef * sigCoefs[l]
            B = matrix(array(self.mData.data[i].values) - array(lCenters[l]))
            C = matrix(lSig[l]).I
            g_EM =  A * np.exp(-0.5 * B * C * B.T )
            g_EM = g_EM[0,0]
            
            if bPPC:
                return g_EM * np.exp(2*ppc_lambda *
                                     np.sum( array(  [ 1,2,3  ] ) ) )

            return g_EM

        # |data| x |centers|
        G = matrix([ [ g()
                       for l in range(len(lCenters)) ]
                     for i in range(nData) ])

        # normalize each row ( over l for each i )
        rowsums = G.sum(axis=1)  # matrix  |row| x 1
        G = [ [ (G[i,l] / rowsums[i,0])
                for l in range(len(lCenters)) ]
              for i in range(nData) ]

        return matrix(G)
    
    def EM(self, numCenters):
        iterBound = 20
        
        # get the initial centers
        if self.lInitialCenters != []:
            if len(self.lInitialCenters) > numCenters:
                self.lInitialCenters = self.lInitialCenters[:numCenters]
                print "NOTE: given list of initial centers is too long, truncating"
            elif len(self.lInitialCenters) < numCenters:
                print "ERROR: provided too few initial centers"
                sys.exit(1)
        else:  # pick centers from data
            self.lInitialCenters = random.sample(self.mData.data, numCenters)
            self.lInitialCenters = [ c.values for c in self.lInitialCenters ]

        # initialization
        nDataDim = len(self.mData.data[0].values)
        nData = len(self.mData.data)
        self.lCenters = self.lInitialCenters[:]
        lPl = [ double(1)/len(self.lCenters) for l in range(len(self.lCenters)) ]
        lSig = [ np.eye(nDataDim) for i in range(len(self.lCenters)) ]
        lXi = [ self.mData.data[i].values for i in range(nData) ]

        iters = 0
        while iters < iterBound and not (iters != 1 and self.convergence()):
            self.lLastCenters = self.lCenters[:]
            print "new iter, last centers are: "
            print self.lLastCenters
            
            # estimate cluster membership
            # gets 2D array of Gamma_i,j
            # the prob that isntance i is in cluster j
            self.mGammas = self.clusterMembership(nData, nDataDim, self.lCenters, lPl, lSig, self.bPPC)
            
            # recalculate parameters
            lNl = [ np.sum(self.mGammas[:,l]) for l in range(len(self.lCenters)) ]
            lPl = [ lNl[l] / nData for l in range(len(self.lCenters)) ]
            self.lCenters = [ np.multiply(1/lNl[l],
                                          reduce(lambda x,y: x + y,
                                                 [ np.multiply(self.mGammas[i,l], lXi[i])
                                                   for i in range(nData) ],
                                                 0 ) )
                              for l in range(len(self.lCenters)) ]
            print "new centers: "
            print self.lCenters
            
            aXMuDiff = [ [ matrix(array(lXi[i]) - array(self.lCenters[l]))
                           for i in range(nData) ]
                         for l in range(len(self.lCenters)) ]

            lSig = [ np.multiply(1/lNl[l],
                                 reduce(lambda x,y: x+ y,
                                        [ np.multiply(self.mGammas[i,l],
                                                      aXMuDiff[l][i].T *
                                                      aXMuDiff[l][i]) 
                                          for i in range(nData) ] ) )
                     for l in range(len(self.lCenters)) ]

            iters += 1


def printDims(v, textv):
    print "dims ", textv, np.size(v,0), np.size(v,1)

def printDim(v, textv):
    print "dim ", textv, np.size(v,0)



if __name__ == "__main__":
    import basis
    D = basis.cData("data/winenorm3_pyre.txt")
    M = EM(D)
    M.EM(3)
    
