import numpy as np
import scipy as sp
from numpy import double
from numpy import matrix
from numpy import array
#from numpy import *
#from scipy import *
from basis import *

# this module contains a class that performs EM
# clustering.  The input data is in the 'machine' format
# given in data.py

class EM:
    lInitialCenters = []  # [ centers as [values] ]
    mData = cData()

    # iteration information
    numSteps = 0
    lLastCenters = []  # [ centers as [values] ]
    lCenters = [] # [ centers as [values] ]
    
    def __init__(self, _mData):
        self.mData = _mData


    # returns true if we have reached convergence criteria
    def convergence(self):
        threshold = 0.5
        diffs = [ scipy.spatial.distance.euclidean(array(lLastCenters[i]),array(lCenters[i]))
                  for i in len(lCenters) ]
        if max(diffs) < threshold:
            return true

    
    # returns an array that is gamma_i,l values,
    # ie the probability that instance i is in cluster l
    # where i and l are indices into the data and centers lists resp.
    # inputs are number of data points, the current centers,
    # the P_l and the Sigma_l
    def clusterMembership(self, nData, nDataDim, lCenters, lPl, lSigs):
        coef = 1/(double(2*np.pi)**(nDataDim/double(2)))
        sigCoefs = [ 1/np.sqrt(np.linalg.det(sig)) for sig in lSigs ]
        
        G = [ [  lPl[l] * coef * sigCoefs[l] *
                 np.exp(-0.5 *
                        matrix(array(self.mData.data[i].values) - array(lCenters[l]))
                        * lSigs[l].I *
                        matrix(array(self.mData.data[i].values) - array(lCenters[l])).T )
                 for i in range(nData) ]
              for l in range(len(lCenters)) ]

        # normalize each row
        rowsums = [ np.sum(R) for R in G ]
        G = [ [ G[i][l] / rowsums[i]
                for i in range(len(lCenters)) ]
              for l in range(len(nData))  ]
        
        return array(G)

    def EM(self, numCenters):
        # get the initial centers
        if lInitialCenters != []:
            if len(lInitialCenters) > numCenters:
                lInitialCenters = lInitialCenters[:numCenters]
                print "NOTE: given list of initial centers is too long, truncating"
            elif len(lInitialCenters < numCenters):
                print "ERROR: provided too few initial centers"
                sys.exit(1)
        else:  # pick centers from data
            lInitialCenters = random.sample(mData.data, numCenters)

        # initialization
        nDataDim = len(mData.data[0].values)
        nData = len(mData.data)
        lCenters = [ c.values for c in lInitialCenters ]
        lPl = [ double(1)/len(lCenters) ]
        lSig = [ eye(nDataDim) for i in range(len(lCenters)) ]

        while not convergence():
            lLastCenters = lCenters[:]
            print "new iter, last centers are: "
            print lLastCenters
            
            # estimate cluster membership
            # gets 2D array of Gamma_i,j
            # the prob that isntance i is in cluster j
            aGammas = clusterMembership(nData, nDataDim, lCenters, lPl, lSigs)
            
            # recalculate parameters
            lNl = [ sum(aGammas[:,l]) for l in range(len(lCenters)) ]
            lPl = [ lNl[l] / nData for l in range(len(lCenters)) ]
            lCenters = [ 1/lNl[l] * sum( [ aGammas[i,l] * array(Xi.values) for Xi in mData.data] )
                         for l in range(len(lCenters)) ]
            aXMuDiff = [ [ matrix(array(self.mData.data[i].values) - array(lCenters[l]))
                           for i in range(nData) ]
                         for l in range(len(lCenters)) ]
            lSig = [ 1/lNl[l] *
                     sum( [ aGammas[i,l] * aXMuDiff[i][l].T * aXMuDiff[i][l]
                            for i in range(nData) ] )
                     for l in range(len(lCenters)) ]
