import numpy as np
import scipy as sp
from numpy import double
from numpy import matrix
from numpy import array
from numpy import mat
import scipy.spatial as sp_s
#from basis import cData
import random
import sys

# this module contains a class that performs EM
# clustering.  The input data is in the 'machine' format
# given in data.py

class cEM:
    def __init__(self, _mData):
        self.mData = _mData

        self.lInitialCenters = []  # [ centers as [values] ]
        
        # iteration information
        self.numSteps = 0
        self.lLastCenters = []  # [ centers as [values] ]
        self.lCenters = [] # [ centers as [values] ]
        self.lUsedCenterInds = [] # indices of data already chosen as centers

        self.bPPC = False # use PPC

        # init Cij to 0's
        self.mCij =  [ [ 0 for i in range(len(_mData.data)) ]
                       for j in range(len(_mData.data)) ]

        # matrix prob/log_likelihood of instance i in cluster j
        self.mGammas = []
        self.mLikelihood_il = []

        self.bVerbose = False

        self.sErrInfo = ""
        #self.saved_handler = np.seterrcall(self)
        #self.save_err = np.seterr(all='log')
        self.numErrs = 0

        self.bEMLikelihoodEachStep = False
        self.dEMLikelihood = 0
        
    # error logging function for numpy (required name)
    def write(self, msg):
        if self.numErrs <2:
            sys.stderr.write("ERROR: %s\n" % msg)
            sys.stderr.write("   { %s }\n" % self.sErrInfo)
        self.numErrs += 1


    # returns true if we have reached convergence criteria
    def convergence(self):
        if len(self.lLastCenters) == 0:
            return False
        
        threshold = 0.00005
            
        diffs = [ sp_s.distance.euclidean(array(self.lLastCenters[i]),array(self.lCenters[i]))
                  for i in range(len(self.lCenters)) ]
        if max(diffs) < threshold:
            sys.stderr.write("EM CONVERGED!!\n")
            return True

    
    # returns an array that is gamma_i,l values,
    # ie the probability that instance i is in cluster l
    # where i and l are indices into the data and centers lists resp.
    # inputs are number of data points, the current centers,
    # the P_l and the Sigma_l
    def clusterMembership(self, nData, nDataDim, lCenters, lPl, lSig, bPPC):
        coef = 1/(double(2*np.pi)**(nDataDim/double(2)))
        if self.bVerbose:
            print "lsig: " , lSig[0].shape
        sigCoefs = [ 1/np.sqrt(np.linalg.det(sig)) for sig in lSig ]

        if self.bVerbose:
            print "lpl: " , lPl

        ppc_lambda = 1
        # compute single gamma value
        #  (return log of that gamma value)
        def g(di, dj, bPPC):
            i = int(di)
            l = int(dj)
            
            # gamma value for standard EM
            A = lPl[l] * coef * sigCoefs[l]
            #self.sErrInfo = "(abits - %f * %f * %f)" % (lPl[l], coef, sigCoefs[l])
            B = mat(self.mData.data[i].values) - mat(lCenters[l])
            C = mat(lSig[l])
            if C.min() == 0:
                 cis = np.where(C == 0)
                 C[cis] = np.exp(-745)

            try:
                C = C.I
            except:
                sys.stderr.write("Singular Matrix C, moving on\n")
                
            #self.sErrInfo = self.sErrInfo + "  A: %f, D: %f" % (A, -0.5 * B * C * B.T)
            #self.sErrInfo = self.sErrInfo + "um1"
            #g_EM =  A * np.exp(-0.5 * B * C * B.T )
            
            g_EM =  np.log(A) + B * C * B.T * -0.5
            #self.sErrInfo += "um2"
            g_EM = g_EM[0,0]
            
            if bPPC:
                g_PPC = 2*ppc_lambda * \
                           np.sum( array(  [ self.mCij[i][j] * G_old[j,l]
                                             for j in range(nData) if i != j ]
                                           ) )
                return g_PPC + g_EM
            #self.sErrInfo += "um3"
            return g_EM

        def gammaConverge():
            # compare new with old gammas
            threshold = np.exp(-20)
            step = np.max(np.abs(G - G_old))
            if step < threshold and self.bVerbose:
                sys.stderr.write("gammaconverged: " + str(step) + "\n")
            return step < threshold

        # normalize each row ( over l for each i )
        # also apply exp
        def normalize(G):
            #rowmins = G.min(axis=1)
            rowmaxs = G.max(axis=1)
            #rowranges = rowmaxs - rowmins
            # 18 is half difference between underflow and overflow abs values
            G = G - (rowmaxs - 700)
            G = np.exp(G)
            rowsums = G.sum(axis=1)  # matrix  |row| x 1
            
            nG = G / rowsums
            if nG.min() == 0:
                gis = np.where(nG == 0)
                nG[gis] = np.exp(-745)

            return nG

        # |data| x |centers|
        iterBound = 20
        if not bPPC:
            iterBound = 1
        iters = 0
        G_old = self.mGammas

        if self.bVerbose:
            print "ppc converge...",
        
        bGammaConverged = False
        while iters < iterBound and not (iters != 0 and bGammaConverged):
            if self.bVerbose:
                print iters, ",",

            try:
                G = mat([ [ g(i, l, bPPC) if iters > 0 else g(i,l,False)
                            for l in range(len(lCenters)) ]
                          for i in range(nData) ] )
                #G = mat( np.fromfunction( np.vectorize(lambda i, j: g(i, j, bPPC) if iters > 0 else g(i, j, False)),
                #                          (nData, len(lCenters)) ) )
            except numpy.linalg.linalg.LinAlgError:
                print "Singular matrix: moving on"
                break

            self.mLikelihood_il = G.copy()
            G = normalize(G)
            if iters > 0:
                bGammaConverged = gammaConverge()
            G_old = G.copy()

            iters += 1

        if self.bVerbose:
            sys.stderr.write("ppcifinal " + str(iters) + "\n")

        return mat(G)
    
    # in the list of centers _lCenters_, swap out the items at indices in the
    # list _lIndices_ with a randomly chosen data point
    def resetSomeCenters(self, lCenters, lIndices):
        if self.lUsedCenterInds == []:
            self.lUsedCenterInds = range(self.mData.data)
        
        for i in lIndices:
            icenteri = random.sample(range(len(self.lUsedCenterInds)), 1)
            icenter = self.lUsedCenterInds.pop(icenteri)
            lCenters[i] = self.mData.data[icenter].values

    # run the EM algorithm looking for given number of clusters
    def EM(self, numCenters):
        iterBound = 20
        
        # get the initial centers
        if self.lInitialCenters != []:
            if self.bVerbose:
                print "need initicenters ", numCenters, " ", len(self.lInitialCenters)
            if len(self.lInitialCenters) > numCenters:
                self.lInitialCenters = self.lInitialCenters[:numCenters]
                if self.bVerbose:
                    sys.stderr.write("NOTE: given list of initial centers is too long, truncating\n")
            elif len(self.lInitialCenters) < numCenters:
                if self.bVerbose:
                    sys.stderr.write("ERROR: provided too few initial centers\n")
                sys.exit(1)
        else:  # pick centers from data
            self.lInitialCenters = range(numCenters)
            self.resetSomeCenters(self.lInitialCenters, range(numCenters))
            if self.bVerbose:
                print "initcenters ", self.lInitialCenters

        # initialization
        nDataDim = len(self.mData.data[0].values)
        nData = len(self.mData.data)
        self.lCenters = self.lInitialCenters[:]
        lPl = [ double(1)/len(self.lCenters) for l in range(len(self.lCenters)) ]
        lSig = [ np.eye(nDataDim) for i in range(len(self.lCenters)) ]
        lXi = [ self.mData.data[i].values for i in range(nData) ]

        iters = 0
        while iters < iterBound and not (iters != 0 and self.convergence()):
            if self.bVerbose:
                print "aGamma convergence: ", iters
            self.lLastCenters = self.lCenters[:]
            if self.bVerbose:
                print "new iter, last centers are: "
                print self.lLastCenters
            
            # estimate cluster membership
            # gets 2D array of Gamma_i,j
            # the prob that isntance i is in cluster j
            self.mGammas = self.clusterMembership(nData, nDataDim, self.lCenters, lPl, lSig, self.bPPC)
            
            # recalculate parameters
            lNl = array([ np.sum(self.mGammas[:,l]) for l in range(len(self.lCenters)) ])
            lnli = np.where(lNl == 0)
            lNl[lnli] = np.exp(-745)
            lPl = [ lNl[l] / nData for l in range(len(self.lCenters)) ]
            # ****** lcenters as matrix?
            #self.lCenters = [ np.multiply(1/lNl[l],
            #                              reduce(lambda x,y: x + y,
            #                                     [ np.multiply(self.mGammas[i,l], lXi[i])
            #                                       for i in range(nData) ],
            #                                     0 ) )
            #                  for l in range(len(self.lCenters)) ]
            self.lCenters = [ np.multiply(1/lNl[l],
                                          reduce(lambda x,y: x + y,
                                                 [ np.multiply(self.mGammas[i,l], lXi[i])
                                                   for i in range(nData) ],
                                                 0 ) )
                              for l in range(len(self.lCenters)) ]

            if self.bVerbose:
                print "new centers: "
                print self.lCenters
            
            aXMuDiff = [ [ mat(lXi[i]) - mat(self.lCenters[l])
                           for i in range(nData) ]
                         for l in range(len(self.lCenters)) ]

            # ***** OPT... each elem should be a matrix
            lSig = [ np.multiply(1/lNl[l],
                                 reduce(lambda x,y: x+ y,
                                        [ np.multiply(self.mGammas[i,l],
                                                      aXMuDiff[l][i].T *
                                                      aXMuDiff[l][i]) 
                                          for i in range(nData) ] ) )
                     for l in range(len(self.lCenters)) ]
            #lSig = [ np.multiply(1/lNl[l],
            #                     np.add.reduce(
            #                       np.fromfunction( np.vectorize(
            #                         lambda i:
            #                           np.multiply(self.mGammas[int(i),l],
            #                                          aXMuDiff[l][int(i)].T *
            #                                          aXMuDiff[l][int(i)]) ),
            #                         (nData,) ) ))
            #                     for l in range(len(self.lCenters)) ]

            if self.bEMLikelihoodEachStep:
                ll = EMLikelihood()
                if bVerbose:
                    print "likelihood at step " + str(iters) + ": " + str(ll)

            iters += 1

        # calculate the EM Likelihood
        self.EMLikelihood()

    def EMLikelihood(self):
        # formula is sum over n,l of gamma(i,l) * z(i,l) where
        # z = 1 iff i is in cluster k, else 0
        nData = len(self.mData.data)
        LL = 0
        membership = np.ravel(self.mGammas.argmax(1).T)
        for i in range(nData):
            for k in range(len(self.lCenters)):
                if membership[i] == k:
                    LL += self.mLikelihood_il[i,k]

        self.dEMLikelihood = LL
        return LL

def printDims(v, textv):
    print "dims ", textv, np.size(v,0), np.size(v,1)

def printDim(v, textv):
    print "dim ", textv, np.size(v,0)



if __name__ == "__main__":
    import basis
    D = basis.cData("data/winenorm3_pyre.txt")
    M = EM(D)
    M.EM(3)
    
