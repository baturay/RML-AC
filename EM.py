import numpy as np
import scipy as sp
from numpy import double
from numpy import matrix
from numpy import array
from numpy import mat
import scipy.spatial as sp_s
from cData import *
import random
import sys
import copy
import collections

# this module contains a class that performs EM
# clustering.  The input data is in the 'machine' format
# given in data.py

epsilon = 1e-250

class cEM:
    def __init__(self, _mData):
        self.mData = _mData

        self.lInitialCenters = []  # [ centers as [values] ]
        
        # iteration information
        self.numSteps = 0
        self.lLastCenters = []  # [ centers as [values] ]
        self.lCenters = [] # [ centers as [values] ]
        self.lUnUsedCenterInds = [] # indices of data already chosen as centers

        self.bPPC = False # use PPC

        # init Cij to 0's
        self.mCij =  [ [ 0 for i in range(len(_mData.data)) ]
                       for j in range(len(_mData.data)) ]

        # matrix prob/log_likelihood of instance i in cluster j
        self.mGammas = []
        self.mLikelihood_il = []

        self.bVerbose = False

        self.sErrInfo = ""
        self.saved_handler = np.seterrcall(self)
        self.save_err = np.seterr(all='log')
        self.numErrs = 0

        self.bEMLikelihoodEachStep = False
        self.dEMLikelihood = 0
        
    # error logging function for numpy (required name)
    def write(self, msg):
        if self.numErrs <2:
            sys.stderr.write("ERROR: %s" % msg)
            sys.stderr.write(" { %s }\n" % self.sErrInfo)
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
        sigCoefs = [ 1/np.sqrt(np.linalg.det(sig)) for sig in lSig ]

        # if we got inf in sigCoefs (caused by certain matrix conditions of sig)
        # change the formula used to calculate the coefficient from
        # 1/sqrt(det(sig)) to sqrt(10^m) * 1/sqrt(det(10*sig))
        for sigI in range(len(lSig)):
            if sigCoefs[sigI] == inf:
                sig = lSig[sigI]
                trymult = 10.0
                tryval = inf
                while tryval == inf:
                    #if trymult > 1000:
                    #    pdb.set_trace()
                    print "trymult,tryval, pow:",
                    print trymult,
                    print tryval,
                    print np.sqrt(pow(trymult,nDataDim))
                    tryval = np.sqrt(pow(trymult,nDataDim)) * (1/np.sqrt(np.linalg.det(trymult*sig)))
                    trymult *= 10
                sigCoefs[sigI] = tryval

        #if inf in sigCoefs:
        #    pdb.set_trace()

        if self.bVerbose:
            print "lsig: " , lSig[0].shape

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
            B = mat(self.mData.data[i].values) - mat(lCenters[l])
            C = mat(lSig[l])
            if C.min() == 0:
                 cis = np.where(C == 0)
                 C[cis] = np.exp(-745)

            try:
                C = C.I
            except:
                print "singular handled"
                print C
                #pdb.set_trace()
                C = C + epsilon * np.eye(len(C))
                C = C.I
                
            #g_EM =  A * np.exp(-0.5 * B * C * B.T )
            
            g_EM =  np.log(A) + B * C * B.T * -0.5
            g_EM = g_EM[0,0]
            
            if bPPC:
                g_PPC = 2*ppc_lambda * \
                           np.sum( array(  [ self.mCij[i][j] * G_old[j,l]
                                             for j in range(nData) if i != j ]
                                           ) )
                return g_PPC + g_EM
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
            # say M is the max value such that exp(M) < inf
            # we must avoid any items of LL in G being >= M.
            # further, we must avoid M - log(k), since if there
            # are k elements of M, we get a row sum, after the
            # exponent is applied, of k exp(M).  To avoid this,
            # me make sure the max of any row's LL is M - log(k) -1
            # since k exp(M-log(k)) = exp(M)  ( 1 is for rounding error)
            rowmaxs = G.max(axis=1)
            k = G.shape[1]
            maxval = 700 - np.log(k) # was 709

            # fix G before applying exponent (as explained above)
            # * this syntax lets us fix each row and stack them back
            #   into a matrix
            G = np.vstack( [ G[i] - (rowmaxs[i]-(maxval-1))
                             if rowmaxs[i] >= maxval
                             else G[i]
                             for i in range(len(G)) ] )

            G = np.exp(G)
            rowsums = G.sum(axis=1)  # matrix  |row| x 1
            #print "rowsums: ", rowsums
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
            
            G = mat(G_old).copy()
            try:
                G = mat([ [ g(i, l, bPPC) if iters > 0 else g(i,l,False)
                            for l in range(len(lCenters)) ]
                          for i in range(nData) ] )
                #G = mat( np.fromfunction( np.vectorize(lambda i, j: g(i, j, bPPC) if iters > 0 else g(i, j, False)),
                #                          (nData, len(lCenters)) ) )
            except np.linalg.linalg.LinAlgError:
                print "Singular matrix: moving on"
                print "  problem at iter %d of EM" % iters
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
    # Cannot include indices in lExclusions
    def resetSomeCenters(self, lCenters, lIndices, lExclusions):
        if self.lUnUsedCenterInds == []:
            self.lUnUsedCenterInds = range(len(self.mData.data))

        lAvailable = self.lUnUsedCenterInds[:]
        for i in lExclusions:
            if i in lAvailable:
                lAvailable.remove(i)
        
        for i in lIndices:
            icenteri = random.sample(range(len(lAvailable)), 1)
            icenteri = icenteri[0]  # random.sample returns a list
            icenter = lAvailable.pop(icenteri)
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
            self.resetSomeCenters(self.lInitialCenters, range(numCenters), [])
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
            lNl[lnli] = np.exp(-709)
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
                                 reduce(lambda x,y: x + y,
                                        [ np.multiply(self.mGammas[i,l],
                                                      aXMuDiff[l][i].T *
                                                      aXMuDiff[l][i]) 
                                          for i in range(nData) ] ) )
                     for l in range(len(self.lCenters)) ]
            # in all Sig matrices, eliminate diagonal entries that
            # are too low (they cause underflow type errors in inversion)
            # * 1e-308 is the magic number such that if the last row and
            #   column are zeros except this number is on the diagonal
            #   the matrix can be inverted cleanly (any lower and there
            #   are a lot of nan values)
            for S in lSig:
                ww = np.where(S < 1e-308)
                for w in range(ww[0].shape[1]):
                    wi = ww[0][0,w]
                    wj = ww[1][0,w]
                    if wi == wj:
                        S[wi,wj] = 1e-308
                if np.linalg.det(S) == 0:
                    S = S + epsilon * np.eye(len(S))
            
            # supposedly faster alternative... slower in my tests
            #lSig = [ np.multiply(1/lNl[l],
            #                     np.add.reduce(
            #                       np.fromfunction( np.vectorize(
            #                         lambda i:
            #                           np.multiply(self.mGammas[int(i),l],
            #                                          aXMuDiff[l][int(i)].T *
            #                                          aXMuDiff[l][int(i)]) ),
            #                         (nData,) ) ))
            #                     for l in range(len(self.lCenters)) ]

            mems = self.Membership()
            countcol = collections.Counter(mems)
            print "membership: ", [ cnt for elt,cnt in countcol.most_common() ]

            if self.bEMLikelihoodEachStep:
                ll = EMLikelihood()
                if bVerbose:
                    print "likelihood at step " + str(iters) + ": " + str(ll)

            iters += 1

        # calculate the EM Likelihood
        self.EMLikelihood()

    # internal function for saving the LL of total EM
    def EMLikelihood(self):
        # formula is sum over n,l of gamma(i,l) * z(i,l) where
        # z = 1 if i is in cluster k, else 0
        nData = len(self.mData.data)
        LL = 0
        membership = np.ravel(self.mGammas.argmax(1).T)
        for i in range(nData):
            for k in range(len(self.lCenters)):
                if membership[i] == k:
                    LL += self.mLikelihood_il[i,k]

        self.dEMLikelihood = LL
        return LL

    # give the class labels of each data item
    def Membership(self):
        return np.ravel(self.mGammas.argmax(1).T)

def printDims(v, textv):
    print "dims ", textv, np.size(v,0), np.size(v,1)

def printDim(v, textv):
    print "dim ", textv, np.size(v,0)

if __name__ == "__main__":
    import cData
    D = cData.cData("data/winenorm3_pyre.txt")
    M = cEM(D)
    M.EM(3)
    
