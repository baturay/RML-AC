import sys
import basis
import EM
import math
import numpy as np

bV = True

# calculate RSS value needed for BIC
# RSS = sum over i, l: z_il (1- p(i in l))^2
def RSS(mEM):
    # z = 1 iff i is in cluster k, else 0
    nData = len(mEM.mData.data)
    rss = 0
    membership = np.ravel(mEM.mGammas.argmax(1).T)
    for i in range(nData):
        for k in range(len(mEM.lCenters)):
            if membership[i] == k:
                rssp = (1 - np.exp(mEM.mLikelihood_il[i,k]))
                rss += rssp * rssp

    return rss

# return the best k for dataset D chosen
# by running EM with 10 starting points
# for each choice of k and choosing the
# k which minimizes BIC
def bicK(D):
    mEM = EM.EM(D)
    mEM.bPPC = False

    bestk = 0
    bestBic = 0
    for ik, tryk in enumerate(range(2,11)):
        if bV:
            print "testing k = %d ..." % tryk

        bestround = 0
        print "likilihoods: ",
        for iterEM in range(10):
            mEM.lInitialCenters = []
            mEM.EM(tryk)
            print mEM.dEMLikelihood,
            if mEM.dEMLikelihood > bestround or iterEM == 0:
                bestround = mEM.dEMLikelihood
        print ""

        # -2 LL + k ln(n)
        #bicVal = -2 * bestround + tryk * math.log(len(D.data))
        # BIC = N log(RSS/N) + K log(N)
        N = len(D.data)
        rss = RSS(mEM)
        bicVal = N * math.log(rss/N) + tryk * math.log(N)
        print "best round: ", bestround
        print "rss: ", rss
        print "bicVal: ", bicVal
        if bicVal < bestBic or ik == 0:
            bestBic = bicVal
            bestk = tryk
            if bV:
                print "new best: %d at %d" % (bestBic, bestk)

    return bestk


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print "usage: " + sys.argv[0] + " filename"
        sys.exit(1)
    
    D = basis.cData(sys.argv[1])
    print "best k %d:" % bicK(D)
