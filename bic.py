import sys
import basis
import EM
import math

bV = True

# return the best k for dataset D chosen
# by running EM with 10 starting points
# for each choice of k and choosing the
# k which minimizes BIC
def bicK(D):
    mEM = EM.EM(D)
    mEM.bPPC = False

    # -2 LL + k ln (n)
    bestk = 0
    bestBic = 0
    for ik, tryk in enumerate(range(2,11)):
        if bV:
            print "testing k = %d ..." % tryk

        bestround = 0
        for iterEM in range(10):
            mEM.lInitialCenters = []
            mEM.EM(tryk)
            if mEM.dEMLikelihood > bestround or iterEM == 0:
                bestround = mEM.dEMLikelihood

        bicVal = -2 * bestround + tryk * math.log(len(D.data))
        print bestround
        print bicVal
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
