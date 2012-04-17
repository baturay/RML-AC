import numpy
import EM
import cData

D = cData.cData("data/winenorm3_pyre.txt")
M = EM.cEM(D)

M.bPPC = False

# run EM several times and get the likelihood
for iRestart in range(20):
    M.lInitialCenters = []
    M.EM(3)
    print M.dEMLikelihood,
    print " nmi: ",
    estimated = numpy.ravel(M.mGammas.argmax(1).T)
    print cData.nmi(estimated, D.real)

    
