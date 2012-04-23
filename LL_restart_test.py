import numpy
import EM
import cData
import utils

# run EM several times and get the likelihood
for iRestart in range(20):
    D = cData.cData("data/winenorm3_pyre.csv")
    # D = cData.cData("data/normvert.csv")
    M = EM.cEM(D)
    M.bPPC = False

    M.EM(3)
    print M.dEMLikelihood,
    print " nmi: ",
    print utils.evaluateEM_NMI(D, M)

    
