import numpy as np
from NMI import nmi

def evaluateEM_NMI(D,em):
    Estimated = np.ravel(em.mGammas.argmax(1).T)
    return nmi(Estimated,D.classes)
