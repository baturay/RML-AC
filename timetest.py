import EM
import basis

import time
import sys
import pickle
import os
#import random

sfData = ""
lICenters = []
if (len(sys.argv)) < 2:
    sfData = "data/DS3_samp20.csv"
    #sfData = "data/winenorm3_pyre.txt"
else:
    sfData = sys.argv[1]
    
bNewCenters = False
if len(sys.argv) > 2:
    # check for 'new centers' bypass
    if sys.argv[2] == "__NEW__":
        bNewCenters = True
    else:
        lICenters = pickle.load(open(sys.argv[2]))
else:
    # check for default initial centers pickle
    sfDefltPkl = sfData + ".icenters.pickle"
    if os.access(sfDefltPkl, os.R_OK):
        lICenters = pickle.load(open(sys.argv[2]))
    else:
        bNewCenters = True

D = basis.cData(sfData)
PPC = EM.EM(D)
if bNewCenters:
    #PPC = D.emRestarts(10)
    #PPC.EM(len(D.data))
    PPC.EM(3)
    
    print PPC.lInitialCenters
    lICenters = PPC.lInitialCenters[:]

D.poscons = [(i,j) for i in range(len(D.data)) for j in range(len(D.data))]
cons = D.pairCons(500)
for i in cons:
    PPC.mCij[i[0]][i[1]] = i[2]
    PPC.mCij[i[1]][i[0]] = i[2]

print "starting with data file: " + sfData
print "starting with centers: ",
print ",".join( [ str(i) for i in lICenters ])
starttime = time.clock()
for iters in range(2):
    print "iters: ", iters
    PPC.lInitialCenters = lICenters[:]
    PPC.bPPC = True
    #PPC.EM(len(D.classlist))
    PPC.EM(3)
endtime = time.clock()

print "total time: " + str(endtime - starttime)
print "average time: " + str((endtime - starttime)/10.0)
