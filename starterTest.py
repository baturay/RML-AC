# this test is designed to see how effective our
# starting point selection is on its own.
# After several rounds of having the user provide
# feedback on starting points, we proceed with
# random pairwise constraints and watch to see our
# NMI increase

import pickle
import EM
import cData
import sys

if len(sys.argv) < 2:
    print "provide filename and optionally a filename for the pickle of centers"
    exit(1)

# use the same code for getting initial points as baturay did
D = cData.cData(sys.argv[1])
D.setType("2", "random")
EmAlg = EM.cEM(D)
EmAlg.EM(len(D.classlist))
EmAlg.bPPC = True 
#Creates clusters depending on what EM guessed.
D.createClusters(EmAlg)
#Finds the outerpoints and the midpoints and assigns them in emclusters.
D.repPoints(EmAlg)
#This makes the algorithm start with good initial points.
EmAlg = D.goodInitial(EmAlg)

print "pickling starting position to: ",
picklefname = "pickles/"+sys.argv[1].split('/')[-1]+".pickle"
if len(sys.argv) > 2:
    picklefname = sys.argv[2]
f = open(picklefname,"w")
l = EmAlg.lCenters
pickle.dump(l,f)

nmiResult = D.evaluateEM(EmAlg)
print "Initial nmi: ",nmiResult

consperstep = 20
for iters in range(10):
    EmAlg.lInitialCenters = l[:]
    for numCons in range(1,len(D.data), consperstep):
        cons = D.pairCons(consperstep)
        for i in cons:
            EmAlg.mCij[i[0]][i[1]] = i[2]
            EmAlg.mCij[i[1]][i[0]] = i[2]

        EmAlg.EM(len(D.classlist))

        nmiresult = D.evaluateEM(EmAlg)
        print numCons, ",", nmiresult
        if(nmiresult > 0.999):
            break
