from EM import *
from cData import *

if len(sys.argv) > 3:
    print("Error - usage is " + sys.argv[0] + " <data_file>")
    sys.exit(1)
    
    
    m = cData(sys.argv[1])
    if len(sys.argv) == 2:
        m.poscons = [(i,j) for i in range(len(m.data)) for j in range(len(m.data))]
    else:
        m.parseConstraints(sys.argv[2]) 
    iteration = m.emRestarts(10)
    iteration.bPPC = True
    prevCons = 0
    for numCons in [0,25,50,100,200,500,1000,2000,5000,10000]:
    #for numCons in range(0,len(m.poscons)):
        cons=m.makeConst(numCons-prevCons)
        prevCons = numCons
        for i in cons:
            iteration.mCij[i[0]][i[1]] = i[2]
            iteration.mCij[i[1]][i[0]] = i[2]
        nmiresult = m.evaluateEM(iteration)
        print numCons, ",", nmiresult
        if(nmiresult == 1):
            break
