import EM
import cData
D = cData.cData("data/DATASET3_trim.csv")
# load full constraint set into cData

# run PPC with same initial clusters with increasing number constraints
M = EM(D)

numCent = 6

for numCons in [500]:
    
