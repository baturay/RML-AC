import sys
from utils import evaluateEM_NMI
import cons
from cData import *
import RepPoints
import Starts

def parseCommandLine (D,argv):
   if len(argv) != 2 and len(argv) != 3:
      print("Error - usage is " + argv[0] + " <data_file> <startingdatapoints> ")
      sys.exit(1)   

   RP = RepPoints.RepPoints()
   EmAlg = cEM(D)
   if len(argv)>2:
      f = open(argv[2],"r")
      EmAlg = pickle.load(f)
      EmAlg.EM(len(D.classlist))
      EmAlg.bPPC = True
      emclusts = RP.createClusters(EmAlg)
      RP.repPoints(EmAlg, emclusts)
   else:
      EmAlg.EM(len(D.classlist))
      EmAlg.bPPC = True 
      # Creates clusters depending on what EM guessed.
      emclusts = RP.createClusters(EmAlg)
      # Finds the outerpoints and the midpoints and assigns them in emclusters.
      RP.repPoints(EmAlg, emclusts)
      # This makes the algorithm start with good initial points.
      starter = Starts.starts()
      EmAlg = starter.goodInitial(D, EmAlg, emclusts, RP)
      f = open("pickles/"+argv[1].split('/')[-1]+"pickle","w")
      pickle.dump(EmAlg,f)
      EmAlg.EM(len(D.classlist))
      f.close()
   
   return EmAlg
   

def main():
    # Takes in the file and parses into datum's.
    D = cData(sys.argv[1])
    # Starting points from the pickle or not.
    EmAlg = parseCommandLine(D, sys.argv)   
   
    prevCons = 0
    totalcons = 0
    nmiResult = evaluateEM_NMI(D, EmAlg)
    print "Initial nmi: ",nmiResult
    
    consobj = cons.cCons()
    for numCons in range(1,len(m.data)/4,1):
        cons = consobj.tripConsgamma(EmAlg.mGammas,numCons-prevCons)
        prevCons = numCons
        totalcons += len(cons)
        for i in cons:
            EmAlg.mCij[i[0]][i[1]] = i[2]
            EmAlg.mCij[i[1]][i[0]] = i[2] 
        EmAlg.EM(len(D.classlist))
        nmiresult = evaluateEM_NMI(D, EmAlg)
        print numCons, ",", nmiresult, ",", totalcons 
        if(nmiresult > 0.999 or len(D.data)==numCons):
            break

if __name__ == "__main__":
   main()
         
