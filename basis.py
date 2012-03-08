import sys
import random
from EM import *
from NMI import *
class datum:
   values = []
   name = ""
   cl = ""
class cData:
   def __init__(self,filename):
      self.data = []
      self.classlist = dict()
      self.cons = [] #constraints so far
      self.parseCsv(filename)
      self.poscons = [] #possible constraints
      self.consfile = 0
      self.classes = array([ self.classlist[i.cl] for i in self.data])
      self.constype = 0
      self.random = 0
      self.startpoins = []
   def addDatum(self, values, index):
      new_datum = datum()
      new_datum.index = index #Every data's id is its row number
      new_datum.cl = values[0]
      #Classlist keeps all possible classes.
      if(not (new_datum.cl in self.classlist)):
         self.classlist[new_datum.cl] = len(self.classlist)
      new_datum.values = [float(x) for x in values[1:]]
      self.data.append(new_datum)

   def makeConst(self,numC):
      if self.constype == 2:
         return pairCons(numC)
      elif self.constype == 3:
         return tripCons(numC)
      else:
         return []
         
   def pairCons(self,numC):
      for i in range(numC):
         pair = random.choice(self.poscons)
         self.poscons.remove(pair)
         if(not self.consfile):
            link = 0
            if(self.data[pair[0]].cl==self.data[pair[1]].cl):
               link = 1
            else:
               link = -1 
            pair += (link,)
         self.cons.append(pair)
         cons.append(pair)
      return cons
                 
   def parseCsv(self,filename):
      with open(filename,"r") as fin:
         lines = fin.readlines()
      for i in range(1, len(lines)):
         values = lines[i].rstrip().split(",")
         self.addDatum(array(values),i)
         
   def parseConstraints(self,filename):
     self.consfile = 1
     with open(filename,"r") as fin:
        lines = fin.readlines()
     for i in range(0, len(lines)):        
        values = lines[i].rstrip().split(",")
        self.poscons.append(array([ int(v) for v in values] ))
        
   def tripCons(self,gammadiffs,k):
      constraints = []
      link = 0
     
      for i in range(k):
         if(self.random):
            cons = choice(gammadiffs)
            gammadiffs.remove(cons)
         else:
            cons = gammadiffs.pop(0)
         class1 = filter(lambda x: x[1] == cons[1],gammadiffs)
         class2 = filter(lambda x: x[1] == cons[2],gammadiffs)
         class1 = class1[int(np.floor(0.8*len(class1))):]
         class2 = class2[int(np.floor(0.8*len(class2))):]

         if len(class1)>0 and (m.data[cons[3]].cl == str(m.data[class1[-1][3]].cl)):
            link = 1
         elif len(class2)>0 and (m.data[cons[3]].cl == str(m.data[class2[-1][3]].cl)):
            link = -1
       
         if(link != 0):
            for i in class1:
               constraints.append((cons[3],i[3],link))
            for i in class2:
               constraints.append((cons[3],i[3],-1*link))
      return constraints
      
   def findDiffs(self,gammas):
     gammadiffs = []
     maxindices = np.ravel(gammas.argmax(1).T)
     gammas = gammas.tolist()
     for i in range(len(gammas)):
         secondprob= -1
         firstindex = maxindices[i]
         firstprob = gammas[i][firstindex]
         for j in range(len(gammas[i])):
            if (gammas[i][j] > secondprob) and (j!=firstindex):
               secondindex = j
               secondprob = gammas[i][j]
         metric = (firstprob-secondprob)/(firstprob+secondprob)
         gammadiffs.append([metric,firstindex,secondindex,i])
     gammadiffs = sorted(gammadiffs,key=lambda x:x[0])
     return gammadiffs
     
   def emRestarts(self,k):
      centers = []
      maxEM = EM(self)
      maxNMI = self.evaluateEM(maxEM)
      for i in range(k):
         iteration = EM(self)
         nmi = self.evaluateEM(iteration)
         if maxNMI < nmi:
            maxNMI = nmi
            maxEM = iteration
      return maxEM
      
   def evaluateEM(self,em):
      em.EM(len(m.classlist))
      Estimated = np.ravel(em.mGammas.argmax(1).T)
      return nmi(Estimated,self.classes)
      
   def emInitialPoints(self,file):
      emstart = EM(self)
      with open(filename,"r") as fin:
         lines = fin.readlines()
      for i in range(1, len(lines)):
         values = lines[i].rstrip().split(",")
         self.startpoints.append(array([ float(v) for v in values] ))
      return emstart
if __name__ == "__main__":
   if len(sys.argv) != 4 and len(sys.argv) != 5:
      print("Error - usage is " + sys.argv[0] + " <data_file> <constrainttype> <constraintsource> <startingdatapoints>")
      sys.exit(1)   
   m = cData(sys.argv[1])
   
   if sys.argv[2] == "3":
      m.constype = 3
      if sys.argv[3] == "random":
         m.random = 1
         m.poscons = m.data[:]
      elif sys.argv[3] == "metric":
         m.poscons = m.data[:]
      else:
         parseConstraints(sys.argv[3])
   elif sys.argv[2] == "2":
      m.constype = 2
      if sys.argv[3] == "random" or sys.argv[3] == "metric":
         m.poscons = [(i,j) for i in range(len(m.data)) for j in range(len(m.data))]
      else:
         m.parseConstraints(sys.argv[3]) 
   else:
      m.constype = 1
   
   
   if len(sys.argv)>4:
      iteration = m.emInitial(file)
      iteration.lInitialCenters = m.startpoints.pop(0)      
   else:
      iteration = m.emRestarts(9)
   iteration.bPPC = True
   prevCons = 0

   totalcons = 0
   for numCons in range(5,len(m.data),5): 
      gammadiffs = m.findDiffs(iteration.mGammas)
      #The boolean at the end of next line is True for random
      #False for metric constraints.
      cons = m.tripCons(gammadiffs,numCons-prevCons)
      prevCons = numCons
      totalcons += len(cons)
      for i in cons:
        iteration.mCij[i[0]][i[1]] = i[2]
        iteration.mCij[i[1]][i[0]] = i[2]
      
      nmiresult = m.evaluateEM(iteration)
      print numCons, ",", nmiresult, ",", totalcons
      if(nmiresult == 1 or len(m.data)==numCons):
         break
