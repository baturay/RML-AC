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
      cons = []
      for i in range(numC):
         if(len(self.poscons)==1):
            pair = self.poscons.pop()
         else:
            pair = self.poscons.pop(random.randint(0,len(self.poscons)-1))
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
        
   def tripCons(self,gammadiffs,k,rand=False):
      constraints = []
      link = 0
     
      for i in range(k):
         if(rand):
            cons = gammadiffs.pop(random.randint(0,len(gammadiffs)-1))
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
         print nmi
         if maxNMI < nmi:
            maxNMI = nmi
            maxEM = iteration
      print "maxNMI",maxNMI
      return maxEM
      
   def evaluateEM(self,em):
      em.EM(len(m.classlist))
      Estimated = np.ravel(em.mGammas.argmax(1).T)
      return nmi(Estimated,self.classes)
      
      
      
if __name__ == "__main__":
   if len(sys.argv) > 3:
      print("Error - usage is " + sys.argv[0] + " <data_file>")
      sys.exit(1)   
   m = cData(sys.argv[1])
   
   if len(sys.argv) == 2:
      m.poscons = m.data[:]
   else:
      m.parseConstraints(sys.argv[2]) 
      
   iteration = m.emRestarts(9)
   iteration.EM(len(m.classlist))
   #iteration.lInitialCenters = [m.data[5],m.data[76],m.data[160]]
   iteration.bPPC = True
   prevCons = 0

   #for numCons in range(len(m.data)):
      # print i, " ",m.data[i].values," ",m.data[i].cl
   totalcons = 0
   for numCons in range(5,len(m.data),5): 
      gammadiffs = m.findDiffs(iteration.mGammas)
      #The boolean at the end of next line is True for random
      #False for metric constraints.
      cons = m.tripCons(gammadiffs,numCons-prevCons,False)
      prevCons = numCons
      totalcons += len(cons)
      for i in cons:
        iteration.mCij[i[0]][i[1]] = i[2]
        iteration.mCij[i[1]][i[0]] = i[2]
      
      nmiresult = m.evaluateEM(iteration)
      print numCons, ",", nmiresult, ",", totalcons
      if(nmiresult == 1 or len(m.data)==numCons):
         break
