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
      self.means = []
      self.stddevs = []
      self.classes = dict()
      self.cons = [] #constraints so far
      self.parseCsv(filename)
      self.poscons = [] #possible constraints
      self.zvalues()
      
   def addDatum(self, values, index):
      new_datum = datum()
      new_datum.index = index #Every data's id is its row number
      new_datum.cl = values[0]
      if(not (new_datum.cl in self.classes)):
         self.classes[new_datum.cl] = len(self.classes)
      new_datum.values = [float(x) for x in values[1:]]
      self.data.append(new_datum)
      
   def mean(self):
      sum = 0
      attCount = len(self.data[0].values)
      for j in range(0,attCount):
         for i in self.data:
            sum+=i.values[j] 
         self.means.append(sum/len(self.data))
         sum = 0
         
   def makeConst(self,numC):
      cons = []
      for i in range(numC):
         pair = self.poscons.pop(random.randint(0,len(self.poscons)-1))
         link = 0
         if(self.data[pair[0]].cl==self.data[pair[1]].cl):
            link = 1
         else:
            link = -1 
         pair += (link,)
         self.cons.append(pair)
         cons.append(pair)
      return cons
   def stddev(self):
      attCount = len(self.data[0].values)
      for j in range(0,attCount):
         sd = 0
         for i in self.data:
            sd += pow(i.values[j]-self.means[j],2)
         sd /= len(self.data)-1
         self.stddevs.append(pow(sd,0.5)) #Following the formula of standard deviation        
   def zvalues(self):
      self.mean()
      self.stddev()
      attCount = len(self.data[0].values)
      for i in self.data:
         for j in range(0,attCount):
            i.values[j]=(i.values[j]-self.means[j])/self.stddevs[j]
   def parseCsv(self,filename):
      with open(filename,"r") as fin:
         lines = fin.readlines()
      for i in range(1, len(lines)):
         values = lines[i].rstrip().split(",")
         self.addDatum(array(values),i)
   def parseConstraints(self,filename):
     with open(filename,"r") as fin:
        lines = fin.readlines()
     for i in range(0, len(lines)):
        values = lines[i].rstrip().split(",")
        self.poscons.append(array(values))


if __name__ == "__main__":
   if len(sys.argv) > 3:
      print("Error - usage is " + sys.argv[0] + " <data_file>")
      sys.exit(1)
      
   numClusters = 3
   m = cData(sys.argv[1])
   if len(sys.argv) == 2:
      m.poscons = [(i,j) for i in range(len(m.data)) for j in range(len(m.data))]
   else:
      m.parseConstraints(sys.argv[2])
   iteration = EM(m)
   print(len(m.classes))
   #iteration.bPPC = True
   for numCons in [0,1,25,50,100,200,500,1000,2000,5000,10000]:
   #for numCons in [0]:
       cons=m.makeConst(numCons)
       for i in cons:
           iteration.mCij[i[0]][i[1]] = i[2]
           iteration.mCij[i[1]][i[0]] = i[2]
       iteration.EM(len(m.classes))
       Estimated = np.ravel(iteration.mGammas.argmax(1).T)
       Real = array([ m.classes[i.cl] for i in m.data])
       print numCons, ",", nmi(Estimated,Real)
