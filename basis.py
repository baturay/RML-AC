import sys
import numpy
import scipy
import random
from EM import *
class datum:
   values = []
   name = ""
   cl = ""
class cData:
   def __init__(self):
      self.data = []
      self.means = []
      self.stddevs = []
      self.clusters = []
   def addDatum(self, values, index):
      new_datum = datum()
      new_datum.index = index #Every data's id is its row number
      new_datum.cl = values[0]
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
   def randomInit(self, numClusters = 10):
   #Modify it to pick far cluster centers?
      for x in range(numClusters):
         index = random.randint(0,len(self.data)-1)
         center = self.data.pop(index)
         self.clusters.append(center)
if __name__ == "__main__":
   if len(sys.argv) != 2:
      print("Error - usage is " + sys.argv[0] + " <data_file>")
      sys.exit(1)
   m = cData()
   m.parseCsv(sys.argv[1])
   m.zvalues()
   m.randomInit()
   iteration = EM(m)
   iteration.EM(3)
   # for i in m.data:
      # print  str(i.index) + " " + str(i.values) + " " + i.cl
         
         
         
         
         
