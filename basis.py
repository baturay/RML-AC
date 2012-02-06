import sys
import numpy
import scipy
import random
from numpy import *
from scipy import *
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
   def addDatum(self, values):
      new_datum = datum()
      new_datum.name = values[-1]
      new_datum.cl = values[-2]
      new_datum.values = [float(x) for x in values[:-2]]
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
   def parseData(self, filename):
      with open(filename,"r") as fin:
         lines = fin.readlines()
      for i in range(len(lines)):
         parts = lines[i].rstrip().split(" ")
         if parts[0] == "@data":
            data_start = i+1
            break
      for i in range(data_start, len(lines)):
         values = lines[i].rstrip().split(",")
         self.addDatum(array(values))
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
   m = machine()
   m.parseData(sys.argv[1])
   m.zvalues()
   m.randomInit()
   for i in m.data:
      print  i.name + " " + str(i.values) + " " + i.cl
         
         
         
         
         
