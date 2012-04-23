import sys
import random as R
import pickle
from EM import *
from NMI import *

class datum:
   def __init__(self, values=[]):
      self.values = values
      self.name = ""
      self.cl = ""
      
class cData:
   def __init__(self,filename):
      self.data = []
      self.classlist = dict()
      self.parseCsv(filename)
      self.classes = array([ self.classlist[i.cl] for i in self.data])

   # appends self.data with datum created with _values_
   # and self.classlist with numeric class labels
   def addDatum(self, values, index):
      new_datum = datum()
      new_datum.index = index #Every data's id is its row number
      new_datum.cl = values[0]
      #Classlist keeps all possible classes.
      if(not (new_datum.cl in self.classlist)):
         self.classlist[new_datum.cl] = len(self.classlist)
      new_datum.values = [float(x) for x in values[1:]]
      self.data.append(new_datum)

   # data loaded from _filename_ and written to self.data
   def parseCsv(self,filename):
      with open(filename,"r") as fin:
         lines = fin.readlines()
      for i in range(1, len(lines)):
         values = lines[i].rstrip().split(",")
         self.addDatum(array(values),i-1)
      
