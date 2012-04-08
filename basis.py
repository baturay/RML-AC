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
      
class emcluster:
   def __init__(self):
      self.center = []
      self.index = 0
      self.outerpoints = []
      self.midpoints = []
      self.points = []
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
      self.consselect = 0
      self.clusters = []
      
   def setType(self, constype, consselect):
       if constype == "3":
          m.constype = 3
          if consselect == "random":
             m.consselect = 0
          elif consselect == "metric":
             m.consselect = 1
          else:
             parseConstraints(consselect)
       elif constype == "2":
          m.constype = 2
          if consselect == "random" or consselect == "metric":
             m.poscons = [(i,j) for i in range(len(m.data)) for j in range(len(m.data))]
          else:
             m.parseConstraints(consselect) 
       else:
          m.constype = 1
       
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
      cons = []
      for i in range(numC):
         pair = R.choice(self.poscons)
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
        
   def tripCons(self,mGammas,k):
      gammadiffs = m.findDiffs(mGammas)
      constraints = []
      link = 0   
      for i in range(k):
         if(not self.consselect):
            cons = R.choice(gammadiffs)
            gammadiffs.remove(cons)
         else:
            cons = gammadiffs.pop(0)
         class1 = sorted(filter(lambda x: x[1] == cons[1],gammadiffs),key=lambda y: y[3])
         class2 = sorted(filter(lambda x: x[1] == cons[2],gammadiffs),key=lambda y: y[3])
         class1 = class1[int(np.floor(0.8*len(class1))):]
         class2 = class2[int(np.floor(0.8*len(class2))):]
         if len(class1)>0 and (m.data[cons[4]].cl == str(m.data[class1[-1][4]].cl)):
            link = 1
         elif len(class2)>0 and (m.data[cons[4]].cl == str(m.data[class2[-1][4]].cl)):
            link = -1
         print link
         if(link != 0):
            for i in class1:
               constraints.append((cons[4],i[4],link))
            for i in class2:
               constraints.append((cons[4],i[4],-1*link))
      return constraints
      
   def findDiffs(self,gammas):
     gammadiffs = []
     maxindices = np.ravel(gammas.argmax(1).T)
     gammas = gammas.tolist()
     for i,gamma in enumerate(gammas):
         secondprob= -inf
         firstindex = maxindices[i]
         firstprob = gamma[firstindex]
         for j in range(len(gamma)):
            if (gamma[j] > secondprob) and (j!=firstindex):
               secondindex = j
               secondprob = gamma[j]
         metric = (firstprob-secondprob)/(firstprob+secondprob)
         gammadiffs.append([metric,firstindex,secondindex,firstprob,i])
     gammadiffs = sorted(gammadiffs,key=lambda x:x[0])
     return gammadiffs
     
   def findMin(self,sources,cluster):
      distances = []
      for point in cluster.points:
         mindist = inf
         for source in sources:
         #Finds the minimum distance of a point to source points.
            s = source.values
            distance = 0
            for index,value in enumerate(point.values):
               distance += (value - s[index])**2
            #Just distance calculation.
            if distance < mindist:
               mindist = distance
         distances.append([point,mindist])
      #The minimum distance of every point in the cluster
      #to the given source points.
      return distances    
   def createClusters(self,EM):
      maxindices = np.ravel(EM.mLikelihood_il.argmax(1).T)
      #Finds the indices of maximum likelihoods.
      newdatums = []
      for i,cl in enumerate(maxindices):
         newdatum = datum(self.data[i].values)
         newdatum.cl = cl
         newdatum.name = i
         newdatums.append(newdatum)
      #Points turned into datums with Em's guess as their class.
      clusters = []
      for i,center in enumerate(EM.lCenters):
         c = emcluster()
         c.points = filter(lambda x: x.cl==i,newdatums)  
         c.index = i
         #Clusters are formed.
         c.center = datum(center.tolist())
         #Non-point center.
         self.clusters.append(c)
   def repPoints(self,EM):     
      print "Classes of the midpoints of clusters: "
      for cl in self.clusters:
         cdist = sorted(self.findMin([cl.center],cl),key = lambda x : x[1])
         #Compares every point to the calculated center.
         cl.center = cdist[0][0]
         #Finds the closest real point and assigns it.
         cl.points.remove(cdist[0][0])        
         #Removes it from the remaining points.
         cl.outerpoints.append(cdist[-1][0])
         #The largest distance from the center is the first outerpoint.
         cl.points.remove(cdist[-1][0])
         for i in range(5):
            #Other outerpoints are found by finding the maxmin of a point.
            odist = sorted(self.findMin(cl.outerpoints,cl),key = lambda x : x[1])                     
            cl.outerpoints.append(odist[-1][0])
            cl.points.remove(odist[-1][0])
         
         for o in cl.outerpoints:
            midvalues = []
            for index,value in enumerate(o.values):
               midvalues.append((value+cl.center.values[index])/2)
               #A point equidistant to the center and an outerpoint.
            npmidpoint = datum()
            npmidpoint.values = midvalues
            #An imaginary midpoint is found.
            mdist = sorted(self.findMin([npmidpoint],cl),key = lambda x : x[1])
            #The real point closest to the imaginary midpoint is found and added.
            cl.midpoints.append(mdist[0][0])
            cl.points.remove(mdist[0][0])
            
         print [i.cl for i in cl.midpoints]
            
         
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
      print maxNMI
      return maxEM
      
   def evaluateEM(self,em):
      em.EM(len(self.classlist))
      Estimated = np.ravel(em.mGammas.argmax(1).T)
      return nmi(Estimated,self.classes)
      
   def parseCommandLine (self,argv):
      if len(argv) != 2 and len(argv) != 3:
         print("Error - usage is " + argv[0] + " <data_file> <startingdatapoints>")
         sys.exit(1)   
        
      if len(argv)>2:
         EmAlg = EM(self)
         f = open(argv[2],"r")
         centers = pickle.load(f)
         EmAlg.lInitialCenters = centers
         EmAlg.EM(len(self.classlist))
      else:
         EmAlg = self.emRestarts(5)
         f = open("centers.pickle","w")
         l = EmAlg.lCenters
         pickle.dump(l,f)
      return EmAlg
      
if __name__ == "__main__":
   m = cData(sys.argv[1])
   #Takes in the file and parses into datum's.
   m.setType("3","metric")
   #Uses triple constraints and the metric to choose.
   EmAlg = m.parseCommandLine(sys.argv) 
   #Starting points from the pickle or not.
   m.createClusters(EmAlg)
   m.repPoints(EmAlg)
   
   EmAlg.bPPC = True
   
   prevCons = 0
   totalcons = 0
   initial = EmAlg.lInitialCenters
   nmiResult = m.evaluateEM(EmAlg)
   print "Initial nmi: ",nmiResult
   for numCons in range(1,len(m.data),1):     
      cons = m.tripCons(EmAlg.mGammas,numCons-prevCons)
      prevCons = numCons
      totalcons += len(cons)
      for i in cons:
        EmAlg.mCij[i[0]][i[1]] = i[2]
        EmAlg.mCij[i[1]][i[0]] = i[2]    
      nmiresult = m.evaluateEM(EmAlg)
      print numCons, ",", nmiresult, ",", totalcons 
      if(nmiresult == 1 or len(m.data)==numCons):
         break
   print numCons, ",", nmiresult, ",", totalcons      
         