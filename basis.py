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
      self.emclusters = []

      # options for constraint generation (set these with setType)
      self.constype = 3 # triplet or pair constraints (3 or 2)
      self.consselect = 1 # 0 for random, 1 for using metric

      
   # public use: sets options for choosing constraints
   # * constype as string tells triplet ("3") or pair ("2")
   # * consselect tells 'random' or 'metric' selection or a filename
   # if random pairwise constraints, they are calculated
   #  and stored in self.poscons
   def setType(self, constype, consselect):
       if constype == "3":
          self.constype = 3
          if consselect == "random":
             self.consselect = 0
          elif consselect == "metric":
             self.consselect = 1
          else:
             parseConstraints(consselect)
       elif constype == "2":
          self.constype = 2
          if consselect == "random" or consselect == "metric":
             self.poscons = [(i,j) for i in range(len(self.data)) for j in range(len(self.data))]
          else:
             self.parseConstraints(consselect) 
       else:
          self.constype = 1
       
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

   # return list of _numC_ triplet or pairwise constraint
   # depending on self.constype
   def makeConst(self,numC):
      if self.constype == 2:
         return pairCons(numC)
      elif self.constype == 3:
         return tripCons(numC)
      else:
         return []
         
   # return list of _numC_ pairwise constraints from self.poscons
   # * list is also appended to self.cons
   # * constraint format is (x,y,link as {-1,1})
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
                 
   # data loaded from _filename_ and written to self.data
   def parseCsv(self,filename):
      with open(filename,"r") as fin:
         lines = fin.readlines()
      for i in range(1, len(lines)):
         values = lines[i].rstrip().split(",")
         self.addDatum(array(values),i)
         
   # store constraints from _filename_ to
   # self.poscons
   def parseConstraints(self,filename):
     self.consfile = 1
     with open(filename,"r") as fin:
        lines = fin.readlines()
     for i in range(0, len(lines)):        
        values = lines[i].rstrip().split(",")
        self.poscons.append(array([ int(v) for v in values] ))
        
   def tripCons(self,mGammas,k):
      gammadiffs = self.findDiffs(mGammas)
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
         # class1 = self.emclusters[cons[1]].midpoints[:]
         # class2 = self.emclusters[cons[2]].midpoints[:]
         # class1.append(self.emclusters[cons[1]].center)
         # class2.append(self.emclusters[cons[2]].center)
         # if self.data[class1[-1].name].cl == str(self.data[cons[4]].cl):
            # link = 2
         # elif self.data[class2[-1].name].cl == str(self.data[cons[4]].cl):
            # link = -2           
         if(link != 0):
            for i in class1:
               constraints.append((cons[4],i[4],link))
            for i in class2:
               constraints.append((cons[4],i[4],-1*link))
      return constraints
      
   # from an EM.mGammas matrix (normd likelihood of pt i in clust j)
   # determine for each datum, most (CA) and second-most (CB) likely
   # cluster assignment and evaluate metric of (CA-CB)/(CA+CB)
   # return [ [metric,firstindex,secondindex,firstprob,i] ]
   #   with entry for each datum, all sorted by metric
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
         newdatum.index = i
         newdatums.append(newdatum)
      #Points turned into datums with Em's guess as their class.
      self.emclusters = []
      for i,center in enumerate(EM.lCenters):
         c = emcluster()
         c.points = filter(lambda x: x.cl==i,newdatums)  
         c.index = i
         #Clusters are formed.
         c.center = datum(center.tolist())
         #Non-point center.
         self.emclusters.append(c)
   def repPoints(self,EM):     
      print "Classes of the midpoints of clusters: "
      for cl in self.emclusters:
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
            
         print [self.data[i.index].cl for i in cl.midpoints]
            
         
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
   def goodInitial (self,em):
      consistent = 0
      constraints = []
      while consistent != len(self.emclusters):
         consistent = 0
         for cl in self.emclusters:
            realpoints = [self.data[i.index] for i in cl.midpoints]
            realcenter = self.data[cl.center.index]
            rightclass = filter(lambda x: x.cl==realcenter.cl,realpoints)
            rightclass.append(realcenter)
            wrongclass = filter(lambda x: x.cl!=realcenter.cl,realpoints)
            if len(wrongclass) == 0:
               consistent += 1
            for i in rightclass:
               for j in realpoints:
                  if j in wrongclass:
                    constraints.append([i.index,j.index,-1])
                  elif j!= i:
                     constraints.append([i.index,j.index,1])
            for i in constraints:
              em.mCij[i[0]][i[1]] = i[2]
              em.mCij[i[1]][i[0]] = i[2]  
         print consistent
         if consistent != len(self.emclusters):    
            em.lInitialCenters = []
            em.EM(len(self.emclusters))
            self.createClusters(em)
            self.repPoints(em)
         Estimated = np.ravel(em.mGammas.argmax(1).T)
         print "nmi: ",nmi(Estimated,self.classes)
         
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
         EmAlg = self.emRestarts(1)
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
   #Creates clusters depending on what EM guessed.
   m.repPoints(EmAlg)
   EmAlg.bPPC = True 
   #Finds the outerpoints and the midpoints and assigns them in emclusters.
   m.goodInitial(EmAlg)
   
   
   prevCons = 0
   totalcons = 0
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
         
