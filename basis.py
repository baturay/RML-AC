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
        
   # pick A, find B and C and then return a list of corresponding
   # pairwise constraints
   def tripCons(self,mGammas,numTrips):
      gammadiffs = self.findDiffs(mGammas)
      constraints = []
      link = 0   
      for i in range(numTrips):
         if(not self.consselect):  # use random (not metric)
            A = R.choice(gammadiffs)
            gammadiffs.remove(cons)
         else:
            A = gammadiffs.pop(0)
         # class1 = sorted(filter(lambda x: x[1] == cons[1],gammadiffs),key=lambda y: y[3])
         # class2 = sorted(filter(lambda x: x[1] == cons[2],gammadiffs),key=lambda y: y[3])
         # class1 = class1[int(np.floor(0.8*len(class1))):]
         # class2 = class2[int(np.floor(0.8*len(class2))):]
         # if len(class1)>0 and (m.data[cons[4]].cl == str(m.data[class1[-1][4]].cl)):
            # link = 1
         # elif len(class2)>0 and (m.data[cons[4]].cl == str(m.data[class2[-1][4]].cl)):
            # link = -1
         class1mids = self.emclusters[A[1]].midpoints[:]
         class2mids = self.emclusters[A[2]].midpoints[:]
         class1mids.append(self.emclusters[A[1]].center)
         class2mids.append(self.emclusters[A[2]].center)
         if self.data[class1mids[-1].index].cl == self.data[A[4]].cl:
            link = 2
         elif self.data[class2mids[-1].index].cl == self.data[A[4]].cl:
            link = -2           
         if(link != 0):
            for i in class1mids:
               constraints.append((A[4],i.index,link))
            for i in class2mids:
               constraints.append((A[4],i.index,-1*link))
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

   # return a list of distances from each point in the cluster
   # to the closest source
   # **** consider using normalized data when calculating distances
   #      or using a distance metric based on EM parameters
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
      
   # populates self.emclusters with a list of emcluster objects
   # representing the clusters given by the _EM_ argument
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
         
   # updates the emcluster objects in self.emclusters
   # to include the correct center point and the mid-
   # and outer points
   def repPoints(self,EM):     
      print "Classes of the midpoints of clusters: "
      for cl in self.emclusters:
         #Compares every point to the calculated center.
         cdist = sorted(self.findMin([cl.center],cl),key = lambda x : x[1])
         #Finds the closest real point and assigns it.
         cl.center = cdist[0][0]
         #Removes it from the remaining points.
         cl.points.remove(cdist[0][0])        
         #The largest distance from the center is the first outerpoint.
         cl.outerpoints.append(cdist[-1][0])
         cl.points.remove(cdist[-1][0])
         #Other outerpoints are found by finding the maxmin of a point.
         for i in range(5):
            odist = max(self.findMin(cl.outerpoints,cl),key = lambda x : x[1])                     
            cl.outerpoints.append(odist[0])
            cl.points.remove(odist[0])
         
         for o in cl.outerpoints:
            midvalues = []
            # midvalues <- calculated mid point
            for index,value in enumerate(o.values):
               midvalues.append((value+cl.center.values[index])/2)               
            datMidvalues = datum()
            datMidvalues.values = midvalues
            
            #The real point closest to the imaginary midpoint is found and added.
            mdist = min(self.findMin([datMidvalues],cl),key = lambda x : x[1])
            cl.midpoints.append(mdist[0])
            cl.points.remove(mdist[0])
            
         print [self.data[i.index].cl for i in cl.midpoints]
            
         
   def emRestarts(self,k):
      centers = []
      maxEM = cEM(self)
      maxEM.EM(len(self.classlist))
      maxNMI = self.evaluateEM(maxEM)
      for i in range(k):
         iteration = cEM(self)
         iteration.EM(len(self.classlist))
         nmi = self.evaluateEM(iteration)
         if maxNMI < nmi:
            maxNMI = nmi
            maxEM = iteration
      print maxNMI
      return maxEM

   def goodInitial (self,em):
      consistent = 0
      #Consistent means all the midpoints are same with the center.
      constraints = []
      while consistent != len(self.emclusters):
         consistent = 0
         resetCenters = []
         for ind,cl in enumerate(self.emclusters):
            realpoints = [self.data[i.index] for i in cl.midpoints]
            realcenter = self.data[cl.center.index]
            #Gets the real classes so we can do the simulation.
            rightclass = filter(lambda x: x.cl==realcenter.cl,realpoints)
            rightclass.append(realcenter)
            #Filters the midpoints same with center and then adds center.
            wrongclass = filter(lambda x: x.cl!=realcenter.cl,realpoints)
            #All the leftovers...
            if len(wrongclass) == 0:
               consistent += 1
               resetCenters.append(ind)
            #Cross constraints between right and wrong classes.
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
         #If all classes are not right, restart.
         if consistent != len(self.emclusters):    
            em.lInitialCenters = []
            em.EM(len(self.emclusters))
            self.createClusters(em)
            self.repPoints(em)
         Estimated = np.ravel(em.mGammas.argmax(1).T)
         print "nmi: ",nmi(Estimated,self.classes)
      return em   
   def evaluateEM(self,em):
      Estimated = np.ravel(em.mGammas.argmax(1).T)
      return nmi(Estimated,self.classes)
      
   def parseCommandLine (self,argv):
      if len(argv) != 2 and len(argv) != 3:
         print("Error - usage is " + argv[0] + " <data_file> <startingdatapoints>")
         sys.exit(1)   
        
      if len(argv)>2:
         EmAlg = cEM(self)
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
   EmAlg = m.goodInitial(EmAlg)
   #This makes the algorithm start with good initial points.
   
   
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
      EmAlg.EM(len(m.classlist))
      nmiresult = m.evaluateEM(EmAlg)
      print numCons, ",", nmiresult, ",", totalcons 
      if(nmiresult == 1 or len(m.data)==numCons):
         break
   print numCons, ",", nmiresult, ",", totalcons      
         
