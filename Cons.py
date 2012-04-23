import sys
import random as R
import numpy as np

class cCons:
    def __init__(self, _D):
        self.D = _D # link to the data
        self.cons = [] #constraints so far
        self.poscons = [] #possible constraints
        self.consfile = 0
        self.emclusters = []
        # options for constraint generation (set these with setType)
        self.constype = cCons.eConsType.TripCenterChunk
        self.centerChunkSize = 0.20
        self.consselect = 1 # 0 for random, 1 for using metric
        self.ConsStrength = 2
        
        # keep track of A values selected for triplet constraints
        self.AHistory = []

    class cMetricData:
        def __init__(self, metric, firstindex, secondindex, firstprob, i):
            self.metric = metric
            self.firstindex = firstindex
            self.secondindex = secondindex
            self.firstprob = firstprob
            self.index = i

    class eConsType:
        TripMids, TripCenterChunk = range(2)
        
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
             if(self.D.data[pair[0]].cl==self.D.data[pair[1]].cl):
                link = self.ConsStrength
             else:
                link = -self.ConsStrength
             pair += (link,)
          self.cons.append(pair)
          cons.append(pair)
       return cons

    # store constraints from _filename_ to
    # self.poscons
    def parseConstraints(self,filename):
        self.consfile = 1
        with open(filename,"r") as fin:
            lines = fin.readlines()
        for i in range(0, len(lines)):        
            values = lines[i].rstrip().split(",")
            self.poscons.append(array([ int(v) for v in values] ))
        
    def tripCons(self, mGammas, numTrips):
        if self.constype == cCons.eConsType.TripMids:
            return self.tripConsMid(mGammas, numTrips)
        elif self.constype == cCons.eConsType.TripCenterChunk:
            return self.tripConsCenterChunk(mGammas,numTrips)

    # pick A, find B and C and then return a list of corresponding
    # pairwise constraints
    def tripConsMid(self,mGammas,numTrips):
      gammadiffs = self.findDiffs(mGammas)
      constraints = []
      link = 0   
      for i in range(numTrips):
         if(not self.consselect):  # use random (not metric)
            A = R.choice(gammadiffs)
            gammadiffs.remove(A)
         else:
            A = -1
            while A == -1 or A.index in self.AHistory:
               A = gammadiffs.pop(0)
            self.AHistory.append(A.index)
         class1mids = self.emclusters[A.firstindex].midpoints[:]
         class2mids = self.emclusters[A.secondindex].midpoints[:]
         class1mids.append(self.emclusters[A.firstindex].center)
         class2mids.append(self.emclusters[A.secondindex].center)
         if self.D.data[class1mids[-1].index].cl == self.D.data[A.index].cl:
            link = self.ConsStrength
            self.emclusters[A.firstindex].determined.append(A.index)
         elif self.D.data[class2mids[-1].index].cl == self.D.data[A.index].cl:
            link = -self.ConsStrength
            self.emclusters[A.firstindex].determined.append(A.index)
         if(link != 0):
            for i in class1mids:
               constraints.append((A.index,i.cl,link))
            for i in class2mids:
               constraints.append((A.index,i.cl,-1*link))
      return constraints

    # returns a list of pairwise con
    def tripConsCenterChunk(self,mGammas,numTrips):
      gammadiffs = self.findDiffs(mGammas)
      constraints = []
      link = 0   
      for i in range(numTrips):
         if(not self.consselect):  # use random (not metric)
            A = R.choice(gammadiffs)
            gammadiffs.remove(A)
         else:
            A = -1
            while A == -1 or A.index in self.AHistory:
               A = gammadiffs.pop(0)
            self.AHistory.append(A.index)
         class1 = sorted(filter(lambda x: x.firstindex == A.firstindex,gammadiffs),
                         key=lambda y: y.firstprob)
         class2 = sorted(filter(lambda x: x.firstindex == A.secondindex,gammadiffs),
                         key=lambda y: y.firstprob)
         class1 = class1[int((1-self.centerChunkSize)*len(class1)):]
         class2 = class2[int((1-self.centerChunkSize)*len(class2)):]
         if len(class1)>0 and (self.D.data[A.index].cl == str(self.D.data[class1[-1].index].cl)):
            link = self.ConsStrength
         elif len(class2)>0 and (self.D.data[A.index].cl == str(self.D.data[class2[-1].index].cl)):
            link = -self.ConsStrength
         if(link != 0):
            for i in class1:
               constraints.append((A.index,i.index,link))
            for i in class2:
               constraints.append((A.index,i.index,-1*link))
      return constraints   

    # from an EM.mGammas matrix (normd likelihood of pt i in clust j)
    # determine for each datum, most (CA) and second-most (CB) likely
    # cluster assignment and evaluate metric of (CA-CB)/(CA+CB)
    # return [ cMetricData ]
    #   with entry for each datum, all sorted by metric
    def findDiffs(self,gammas):
     gammadiffs = []
     maxindices = np.ravel(gammas.argmax(1).T)
     gammas = gammas.tolist()
     for i,gamma in enumerate(gammas):
         secondprob= -np.inf
         firstindex = maxindices[i]
         firstprob = gamma[firstindex]
         for j in range(len(gamma)):
            if (gamma[j] > secondprob) and (j!=firstindex):
               secondindex = j
               secondprob = gamma[j]
         metric = (firstprob-secondprob)/(firstprob+secondprob)
         gammadiffs.append(cCons.cMetricData(metric,firstindex,secondindex,firstprob,i))
     gammadiffs = sorted(gammadiffs,key=lambda x:x.metric)
     return gammadiffs

            
