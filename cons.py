import sys
import random as R
import numpy as np

class cCons:
    def __init__(self):
        self.cons = [] #constraints so far
        self.poscons = [] #possible constraints
        self.consfile = 0
        self.emclusters = []
        # options for constraint generation (set these with setType)
        self.constype = 3 # triplet or pair constraints (3 or 2)
        self.consselect = 1 # 0 for random, 1 for using metric
        
        # keep track of A values selected for triplet constraints
        self.AHistory = []
        
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
    def tripConsmid(self,mGammas,numTrips):
      gammadiffs = self.findDiffs(mGammas)
      constraints = []
      link = 0   
      for i in range(numTrips):
         if(not self.consselect):  # use random (not metric)
            A = R.choice(gammadiffs)
            gammadiffs.remove(A)
         else:
            A = -1
            while A == -1 or A[4] in self.AHistory:
               A = gammadiffs.pop(0)
            self.AHistory.append(A[4])
         class1mids = self.emclusters[A[1]].midpoints[:]
         class2mids = self.emclusters[A[2]].midpoints[:]
         class1mids.append(self.emclusters[A[1]].center)
         class2mids.append(self.emclusters[A[2]].center)
         if self.data[class1mids[-1].index].cl == self.data[A[4]].cl:
            link = 1
            self.emclusters[A[1]].determined.append(A[4])
         elif self.data[class2mids[-1].index].cl == self.data[A[4]].cl:
            link = -1
            self.emclusters[A[1]].determined.append(A[4])
         if(link != 0):
            for i in class1mids:
               constraints.append((A[4],i.cl,link))
            for i in class2mids:
               constraints.append((A[4],i.cl,-1*link))
      return constraints

    def tripConsgamma(self,mGammas,numTrips):
      gammadiffs = self.findDiffs(mGammas)
      constraints = []
      link = 0   
      for i in range(numTrips):
         if(not self.consselect):  # use random (not metric)
            A = R.choice(gammadiffs)
            gammadiffs.remove(A)
         else:
            A = -1
            while A == -1 or A[4] in self.AHistory:
               A = gammadiffs.pop(0)
            self.AHistory.append(A[4])
         class1 = sorted(filter(lambda x: x[1] == A[1],gammadiffs),key=lambda y: y[3])
         class2 = sorted(filter(lambda x: x[1] == A[2],gammadiffs),key=lambda y: y[3])
         class1 = class1[int(np.floor(0.8*len(class1))):]
         class2 = class2[int(np.floor(0.8*len(class2))):]
         if len(class1)>0 and (m.data[A[4]].cl == str(m.data[class1[-1][4]].cl)):
            link = 1
         elif len(class2)>0 and (m.data[A[4]].cl == str(m.data[class2[-1][4]].cl)):
            link = -1        
         if(link != 0):
            for i in class1:
               constraints.append((A[4],i[4],link))
            for i in class2:
               constraints.append((A[4],i[4],-1*link))
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

            
