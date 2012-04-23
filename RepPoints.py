import numpy as np
from cData import datum

class emcluster:
    def __init__(self):
        self.center = []
        self.index = 0
        self.outerpoints = []
        self.midpoints = []
        self.points = []
        self.determined = []

class RepPoints:
    def __init__(self):
        self.bOptions = True
        self.numMidpoints = 4
        
    # return list of emcluster objects 
    # representing the clusters given by the _EM_ argument
    def createClusters(self,EM):
        D = EM.mData
        # Finds the indices of maximum likelihoods.
        maxindices = np.ravel(EM.mLikelihood_il.argmax(1).T)
        newdatums = []
        for i,cl in enumerate(maxindices):
            newdatum = datum(D.data[i].values)
            newdatum.cl = cl
            newdatum.index = i
            newdatums.append(newdatum)
        # Points turned into datums with Em's guess as their class.
        emclusters = []
        for i,center in enumerate(EM.lCenters):
            c = emcluster()
            c.points = filter(lambda x: x.cl==i,newdatums)
            c.index = i
            #Clusters are formed.
            c.center = datum(center.tolist())
            #Non-point center.
            emclusters.append(c)
        return emclusters
           
    # updates the emcluster objects in self.emclusters
    # to include the correct center point and the mid-
    # and outer points
    def repPoints(self, EM, emclusters):     
        for cl in emclusters:
            # Compares every point to the calculated center.
            cdist = sorted(self.findMin([cl.center],cl),key = lambda x : x[1])
            if (len (cl.points) <= 1): 
               continue
            # Finds the closest real point and assigns it.
            cl.center = cdist[0][0]
            # Removes it from the remaining points.
            cl.points.remove(cdist[0][0])        
            # The largest distance from the center is the first outerpoint.
            cl.outerpoints.append(cdist[-1][0])
            cl.points.remove(cdist[-1][0])
            # Other outerpoints are found by finding the maxmin of a point.
            if(len(cl.points) <= self.numMidpoints):
                cl.midpoints = cl.points[:]
                print [EM.mData.data[i.index].cl for i in cl.midpoints],EM.mData.data[cl.center.index].cl
                continue
            for i in range(self.numMidpoints-1):
                odist = max(self.findMin(cl.outerpoints,cl),key = lambda x : x[1])                     
                cl.outerpoints.append(odist[0])
                cl.points.remove(odist[0])
            if(len(cl.points) <= self.numMidpoints):
                cl.midpoints = cl.points[:]
                print [EM.mData.data[i.index].cl for i in cl.midpoints],EM.mData.data[cl.center.index].cl
                continue          
            for o in cl.outerpoints:
                midvalues = []
                # midvalues <- calculated mid point
                for index,value in enumerate(o.values):
                    midvalues.append((value+cl.center.values[index])/2)               
                datMidvalues = datum()
                datMidvalues.values = midvalues
                # The real point closest to the imaginary midpoint is found and added.
                mdist = min(self.findMin([datMidvalues],cl),key = lambda x : x[1])
                cl.midpoints.append(mdist[0])
                cl.points.remove(mdist[0])
            
            print [EM.mData.data[i.index].cl for i in cl.midpoints],EM.mData.data[cl.center.index].cl

    # return a list of distances from each point in the cluster
    # to the closest source
    # **** consider using normalized data when calculating distances
    #      or using a distance metric based on EM parameters
    def findMin(self,sources,cluster):
        distances = []
        for point in cluster.points:
            mindist = np.inf
            for source in sources:
                # Finds the minimum distance of a point to source points.
                s = source.values
                distance = 0
                for index,value in enumerate(point.values):
                    distance += (value - s[index])**2
                # Just distance calculation.
                if distance < mindist:
                    mindist = distance
            distances.append([point,mindist])
        # The minimum distance of every point in the cluster
        # to the given source points.
        return distances    
