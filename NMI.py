import numpy
import scipy
from numpy import *
from scipy import *

#Mutual information#Calculate mutual informaion of real num vector x and y
def mi (x,y):
    N=double(x.size)
    I=0.0
    eps = numpy.finfo(float).eps
    for l1 in unique(x):
        for l2 in unique(y):
            #Find intersections
            l1_ids=nonzero(x==l1)[0]
            l2_ids=nonzero(y==l2)[0]
            pxy=(double(intersect1d(l1_ids,l2_ids).size)/N)
            if pxy>0:
                I+=pxy*log2(pxy/((double(l1_ids.size)/N)*(double(l2_ids.size)/N)))
    return I

#Normalized mutual information
def nmi (x,y):#calculate NMI value of vector x and y
    N=x.size
    I=mi(x,y)
    Hx=0.0
    for l1 in unique(x):
        l1_count=nonzero(x==l1)[0].size
        Hx+=-(double(l1_count)/N)*log2(double(l1_count)/N)
        
    Hy=0.0
    for l2 in unique(y):
        l2_count=nonzero(y==l2)[0].size
        Hy+=-(double(l2_count)/N)*log2(double(l2_count)/N)
            
    return 2*I/(Hx+Hy)
        
if __name__=="__main__":
    print nmi(array([1,2,2]),array([1,1,2]))
