from numpy import *
import sys
if __name__ == "__main__":
   f=open(sys.argv[1])
   g=open("norm"+sys.argv[1],"w")
   lines = f.readlines()
   g.write(lines[0])
   lines.pop(0)
   linevalues = [i.rstrip().split(",")  for i in lines]
   minarray = [inf for i in linevalues[1]]
   maxarray = [-inf for i in linevalues[1]]
   for i in range(len(linevalues)):
      linevalues[i] = [linevalues[i][0]] + [float(s) for s in linevalues[i][1:]] 
      for j in range(1,len(linevalues[i])):
         if linevalues[i][j] < minarray[j]:
           minarray[j] = linevalues[i][j]
         if linevalues[i][j] > maxarray[j]:
           maxarray[j] = linevalues[i][j]
   for i in range(len(maxarray)):
      maxarray[i] -= minarray[i]
       
   for i in linevalues:
       for j in range(1,len(i)):
          i[j] = i[j]-minarray[j]
          if maxarray[j] != 0:
             i[j] /= maxarray[j]
       g.write(",".join(map(str,i))+"\n")
