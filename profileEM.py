import EM
import cData

def run():
    D = cData.cData("data/winenorm3_pyre.csv")
    E = EM.cEM(D)
    E.EM(3)
    
