import EM
import cData

def run():
    D = cData.cData("data/DS3_samp20.csv")
    E = EM.EM(D)
    E.EM(3)
    
