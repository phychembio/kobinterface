import math

class KobInfo:
    
    stepskipped=0    
    smalleststep=0
    noofatom=0
    nooftype=0                
    binsize=0.05
    skipstep=0
    trajext='.lammpstrj'    
    xl=2    
    statesL, Typemap = ([] for i in range(2))
