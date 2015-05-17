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
    
    MSDtypes,FPTpairL,bigFPTL,bigtrackingFPTL, angletypes,pairL,CECorderL,bigMSDL,bigangleL,bigRDFL,watersizeclustersL,nowaterclustersL,statesL, Typemap = ([] for i in range(14))
 
# essential lists 
   #statesL
   #MSDtypes
   #pairL
   #CECorderL
   #bigMSDL
   #bigRDFL
   #Typemap
#not essential lists  
    #respairL,closestL,conditionL,notpolymerL, atommapC1S,OH_L,counthopsL,
    #bigClosestL,bigConditionalL,bigtrackingresL,bigresL,hopsL,C1SL,statesL,waterorderL,molL, bigciL,NstateL