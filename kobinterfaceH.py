#modified on 4/12/2015 
#only for interfaces
import math
import sys
import datetime
from koblib import *
from kobcompute import *
from kobinfo import*

starttime=datetime.datetime.now()
args=sys.argv[1:]
if len(args)==0:
    args.append("")

q=KobInfo()


OWID=1
q.OWID=OWID
HWID=2
q.HWID=HWID
OHID=3
q.OHID=OHID
HHID=4
q.HHID=HHID

prefix="H-interface/short4"

q.skipstep=0
q.trajfile=prefix+''+args[0]

q.calcSurf=False
q.studysurf=False
q.calcCN=True
q.calcenergies=False
q.surflevel=0.015
q.surftol=0.0001
q.cleanevbfilename=prefix+".cleanevb"
q.massL=[16,1,16,1]
q.sigmatypeL=[[3.165492,1.582746,3.142,1.582746],[0.0,0.0,1.6,0.0],[3.142,1.6,3.118508,1.559254],[1.582746,0.0,1.559254,0.0]]
q.epstypeL=[[0.1554253,0.0,0.1238,0.0025076 ],[0.0,0.0,3.0,0.0],[0.1238,3.0,0.098609686,0.0019974],[0.0025076,0.0,0.0019974,0.000040458]]
q.radius=3.2
q.CNinterfacez=13.5
q.prefix="full"
q.suffix=str(q.CNinterfacez)
q.chargetypeL=[-0.82,0.41,-0.32,0.44]
q.HorOH=1 #1=H, 2=OH

#q.debugf=open("zpos.txt","w")

initialize(q)


stepskipped=0 #time steps already skipped, no need to change


trajf=open(q.trajfile+q.trajext,'rU')             
currenttime=trajtime(trajf,q)  
####End initialize#####    

tempstate=inputstate(trajf,q,currenttime)

                   
fixpbc(tempstate,q)
q.statesL.append(tempstate)


while True:
    try:
        currenttime=q.statesL[-1][0]
        q.time=currenttime
        if currenttime%10000==0:
            print currenttime
            sys.stdout.flush()    
        if currenttime%50000==0 and currenttime!=q.statesL[0][0]:
            printresults(args[0],q)
        compute(q) #calculate quantities from statesL        
       

        uppert=upper(currenttime)
        beforesize=len(q.statesL)
        kobclean(q.statesL,currenttime)  #clean up states
        
        
        nexttime=trajtime(trajf,q)
        tempstate=inputstate(trajf,q,nexttime)                                             
        if q.statesL[-1][0]<nexttime:                             
            fixpbc(tempstate,q)
            q.statesL.append(tempstate) 
     
    except StopIteration:           
        break
    
    
endprogram(args[0],q)

endtime=datetime.datetime.now()
print endtime-starttime

