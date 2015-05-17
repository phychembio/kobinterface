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

prefix="OH-interface/"

q.trajfile=prefix+'short'+args[0]

q.calculate_interface=True


q.cleanevbfilename=prefix+"short.cleanevb"
q.massL=[16,1,16,1]
q.radius=3.5
q.sigmatypeL=[[3.165492,0.0,3.165492,1.605999],[0.0,0.0,1.619017,2.670681],[3.165492,1.619017,3.165492,0.0],[1.605999,2.670681,0.0,0.0]]
q.epstypeL=[[0.1554253,0.0,0.1554253,7.664461],[0.0,0.0,0.094901,0.011165],[0.1554253,0.094901,0.155425,0.0],[7.664461,0.0111653,0.0,0.0]]
q.chargetypeL=[-0.82,0.41,-1.13,0.13]
q.initialzcenterofmass=-1 #find out what this is.
q.HorOH=2 #1=H, 2=OH



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
    
    

printresults(args[0],q)


endtime=datetime.datetime.now()
print endtime-starttime


