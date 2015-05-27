import re
import sys
import math
import kobcompute
import gzip
import numpy as np

def initialize(q):
    gettrajinfo(q)
    q.cleanevbf=open(q.cleanevbfilename,"r")
    q.moltoIDL=[]
    q.IDtomolL=[-1] #molID,row index in q.statesL    
    if q.calcenergies:
        q.dataoutL=open(q.prefix+"infoperframe"+q.suffix+".txt","w")
        q.e_iwcsum=0
        q.e_iwisum=0
        q.e_iwbsum=0
        q.e_wwcsum=0
        q.e_wwisum=0
        q.e_wwbsum=0
        q.cwNtotal=0
        q.iwNtotal=0
        q.bwNtotal=0
        q.ecounter=0
    if q.studysurf:
        q.calcSurf=False
        q.surff=gzip.open("surf.pickle","rb")
        q.surfdiffL=[]
        
    if q.calcSurf:
        q.surff=gzip.open("surf.pickle","wb", compresslevel=5)
    if q.calcCN:
        q.OCN_L=[]
    if q.calcAngle:        
        q.anglebinsize=1.0
        q.angleL=[]
        Nbin=int(180./q.anglebinsize)
        for i in range(Nbin):
            q.angleL.append([q.anglebinsize/2+i*q.anglebinsize,0])
    if q.calcGDS:
        q.GDSbinL=[0]
            
    

def readcleanevb(q,time,cleanevbf):
    evbtime=None
    for line in cleanevbf:
        match = re.match("--",line)
        if match:
            evbtime = int(cleanevbf.next())
            if time == evbtime:                
                q.evbciL = [float(x) for x in cleanevbf.next().split()]
                q.evbmolBL = [int(x) for x in cleanevbf.next().split()]
                q.evbmolAL = [int(x) for x in cleanevbf.next().split()]
                q.evbshellL = [int(x) for x in cleanevbf.next().split()]
                q.evbrepL = [float(x) for x in cleanevbf.next().split()]
                q.evboffdiagL = [float(x) for x in cleanevbf.next().split()]
                q.evbcijL = [float(x) for x in cleanevbf.next().split()]
                q.CEC = [float(x) for x in cleanevbf.next().split()]
                return evbtime
    sys.exit("At the end of cleanevb. Quitting")
     

def gettrajinfo(q):
    trajfile=q.trajfile+q.trajext
    trajf=open(trajfile,"r")
    line=gotoline(trajf,'ITEM: TIMESTEP\n',0)
    time1=int(trajf.next())  
    
    trajf.next()
    
    q.noofatom=int(trajf.next())
    
    trajf.next()
    info=trajf.next().split()
    q.xlo=float(info[0])
    q.xhi=float(info[1])
    info=trajf.next().split()
    q.ylo=float(info[0])
    q.yhi=float(info[1])
    info=trajf.next().split()
    q.zlo=float(info[0])
    q.zhi=float(info[1])
    q.boxlengthL=[]
    q.boxlengthL.append(q.xhi-q.xlo)    
    q.boxlengthL.append(q.yhi-q.ylo)
    q.boxlengthL.append(q.zhi-q.zlo)
    q.xboxlength=q.boxlengthL[0]
    q.yboxlength=q.boxlengthL[1]
    q.zboxlength=q.boxlengthL[2]
    q.boxlength=q.xhi-q.xlo
    labelL=trajf.next().split()[2:] #ITEM: ATOMS id type....
    typelabelindex=None
    for i in range(len(labelL)):
        label=labelL[i]
        if label=="type":
            typelabelindex=i
            break
    maxtypeid=0
    for line in range(q.noofatom):
        line=trajf.next()        
        typeid=int(line.split()[typelabelindex])
        if typeid>maxtypeid:
            maxtypeid=typeid    
    q.nooftype=maxtypeid
    line=gotoline(trajf,'ITEM: TIMESTEP\n',0)
    time2=int(trajf.next())
    trajf.close()    
    q.smalleststep=time2-time1
  
    


def histo(binL,binsize,whichbin,increment):
     maxbin=len(binL)-1 
     finalbin=whichbin
     if whichbin>maxbin:                      
        while maxbin<whichbin: 
            if len(binL)==1:
                binL.append([binsize/2,0,0])
            else:
                binL.append([binL[-1][0]+binsize,0,0])
            maxbin+=1
        finalbin=-1
     binL[finalbin][1]+=increment #last entry
     binL[finalbin][2]+=1
     
def setupTypemap(state,q):    
    q.Typemap=[]
    Typemap=q.Typemap
    for i in range(q.nooftype+1): #+1 is for the spaceholder to match the atomtype IDs
        Typemap.append([])
    index=3
    for item in state[3:]:    
        typeID=item[1]
        atomID=item[0]
        Typemap[typeID].append(index)       
        index+=1



def gotoline(fstream, text, after):
    
    for line in fstream:
        match=re.match(text,line)
        if match:
           found= True
           while after > 0:        
              line=fstream.next()
              after=after -1
           return line
    return 'String Not Found'   




def inputstate(trajf,q,timestep):
    q.typeL=[-1]    
    q.coordsL=[[-999,-999,-999]]
    state=[]  #(timestep, q.noofatom, halfwidth, (atom_1)...(atom_N)(CEC_1)..(CEC_N))
    #atom/CEC info=(atomid, atom type, x, y, z)
    
    state.append(timestep)
    gotoline(trajf,'ITEM: NUMBER',0)
    
    state.append(q.noofatom)
    gotoline(trajf,'ITEM: BOX',0)
    
    state.append("placeholder") #will be used to store Typemap
    line=gotoline(trajf,'ITEM: ATOMS',0)    
    labelL=line.split()[2:]
    xorxs=None
    xindex=None
    molindex=None
    idindex=None
    typeindex=None
    for i in range(len(labelL)):
        label=labelL[i]
        if label=="x":            
            xorxs=label    
            xindex=i        
        if label=="id":
            idindex=i
        if label=="mol":
            molindex=i
        if label=="type":
            typeindex=i            
    tempL=[]
    for i in range(0,q.noofatom):
        atomdata=trajf.next().split()                       
        x=float(atomdata[xindex])
        y=float(atomdata[xindex+1])
        z=float(atomdata[xindex+2])
        atomID=int(atomdata[idindex])
        atomType=int(atomdata[typeindex])
        tempL.append((atomID,atomType,x,y,z))
        #q.coordsL.append([x,y,z])
        if i==0:
            q.moltoIDL=[]            
        molID=int(atomdata[molindex])

        if len(q.moltoIDL)<= molID:
            while len(q.moltoIDL)<=molID:
                q.moltoIDL.append([])    
        if atomType==q.OHID or atomType==q.OWID:           
            q.moltoIDL[molID].insert(0,[atomID,atomType])
        else:
            q.moltoIDL[molID].append([atomID,atomType])

        if atomType==q.OHID:
            q.moltoIDL[0].insert(0,[atomID,atomType])
        elif atomType==q.HHID:
            q.moltoIDL[0].append([atomID,atomType])
        if len(q.IDtomolL)<= atomID:
            while len(q.IDtomolL)<=atomID:
                q.IDtomolL.append([-1,-1]) 
        q.IDtomolL[atomID]=molID

    for i in sorted(tempL,key=lambda k: k[0]):
        state.append(i)       
        q.typeL.append(i[1])   
        q.coordsL.append(i[2:5])
    setupTypemap(state,q)   
    state[2]=q.Typemap



    return state

def fixpbc(tempstate, q):
      for i in range(3,len(tempstate)):
               if len(q.statesL)==0:
                  tempstate[i]=tempstate[i]+(0,0,0) 
               else:                            
                  crossxdist=tempstate[i][2]-q.statesL[-1][i][2]
                  crossydist=tempstate[i][3]-q.statesL[-1][i][3]
                  crosszdist=tempstate[i][4]-q.statesL[-1][i][4]
                  crossx=q.statesL[-1][i][5]
                  crossy=q.statesL[-1][i][6]                  
                  crossz=q.statesL[-1][i][7]                  
                  if abs(crossxdist) >= q.xboxlength/2:
                     if crossxdist< 0:
                        crossx = crossx + 1
                     else:
                        crossx = crossx - 1                  
                  if abs(crossydist) >= q.yboxlength/2:
                     if crossydist< 0:
                        crossy = crossy + 1
                     else:
                        crossy = crossy - 1
                  if abs(crosszdist) >= q.zboxlength/2:
                     if crosszdist< 0:
                        crossz = crossz + 1
                     else:
                        crossz = crossz - 1            		                                     		                                     		                                     		             
                  tempstate[i]=tempstate[i]+(crossx,crossy,crossz)     

   

def trajtime(trajf,q):                       
    line=gotoline(trajf,'ITEM: TIMESTEP\n',0)
    time=int(trajf.next()) 
    
    while q.stepskipped<q.skipstep:
        line=gotoline(trajf,'ITEM: TIMESTEP\n',0);
        time=int(trajf.next()) ;   
        q.stepskipped+=1

    if line!="String Not Found":     
        return time 
    else:
        print "End of file: trajf"
        while True:
            trajf.next()
        


def countatoms(map, atomIDlist):
    sum=0
    for atomID in atomIDlist:
        sum+=len(map[atomID])
    return sum

def endprogram(suffix,q):
    printresults(suffix,q)    
    if q.calcenergies:
        q.dataoutL.close()
    if q.calcSurf:
        q.surff.close()
    #q.debugf.close()


def printresults(suffix,q):
    #volume=q.xboxlength*q.yboxlength*q.zboxlength
    if q.calcenergies:
        inf=open(q.prefix+"localenergies"+q.suffix+".txt","w")
        inf.write("%f\n"%(q.e_iwcsum/q.ecounter))
        inf.write("%f\n"%(q.e_iwisum/q.ecounter))
        inf.write("%f\n"%(q.e_iwbsum/q.ecounter))
        inf.write("%f\n"%(q.e_wwcsum/q.ecounter))
        inf.write("%f\n"%(q.e_wwisum/q.ecounter))
        inf.write("%f\n"%(q.e_wwbsum/q.ecounter))
        inf.write("%f\n"%(q.cwNtotal*1./q.ecounter))
        inf.write("%f\n"%(q.iwNtotal*1./q.ecounter))
        inf.write("%f\n"%(q.bwNtotal*1./q.ecounter))
        inf.close()
    if q.studysurf:
        outf=open("surflateralstats.txt","w")
        #for bin in q.surfdiffL:
        #    if len(bin[1])!=0:
        #        outf.write("%f %f\n"%(bin[0], np.std(bin[1],ddof=1)))
        for bin in q.surfdiffL:
            if bin[1][0]>=2:
                outf.write("%f %f\n"%(bin[0], bin[1][1]*1./(bin[1][0]-1.0)))   
        outf.close()
    if q.calcCN:
        outf=open("OCNstats3.2.txt","w")
        #for bin in q.surfdiffL:
        #    if len(bin[1])!=0:
        #        outf.write("%f %f\n"%(bin[0], np.std(bin[1],ddof=1)))
        for bin in q.OCN_L:
            if bin[1][0]>=2:
                outf.write("%f %f %f\n"%(bin[0], bin[1][2],bin[1][1]*1./(bin[1][0]-1.0)))  
        outf.close()
    if q.calcAngle:
        outf=open("anglestats.txt","w")
        for bin in q.angleL:
             outf.write("%f %d\n"%(bin[0], bin[1]))
        outf.close()
    
    
    outf=open('GDSdensity.txt','w')
    for bin in q.GDSbinL[1:]:
       outf.write("%f  %f \n" %(bin[0] ,bin[1]*1./q.GDSbinL[0]))
    outf.close()
    
