import datetime
import re
import sys
import math
import koblib
import vecmath as vec
import numpy as np
import os
import pyximport
import cPickle as pickle
import gzip

if os.name == 'nt':    
    
    mingw_setup_args = { 'options': { 'build_ext': { 'compiler': 'msvc','include_dirs': np.get_include()}}}
    pyximport.install(setup_args=mingw_setup_args)
#else:
#    pyximport.install()

import kobpairwise as kobpair

def upper(s):
    multipleofsmallest = s
    done = False
    counter = 0
    temp = multipleofsmallest
    while not done and s != 0:
        if temp % 10 != 0:
            done = True
        else:
            temp = temp / 10
            counter+=1 
    return multipleofsmallest - temp % 10 * (10 ** counter) + 10 ** (counter + 1)

def kobclean(statesL,currenttime):
    for state in statesL[1:]:
            s = state[0]
            if upper(s) == currenttime:
                statesL.remove(state)                                     


 
def minimagedist(r2,r1,boxlength):
    diff = r2 - r1    
    increment = 0
    if diff / boxlength > 0.5:
        increment = 1
    elif diff / boxlength < -0.5:
        increment = -1      
    return diff - increment * boxlength


def minimage3Dsqdist(r2,r1,boxlengthL):
    sumsq = 0
    
    for i in range(3):
        boxlength = boxlengthL[i]
        increment = 0
        diff = r2[i] - r1[i]
        if diff / boxlength > 0.5:
            increment = 1
        elif diff / boxlength < -0.5:
            increment = -1      
        diff = diff - increment * boxlength
        sumsq+=diff * diff
    return sumsq


def findzcenter2(q):  # using z center of mass
    noofatom = q.statesL[-1][1]
    sum = 0
    masssum = 0
    for i in range(noofatom):
        type = q.statesL[-1][i + 3][1]
        zpos = q.statesL[-1][i + 3][4]+q.statesL[-1][i + 3][-1]*q.zboxlength
        mass = q.massL[type - 1]
        #if abs(zpos - q.initialzcenterofmass) < 50:  #not letting gas molecules affect this
        sum += q.massL[type - 1] * zpos
        masssum += mass
    return sum / masssum



def backupevbstate(q):
    for i in range(4):
        molBID = q.evbmolBL[i]
        tempL = []
        for item in q.moltoIDL[molBID]:
            tempL.append(item[:])
        q.backupmoltoIDL.append([molBID,tempL])
    
    tempL = []
    for item in q.moltoIDL[0]:
           tempL.append(item[:])
    q.backupmoltoIDL.append([0,tempL])
    for i in range(4):
        molBID = q.evbmolBL[i]
        for atomL in q.moltoIDL[molBID][1:]:
            atomID = atomL[0]            
            q.backupidtomolL.append([atomID,q.IDtomolL[atomID]])
    
        
def changeEVBstateH(q,whichstate):
    if q.currentEVBstate == whichstate:
        return
    nobackup = False
    previousstate = 0
    if len(q.backupmoltoIDL) == 0:
        nobackup = True
        previousstate = q.currentEVBstate
    if whichstate > len(q.evbmolBL) - 1:
        sys.exit("whichstate (%d) out of bound\n" % len(q.evbmolBL))

    reactionshell = q.evbshellL[whichstate]
    fromshell = q.evbshellL[q.currentEVBstate]
    
    if not nobackup:
        if (reactionshell == 0 or reactionshell == 1):           
            for item in q.backupmoltoIDL:
                tempL = []
                for subitem in item[1]:
                    tempL.append(subitem[:])
                    q.typeL[subitem[0]] = subitem[1]
                q.moltoIDL[item[0]] = tempL

            for item in q.backupidtomolL:
                q.IDtomolL[item[0]] = item[1]
                

            if reactionshell == 0:
                return
    else: #there is a reaction in this step.
        if reactionshell == 0:
            previousstate = q.pivotstate
            
    if reactionshell == 2:       
        molA = q.evbmolAL[whichstate]
        previousstate = q.evbmolBL.index(molA)         
        changeEVBstateH(q,previousstate)
    elif reactionshell > 2:
        sys.exit("Reaction shell >2!\n")                   
    
    oldwatermolID = q.evbmolBL[whichstate]
    oldwaterOID = q.moltoIDL[oldwatermolID][0][0]  
    oldwaterH1ID = q.moltoIDL[oldwatermolID][1][0]  
    oldwaterH2ID = q.moltoIDL[oldwatermolID][2][0]  
    r_oldwaterO = q.coordsL[oldwaterOID]

  
    oldOstarmolID = q.evbmolBL[previousstate]        
    oldOstarID = q.moltoIDL[oldOstarmolID][0][0] 
    oldH1starID = q.moltoIDL[oldOstarmolID][1][0]           
    oldH2starID = q.moltoIDL[oldOstarmolID][2][0]
    oldH3starID = q.moltoIDL[oldOstarmolID][3][0]
    r_oldH1star = q.coordsL[oldH1starID]
    r_oldH2star = q.coordsL[oldH2starID]
    r_oldH3star = q.coordsL[oldH3starID]
    
    sqdist1 = minimage3Dsqdist(r_oldwaterO,r_oldH1star,q.boxlengthL)
    sqdist2 = minimage3Dsqdist(r_oldwaterO,r_oldH2star,q.boxlengthL)
    sqdist3 = minimage3Dsqdist(r_oldwaterO,r_oldH3star,q.boxlengthL)

    
        
    q.moltoIDL[oldwatermolID][0][1] = 3
    q.moltoIDL[oldwatermolID][1][1] = 4
    q.moltoIDL[oldwatermolID][2][1] = 4
    q.typeL[oldwaterOID] = 3
    q.typeL[oldwaterH1ID] = 4
    q.typeL[oldwaterH2ID] = 4
    
    
    
    q.moltoIDL[oldOstarmolID][0][1] = 1
    q.typeL[oldOstarID] = 1

    minsqdist = min([sqdist1,sqdist2,sqdist3])
    whichH = None
    if minsqdist == sqdist1:
        q.moltoIDL[oldOstarmolID][0][1] = 1
        q.moltoIDL[oldOstarmolID][2][1] = 2 
        q.moltoIDL[oldOstarmolID][3][1] = 2 
        del q.moltoIDL[oldOstarmolID][1]
        q.moltoIDL[oldwatermolID].append([oldH1starID,4])

        q.typeL[oldH2starID] = 2
        q.typeL[oldH3starID] = 2
        
        q.moltoIDL[0][3] = [oldH1starID,4]
        
    elif minsqdist == sqdist2:
        q.moltoIDL[oldOstarmolID][0][1] = 1
        q.moltoIDL[oldOstarmolID][1][1] = 2 
        q.moltoIDL[oldOstarmolID][3][1] = 2 
        del q.moltoIDL[oldOstarmolID][2]
        q.moltoIDL[oldwatermolID].append([oldH2starID,4])

        q.typeL[oldH1starID] = 2
        q.typeL[oldH3starID] = 2
        
        q.moltoIDL[0][3] = [oldH2starID,4]
    else:
        q.moltoIDL[oldOstarmolID][0][1] = 1
        q.moltoIDL[oldOstarmolID][1][1] = 2 
        q.moltoIDL[oldOstarmolID][2][1] = 2 
        del q.moltoIDL[oldOstarmolID][3]
        q.moltoIDL[oldwatermolID].append([oldH3starID,4])

        q.typeL[oldH1starID] = 2
        q.typeL[oldH2starID] = 2
        
        q.moltoIDL[0][3] = [oldH3starID,4]



    q.moltoIDL[0][0] = [oldwaterOID,3]
    q.moltoIDL[0][1] = [oldwaterH1ID,4]
    q.moltoIDL[0][2] = [oldwaterH2ID,4]
    
    

    q.currentEVBstate = whichstate                
                
def changeEVBstateOH(q,whichstate):
    if q.currentEVBstate == whichstate:
        return
    nobackup = False
    previousstate = 0
    if len(q.backupmoltoIDL) == 0:
        nobackup = True
        previousstate = q.currentEVBstate
    if whichstate > len(q.evbmolBL) - 1:
        sys.exit("whichstate (%d) out of bound\n" % len(q.evbmolBL))

    reactionshell = q.evbshellL[whichstate]
    fromshell = q.evbshellL[q.currentEVBstate]
    
    if not nobackup:
        if (reactionshell == 0 or reactionshell == 1):           
            for item in q.backupmoltoIDL:
                tempL = []
                for subitem in item[1]:
                    tempL.append(subitem[:])
                    q.typeL[subitem[0]] = subitem[1]
                q.moltoIDL[item[0]] = tempL

            for item in q.backupidtomolL:
                q.IDtomolL[item[0]] = item[1]
                

            if reactionshell == 0:
                return
    else: #there is a reaction in this step.
        if reactionshell == 0:
            previousstate = q.pivotstate
            
    if reactionshell >= 2:       
        molA = q.evbmolAL[whichstate]
        previousstate = q.evbmolBL.index(molA)         
        changeEVBstateOH(q,previousstate)
    #elif reactionshell > 2:
    #    sys.exit("Reaction shell >2!\n")                   
    
    oldwatermolID = q.evbmolBL[whichstate]
    oldwaterID = q.moltoIDL[oldwatermolID][0][0]
  
    oldwaterH1ID = q.moltoIDL[oldwatermolID][1][0]
    r_oldwaterH1ID = q.coordsL[oldwaterH1ID]
    
    oldwaterH2ID = q.moltoIDL[oldwatermolID][2][0]
    r_oldwaterH2ID = q.coordsL[oldwaterH2ID]
  
    oldOstarmolID = q.evbmolBL[previousstate]    
    oldHstarID = q.moltoIDL[oldOstarmolID][1][0]    
    oldOstarID = q.moltoIDL[oldOstarmolID][0][0]        
    r_oldOstar = q.coordsL[oldOstarID]
    
    sqdist1 = minimage3Dsqdist(r_oldOstar,r_oldwaterH1ID,q.boxlengthL)
    sqdist2 = minimage3Dsqdist(r_oldOstar,r_oldwaterH2ID,q.boxlengthL)

    q.moltoIDL[oldOstarmolID][0][1] = 1
    q.moltoIDL[oldOstarmolID][1][1] = 2
    q.typeL[oldOstarID] = 1
    q.typeL[oldHstarID] = 2
    
    q.moltoIDL[0][0] = [oldwaterID,3]
    q.typeL[oldwaterID] = 3
    if sqdist1 < sqdist2:
        q.moltoIDL[oldOstarmolID].append([oldwaterH1ID,2])
        q.IDtomolL[oldwaterH1ID] = oldOstarmolID        
        del q.moltoIDL[oldwatermolID][1]
        q.moltoIDL[0][1] = [oldwaterH2ID,4]    
        q.typeL[oldwaterH2ID] = 4
    else:
        q.moltoIDL[oldOstarmolID].append([oldwaterH2ID,2])
        q.IDtomolL[oldwaterH2ID] = oldOstarmolID
        del q.moltoIDL[oldwatermolID][2]    
        q.moltoIDL[0][1] = [oldwaterH1ID,4]   
        q.typeL[oldwaterH1ID] = 4 
    q.moltoIDL[oldwatermolID][0][1] = 3    
    q.moltoIDL[oldwatermolID][1][1] = 4    

    q.currentEVBstate = whichstate

def findwatertypes(q,watertypeL):
    intwaternum = 0
    bulkwaternum = 0
    coordwaternum = 0    
    OstarID = q.moltoIDL[0][0][0]    
    Ostarcoords=q.statesL[-1][OstarID+2][q.xl:q.xl+3]
    #Ostarcoords = q.CEC
    waterIDL = [subL[0][0] for subL in q.moltoIDL[1:]]
    radius = q.radius
    zcom = q.zcom
    watertypeL[0] = -1
    CNinterfacez=q.CNinterfacez
    for ID in waterIDL:    
        if ID == OstarID:
            continue
        molID = q.IDtomolL[ID]
        watercoords = q.statesL[-1][ID + 2][q.xl:q.xl + 3]
        #abswaterz=q.statesL[-1][ID + 2][q.xl+2]+q.statesL[-1][ID + 2][-1]*q.zboxlength
        waterz=q.statesL[-1][ID + 2][q.xl+2]
        xdist = minimagedist(Ostarcoords[0],watercoords[0],q.xboxlength)
        ydist = minimagedist(Ostarcoords[1],watercoords[1],q.yboxlength)
        zdist = minimagedist(Ostarcoords[2],watercoords[2],q.zboxlength)        
        sqdist = xdist ** 2 + ydist ** 2 + zdist ** 2
        if sqdist < radius * radius:         
            coordwaternum+=1
            watertypeL[molID] = 1
        else:            
            if math.fabs(waterz - zcom) > CNinterfacez:
                #print waterz - zcom, molID
                intwaternum+=1
                watertypeL[molID] = 2
            else:                
                bulkwaternum+=1
                watertypeL[molID] = 3
    return [coordwaternum,intwaternum,bulkwaternum]

def findsurf(q,coordsL,boxlengthL):
    xmin=ymin=0
    xmax=ymax=q.yboxlength
    xybin=0.5
    zmin=10 
    zmax=20
    zbin=0.1
    Nxbin=int((xmax-xmin)/xybin)
    Nybin=int((ymax-ymin)/xybin)
    Nzbin=int((zmax-zmin)/zbin)
    #gaussmatrix=np.zeros((Nxbin,Nybin,Nzbin), dtype=np.float64)
    surfacepoint=np.zeros((Nxbin*Nybin,3), dtype=np.float64)           

    nptypeL=np.array(q.typeL,dtype=np.intc)
    starttime = datetime.datetime.now()    
    kobpair.calcSurf2(surfacepoint,coordsL,nptypeL,xmin,xmax,ymin,ymax,zmin,zmax,q.zcom,xybin,zbin,boxlengthL,q.surflevel,q.surftol)                        
    pickle.dump("%f\n"%q.time, q.surff)
            
    str=""
    for index,item in enumerate(surfacepoint):
        x=item[0];y=item[1];z=item[2]                
                
        if x!=0 or y!=0 or z!=0:
            str=str+"%f %f %f "%(x,y,z)                    
    str=str+"\n"
    pickle.dump(str, q.surff)
    endtime = datetime.datetime.now()                                
    print endtime - starttime


def findenergies(q,coordsL,boxlengthL,sqdistL):
    changeEVBfct = None
    if q.HorOH == 1:
        changeEVBfct = changeEVBstateH
    else:
        changeEVBfct = changeEVBstateOH

    
    if q.evbtime == None:
        sys.exit("evbtime is none. Cleanevb is probably at the end or not matched") 
    ####setup intial state####
    starmolID=q.IDtomolL[q.moltoIDL[0][0][0]]
    q.currentEVBstate = q.evbmolBL.index(starmolID)        
    q.pivotstate = q.currentEVBstate
    q.backupmoltoIDL = []        
    q.backupidtomolL = []               

    if q.currentEVBstate != 0:
        changeEVBfct(q,0)
    backupevbstate(q)
               
    sqdistL=np.zeros(q.noofatom*q.noofatom)               
    sqdistL=kobpair.sqdistmatrix(sqdistL,coordsL,boxlengthL) 
    chargetypeL = np.array(q.chargetypeL)
    esptypeL = np.array(q.epstypeL)
    sigmatypeL = np.array(q.sigmatypeL)
           
    
    normalization=1./sum([x*x for x in q.evbciL])   


    q.dataoutL.write("%d "%q.time)
    for i in xrange(4):
        changeEVBfct(q,i)        
        waterindexL = []
        starIDL = np.array([subL[0] for subL in q.moltoIDL[0]])
        starindexL = np.array([ID - 1 for ID in starIDL],dtype=np.intc)            
        for mol in q.moltoIDL[1:]:
            tempL = [subL[0]-1 for subL in mol]
            if starIDL[0]-1 not in tempL:
                waterindexL=waterindexL+tempL
        waterindexL=np.array(waterindexL,dtype=np.intc)
        nptypeL=np.array(q.typeL,dtype=np.intc)
        npidtomolL=np.array(q.IDtomolL,dtype=np.intc)
        watertypeL = np.zeros(len(q.moltoIDL),dtype=np.intc) 

        [coordwaternum,intwaternum,bulkwaternum] = findwatertypes(q,watertypeL)                       
        #print coordwaternum,intwaternum,bulkwaternum,q.evbmolBL[0],zOstar,q.CEC[2], zcom
        #[eng_iwc,eng_iwi,eng_iwb] = jitIWenergy(watertypeL,starindexL,waterindexL,np.array(q.typeL),np.array(q.IDtomolL),coordsL,chargetypeL,esptypeL,sigmatypeL,boxlengthL)   
        #[eng_wwc,eng_wwi,eng_wwb]=jitWWenergy(watertypeL,waterindexL,np.array(q.typeL),np.array(q.IDtomolL),coordsL,chargetypeL,esptypeL,sigmatypeL,boxlengthL)
        #[eng_iwc,eng_iwi,eng_iwb] = jitIWenergybymol(watertypeL,starindexL,waterindexL,np.array(q.typeL),np.array(q.IDtomolL),sqdistL,chargetypeL,esptypeL,sigmatypeL,boxlengthL)   
        #[eng_wwc,eng_wwi,eng_wwb]= jitWWenergybymol(watertypeL,waterindexL,np.array(q.typeL),np.array(q.IDtomolL),sqdistL,chargetypeL,esptypeL,sigmatypeL,boxlengthL)
        [eng_iwc,eng_iwi,eng_iwb] = kobpair.fullIWenergy(watertypeL,starindexL,waterindexL,nptypeL,npidtomolL,sqdistL,chargetypeL,esptypeL,sigmatypeL,boxlengthL,q.noofatom)            
        [eng_wwc,eng_wwi,eng_wwb] = kobpair.fullWWenergy(watertypeL,waterindexL,nptypeL,npidtomolL,sqdistL,chargetypeL,esptypeL,sigmatypeL,boxlengthL,q.noofatom)                        
        weight = q.evbciL[i] * q.evbciL[i]*normalization          
        if coordwaternum>0: q.e_iwcsum+=(eng_iwc + q.evbrepL[i]) * weight / coordwaternum
        q.e_iwisum+=eng_iwi * weight / intwaternum
        q.e_iwbsum+=eng_iwb * weight / bulkwaternum
        if coordwaternum>0: q.e_wwcsum+=eng_wwc * weight / coordwaternum
        q.e_wwisum+=eng_wwi * weight / intwaternum
        q.e_wwbsum+=eng_wwb * weight / bulkwaternum   
        q.dataoutL.write("%f %f %f %f %f %f "%(eng_iwc,eng_iwi,eng_iwb,eng_wwc,eng_wwi,eng_wwb))
        q.cwNtotal+=coordwaternum*weight
        q.iwNtotal+=intwaternum*weight
        q.bwNtotal+=bulkwaternum*weight
        q.dataoutL.write("%d %d %d "%(coordwaternum,intwaternum,bulkwaternum))
                
        #q.debugf.write("%d %d\n"%(q.time,bulkwaternum))            
        #waterindexL = np.array(waterindexL)
        #[eng_wwc,eng_wwi,eng_wwb] = jitWWenergybymol(watertypeL,waterindexL,np.array(q.typeL),np.array(q.IDtomolL),coordsL,chargetypeL,esptypeL,sigmatypeL,boxlengthL)  
    q.dataoutL.write("\n")
    q.ecounter+=1

def studysurf(q,boxlengthL):
    binsize=0.5
    binL=q.surfdiffL
    CECcoords=q.CEC[:]
    CECcoords[2]-=q.zcom
    q.time=pickle.load(q.surff)
    surfdata=pickle.load(q.surff)
    info=[float(x) for x in surfdata.split()]
    surfpts=[]    

    for index,i in enumerate(info):        
        if index%3==0:
            x,y,z=info[index],info[index+1],info[index+2]
            surfpts.append([x,y,z])

    npsurfpts=np.array(surfpts)
    Ndata=len(npsurfpts)
    for i in xrange(Ndata):
        coords=npsurfpts[i]
        sqdist=minimage3Dsqdist(CECcoords,coords,boxlengthL)
        whichbin=int(math.sqrt(sqdist)/binsize)        
        maxbin=len(binL)-1
        finalbin=whichbin
        if whichbin>maxbin:                      
            while maxbin<whichbin: 
                if len(binL)==0:
                    binL.append([binsize/2,[0,0,0]])
                else:
                    binL.append([binL[-1][0]+binsize,[0,0,0]])
                maxbin+=1
            finalbin=-1
        #print i,finalbin,len(npsurfpts),len(binL)   
        bin= binL[finalbin][1]   
        bin[0]+=1        
        m_prev = bin[2]
        x_i=npsurfpts[i][2]-CECcoords[2]

        bin[2] += (x_i - bin[2]) / bin[0]
        bin[1] += (x_i - bin[2]) * (x_i - m_prev)
        

def studysurf2(q,boxlengthL):
    binsize=0.5
    binL=q.surfdiffL
    CECcoords=q.CEC[:]
    CECcoords[2]-=q.zcom
    q.time=pickle.load(q.surff)
    surfdata=pickle.load(q.surff)
    info=[float(x) for x in surfdata.split()]
    surfpts=[]    

    for index,i in enumerate(info):        
        if index%3==0:
            x,y,z=info[index],info[index+1],info[index+2]
            surfpts.append([x,y,z])

    npsurfpts=np.array(surfpts)
    Ndata=len(npsurfpts)
    for i in xrange(Ndata):
        coords=[npsurfpts[i][0],npsurfpts[i][1],CECcoords[2]]
        sqdist=minimage3Dsqdist(CECcoords,coords,boxlengthL)
        whichbin=int(math.sqrt(sqdist)/binsize)        
        maxbin=len(binL)-1
        finalbin=whichbin
        if whichbin>maxbin:                      
            while maxbin<whichbin: 
                if len(binL)==0:
                    binL.append([binsize/2,[0,0,0]])
                else:
                    binL.append([binL[-1][0]+binsize,[0,0,0]])
                maxbin+=1
            finalbin=-1
        #print i,finalbin,len(npsurfpts),len(binL)   
        bin= binL[finalbin][1]   
        bin[0]+=1        
        m_prev = bin[2]
        x_i=npsurfpts[i][2]-CECcoords[2]

        bin[2] += (x_i - bin[2]) / bin[0]
        bin[1] += (x_i - bin[2]) * (x_i - m_prev)        
        

    #sys.exit(0)
def findCN(q,coordsL,boxlengthL,binL):
    nptypeL=np.array(q.typeL,dtype=np.intc)
    sqradius=q.radius*q.radius
    binsize=0.5
    
    OL=np.sort(np.array(q.Typemap[q.OWID]+q.Typemap[q.OHID],dtype=np.intc))
    N=len(OL)
    OIDL=np.array(OL)-2
    sqOOdistL=np.zeros((N,N),dtype=np.float64)             
    sqOOdistL=kobpair.sqdistOOmatrix(sqOOdistL,coordsL,OIDL,boxlengthL) 
    
    CECcoords=q.CEC[:]    
    for i in xrange(N):
        CNcount=0
        sortedOindexL=np.argsort(sqOOdistL[i],kind='quicksort')        
        for j in xrange(1,10):            
            if sqOOdistL[i][sortedOindexL[j]]<sqradius:
                #print sqOOdistL[i][sortedOindexL[j]]
                CNcount+=1
            else:
                break
        CECOdist=math.sqrt(minimage3Dsqdist(CECcoords,coordsL[OIDL[i]],boxlengthL))
        if CECOdist<15.0:
            whichbin=int(CECOdist/binsize)        
            maxbin=len(binL)-1
            finalbin=whichbin
            if whichbin>maxbin:                      
                while maxbin<whichbin: 
                    if len(binL)==0:
                        binL.append([binsize/2,[0,0,0]])
                    else:
                        binL.append([binL[-1][0]+binsize,[0,0,0]])
                    maxbin+=1
                finalbin=-1
            #print i,finalbin,len(npsurfpts),len(binL)   
            bin= binL[finalbin][1]   
            bin[0]+=1        
            m_prev = bin[2]
            x_i=CNcount
            bin[2] += (x_i - bin[2]) *1./ bin[0]
            bin[1] += (x_i - bin[2]) * (x_i - m_prev) 
        

                
def compute(q):                 
        print "time= ",q.time    
        starttime = datetime.datetime.now() 
        q.zcom = findzcenter2(q)
        coordsL = np.array(q.coordsL)
        boxlengthL = np.array(q.boxlengthL)
        q.evbtime = koblib.readcleanevb(q,q.time,q.cleanevbf)          
        if q.calcSurf: findsurf(q,coordsL,boxlengthL)            

        if q.studysurf: 
            studysurf2(q,boxlengthL)
        
        sqdistL=None            
        if q.calcenergies: findenergies(q,coordsL,boxlengthL,sqdistL)
        if q.calcCN: findCN(q,coordsL,boxlengthL,q.OCN_L)


        endtime = datetime.datetime.now()        
        print endtime - starttime
        
        

    




