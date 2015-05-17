import datetime
import re
import sys
import math
import koblib
import vecmath as vec
import numba as nb
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

@nb.jit 
def jitIWenergy(watertypeL,starindexL, waterindexL,typeL,IDtomolL,coordsL,chargetypeL,epstypeL,sigmatypeL,boxlengthL): 
    
    eng_iwc = 0
    eng_iwi = 0
    eng_iwb = 0
    noofatom = len(coordsL) - 1    
    for i in starindexL:        
        typei = typeL[i + 1]
        molIDi = IDtomolL[i + 1]
        coordsi = coordsL[i + 1]
        chargei = chargetypeL[typei - 1]                
        epsL = epstypeL[typei - 1]
        sigmaL = sigmatypeL[typei - 1]
        for j in waterindexL:          
               typej = typeL[j + 1]
               eps = epsL[typej - 1]
               sigma = sigmaL[typej - 1]
               molIDj = IDtomolL[j + 1]
               watertype = watertypeL[molIDj]               
               coordsj = coordsL[j + 1]
               chargej = chargetypeL[typej - 1]                       
               if molIDi != molIDj:                       
                    sqdistij = 0                                      
                    for k in range(3):
                        increment = 0                            
                        diff = coordsi[k] - coordsj[k]
                        if diff / boxlengthL[k] > 0.5:increment = 1
                        elif diff / boxlengthL[k] < -0.5:increment = -1      
                        diff = diff - increment * boxlengthL[k]                            
                        sqdistij+=diff * diff                         
                    distij = np.sqrt(sqdistij)                   
                    xsq = sigma * sigma / sqdistij
                    x6p = xsq * xsq * xsq
                    e1 = 4 * eps * (x6p * x6p - x6p)                    
                    e2 = 332.0636 * chargei * chargej / distij

                    if watertype == 1:
                        eng_iwc+=e1 + e2                          
                    elif watertype == 2:
                        eng_iwi+=e1 + e2                           
                    elif watertype == 3:                        
                        eng_iwb+=e1 + e2                        


    return eng_iwc,eng_iwi,eng_iwb


@nb.jit
def jitWWenergy(watertypeL,waterindexL,typeL,IDtomolL,coordsL,chargetypeL,epstypeL,sigmatypeL,boxlengthL):  
    
    eng_wwc = 0
    eng_wwi = 0    
    eng_wwb = 0
    
    noofatom = len(coordsL) - 1    
    for i in waterindexL:
        typei = typeL[i + 1]
        molIDi = IDtomolL[i + 1]
        coordsi = coordsL[i + 1]
        chargei = chargetypeL[typei - 1]                
        epsL = epstypeL[typei - 1]
        sigmaL = sigmatypeL[typei - 1]
        watertype1 = watertypeL[molIDi]
        for j in waterindexL:          
               typej = typeL[j + 1]
               eps = epsL[typej - 1]
               sigma = sigmaL[typej - 1]
               molIDj = IDtomolL[j + 1]
               watertype2 = watertypeL[molIDj]
               if j > i:
                   coordsj = coordsL[j + 1]
                   chargej = chargetypeL[typej - 1]                       
                   if molIDi != molIDj:                       
                        sqdistij = 0                                      
                        for k in range(3):
                            increment = 0                            
                            diff = coordsi[k] - coordsj[k]
                            if diff / boxlengthL[k] > 0.5:increment = 1
                            elif diff / boxlengthL[k] < -0.5:increment = -1      
                            diff = diff - increment * boxlengthL[k]                            
                            sqdistij+=diff * diff                         
                        distij = np.sqrt(sqdistij)                   
                        xsq = sigma * sigma / sqdistij
                        x6p = xsq * xsq * xsq
                        e1 = 4 * eps * (x6p * x6p - x6p)
                        e2 = 332.0636 * chargei * chargej / distij

                        if watertype1 == 1:                            
                            eng_wwc+=e1 + e2                            
                        elif watertype1 == 2:                            
                            eng_wwi+=e1 + e2                            
                        elif watertype1 == 3:                            
                            eng_wwb+=e1 + e2
                            
                        if watertype2 == 1:                            
                            eng_wwc+=e1 + e2                            
                        elif watertype2 == 2:                            
                            eng_wwi+=e1 + e2                            
                        elif watertype2 == 3:                            
                            eng_wwb+=e1 + e2                                    


    return eng_wwc,eng_wwi,eng_wwb

@nb.jit(nopython=True) 
def jitIWenergybymol(watertypeL,starindexL, waterindexL,typeL,IDtomolL,sqdistL,chargetypeL,epstypeL,sigmatypeL,boxlengthL):     
    eng_iwc = 0
    eng_iwi = 0
    eng_iwb = 0
    N=len(typeL)-1
    Ostarindex=starindexL[0]
    for i in starindexL:
        typei = typeL[i + 1]
        molIDi = IDtomolL[i + 1]        
        chargei = chargetypeL[typei - 1]                              
        for indexj,j in enumerate(waterindexL):      
               Oindex= waterindexL[indexj/3*3]
               typej = typeL[j + 1]
               eps = epstypeL[typei - 1][typej - 1]
               sigma = sigmatypeL[typei - 1][typej - 1]
               molIDj = IDtomolL[j + 1]
               watertype = watertypeL[molIDj]                              
               chargej = chargetypeL[typej - 1]                       
               if molIDi != molIDj:                       
                    sqdistij = sqdistL[i*N+j]                    
                    OOsqdist = sqdistL[Ostarindex*N+Oindex]
                    e1=0
                    if OOsqdist<9*9:
                        xsq = sigma * sigma / sqdistij
                        x6p = xsq * xsq * xsq
                        e1 = 4 * eps * (x6p * x6p - x6p)                    
                    e2=0
                    if OOsqdist<15*15:
                        distij=np.sqrt(sqdistij)
                        e2 = 332.0636 * chargei * chargej / distij
                    if watertype == 1:
                        eng_iwc+=e1 + e2                          
                    elif watertype == 2:
                        eng_iwi+=e1 + e2                           
                    elif watertype == 3:                        
                        eng_iwb+=e1 + e2                                


    return eng_iwc,eng_iwi,eng_iwb

@nb.jit(nopython=True) 
def jitWWenergybymol(watertypeL,waterindexL,typeL,IDtomolL,sqdistL,chargetypeL,epstypeL,sigmatypeL,boxlengthL):  
    
    eng_wwc = 0
    eng_wwi = 0    
    eng_wwb = 0
    
    N = len(typeL) - 1    
    for indexi,i in enumerate(waterindexL):
        Oiindex=waterindexL[indexi/3*3]
        typei = typeL[i + 1]
        molIDi = IDtomolL[i + 1]        
        chargei = chargetypeL[typei - 1]               
        watertype1 = watertypeL[molIDi]
        for indexj,j in enumerate(waterindexL):      
               Ojindex=waterindexL[indexj/3*3]           
               typej = typeL[j + 1]
               eps = epstypeL[typei - 1][typej - 1]
               sigma = sigmatypeL[typei - 1][typej - 1]
               molIDj = IDtomolL[j + 1]
               watertype2 = watertypeL[molIDj]
               if j > i:                   
                   chargej = chargetypeL[typej - 1]                       
                   if molIDi != molIDj:                       
                        sqdistij = sqdistL[N*i+j]                                
                        OOsqdistij=sqdistL[Oiindex*N+Ojindex]
                        e1=0
                        if  OOsqdistij<9*9:          
                            xsq = sigma * sigma / sqdistij
                            x6p = xsq * xsq * xsq
                            e1 = 4 * eps * (x6p * x6p - x6p)
                        e2=0
                        if  OOsqdistij<15*15:   
                           distij = np.sqrt(sqdistij)         
                           e2 = 332.0636 * chargei * chargej / distij

                        if watertype1 == 1:                            
                            eng_wwc+=e1 + e2                            
                        elif watertype1 == 2:                            
                            eng_wwi+=e1 + e2                            
                        elif watertype1 == 3:                            
                            eng_wwb+=e1 + e2
                            
                        if watertype2 == 1:                            
                            eng_wwc+=e1 + e2                            
                        elif watertype2 == 2:                            
                            eng_wwi+=e1 + e2                            
                        elif watertype2 == 3:                            
                            eng_wwb+=e1 + e2                                    


    return eng_wwc,eng_wwi,eng_wwb

def calcSurf(coordsL,typeL,xmin,xmax,ymin,ymax,zmin,zmax,zcom,bin,boxlengthL):
    Nxbin=int((xmax-xmin)/bin)
    Nybin=int((ymax-ymin)/bin)
    Nzbin=int((zmax-zmin)/bin)
    Natom=len(coordsL)-1
    for i in range(Nxbin):
        x=xmin+i*bin+bin/2        
        starttime = datetime.datetime.now()    
        for j in range(Nybin):
            y=ymin+j*bin+bin/2            
            for k in range(Nzbin):
                z=zmin+k*bin+bin/2                
                for l in range(Natom):
                    if typeL[l]==1 or typeL[l]==3:
                        sqdist=kobpair.minimage3Dsqdist(coordsL[l][0],coordsL[l][1],coordsL[l][2],x,y,z,boxlengthL[0],boxlengthL[1],boxlengthL[2])
                        phi=(2*3.141592654*3*3)**(1.5)*math.exp(-sqdist/(2*3*3))
        endtime = datetime.datetime.now()        
        print endtime - starttime
                        
                
    


def compute(q):                     
        changeEVBfct = None
        if q.HorOH == 1:
            changeEVBfct = changeEVBstateH
        else:
            changeEVBfct = changeEVBstateOH
        evbtime = koblib.readcleanevb(q,q.time,q.cleanevbf)          
        if evbtime == None:
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
        






  
        print "time= ",q.time

        starttime = datetime.datetime.now()        
        coordsL = np.array(q.coordsL)
        chargetypeL = np.array(q.chargetypeL)
        esptypeL = np.array(q.epstypeL)
        sigmatypeL = np.array(q.sigmatypeL)
        boxlengthL = np.array(q.boxlengthL)
        
        
        q.zcom = findzcenter2(q)
        #print zcom
        normalization=1./sum([x*x for x in q.evbciL])
        
        
        sqdistL=np.zeros(q.noofatom*q.noofatom)        
        
        sqdistL=kobpair.sqdistmatrix(sqdistL,coordsL,boxlengthL)
        
        if q.calcSurf:
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
            
            #q.surff.close()
            #f=gzip.open("surf.pickle","rb")
            #outf=open("test.xyz","w")
            #while True:
            #    try:
            #        value1=pickle.load(f)
            #        value2=pickle.load(f)
            #        print value1
            #        #info=[int(x) for x in value2.split()]
            #        info=[float(x) for x in value2.split()]
            #        outf.write("%d\n"%(len(info)/3))
            #        outf.write("Comments\n")
            #        for index,i in enumerate(info):
            #            if index%3==0:                            
            #                x=info[index]
            #                y=info[index+1]
            #                z=info[index+2]                            
            #                outf.write(" 1 %f %f %f\n"%(x,y,z))                    
            #    except EOFError:           
            #        break
                    
            #sys.exit(0)
            endtime = datetime.datetime.now()                                
            print endtime - starttime

        if q.calcenergies:
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
        endtime = datetime.datetime.now()        
        print endtime - starttime
        
        

    




