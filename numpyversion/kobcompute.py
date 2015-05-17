import datetime
import re
import sys
import math
import koblib
import vecmath as vec
import numpy

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

def countatoms(map, atomIDlist):
    sum = 0
    for atomID in atomIDlist:
        sum+=len(map[atomID])    
    return sum       


 
def minimagedist(r2,r1,boxlength):
    diff = r2 - r1    
    increment=0
    if diff/boxlength >0.5:
        increment=1
    elif diff/boxlength <-0.5:
        increment=-1      
    return diff - increment*boxlength


def minimage3Dsqdist(r2,r1,boxlengthL):
    sumsq=0
    for i in range(3):
        increment=0
        diff=r2[i]-r1[i]
        if diff/boxlengthL[i]>0.5:
            increment=1
        elif diff/boxlengthL[i] <-0.5:
            increment=-1      
        diff=diff-increment*boxlengthL[i]
        sumsq+=diff*diff
    return sumsq





def findzcenter1(q):  # gibbs dividing surface across all xy-plane
    noofatom = q.statesL[-1][1]
    binsize = 0.5
    bulkdensity = 0.63
    binvol = q.xboxlength * q.yboxlength * binsize

    zhisto = ["local"]
    for i in range(noofatom):
        z = q.statesL[-1][3 + i][q.xl + 2]
        if z > 0:
            typeIndex = q.statesL[-1][3 + i][1] - 1
            whichbin = int(z / binsize) + 1
            koblib.histo(zhisto, binsize, whichbin, q.massL[typeIndex] / binvol)
    firstboundary = -1
    secondboundary = -1
    count = 1
    for item in zhisto[1:]:
        density = 0
        if item[2] > 0:
            density = item[1]
            if firstboundary == -1:
                if density >= bulkdensity / 2:
                    deltadensity = bulkdensity - item[1]
                    slope = (zhisto[count][1] - zhisto[count - 1][1]) / (zhisto[count][0] - zhisto[count - 1][0])
                    firstboundary = item[0] - slope * deltadensity
            elif secondboundary == -1:
                if density <= bulkdensity / 2 and item[0] - firstboundary > 5:
                    deltadensity = bulkdensity - item[1]
                    slope = (zhisto[count][1] - zhisto[count - 1][1]) / (zhisto[count][0] - zhisto[count - 1][0])
                    secondboundary = item[0] + slope * deltadensity
        count += 1

    return (secondboundary + firstboundary) / 2


def findzcenter2(q):  # using z center of mass
    noofatom = q.statesL[-1][1]
    sum = 0
    masssum = 0
    for i in range(noofatom):
        type = q.statesL[-1][i + 3][1]
        zpos = q.statesL[-1][i + 3][4]
        mass = q.massL[type - 1]
        if q.initialzcenterofmass != -1:
            if abs(zpos - q.initialzcenterofmass) < 50:  #not letting gas molecules affect this
                sum += q.massL[type - 1] * zpos
                masssum += mass
        else:
            sum += q.massL[type - 1] * zpos
            masssum += mass
    if (q.initialzcenterofmass) == -1:
        q.initialzcenterofmass = sum / masssum
    return sum / masssum


def findzcenter3(q):  # using z center of mass for each \Delta x \Delta y
    noofatom = q.statesL[-1][1]

    xwidth = q.xhi - q.xlo
    ywidth = q.yhi - q.ylo
    noofxybins = 10
    xbinsize = xwidth / noofxybins
    ybinsize = ywidth / noofxybins

    sum = []
    for i in range(noofxybins):
        L = []
        for j in range(noofxybins):
            L.append(0)
        sum.append(L)

    masssum = []
    for i in range(noofxybins):
        L = []
        for j in range(noofxybins):
            L.append(0)
        masssum.append(L)
    zposL = []
    for i in range(noofxybins):
        L = []
        for j in range(noofxybins):
            L.append(0)
        zposL.append(L)

    for i in range(noofatom):
        type = q.statesL[-1][i + 3][1]
        xpos = q.statesL[-1][i + 3][2]
        if xpos <= q.xlo:
            while xpos <= q.xlo:
                xpos += xwidth
        elif xpos >= q.xhi:
            while xpos >= q.xhi:
                xpos -= xwidth
        ypos = q.statesL[-1][i + 3][3]
        if ypos <= q.ylo:
            while ypos <= q.ylo:
                ypos += ywidth
        elif ypos >= q.yhi:
            while ypos >= q.yhi:
                ypos -= ywidth
        whichxbin = int(xpos / xbinsize)
        whichybin = int(ypos / ybinsize)
        zpos = q.statesL[-1][i + 3][4]
        mass = q.massL[type - 1]

        if (q.initialzcenterofmass) != -1:
            if abs(zpos - q.initialzcenterofmass) < 50:
                sum[whichxbin][whichybin] += mass * zpos
                masssum[whichxbin][whichybin] += mass
        else:
            sum[whichxbin][whichybin] += mass * zpos
            masssum[whichxbin][whichybin] += mass

    if (q.initialzcenterofmass) == -1:
        tempsum = 0
        tempmasssum = 0
        for i in range(noofxybins):
            for j in range(noofxybins):
                tempsum += sum[i][j]
                tempmasssum += masssum[i][j]
        q.initialzcenterofmass = tempsum / tempmasssum
    for i in range(noofxybins):
        for j in range(noofxybins):
            zposL[i][j] = sum[i][j] / masssum[i][j]

    return zposL

    return (secondboundary + firstboundary) / 2


def interfacelocalcompute(q):
    noofatom = q.statesL[-1][1]

    variance = 0
    xwidth = q.xhi - q.xlo
    ywidth = q.yhi - q.ylo
    noofxybins = 10
    xbinsize = xwidth / noofxybins
    ybinsize = ywidth / noofxybins

    zcenterL = findzcenter3(q)

    q.bigionL[0] += 1
    ionbinsize = 0.1

    for type in q.iontype:
        for index in q.Typemap[type]:
            if type != 4 or (type == 4 and (index - q.Typemap[4]) % 6 == 0):
                for L in q.bigionL[1:]:
                    if type == L[0]:
                        xpos = q.statesL[-1][index][2]
                        if xpos <= q.xlo:
                            while xpos <= q.xlo:
                                xpos += xwidth
                        elif xpos >= q.xhi:
                            while xpos >= q.xhi:
                                xpos -= xwidth
                        ypos = q.statesL[-1][index][3]
                        if ypos <= q.ylo:
                            while ypos <= q.ylo:
                                ypos += ywidth
                        elif ypos >= q.yhi:
                            while ypos >= q.yhi:
                                ypos -= ywidth
                        whichxbin = int(xpos / xbinsize)
                        whichybin = int(ypos / ybinsize)
                        zcenter = zcenterL[whichxbin][whichybin]
                        zpos = q.statesL[-1][index][q.xl + 2] + q.statesL[-1][index][q.xl + 5] * q.zboxlength
                        # zcenter=zcenterL[whichxbin][whichybin]
                        if zpos - zcenter > 0:
                            whichbin = int((zpos - zcenter) / ionbinsize) + 1
                        else:
                            whichbin = int((zcenter - zpos) / ionbinsize) + 1
                        koblib.histo(L, ionbinsize, whichbin, 1)


def findwatertypes(q):
    intwaternum = 0
    bulkwaternum = 0
    coordwaternum = 0    
    OstarID = q.molL[0][0][0]    
    Ostarcoords=q.statesL[-1][OstarID+2][q.xl:q.xl+3]
    waterindexL = q.Typemap[q.OWID]
    zcom = findzcenter2(q)
    for index in waterindexL:
        watercoords = q.statesL[-1][index][q.xl:q.xl + 3]
        xdist = minimagedist(Ostarcoords[0],watercoords[0],q.xboxlength)
        ydist = minimagedist(Ostarcoords[1],watercoords[1],q.yboxlength)
        zdist = Ostarcoords[2] - watercoords[2]
        sqdist = xdist ** 2 + ydist ** 2 + zdist ** 2
        if sqdist < 12.25: #3.5*3.5
            coordwaternum+=1
        else:
            if math.fabs(watercoords[2] - zcom) > 13.5:
                intwaternum+=1
            else:
                bulkwaternum+=1
    print coordwaternum,intwaternum,bulkwaternum

                
def readcleanevb(q,time,cleanevbf):
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
                break

def backupevbstate(q):
    for i in range(4):
        molBID = q.evbmolBL[i]
        tempL=[]
        for item in q.molL[molBID]:
            tempL.append(item[:])
        q.backupmolL.append([molBID,tempL])
    for i in range(4):
        molBID = q.evbmolBL[i]
        for atomL in q.molL[molBID][1:]:
            atomID = atomL[0]
            tempL=[]
            q.backupidtomolL.append([atomID,q.idtomolL[atomID]])
    
        
                
                
def changeEVBstateOH(q,whichstate):
    first=None
    if len(q.backupmolL)==0:
        first=True
    reactionshell = q.evbshellL[whichstate]
    previousstate=0
    if reactionshell == 1 and q.currentEVBstate != 0 and not first:        
        for item in q.backupmolL:
            tempL=[]
            for subitem in item[1]:
                tempL.append(subitem[:])
            q.molL[item[0]] = tempL
        for item in q.backupidtomolL:
            q.idtomolL[item[0]] = item[1]
    if first:
        previousstate=q.currentEVBstate
    if reactionshell == 2:       
        molA = q.evbmolAL[whichstate]
        previousstate = q.evbmolBL.index(molA)         
        changeEVBstateOH(q,previousstate)
    elif reactionshell > 2:
        sys.exit("Reaction shell >2!\n")                   

    oldwatermolID = q.evbmolBL[whichstate]
    oldwaterID = q.molL[oldwatermolID][0][0]
    oldwaterH1ID = q.molL[oldwatermolID][1][0]
    #print q.time,q.evbmolBL, q.molL[oldwatermolID]
    oldwaterH2ID = q.molL[oldwatermolID][2][0]
    oldOstarmolID = q.evbmolBL[previousstate]
    oldOstarID = q.molL[oldOstarmolID][0][0]    
    oldHstarID = q.molL[oldOstarmolID][1][0]    
        
    
    sqdist1 = q.distsq2D[oldwaterH1ID-1][oldOstarID-1]
    sqdist2 = q.distsq2D[oldwaterH2ID-1][oldOstarID-1]

    q.molL[oldOstarmolID][0][1] = 1
    q.molL[oldOstarmolID][1][1] = 2
    
    q.molL[0][0] = [oldwaterID,3]
    if sqdist1 < sqdist2:
        q.molL[oldOstarmolID].append([oldwaterH1ID,2])
        q.idtomolL[oldwaterH1ID]= oldOstarmolID
        del q.molL[oldwatermolID][1]
        q.molL[0][1] = [oldwaterH2ID,4]    
    else:
        q.molL[oldOstarmolID].append([oldwaterH2ID,2])
        q.idtomolL[oldwaterH2ID]= oldOstarmolID
        del q.molL[oldwatermolID][2]    
        q.molL[0][1] = [oldwaterH1ID,4]    
    q.molL[oldwatermolID][0][1] = 3
    q.molL[oldwatermolID][1][1] = 4
    
    if not first:
        newOstarID=q.molL[oldwatermolID][0][0]
        newHstarID=q.molL[oldwatermolID][1][0]
        for i in xrange(q.noofatom):
            typei=q.statesL[-1][i+3][1]-1
            chargei=q.chargetypeL[typei]            
            q.npcharge2D[oldOstarID-1][i]=chargei*q.chargetypeL[0]
            q.npcharge2D[oldHstarID-1][i]=chargei*q.chargetypeL[1]
            q.npcharge2D[newOstarID-1][i]=chargei*q.chargetypeL[2]
            q.npcharge2D[newHstarID-1][i]=chargei*q.chargetypeL[3]
            q.npcharge2D[i][oldOstarID-1]=chargei*q.chargetypeL[0]
            q.npcharge2D[i][oldHstarID-1]=chargei*q.chargetypeL[1]
            q.npcharge2D[i][newOstarID-1]=chargei*q.chargetypeL[2]
            q.npcharge2D[i][newHstarID-1]=chargei*q.chargetypeL[3]
            q.eps2DL[oldOstarID-1][i]=q.epstypeL[typei][0]
            q.eps2DL[oldHstarID-1][i]=q.epstypeL[typei][1]
            q.eps2DL[newOstarID-1][i]=q.epstypeL[typei][2]
            q.eps2DL[newHstarID-1][i]=q.epstypeL[typei][3]
            q.eps2DL[i][oldOstarID-1]=q.epstypeL[typei][0]
            q.eps2DL[i][oldHstarID-1]=q.epstypeL[typei][1]
            q.eps2DL[i][newOstarID-1]=q.epstypeL[typei][2]
            q.eps2DL[i][newHstarID-1]=q.epstypeL[typei][3]
            q.sigma2DL[oldOstarID-1][i]=q.sigmatypeL[typei][0]
            q.sigma2DL[oldHstarID-1][i]=q.sigmatypeL[typei][1]
            q.sigma2DL[newOstarID-1][i]=q.sigmatypeL[typei][2]
            q.sigma2DL[newHstarID-1][i]=q.sigmatypeL[typei][3]
            q.sigma2DL[i][oldOstarID-1]=q.sigmatypeL[typei][0]
            q.sigma2DL[i][oldHstarID-1]=q.sigmatypeL[typei][1]
            q.sigma2DL[i][newOstarID-1]=q.sigmatypeL[typei][2]
            q.sigma2DL[i][newHstarID-1]=q.sigmatypeL[typei][3]
           
            
    
    q.currentEVBstate = whichstate

def oldLJenergy(q):  
    LJsum=0
    coulsum=0    
    for i in xrange(q.noofatom):
        typei=q.statesL[-1][i+3][1]
        molIDi=q.idtomolL[i+1]
        coordsi=q.statesL[-1][i+3][q.xl:q.xl+3]
        chargei=q.chargetypeL[typei-1]                
        for j in xrange(q.noofatom):          
               typej=q.statesL[-1][j+3][1]
               eps=q.epstypeL[typei-1][typej-1]
               sigma=q.sigmatypeL[typei-1][typej-1]
               molIDj=q.idtomolL[j+1]
               if j>i:
                   coordsj=q.statesL[-1][j+3][q.xl:q.xl+3]
                   chargej=q.chargetypeL[typej-1]
                   #sqdistij=minimage3Dsqdist(coordsi,coordsj,q.boxlengthL)
                   sqdistij=q.distsq2D[i][j]
                   distij=math.sqrt(sqdistij)                   
                   if molIDi!=molIDj:
                       e1=4*eps*((sigma/distij)**12-(sigma/distij)**6)
                       LJsum+=e1    
                       e2=332.0636*chargei*chargej/distij
                       coulsum+=e2                       
                       #if i==0 and j==2981:
                       #    print sigma,eps,distij,e1                          
    print "old:", coulsum, LJsum
    return coulsum, LJsum
def LJcoulenergy(q):
    q.coulenergy_mat=q.npcharge2D*q.inv1p2D
    q.LJenergy_mat=4*q.npeps2D*(numpy.power(q.npsigma2D*q.inv1p2D,12)-numpy.power(q.npsigma2D*q.inv1p2D,6))
    
    coulenergy=q.coulenergy_mat.sum()
    LJenergy=q.LJenergy_mat.sum()
    print "numpy", 0.5*332.0636*coulenergy, 0.5*LJenergy
    return 0.5*(332.0636*coulenergy+LJenergy)

def initparam2D(q):
    q.chargeL=[]
    q.eps2DL=[]      
    q.sigma2DL=[]
    for i in xrange(q.noofatom):
        itype=q.statesL[-1][i+3][1]-1
        chargetempL=[]
        epstempL=[]
        sigmatempL=[]  
        q.chargeL.append(q.chargetypeL[itype])
        for j in xrange(q.noofatom):
            jtype=q.statesL[-1][j+3][1]-1            
            epstempL.append(q.epstypeL[itype][jtype])
            sigmatempL.append(q.sigmatypeL[itype][jtype])
        q.eps2DL.append(epstempL[:])
        q.sigma2DL.append(sigmatempL[:])

def calcsqdistmatrixnumpy(q):
    r=q.npcoordsL    
    diff=r[:,:,None] - r[:,:,None].T
    diff-=numpy.round(diff/q.npboxlength)*q.npboxlength      
    ssdiff=numpy.multiply(diff,diff,out=diff).sum(1)
    
    return ssdiff

def compute(q):                     
        starttime=datetime.datetime.now()
        #if len(q.chargeL)==0:            
        initparam2D(q)            
        readcleanevb(q,q.time,q.cleanevbf)          
        q.npcoordsL=numpy.array(q.coordsL)        
        q.distsq2D=calcsqdistmatrixnumpy(q)        
        ####setup intial state####       
        q.currentEVBstate=q.evbciL.index(max(q.evbciL))        
        q.backupmolL = []        
        q.backupidtomolL = []               
        if q.currentEVBstate!=0:
            changeEVBstateOH(q,0)
        q.npcharge=numpy.array(q.chargeL)
        q.npcharge2D=q.npcharge[:,None]*q.npcharge[:,None].T          
        q.npsigma2D=numpy.array(q.sigma2DL)      
        q.npeps2D=numpy.array(q.eps2DL)      
        backupevbstate(q)
        ###setup intial state#### 
        q.npidtomol=numpy.array(q.idtomolL[1:])
        q.npboolmol2D=q.npidtomol[:,None]!=q.npidtomol[:,None].T   #same mol ID gets 0          
        q.inv1p2D=numpy.where(q.npboolmol2D==True,numpy.sqrt(1./q.distsq2D),0)    
        print "time= ",q.time
        energy=LJcoulenergy(q)
        #oldLJenergy(q)
        anothersum=0
        anothersum2=0
        diffenergysum=0
        #for i in range(q.noofatom):
        #    molIDi=q.idtomolL[i+1]
        #    for j in range(q.noofatom):
        #        molIDj=q.idtomolL[j+1]
        #        if j>i and molIDi!=molIDj:
        #            diff=q.LJenergy_mat[i][j]-q.LJenergy_mat[j][i]
        #            if diff>0.00000001:
        #                print i,j, q.LJenergy_mat[i][j],q.LJenergy_mat[j][i],diff
        #                sys.exit()
                    
        #for i in range(q.noofatom):
        #    molIDi=q.idtomolL[i+1]
        #    print i, q.LJenergy_mat[i][i]
        #    for j in range(q.noofatom):
        #        molIDj=q.idtomolL[j+1]
        #        if j>i: 
        #             if molIDi==molIDj:
        #                if q.LJenergy_mat[i][j]>0:
        #                    print i,j,q.LJenergy_mat[i][j]
        #             else:                
        #                diffenergy=q.LJenergy_mat[i][j]-q.LJij[i][j]
        #                anothersum+=q.LJenergy_mat[i][j]
        #                anothersum2+=q.LJij[i][j]
        #                if diffenergy>1e-10:
        #                    print i,j, diffenergy
        #                    sys.exit()
               
                    
        #            #print("i=%d j=%d, %f"%(i,j,diffenergy))
        #print anothersum, anothersum2,q.LJenergy_mat.sum()*0.5
        #sys.exit()
        #i=0;j=2981
        #print("i=%d j=%d, diff=%f"%(i,j,q.LJenergy_mat[i][j]-q.LJij[i][j]))
        #print i,j,q.npsigma2D[i][j],q.npeps2D[i][j],math.sqrt(q.distsq2D[i][j]),q.LJenergy_mat[i][j]
        #print("i=%d j=%d, %s,%s"%(i,j,q.statesL[-1][i+3][q.xl:q.xl+3],q.statesL[-1][j+3][q.xl:q.xl+3]))
       # sys.exit()

        #if q.HorOH == 1:
        #    changeEVBstateH(q,2)
        #else:
        #    changeEVBstateOH(q,1)
        #    changeEVBstateOH(q,2)
        #    changeEVBstateOH(q,3)
        #    changeEVBstateOH(q,0)
        endtime=datetime.datetime.now()
        print endtime-starttime
        

    



