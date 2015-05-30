import cython
import sys
from libc.math cimport sqrt
from libc.math cimport exp
from libc.math cimport pow
from libc.math cimport abs
import numpy as np
cimport numpy as np
from cython.parallel import prange

ctypedef np.float64_t DTYPE_t
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
def sqdistmatrix(np.ndarray[DTYPE_t, ndim=1,mode="c"]sqdistL,double[:,::1] coordsL,double[:] boxlengthL):
    cdef int i,j,k
    cdef double sqdistij,diff
    cdef int N=coordsL.shape[0]-1
    for i in range(N):
        for j in range(N):
            if j>i:
                sqdistij = 0                                      
                for k in range(3):
                    diff = coordsL[i+1,k] - coordsL[j+1,k]
                    if diff / boxlengthL[k] > 0.5:
                        diff-= boxlengthL[k]
                    elif diff / boxlengthL[k] < -0.5:
                        diff+= boxlengthL[k]                                         
                    sqdistij+=diff * diff 
                sqdistL[i*N+j]= sqdistij      
                sqdistL[j*N+i]= sqdistij 
    return sqdistL

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
def sqdistOOmatrix(double[:,::1] sqOOdistL,double[:,::1] coordsL,int[:] OIDL, double[:] boxlengthL):
    cdef int i,j,k
    cdef double sqdistij,diff
    cdef int N=OIDL.shape[0]    
    
    counti=0    
    for ii in range(N):                            
       countj=0 
       i=OIDL[ii]
       for jj in range(N):                                   
          j=OIDL[jj]
          if j>i:                
            sqdistij = 0                                      
            for k in range(3):
                diff = coordsL[i,k] - coordsL[j,k]
                if diff / boxlengthL[k] > 0.5:
                    diff-= boxlengthL[k]
                elif diff / boxlengthL[k] < -0.5:
                    diff+= boxlengthL[k]                                         
                sqdistij+=diff * diff 
            sqOOdistL[ii,jj]= sqdistij 
            sqOOdistL[jj,ii]= sqdistij 
              

         
    return sqOOdistL

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
def WWenergy(int[:] watertypeL, int[:] waterindexL,int[:] typeL,int[:] idtomolL,double[:] sqdistL,double[:] chargetypeL,double[:,::1] epstypeL,double[:,::1]sigmatypeL,double[:]boxlengthL,int N):
    cdef double eng_wwc = 0
    cdef double eng_wwi = 0
    cdef double eng_wwb = 0
    cdef int typei, molIDi,typej,molIDj,Oiindex,Ojindex,ii,i,indexi,j,jj,indexj,watertype1,watertype2
    cdef double chargei,chargej,eps,sigma,distij,OOsqdist, sqdistij,e1,e2,xsq,x6p
    #cdef int wrange=len(waterindexL)
    indexi=0
    for ii in range(waterindexL.shape[0]):
        i=waterindexL[ii]
        Oiindex=waterindexL[indexi/3*3]
        typei = typeL[i + 1]
        molIDi = idtomolL[i + 1]        
        chargei = chargetypeL[typei - 1]                
        watertype1 = watertypeL[molIDi]
        indexj=0
        for jj in range(waterindexL.shape[0]):      
               j=waterindexL[jj]
               Ojindex=waterindexL[indexj/3*3]           
               typej = typeL[j + 1]
               molIDj = idtomolL[j + 1]
               watertype2 = watertypeL[molIDj]
               if j > i:                   
                   chargej = chargetypeL[typej - 1]                       
                   if molIDi != molIDj:                       
                        sqdistij = sqdistL[N*i+j]                                
                        OOsqdist=sqdistL[Oiindex*N+Ojindex]
                        e1=0
                        if  OOsqdist<9*9:          
                            sigma=sigmatypeL[typei-1,typej-1]
                            eps=epstypeL[typei-1,typej-1]
                            xsq = sigma * sigma / sqdistij
                            x6p = xsq * xsq * xsq
                            e1 = 4 * eps * (x6p * x6p - x6p)
                        e2=0
                        if  OOsqdist<15*15:   
                           distij = sqrt(sqdistij)         
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
               indexj+=1
        indexi+=1
    return eng_wwc,eng_wwi,eng_wwb

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
def IWenergy(int[:] watertypeL,int[:] starindexL, int[:] waterindexL,int[:] typeL,int[:] idtomolL,double[:] sqdistL,double[:] chargetypeL,double[:,::1] epstypeL,double[:,::1]sigmatypeL,double[:]boxlengthL,int N):     
    cdef double eng_iwc = 0
    cdef double eng_iwi = 0
    cdef double eng_iwb = 0
    cdef int Ostarindex=starindexL[0]
    cdef int typei, molIDi,typej,molIDj,Oindex,i,j,indexj,watertype
    cdef double chargei,chargej,eps,sigma,distij,OOsqdist, sqdistij,e1,e2,xsq,x6p
    cdef int ii,jj
    #cdef int starrange=len(starindexL)
    #cdef int wrange=len(waterindexL)

    for ii in range(starindexL.shape[0]):
        i=starindexL[ii]
        typei = typeL[i + 1]
        molIDi = idtomolL[i + 1]        
        chargei = chargetypeL[typei - 1]                
        indexj=0
        for jj in range(waterindexL.shape[0]):
               j=waterindexL[jj]      
               Oindex= waterindexL[indexj/3*3]
               typej = typeL[j + 1]
               molIDj = idtomolL[j + 1]
               watertype = watertypeL[molIDj]                              
               chargej = chargetypeL[typej - 1]                       
               if molIDi != molIDj:                       
                    sqdistij = sqdistL[i*N+j]                    
                    OOsqdist = sqdistL[Ostarindex*N+Oindex]
                    e1=0
                    if OOsqdist<9*9:
                        sigma=sigmatypeL[typei-1,typej-1]
                        eps=epstypeL[typei-1,typej-1]
                        xsq = sigma * sigma / sqdistij
                        x6p = xsq * xsq * xsq
                        e1 = 4 * eps * (x6p * x6p - x6p)                    
                    e2=0
                    if OOsqdist<15*15:
                        distij=sqrt(sqdistij)
                        e2 = 332.0636 * chargei * chargej / distij
                    if watertype == 1:
                        eng_iwc+=e1 + e2                          
                    elif watertype == 2:
                        eng_iwi+=e1 + e2                           
                    elif watertype == 3:                        
                        eng_iwb+=e1 + e2                                
               indexj+=1

    return eng_iwc,eng_iwi,eng_iwb

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
def fullWWenergy(int[:] watertypeL, int[:] waterindexL,int[:] typeL,int[:] idtomolL,double[:] sqdistL,double[:] chargetypeL,double[:,::1] epstypeL,double[:,::1]sigmatypeL,double[:]boxlengthL,int N):
    cdef double eng_wwc = 0
    cdef double eng_wwi = 0
    cdef double eng_wwb = 0
    cdef int typei, molIDi,typej,molIDj,Oiindex,Ojindex,ii,i,indexi,j,jj,indexj,watertype1,watertype2
    cdef double chargei,chargej,eps,sigma,distij,OOsqdist, sqdistij,e1,e2,xsq,x6p
    #cdef int wrange=len(waterindexL)
    indexi=0
    for ii in range(waterindexL.shape[0]):
        i=waterindexL[ii]
        Oiindex=waterindexL[indexi/3*3]
        typei = typeL[i + 1]
        molIDi = idtomolL[i + 1]        
        chargei = chargetypeL[typei - 1]                
        watertype1 = watertypeL[molIDi]
        indexj=0
        for jj in range(waterindexL.shape[0]):      
               j=waterindexL[jj]
               Ojindex=waterindexL[indexj/3*3]           
               typej = typeL[j + 1]
               molIDj = idtomolL[j + 1]
               watertype2 = watertypeL[molIDj]
               if j > i:                   
                   chargej = chargetypeL[typej - 1]                       
                   if molIDi != molIDj:                       
                        sqdistij = sqdistL[N*i+j]                                
                        OOsqdist=sqdistL[Oiindex*N+Ojindex]
                        sigma=sigmatypeL[typei-1,typej-1]
                        eps=epstypeL[typei-1,typej-1]
                        xsq = sigma * sigma / sqdistij
                        x6p = xsq * xsq * xsq
                        e1 = 4 * eps * (x6p * x6p - x6p)
                        distij = sqrt(sqdistij)         
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
               indexj+=1
        indexi+=1
    return eng_wwc,eng_wwi,eng_wwb



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
def fullIWenergy(int[:] watertypeL,int[:] starindexL, int[:] waterindexL,int[:] typeL,int[:] idtomolL,double[:] sqdistL,double[:] chargetypeL,double[:,::1] epstypeL,double[:,::1]sigmatypeL,double[:]boxlengthL,int N):     
    cdef double eng_iwc = 0
    cdef double eng_iwi = 0
    cdef double eng_iwb = 0
    cdef int Ostarindex=starindexL[0]
    cdef int typei, molIDi,typej,molIDj,Oindex,i,j,indexj,watertype
    cdef double chargei,chargej,eps,sigma,distij,OOsqdist, sqdistij,e1,e2,xsq,x6p
    cdef int ii,jj
    #cdef int starrange=len(starindexL)
    #cdef int wrange=len(waterindexL)

    for ii in range(starindexL.shape[0]):
        i=starindexL[ii]
        typei = typeL[i + 1]
        molIDi = idtomolL[i + 1]        
        chargei = chargetypeL[typei - 1]                
        indexj=0
        for jj in range(waterindexL.shape[0]):
               j=waterindexL[jj]      
               Oindex= waterindexL[indexj/3*3]
               typej = typeL[j + 1]
               molIDj = idtomolL[j + 1]
               watertype = watertypeL[molIDj]                              
               chargej = chargetypeL[typej - 1]                       
               if molIDi != molIDj:                       
                    sqdistij = sqdistL[i*N+j]                    
                    OOsqdist = sqdistL[Ostarindex*N+Oindex]                    
                    sigma=sigmatypeL[typei-1,typej-1]
                    eps=epstypeL[typei-1,typej-1]
                    xsq = sigma * sigma / sqdistij
                    x6p = xsq * xsq * xsq
                    e1 = 4 * eps * (x6p * x6p - x6p)                    
                    
                    distij=sqrt(sqdistij)
                    e2 = 332.0636 * chargei * chargej / distij
                    if watertype == 1:
                        eng_iwc+=e1 + e2                          
                    elif watertype == 2:
                        eng_iwi+=e1 + e2                           
                    elif watertype == 3:                        
                        eng_iwb+=e1 + e2                                
               indexj+=1

    return eng_iwc,eng_iwi,eng_iwb


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cdef double minimage3Dsqdist(double r1x,double r1y,double r1z,double r2x, double r2y, double r2z,double boxx,double boxy, double boxz):
    cdef double sum = 0
    cdef double diff
    cdef int k
    cdef double *r1 = [r1x,r1y,r1z]    
    cdef double *r2 = [r2x,r2y,r2z]    
    cdef double *boxlengthL = [boxx,boxy,boxz]    
    for k in range(3):
        diff = r2[k] - r1[k]
        if diff / boxlengthL[k] > 0.5:
            diff-= boxlengthL[k]
        elif diff / boxlengthL[k] < -0.5:
            diff+= boxlengthL[k]                                         
        sum+=diff*diff 
    return sum

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
def calcSurf(double[:,:,::1] gaussmatrix,double[:,::1] surfacepoint,double[:,::1] coordsL,int[:] typeL,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,double zcom,double bin,double[:] boxlengthL):    
    cdef int Nxbin,Nybin,Nzbin, Natom,i,j,k,l,count
    cdef double x,y,z,sum,sqdist,diff,phi,value
    cdef double boxx=boxlengthL[0]
    cdef double boxy=boxlengthL[1]
    cdef double boxz=boxlengthL[2]
    cdef double xi=2.4
    cdef double factor=pow(2*3.141592654*xi*xi,-1.5)
    #cdef double* r1=[0,0,0]
    cdef double* r2=[0,0,0]
    cdef double* a=[0.0,0.0,0.0]
    cdef int index =0
    Nxbin=<int>((xmax-xmin)/bin)
    Nybin=<int>((ymax-ymin)/bin)
    Nzbin=<int>((zmax-zmin)/bin)
    Natom=coordsL.shape[0]-1       
    for l in range(Natom):                
        if (typeL[l+1]==1 or typeL[l+1]==3):# and coordsL[l,2]-zcom>0
            sum=0.0                
            for k in range(Nzbin):
                z=zmin+k*bin+bin/2          
                for i in range(Nxbin):
                    x=xmin+i*bin+bin/2                
                    for j in range(Nybin):
                        #print i,j,k,x,y,z
                        y=ymin+j*bin+bin/2                
                        sqdist=0.0
                        r2[0]=x;r2[1]=y;r2[2]=z+zcom						
                        for m in range(3):
                            diff = r2[m] - coordsL[l+1,m]
                            if diff / boxlengthL[m] > 0.5:
                                diff-= boxlengthL[m]
                            elif diff / boxlengthL[m] < -0.5:
                                diff+= boxlengthL[m]                                         
                            sqdist+=diff*diff
                        #print sqdist ,factor					                 	                                         
                        gaussmatrix[i,j,k]+=factor*exp(-sqdist/(2*xi*xi))  
    count=0                                       
    for k in range(Nzbin):
       z=zmin+k*bin+bin/2          
       for i in range(Nxbin):
          x=xmin+i*bin+bin/2                
          for j in range(Nybin):
             y=ymin+j*bin+bin/2 
             value=gaussmatrix[i,j,k]
             if value>0.0155 and value<0.016:
                 #print count
                 surfacepoint[count,0]=x
                 surfacepoint[count,1]=y
                 surfacepoint[count,2]=z
                 count+=1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
def calcSurf2(double[:,::1] surfacepoint,double[:,::1] coordsL,int[:] typeL,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,double zcom,double xybin,double zbin,double[:] boxlengthL,double level,double tol):    
    cdef int lowbin,highbin,Nxbin,Nybin,Nzbin, Natom,i,j,k,l,searchcount
    cdef double x,y,z,sum,sqdist,diff,phi,value
    cdef double boxx=boxlengthL[0]
    cdef double boxy=boxlengthL[1]
    cdef double boxz=boxlengthL[2]
    cdef double xi=2.4
    cdef double factor=pow(2*3.141592654*xi*xi,-1.5)
    #cdef double* r1=[0,0,0]
    cdef double* r2=[0,0,0]
    cdef double* a=[0.0,0.0,0.0]
    cdef int index =0
    Nxbin=<int>((xmax-xmin)/xybin)
    Nybin=<int>((ymax-ymin)/xybin)
    Nzbin=<int>((zmax-zmin)/zbin)
    Natom=coordsL.shape[0]-1      
    for i in range(Nxbin):
        x=xmin+i*xybin+xybin/2                
        for j in range(Nybin):
            #print i,j,k,x,y,z
            y=ymin+j*xybin+xybin/2        
            lowbin=0
            highbin=Nzbin-1    
            searchcount=0          
            while lowbin<=highbin:
                k=(lowbin+highbin)/2
                z=zmin+k*zbin+zbin/2                                 
                sum=0.0 
                r2[0]=x;r2[1]=y;r2[2]=z+zcom	                  
                for l in range(Natom):                
                    if (typeL[l+1]==1 or typeL[l+1]==3):                        					
                        sqdist=0.0
                        for m in range(3):
                            diff = r2[m] - coordsL[l+1,m]
                            if diff / boxlengthL[m] > 0.5:
                                diff-= boxlengthL[m]
                            elif diff / boxlengthL[m] < -0.5:
                                diff+= boxlengthL[m]                                         
                            sqdist+=diff*diff                        			                 	                                         
                        sum+=factor*exp(-sqdist/(2*xi*xi))                    
                if abs(sum-level)<tol:             
                    break
                elif sum>level:                                          
                     lowbin=k+1                     
                else:                     
                     highbin=k-1   
                searchcount+=1         
            surfacepoint[i*Nxbin+j,0]=x
            surfacepoint[i*Nxbin+j,1]=y
            surfacepoint[i*Nxbin+j,2]=z
            #print z,sum
            #print searchcount
            
            




 

