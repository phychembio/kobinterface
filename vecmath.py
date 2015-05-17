import math

def minimagedist(r2,r1,boxlength):
    diff = r2 - r1
    while math.fabs(diff) > boxlength / 2:
        if diff > 0:
            diff=diff-boxlength
        else:
            diff=diff+boxlength
    return diff

def dotvec(v1,v2):
    sum=0
    for i in range(len(v1)):
       sum+=v1[i]*v2[i]
    return sum

def lenvec(v):
    sum=0
    for i in range(len(v)):
        sum+=v[i]*v[i]
    return math.sqrt(sum)

def addvec(v1,v2):
    v=[]
    for i in range(len(v1)):
       v.append(v1[i]+v2[i])    
    return v

def subtractvec(v1,v2):
    v=[]
    for i in range(len(v1)):
       v.append(v1[i]-v2[i])    
    return v

def crossvec(v1,v2):
    v=[0,0,0]
    v[0]=v1[1]*v2[2]-v1[2]*v2[1]
    v[1]=v1[2]*v2[0]-v1[0]*v2[2]
    v[2]=v1[0]*v2[1]-v1[1]*v2[0]
    return v

def scalermultiply(c,v):
    v2=[]
    for i in range(len(v)):
       v2.append(c*v[i])    
    return v2

def uniterize(v):
    length=lenvec(v)
    return scalermultiply(1/length,v)

def anglev1v2(v1,v2):
    dotp=0
    lenv1sq=0
    lenv2sq=0
    for i in range(3):                
        dotp+=v1[i]*v2[i]
    crossp=crossvec(v1,v2)
    angle=math.atan2(lenvec(crossp),dotp)  #atan2(y,x) not (x,y)
    if crossp[2]<0:
        angle=2*3.141592654-angle    
    return angle

def angle(A,B,C,boxlength): #angle between B->A and B->C
    AB=[]
    CB=[]
    for i in range(3):
        AB.append(minimagedist(A[i],B[i],boxlength))
        CB.append(minimagedist(C[i],B[i],boxlength))
    lenAB=0
    lenCB=0
    cdot=0
    for i in range(3):
        lenAB+=AB[i]**2
        lenCB+=CB[i]**2
        cdot+=AB[i]*CB[i]
    lenAB=math.sqrt(lenAB)
    lenCB=math.sqrt(lenCB)

    product=cdot/lenAB/lenCB
    if product>1:
        product=1
    elif product<-1:
        product=-1
    return math.acos(product)


def formminimagevec(A,B,boxlength): #form a vector from point B to point A.
    tempv=[]
    for i in range(3):
       tempv.append(minimagedist(A[i],B[i],boxlength))
    return tempv


def minimagevec(v,boxlength):
    tempv=[]
    for i in range(3):
       tempv.append(minimagedist(v,[0,0,0],boxlength))
    return tempv

