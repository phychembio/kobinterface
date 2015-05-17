import re

basedir="OH-interface/"
filename="out"
inf=open(basedir+filename+".evb","r")
outf=open(basedir+filename+".cleanevb","w")
Nstates=None
evbinfoL=[]
CECcoords=[]
timestep=-1
laststep=-1
for line in inf:
    #print line,
    match1=re.search("TIMESTEP (\d+)",line)
    match2=re.search("STATES",line)
    match3=re.search("EIGEN_VECTOR",line)
    match4=re.search("CEC\_COORDINATE",line)
    match5=re.search('COMPLEX 1: (\d+)',line)
    match6=re.search('^DIAGONAL ',line)
    match7=re.search('OFF-DIAGONAL',line)
    if match1:
        timestep=int(match1.group(1))
        if timestep%10000==0:
           print timestep
        evbinfoL=[]
        CECcoords=[]        
    if match5:
        Nstates=int(match5.group(1))        
    if match2:
        for i in range(Nstates):
            info=inf.next().split()
            evbinfoL.append([-1,int(info[4]),int(info[3]),int(info[2]),-1,-1]) #ci,molB,molA,rep,offdiagonal
    if match3:
        cL=[float(x) for x in inf.next().split()]
        for i in range(Nstates):
            evbinfoL[i][0]=cL[i]        
    if match6:         
        
        for i in range(Nstates):
             info=inf.next().split()                          
             evbinfoL[i][4]=float(info[-1])
    if match7:        
        for i in range(Nstates-1):
            info=inf.next().split()
            evbinfoL[i+1][5]=float(info[1])
    if match4 and timestep%100==0 and laststep!=timestep:
        evbinfoL=[evbinfoL[0]]+sorted(evbinfoL[1:], reverse=True)
        CECcoords=[float(x) for x in inf.next().split()]
        outf.write("--\n")
        outf.write("%d \n"%timestep)        
        for item in evbinfoL[:4]:
            outf.write("%f "%item[0])
        outf.write("\n")
        for item in evbinfoL[:4]:
            outf.write("%d "%item[1])
        outf.write("\n")
        for item in evbinfoL[:4]:
            outf.write("%d "%item[2])
        outf.write("\n")
        for item in evbinfoL[:4]:
            outf.write("%d "%item[3])
        outf.write("\n")
        for item in evbinfoL[:4]:
            outf.write("%f "%item[4])
        outf.write("\n")
        for item in evbinfoL[1:4]:
            outf.write("%f "%item[5])
        outf.write("\n")
        for item in evbinfoL[1:4]:
            ci=item[0]
            weight=0
            molA=item[2]
            for item in evbinfoL:
                if item[1]==molA:
                    weight=ci*item[0]
                    break
            outf.write("%f "%weight)
        outf.write("\n")
        outf.write("%s\n"%" ".join(map(str,CECcoords)))
        laststep=timestep
     

