import re

def findstart(trajf):
    trajf.next()
    return int(trajf.next())

inf=open("OH-interface/out.cleanevb","r")
trajf=open("OH-interface/out.lammpstrj","r")
outf=open("offdiag.txt","w")
mintime=findstart(trajf)
sum=0
count=0
for line in inf:
     match = re.match("--",line)
     if match:
        evbtime = int(inf.next())        
        evbciL = [float(x) for x in inf.next().split()]
        evbmolBL = [int(x) for x in inf.next().split()]
        evbmolAL = [int(x) for x in inf.next().split()]
        evbshellL = [int(x) for x in inf.next().split()]
        evbrepL = [float(x) for x in inf.next().split()]
        evboffdiagL = [float(x) for x in inf.next().split()]
        evbcijL = [float(x) for x in inf.next().split()]
        CEC = [float(x) for x in inf.next().split()]
        if evbtime>=mintime:
            for i in range(len(evboffdiagL)):
                sum+=evboffdiagL[i]*evbcijL[i]
            count+=1
outf.write("%d %f\n"%(count, sum/count))
inf.close()
outf.close()
trajf.close()