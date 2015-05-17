filename="debugcoul.txt"
inf=open(filename,"r")
outf=open("sort"+filename,"w")

dataL=[]
for line in inf:
    info=[float(a) for a in line.split()]	
    dataL.append(info)

for item in sorted(dataL):   
   outf.write("%f %f\n"%(item[0],item[1]))
inf.close()
outf.close()