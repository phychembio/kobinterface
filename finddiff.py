inf1=open("sortdebugcoul.txt","r")
dataL1=[]
dataL2=[]
for line in inf1:
    dataL1.append(line.split()[0])

inf2=open("sortdebugcoulf90.txt","r")
for line in inf2:
    dataL2.append(line.split()[0])
	
for i in range(len(dataL1)):
    if dataL1[i]!= dataL2[i]:
	   print i, dataL1[i],dataL2[i]
	   break

