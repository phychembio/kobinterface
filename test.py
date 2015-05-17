#5/16/2015
import cPickle as pickle
import gzip
import sys
#str="3493493284903820583258032532"

#for i in range(1000000):
#    str=str+"3493493284903820583258032532"

#f = gzip.open("example", "wb", compresslevel=3)

#f.close()
f=gzip.open("surf.pickle","rb")

#xmin=ymin=0
#xmax=ymax=31.1
#zmin=14 
#zmax=30
#bin=0.5
#Nxbin=int((xmax-xmin)/bin)
#Nybin=int((ymax-ymin)/bin)
#Nzbin=int((zmax-zmin)/bin)
outf=open("test.xyz","w")
while True:
    try:
        value1=pickle.load(f)
        value2=pickle.load(f)
        print value1
        #info=[int(x) for x in value2.split()]
        info=[float(x) for x in value2.split()]
        outf.write("%d\n"%(len(info)/3))
        outf.write("Comments\n")
        for index,i in enumerate(info):
            if index%3==0:
                #x=xmin+ info[index]*0.5
                #y=ymin + info[index+1]*0.5
                #z=zmin+info[index+2]*0.5
                x=info[index]
                y=info[index+1]
                z=info[index+2]
                #if x==7 and y==3 and z==15:
                #    print "here"
                outf.write(" 1 %f %f %f\n"%(x,y,z))
        sys.exit()
    except EOFError:           
        break
    info=value2.split()

