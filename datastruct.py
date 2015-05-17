#!/usr/bin/python

class hist2d:
    xbin = 0.1
    ybin = 0.1
    xmin = 0.0
    ymin = 0.0
    inc = 1.0
    dataL = []

    def __init__(self):
        print "Hist2d instance started"

    def adddata(self, x, y):
        L = self.dataL
        if x < self.xmin or y < self.ymin:
            exit("x or y is out of bound. x is %f and y is %f" % (x, y))
        else:
            whichxbin = int((x - self.xmin) / self.xbin)
            whichybin = int((y - self.ymin) / self.ybin)
            while len(L) <= whichxbin:
                if not L:
                    L.append([self.xmin + self.xbin / 2, []])
                else:
                    L.append([L[-1][0] + self.xbin, []])

            yL = L[whichxbin][1]
            while len(yL) <= whichybin:
                if not yL:
                    yL.append([self.ymin + self.ybin / 2, 0])
                else:
                    yL.append([yL[-1][0] + self.ybin, 0])
            L[whichxbin ][1][whichybin][1] += self.inc


    def output(self, filename):
        outf = open(filename, "w")
        L = self.dataL
        maxy = 0
        whichxindex = 0
        count = 0
        for item in L:
            if len(item[1]) > maxy:
                maxy = len(item[1])
                whichxindex = count
            count += 1
        labelL=[]
        for item in L[whichxindex][1]:
            labelL.append(item[0])
        outf.write("x, %s\n" %" ,".join(map(str,labelL)))
        for item in L:
            tempdataL = []
            for subitem in item[1]:
                tempdataL.append(subitem[1])
            while len(tempdataL)<maxy:
                tempdataL.append(0)
            outf.write("%f, %s\n" % (item[0], " ,".join(map(str,tempdataL))))

