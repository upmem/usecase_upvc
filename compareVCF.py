#!/usr/bin/env python

import sys

def computeTPFN (V1,V2):
    XX = 10
    tp = 0
    for v in V2 :
        X = False
        for x in range(-XX,XX,1):
            if v+x in V1 :
                tp = tp+1
                X = True
                break;
        #if X == False:
        #    print ("FP:",v,V2[v])

    fn = 0
    for v in V1 :
        X = False
        for x in range(-XX,XX,1):
            if v+x in V2:
                X = True
                break
        if X == False:
            fn = fn+1
            #print ("FN:",v,V1[v][0:2])

    print ("tp :",str((tp*100)/len(V2))[0:5]+"%",tp,"/",len(V2))
    print ("fp :",str(100-(tp*100)/len(V2))[0:5]+"%",len(V2)-tp,"/",len(V2))
    print ("fn :",str((fn*100)/len(V1))[0:5]+"%",fn,"/",len(V1))
    print ("common match :",tp,"i.e",str((tp*100)/len(V1))[0:5]+"%","with reference")

def getData(filename):
    SUB = {}
    INS = {}
    DEL = {}
    ff = open(filename)
    line = ff.readline()
    while line != "":
        if line[0] != "#":
            L = line.split()
            pos = int(L[1])
            ref = L[3]
            alt = L[4]
            info = L[7]
            if len(ref) > 1 and len(ref) <= 5:
                DEL[pos] = [ref,alt,info]
            elif len(alt) > 1 and len(alt) <= 5:
                INS[pos] = [ref,alt,info]
            elif len(ref) == 1 and len(alt) == 1:
                SUB[pos] = [ref,alt,info]
        line = ff.readline()
    ff.close()
    return SUB, INS, DEL

print ("read",sys.argv[1])
SUB1, INS1, DEL1 = getData(sys.argv[1])
print (len(SUB1),len(INS1),len(DEL1))
print ("read",sys.argv[2])
SUB2, INS2, DEL2 = getData(sys.argv[2])
print (len(SUB2),len(INS2),len(DEL2))



print ("\nsubstitution")
computeTPFN(SUB1,SUB2)

print ("\ninsertions")
computeTPFN(INS1,INS2)

print ("\ndeletion")
computeTPFN(DEL1,DEL2)


