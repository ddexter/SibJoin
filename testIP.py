import gc
import math
import os
import pickle
import evaluationTools
import random
import re
import time

from IPSolver import IPSolver
from Population import Population
from SibJoin import SibJoin

def generateTests():
    # Generate populations of 10 - 100 individuals in increments of 10
    # 6 alleles, 6 loci, popSize / 3 mothers and fathers, 5 children
    for i in range(10, 101, 10):
        rSeed = 37 * i + 931

        pop = Population(6, 6, i / 3, i / 3, 5, 5, randomSeed=rSeed)
        pop.genPopulation(i)
        random.shuffle(pop.individuals)

        if i < 100:
            nI = "0%d" % (i)
        else:
            nI = str(i)
        pop.saveState("../tests/IP/%s_6_6.pkl" % (nI))

def test():
    out = open("results/ip.txt", 'w')
    path = '../tests/IP/'

    fileList = os.listdir(path)
    fileList.sort()
    for test in fileList:
        gc.collect()

        f = test.split('.')

        sj = SibJoin("pkl", fn='../tests/%s%s.pkl' % ('IP/', f[0]))
        mHS, pHS = sj.getClusterings()
        ip = IPSolver("pkl", fn='../tests/%s%s.pkl' % ('IP/', f[0]),\
            guesses=[len(mHS) + 2, len(pHS) + 2])

        # File name, number wrong, VI, time
        sjRes = sj.getResults()
        out.write("SibJoin: %s %d %f %f"%\
            (f[0], sjRes[0], sjRes[2], sjRes[4]))
        print("SibJoin: %s %d %f %f"%\
            (f[0], sjRes[0], sjRes[2], sjRes[4]))

        ipRes = ip.getResults()
        out.write("IP: %s %d %f %f" %\
            (f[0], ipRes[0], ipRes[1], ipRes[2]))
        print("IP: %s %d %f %f" %\
            (f[0], ipRes[0], ipRes[1], ipRes[2]))
        
    close(out)
    
if __name__ == "__main__":
    test()
