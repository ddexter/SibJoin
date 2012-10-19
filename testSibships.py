import matplotlib
import math
import os
import pickle
import evaluationTools
import random
import re
import time

from EvaluationToolkit import EvaluationToolkit
from IPSolver import IPSolver
from Population import Population
from SibJoin import SibJoin

def generateTests():
    # Test IP
    # 6 alleles, 6 loci, popSize /  3 mothers and fathers, 5 children
    for i in range(10, 26, 5):
        for j in range(10):
            rSeed = 1399 * i + 931 * j
            pop = Population(6, 6, i / 3, i / 3, 5, 5, randomSeed=rSeed)
            pop.genPopulation(i)
            random.shuffle(pop.individuals)

            pop.saveState("tests/ip/%d_6_6_%d.pkl" % (i, j))

    # Test number of alleles
    # 6 loci, 7 mothers, 7 fathers, 5 children, 40 individuals
    alleles = [2, 5, 10, 15, 20]
    for i in alleles:
        for j in range(10):
            rSeed = 293 * i + 13 * j
            pop = Population(i, 6, 7, 7, 5, 5, randomSeed=rSeed)
            pop.genPopulation(40)
            random.shuffle(pop.individuals)

            if i < 10:
                nA = "0%d" % (i)
            else:
                nA = str(i)
            pop.saveState("tests/alleles/40_%s_6_%s.pkl" % (nA, str(j)))

    # Test number of loci
    # 6 alleles, 7 mothers, 7 fathers, 5 children, 40 individuals
    loci = [2, 5, 10, 15, 20]
    for i in loci:
        for j in range(10):
            rSeed = 631 * i + 17 * j
            pop = Population(6, i, 7, 7, 5, 5, randomSeed=rSeed)
            pop.genPopulation(40)
            random.shuffle(pop.individuals)

            if i < 10:
                nL = "0%d" % (i)
            else:
                nL = str(i)
            pop.saveState("tests/loci/40_6_%s_%d.pkl" % (nL, j))

    # Test number of individuals
    # 6 alleles, 6 loci, popSize /  3 mothers and fathers, 5 children
    popSizes = [10, 50, 100, 200]
    for i in popSizes:
        for j in range(10):
            rSeed = 37 * i + 931 * j
            pop = Population(6, 6, i / 3, i / 3, 5, 5, randomSeed=rSeed)
            pop.genPopulation(i)
            random.shuffle(pop.individuals)

            if i < 10:
                nI = "00%d" %False
            elif i < 100:
                nI = "0%d" % (i)
            else:
                nI = str(i)
            pop.saveState("tests/indivs/%s_6_6_%d.pkl" % (nI, j))

    # Test family size
    # 80 individuals, 6 alleles, 6 loci, 7 mothers and fathers
    famSizes = [1, 5, 10, 20]
    for i in famSizes:
        for j in range(10):
            rSeed = 31 * i + 1931 * j
            pop = Population(6, 6, 20, 20, i, i, randomSeed=rSeed)
            pop.genPopulation(80)
            random.shuffle(pop.individuals)

            if i < 10:
                nI = "0%d" % (i)
            else:
                nI = str(i)
            pop.saveState("tests/famsize/80_6_6_%s_%d.pkl" % (nI, j))

    # Test included families
    # 40 individuals, 6 alleles, 6 loci, 7 mothers and fathers
    for i in range(10):
        rSeed = 329019 * i
        pop = Population(6, 6, 8, 8, 5, 5, randomSeed=rSeed)
        pop.genPopulation(40)
        random.shuffle(pop.individuals)

        pop.saveState("tests/big/40_6_6_%d.pkl" % (i))

    popSizes = [500, 1000, 2000]
    for i in popSizes:
        for j in range(10):
            rSeed = 371 * i + 1931 * j
            pop = Population(6, 10, i / 20, i / 20, 20, 20, randomSeed=rSeed)
            pop.genPopulation(i)
            random.shuffle(pop.individuals)

            nI = ''
            if i < 1000:
                nI = "0%d" % (i)
            else:
                nI = str(i)
            pop.saveState("tests/big/%s_6_10_%d.pkl" % (nI, j))

def computeNumClusters(clusters):
    nClusters = 0
    for cluster in clusters:
        if cluster != []:
            nClusters += 1

    return nClusters

def computeSqClustSize(clusters):
    clustSum = 0
    for cluster in clusters:
        clustSum += len(cluster) * len(cluster)

    return clustSum

def computeInterClustDist(clusters, d):
    dist = 0.0
    nClust = 0
    for cluster in clusters:
        if len(cluster) == 0:
            continue

        nClust += 1

        clustAvg = 0.0
        nI = 0
        for i in range(len(cluster)):
            for j in range(i + 1, len(cluster)):
                nI += 1
                clustAvg += d[cluster[i]][cluster[j]]
        if nI > 0:
            clustAvg = float(clustAvg) / float(nI)
        dist += clustAvg
    
    return dist / float(nClust)

def computeLogClustSize(clusters):
    clustSum = 0
    for cluster in clusters:
        if len(cluster) > 0:
            clustSum += len(cluster) * math.log(len(cluster))

    return clustSum

def computeOptFn(clusters, d):
    return [computeInterClustDist(clusters, d),\
        computeNumClusters(clusters), computeSqClustSize(clusters),\
        computeLogClustSize(clusters)]

def test(fn, vi, verbose=False, outfile=''):

    vi = []
    time = []

    numTests = 10

    if outfile != '':
        out = open(outfile, 'a')

    path = 'tests/' + fn
    fileList = os.listdir(path)
    fileList.sort()
    res = []
    avgCnt = 0
    eRDiv = 0
    for i,test in enumerate(fileList):
        f = test.split('.')

        sj = SibJoin("pkl", fn='tests/%s%s.pkl' % (fn, f[0]))

        evaluator = EvaluationToolkit(sj.builder.pop)
        viRun = evaluator.computeVI()
        vi.append(viRun)
        time.append(sj.runTime)

        if out != '':
            out.write("%s%s %f %f\n" % (fn, f[0], viRun, sj.runTime))

            if (i + 1) % numTests == 0:
                aVI = sum(vi) / float(numTests)
                aTime = sum(time) / float(numTests)
                out.write("***%s%s %f %f\n" % (fn, f[0], aVI, aTime))

                vi = []
                time = []

        if verbose:
            print('%s%s: %f %f' %\
                (fn, f[0], viRun, sj.runTime))

    if verbose:
        print('')
        print("---------------")
        print('')

def testCOLONY():
    #testTypes = ['alleles/', 'loci/', 'indivs/', 'famsize/']
    testTypes = ['indivs/']

    ola = open("results/colonyResults.txt", "w")
    for testType in testTypes:
        path = 'tests/' + testType
        fileList = os.listdir(path)
        fileList.sort()
        t = []
        vi = []
        match = []
        numTests = 10

        for i,test in enumerate(fileList):
            '''
            TEST COLONY
            '''

            f2 = open("%s%s" % (path, test), "r")
            pop = pickle.load(f2)
            f2.close()
            if len(pop.individuals) <= 400:
                pop.writeC2('test')
                startTime = time.time()
                os.system('Colony2')
                stopTime = time.time()
            else:
                continue

            fin = open("test.BestConfig", "r")
            fin.readline()
            nIndvs = 0
            lines = []
            for line in fin.readlines():
                nIndvs += 1
                lines.append([x for x in re.split(r"[#\*\s]", line) if x])
            fin.close()

            mC = [[] for x in range(nIndvs + 1)]
            pC = [[] for x in range(nIndvs + 1)]

            for line in lines:
                pC[int(line[1])].append(int(line[0]))
                mC[int(line[2])].append(int(line[0]))

            clusters = []
            for cluster in mC:
                clusters.append(cluster)
            for cluster in pC:
                clusters.append(cluster)

            cPos = [[-1, -1] for x in range(nIndvs)]
            for i, cluster in enumerate(clusters):
                for ind in cluster:
                    if cPos[ind][0] == -1:
                        cPos[ind][0] = i
                    else:
                        cPos[ind][1] = i

            eT = evaluationTools.EvaluationTools(pop)
            correct = eT.computeMaxMatchings(clusters)
            wrong = 2 * nIndvs - correct
            pdNorm = float(wrong) / float(2 * nIndvs)
            viNorm = eT.compareResults3(clusters, cPos)

            vi.append(viNorm)
            match.append(pdNorm)
            t.append(stopTime - startTime)

            f = test.split('.')
            print("(COLONY)%s VI: %f PD: %f, Time: %f" % (f[0], viNorm, pdNorm, stopTime - startTime))
            ola.write("%s%s %f %f %f\n" % (testType, f[0], viNorm,\
                pdNorm, stopTime - startTime))
            ola.flush()
            print("here")
            if (i + 1) % numTests == 0:
                aVI = sum(vi) / float(numTests)
                aMatch = sum(match) / float(numTests)
                aTime = sum(t) / float(numTests)
                ola.write("***%s%s %f %f %f\n" % (testType, f[0], aVI,\
                    aMatch, aTime))
                ola.flush()
                vi = []
                match = []
                t = []
            
def testPaper():
    vi = []
    testDirs = ['alleles/', 'loci/', 'indivs/', 'famsize/', 'big/']
    #testDirs = ['alleles/', 'loci/']
    #testDirs = ['big/']
    
    for fn in testDirs:
        test(fn, vi, outfile='results/results.txt', verbose=True)

def testIP():
    path = "tests/ip"
    fileList = os.listdir(path)
    fileList.sort()
    t = []
    viSJ = []
    viIP = []
    out = open("results/ipResults.txt", "a")
    numTests = 10

    for i,test in enumerate(fileList):
        f = test.split('.')

        '''
        sj = SibJoin("pkl", fn='tests/%s%s.pkl' % ("ip/", f[0]))
        res = sj.getResults()
        viSJ.append(res[2])
        clusterings = sj.getClusterings()
        '''

        '''
        ipGuesses = [int(len(clusterings[0]) * 1.1),\
            int(len(clusterings[1]) * 1.1)]
        '''
        ip = IPSolver("pkl", fn="tests/%s%s.pkl" % ("ip/", f[0]),
            guesses=[0,0])
        res2 = ip.getResults()
        viIP.append(res2[1])
        t.append(res2[2])
        '''
        out.write("%s: %f %f %f\n" %\
            (f[0], res[2], res2[1], res2[2]))
        '''
        out.write("%s: %f %f\n" %\
            (f[0], res2[1], res2[2]))

        #print("%s: %f %f %f" % (f[0], res[2], res2[1], res2[2]))
        print("%s: %f %f" % (f[0], res2[1], res2[2]))
        
        if (i + 1) % 10 == 0:
            #sjVIMean = sum(viSJ) / 10.0
            ipVIMean = sum(viIP) / 10.0
            ipTimeMean = sum(t) / 10.0
            '''
            out.write("**%s: %f %f %f\n" %\
                (f[0], sjVIMean, ipVIMean, ipTimeMean))
            '''
            out.write("**%s: %f %f\n" %\
                (f[0], ipVIMean, ipTimeMean))

            t = []
            viSJ = []
            viIP = []

    out.close()

if __name__ == '__main__':
    #generateTests()
    #testCOLONY()
    testPaper()
    #testIP()

