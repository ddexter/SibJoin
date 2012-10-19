import copy
import evaluationTools
import pickle
import sys
import time

from containers.Cluster import Cluster
from containers.Clusters import Clusters
from containers.Individual import Individual
from OriginalIP import OriginalIP
from Population import Population
from SJGlobals import SJGlobals

class IPSolver:
    def __init__(self, filetype, fn=None, pop=None, guesses=[0,0]):
        if (fn != None and pop != None) or (fn == None and pop == None):
            print("Error: Expected exactly one of (fn, pop) to be None")

        self.fn = fn
        SJGlobals.joinHistory = []
        # Load test file
        if filetype == "pkl":
            data = self.readPkl(fn)
        elif filetype == "txt":
            data = self.readText(fn)
        else:
            data = self.readPop(pop)

        self.startTime = time.time()

        SJGlobals.nLoci = data[0]
        individuals = self.createIndvs(data[1])
        SJGlobals.individuals = individuals
        SJGlobals.nIndvs = len(individuals)

        # Create a list of all alleles found at each locus in the pop.
        allAlleles = [set() for x in range(SJGlobals.nLoci)]
        for i in range(SJGlobals.nIndvs):
            ind = individuals[i]
            for l in range(SJGlobals.nLoci):
                allAlleles[l].add(ind.loci[l][0])
                allAlleles[l].add(ind.loci[l][1])
        SJGlobals.allAlleles = allAlleles

        # Guesses are the the number of clusters the IP should create in
        # each of the two partitionings
        for i in range(len(guesses)):
            if guesses[i] == 0:
                guesses[i] = len(individuals)

        self.ip = OriginalIP(individuals, guesses=guesses)
        self.stopTime = time.time()

    def createIndvs(self, inds):
        individuals = []
        nLoci = SJGlobals.nLoci

        for i in range(len(inds)):
            ind = inds[i]
            loci = []

            for j in range(nLoci):
                loci.append([ind[2 * j], ind[2 * j + 1]])

            individuals.append(Individual(-1, [-1, -1], i, loci))

        return individuals

    def readPkl(self, fn):
        f = open(fn, 'r')
        pop = pickle.load(f)
        self.pop = pop
        nLoci = len(pop.individuals[0]) / 2
        
        return [nLoci, pop.individuals, []]
   
    def readPop(self, pop):
        ret = []
        ret.append(pop.nLoci)
        ret.append(pop.individuals)
        ret.append([])

        return ret

    def readTxt(self, fn):
        f = open(fn, 'r')
        lines = f.readlines()

        nIndvs = int(lines[0])
        nLoci = int(lines[1])
        individuals = lines[2:nIndvs+3]
        candidateParents = lines[nIndvs+3:]

        return [nLoci, individuals, candidateParents]

    def getResults(self):
        eT = evaluationTools.EvaluationTools(self.pop)
        individuals = SJGlobals.individuals
        clusters = self.ip.getClusters()

        cHS = []
        for cluster in clusters.hsClusters:
            inds = [ind.index for ind in cluster.individuals]
            cHS.append(inds)
        cPos = [ind.hsClusters for ind in individuals]

        correct = eT.computeMaxMatchings(cHS)
        wrong = 2 * SJGlobals.nIndvs - correct
        vi = eT.compareResults3(cHS, cPos)
        timing = self.getTiming()

        return [wrong, vi, timing]

    def getTiming(self):
        return self.stopTime - self.startTime

    def getClusters(self):
        clusters = self.ip.getClusters()

        return clusters.sortMaternalPaternal()

if __name__ == "__main__":
    ips = IPSolver("pkl", fn="../tests/indivs/010_6_6_0.pkl")

    print ips.getResults()
