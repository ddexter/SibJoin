import copy
import evaluationTools
import math
import pickle
import Queue
import random
import sys
import time

import matplotlib.pyplot as plt

from containers.Bipartite import Bipartite
from containers.Cluster import Cluster
from containers.Clusters import Clusters
from containers.Individual import Individual
from JoinTests import JoinTests
from MinRemovals import MinRemovals
from Population import Population
from SJGlobals import SJGlobals

class SibJoin:
    def __init__(self, filetype, fn=None, pop=None, scoreSensitivity=4.0,\
        rSeed=289):
        if (fn != None and pop != None) or (fn == None and pop == None):
            print("Error: Expected exactly one of (fn, pop) to be None")

        random.seed(rSeed)

        self.fn = fn
        self.scoreSensitivity = scoreSensitivity
        relatedScores = []

        # Load test file
        if filetype == "pkl":
            data = self.readPkl(fn)
        elif filetype == "txt":
            data = self.readText(fn)
        else:
            data = self.readPop(pop)

        for its in range(100):
            print its
            self.startTime = time.time()

            SJGlobals.nLoci = data[0]
            individuals = self.createIndvs(data[1])
            candidateParents = data[2]

            if relatedScores == []:
                relatedScores = [[0] * len(individuals)\
                    for i in range(len(individuals))]

            self.setup(individuals, candidateParents = candidateParents,\
                strictAlleles=False)

            self.run()

            self.stopTime = time.time()
            self.runTime = self.stopTime - self.startTime

            for cluster in SJGlobals.clusters.hsClusters:
                for i in range(len(cluster.individuals)):
                    for j in range(i + 1, len(cluster.individuals)):
                        ind0 = cluster.individuals[i].index
                        ind1 = cluster.individuals[j].index
                        relatedScores[ind0][ind1] += 1

            #self.resetSJGlobals()

        for row in relatedScores:
            print row


    def resetSJGlobals(self):
        SJGlobals.avgLinkage = False
        SJGlobals.candidateParents = []
        SJGlobals.isFS = []
        SJGlobals.isHS = []
        SJGlobals.joinHistory = []
        SJGlobals.nIndvs = -1
        SJGlobals.nLoci = -1
        SJGlobals.strictAlleles = False

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
        
    def setup(self, individuals, candidateParents=[], strictAlleles=False,\
        excludedIndvs=[]):

        nLoci = SJGlobals.nLoci
        SJGlobals.strictAlleles = strictAlleles
        if len(candidateParents) > 0:
            SJGlobals.candidateParents = True
        else:
            SJGlobals.candidateParents = False

        SJGlobals.individuals = individuals
        SJGlobals.nIndvs = len(individuals)

        # Create a list of all alleles found at each locus in the pop.
        allAlleles = [set() for x in range(nLoci)]
        for i in range(SJGlobals.nIndvs):
            ind = individuals[i]
            for l in range(nLoci):
                allAlleles[l].add(ind.loci[l][0])
                allAlleles[l].add(ind.loci[l][1])
        SJGlobals.allAlleles = allAlleles

        # Initialize the clusters
        fsClusters = []
        hsClusters = []
        for i in range(SJGlobals.nIndvs):
            ind = individuals[i]
            fsClusters.append(Cluster(i, [ind], fsCluster=True))
            hsClusters.append(Cluster(2 * i, [ind]))
            hsClusters.append(Cluster(2 * i + 1, [ind]))
        SJGlobals.clusters = Clusters(fsClusters, hsClusters)

        # Create the dissimilarity mtx and allow all initial joins
        self.d, self.threshold = self.createDissimilarityMtx(threshold=True)
        SJGlobals.allowableJoins = self.createAllowableJoins([])

        # Create a list of scores for each possible pairwise join
        d = self.d
        self.scores = []
        self.sumOfScores = 0.0
        maxScore = 2 * nLoci
        for i in range(SJGlobals.nIndvs):
            for j in range(i+1, SJGlobals.nIndvs):
                score = math.pow(math.e,\
                    self.scoreSensitivity * (maxScore - d[i][j]))
                self.sumOfScores += score
                self.scores.append([score, [individuals[i], individuals[j]]])
        # We want higher scores toward the end for efficiency since they are
        # more likely to be removed first.
        self.scores.sort()

    def setupRelationships(fs):
        nIndvs = SJGlobals.nIndvs

        famMtx = [[False] * nIndvs for i in range(nIndvs)]
        if fs:
            clusters = SJGlobals.clusters.fsClusters
        else:
            clusters = SJGlobals.clusters.hsClusters

        for cluster in clusters:
            for i in range(len(clusters.individuals)):
                for j in range(len(clusters.individuals)):
                    ind0 = clusters.individuals[i].index
                    ind1 = clusters.individuals[j].index

                    famMtx[ind0][ind1] = True
                    famMtx[ind1][ind0] = True

        if fs:
            SJGlobals.isFS = famMtx
        else:
            SJGlobals.isHS = famMtx

    def createAllowableJoins(self, excludedIndvs):
        nIndvs = SJGlobals.nIndvs
        allowableJoins = [[True] * nIndvs for i in range(nIndvs)]
        clusters = SJGlobals.clusters.hsClusters

        for i in range(nIndvs):
            allowableJoins[i][i] = False

        for i in excludedIndvs:
            for j in range(nIndvs):
                allowableJoins[i][j] = False
                allowableJoins[j][i] = False

        for cluster in clusters:
            for i in range(len(cluster.individuals)):
                for j in range(i + 1, len(cluster.individuals)):
                    ind0 = cluster.individuals[i]
                    ind1 = cluster.individuals[j]

                    allowableJoins[ind0.index][ind1.index] = False
                    allowableJoins[ind1.index][ind0.index] = False

        return allowableJoins

    def createDissimilarityMtx(self, threshold=False):
        d = []
        threshold = -1
        nIndvs = SJGlobals.nIndvs
        individuals = SJGlobals.individuals

        for i in range(nIndvs):
            tmp = []
            ind0 = individuals[i]

            for j in range(nIndvs):
                if i == j:
                    tmp.append(-1)
                else:
                    dis = ind0.alleleDissimilarity(individuals[j])
                    tmp.append(dis)

                    if dis > threshold:
                        threshold = dis

            d.append(tmp)

        if threshold:
            return d, threshold
        else:
            return d

    def getResults(self):
        return [self.wrong, 2 * SJGlobals.nIndvs, self.viNorm,\
            self.matchNorm, self.runTime]

    def removeScores(self, scores, sumOfScores, pairsToRemove):
        for pair in pairsToRemove:
            inds = sorted([pair[0].index, pair[1].index])
            for i, entry in enumerate(scores):
                if inds[0] == entry[1][0].index and\
                    inds[1] == entry[1][1].index:

                    sumOfScores -= entry[0]
                    scores.pop(i)
                    break

        return scores, sumOfScores

    def run(self):
        allowable = SJGlobals.allowableJoins
        d = self.d
        clusters = SJGlobals.clusters
        individuals = SJGlobals.individuals
        nIndvs = SJGlobals.nIndvs
        scores = self.scores
        sumOfScores = float(self.sumOfScores)

        while len(scores) > 0:
            candidateJoin = random.randint(0, len(scores) - 1)
            tmp = scores[candidateJoin][0] / sumOfScores
            if random.random() < scores[candidateJoin][0] / sumOfScores:
                ind0 = scores[candidateJoin][1][0]
                ind1 = scores[candidateJoin][1][1]

                pairsToRemove = set()

                # Try joining full-siblings first
                vFS = JoinTests.isValidFSWithHS(ind0, ind1)
                if vFS[0]:
                    # Join initial half-sibs
                    clusters.joinFS(ind0.fsCluster, ind1.fsCluster, vFS[1])

                    # Make any other full-sib joins that were made necessary
                    # by the half-sib joins
                    if len(vFS) == 3:
                        for fsJoin in vFS[2]:
                            pairsToRemove = clusters.join(\
                                fsJoin[1], fsJoin[0], fs=True)
                            scores, sumOfScores = self.removeScores(\
                                scores, sumOfScores, pairsToRemove)
                # If FS doesn't work, try HS
                else:
                    hs = JoinTests.getCandidateHS(ind0, ind1)
                    hsJoin = []
                    if len(hs) > 0:
                        hsJoin = hs[random.randint(0, len(hs) - 1)]

                        clusterID0 = hsJoin[0]
                        clusterID1 = hsJoin[1]
                        vHS = JoinTests.isValidHS(clusterID0, clusterID1)
                        if vHS[0]:
                            tmp = sorted([clusterID0, clusterID1])
                            clusters.join(tmp[0], tmp[1])
                            for fsJoin in vHS[1]:
                                # fsJoin[1] is actually the smaller index
                                pairsToRemove = clusters.join(\
                                    fsJoin[1], fsJoin[0], fs=True)
                                scores, sumOfScores = self.removeScores(\
                                    scores, sumOfScores, pairsToRemove)

                scores, sumOfScores = self.removeScores(\
                    scores, sumOfScores, set([tuple([ind0, ind1])]))

        self.ca = []
        for cluster in clusters.hsClusters:
            tmp = sorted([ind.index for ind in cluster.individuals])
            self.ca.append(tmp)

if __name__ == '__main__':
    #sj = SibJoin("pkl", fn="../tests/indivs/050_6_6_7.pkl")
    #sj = SibJoin("pkl", fn="../tests/alleles/40_20_6_0.pkl")
    sj = SibJoin("pkl", fn="../tests/loci/40_6_05_5.pkl")

