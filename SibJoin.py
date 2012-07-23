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
    def __init__(self, filetype, fn=None, pop=None, scoreSensitivity=2.0):
        if (fn != None and pop != None) or (fn == None and pop == None):
            print("Error: Expected exactly one of (fn, pop) to be None")

        self.fn = fn
        self.scoreSensitivity = scoreSensitivity
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
        candidateParents = data[2]

        self.setup(individuals, candidateParents = candidateParents,\
            strictAlleles=False)

        self.run()

        self.stopTime = time.time()
        self.runTime = self.stopTime - self.startTime

        self.calculateStats()

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

    def calculateStats(self):
        clusters = SJGlobals.clusters
        individuals = SJGlobals.individuals
        nIndvs = SJGlobals.nIndvs
        joinHistory = SJGlobals.joinHistory

        cHS = []
        for cluster in clusters.hsClusters:
            inds = [ind.index for ind in cluster.individuals]
            cHS.append(inds)

        eT = evaluationTools.EvaluationTools(self.pop)
        correct = eT.computeMaxMatchings(cHS)
        self.wrong = 2 * nIndvs - correct
        self.matchNorm = float(self.wrong) / float(2 * nIndvs)

        # Calculate variation of information
        cPos = [ind.hsClusters for ind in individuals]
        self.viNorm = eT.compareResults3(cHS, cPos)
        print self.viNorm

        x = [0]
        y = [0]
        cnt = 0
        nJoins = len(joinHistory)
        for i, join in enumerate(joinHistory):
            cID = join[2]
            if join[1] == "fs":
                continue
            nFP, a, b, bM, d = eT.findBestMatch(join[0])
            # Account for individuals which have already been involved in a 
            # bad join
            for ind in join[0]:
                if ind not in bM:
                    # Part of bad join, mark that it's already been involved
                    if not individuals[ind].hasFalseJoin(cID):
                        individuals[ind].markFalseJoin(cID)
                    # Already part of a bad join, don't count it
                    else:
                        nFP -= 1

            if nFP > 0:
                cnt += 1

            x.append(float(i + 1) / float(nJoins))
            y.append(cnt)

        '''
        fParts = self.fn.split("/")
        fPrefix = fParts[len(fParts) - 1].split(".")
        graphOut = "results/graphs/" + fPrefix[0] + ".png"
        plt.clf()
        plt.plot(x, y)
        plt.xlabel("% complete - HS joins")
        plt.ylabel("Total incorrect joins")
        plt.savefig(graphOut, format='png')
        '''

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
        nIndvs = SJGlobals.nIndvs

        return [self.wrong, 2 * nIndvs, self.viNorm, self.matchNorm, self.runTime]

    def removeScores(self, scores, sumOfScores, pairsToRemove):
        for pair in pairsToRemove:
            inds = sorted(pair)
            for i, entry in enumerate(scores):
                if inds[0].index == entry[1][0].index and\
                    inds[1].index == entry[1][1].index:

                    sumOfScores -= entry[0]
                    scores.pop(i)
                    break

        return scores, sumOfScores

    def runAnalytics(self, indsToRemove):
        allowableJoins = SJGlobals.allowableJoins
        clusters = SJGlobals.clusters
        nIndvs = SJGlobals.nIndvs

        eT = evaluationTools.EvaluationTools(self.pop)

        for ind in indsToRemove:
            fsCluster = ind.fsCluster
            clusters.fsClusters[fsCluster].remove([ind.index])
            if len(clusters.fsClusters[fsCluster]) == 0:
                clusters.remove(fsCluster, True)

            for i in range(2):
                hsCluster = ind.hsClusters[i]
                clusters.hsClusters[hsCluster].remove([ind.index])
                # Remove empty clusters
                if len(clusters.hsClusters[hsCluster]) == 0:
                    clusters.remove(hsCluster, False)

        newFSNum = len(clusters.fsClusters)
        newHSNum = len(clusters.hsClusters)
        for ind in indsToRemove:
            best = []
            tmpCluster = Cluster(newHSNum, [ind])
            for cluster in clusters.hsClusters:
                best.append([cluster.avgLinkageDissimilarity(tmpCluster),
                    cluster.clusterID])
            best.sort()

            print(ind.index, best[0][0],\
                sorted([s.index for s in\
                clusters.hsClusters[best[0][1]].individuals]))

            print(ind.index, best[1][0],\
                sorted([s.index for s in\
                clusters.hsClusters[best[1][1]].individuals]))

            eT.printRealClusters(ind.index)
            print("")


    def run(self):
        allowable = SJGlobals.allowableJoins
        d = self.d
        clusters = SJGlobals.clusters
        individuals = SJGlobals.individuals
        nIndvs = SJGlobals.nIndvs
        scores = self.scores
        sumOfScores = float(self.sumOfScores)

        random.seed(289)

        while len(scores) > 0:
            candidateJoin = random.randint(0, len(scores) - 1)
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

        ca = []
        for cluster in clusters.hsClusters:
            print sorted([ind.index for ind in cluster.individuals])
            ca.append(sorted([ind.index for ind in cluster.individuals]))
        print("")

if __name__ == '__main__':
    #sj = SibJoin("pkl", fn="../tests/indivs/050_6_6_7.pkl")
    #sj = SibJoin("pkl", fn="../tests/alleles/40_20_6_0.pkl")
    sj = SibJoin("pkl", fn="../tests/loci/40_6_05_5.pkl")

