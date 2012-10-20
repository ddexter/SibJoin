import copy
import matplotlib.pyplot as plt
import pickle
import Queue
import sys
import time

from containers.Bipartite import Bipartite
from containers.Cluster import Cluster
from containers.Clusters import Clusters
from containers.Individual import Individual
from EvaluationToolkit import EvaluationToolkit
from IPSolver import IPSolver
from JoinTests import JoinTests
from MinRemovals import MinRemovals
from OriginalIP import OriginalIP
from Population import Population
from SibJoinBuilder import SibJoinBuilder
from SJGlobals import SJGlobals

class SibJoin:
    def __init__(self, filetype, fn=None, pop=None):
        self.builder = SibJoinBuilder(filetype, fn=fn, pop=pop)
        self.d, self.threshold = self.builder.setupResults()

        self.startTime = time.time()

        self.run()

        self.stopTime = time.time()
        self.runTime = self.stopTime - self.startTime

        for i in SJGlobals.clusters.hsClusters.keys():
            print sorted([ind.index for ind in SJGlobals.clusters.hsClusters[i].individuals])

    def run(self):
        allowable = SJGlobals.allowableJoins
        d = self.d
        clusters = SJGlobals.clusters
        individuals = SJGlobals.individuals
        nIndvs = SJGlobals.nIndvs

        for t in range(self.threshold + 1):
            largestFirst = []
            for i in range(nIndvs):
                for j in range(i + 1, nIndvs):
                    if d[i][j] == t and allowable[i][j]:
                        ind0 = individuals[i]
                        ind1 = individuals[j]

                        # Join full-siblings first
                        vFS = JoinTests.isValidFSWithHS(ind0, ind1)
                        if vFS[0]:
                            clusters.joinFS(ind0.fsCluster, ind1.fsCluster,\
                                vFS[1])
                            if len(vFS) == 3:
                                for fsJoin in vFS[2]:
                                    clusters.join(fsJoin[1], fsJoin[0], fs=True)
                        # Construct half-sibs.  Prefer larger families first
                        else:
                            for k in range(2):
                                for l in range(2):
                                    cSize = min(\
                                        len(clusters.hsClusters[ind0.hsClusters[k]].individuals),\
                                        len(clusters.hsClusters[ind1.hsClusters[l]])\
                                    )
                                    largestFirst.append([cSize, i, j, k, l])
            largestFirst.sort(reverse=True)
            for fam in largestFirst:
                # Don't merge families which have already been merged or are
                # invalid
                if d[fam[1]][fam[2]] == -1 and allowable[fam[1]][fam[2]]:
                    continue

                ind0 = individuals[fam[1]]
                ind1 = individuals[fam[2]]
                clusterID0 = ind0.hsClusters[fam[3]]
                clusterID1 = ind1.hsClusters[fam[4]]
                vHS = JoinTests.isValidHS(clusterID0, clusterID1)
                if vHS[0]:
                    tmp = sorted([clusterID0, clusterID1])
                    clusters.join(tmp[0], tmp[1])
                    for fsJoin in vHS[1]:
                        # fsJoin[1] is actually the smaller index
                        clusters.join(fsJoin[1], fsJoin[0], fs=True)

    def getClusterings(self):
        return SJGlobals.clusters.sortMaternalPaternal()

if __name__ == '__main__':
    if len(sys.argv) > 1:
        sj = SibJoin("txt", sys.argv[1])
    #sj = SibJoin("pkl", fn="../tests/big/2000_6_10_7.pkl")
    #sj = SibJoin("pkl", fn="../tests/indivs/200_6_6_7.pkl")
    #sj = SibJoin("pkl", fn="../tests/alleles/40_20_6_0.pkl")
    #sj = SibJoin("pkl", fn="../tests/loci/40_6_05_5.pkl")
    #sj = SibJoin("pkl", fn="../tests/indivs/010_6_6_5.pkl")

