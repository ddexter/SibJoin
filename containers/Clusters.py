import copy
import Queue
import networkx as nx

from SJGlobals import SJGlobals
from containers.Bipartite import Bipartite
from containers.Cluster import Cluster

class Clusters:
    def __init__(self, fsClusters, hsClusters):
        self.fsClusters = fsClusters
        self.hsClusters = hsClusters

        self.biGraph = Bipartite(hsClusters)

    # In each case, destruction of old cluster must occur first
    def join(self, clusterID0, clusterID1, fs=False):
        joinHistory = SJGlobals.joinHistory
        allowableJoins = SJGlobals.allowableJoins
        allowableClusterJoins = SJGlobals.allowableClusterJoins

        clusters = self.hsClusters
        jType = "hs"
        if fs:
            clusters = self.fsClusters
            jType = "fs"

            # None of the individuals can be joined again.  They're already
            # full-sibs
            for ind0 in clusters[clusterID0].individuals:
                for ind1 in clusters[clusterID1].individuals:
                    allowableJoins[ind0.index][ind1.index] = False
                    allowableJoins[ind1.index][ind0.index] = False
        else:
            self.biGraph.combineNodes(clusterID0, clusterID1)

            for i in range(2 * SJGlobals.nIndvs):
                allowableClusterJoins[i][clusterID1] = False
                if not allowableClusterJoins[clusterID1][i]:
                    allowableClusterJoins[clusterID0][i] = False
            
        clusters[clusterID1].deleteCluster()

        clusters[clusterID0].add(clusters[clusterID1])
        del clusters[clusterID1]

    def joinFS(self, fsClusterID0, fsClusterID1, hsClusterIDs):
        if len(hsClusterIDs) == 2:
            hsClusterIDs[0].sort()
            hsClusterIDs[1].sort()

            self.join(hsClusterIDs[0][0], hsClusterIDs[0][1])
            self.join(hsClusterIDs[1][0], hsClusterIDs[1][1])
        else:
            hsClusterIDs[0].sort()
            self.join(hsClusterIDs[0][0], hsClusterIDs[0][1])

    def sortMaternalPaternal(self):
        clusters = self.hsClusters

        '''
        Construct the algorithm's half-sibs.  We partition the clusters in
        to a maternal and paternal side 0 = Maternal, 1 = Paternal.  This
        part of the algorithm iterates through each cluster and adds the
        clusters to the respective 0 or 1 side based on the assumption that
        each indiv. will have a mother and a father.  Once a cluster has
        been added, the members of that cluster are queued up so that their 
        cluster representing the other parent can be added to the opposite
        partition.
        '''

        aHS = []
        addedClusters = [False] * len(clusters)
        numClusters = len(clusters)
        for i in clusters.keys():
            # Only process clusters which have not yet been assigned
            if not addedClusters[i]:
                q = Queue.Queue()
                tmp = [[],[]]

                tmp[0].append(copy.deepcopy(clusters[i]))
                addedClusters[i] = True

                # Stores a list of found individuals
                for ind in clusters[i].individuals:
                    q.put([ind, i, 0])

                # Add the other sex parents to the opposite partition
                while not q.empty():
                    cur = q.get()

                    # Find the cluster which has not just been added
                    if cur[1] == cur[0].hsClusters[0]:
                        j = cur[0].hsClusters[1]
                    else:
                        j = cur[0].hsClusters[0]

                    # Add the cluster and mark as added
                    pMIdx = (cur[2] + 1) % 2
                    if not addedClusters[j]:
                        tmp[pMIdx].append(copy.deepcopy(clusters[j]))
                        addedClusters[j] = True
                        
                        for ind in clusters[j].individuals:
                            q.put([ind, j, pMIdx])

                aHS.append(tmp)

        self.aHS = aHS

        mHS = []
        pHS = []
        for s in aHS:
            for i in s[0]:
                mHS.append(i)
            for i in s[1]:
                pHS.append(i)
        
        return mHS, pHS
    
    def matchingClusters(self, ind0, ind1):
        matchingClusters = []
        if ind0.hsClusters[0] == ind0.hsClusters[0]:
            matchingClusters.append([0, 0])
        if ind0.hsClusters[0] == ind0.hsClusters[1]:
            matchingClusters.append([0, 1])
        if ind0.hsClusters[1] == ind0.hsClusters[0]:
            matchingClusters.append([1, 0])
        if ind0.hsClusters[1] == ind0.hsClusters[1]:
            matchingClusters.append([1, 1])

        return matchingClusters

    def findFSInd(self, clusterID0, clusterID1):
        c0 = self.hsClusters[clusterID0]
        c1 = self.hsClusters[clusterID1]

        for ind0 in c0.individuals:
            for ind1 in c1.individuals:
                if ind0.index == ind1.index:
                    return ind0

        return None
