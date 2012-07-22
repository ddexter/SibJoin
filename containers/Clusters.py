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

        if SJGlobals.avgLinkage:
            self.createAvgLinkageMtx()
    
    # In each case, destruction of old cluster must occur first
    def join(self, clusterID0, clusterID1, fs=False):
        joinHistory = SJGlobals.joinHistory
        allowableJoins = SJGlobals.allowableJoins

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

        if SJGlobals.avgLinkage:
            self.updateAvgLinkageMtx(clusterID0, clusterID1, fs=fs)
            #print sorted([s.index for s in clusters[clusterID0].individuals])

        # Remove the old cluster
        # Save off indices to add to first cluster
        tmpCluster = clusters[clusterID1]
        self.remove(clusterID1, fs)

        clusters[clusterID0].add(tmpCluster)

        joinHistory.append(\
            [[ind.index for ind in clusters[clusterID0].individuals],\
            jType, clusterID0])

    def joinFS(self, fsClusterID0, fsClusterID1, hsClusterIDs):
        '''
        tmp = sorted([fsClusterID0, fsClusterID1])
        self.join(tmp[0], tmp[1], fs=True)
        '''

        # We need to account for the removal of a HS family in the indexing
        # e.g. if cluster 45 is removed and hsClusters[1][0] == 46, it will
        # be shifted to 45 by the automatic renumbering process
        if len(hsClusterIDs) == 2:
            hsClusterIDs[0].sort()
            hsClusterIDs[1].sort()
            if hsClusterIDs[1][0] > hsClusterIDs[0][1]:
                hsClusterIDs[1][0] -= 1
            if hsClusterIDs[1][1] > hsClusterIDs[0][1]:
                hsClusterIDs[1][1] -= 1

            hsClusterIDs[0].sort()
            hsClusterIDs[1].sort()
            self.join(hsClusterIDs[0][0], hsClusterIDs[0][1])
            self.join(hsClusterIDs[1][0], hsClusterIDs[1][1])
        else:
            hsClusterIDs[0].sort()
            self.join(hsClusterIDs[0][0], hsClusterIDs[0][1])

    def remove(self, clusterID, fs):
        clusters = self.hsClusters
        if fs:
            clusters = self.fsClusters

        clusters[clusterID].deleteCluster()

        # Shift each cluster number down by 1
        # Note: In the first draft, I wanted to avoid shifting after every
        #   join because it's super inefficient, but it makes the code
        #   illegible and fragile otherwise
        nodeNameMapping = {}
        for cID in range(clusterID + 1, len(clusters)):
            clusters[cID].updateClusterID(cID - 1)
            nodeNameMapping[cID] = cID - 1
        # Relabel nodes in bipartite graph
        if not fs:
            self.biGraph.relabelNodes(nodeNameMapping)

        # Remove the cID1 index from our distance matrix
        if SJGlobals.avgLinkage:
            d = self.hsAD
            if fs:
                d = self.fsAD

            for row in d:
                row.pop(clusterID)
            d.pop(clusterID)

        clusters.pop(clusterID)

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
            for i in range(numClusters):
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
    
    def createAvgLinkageMtx(self):
        self.fsAD = []
        self.hsAD = []

        # Create average linkage dissimilarity matrix for full-sibs
        for c0 in self.fsClusters:
            tmp = []
            for c1 in self.fsClusters:
                if c0.clusterID == c1.clusterID:
                    tmp.append(-1)
                else:
                    tmp.append(c0.avgLinkageDissimilarity(c1))
            self.fsAD.append(tmp)

        # Create average linkage dissimilarity matrix for half-sibs
        for c0 in self.hsClusters:
            tmp = []
            for c1 in self.hsClusters:
                if c0.clusterID == c1.clusterID:
                    tmp.append(-1)
                else:
                    tmp.append(c0.avgLinkageDissimilarity(c1))
            self.hsAD.append(tmp)

    '''
    Update the average linkage distances during a join.  cID1 will be
    removed from the matrix and all higher ID's will be shifted down
    '''
    def updateAvgLinkageMtx(self, cID0, cID1, fs=False):
        clusters = self.hsClusters
        d = self.hsAD
        if fs:
            clusters = self.fsClusters
            d = self.fsAD

        # Update the distances from each cluster to cID0
        for cluster in clusters:
            if cluster.clusterID == cID0 or cluster.clusterID == cID1:
                continue

            c0Len = float(len(clusters[cID0].individuals))
            c1Len = float(len(clusters[cID1].individuals))
            totIndvs = float(c0Len + c1Len)

            d[cID0][cluster.clusterID] =\
                (c0Len / totIndvs) * d[cID0][cluster.clusterID] +\
                (c1Len / totIndvs) * d[cID1][cluster.clusterID]
            d[cluster.clusterID][cID0] = d[cID0][cluster.clusterID]

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
