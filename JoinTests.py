import itertools

import networkx as nx

from containers.Clusters import Clusters
from SJGlobals import SJGlobals

# Tests for decisions about joining individuals.
class JoinTests:
    @classmethod
    def isValidHS(self, clusterID0, clusterID1):
        allowable = SJGlobals.allowableJoins
        candidateParents = SJGlobals.candidateParents
        clusters = SJGlobals.clusters
        nLoci = SJGlobals.nLoci
        strictAlleles = SJGlobals.strictAlleles

        cluster0 = clusters.hsClusters[clusterID0]
        cluster1 = clusters.hsClusters[clusterID1]

        # Test that we're not combining the same individual
        for ind0 in cluster0.individuals:
            for ind1 in cluster1.individuals:
                if ind0.index == ind1.index:
                    return [False]

        # Test that any newly forced FS groups are valid FS
        fsToJoin = set()
        for ind0 in cluster0.individuals:
            for ind1 in cluster1.individuals:
                if ind0.getOppositeHSFam(clusterID0) ==\
                    ind1.getOppositeHSFam(clusterID1) and\
                    ind0.fsCluster != ind1.fsCluster:
                    
                    candidateCluster =\
                        clusters.fsClusters[ind0.fsCluster].individuals[:] +\
                        clusters.fsClusters[ind1.fsCluster].individuals[:]

                    if self.isValidFS(candidateCluster):
                        fsToJoin.add(tuple(sorted([\
                            ind0.fsCluster, ind1.fsCluster], reverse=True)))
                    else:
                        return [False]

        # Test that all individuals can be explained by at least one parent
        #   allele set at each locus
        parents = []
        if candidateParents:
            parents = cluster0.parentList & cluster1.parentList
            if len(parents) == 0:
                return [False]
        else:
            for l in range(nLoci):
                parents = cluster0.allParents[l] & cluster1.allParents[l]
                if len(parents) == 0:
                    return [False]

        # Test that we maintain bipartiteness
        if not clusters.biGraph.isValidCombining([[clusterID0, clusterID1]]):
            return [False]

        if strictAlleles and\
            not self.isHSAlleleCompat(clusterID0, clusterID1):

            return [False]

        return [True, sorted(fsToJoin, reverse = True)]

    @classmethod
    def isValidFS(self, candidateCluster):
        nLoci = SJGlobals.nLoci
        clusters = SJGlobals.clusters
        candidateParents = SJGlobals.candidateParents
        strictAlleles = SJGlobals.strictAlleles

        for l in range(nLoci):
            # Graph to test for full sibling compatibility on allele
            # graph (allele nodes share an edge if they share an indiv.)
            # Incompatible FS group if:
            #   1: # connected components > 2
            #   2: degree of any node > 3
            
            alleles = set()
            for ind in candidateCluster:
                alleles.add(ind.loci[l][0])
                alleles.add(ind.loci[l][1])

            # All clusters with 2 or less alleles are valid FS
            if len(alleles) < 3:
                continue
            if len(alleles) > 4:
                return False

            # Check all allele degrees < 3
            testG = nx.Graph()
            for ind in candidateCluster:
                if ind.loci[l][0] == ind.loci[l][1]:
                    testG.add_edge(str(ind.loci[l][0]),\
                        str(ind.loci[l][1]) + 'p')
                else:
                    testG.add_edge(str(ind.loci[l][0]), str(ind.loci[l][1]))

            for p0 in itertools.permutations(list(alleles), 2):
                for p1 in itertools.permutations(list(alleles), 2):
                    parentsFound = True
                    for i in candidateCluster:
                        # (0,0)
                        if (i.loci[l][0] == p0[0] and i.loci[l][1] == p1[0]) or\
                            (i.loci[l][0] == p1[0] and i.loci[l][1] == p0[0]):
                            continue
                        # (0,1)
                        elif (i.loci[l][0] == p0[0] and i.loci[l][1] == p1[1]) or\
                            (i.loci[l][0] == p1[1] and i.loci[l][1] == p0[0]):
                            continue
                        # (1,0)
                        elif (i.loci[l][0] == p0[1] and i.loci[l][1] == p1[0]) or\
                            (i.loci[l][0] == p1[0] and i.loci[l][1] == p0[1]):
                            continue
                        # (1,1)
                        elif (i.loci[l][0] == p0[1] and i.loci[l][1] == p1[1]) or\
                            (i.loci[l][0] == p1[1] and i.loci[l][1] == p0[1]):
                            continue
                        else:
                            parentsFound = False
                            break
                    if parentsFound:
                        break
                if parentsFound:
                    break
            if not parentsFound:
                return False

            # More than 2 connected components means more than 2 alleles
            # to explain
            if nx.number_connected_components(testG) > 2:
                return False

        return True
         
    @classmethod
    def isValidFSWithHS(self, ind0, ind1):
        clusters = SJGlobals.clusters

        # Test 1: Forms valid FS Group
        candidateCluster =\
            clusters.fsClusters[ind0.fsCluster].individuals[:] +\
            clusters.fsClusters[ind1.fsCluster].individuals[:]

        if not self.isValidFS(candidateCluster):
            return [False]

        origFS = [tuple(sorted([ind0.fsCluster, ind1.fsCluster], reverse=True))]

        # Test 2: The maternal and paternal clusters from both groups are
        #   each mutually compatible.  Either fams 0,0 and 1,1 can join or
        #   0,1 and 1,0 can join.  Both need to be checked.
        hs00 = ind0.hsClusters[0]
        hs01 = ind0.hsClusters[1]
        hs10 = ind1.hsClusters[0]
        hs11 = ind1.hsClusters[1]
        sharedCluster = ind0.sharesHSCluster(ind1)
        if sharedCluster[0]:
            clusterID0 = ind0.getOppositeHSFam(sharedCluster[1])
            clusterID1 = ind1.getOppositeHSFam(sharedCluster[1])


            vHS = self.isValidHS(clusterID0, clusterID1)
            if vHS[0]:
                fsToJoin = sorted(set(vHS[1] + origFS), reverse=True)
                return [True, [[clusterID0, clusterID1]], fsToJoin]

        vHS0 = self.isValidHS(hs00, hs10)
        vHS1 = self.isValidHS(hs01, hs11)
        if vHS0[0] and vHS1[0] and\
            clusters.biGraph.isValidCombining([[hs00, hs10],\
                [hs01, hs11]]):

            fsToJoin = sorted(set(vHS0[1] + vHS1[1] + origFS), reverse=True)

            # Check to make sure we're not joining a FS fam more than once
            if len(fsToJoin) > 1:
                for i in range(len(fsToJoin)):
                    for j in range(i + 1, len(fsToJoin)):
                        if fsToJoin[i][0] == fsToJoin[j][0] or\
                            fsToJoin[i][0] == fsToJoin[j][1] or\
                            fsToJoin[i][1] == fsToJoin[j][0] or\
                            fsToJoin[i][1] == fsToJoin[j][1]:

                            return [False]

            return [True, [[hs00, hs10], [hs01, hs11]], fsToJoin]

        # And the second option...
        vHS0 = self.isValidHS(hs00, hs11)
        vHS1 = self.isValidHS(hs01, hs10)
        if vHS0[0] and vHS1[0] and\
            clusters.biGraph.isValidCombining([[hs00, hs11],\
                [hs01, hs10]]):

            fsToJoin = sorted(set(vHS0[1] + vHS1[1] + origFS), reverse=True)
            
            # Check to make sure we're not joining a FS fam more than once
            if len(fsToJoin) > 1:
                for i in range(len(fsToJoin)):
                    for j in range(i + 1, len(fsToJoin)):
                        if fsToJoin[i][0] == fsToJoin[j][0] or\
                            fsToJoin[i][0] == fsToJoin[j][1] or\
                            fsToJoin[i][1] == fsToJoin[j][0] or\
                            fsToJoin[i][1] == fsToJoin[j][1]:

                            return [False]

            return [True, [[hs00, hs11], [hs01, hs10]], fsToJoin]

        return [False]

    # Given two individuals, return a list of tuples where each tuple
    #   explains a valid HS to join, 1 for each individual
    @classmethod
    def getCandidateHS(self, ind0, ind1):
        ret = []
        for clusterID0 in ind0.hsClusters:
            for clusterID1 in ind1.hsClusters:
                vHS = self.isValidHS(clusterID0, clusterID1)
                if vHS[0]:
                    ret.append([clusterID0, clusterID1])

        return ret

    '''
    # Tests to make sure that the proposed FS join would not force any
    #   allele incompatibilities.
    @classmethod
    def isFSAlleleCompat(self, cluster, parents0, parents1):
        nLoci = SJGlobals.nLoci
        candidateParents = SJGlobals.candidateParents

        for ind in cluster:
            for l in range(nLoci):
                allele = -1
                validAlleles0 = set()
                validAlleles1 = set()

                # Ignore homozygotes and allelic dropouts
                if ind.loci[l][0] == -1 or ind.loci[l][1] == -1 or\
                    ind.loci[l][0] == ind.loci[l][1]:

                    continue

                if candidateParents:
                    for p in parents0:
                        if ind.loci[l][0] == p.loci[l][0] or\
                            ind.loci[l][1] == p.loci[l][0]:

                            validAlleles0.add(p.loci[l][0])
                    
                        if ind.loci[l][0] == p.loci[l][1] or\
                            ind.loci[l][1] == p.loci[l][1]:

                            validAlleles0.add(p.loci[l][1])

                    for p in parents1:
                        if ind.loci[l][0] == p.loci[l][0] or\
                            ind.loci[l][1] == p.loci[l][0]:

                            validAlleles1.add(p.loci[l][0])
                    
                        if ind.loci[l][0] == p.loci[l][1] or\
                            ind.loci[l][1] == p.loci[l][1]:

                            validAlleles1.add(p.loci[l][1])
                else:
                    for p in parents0:
                        if ind.loci[l][0] == p[0] or ind.loci[l][1] == p[0]:
                            validAlleles0.add(p[0])
                        
                        if ind.loci[l][0] == p[1] or ind.loci[l][1] == p[1]:
                            validAlleles0.add(p[1])

                    for p in parents1:
                        if ind.loci[l][0] == p[0] or ind.loci[l][1] == p[0]:
                            validAlleles1.add(p[0])
                        
                        if ind.loci[l][0] == p[1] or ind.loci[l][1] == p[1]:
                            validAlleles1.add(p[1])

                # If both force the same allele, there is an incompat
                if len(validAlleles0) == 1 and\
                    len(validAlleles1) == 1 and\
                    validAlleles0.pop() == validAlleles1.pop():

                    return False

        return True
    '''

    @classmethod
    def isHSAlleleCompat(self, clusterID0, clusterID1):
        candidateParents = SJGlobals.candidateParents
        nLoci = SJGlobals.nLoci
        clusters = SJGlobals.clusters

        cluster0 = clusters.hsClusters[clusterID0]
        cluster1 = clusters.hsClusters[clusterID1]
        cluster = cluster0.individuals[:] + cluster1.individuals[:]
        parents = []
        if candidateParents:
            parents = cluster0.parentList & cluster1.parentList
        else:
            for l in range(nLoci):
                parents.append(cluster0.allParents[l] &\
                    cluster1.allParents[l])

        for ind in cluster:
            for l in range(nLoci):
                validAlleles = set()

                # Ignore homozygotes and allelic dropouts
                if ind.loci[l][0] == -1 or ind.loci[l][1] == -1 or\
                    ind.loci[l][0] == ind.loci[l][1]:

                    continue

                if candidateParents:
                    for p in parents:
                        if ind.loci[l][0] == p.loci[l][0] or\
                            ind.loci[l][1] == p.loci[l][0]:

                            validAlleles.add(p.loci[l][0])
                    
                        if ind.loci[l][0] == p.loci[l][1] or\
                            ind.loci[l][1] == p.loci[l][1]:

                            validAlleles.add(p.loci[l][1])
                else:
                    for p in parents[l]:
                        if ind.loci[l][0] == p[0] or ind.loci[l][1] == p[0]:
                            validAlleles.add(p[0])
                        
                        if ind.loci[l][0] == p[1] or ind.loci[l][1] == p[1]:
                            validAlleles.add(p[1])

                # Get the forced allele from the other family
                otherFam = 0
                if ind.hsClusters[0] == clusterID0 or\
                    ind.hsClusters[0] == clusterID1:
                    
                    otherFam = 1
                otherAllele = ind.forcedAlleles[l][otherFam]

                if otherAllele != -1 and len(validAlleles) == 1 and\
                    validAlleles.pop() == otherAllele:
                        
                    return False

        return True

