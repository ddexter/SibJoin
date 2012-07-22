import itertools
import networkx as nx

from SJGlobals import SJGlobals

'''
Keeps track of information about a cluster of individuals
Using candidate parents and enforcing strict alleles are optional
'''
class Cluster:
    '''
    Mandatory:
    clusterID is the integer index of the cluster in the clusters struct
    individuals is a list of individuals belonging to the cluster
    fsCluster specifies whether or not this cluster is a FS cluster

    Optional:
    parentList is a list of sampled candidate parents from a population
    '''
    def __init__(self, clusterID, individuals, candidateParentList=[],\
        fsCluster=False):

        self.clusterID = clusterID
        self.individuals = individuals
        self.fsCluster = fsCluster
        # self.parentList: Fully defined parents
        if not fsCluster:
            # Store the candidate parent list for future recomputations of
            # valid parents
            self.candidateParentList = candidateParentList
            self.initCompatibleParents(candidateParentList)

        self.setIndClusterID(self.individuals)
    
    def __len__(self):
        return len(self.individuals)

    def initCompatibleParents(self, candidateParentList):
        nLoci = SJGlobals.nLoci
        candidateParents = SJGlobals.candidateParents

        if candidateParents:
            # Find valid parents from candidateParents
            self.parentList = set()
            for p in candidateParentList:
                validParent = True
                for ind in self.individuals:
                    for l in range(nLoci):
                        if ind.loci[l][0] != p.loci[l][0] and\
                            ind.loci[l][0] != p.loci[l][1] and\
                            ind.loci[l][1] != p.loci[l][0] and\
                            ind.loci[l][1] != p.loci[l][1]:

                            validParent = False
                            break
                    if not validParent:
                        break
                if validParent:
                    self.parentList.add(p)

        self.allParents = []
        allAlleles = SJGlobals.allAlleles
        for l in range(nLoci):
            lociParents = set()
            alleleSet = set()

            # Add alleles for each individual in the cluster at given locus
            for ind in self.individuals:
                locus = ind.loci[l]

                # Do not add sites with allelic dropout
                if locus[0] != -1:
                    alleleSet.add(locus[0])
                if locus[1] != -1:
                    alleleSet.add(locus[1])

            otherAlleles = allAlleles[l] - alleleSet

            # Create list of all possible parents (all 2-permutations of
            # alleles)
            parents = [sorted(list(x))\
                for x in itertools.combinations(alleleSet,2)]
            for i in alleleSet:
                parents.append([i, i])
                for j in otherAlleles:
                    parents.append(sorted([i, j]))

            for p in parents:
                validParent = True
                for ind in self.individuals:
                    locus = ind.loci[l]

                    # Found incompatible parent
                    if len(set(locus) & set(p)) == 0:
                        validParent = False
                        break

                if validParent:
                    lociParents.add(tuple(p))

            self.allParents.append(lociParents)

    def computeForcedAlleles(self, ind):
        candidateParents = SJGlobals.candidateParents
        nLoci = SJGlobals.nLoci
        forcedAlleles = []

        for l in range(nLoci):
            allele = -1
            validAlleles = set()

            # Ignore homozygotes and allelic dropouts
            if ind.loci[l][0] == -1 or ind.loci[l][1] == -1 or\
                ind.loci[l][0] == ind.loci[l][1]:

                forcedAlleles.append(-1)
                continue

            # Check for forced allele from candidate parent set
            if candidateParents:
                for p in self.parentList:
                    if ind.loci[l][0] == p.loci[l][0] or\
                        ind.loci[l][1] == p.loci[l][0]:

                        validAlleles.add(p.loci[l][0])
                    
                    if ind.loci[l][0] == p.loci[l][1] or\
                        ind.loci[l][1] == p.loci[l][1]:

                        validAlleles.add(p.loci[l][1])
            # Check for forced allele from all parents set
            else:
                for p in self.allParents[l]:
                    if ind.loci[l][0] == p[0] or ind.loci[l][1] == p[0]:
                        validAlleles.add(p[0])
                    
                    if ind.loci[l][0] == p[1] or ind.loci[l][1] == p[1]:
                        validAlleles.add(p[1])

            if len(validAlleles) == 1:
                allele = validAlleles.pop()
            forcedAlleles.append(allele)

        return forcedAlleles

    '''
    Precondition: each individual in cluster must currently belong to at
       most 1 family.  Take care of destructing the old cluster first
    cluster must be a list
    '''
    def add(self, cluster, fs=False):
        nLoci = self.individuals[0].nLoci
        candidateParents = SJGlobals.candidateParents
        strictAlleles = SJGlobals.strictAlleles

        # Update compatible parent lists if not full-sibling clusters
        # (we don't keep track of FS potential parents)
        if not cluster.fsCluster:
            if candidateParents:
                self.parentList = self.parentList & cluster.parentList
            else:
                for l in range(nLoci):
                    self.allParents[l] = self.allParents[l] & cluster.allParents[l]

        # Update the cluster ID's of the individuals in the joined cluster
        self.setIndClusterID(cluster.individuals)

        self.individuals.extend(cluster.individuals)

    # indices must be a list
    def remove(self, indices):
        for i in indices:
            for ind in self.individuals:
                if ind.index == i:
                    self.individuals.remove(ind)

                    ind.removeCluster(self.clusterID, fs=self.fsCluster)

        # We must redefine compatible parents if this is a hs cluster
        if not self.fsCluster:
            self.initCompatibleParents(self.candidateParentList)

    def deleteCluster(self):
        for ind in self.individuals:
            ind.removeCluster(self.clusterID, fs=self.fsCluster)

    def updateClusterID(self, newClusterID):
        for ind in self.individuals:
            ind.updateCluster(self.clusterID, newClusterID,\
                fs=self.fsCluster)
        self.clusterID = newClusterID

    def setIndClusterID(self, individuals):
        strictAlleles = SJGlobals.strictAlleles

        for ind in individuals:
            if not self.fsCluster:
                forcedAlleles = []
                if strictAlleles:
                    forcedAlleles = self.computeForcedAlleles(ind)
                ind.setCluster(self.clusterID,\
                    forcedAlleles=forcedAlleles)
            else:
                ind.setCluster(self.clusterID, fs=True)

    def closenessScore(self, ind):
        nLoci = SJGlobals.nLoci
        
        score = 0.0
        for l in range(nLoci):
            locusScore = 0.0
            for p in list(self.allParents[l]):
                locusCounter = {}
                locusCounter.setdefault(ind.loci[l][0], 0)
                locusCounter.setdefault(ind.loci[l][1], 0)
                for i in self.individuals:
                    locusCounter.setdefault(i.loci[l][0], 0)
                    locusCounter.setdefault(i.loci[l][1], 0)

                    if p[0] == i.loci[l][0]:
                        if p[1] == i.loci[l][1]:
                            locusCounter[i.loci[l][0]] += 0.5
                            locusCounter[i.loci[l][1]] += 0.5
                        else:
                            locusCounter[i.loci[l][1]] += 1
                    elif p[0] == i.loci[l][1]:
                        if p[1] == i.loci[l][0]:
                            locusCounter[i.loci[l][0]] += 0.5
                            locusCounter[i.loci[l][1]] += 0.5
                        else:
                            locusCounter[i.loci[l][1]] += 1
                    elif p[1] == i.loci[l][0]:
                        locusCounter[i.loci[l][1]] += 1
                    else:
                        locusCounter[i.loci[l][0]] += 1

                if p[0] == ind.loci[l][0]:
                    if p[1] == ind.loci[l][1]:
                        locusScore += 0.5 * locusCounter[ind.loci[l][0]]
                        locusScore += 0.5 * locusCounter[ind.loci[l][1]]
                    else:
                        locusScore += locusCounter[ind.loci[l][1]]
                elif p[0] == ind.loci[l][1]:
                    if p[1] == ind.loci[l][0]:
                        locusScore += 0.5 * locusCounter[ind.loci[l][0]]
                        locusScore += 0.5 * locusCounter[ind.loci[l][1]]
                    else:
                        locusScore += locusCounter[ind.loci[l][0]]
                elif p[1] == ind.loci[l][0]:
                    locusScore += locusCounter[ind.loci[l][1]]
                else:
                    locusScore += locusCounter[ind.loci[l][0]]

            score += (float(locusScore) / float(len(self.allParents[l])))

        return score

    def avgLinkageDissimilarity(self, cluster):
        avg = 0.0
        for ind0 in self.individuals:
            for ind1 in cluster.individuals:
                avg += ind0.alleleDissimilarity(ind1)
        avg /= float(len(cluster.individuals) * len(self.individuals))

        return avg


