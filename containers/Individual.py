'''
Keeps track of individual HS instances in a cluster.
'''
class Individual:
    '''
    Mandatory:
    hsClusters specifies which two half-sib clusters the individual is in
        by [cluster0, cluster1]
    fsCluster specifies which full-sib cluster the individual is in as int

    index is the given index number assigned to an individual
    loci is of the form [[a0,0, a0,1], [a1,0, a1,1] ...]

    Optional:
    forcedAlleles option if inhereted allele is forced.  Specified across
        all loci for each parent cluster.
        E.g. [[-1, -1], [-1, 2], [-1, -1]]
    '''
    def __init__(self, fsCluster, hsClusters, index, loci):
        self.fsCluster = fsCluster
        self.hsClusters = hsClusters
        self.index = index
        self.fsFalseJoin = False
        self.hsFalseJoins = [False, False]
        self.loci = loci
        self.nLoci = len(loci)

        self.forcedAlleles = []
        for l in range(self.nLoci):
            self.forcedAlleles.append([-1, -1])

    '''
    Compute the allele dissimilarity between self and another individual
    '''
    def alleleDissimilarity(self, ind):
        similarity = 0
        # Iterate over each locus
        for l in range(self.nLoci):
            s0 = 0
            if self.loci[l][0] == ind.loci[l][0] and self.loci[l][0] != -1:
                s0 += 1
            if self.loci[l][1] == ind.loci[l][1] and self.loci[l][1] != -1:
                s0 += 1
            
            s1 = 0
            if self.loci[l][0] == ind.loci[l][1] and self.loci[l][0] != -1:
                s1 += 1
            if self.loci[l][1] == ind.loci[l][0] and self.loci[l][1] != -1:
                s1 += 1

            similarity += max(s0, s1)

        return 2 * self.nLoci - similarity

    def hasForcedAlleles(self, locus, clusterID):
        # Determine the correct cluster
        c = 0
        if self.hsClusters[1] == clusterID:
            c = 1

        if self.forcedAlleles[locus][c] == -1:
            return False
        return True

    def getForcedAllele(self, locus, clusterID):
        c = 0
        if self.hsClusters[1] == clusterID:
            c = 1

        return self.forcedAlleles[locus][c]

    # Gets the forced allele from the opposite family
    def getOppositeForcedAllele(self, locus, clusterID):
        c = 0
        if self.hsClusters[0] == clusterID:
            c = 1

        return self.forcedAlleles[locus][c]

    def getOppositeHSFam(self, clusterID):
        if self.hsClusters[0] == clusterID:
            return self.hsClusters[1]
        else:
            return self.hsClusters[0]
   
    def setForcedAllele(self, allele, locus, parentCluster):
        c = 0
        if self.hsClusters[1] == parentCluster:
            c = 1

        self.forcedAlleles[locus][c] = allele

    def removeCluster(self, clusterID, fs=False):
        nLoci = self.nLoci

        if fs:
            self.fsCluster = -1
        else:
            c = 0
            if self.hsClusters[1] == clusterID:
                c = 1

            self.hsClusters[c] = -1

            for locus in range(nLoci):
                self.forcedAlleles[locus][c] = -1

    def setCluster(self, clusterID, forcedAlleles=[], fs=False):
        if fs:
            self.fsCluster = clusterID
        else:
            c = 0
            if self.hsClusters[1] == -1:
                c = 1

            self.hsClusters[c] = clusterID

            if len(forcedAlleles) > 0:
                for locus in range(self.nLoci):
                    self.forcedAlleles[locus][c] = forcedAlleles[locus]
   
    # Used to update the clusterID if a cluster is removed from the list
    # e.g. a merge leaves the second cluster empty
    def updateCluster(self, oldClusterID, newClusterID, fs=False):
        if fs:
            self.fsCluster = newClusterID
        else:
            c = 0
            if self.hsClusters[1] == oldClusterID:
                c = 1

            self.hsClusters[c] = newClusterID

    '''
    Used to determine if individual has already been joined to the wrong
    cluster.  This is mostly for testing purposes (learning fraction of
    false joins) rather than a part of the cluster reconstruct alg.
    '''
    def markFalseJoin(self, clusterID, fs=False):
        if fs:
            self.fsFalseJoin = True
        else:
            c = 0
            if self.hsClusters[1] == clusterID:
                c = 1

            self.hsFalseJoins[c] = True

    def hasFalseJoin(self, clusterID, fs=False):
        if fs:
            return self.fsFalseJoin
        else:
            c = 0
            if self.hsClusters[1] == clusterID:
                c = 1

            return self.hsFalseJoins[c]

    def sharesHSCluster(self, ind):
        if self.hsClusters[0] == ind.hsClusters[0] or\
            self.hsClusters[0] == ind.hsClusters[1]:

            return [True, self.hsClusters[0]]
        elif self.hsClusters[1] == ind.hsClusters[0] or\
            self.hsClusters[1] == ind.hsClusters[1]:

            return [True, self.hsClusters[1]]
        else:
            return [False]
