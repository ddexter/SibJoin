import pickle

from containers.Cluster import Cluster
from containers.Clusters import Clusters
from containers.Individual import Individual
from SJGlobals import SJGlobals

class SibJoinBuilder:
    def __init__(self, filetype, fn=None, pop=None, strictAlleles=False):
        if (fn != None and pop != None) or (fn == None and pop == None):
            print("Error: Expected exactly one of (fn, pop) to be None")

        SJGlobals.clear()

        # Load test file
        named = False
        if filetype == "pkl":
            data = self.readPkl(fn)
        elif filetype == "txt":
            data = self.readText(fn)
            named = True
        else:
            data = self.readPop(fn)
        
        SJGlobals.nLoci = data[0]
        individuals = self.createIndvs(data[1], named=named)
        candidateParents = data[2]

        nLoci = data[0]

        # strictAlleles tells us whether to enforce that the same allele
        # can't come from 2 different parents.  Setting to true degrades
        # performance in practice
        SJGlobals.strictAlleles = strictAlleles

        # Set whether or not a candidate parent list is provided
        if len(candidateParents) > 0:
            SJGlobals.hasCandidateParents = True
            SJGlobals.candidateParents = candidateParents
        else:
            SJGlobals.candidateParents = False
            SJGlobals.candidateParents = []

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
        fsClusters = dict()
        hsClusters = dict()
        for i in range(SJGlobals.nIndvs):
            ind = individuals[i]
            fsClusters[i] = Cluster(i, [ind], fsCluster=True)
            hsClusters[2*i] = Cluster(2 * i, [ind])
            hsClusters[2*i+1] = Cluster(2 * i + 1, [ind])
        SJGlobals.clusters = Clusters(fsClusters, hsClusters)

        # Create the dissimilarity mtx and allow all initial joins
        self.d, self.threshold = self.createDissimilarityMtx(threshold=True)
        SJGlobals.allowableJoins = self.createAllowableJoins([])
        SJGlobals.allowableClusterJoins = self.createAllowableClusterJoins()

    def createAllowableClusterJoins(self):
        nIndvs = SJGlobals.nIndvs

        allowableClusterJoins = [[True]*2*nIndvs for i in range(2*nIndvs)]
        for i in range(0, 2 * nIndvs, 2):
            allowableClusterJoins[i][i+1] = False
            allowableClusterJoins[i+1][i] = False
            
            allowableClusterJoins[i][i] = False
            allowableClusterJoins[i+1][i+1] = False

        return allowableClusterJoins

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

        for cluster in clusters.items():
            for i in range(len(cluster[1].individuals)):
                for j in range(i + 1, len(cluster[1].individuals)):
                    ind0 = cluster[1].individuals[i]
                    ind1 = cluster[1].individuals[j]

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

    # The KINALYZER and COLONY formats accept names for each individual
    # If this is the case, the first column is the name, set named=True
    def createIndvs(self, inds, named=False):
        individuals = []
        nLoci = SJGlobals.nLoci

        for i in range(len(inds)):
            ind = []
            name = ''
            if named:
                name = inds[i][0]
                ind = inds[i][1:]
            else:
                ind = inds[i]
                
            loci = []

            for j in range(nLoci):
                loci.append([int(ind[2 * j]), int(ind[2 * j + 1]]))

            individuals.append(Individual(-1, [-1, -1], i, loci), name=name)

        return individuals

    def readPkl(self, fn):
        f = open(fn, 'r')
        pop = pickle.load(f)
        self.pop = pop
        nLoci = len(pop.individuals[0]) / 2
        
        return [nLoci, pop.individuals, []]
   
    def readPop(self, pop):
        ret = []
        self.pop = pop
        ret.append(pop.nLoci)
        ret.append(pop.individuals)
        ret.append([])

        return ret

    def readTxt(self, fn):
        f = open(fn, 'r')
        lines = f.readlines()
        
        individuals = [line.split(',') for line in lines]
        nLoci = len(individuals[0]) / 2

        return [nLoci, individuals, []]
    
    def setupResults(self):
        return self.d, self.threshold

