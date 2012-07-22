import itertools
import pickle
import random

from copy import deepcopy

def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

# Data structure for storing families
class Family:
    def __init__(self, children, father, mother):
        self.children = children
        self.father = father
        self.mother = mother

# Use Population to generate simulated data for HS families
class Population:
    def __init__(self, numAlleles = 0, numLoci=0, numMothers=0,\
        numFathers=0, minFamilySize=1, maxFamilySize=1, randomSeed=289):

        self.randomSeed = randomSeed
        random.seed(randomSeed)

        self.usedParents = {}

        self.numAlleles = numAlleles
        self.numLoci = numLoci

        self.numFathers = numFathers
        self.numMothers = numMothers

        self.minFamilySize = minFamilySize
        self.maxFamilySize = maxFamilySize

        self.families = []
        self.individuals = []
        self.fathers = self.genParents(numFathers)
        self.mothers = self.genParents(numMothers)

        # Stores [father, mother] pair for each child
        # Key: tuple(child)
        # Value: [father, mother]
        self.parentDict = {}

        # Stores all the children for a given parent
        # Key: tuple(parent)
        # Value: list(children)
        self.childDictFKey = {}
        self.childDictMKey = {}

    def readPopulationFile(self, fn, nLoci):
        self.numLoci = nLoci

        f = open(fn, 'r')
        for line in f:
            tmp = line.split(',')

            self.individuals.append(\
                [int(tmp[x]) for x in range(1, len(tmp))])
        
        for ind in self.individuals:
            print ind

    # Add a family to the population
    def addFamily(self, children, father, mother, unique=False):
        # Do not allow duplicates
        if unique:
            children = self.uniqify(children)

        for child in children:
            # Assume duplicate individuals belong to the same father,mother
            if tuple(child) not in self.parentDict:
                self.parentDict[tuple(child)] = [father, mother]

        # Do not add empty families.  This can occur if all the children are
        # already in the population.
        if len(children) == 0:
            return

        # Update father and mother lookup dictionaries (for easy HS
        # family lookups)
        if tuple(father) not in self.childDictFKey:
            self.childDictFKey[tuple(father)] = children[:]
        else:
            self.childDictFKey[tuple(father)].extend(children[:])
        if tuple(mother) not in self.childDictMKey:
            self.childDictMKey[tuple(mother)] = children[:]
        else:
            self.childDictMKey[tuple(mother)].extend(children[:])
        
        # If the family already exists, just append the children
        for fam in self.families:
            if fam.father == father and fam.mother == mother:
                fam.children.extend(children[:])
                self.individuals.extend(children[:])

                return

        # If the family does not yet exist, add the children
        family = Family(children, father, mother)
        self.families.append(family)
        self.individuals.extend(children[:])

    # Generate parents for the population:
    # numParents fathers and numParents mothers
    def genParents(self, numParents):
        parents = {}

        # Make sure we can generate enough parents
        possibleParents = choose(self.numAlleles + 1, 2) ** self.numLoci
        if possibleParents < numParents:
            numParents = possibleParents

        i = 0
        while i < numParents - 1:
            p = []
            for j in range(self.numLoci):
                tmp = []
                tmp.append(random.randint(0, self.numAlleles - 1))
                tmp.append(random.randint(0, self.numAlleles - 1))
                tmp.sort()
                p.extend(tmp)

            # Make sure the parent does not already exist
            if not(tuple(p) in parents):
                parents[tuple(p)] = 1
                i += 1

        # Return a list of parents
        return map(lambda x: list(x), parents)

    # Generates children from the father and mother alleles
    # Number of children between minFamilySize and maxFamilySize
    def genChildren(self, numChildren, father, mother):
        children = []
        for i in range(0, numChildren):
            child = []
            for j in range(0, len(father), 2):

                fAllele = random.randint(0,1)
                mAllele = random.randint(0,1)

                child.append(father[j + fAllele])
                child.append(mother[j + mAllele])

            children.append(child)
        
        return children

    # Generates a family with numChildren children.  If father and mother
    # are not specified, a random father and mother will be selected, though
    # not recorded
    def genFamily(self, numChildren, father = [], mother = [],\
        oneSexMonog=False):
        fatherAvailable = [True] * len(self.fathers)

        invalid = True
        while invalid:
            if len(father) == 0:
                r = random.randint(0, len(self.fathers) - 1)
                if fatherAvailable[r]:
                    f = self.fathers[r]
                    fatherAvailable[r] = False
            else:
                f = father
            if len(mother) == 0:
                m = self.mothers[random.randint(0, len(self.mothers) - 1)]
            else:
                m = mother
        
            # Don't regenerate the same two parents
            if tuple([tuple(f), tuple(m)]) not in self.usedParents:
                self.usedParents[tuple([tuple(f), tuple(m)])] = True
                invalid = False
                break

        children = self.genChildren(numChildren, f, m)

        self.addFamily(children, f, m)

    # Creates a population with numIndiv distinct individuals.
    def genPopulation(self, numIndivs, oneSexMonog=False):
        self.mothers = self.genParents(self.numMothers)
        self.fathers = self.genParents(self.numFathers)

        # Generate the individuals
        while len(self.individuals) < numIndivs:
            # Random number of children within the specified range
            numChildren = random.randint(self.minFamilySize,\
                self.maxFamilySize)
            # Assure that we don't go over the number of children specified
            if numChildren + len(self.individuals) > numIndivs:
                numChildren = numIndivs - len(self.individuals)
            
            self.genFamily(numChildren, oneSexMonog=oneSexMonog)

    def printFSFams(self):
        for fam in self.families:
            tmp = ''
            for child in fam.children:
                tmp += "%4d" % self.individuals.index(child)
            print(tmp)

    def printHSFams(self):
        print("Parent Sex A:")
        for father in self.childDictFKey:
            tmp = ''
            for child in self.childDictFKey[father]:
                tmp += "%3d" % self.individuals.index(child)
            print(tmp)

        print('')

        print("Parent Sex B:")
        for mother in self.childDictMKey:
            tmp = ''
            for child in self.childDictMKey[mother]:
                tmp += "%3d" % self.individuals.index(child)
            print(tmp)

    def saveState(self, outfile):
        f = open(outfile, 'w')
        pickle.dump(self, f)

    # Removes duplicates from a list
    def uniqify(self, seq):
        seen = {} 
        result = [] 
        for item in seq: 
            if tuple(item) in seen: continue 
            seen[tuple(item)] = 1 
            result.append(item) 
        return result

    # True HS adjacency matrix
    def getAdjMtx(self):
        adjMtx = []
        for c0 in self.individuals:
            adjRow = []
            for c1 in self.individuals:
                p0 = self.parentDict[tuple(c0)]
                p1 = self.parentDict[tuple(c1)]
                if p0[0] == p1[0] or p0[1] == p1[1]:
                    adjRow.append(1)
                else:
                    adjRow.append(0)
            adjMtx.append(adjRow)

        return adjMtx

    # HS adjacency matrix with added random noise
    def getAdjMtxNoise(self, noiseLevel):
        adjMtx = []
        for c0 in self.individuals:
            adjRow = []
            for c1 in self.individuals:
                p0 = self.parentDict[tuple(c0)]
                p1 = self.parentDict[tuple(c1)]
                if p0[0] == p1[0] or p0[1] == p1[1]:
                    adjRow.append(1)
                elif random.random() < noiseLevel:
                    adjRow.append(1)
                else:
                    adjRow.append(0)
            adjMtx.append(adjRow)

        return adjMtx

    def isFS(self, c0, c1):
        if isinstance(c0, int):
            c0 = self.individuals[c0]
        if isinstance(c1, int):
            c1 = self.individuals[c1]

        p0 = self.parentDict[tuple(c0)]
        p1 = self.parentDict[tuple(c1)]

        return p0[0] == p1[0] and p0[1] == p1[1]

    # Parent is either 'P' for paternal or 'M' for maternal
    def isHS(self, c0, c1, parent=None):
        if isinstance(c0, int):
            c0 = self.individuals[c0]
        if isinstance(c1, int):
            c1 = self.individuals[c1]

        p0 = self.parentDict[tuple(c0)]
        p1 = self.parentDict[tuple(c1)]

        if parent is None:
            return p0[0] == p1[0] or p0[1] == p1[1]
        elif parent == 'P':
            return p0[0] == p1[0]
        elif parent == 'M':
            return p0[1] == p1[1]
        else:
            return False

    # Write a Colony 2 file
    def writeC2(self, title, mFS = False, pFS = False, randSeed=289,\
        updatingAlleleFreq=False, monoecious=False, diploid=True,\
        alleleFreq=[], sibSizePrior = False, nRuns = 1, lRuns = 2,\
        fullLikelihood=1, likelihoodPrecision = 3):

        f = open('Colony2.DAT', 'w')
        f.write("'%s'\n" % (title))
        f.write("'%s'\n" % (title))
        f.write("%d\n" % (len(self.individuals)))
        f.write("%d\n" % (self.numLoci))
        f.write("%d\n" % (randSeed))
        f.write("%d\n" % (int(updatingAlleleFreq)))

        if monoecious:
            f.write('1\n')
        else:
            f.write('2\n')

        d = 1
        if diploid:
            d = 0
        f.write("%d\n" % (d))
        # Maternal and paternal full = 1 or half = 0 sib.
        f.write("%d %d\n" % (int(pFS), int(mFS)))
        f.write("%d\n" % (int(sibSizePrior)))

        if alleleFreq == []:
            f.write("0\n")
        # Add allele frequency handling
        else:
            f.write("0\n")

        f.write("%d\n" % (nRuns))
        f.write("%d\n" % (lRuns))
        f.write("0\n")
        f.write("1000000\n")
        f.write("0\n")
        f.write("%d\n" % (fullLikelihood))
        f.write("%d\n" % (likelihoodPrecision))

        for i in range(self.numLoci):
            f.write("\tmk%d" % (i + 1))
        f.write("\n")

        for i in range(2):
            for j in range(self.numLoci):
                f.write("\t0")
            f.write("\n")
        for i in range(self.numLoci):
            f.write("\t0")
        f.write("\n")

        # Write individuals.  Lost alleles denoted as 0
        for i,ind in enumerate(self.individuals):
            f.write("%d " % (i))
            for allele in ind:
                f.write("%d " % (allele + 1))
            f.write("\n")

        f.write("0.0 0.0\n")
        f.write("0 0\n\n")
        f.write("0\n\n")
        f.write("0\n\n")
        f.write("0\n\n")
        f.write("0\n\n")
        f.write("0\n\n")
        f.write("0\n\n")
        f.write("0\n\n")
        f.write("0\n")
        f.close()

