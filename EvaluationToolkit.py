import math
from Population import Population
from SJGlobals import SJGlobals

class EvaluationToolkit:
    def __init__(self, pop=[], parentInfo=True):
        if pop == [] or pop == None:
            return

        self.pop = pop
        self.parentInfo = parentInfo
        pop = self.pop

        # Create hash table for easy individual lookup
        self.indivIndices = {}
        for (i, indiv) in enumerate(pop.individuals):
            self.indivIndices[tuple(indiv)] = i

        if parentInfo:
            # Populate the real M and P HS arrays with child indices.
            self.realHS = []
            for i, mother in enumerate(pop.childDictMKey):
                self.realHS.append([[], mother])
                for child in pop.childDictMKey[mother]:
                    self.realHS[i][0].append(self.indivIndices[tuple(child)])
                self.realHS[i][0].sort()

            self.nMothers = len(self.realHS)
            for i, father in enumerate(pop.childDictFKey):
                self.realHS.append([[], father])
                for child in pop.childDictFKey[father]:
                    self.realHS[i + self.nMothers][0].append(\
                        self.indivIndices[tuple(child)])
                self.realHS[i + self.nMothers][0].sort()

    def vI(self, algCs, realCs, numIndvs):
        n = numIndvs

        # Compute the entropy of each algorithm cluster H(C)
        HC = 0.0
        for cluster in algCs:
            P_k = float(len(cluster)) / float(n)
            HC += (P_k * math.log(P_k, 2))
        HC *= -1.0

        HCp = 0
        for cluster in realCs:
            P_kp = float(len(cluster)) / float(n)
            HCp += (P_kp * math.log(P_kp, 2))
        HCp *= -1.0
        
        # Compute the mutual information I(C,C')
        I = 0
        for aCluster in algCs:
            for rCluster in realCs:
                aC = set(tuple(aCluster))
                rC = set(tuple(rCluster))
                P_k = float(len(aC)) / float(n)
                P_kp = float(len(rC)) / float(n)
                P_k_kp = float(len(aC & rC)) / float(n)

                if P_k_kp > 0:
                    I += P_k_kp * math.log(P_k_kp/(P_k * P_kp), 2)

        VI = HC + HCp - 2.0 * I

        return VI

    def computeVI(self, excludeInds = [], verbose=False):
        clusters = SJGlobals.clusters
        nIndvs = SJGlobals.nIndvs
        pop = self.pop
        indivIndices = self.indivIndices

        # Populate the real M and P HS arrays with child indices.
        rMHS = []
        for i, mother in enumerate(pop.childDictMKey):
            tmp = []
            for child in pop.childDictMKey[mother]:
                if indivIndices[tuple(child)] in excludeInds:
                    continue
                tmp.append(indivIndices[tuple(child)])
            if len(tmp) > 0:
                rMHS.append(tmp)
        numMothers = len(rMHS)

        rPHS = []
        for i, father in enumerate(pop.childDictFKey):
            tmp = []
            for child in pop.childDictFKey[father]:
                if indivIndices[tuple(child)] in excludeInds:
                    continue
                tmp.append(indivIndices[tuple(child)])
            if len(tmp) > 0:
                rPHS.append(tmp)
        numFathers = len(rPHS)

        partitionings = clusters.biGraph.getPartitionings()
        finalAlgM = []
        finalAlgP = []
        for part in partitionings:
            left = list(part[0])
            right = list(part[1])

            algM = []
            algP = []
            for i in left:
                algM.append([ind.index\
                    for ind in clusters.hsClusters[i].individuals])
            for i in right:
                algP.append([ind.index\
                    for ind in clusters.hsClusters[i].individuals])

            # Determine the partitioning assignment (maternal or paternal)
            # that minimizes the VI
            vi0 = self.vI(algM, rMHS, nIndvs) + self.vI(algP, rPHS, nIndvs)
            vi1 = self.vI(algP, rMHS, nIndvs) + self.vI(algM, rPHS, nIndvs)

            if vi0 < vi1:
                finalAlgM.extend(algM)
                finalAlgP.extend(algP)
            else:
                finalAlgM.extend(algP)
                finalAlgP.extend(algM)

        vi = self.vI(finalAlgM, rMHS, nIndvs) + self.vI(finalAlgP, rPHS, nIndvs)
        # 'normalize' the VI.  There are 2 copies of each individual.
        return vi / (2 * math.log(nIndvs, 2))
