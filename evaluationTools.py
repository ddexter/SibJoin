import copy
import itertools
import matching
import math
import pydot
import Queue

import matplotlib.pyplot as plt

class EvaluationTools:
    def __init__(self, pop=[], parentInfo=True):
        if pop == []:
            return

        self.pop = pop

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
    
    # Preference least # false positives
    def findBestMatch(self, cluster):
        pop = self.pop
        indivIndices = self.indivIndices
        realHS = copy.deepcopy(self.realHS)

        # Total individuals wrong
        totalWrong = 0

        # Total missing individuals
        totalMiss = 0

        bestMatch = -1
        bMCount = -1
        nFalse = -1
        nMatch = -1
        nMiss = -1
        parent = []

        algMap = [False] * len(pop.individuals)

        for child in cluster:
            algMap[child] = True
       
        for i, realFam in enumerate(realHS):
            realMap = [False] * len(pop.individuals)

            for child in realFam[0]:
                realMap[child] = True

            # Count the number of children found by the alg, in the real
            # family
            numMatch = 0
            for child in cluster:
                if realMap[child]:
                    numMatch += 1

            # Count the number of children found by the alg, but not in
            # a real family
            numFalse = 0
            for child in cluster:
                if not realMap[child]:
                    numFalse += 1

            # Count the number of children excluded by the algorithm
            numMissing = 0
            for child in realFam[0]:
                if not algMap[child]:
                    numMissing += 1
            
            if numMatch > bMCount:
                nFalse = numFalse
                nMatch = numMatch
                nMiss = numMissing
                bMCount = numMatch
                bestMatch = realFam[0]
                parent = realFam[1]

            # Break ties with number missing
            if numMatch == bMCount and nMiss > numMissing:
                nFalse = numFalse
                nMatch = numMatch
                nMiss = numMissing
                bMCount = numMatch
                bestMatch = realFam[0]
                parent = realFam[1]

        bestMatch.sort()

        return nFalse, nMatch, nMiss, bestMatch, parent

    def compareResults(self, clusters, verbose=False, special=[]):
        """
        Measure different.  The most different are calculated as:

        numMatches = # children in algorithm's HS fam that match the real
            HS family
        """
        pop = self.pop
        nIndvs = len(pop.individuals)
        indWrong = [0] * nIndvs

        totalWrong = 0
        totalMiss = 0
        totalMatch = 0
        for algFam in clusters:
            if len(algFam) < 1:
                continue

            ret = self.findBestMatch(algFam)
            totalWrong += ret[0]
            totalMiss += ret[2]
            realFam = ret[3]
            for ind in algFam:
                if ind not in realFam:
                    indWrong[ind] += 1
            if verbose:
                print(ret[0], ret[1], ret[2])
                algFam.sort()
                print(algFam)
                print(ret[3])
                print('')

        res = [0, 0, 0 ,0]
        for ind in special:
            if indWrong[ind] == 0:
                res[0] += 1
                continue
            elif indWrong[ind] == 1:
                res[3] += 1
                res[1] += 1
            elif indWrong[ind] == 2:
                res[3] += 1
                res[2] += 1
        res.append(len(special))
        res.append(totalWrong)

        if verbose:
            print("Total Wrong: {0}".format(totalWrong))
            print("Total Missing: {0}".format(totalMiss))
            
            for ind in special:
                print("{0}: {1}".format(ind, indWrong[ind]))
                
        return res

    def drawResultsGraphError(self, groupings, fn, fs=False):
        pop = self.pop
        g = pydot.graph_from_dot_file('graphs/%s.dot'%(fn))

        for group in groupings:
            for i in range(0, len(group)):
                for j in range(i+1, len(group)):
                    color = 'black'
                    if pop.isFS(group[i], group[j]):
                        color = 'black'
                    elif pop.isHS(group[i], group[j], parent='P'):
                        color = 'blue'
                    elif pop.isHS(group[i], group[j], parent='M'):
                        color = 'green'
                    else:
                        color = 'red'

                    if color == 'red':
                        g.add_edge(pydot.Edge(str(group[i]), str(group[j]),\
                            color=color))

        if fs:
            graphPath = 'results/%s_FS_Error.svg'%(fn)
        else:
            graphPath = 'results/%s_Error.svg'%(fn)
        g.write_svg(graphPath, prog=['neato', '-n'])

    def printDistanceMatrix(self, d):
        tmp = '%3c' % ' '
        for i in range(0, len(d)):
            tmp += "%3d" % i
        print(tmp)

        for i in range(0, len(d)):
            tmp = "%3d" % i
            for j in range(0, len(d[i])):
                tmp += "%3d" % d[i][j]
            print(tmp)

    def drawResultsGraph(self, groupings, pop, fn, fs=False, altPath = ''):
        numIndiv = len(pop.individuals)
        g = pydot.graph_from_dot_file('graphs/%s.dot'%(fn))
        g.set_strict(True)

        for group in groupings:
            for i in range(0, len(group)):
                for j in range(i+1, len(group)):
                    color = 'black'
                    if pop.isFS(group[i], group[j]):
                        color = 'black'
                    elif pop.isHS(group[i], group[j], parent='P'):
                        color = 'blue'
                    elif pop.isHS(group[i], group[j], parent='M'):
                        color = 'green'
                    else:
                        color = 'red'

                    g.add_edge(pydot.Edge(str(group[i]), str(group[j]),\
                        color=color))

        if fs:
            graphPath = 'results/%s_FS.svg'%(fn)
        elif altPath == '':
            graphPath = 'results/%s.svg'%(fn)
        else:
            graphPath = altPath

        g.write_svg(graphPath, prog=['neato', '-n'])

    def printDistanceMatrix(self, d):
        tmp = '%3c' % ' '
        for i in range(0, len(d)):
            tmp += "%3d" % i
        print(tmp)

        for i in range(0, len(d)):
            tmp = "%3d" % i
            for j in range(0, len(d[i])):
                tmp += "%3d" % d[i][j]
            print(tmp)

    def printDistanceMatrixToFile(self, d, f):
        tmp = '%3c' % ' '
        for i in range(0, len(d)):
            tmp += "%3d" % i
        f.write(tmp + "\n")

        for i in range(0, len(d)):
            tmp = "%3d" % i
            for j in range(0, len(d[i])):
                tmp += "%3d" % d[i][j]
            f.write(tmp + "\n")

    def computeMaxMatchings(self, clusters):
        pop = self.pop
        tmp = copy.deepcopy(self.realHS)
        indivIndices = self.indivIndices
        realHS = [x[0] for x in tmp]

        algHS = []
        for cluster in clusters:
            if len(cluster) > 0:
                algHS.append(cluster)

        # Create hash table for easy individual lookup
        indivIndices = {}
        for (i, indiv) in enumerate(pop.individuals):
            indivIndices[tuple(indiv)] = i

        correct = 0
        for i in algHS:
            i.sort()
            if i in realHS:
                correct += len(i)
                algHS.remove(i)
                realHS.remove(i)

        edges = []
        lClusters = len(algHS)
        lRealHS = len(realHS)
        for i in range(lClusters):
            for j in range(lRealHS):
                c = 0
                for p in algHS[i]:
                    if p in realHS[j]:
                        c += 1
                edges.append((i, j + lClusters, c))
        mate = matching.maxWeightMatching(edges)
        for i in range(lClusters):
            if mate[i] >= 0:
                for p in algHS[i]:
                    if p in realHS[mate[i] - lClusters]:
                        correct += 1

        return correct

    def computeMaxMatchingsFS(self, clusters):
        self.pop = pop
        algFS = []
        for cluster in clusters:
            if len(cluster) > 0:
                algFS.append(cluster)

        # Create hash table for easy individual lookup
        indivIndices = {}
        for (i, indiv) in enumerate(pop.individuals):
            indivIndices[tuple(indiv)] = i

        realFS = []
        for fam in pop.families:
            realFS.append([indivIndices[tuple(child)]\
                for child in fam.children])

        correct = 0
        for i in algFS:
            if i in realFS:
                correct += len(i)
                algFS.remove(i)
                realFS.remove(i)

        edges = []
        lClusters = len(algFS)
        lRealFS = len(realFS)
        for i in range(lClusters):
            for j in range(lRealFS):
                c = 0
                for p in algFS[i]:
                    if p in realFS[j]:
                        c += 1
                edges.append((i, j + lClusters, c))
        mate = matching.maxWeightMatching(edges)
        for i in range(lClusters):
            if mate[i] >= 0:
                for p in algFS[i]:
                    if p in realFS[mate[i] - lClusters]:
                        correct += 1

        return correct

    # x and y are lists of lists.  Each sublist contains the x and y
    # coordinates respectively that should be printed on a single graph
    def plotResults(self, x, y, destination, labels, title, xLabel, yLabel,\
        errorBars=False, lines=True):
        plt.clf()
        ax = plt.axes()

        l = [0] * len(x)
        for i in range(len(x)):
            if lines:
                l[i], = ax.plot(x[i], y[i], label=labels[i])
            else:
                l[i], = ax.plot(x[i], y[i], 'o', label=labels[i])


        handles, labels = ax.get_legend_handles_labels()

        ax.legend(handles, labels)
        plt.xlabel(xLabel)
        plt.ylabel(yLabel)
        plt.title(title)

        plt.savefig(destination, format='png')

    def computeVI(self, algCs, realCs, numIndvs):
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

# I need a rewrite, inefficient, but works
    def compareResults3(self, cS, locations, excludeInds = [], verbose=False):
        pop = self.pop
        indivIndices = self.indivIndices
        clusters = copy.deepcopy(cS)
        numIndvs = len(pop.individuals)

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

        # Construct the algorithm's half-sibs.  We partition the clusters into
        # a maternal and paternal side 0 = Maternal, 1 = Paternal.  This part of
        # the algorithm iterates through each cluster and adds the clusters to
        # the respective 0 or 1 side based on the assumption that each indiv.
        # will have a mother and a father.  Once a cluster has been added, the
        # members of that cluster are queued up so that their cluster
        # representing the other parent can be added to the opposite partition.
        aHS = []
        addedClusters = [False] * len(clusters)
        # Don't process the empty clusters
        for i,cluster in enumerate(clusters):
            if len(cluster) == 0:
                addedClusters[i] = True

        numClusters = len(clusters)
        for i in range(numClusters):
            # Only bother processing clusters which have not yet been assigned
            if not addedClusters[i]:
                q = Queue.Queue()
                tmp = [[],[]]

                tmp[0].append(copy.deepcopy(clusters[i]))
                addedClusters[i] = True

                # Stores a list of found individuals
                for ind in clusters[i]:
                    q.put([ind, i, 0])

                # Add the other sex parents to the opposite partition
                while not q.empty():
                    cur = q.get()

                    # Find the cluster which has not just been added
                    if cur[1] == locations[cur[0]][0]:
                        j = locations[cur[0]][1]
                    else:
                        j = locations[cur[0]][0]

                    # Add the cluster and mark as added
                    pMIdx = (cur[2] + 1) % 2
                    if not addedClusters[j]:
                        tmp[pMIdx].append(copy.deepcopy(clusters[j]))
                        addedClusters[j] = True
                        
                        for ind in clusters[j]:
                            q.put([ind, j, pMIdx])

                aHS.append(tmp)

        # Fix duplicates by re-adding them as individual clusters
        for s in aHS:
            # s contains [[paternal],[maternal]] half-sib groups
            for si in s:
                duplicate = [False] * len(pop.individuals)
                for sj in si:
                    for ind in sj:
                        if duplicate[ind]:
                            # Remove the duplicate entry
                            clusters[locations[ind][0]].remove(ind)
                            clusters[locations[ind][1]].remove(ind)

                            # Add the duplicate, 1 to each side
                            clusters.append([ind])
                            clusters.append([ind])

                            # Update the reference
                            locations[ind][0] = len(clusters) - 2
                            locations[ind][0] = len(clusters) - 1
                                
                            duplicate[ind] = False
                        else:
                            duplicate[ind] = True

        # Construct the algorithm's half-sibs.  We partition the clusters into
        # a maternal and paternal side 0 = Maternal, 1 = Paternal.  This part of
        # the algorithm iterates through each cluster and adds the clusters to
        # the respective 0 or 1 side based on the assumption that each indiv.
        # will have a mother and a father.  Once a cluster has been added, the
        # members of that cluster are queued up so that their cluster
        # representing the other parent can be added to the opposite partition.
        aHS = []
        addedClusters = [False] * len(clusters)
        # Don't process the empty clusters
        for i,cluster in enumerate(clusters):
            if len(cluster) == 0:
                addedClusters[i] = True

        numClusters = len(clusters)
        for i in range(numClusters):
            # Only bother processing clusters which have not yet been assigned
            if not addedClusters[i]:
                q = Queue.Queue()
                tmp = [[],[]]

                tmp[0].append(copy.deepcopy(clusters[i]))
                addedClusters[i] = True

                # Stores a list of found individuals
                for ind in clusters[i]:
                    q.put([ind, i, 0])

                # Add the other sex parents to the opposite partition
                while not q.empty():
                    cur = q.get()

                    # Find the cluster which has not just been added
                    if cur[1] == locations[cur[0]][0]:
                        j = locations[cur[0]][1]
                    else:
                        j = locations[cur[0]][0]

                    # Add the cluster and mark as added
                    pMIdx = (cur[2] + 1) % 2
                    if not addedClusters[j]:
                        tmp[pMIdx].append(copy.deepcopy(clusters[j]))
                        addedClusters[j] = True
                        
                        for ind in clusters[j]:
                            q.put([ind, j, pMIdx])

                aHS.append(tmp)

        aPHS = []
        aMHS = []
        for s in aHS:
            vi0 = [self.computeVI(s[0], rPHS, numIndvs),\
                self.computeVI(s[0], rMHS, numIndvs)]
            vi1 = [self.computeVI(s[1], rPHS, numIndvs),\
                self.computeVI(s[1], rMHS, numIndvs)]
            # Add clusters to paternal and maternal HS sets based on which 
            # configuration minimizes the variation of information
            if vi0[0] + vi1[1] < vi0[1] + vi1[0]:
                aPHS.extend(s[0])
                aMHS.extend(s[1])
            else:
                aPHS.extend(s[1])
                aMHS.extend(s[0])

        vi = (self.computeVI(aPHS, rPHS, numIndvs) +\
            self.computeVI(aMHS, rMHS, numIndvs))/2.0
        viNorm = vi / math.log(numIndvs, 2)

        return viNorm

    def isBadJoin(self, cluster):
        pop = self.pop
        realHS = copy.deepcopy(self.realHS)

        # Check for unrepresented family
        for fam in realHS:
            # Check if the cluster is a subset of fam
            if set(cluster) <= set(fam[0]):
                return False

        return True

    def calculateLikelihood(self, cluster, pList, parents, verbose=False):
        nChildren = len(cluster)
        nLoci = self.pop.numLoci
        nAlleles = self.pop.numAlleles

        lpChildren = [0.0] * nChildren
        for l in range(nLoci):
            lpl = [0.0] * nChildren
            nValid = 0
            for i, validP in enumerate(pList[l]):
                if not validP:
                    continue

                nValid += 1
                counts = {}
                parent = list(parents[l][i])
                parent.sort()

                # Case 0: Parent is homozygotic
                if parent[0] == parent[1]:
                    for child in cluster:
                        ind = self.pop.individuals[child][2 * l : 2 * l + 2]
                        ind.sort()
                        # Identify which allele should contribute
                        if ind[0] == parent[0]:
                            counts[ind[1]] = counts.get(ind[1], 0) + 1
                        else:
                            counts[ind[0]] = counts.get(ind[0], 0) + 1
                # Case 1: Parent is heterozygotic
                else:
                    for child in cluster:
                        ind = self.pop.individuals[child][2 * l : 2 * l + 2]
                        ind.sort()
                        # Case A: Child looks identical to parent
                        # Add the expected weight (0.5) to each allele count
                        if parent == ind:
                            counts[ind[0]] = counts.get(ind[0], 0) + 0.5
                            counts[ind[1]] = counts.get(ind[1], 0) + 0.5
                        # Case B: Child is not identical
                        # Identify which allele should contribute
                        elif parent[0] == ind[0] or parent[1] == ind[0]:
                            counts[ind[1]] = counts.get(ind[1], 0) + 1
                        else:
                            counts[ind[0]] = counts.get(ind[0], 0) + 1
                
                if parent[0] == parent[1]:
                    for j, child in enumerate(cluster):
                        ind = self.pop.individuals[child][2 * l : 2 * l + 2]
                        ind.sort()
                        if ind[0] == parent[0]:
                            lpl[j] += -1.0 * math.log(\
                                float(counts[ind[1]] / float(nChildren)), 2)
                        else:
                            lpl[j] += -1.0 * math.log(\
                                float(counts[ind[0]] / float(nChildren)), 2)
                else:
                    for j, child in enumerate(cluster):
                        ind = self.pop.individuals[child][2 * l : 2 * l + 2]
                        ind.sort()
                        if parent == ind:
                            lpl[j] += -1.0 * math.log(\
                                (float(counts[ind[0]]) / float(nChildren) +\
                                float(counts[ind[1]]) / float(nChildren)) /\
                                2.0, 2)
                        elif parent[0] == ind[0] or parent[1] == ind[0]:
                            lpl[j] += -1.0 * math.log(\
                                float(counts[ind[1]] / float(nChildren)), 2)
                        else:
                            lpl[j] += -1.0 * math.log(\
                                float(counts[ind[0]] / float(nChildren)), 2)

            
            for j in range(nChildren):
                lpChildren[j] += lpl[j] / (float(nValid) * float(len(cluster)))
    
        avg = 0.0
        for i in lpChildren:
            avg += i
        avg /= float(len(lpChildren))

        for i, c in enumerate(cluster):
            lpChildren[i] = [math.fabs(lpChildren[i] - avg), c]
        lpChildren.sort()
        lpChildren.reverse()
        
        if verbose:
            print(len(cluster))
            for child in lpChildren:
                print(child)
            print('')

        return lpChildren

    def computeFPHelper(self, algCs, realCs):
        ret = []
        for aC in algCs:
            best = 10000
            aCS = set(tuple(aC))
            for rC in realCs:
                rCS = set(tuple(rC))
                tot = len(aCS) - len(aCS & rCS)
                if tot < best:
                    best = tot
                    tmpCS = rCS
            ret.append(list(aCS - (aCS & tmpCS)))

        return ret

    def computeFP(self, algCs):
        nIndvs = len(self.pop.individuals)

        rhs = [c[0] for c in self.realHS]
        res = self.computeFPHelper(algCs, rhs)

        indCounter = [0] * nIndvs
        for c in res:
            for ind in c:
                indCounter[ind] += 1

        ret = []
        for i in range(len(indCounter)):
            if indCounter[i] > 0:
                ret.append([i, indCounter[i]])

        return ret

    def printRealClusters(self, ind):
        for fam in self.realHS:
            if ind in fam[0]:
                print fam[0]
        
