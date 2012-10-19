#/usr/bin/python
import cplex
import random

from SJGlobals import SJGlobals
from containers.Cluster import Cluster
from containers.Clusters import Clusters
from containers.Individual import Individual

# z = cluster objective variable
# a = allele variable
# p = parent allele storage variable
# y = inherited allele variable
# x = individual to cluster assignment variable
# i = individual
# k = allele
# j = family


class OriginalIP():
    def __init__(self, individuals, guesses=[0, 0]):
        self.ip = cplex.Cplex()
        self.ip.set_results_stream(None)
        self.ip.set_log_stream(None)

        nClusters = guesses
        for i in range(len(nClusters)):
            if nClusters[i] == 0:
                nClusters[i] = len(individuals)

        nIndvs = SJGlobals.nIndvs
        nLoci = SJGlobals.nLoci
        self.individuals = individuals
        self.nClusters = nClusters
        self.nIndvs = nIndvs
        self.nLoci = nLoci

        constraints= []
        rhs = []
        rVars = []
        sense = ""
    
        self.ip.objective.set_sense(self.ip.objective.sense.minimize)
        self.ip.parameters.workdir.set("/Users/ddexter/cplex/")
        # Tree limit 3 GB
        self.ip.parameters.mip.limits.treememory.set(3000)
        # Strong branching
        self.ip.parameters.mip.strategy.variableselect.set(3)
        # 1 hour cutoff for finding an optimal solution
        #self.ip.parameters.timelimit.set(3600)

        # z objective variables
        z = map(lambda v: "z0," + str(v), range(nClusters[0]))
        z.extend(map(lambda v: "z1," + str(v), range(nClusters[1])))

        # x_i,j variables
        x = []
        for s in range(2):
            for ind in individuals:
                # x_i,j variables
                x += map(lambda v: "x" + str(s) + "," + str(ind.index) +\
                    "," + str(v), range(nClusters[s]))
        rVars.extend(x)

        zObj = [1] * sum(nClusters)
        xObj = []
        numObj = sum(nClusters) + len(x)
        # Random perturbations to force a family
        for i in range(len(x)):
            xObj.append(random.uniform(0.0000000001, 0.000000001))

        self.ip.variables.add(obj = zObj + xObj, names=z+x,\
            types='B' * numObj, lb=[0] * numObj, ub=[1] * numObj)

        # Find list of alleles at each locus
        alleles = [set()] * nLoci
        for ind in individuals:
            for l in range(nLoci):
                alleles[l].add(ind.loci[l][0])
                alleles[l].add(ind.loci[l][1])

        # sum z_{j+1} - z_j <= 0
        # make the families contiguous to cut down on the k! possible families
        for s in range(2):
            for j in range(nClusters[s] - 1):
                zJ = "z" + str(s) + "," + str(j)
                zJP = "z" + str(s) + "," + str(j+1)
                constraints.append([[zJP, zJ], [1,-1]])
                rhs.append(0)
                sense += "L"

        for s in range(2):
            for ind in individuals:
                xI = map(lambda v: "x" + str(s) + "," + str(ind.index) +\
                    "," + str(v), range(nClusters[s]))
                # sum(x_i,j) for all j = 1
                constraints.append([xI, [1] * len(xI)])
                rhs.append(1)
                sense += "E"

        for s in range(2):
            for ind in individuals:
                for j in range(nClusters[s]):
                    # z_j - x_i,j >= 0
                    xS = "x" + str(s) + "," + str(ind.index) +\
                        "," + str(j)
                    zS = "z" + str(s) + "," + str(j)
                    constraints.append([[zS, xS], [1, -1]])
                    rhs.append(0)
                    sense += "G"

                for l in range(nLoci):
                    # y_i,k,l variables (keep track of 2 alleles for ind)
                    y0 = "y" + str(s) + "," + str(ind.index) + "," +\
                        str(ind.loci[l][0]) + "," + str(l)
                    y1 = "y" + str(s) + "," + str(ind.index) + "," +\
                        str(ind.loci[l][1]) + "," + str(l)
                    # heterozygotic
                    if ind.loci[l][0] != ind.loci[l][1]:
                        y = [y0, y1]
                        rVars.extend(y)

                        # sum(y_i,k,l) for k_0 and k_1 = 1
                        constraints.append([y, [1, 1]])
                        rhs.append(1)
                        sense += "E"

                        if s == 1:
                            y00 = "y0" + "," + str(ind.index) + "," +\
                                str(ind.loci[l][0]) + "," + str(l)
                            y01 = "y0" + "," + str(ind.index) + "," +\
                                str(ind.loci[l][1]) + "," + str(l)

                            # sum(y0_i,k,l + y1_i,k,l) = 1
                            # Each partition must cover a unique allele
                            constraints.append([[y0, y00], [1,1]])
                            rhs.append(1)
                            sense += "E"
                            constraints.append([[y1, y01], [1,1]])
                            rhs.append(1)
                            sense += "E"

                    # homozygotic -- don't double add variables
                    else:
                        rVars.append(y0)
                        constraints.append([[y0], [1]])
                        rhs.append(1)
                        sense += "E"

            for j in range(nClusters[s]):
                for l in range(nLoci):
                    # p_j,k,l variables
                    p = map(lambda v: "p" + str(s) + "," + str(j) + "," +\
                        str(v) + "," + str(l), alleles[l])
                    rVars.extend(p)

                    # sum(p_j,k,l) for all k <= 2
                    constraints.append([p, [1] * len(alleles[l])])
                    rhs.append(2)
                    sense += "L"

                    for k in alleles[l]:
                        for ind in individuals:
                            if ind.loci[l][0] == k or ind.loci[l][1] == k:
                                # p_j,k,l - x_i,j - y_i,k,l >= -1
                                pS = "p" + str(s) + "," + str(j) + "," +\
                                    str(k) + "," + str(l)
                                xS = "x" + str(s) + "," + str(ind.index) +\
                                    "," + str(j)
                                yS = "y" + str(s) + "," + str(ind.index) +\
                                    "," + str(k) + "," + str(l)
                                constraints.append(\
                                    [[pS, xS, yS], [1, -1, -1]])
                                rhs.append(-1)
                                sense += "G"
        
        self.ip.variables.add(names=rVars, lb=[0] * len(rVars),\
            ub=[1] * len(rVars), types='B' * len(rVars))
        self.ip.linear_constraints.add(lin_expr=constraints,\
            senses=sense, rhs=rhs)

        self.ip.solve()

    def getClusters(self):
        labels = []
        variables = []
        nIndvs = self.nIndvs
        mHS = [[] for x in range(self.nClusters[0])]
        pHS = [[] for x in range(self.nClusters[1])]

        for s in range(2):
            nClusters = self.nClusters[s]
            varStart = sum(self.nClusters) + s * self.nClusters[0] * nIndvs
            varEnd = varStart + nIndvs * nClusters - 1

            variables.extend(self.ip.solution.get_values(varStart, varEnd))
            labels.extend(self.ip.variables.get_names(varStart, varEnd))

            for ind in self.individuals:
                for j in range(self.nClusters[s]):
                    if variables[ind.index * nClusters + j +\
                        s * self.nClusters[0] * nIndvs] == 1:
                        if s == 0:
                            mHS[j].append(ind)
                        else:
                            pHS[j].append(ind)

        hsClusters = []
        clustNum = 0
        for cluster in mHS:
            if len(cluster) > 0:
                hsClusters.append(Cluster(clustNum, cluster))
                clustNum += 1
        for cluster in pHS:
            if len(cluster) > 0:
                hsClusters.append(Cluster(clustNum, cluster))
                clustNum += 1
        clusters = Clusters([], hsClusters)
        '''
        for cluster in hsClusters:
            print [ind.index for ind in cluster.individuals]
        '''

        return clusters
        
