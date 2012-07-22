#!/usr/bin/python
import cplex

from SJGlobals import SJGlobals

'''
IP for identifying the minimum number of individuals that need to be
removed in order to assign a valid set of parents to the families
'''

class MinRemovals():
    def __init__(self, individuals, mHS, pHS):
        self.ip = cplex.Cplex()
        self.ip.set_results_stream(None)
        self.ip.set_log_stream(None)

        # Build a reverse lookup for the parent clusters of an indiv.
        self.parentLookup = []
        for ind in individuals:
            cnt = 0
            tmp = []
            for j,cluster in enumerate(mHS):
                indices = [s.index for s in cluster.individuals]
                if ind.index in indices:
                    tmp.append(j)
                    break
            for j,cluster in enumerate(pHS):
                indices = [s.index for s in cluster.individuals]
                if ind.index in indices:
                    tmp.append(j)
                    break

            self.parentLookup.append(tmp)

        self.setup(individuals, mHS, pHS)

    def setup(self, individuals, mHS, pHS):
        nIndvs = SJGlobals.nIndvs
        nLoci = SJGlobals.nLoci

        # Set as a maximization problem with low difference stopping
        # threshold
        self.ip.objective.set_sense(self.ip.objective.sense.maximize)
        # Strong branching
        self.ip.parameters.mip.strategy.variableselect.set(3)
        # Set a 5 minute cutoff for finding an optimal solution
        self.ip.parameters.timelimit.set(300)

        # x objective variables
        z = map(lambda x: "x" + str(x), range(nIndvs))
        self.ip.variables.add(obj = [1] * nIndvs, names = z,\
            types='B' * nIndvs, lb=[0] * nIndvs, ub=[1] * nIndvs)

        nM = len(mHS)
        nP = len(pHS)
        relateVars = []
        constraints = []
        rhs = []
        sense = []
        for l in range(nLoci):
            # find list of all alleles at a locus
            alleles = set()
            for cluster in mHS:
                for ind in cluster.individuals:
                    alleles.add(ind.loci[l][0])
                    alleles.add(ind.loci[l][1])

            # f_j,k,l variables
            # sum(f_j,k,l) for all k = 2
            # note: we can generalize a homozygotic parent to anything else
            for j in range(nM):
                pAlleles = map(lambda a: "fm" + str(j) + "," +\
                    str(a) + "," + str(l), alleles)

                # fm_j,k,l variables
                relateVars.extend(pAlleles)

                # sum(fm_j,k,l) for all k = 2
                constraints.append([pAlleles, [1] * len(pAlleles)])
                rhs.append(2)
                sense += "L"
            # Similarly for the fathers
            for j in range(nP):
                pAlleles = map(lambda a: "fp" + str(j) + "," +\
                    str(a) + "," + str(l), alleles)

                relateVars.extend(pAlleles)

                constraints.append([pAlleles, [1] * len(pAlleles)])
                rhs.append(2)
                sense += "L"
            
            # x_i,l variables for heterozygotes
            for i in range(nIndvs):
                ind = individuals[i]
                x = "x" + str(i) + "," + str(l)
                relateVars.append(x)
                parents = self.parentLookup[i]

                # x_i,0,l and x_i,1,l variables for heterozygotes
                if ind.loci[l][0] != ind.loci[l][1]:
                    x0 = "x0" + "," + str(i) + "," + str(l)
                    x1 = "x1" + "," + str(i) + "," + str(l)
                    relateVars.append(x0)
                    relateVars.append(x1)

                    # Constraint x_i,0,l <= 0.5(f_m,a,l + f_p,b,l)
                    m0 = "fm" + str(parents[0]) + "," +\
                        str(ind.loci[l][0]) + "," + str(l)
                    p0 = "fp" + str(parents[1]) + "," +\
                        str(ind.loci[l][1]) + "," + str(l)
                    constraints.append([[m0, p0, x0], [0.5, 0.5, -1]])
                    rhs.append(0)
                    sense += "G"

                    # Constraint x_i,1,l <= 0.5(f_m,b,l + f_p,a,l)
                    m1 = "fm" + str(parents[0]) + "," +\
                        str(ind.loci[l][1]) + "," + str(l)
                    p1 = "fp" + str(parents[1]) + "," +\
                        str(ind.loci[l][0]) + "," + str(l)
                    constraints.append([[m1, p1, x1], [0.5, 0.5, -1]])
                    rhs.append(0)
                    sense += "G"

                    # Constraint x_i,l <= x_i,0,l + x_i,1,l
                    constraints.append([[x0, x1, x], [1, 1, -1]])
                    rhs.append(0)
                    sense += "G"
                else:
                    # Constraint x_i,1,l <= 0.5(f_m,a,l + f_p,a,l)
                    m0 = "fm" + str(parents[0]) + "," +\
                        str(ind.loci[l][0]) + "," + str(l)
                    p0 = "fp" + str(parents[1]) + "," +\
                        str(ind.loci[l][0]) + "," + str(l)
                    constraints.append([[m0, p0, x], [0.5, 0.5, -1]])
                    rhs.append(0)
                    sense += "G"

        # Constraint x_i <= (1/|L|) * sum(x_i,l) for all l
        for i in range(nIndvs):
            x = "x" + str(i)

            lhs = map(lambda l: "x" + str(i) + "," + str(l), range(nLoci))
            lhs.append(x)

            right = [1.0/nLoci] * nLoci
            right.append(-1)

            constraints.append([lhs, right])
            rhs.append(0)
            sense += "G"

        self.ip.variables.add(names = relateVars, lb=[0] * len(relateVars),\
            ub = [1] * len(relateVars), types='B' * len(relateVars))
        self.ip.linear_constraints.add(lin_expr = constraints,\
            senses=sense, rhs=rhs)

    def solve(self):
        nIndvs = SJGlobals.nIndvs
        individuals = SJGlobals.individuals

        self.ip.solve()

        '''
        print("Solution status: {0}".format(self.ip.solution.get_status()))
        print("Solution value: {0}".format(\
            self.ip.solution.get_objective_value()))
        '''

        names = self.ip.variables.get_names(0, nIndvs - 1)
        variables = self.ip.solution.get_values(0, nIndvs - 1)
        ret = []
        for i in range(len(variables)):
            if variables[i] < 0.5:
                ret.append(individuals[i])
                #print names[i], variables[i]

        return ret
