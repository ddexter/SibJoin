import cplex

class MinRemovals:
    def __init__(self, individuals, mHS, pHS):
        self.ip = cplex.Cplex()

        # Build a reverse lookup for the parent clusters of an indiv.
        # parentLookup[ind] --> [mother cluster #, father cluster #]
        self.parentLookup = []
        for i in range(len(individuals)):
            tmp = []
            for j,cluster in enumerate(mHS):
                if i in cluster:
                    tmp.append(j)
                    break
            for j,cluster in enumerate(pHS):
                if i in cluster:
                    tmp.append(j)
                    break

            self.parentLookup.append(tmp)

        self.setup(individuals, mHS, pHS)
    
        self.solve(len(individuals))

    def setup(self, individuals, mHS, pHS):
        nIndvs = len(individuals)
        nLoci = len(individuals[0]) / 2

        # Minimization problem
        self.ip.objective.set_sense(self.ip.objective.sense.minimize)
        # Strong branching
        self.ip.parameters.mip.strategy.variableselect.set(3)
        # Set a 5 min. cutoff for opt.
        self.ip.parameters.tuning.timelimit.set(300)

        # x_i objective variables
        z = map(lambda x: 'x' + str(x), range(nIndvs))
        self.ip.variables.add(obj=[1]*nIndvs, names=z, types='B'*nIndvs,\
            lb=[0]*nIndvs, ub=[1]*nIndvs)

        nM = len(mHS)
        nP = len(pHS)
        relateVars = []
        constraints = []
        rhs = []
        sense = []

        for l in range(nLoci):
            # Find list of all alleles at locus
            alleles = set()
            for cluster in mHS:
                for idx in cluster:
                    ind = individuals[idx]

                    alleles.add(ind[2*l])
                    alleles.add(ind[2*l+1])

            # y_j,k,l variables
            #** sum(y_j,k,l) <= 2 constraints**
            for j in range(nM):
                # ym_j,k,l variables
                pAlleles = map(lambda a: "ym" + str(j) + ',' +\
                    str(a) + ',' + str(l), alleles)
                relateVars.extend(pAlleles)

                # sum(ym_j,k,l) <= 2 constraint
                constraints.append([pAlleles, [1] * len(pAlleles)])
                rhs.append(2)
                sense += 'L'

            # Do the same for the fathers
            for j in range(nP):
                # yp_j,k,l variables
                pAlleles = map(lambda a: "yp" + str(j) + ',' +\
                    str(a) + ',' + str(l), alleles)
                relateVars.extend(pAlleles)

                # sum(yp_j,k,l) <= 2 constraint
                constraints.append([pAlleles, [1] * len(pAlleles)])
                rhs.append(2)
                sense += "L"

            # x_i,l variables
            for i in range(nIndvs):
                ind = individuals[i]
                x = 'x' + str(i) + ',' + str(l)
                relateVars.append(x)
                parents = self.parentLookup[i]

                x0 = "x0" + ',' + str(i) + ',' + str(l)
                x1 = "x1" + ',' + str(i) + ',' + str(l)
                relateVars.append(x0)
                relateVars.append(x1)

                # For the right ym and yp (parents)
                # x0_i,l + 0.5(ym_a,l + yp_b,l) >= 1 constraints
                m0 = "ym" + str(parents[0]) + ',' + str(ind[2*l]) +\
                    ',' + str(l)
                p0 = "yp" + str(parents[1]) + ',' + str(ind[2*l+1]) +\
                    ',' + str(l)
                constraints.append([[m0, p0, x0], [0.5, 0.5, 1]])
                rhs.append(1)
                sense += "G"

                # Alternate case where 'b' allele is received from mother
                # x1_i,l + 0.5(ym_b,l + yp_a,l) >= 1 constraints
                m1 = "ym" + str(parents[0]) + ',' + str(ind[2*l+1]) +\
                    ',' + str(l)
                p1 = "yp" + str(parents[1]) + ',' + str(ind[2*l]) +\
                    ',' + str(l)
                constraints.append([[m1, p1, x1], [0.5, 0.5, 1]])
                rhs.append(1)
                sense += "G"
                
                # (x0_i,l + x1_i,l) - x_i,l <= 1 constraint
                constraints.append([[x0, x1, x], [1, 1, -1]])
                rhs.append(1)
                sense += "L"

                # x_i - x_i,l >= 0 constraint
                xi = 'x' + str(i)
                constraints.append([[xi, x], [1, -1]])
                rhs.append(0)
                sense += "G"

        # Add all variables and constraints to the IP
        self.ip.variables.add(names=relateVars, lb=[0] * len(relateVars),\
            ub=[1]*len(relateVars), types='B'*len(relateVars))
        self.ip.linear_constraints.add(lin_expr=constraints,\
            senses=sense, rhs=rhs)

    def solve(self, nIndvs):
        self.ip.solve()

        print("Solution Status: {0}".format(self.ip.solution.get_status()))
        print("Solution value: {0}".format(\
            self.ip.solution.get_objective_value()))
        
        names = self.ip.variables.get_names(0, nIndvs-1)
        variables = self.ip.solution.get_values(0, nIndvs-1)

        totalRemoved = 0
        for i in range(len(variables)):
            if variables[i] > 0.5:
                totalRemoved += 1
                print names[i], variables[i]

        print("Total Removed: {0}".format(totalRemoved))
        print("")

