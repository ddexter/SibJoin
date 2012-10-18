'''
SJGlobals is meant to be used as a "namespace" -- as close to a namespace as
python gets -- for variables which need to be shared across many SibJoin
sub-packages
'''
class SJGlobals:
    hasCandidateParents = False
    allowableClusterJoins = []
    allowableJoins = []
    candidateParents = []
    nIndvs = -1
    nLoci = -1
    strictAlleles = False

    def clear(self):
        self.allowableJoins = []
        self.allowableClusterJoins = []
        self.candidateParents = []
        self.hasCandidateParents = False
        self.nIndvs = -1
        self.nLoci = -1
        self.strictAlleles = False
