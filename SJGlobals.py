'''
SJGlobals is meant to be used as a "namespace" -- as close to a namespace as
python gets -- for variables which need to be shared across many SibJoin
sub-packages
'''
class SJGlobals:
    allowableClusterJoins = []
    allowableJoins = []
    candidateParents = []
    hasCandidateParents = False
    joinHistory = []
    nIndvs = -1
    nLoci = -1
    strictAlleles = False

    @classmethod
    def clear(cls):
        cls.allowableJoins = []
        cls.allowableClusterJoins = []
        cls.candidateParents = []
        cls.hasCandidateParents = False
        cls.joinHistory = []
        cls.nIndvs = -1
        cls.nLoci = -1
        cls.strictAlleles = False
