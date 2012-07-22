'''
SJGlobals is meant to be used as a "namespace" -- as close to a namespace as
python gets -- for variables which need to be shared across many SibJoin
sub-packages
'''
class SJGlobals:
    avgLinkage = False
    candidateParents = []
    isFS = []
    isHS = []
    joinHistory = []
    nIndvs = -1
    nLoci = -1
    strictAlleles = False

