from SJGlobals import SJGlobals
from Population import Population
from MinRemovals import MinRemovals

"""""""""""""""""""""""""""""""
Runs random SibJoin r times and calculates the pairwise co-clustered count
for each pair of individuals.
"""""""""""""""""""""""""""""""
class SJEC:
    def __init__(self, filetype, r, fn=None, pop=None,\
        scoreSensitivity = 4.0, rSeed = 289):

        # Run SibJoin to get the deterministic solution
        sj = SibJoin(self, filetype, fn=fn, pop=pop,\
            scoreSensitivity=scoreSensitivity, rSeed=rSeed):

        '''
        Run Random SibJoin and count how often pairs of indivs. share a
        family.  We use this to determine replacement of incompatible
        indivs. in the original SibJoin instance.
        '''
        ccMtx = []
        for i in range(r):
            sjr = SibJoinRandom(self, filetype, fn=fn, pop=pop,\
                scoreSensitivity=scoreSensitivity, rSeed=rSeed)

            '''
            Now that we have a count of the number of indivs., we can
            construct the co-clustered matrix if it hasn't already been
            constructed.
            '''
            if ccMtx == []:
                ccMtxt = [[0] * SJGlobals.nIndvs\
                    for j in range(SJGlobals.nIndvs)]
            


        
        

