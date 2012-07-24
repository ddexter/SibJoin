import evaluationTools

from SibJoin import SibJoin

sjCounts = [[0] * 40 for x in range(40)]
pop = []
for i in range(100):
    print(i)
    sj = SibJoin("pkl", fn="../tests/loci/40_6_05_5.pkl", rSeed=289 + i)

    res = sj.ca
    pop = sj.pop
    for tmp in res:
        for j in range(len(tmp)):
            for k in range(j + 1, len(tmp)):
                sjCounts[tmp[j]][tmp[k]] += 1
                sjCounts[tmp[k]][tmp[j]] += 1

    
eT = evaluationTools.EvaluationTools(pop=pop)
eT.printDistanceMatrix(sjCounts)
