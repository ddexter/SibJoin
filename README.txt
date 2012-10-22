SibJoin is a heuristic for rapidly reconstructing half-sibling families from microsatellite DNA information of diploid (and some haplo-diploid) individuals.  On simulated and real population sets, we achieve accuracy similar to COLONY 2.0, but thousands to tens of thousands of times faster.

SibJoin is free to use under the MIT License (see below).  However, please cite our work when using the SibJoin algorithm:
    http://link.springer.com/chapter/10.1007/978-3-642-33122-0_4

Copyright (C) <2012> <Daniel G. Brown, Daniel C. Dexter>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

SibJoin architecture:
The current deployment is written in Python and requires at least Python 2.7

3rd Party Packages:
    matplotlib
    pickle
    networkx

Entry Point:
    SibJoin.py: Takes a file and calls SibJoinBuilder object to prepare
	the algorithm.  The run function handles searching at different
	thresholds.
    
    Example (Run with CSV file):
	python SibJoin.py individuals.txt

    CSV Input Format (denote missing alleles with -1):
	Animal0,13,17,-1,-1,-1,-1
	Animal1,13,19,10,3,7,7
	.
	.
	.
	AnimalN,18,17,3,4,8,2

Container Classes:
    Containers hold information about individuals and clusters

    containers/Cluster.py: An individual cluster object which holds a group
	of individuals.
    containers/Clusters.py: Holds many Cluster objects. Clusters is also
	where the join algorithm is implemented
    containers/Individual.py: Stores relevant information about an
	individual, including alleles and loci, half-sibling cluster
	location, and full sibling cluster location.
    containers/Bipartite.py: Structure for the internal bipartite graph
	used to enforce that no individual ends up with two parents of the
	same sex.

Evaluation:
    SibJoin uses variation of information to assess accuracy, but may only
    be used when the actual population is known.  Therefore, only for
    benchmarking the algorithm.

    EvaluationToolkit.py:  To run this, after SibJoin completes, an
	evaluation instance must be established by passing the population
	structure of the true population e.g.:
	    eval = EvaluationToolkit(sjInstance.builder.pop)

Contact:
    ddexter@uwaterloo.ca with questions

