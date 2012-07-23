import networkx as nx

from networkx.algorithms import bipartite
from SJGlobals import SJGlobals

class Bipartite:
    def __init__(self, clusters):
        self.g = nx.Graph()
        individuals = SJGlobals.individuals

        for cluster in clusters:
            self.g.add_node(cluster.clusterID)

        for cluster in clusters:
            for ind in cluster.individuals:
                self.g.add_edge(ind.hsClusters[0], ind.hsClusters[1])

    def combineNodes(self, clusterID0, clusterID1):
        # Get edge list of n1
        neighbors = self.g.neighbors(clusterID1)
        # Set the neighbors of n1 as adjacent to n0
        newEdges = []
        for neighbor in neighbors:
            if not self.g.has_edge(clusterID0, neighbor):
                newEdges.append([clusterID0, neighbor])
                self.g.add_edge(clusterID0, neighbor)

        self.g.remove_node(clusterID1)

        return [clusterID1, newEdges, neighbors]

    def undoCombineNodes(self, clusterID, newEdges, neighbors):
        self.g.add_node(clusterID)

        for e in newEdges:
            self.g.remove_edge(e[0], e[1])

        for n in neighbors:
            self.g.add_edge(clusterID, n)

    # joins expects a list of pairs of elements to join.  Testing for valid
    # full-sib joins requires testing more than one join at a time, so
    # always make joins a list
    def isValidCombining(self, joins):
        rollback = []
        isBipartite = False

        # Test the joins
        for i,join in enumerate(joins):
            join.sort()
            rollback.append(self.combineNodes(join[0], join[1]))

        isBipartite = nx.bipartite.is_bipartite(self.g)
        
        # Roll back the joins
        for join in reversed(rollback):
            self.undoCombineNodes(join[0], join[1], join[2])

        return isBipartite

    def relabelNodes(self, nodeNameMapping):
        nx.relabel_nodes(self.g, nodeNameMapping, copy=False)

    def getPartitionings(self):
        partitionings = [list(bipartite.sets(subG))\
            for subG in nx.connected_component_subgraphs(self.g)]

        return partitionings

