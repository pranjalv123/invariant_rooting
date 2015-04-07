import dendropy
import InvariantScores
import time

S = dendropy.Tree.get_from_path('/Users/ruthdavidson/code/Outroup_Rooting_Project/julies.genes.4.3.15.astral.tre','newick')

Taxa = [n.taxon.label for n in S.leaf_nodes()]

edgelist = [e for e in S.postorder_edge_iter()]

print 'taxon set is'

print Taxa

treelist = dendropy.TreeList.get_from_path('/Users/ruthdavidson/code/Outroup_Rooting_Project/julies.genes.4.3.15','newick', taxon_set = S.taxon_set)

l = [Taxa[i] for i in range(5)]

d = InvariantScores.get_dist(l,treelist)

print 'quintet dist for first 5 taxa is'

print d

print 'calculating score for 4th edge'

t0 = time.clock()

score4 = InvariantScores.total_quintet_score(S,4,treelist)

t1 = time.clock()

print 'the score for the 4th edge is' 

print score4

print 'it took processor time'

print t1-t0

print 'to score the 4th edge'

edgescores = []

print 'there are this many edges in the species tree:'

print len(edgelist)

#print 'edge 0'
for i in range(len(edgelist)):
    edgescores.append(InvariantScores.total_quintet_score(S,i,treelist))

print 'edgescores are'

print edgescores

#T = InvariantScores.find_best_edge_by_total_quintet_score(S,treelist)

#print 'finding best rooting for species tree julies.genes.4.3.15.astral.tre from treelist julies.genes.4.3.15'

#:print T.as_newick_string()
