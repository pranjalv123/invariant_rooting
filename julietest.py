import dendropy
import InvariantScores

S = dendropy.Tree.get_from_path('/Users/ruthdavidson/code/Outroup_Rooting_Project/julies.genes.4.3.15.astral.tre','newick')

Taxa = [n.taxon.label for n in S.leaf_nodes()]

print 'taxon set is'

print Taxa

treelist = dendropy.TreeList.get_from_path('/Users/ruthdavidson/code/Outroup_Rooting_Project/julies.genes.4.3.15','newick', taxon_set = S.taxon_set)

l = [Taxa[i] for i in range(5)]

d = InvariantScores.get_dist(l,treelist)

print 'quintet dist for first 5 taxa is'

print d

#T = InvariantScores.find_best_edge_by_total_quintet_score(S,treelist)

#print 'finding best rooting for species tree julies.genes.4.3.15.astral.tre from treelist julies.genes.4.3.15'

#:print T.as_newick_string()
