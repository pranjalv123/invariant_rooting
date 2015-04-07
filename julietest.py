import dendropy
import InvariantScores
import time

S = dendropy.Tree.get_from_path('/Users/ruthdavidson/code/Outroup_Rooting_Project/julies.genes.4.3.15.astral.tre','newick')

S.deroot()

S.encode_splits()
S.update_splits()

DS = S.split_edges

print 'split_bitmasks are'

print DS

print "species tree is"

print S.as_newick_string()

Taxa = [n.taxon.label for n in S.leaf_nodes()]

edgelist = [e for e in S.postorder_edge_iter()]

INVDS = {v: k for k, v in DS.items()}

print INVDS

ESBITS = [INVDS[ES[i]] for i in range(len(ES))]

print 'dict that gives bitmasks in same order as edge list should be'

print ESBITS


print 'there are'

print len(edgelist)

print 'edges in the species tree'

H = open('log.julietest','w')

print 'taxon set is'

print Taxa

Splits = [S.taxon_set.split_taxa_list(ESBITS[i]) for i in range(len(ESBITS))]

print 'splits in order of edges are'

print Splits

treelist = dendropy.TreeList.get_from_path('/Users/ruthdavidson/code/Outroup_Rooting_Project/julies.genes.4.3.15','newick', taxon_set = S.taxon_set)

for i in range(len(treelist)):
    treelist[i].deroot()

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

edgescores = [0.0 for i in range(len(edgelist)-1)]

print 'there are this many edges in the species tree, because one is a dendropy artifact:'

print len(edgelist)

#print 'edge 0'
#for i in range(len(edgelist)):

#I'm stopping at the second to last edge because I think in the postorder traversal the final edge is redundant due to Dendropy's "seed node" structure
for i in range(len(edgelist)-1):
    se = InvariantScores.total_quintet_score(S,i,treelist)
    #edgescores.append(se)
    edgescores[i] = se
    s = 'edgescores after edge ' + str(i) + 'are'
    H.write(s)
    H.write(str(edgescores))
    print 'edgescores are'
    print edgescores


print 'edgescores are'

print edgescores

print 'edgescores are indexed by the list'

print edgelist

print 'minimum scoring edge is has index'

minscore = edgescores.index(min(edgescores))

H.write('minimum scoring edge is has index')
H.write(str(minscore))

print minscore

T = dendropy.Tree(S)

T.reroot_at_edge(edgelist[minscore])

print 'tree rooted at minscoring edge is'

print T.as_ascii_plot()
H.write(str(T.as_ascii_plot()))

H.close()
#T = InvariantScores.find_best_edge_by_total_quintet_score(S,treelist)

#print 'finding best rooting for species tree julies.genes.4.3.15.astral.tre from treelist julies.genes.4.3.15'

#:print T.as_newick_string()
