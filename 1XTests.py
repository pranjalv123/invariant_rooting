import dendropy
import InvariantScores

treelist = dendropy.TreeList.get_from_path("true-trees-1X.tre", "newick")

S = dendropy.Tree.get_from_path("species-1X.tre", "newick", taxon_set = treelist.taxon_set)

print S.taxon_set == treelist.taxon_set
print len(treelist)

leaves = [n.taxon.label for n in S.leaf_nodes()] 

L = [leaves[i] for i in range(5)]

M = [leaves[m] for m in range(4,9)]

A = InvariantScores.get_dist(L,treelist)
print A[0]
print sum(A[0][1])

B = InvariantScores.get_dist(M,treelist)
print B[0]
print sum(B[0][1])
