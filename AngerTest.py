import dendropy
import InvariantScores
treelist = dendropy.TreeList.get_from_path("nardtest_song_mammals.424.gene.tre", "newick")
print len(treelist)

S = dendropy.Tree.get_from_path("nardtest_song_mammals.species.tre", "newick", taxon_set = treelist.taxon_set)
#print S.taxon_set == treelist.taxon_set

leaves = [n.taxon.label for n in S.leaf_nodes()]
#print leaves

print len(leaves)
#LL = [[leaves[i] for i in range(0 + j, 5 + j)] for j in range(len(leaves))]

L = [leaves[i] for i in range(5)]
M = [leaves[i] for i in range(3,8)]
K = [leaves[i] for i in range(10,15)]
#A = InvariantScores.get_dist(L,treelist)
#print A[0]
#print errors
#K = [leaves[i] for i in range(3,6)]
#
B = InvariantScores.get_dist(M,treelist)
print B[0]
print sum(B[0][1])
#print [B[1][i] for i in range(10)]

C = InvariantScores.get_dist(L,treelist)
print C[0]
print sum(C[0][1])

D = InvariantScores.get_dist(K,treelist)
print D[0]
print sum(D[0][1])
