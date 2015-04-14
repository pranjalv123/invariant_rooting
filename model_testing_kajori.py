import dendropy
import InvariantScores
import time


# STRAT Module testing for InvariantScores.taxon_with_split_edge
S = dendropy.Tree.get_from_path('/Users/kajori/Box Sync/UIUC/Tandy/invariant_rooting_v2/10-taxon/model.10.5400000.0.000000037/01/s_tree.trees','newick')
S.deroot()
S.encode_splits()
S.update_splits()
DS = S.split_edges
print 'split_bitmasks are'

print DS
ES = [e for e in S.postorder_edge_iter()]
INVDS = {v: k for k, v in DS.items()}
print INVDS
ESBITS = [INVDS[ES[i]] for i in range(len(ES))]
print 'dict that gives bitmasks in same order as edge list should be'
print ESBITS
S.print_plot()

for k in  ESBITS:
    print InvariantScores.taxon_with_split_edge(S,k)
# End Module testing for InvariantScores.taxon_with_split_edge