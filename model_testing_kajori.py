import dendropy
import InvariantScores
import time
import copy 

#in progress
#formatted output
S = dendropy.Tree.get_from_path('/Users/kajori/Box Sync/UIUC/Tandy/invariant_rooting_v2/10-taxon/model.10.5400000.0.000000037/01/s_tree.trees','newick') #kajori
S.deroot()
S.encode_splits()
S.update_splits()
DS = S.split_edges

print 'split_bitmasks are'
print DS

print "species tree is"
print S.as_newick_string()

ES = [e for e in S.postorder_edge_iter()]
INVDS = {v: k for k, v in DS.items()}
ESBITS = [INVDS[ES[i]] for i in range(len(ES))]
print 'dict that gives bitmasks in same order as edge list should be'
print ESBITS


print 'there are' ,len(ES),'edges'

print 'edges in the species tree'
Splits = [S.taxon_set.split_taxa_list(ESBITS[i]) for i in range(len(ESBITS))]
print 'splits in order of edges are'
print Splits

#treelist = dendropy.TreeList.get_from_path('/Users/ruthdavidson/Downloads/10-taxon/model.10.5400000.0.000000037/01/truegenetrees','newick', taxon_set = S.taxon_set)
treelist = dendropy.Tree.get_from_path('/Users/kajori/Box Sync/UIUC/Tandy/invariant_rooting_v2/10-taxon/model.10.5400000.0.000000037/01/truegenetrees','newick', taxon_set = S.taxon_set) #kajori
print 'type(treelist)',type(treelist)
for i in range(len(treelist)):
    treelist[i].deroot()

score=[]
#I'm stopping at the second to last edge because I think in the postorder traversal the final edge is redundant due to Dendropy's "seed node" structure
for i in range(len(ES)-1):
    se = InvariantScores.total_quintet_score_kajori(S,i,treelist)
    s = 'edgescores after edge ' + str(i) + 'are'
    score.append(se)

### FORMATTED OUTPUT #########
ES_new=[e for e in S.postorder_edge_iter()]

print '\n FINAL OUTPUT'
print 'FORMAT edge.oid,score,split '
for index in range(len(ES_new)):
        taxon=InvariantScores.taxon_with_split_edge(S,INVDS_2[index])
        print index, ', score[',index,']=',score[index], ', split=',INVDS_2[index], taxon


print '\n score ='
print score

#T = InvariantScores.find_best_edge_by_total_quintet_score(S,treelist)

#print 'finding best rooting for species tree julies.genes.4.3.15.astral.tre from treelist julies.genes.4.3.15'

#:print T.as_newick_string()

'''
# START Unit testing for checking the following :-
#  1) Check if split_hash_bitmask_changes with deep copying -PASSED
#  2) Check if edge.oid with copying -FAILED
#  3)checking the taxons associated with the splits_hash_bitmask do not change by copying -PASSED
#  4) the splits_hash_bitmask associated with an edge donot change with copying
  
S_old= dendropy.Tree.get_from_path('/Users/kajori/Box Sync/UIUC/Tandy/invariant_rooting_v2/10-taxon/model.10.5400000.0.000000037/01/s_tree.trees','newick')
S_old.deroot()
S_old.encode_splits()
S_old.update_splits()
DS_old = S_old.split_edges
print 'OLD split_bitmasks are'
print DS_old.keys()

S_new=copy.deepcopy(S_old)
S_new.deroot()
S_new.encode_splits()
S_new.update_splits()
DS_new = S_new.split_edges
print 'NEW split_bitmasks are'
print DS_new.keys()



print '\n checking assertions'
assert len(set(DS_new.keys())-set(DS_old.keys()))==0
assert len(set(DS_old.keys())-set(DS_new.keys()))==0

#convert the keys to a list 
key_old=[ k for k in DS_old.keys()]
key_new=[ k for k in DS_new.keys()]
assert len(key_old)==len(key_new)
print '\n checking the taxons associated with the splits_hash_bitmask do not change by copying'
for index in range(len(key_old)):
    #print index
    taxon_set_old=InvariantScores.taxon_with_split_edge(S_old,key_old[index])
    #print taxon_set_old
    taxon_set_new=InvariantScores.taxon_with_split_edge(S_new,key_new[index])
    #print taxon_set_new
    assert (set(taxon_set_old) == set(taxon_set_new) and len(taxon_set_old) == len(taxon_set_new))
    assert key_old[index]==key_new[index] 
    #assert DS_old[key_old[index]].oid==DS_new[key_new[index]].oid #this assertion fails
    #print key_old[index],key_new[index],' = ',DS_old[key_old[index]].oid,DS_new[key_new[index]].oid
print 'assertions passed'


edgeset_old= [e for e in S_old.preorder_edge_iter()]
edgeset_new= [e for e in S_new.preorder_edge_iter()]
for index in range(len(edgeset_old)):
    #assert edgeset_old[index].oid==edgeset_new[index].oid #this assertion fails
    print index,edgeset_old[index],edgeset_new[index]
#END Unit testing
'''
'''
# START Module testing for InvariantScores.taxon_with_split_edge
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
'''