'''
It runs tests on the 10-taxon dataset  on 4 mode conditions, 3 replicates
using MORE LOCAL quintets and reports the scores
'''

import dendropy
import InvariantScores
import time
import copy 


# the model_coditions are arranged from lowest to highest ILS
model_coditions=['model.10.5400000.0.000000037','model.10.1800000.0.000000111','model.10.600000.0.000000333','model.10.200000.0.000001000']

for mc in model_coditions:
    for replicate in range(1,4):
        S = dendropy.Tree.get_from_path('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/10-taxon/'+mc+'/0'+str(replicate)+'/s_tree.trees','newick') #kajori
        print '\n ********** '+mc+'/0'+str(replicate)+ '/ ********** '
        S.deroot()
        S.encode_splits()
        S.update_splits()
        DS = S.split_edges

        ES = [e for e in S.postorder_edge_iter()]
        INVDS = {v: k for k, v in DS.items()}
        #dict that gives bitmasks in same order as edge list should be
        ESBITS = [INVDS[ES[i]] for i in range(len(ES))]
    
        treelist =  dendropy.TreeList.get_from_path('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/10-taxon/'+mc+'/0'+str(replicate)+'/truegenetrees','newick', taxon_set = S.taxon_set) 
        for i in range(len(treelist)):
            treelist[i].deroot()


        print(S.as_ascii_plot())
        #I'm stopping at the second to last edge because I think in the postorder traversal the final edge is redundant due to Dendropy's "seed node" structure
        ### FORMATTED OUTPUT ########
        for i in range(len(ES)-1):
            (se,U,quintet)= InvariantScores.total_quintet_score_distance_kajori(S,i,treelist)
            taxon=InvariantScores.taxon_with_split_edge(S,ESBITS[i])
            # here i = index of edge in postorder_iteration
            print 'index =',i,' score = ',se, 'quintet = ',quintet,' taxon=',taxon, ' U =',U
        


