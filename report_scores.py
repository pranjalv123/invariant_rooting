import dendropy
import InvariantScores
import time
import copy 

model_coditions=['model.10.5400000.0.000000037','model.10.1800000.0.000000111','model.10.600000.0.000000333','model.10.200000.0.000001000']

#It runs tests on the Avian dataset 
#select a quintet on the basis of the topological distance from the root
#and score all the edges on that quintet

S = dendropy.Tree.get_from_path('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/Avian_0.5X/avian_species_tree.trees','newick') #kajori
for replicate in range(1,4): #replicate
    #combine the gene trees from the 1000 sub-dir
    f_write = open('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/Avian_0.5X/0.5X-1000-true/combined_truegenetrees/truegenetrees_R'+str(replicate), 'w')
    for gt in range((replicate-1)*1000+1,replicate*1000):
        f_read = open('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/Avian_0.5X/0.5X-1000-true/R'+str(replicate)+'/'+str(gt)+'/true.gt','r')
        f_write.write(f_read.read())
        f_read.close()
    f_write.close()
    treelist =  dendropy.TreeList.get_from_path('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/Avian_0.5X/0.5X-1000-true/combined_truegenetrees/truegenetrees_R'+str(replicate),'newick') #kajori

    print '\n ********** /0'+str(replicate)+ '/ ********** '
    print(S.as_ascii_plot())
    
    quintet = ['GEOFO','APTFO','HALAL','CHAPE','PODCR']
    print 'quintet',quintet
    
    score=InvariantScores.edge_score_on_quintet(S,quintet,treelist)
    print 'score = ',score
    
'''
#It runs tests on the 10-taxon dataset  on 4 mode conditions, 3 replicates
#it selects a quintet and scores the edges on the quintet and checks if the method rules out leaf edges
for mc in model_coditions:
    for replicate in range(1,4):
        S = dendropy.Tree.get_from_path('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/10-taxon/'+mc+'/0'+str(replicate)+'/s_tree.trees','newick') #kajori
        treelist =  dendropy.TreeList.get_from_path('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/10-taxon/'+mc+'/0'+str(replicate)+'/truegenetrees','newick', taxon_set = S.taxon_set) 
        
        print '\n ********** '+mc+'/0'+str(replicate)+ '/ ********** '
        print(S.as_ascii_plot())
        print 'quintet',['1','3','5','7','8']
    
        score=InvariantScores.edge_score_on_quintet(S,['1','3','5','7','8'],treelist)
        print 'score = ',score

'''
'''
#It runs tests on the 10-taxon dataset  on 4 mode conditions, 3 replicates
#using MORE LOCAL quintets and reports the scores
# the model_coditions are arranged from lowest to highest ILS

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
            print 'index =',i,' score = ',se, 'quintet = ',quintet,' taxon=',taxon, ' U =',U,'\n'

'''