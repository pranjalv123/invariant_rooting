import dendropy
import InvariantScores
import time
import copy 

model_coditions=['model.10.5400000.0.000000037' ,'model.10.1800000.0.000000111','model.10.600000.0.000000333','model.10.200000.0.000001000']
# the model_coditions are arranged from lowest to highest ILS


#It runs tests on the Avian dataset 
#select a quintet on the basis of the topological distance from the root
#and score all the edges on that quintet
'''

#PART 1
S = dendropy.Tree.get_from_path('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/Avian_0.5X/avian_species_tree.trees','newick') #kajori
S.deroot()
for replicate in range(4,5): #replicate
    print '\n ********** '+str(replicate)+ '/ ********** '
    
    #combine the gene trees from the 1000 sub-dir
    f_write = open('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/Avian_0.5X/0.5X-1000-true/combined_truegenetrees/truegenetrees_R'+str(replicate), 'w')
    for gt in range((replicate-1)*1000+1,replicate*1000):
        f_read = open('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/Avian_0.5X/0.5X-1000-true/R'+str(replicate)+'/'+str(gt)+'/true.gt','r')
        f_write.write(f_read.read())
        f_read.close()
    f_write.close()
   
    treelist =  dendropy.TreeList.get_from_path('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/Avian_0.5X/0.5X-1000-true/combined_truegenetrees/truegenetrees_R'+str(replicate),'newick') #kajori
    
    quintet = ['GEOFO','APTFO','HALAL','CHAPE','PODCR']
    score,U=InvariantScores.edge_score_on_quintet(copy.deepcopy(S),quintet,treelist)
    print 'score = ',score
    print 'U =',U
    
    output_filename='/Users/kajori/Box Sync/UIUC/Tandy/invariant_rooting/output_trees/'+str(replicate)+'_avian_species_tree_wih_score.trees'
    InvariantScores.my_print_tree(copy.deepcopy(S),score,output_filename)
    
    f_write=open('/Users/kajori/Box Sync/UIUC/Tandy/invariant_rooting/output_trees/U_distribution.txt','a')
    f_write.write('\n\n Avian_0.5X/0.5X-1000-true/R'+str(replicate)+'/ \n U = ')
    f_write.write(str(U))
    print 'str(U)',str(U), type(str(U))
    f_write.close()
   
    print 'quintet',quintet
    
'''


#PART 3
#It runs tests on the 10-taxon dataset  on 4 mode conditions, 3 replicates
#it selects a quintet and scores the edges on the quintet and checks if the method rules out leaf edges
#f_write=open('/Users/kajori/Box Sync/UIUC/Tandy/invariant_rooting/output_trees/U_distribution_10.txt','a')
for mc in model_coditions:
    for replicate in range(1,2):
        S = dendropy.Tree.get_from_path('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/10-taxon/'+mc+'/0'+str(replicate)+'/s_tree.trees','newick') #kajori
        S.deroot()
        treelist =  dendropy.TreeList.get_from_path('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/10-taxon/'+mc+'/0'+str(replicate)+'/truegenetrees','newick', taxon_set = S.taxon_set) 
        
        print '\n ********** Model Condition -'+mc+'/ Replicate -'+str(replicate)+ '/ ********** '
        #print(S.as_ascii_plot())
        #quintet=['1','3','5','7','8']
        quintet=['0','3','5','7','8']
        print 'quintet',quintet
    
        
        score,U=InvariantScores.edge_score_on_quintet(copy.deepcopy(S),quintet,treelist)
        print 'original tree for debugging purpose'
        InvariantScores.my_print_tree(S,score,'dump.trees')
        print ' Corresponding 5-taxon tree'
        S = dendropy.Tree.get_from_path('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/10-taxon/'+mc+'/0'+str(replicate)+'/s_tree.trees','newick') #kajori
        S.deroot()
        InvariantScores.format_tree(S,quintet,score,'dump.trees')
        print 'U = ',U
        
        
        print '\n Analysis:- \n 1)best score on the dataset - ',min(score),'\n 2) # edges that have the best score - ',score.count(min(score))
        #S_new=dendropy.Tree.clone_from(S)
    

        #Clones the structure and properties of Tree object other
        #output_filename='/Users/kajori/Box Sync/UIUC/Tandy/invariant_rooting/output_trees/'+str(mc)+'_'+str(replicate)+'_five_taxon_tree_10_taxon_dataset.trees'
        
        #output_filename='/Users/kajori/Box Sync/UIUC/Tandy/invariant_rooting/output_trees/'+str(mc)+'_'+str(replicate)+'_10_taxon_tree_with_score.trees'
        #InvariantScores.my_print_tree(copy.deepcopy(S),score,output_filename)
        
       
        #f_write.write('\n\n 10-taxon datatset/ model condition - '+str(mc)+'/R'+str(replicate)+'/ \n U = ')
        #f_write.write(str(U))
        #print 'writted to file'
        
#f_write.close()
   
  

'''

#PART 4
#It runs tests on the 10-taxon dataset  on 4 mode conditions, 3 replicates
#using MORE LOCAL quintets and reports the scores


for mc in model_coditions:
    for replicate in range(1,6):
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


        #print(S.as_ascii_plot())
        #I'm stopping at the second to last edge because I think in the postorder traversal the final edge is redundant due to Dendropy's "seed node" structure
        ### FORMATTED OUTPUT ########
        score=[]
        for i in range(len(ES)-1):
            (se,U,quintet)= InvariantScores.total_quintet_score_distance_kajori(S,i,treelist)
            taxon=InvariantScores.taxon_with_split_edge(S,ESBITS[i])
            # here i = index of edge in postorder_iteration
            score.append(se)
            print 'index =',i,' score = ',se, 'quintet = ',quintet,' taxon=',taxon, ' U =',U,'\n'

        output_filename='/Users/kajori/Box Sync/UIUC/Tandy/invariant_rooting/output_trees/'+mc+'_'+str(replicate)+'_10_taxon_short_quintet_tree_wih_score.trees'
        InvariantScores.my_print_tree(copy.deepcopy(S),score,output_filename)
'''
'''
#PART 5
#It runs tests on the Avian dataset  on 5 replicates
#using MORE LOCAL quintets and reports the scores
S = dendropy.Tree.get_from_path('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/Avian_0.5X/avian_species_tree.trees','newick') #kajori
S.deroot()
S.encode_splits()
S.update_splits()
DS = S.split_edges
ES = [e for e in S.postorder_edge_iter()]
INVDS = {v: k for k, v in DS.items()}
#dict that gives bitmasks in same order as edge list should be
ESBITS = [INVDS[ES[i]] for i in range(len(ES))]
for replicate in range(3,5): #replicate
    print '\n ********** '+str(replicate)+ '/ ********** '
    
    
    #combine the gene trees from the 1000 sub-dir
    #f_write = open('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/Avian_0.5X/0.5X-1000-true/combined_truegenetrees/truegenetrees_R'+str(replicate), 'w')
    #for gt in range((replicate-1)*1000+1,replicate*1000):
    #    f_read = open('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/Avian_0.5X/0.5X-1000-true/R'+str(replicate)+'/'+str(gt)+'/true.gt','r')
    #    f_write.write(f_read.read())
    #    f_read.close()
    #f_write.close()
    
    treelist =  dendropy.TreeList.get_from_path('/Users/kajori/Box Sync/UIUC/Tandy/Data_Set/Avian_0.5X/0.5X-1000-true/combined_truegenetrees/truegenetrees_R'+str(replicate),'newick') #kajori
    for i in range(len(treelist)):
            treelist[i].deroot()
    
    score=[]
    for i in range(len(ES)-1):
        (se,U,quintet)= InvariantScores.total_quintet_score_distance_kajori(copy.deepcopy(S),i,treelist)
        taxon=InvariantScores.taxon_with_split_edge(copy.deepcopy(S),ESBITS[i])
        # here i = index of edge in postorder_iteration
        score.append(se)
        print 'index =',i,' score = ',se, 'quintet = ',quintet,' taxon=',taxon, ' U =',U,'\n'


    
    output_filename='/Users/kajori/Box Sync/UIUC/Tandy/invariant_rooting/output_trees/'+str(replicate)+'_avian_short_quintet_tree_wih_score.trees'
    InvariantScores.my_print_tree(copy.deepcopy(S),score,output_filename)

'''