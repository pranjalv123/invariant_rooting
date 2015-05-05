#!/bin/bash

#Note that in this file, we are not making quartets, trees or astral clades because those have already been made. 

#$1 is the input file.  First make the quartets for the experiment:
#sh ~/code/quartet_scores/quartet-controller.sh ~/code/11_taxon_experiments/gene_tree_files/$1 ~/code/11_taxon_experiments/quartet_files/$1.quartets

#echo made quartets

#then run astral and capture its search space:
#java -jar ~/code/Astral\ 4.7.6/astral.4.7.6.jar -i ~/code/11_taxon_experiments/gene_tree_files/$1 -k searchspace >  ~/code/11_taxon_experiments/astral_clades/$1.astralclades -o ~/code/11_taxon_experiments/astral_species_trees/$1.astral_species_tree

#echo made astral tree and clades

#Note currently logs are empty.  Need to figure this out later
#java -jar ~/code/Astral\ 4.7.5/astral.4.7.5.jar -i ~/code/11_taxon_experiments/gene_tree_files/$1 -o ~/code/11_taxon_experiments/astral_species_trees/$1.astral_species_tree > ~/code/11_taxon_experiments/astral_species_logs/log.$1
echo making inequalities-only invariant tree from file $1

#Now run the invariant species tree algorithm
python ~/code/invariant_rooting/11TaxonRunIneqsOnly.py ~/code/11_taxon_experiments/quartet_files/$1.quartets ~/code/11_taxon_experiments/astral_clades/$1.astralclades  ~/code/11_taxon_experiments/invariant_species_trees_inequalities_only/$1.invariant_species_tree > ~/code/11_taxon_experiments/invariant_species_trees__inequalities_only_logs/log.$1

echo made invariant tree
#now compare missing branches for invariant trees
/Users/ruthdavidson/code/compare_trees/compareTrees.missingBranch /Users/ruthdavidson/code/11_taxon_experiments/11modeltree ~/code/11_taxon_experiments/invariant_species_trees_inequalities_only/$1.invariant_species_tree > ~/code/11_taxon_experiments/missing_branch_invariant_inequalities_only_logs/log.missingBranchInv.$1

echo comparing branches now
#now compare missing branches for astral trees
#/Users/ruthdavidson/code/compare_trees/compareTrees.missingBranch /Users/ruthdavidson/code/11_taxon_experiments/11modeltree ~/code/11_taxon_experiments/astral_species_trees/$1.astral_species_tree > ~/code/11_taxon_experiments/missing_branch_astral_logs/log.missingBranchAstral.$1

#Now compare the Astral and invariant trees

/Users/ruthdavidson/code/compare_trees/compareTrees.missingBranch ~/code/11_taxon_experiments/astral_species_trees/$1.astral_species_tree  ~/code/11_taxon_experiments/invariant_species_trees_inequalities_only/$1.invariant_species_tree > ~/code/11_taxon_experiments/missing_branch_astral_vs_inv_ineqs_only_logs/log.missingBranchAstral.$1

echo compared branches, now making quartets

#Finally get the quartets for the species tree

sh ~/code/quartet_scores/quartet-controller.sh ~/code/11_taxon_experiments/invariant_species_trees_inequalities_only/$1.invariant_species_tree ~/code/11_taxon_experiments/invtree_ineqs_only_quartet_files/$1.species.tree.quartets

#echo done with invariant species tree quartets

echo done with ineqs only  invariant species tree quartets, scoring tree

python ~/code/invariant_rooting/CheckScoresIneqsOnly.py ~/code/11_taxon_experiments/quartet_files/$1.quartets  ~/code/11_taxon_experiments/invtree_ineqs_only_quartet_files/$1.species.tree.quartets > ~/code/11_taxon_experiments/score_logs_inv_ineqs_only/log.$1

echo finished scoring tree

