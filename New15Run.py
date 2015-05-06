import dendropy
import matrixmaker
import InvariantScores
import time
import sys


quartetsfile = str(sys.argv[1])

cladesfile = str(sys.argv[2])

treefilename = str(sys.argv[3])

L = ['I', 'H', 'F', 'G', 'E', 'A', 'B', 'D', 'C', 'J', 'K','L','M','N','O']

print 'taxon set is'
print L


F = InvariantScores.FileClades(cladesfile)

print 'using input clades from file'
print F.filename

Subsets = F.clades()

print len(F.clades())
print 'is number of clades from input file'

#java -jar ~/code/Astral\ 4.7.5/astral.4.7.5.jar -i ~/code/Astral\ 4.7.5/true-trees-1X.tre -k searchspace_norun -o dud > Subsets

print len(Subsets)

print Subsets

start = time.time()

print 'Making Class instance SubsetPenalties'

A = InvariantScores.SubsetPenalties(L, Subsets, quartetsfile)

print 'time for class instance is'

#print A.setlist
end0 = time.time()

print end0-start

print 'Making scored matrix'

print A.scored_matrix()

end1 = time.time()

print 'time to make scored matrix is'

print end1- end0

print 'making clades'

print A.clades()

end2 = time.time()

print 'time to make clades was'

print end2-end1

print 'making tree'

speciestree =  InvariantScores.dendropy_clades_tree(L,A.clades())

print speciestree.as_newick_string()

#treefilename = quartetsfile +'.and.'+cladesfile+'.tre'

H = open(treefilename, 'w')

H.write(speciestree.as_newick_string()+';')

H.close()

end3 = time.time()

print 'time to make tree was'

print end3-end2
