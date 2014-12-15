import dendropy
import matrixmaker
import InvariantScores
import time

L = ['ECH', 'GAL', 'MAC', 'MON', 'ORN']

print L

Subsets = InvariantScores.powerset(L)

print len(Subsets)

start = time.time()

print 'Making Class instance SubsetPenalties'

A = InvariantScores.SubsetPenalties(L, Subsets, 'true-trees-quartets')

print 'setlist is'

print A.setlist
end0 = time.time()

print end0-start

print 'Making scored matrix'

print A.scored_matrix()

end1 = time.time()

print end1- end0

print 'making clades'

print A.clades()

end2 = time.time()

print end2-end1

print 'making tree'

print A.tree()

end3 = time.time()

print end3-end2
