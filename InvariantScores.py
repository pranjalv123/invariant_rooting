#inv3 gets the score for 3-taxon trees or 4-taxon trees.  Use for u1 the
#probability that the gene tree matches the species tree, then u2, u3 those that don't
def inv3(u1, u2, u3):
    a12 = min(u1-u2,0)
    a13 = min(u1-u3,0)
    score = abs(u2-u3) + (-1)*a12 + (-1)*a13
    return score


#inv51 is for the balanced species tree 5 leaves (((a,b),c),(d,e))
def inv51(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15):
    a12 = min(u1-u2,0)
    a14 = min(u1-u4,0)
    a25 = min(u2-u5,0)
    a45 = min(u4-u5,0)
    a57 = min(u5-u7,0)
    score1 = (-1)*a12 + (-1)*a14 + (-1)*a25 + (-1)*a45 + (-1)*a57
    score2 = abs(u14 - u15) + abs(u11-u15) + abs(u10-u15) + abs(u9-u12) + abs(u8-u15) + abs(u7-u15) + abs(u6-u12) + abs(u5-u12) + abs(u4-u13) + abs(u2-u3)
    score = score1 + score2
    return score

#inv52 is for the caterpillar tree 5 leaves (((a,b),c),d),e)
def inv52(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15): 
    a12 = min(u1-u2,0)
    a14 = min(u1-u4,0)
    a25 = min(u2-u5,0)
    a45 = min(u4-u5,0)
    a57 = min(u5-u7,0) 
    a32 = min(u3-u2,0)
    a36 = min(u3-u6,0)
    a65 = min(u6-u5,0)
    score1 = (-1)*a12 + (-1)*a14 + (-1)*a25 + (-1)*a45 + (-1)*a57 + (-1)*a32 + (-1)*a36 + (-1)*a65
    score2 = abs(u14-u15) + abs(u11-u15) + abs(u10-u15) + abs(u8-u15) + abs(u7-u15) + abs(u6-u9) + abs(u5-u12) + abs(u4-u13) + abs(u2-u3 + u9-u12)
    score = score1 + score2
    return score

#inv53 is for the pseudocaterpillar tree (((a,b),(d,e)),c)
def inv53(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15):
    a12 = min(u1-u2,0)
    a14 = min(u1-u4,0)
    a18 = min(u1-u8,0)
    a25 = min(u2-u5,0)
    a45 = min(u4-u5,0) 
    a85 = min(u8-u5,0)
    score1 = (-1)*a12 + (-1)*a14 + (-1)*a18 + (-1)*a25 + (-1)*a45 + (-1)*a85
    score2 = abs(u14-u15) + abs(u12-u15) + abs(u10-u15) + abs(u9-u15) + abs(u8-u11) + abs(u7-u15) + abs(u6-u15) + abs(u5-u15) + abs(u4-u13) + abs(u2-u3)
    score = score1 + score2
    return score


