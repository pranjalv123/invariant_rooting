import dendropy
import itertools
import collections
import pprint
import string
from dendropy import treecalc
import matrixmaker
import numpy as np
import copy #by kajori
import Queue #by kajori
import re  #by kajori
import math  #by kajori

#from dendropy import *
#nuclear option to not type dendropy. I think this is frowned upon
#pp = pprint.PrettyPrinter(indent=4)
#inv3 gets the score for 3-taxon trees or 4-taxon trees.  Use for u1 the
#probability that the gene tree matches the species tree, then u2, u3 those that don't

#this function penalty takes the values u_i and u_j and calculates the penalty score depending on the user_function
# The error str is printed which shows the invariants that are violated
def penalty(str_x,x,str_y,y,user_function):
    if (user_function=='diff'): score= (-1)*min(x-y,0)
    elif (user_function=='ratio'):score=max(float(y+1)/float(x+1),1)
    elif (user_function=='stat' and x>=y): score=0
    elif  (user_function=='stat' and x<y): 
        temp=(float(x+1)/float(y+1+x+1))
        #print 'temp ', temp, ' x =',x,'y =',y
        #temp_1=(float(y+1)/float(y+1+x+1))
        score = (float(-1*(x+1))* math.log(temp))# + (float(-1*(y+1))* math.log(temp_1)) #+(float(y)* math.log(float(y)/float(y+x)))+ 
    if (score>0 and user_function=='diff'): error_str= str_y +' > '+str_x +', '+ str(y) + '>' + str(x)+ ', penalty  = ' + str(score)
    elif (score>1 and user_function=='ratio'):error_str= str_y +' > '+str_x +', '+ str(y) + '>' + str(x)+ ', penalty  = ' + str(score)
    elif (score>0 and user_function=='stat'):error_str= str_y +' > '+str_x +', '+ str(y) + '>' + str(x)+ ', penalty  = ' + str(score)
    else:error_str=''
    #print user_function, score, x, y, error_str
    return score,error_str    

def inv3(u1, u2, u3):
    a12 = min(u1-u2,0)
    a13 = min(u1-u3,0)
    score = abs(u2-u3) + (-1)*a12 + (-1)*a13
    return score

# a version of the penalty without equality....
def inv3I(u1, u2, u3):
    a12 = min(u1-u2,0)
    a13 = min(u1-u3,0)
    score (-1)*a12 + (-1)*a13
    return score

'''
#inv51 is for the balanced species tree 5 leaves (((a,b),c),(d,e))
def inv51(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15):
    a12 = min(u1-u2,0)
    a14 = min(u1-u4,0)
    a25 = min(u2-u5,0)
    a45 = min(u4-u5,0)
    a57 = min(u5-u7,0)
    score1 = (-1)*a12 + (-1)*a14 + (-1)*a25 + (-1)*a45 + (-1)*a57
    #score2 = abs(u14 - u15) + abs(u11-u15) + abs(u10-u15) + abs(u9-u12) + abs(u8-u15) + abs(u7-u15) + abs(u6-u12) + abs(u5-u12) + abs(u4-u13) + abs(u2-u3)
    score = score1 #+ score2
   
    error_str = ''
    if score> 0:  error_str =  error_str+ 'violation  inv51 - balanced species tree '
    if (u2>u1):  error_str =  error_str+  '\n  u2 > u1 ,'+ str(u2)+'>' + str(u1)+ ', penalty  = ' + str(u2-u1)
    if (u4>u1):  error_str =  error_str+  '\n  u4 > u1 ,'+ str(u4)+'>'+ str(u1)+ ', penalty  = ' + str(u4-u1)
    if (u5>u2):  error_str =  error_str+  '\n  u5 > u2 ,'+ str(u5)+'>'+ str(u2)+ ', penalty  = ' + str(u5-u2)
    if (u5>u4):  error_str =  error_str+  '\n  u5 > u4 ,'+ str(u5)+'>'+ str(u4)+ ', penalty  = ' + str(u5-u4)
    if (u5>u7):  error_str =  error_str+  '\n  u7 > u5 ,'+ str(u5)+'>'+ str(u7)+ ', penalty  = ' + str(u5-u7)
    
    return score,error_str 
'''


#inv51_func_part2 is for the balanced species tree 5 leaves (((a,b),c),(d,e))
#inv51_func_part2 is used in additon to inv51_func
#it has additional equalities turned into inequalities
def inv51_func_part2(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15,user_function):
    a5_15,error_str_5_15 = penalty('u5',u5,'u15',u15,user_function) # abs(u7-u15) and u5>u7 therefore u5>15
    a2_15,error_str_2_15 = penalty('u2',u2,'u15',u15,user_function) # abs(u7-u15) and u5>u7 and u2>u5 therefore u2>15
    a4_15,error_str_4_15 = penalty('u4',u4,'u15',u15,user_function) # abs(u7-u15) and u5>u7 and u4>u5 therefore u4>15
    
    a2_12,error_str_2_12 = penalty('u2',u2,'u12',u12,user_function) # abs(u5-u12) and u2>u5 therefore u2>u12
    a4_12,error_str_4_12 = penalty('u4',u4,'u12',u12,user_function) # abs(u5-u12) and u4>u5 therefore u4>u12
    a1_13,error_str_1_13 = penalty('u1',u1,'u13',u13,user_function) # abs(u4-u13) and u1>u4 therefore u1>u13
    a13,error_str_13     = penalty('u1',u1,'u3',u3,user_function) # abs(u2-u3) and u1>u2 therefore u1>u3
    
    score = a5_15 +  a2_15 + a4_15 +a2_12 + a4_12 + a1_13 + a13

    error_str = error_str_5_15 +  error_str_2_15 + error_str_4_15 +error_str_2_12 + error_str_4_12 + error_str_1_13 + error_str_13 
    return score,error_str 

#inv51_func is for the balanced species tree 5 leaves (((a,b),c),(d,e))
#it is just like inv51 but can implement more more penalty function other than just 'diff'
def inv51_func(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15,user_function):
    a12,error_str_12 = penalty('u1',u1,'u2',u2,user_function)
    a14,error_str_14 = penalty('u1',u1,'u4',u4,user_function)
    a25,error_str_25 = penalty('u2',u2,'u5',u5,user_function)
    a45,error_str_45 = penalty('u4',u4,'u5',u5,user_function)
    a57,error_str_57 = penalty('u5',u5,'u7',u7,user_function)
    score1 = a12 + a14 + a25 + a45 + a57
    #score2 = abs(u14 - u15) + abs(u11-u15) + abs(u10-u15) + abs(u9-u12) + abs(u8-u15) + abs(u7-u15) + abs(u6-u12) + abs(u5-u12) + abs(u4-u13) + abs(u2-u3)
    score3,err_str3=inv51_func_part2(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15,user_function)
    score = score1 + score3 #+ score2
    error_str = error_str_12 +error_str_14+error_str_25+error_str_45+error_str_57 + err_str3
    if bool(error_str):error_str = 'violation  inv51 - balanced tree \n'+ error_str
    return score,error_str 

'''   
#inv52 is for the caterpillar tree 5 leaves ((((a,b),c),d),e)
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
    #score2 = abs(u14-u15) + abs(u11-u15) + abs(u10-u15) + abs(u8-u15) + abs(u7-u15) + abs(u6-u9) + abs(u5-u12) + abs(u4-u13) + abs(u2-u3 + u9-u12)
    score = score1 #+ score2
    
    error_str = ''
    if score> 0:  error_str =  error_str+ 'violation  inv52 - caterpillar tree'
    if (u2>u1):  error_str =  error_str+ '\n  u2 > u1 ,' + str(u2)+'>'+ str(u1) + ', penalty  = ' + str(u2-u1)
    if (u4>u1):  error_str =  error_str+ '\n  u4 > u1 ,' + str(u4)+'>'+ str(u1)+ ', penalty  = ' + str(u4-u1)
    if (u5>u2):  error_str =  error_str+ '\n  u5 > u2 ,' + str(u5)+'>'+ str(u2)+ ', penalty  = ' + str(u5-u2)
    if (u5>u4):  error_str =  error_str+ '\n  u5 > u4 ,' + str(u5)+'>'+ str(u4)+ ', penalty  = ' + str(u5-u4)
    if (u5>u7):  error_str =  error_str+ '\n  u5 > u7 ,' + str(u5)+'>'+ str(u7)+ ', penalty  = ' + str(u5-u7)
    if (u2>u3):  error_str =  error_str+ '\n  u2 > u3 ,' + str(u2)+'>'+ str(u3)+ ', penalty  = ' + str(u2-u3)
    if (u6>u3):  error_str =  error_str+ '\n  u6 > u3 ,' + str(u6)+'>'+ str(u3)+ ', penalty  = ' + str(u6-u3)
    if (u6>u5):  error_str =  error_str+ '\n  u6 > u5 ,' + str(u6)+'>'+ str(u5)+ ', penalty  = ' + str(u6-u5)
    
    return score,error_str 
'''

#inv52_func_part2 is for the balanced species tree 5 leaves (((a,b),c),(d,e))
#inv52_func_part2 is used in additon to inv52_func
#it has additional equalities turned into inequalities
def inv52_func_part2(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15,user_function):
    
    a5_15,error_str_5_15 = penalty('u5',u5,'u15',u15,user_function) # abs(u7-u15) and u5>u7 therefore u5>15
    a2_15,error_str_2_15 = penalty('u2',u2,'u15',u15,user_function) # abs(u7-u15) and u5>u7 and u2>u5 therefore u2>u15
    a1_15,error_str_1_15 = penalty('u1',u1,'u15',u15,user_function) # abs(u7-u15) and u5>u7 and u2>u5 and u1>u2 therefore u1>u15
    a3_15,error_str_3_15 = penalty('u3',u3,'u15',u15,user_function) # abs(u7-u15) and u5>u7 and u2>u5 and u3>u2 therefore u3>u15
    a4_15,error_str_4_15 = penalty('u4',u4,'u15',u15,user_function) # abs(u7-u15) and u5>u7 and u4>u5 therefore u4>u15
    a6_15,error_str_6_15 = penalty('u6',u6,'u15',u15,user_function) # abs(u7-u15) and u5>u7 and u6>u5 therefore u6>u15
    
    a39,error_str_39 = penalty('u3',u3,'u9',u9,user_function) # abs(u6-u9) and u3>u6 therefore u3>u9
    
    a6_12,error_str_6_12 = penalty('u6',u6,'u12',u12,user_function) #abs(u5-u12)  and u6>u5 therefore  u6>u12
    a2_12,error_str_2_12 = penalty('u2',u2,'u12',u12,user_function) #abs(u5-u12)  and u2>u5 therefore  u2>u12
    a4_12,error_str_4_12 = penalty('u4',u4,'u12',u12,user_function) #abs(u5-u12)  and u4>u5 therefore  u4>u12
    a1_12,error_str_1_12 = penalty('u1',u1,'u12',u12,user_function) #abs(u5-u12)  and u2>u5 and u1>u2 therefore  u1>u12
    a3_12,error_str_3_12 = penalty('u3',u3,'u12',u12,user_function) #abs(u5-u12)  and u2>u5 and u3>u2 therefore  u3>u12
    a1_13,error_str_1_13 = penalty('u1',u1,'u13',u13,user_function) # abs(u4-u13) and u1>u4 therefore u1>u13
    
    score = a5_15 +  a2_15 + a1_15+a3_15+ a4_15 +a6_15+ a39+ a6_12+ a2_12+ a4_12+ a1_12+ a3_12  + a1_13 
    error_str= error_str_5_15 +  error_str_2_15 + error_str_1_15+ error_str_3_15+ error_str_4_15 + error_str_6_15+ error_str_39+ error_str_6_12+ error_str_2_12+ error_str_4_12+ error_str_1_12+ error_str_3_12  + error_str_1_13 
    return score,error_str 
    
#inv52_func is for the balanced species tree 5 leaves (((a,b),c),(d,e))
#it is just like inv52 but can implement more more penalty function other than just 'diff'
def inv52_func(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15,user_function): 
    a12,error_str_12 = penalty('u1',u1,'u2',u2,user_function)
    a14,error_str_14 = penalty('u1',u1,'u4',u4,user_function)
    a25,error_str_25 = penalty('u2',u2,'u5',u5,user_function)
    a45,error_str_45 = penalty('u4',u4,'u5',u5,user_function)
    a57,error_str_57 = penalty('u5',u5,'u7',u7,user_function)
    a32,error_str_32 = penalty('u3',u3,'u2',u2,user_function)
    a36,error_str_36 = penalty('u3',u3,'u6',u6,user_function)
    a65,error_str_65 = penalty('u6',u6,'u5',u5,user_function)
    score1 = a12 + a14 + a25 + a45 + a57 + a32 + a36 + a65
    #score2 = abs(u14-u15) + abs(u11-u15) + abs(u10-u15) + abs(u8-u15) + abs(u7-u15) + abs(u6-u9) + abs(u5-u12) + abs(u4-u13) + abs(u2-u3 + u9-u12)
    score3,err_str3=inv52_func_part2(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15,user_function)
    score = score1 + score3 #+ score2
    error_str = error_str_12 +error_str_14+error_str_25+error_str_45+error_str_57+ error_str_32 + error_str_36 + error_str_65 + err_str3
    if bool(error_str):error_str = 'violation  inv52 - caterpillar tree \n'+ error_str
    return score,error_str 

'''
#inv53 is for the pseudocaterpillar tree (((a,b),(d,e)),c)
def inv53(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15):
    a12 = min(u1-u2,0)
    a14 = min(u1-u4,0)
    a18 = min(u1-u8,0)
    a25 = min(u2-u5,0)
    a45 = min(u4-u5,0) 
    a85 = min(u8-u5,0)
    score1 = (-1)*a12 + (-1)*a14 + (-1)*a18 + (-1)*a25 + (-1)*a45 + (-1)*a85
    #score2 = abs(u14-u15) + abs(u12-u15) + abs(u10-u15) + abs(u9-u15) + abs(u8-u11) + abs(u7-u15) + abs(u6-u15) + abs(u5-u15) + abs(u4-u13) + abs(u2-u3)
    score = score1 #+ score2
    
    error_str = ''
    if score> 0: error_str =  error_str+ 'violation  inv53 -  pseudocaterpillar tree'
    if (u2>u1):  error_str =  error_str+ '  u2 > u1 ,'+ str(u2)+'>'+ str(u1)+ ', penalty  = ' + str(u2-u1)
    if (u4>u1):error_str =  error_str+ '  u4 > u1 ,'+ str(u4)+'>'+ str(u1)+ ', penalty  = ' + str(u4-u1)
    if (u8>u1):error_str =  error_str+ '  u8 > u1 ,'+ str(u8)+'>'+ str(u1)+ ', penalty  = ' + str(u8-u1)
    if (u5>u2):error_str =  error_str+ '  u5 > u2 ,'+ str(u5)+'>'+ str(u2)+ ', penalty  = ' + str(u5-u2)
    if (u5>u4):error_str =  error_str+ '  u5 > u4 ,'+ str(u5)+'>'+ str(u4)+ ', penalty  = ' + str(u5-u4)
    if (u5>u8):error_str =  error_str+ '  u5 > u8 ,'+ str(u5)+'>'+ str(u8)+ ', penalty  = ' + str(u5-u8)
    
    return score,error_str    
'''

#inv53_func_part2 is for the pseudocaterpillar tree (((a,b),(d,e)),c)
#inv53_func_part2 is used in additon to inv53_func
#it has additional equalities turned into inequalities
def inv53_func_part2(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15,user_function):
    
    a1_11,error_str_1_11 = penalty('u1',u1,'u11',u11,user_function) # abs(u8-u11) and u1>u8 therefore u1>u11
    a2_15,error_str_2_15 = penalty('u2',u2,'u15',u15,user_function) # abs(u5-u15) and u2>u5 therefore u2>u15
    a4_15,error_str_4_15 = penalty('u4',u4,'u15',u15,user_function) # abs(u5-u15) and u4>u5 therefore u4>u15
    a8_15,error_str_8_15 = penalty('u8',u8,'u15',u15,user_function) # abs(u5-u15) and u8>u5 therefore u8>u15
    a1_15,error_str_1_15 = penalty('u1',u1,'u15',u15,user_function) # abs(u5-u15) and u2>u5 and u1>u2 therefore u1>u15
    
    a1_13,error_str_1_13 = penalty('u1',u1,'u13',u13,user_function) # abs(u4-u13) and u1>u4 therefore u1>u13
    a13,error_str_13 = penalty('u1',u1,'u3',u3,user_function) # abs(u2-u3) and u1>u2therefore u1>u3

    score =a1_11 +a2_15 + a4_15 +  a8_15 + a1_15+ a1_13+ a13 
    error_str_= error_str_1_11 +error_str_2_15 + error_str_4_15 +  error_str_8_15 + error_str_1_15+ error_str_1_13+ error_str_13 
    return score,error_str 
    
#inv53_func is for the pseudocaterpillar tree (((a,b),(d,e)),c)
def inv53_func(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15,user_function):
    a12,error_str_12 = penalty('u1',u1,'u2',u2,user_function)
    a14,error_str_14 = penalty('u1',u1,'u4',u4,user_function)
    a18,error_str_18 = penalty('u1',u1,'u8',u8,user_function)
    a25,error_str_25 = penalty('u2',u2,'u5',u5,user_function)
    a45,error_str_45 = penalty('u4',u4,'u5',u5,user_function)
    a85,error_str_85 = penalty('u8',u8,'u5',u5,user_function)
    score1 = a12 + a14 + a18 + a25 + a45 + a85
    #score2 = abs(u14-u15) + abs(u12-u15) + abs(u10-u15) + abs(u9-u15) + abs(u8-u11) + abs(u7-u15) + abs(u6-u15) + abs(u5-u15) + abs(u4-u13) + abs(u2-u3)
    score3,err_str3=inv52_func_part2(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15,user_function)
    score = score1 + score3  #+ score2
    error_str = error_str_12 + error_str_14 + error_str_18 + error_str_25 + error_str_45 + error_str_85+ err_str3
    if bool(error_str):error_str = 'violation  inv53 - pseudocaterpillar tree \n'+ error_str
    return score,error_str 
    
#re-verified 4/15/15 that t_1 are in the same order as in appendix of allman-rhodes-degnan
#l=list of taxon labels as strings from listofgenetreeschoices,listofgenetrees=dendropy TreeList
def dist_counter(l,listofgenetrees):
    M = listofgenetrees.taxon_set
    t1 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[1]+'),' +l[2]+ ',(' +l[3]+ ',' +l[4]+ '));\n', schema = 'newick', taxon_set = M) 
    t2 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[1]+'),' +l[3]+ ',(' +l[2]+ ',' +l[4]+ '));\n', schema = 'newick', taxon_set = M)
    t3 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[1]+'),' +l[4]+ ',(' +l[2]+ ',' +l[3]+ '));\n', schema = 'newick', taxon_set = M)
    t4 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[2]+'),' +l[1]+ ',(' +l[3]+ ',' +l[4]+ '));\n', schema = 'newick', taxon_set = M)
    t5 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[2]+'),' +l[3]+ ',(' +l[1]+ ',' +l[4]+ '));\n', schema = 'newick', taxon_set = M)
    t6 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[2]+'),' +l[4]+ ',(' +l[1]+ ',' +l[3]+ '));\n', schema = 'newick', taxon_set = M)
    t7 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[3]+'),' +l[1]+ ',(' +l[2]+ ',' +l[4]+ '));\n', schema = 'newick', taxon_set = M)
    t8 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[3]+'),' +l[2]+ ',(' +l[1]+ ',' +l[4]+ '));\n', schema = 'newick', taxon_set = M)
    t9 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[3]+'),' +l[4]+ ',(' +l[1]+ ',' +l[2]+ '));\n', schema = 'newick', taxon_set = M)
    t10 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[4]+'),' +l[1]+ ',(' +l[2]+ ',' +l[3]+ '));\n', schema = 'newick', taxon_set = M)
    t11 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[4]+'),' +l[2]+ ',(' +l[1]+ ',' +l[3]+ '));\n', schema = 'newick', taxon_set = M)
    t12 = dendropy.Tree.get_from_string('[&U] ((' +l[0]+ ',' +l[4]+'),' +l[3]+ ',(' +l[1]+ ',' +l[2]+ '));\n', schema = 'newick', taxon_set = M)
    t13 = dendropy.Tree.get_from_string('[&U] ((' +l[1]+ ',' +l[2]+'),' +l[0]+ ',(' +l[3]+ ',' +l[4]+ '));\n', schema = 'newick', taxon_set = M)
    t14 = dendropy.Tree.get_from_string('[&U] ((' +l[1]+ ',' +l[3]+'),' +l[0]+ ',(' +l[2]+ ',' +l[4]+ '));\n', schema = 'newick', taxon_set = M)
    t15 = dendropy.Tree.get_from_string('[&U] ((' +l[1]+ ',' +l[4]+'),' +l[0]+ ',(' +l[2]+ ',' +l[3]+ '));\n', schema = 'newick', taxon_set = M)
    counter1 = [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15]
    counter2 = [0 for i in range(15)]
    return [counter1, counter2]

# I should merge this next with dist_counter after I test them both
def get_dist(l,treelist):
    #treelist = dendropy.TreeList(inputtreelist)
    dist = dist_counter(l, treelist)
    #pp.pprint(dist)
    #viztrees = [dist[0][k].as_string("newick") for k in range(15)]
    #print viztrees
    garbage = []
    incompletes = []
    labels = [n.label for n in treelist.taxon_set]
    for i in range(len(treelist)):
        smalltree = copy.deepcopy(treelist[i])
        if len(list(smalltree.leaf_nodes())) == len(labels):
            smalltree.retain_taxa_with_labels(l)
            smalltree.deroot()
            smalltree.update_splits()
            trace = [dendropy.treecalc.symmetric_difference(smalltree,dist[0][j]) for j in range(15)]
            if trace.count(0) > 1:
                garbage.append(['error: too many zeroes', smalltree.as_string("newick"), smalltree,  treelist[i], trace])
            elif 0 in trace:
                dist[1][trace.index(0)] = dist[1][trace.index(0)] + 1
            else:
                garbage.append([smalltree.as_newick_string(), treelist[i].as_newick_string(), min(trace)])
        else:
            incompletes.append(str(i)+'th gene tree incomplete')
    return [dist, garbage, incompletes]

#S is a rooted species tree, l is a list of five taxa labels like in get_dist, etc.Picking the edge and the quintet are unresolved
def basic_score_quintet(S,l,treelist,user_function):
    T = dendropy.Tree(S)
    T.retain_taxa_with_labels(l)
    
    T.ladderize(ascending=False)
    #here change the order of the entries of l to match what happened to them in ladderize
    rooted_shape_str = T.as_newick_string()
    rooted_shape = ''
    y = ['(', ')', ',']
    for i in range(len(rooted_shape_str)):
        if rooted_shape_str[i] in y:
            rooted_shape = rooted_shape + rooted_shape_str[i]
    shapes = ['(((,),),(,))', '((((,),),),)', '(((,),(,)),)']
    
    #Here is where l needs to be reordered to match order of T after ladderizing and retain_taxa blah blah 
    T_temp=copy.deepcopy(T.as_newick_string())
    #re.sub(':[^,]+,', '', T_temp)
    #print '\n re.sub =',re.sub(':[^,]+,', ',', T_temp)
    temp=re.sub(':[^\)]+', '',re.sub(':[^,]+,', ',', T_temp)) 
    temp=re.sub('\(','',re.sub('\)', '',temp))
    #print '\n ************************************ \n T_temp =',T_temp
    #print  '\n re    =',temp
    ordered_quintet=[i for i in re.split(',',temp)]
    #print  '\n ordered_quintet    =',ordered_quintet
    assert((set(l) | set(ordered_quintet))==set(l))
    l=ordered_quintet
   
    U = get_dist(l,treelist)
    [u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15] = U[0][1] 
    scorefuncs= [inv51_func, inv52_func, inv53_func]
  
    if (user_function!='diff'):
        if rooted_shape == shapes[0]:
            score,error_str = inv51_func(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15,user_function)
        elif rooted_shape == shapes[1]:
            score,error_str = inv52_func(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15,user_function)
        elif rooted_shape == shapes[2]:
            score,error_str = inv53_func(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15,user_function)
        #score_function = scorefuncs[shapes.index(rooted_shape)]
        #score = score_function(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
        else:
            print 'error: quintet topology ' + T.as_newick_string() + '  not in list'
            score = 0
            error_str=''
    
    elif(user_function=='diff'):
        if rooted_shape == shapes[0]:
            score,error_str = inv51(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
        elif rooted_shape == shapes[1]:
            #print 'u2 =',u2,'u3=',u3
            score,error_str = inv52(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
        elif rooted_shape == shapes[2]:
            score,error_str= inv53(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
        #score_function = scorefuncs[shapes.index(rooted_shape)]
        #score = score_function(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
        else:
            print 'error: quintet topology ' + T.as_newick_string() + '  not in list'
            score = 0
            error_str=''
    #print ' 1. l= ',l,' U=',U[0][1] #RUTH - U changes with l
    return score,U[0][1],error_str 

#S is a rooted species tree, l is a list of five taxa labels like in get_dist, etc.Picking the edge and the quintet are unresolved
# the only difference between basic_score_quintet_kajori and basic_score_quintet is the last return sentence. I have added it to associated the
#invariant scores with the edges
'''
def basic_score_quintet_kajori(S,l, treelist):
    T = dendropy.Tree(S)
    T.retain_taxa_with_labels(l)
    T.ladderize(ascending=False)
    rooted_shape_str = T.as_newick_string()
    rooted_shape = ''
    y = ['(', ')', ',']
    for i in range(len(rooted_shape_str)):
        if rooted_shape_str[i] in y:
            rooted_shape = rooted_shape + rooted_shape_str[i]
    shapes = ['(((,),),(,))', '((((,),),),)', '(((,),(,)),)']
    U = get_dist(l,treelist)
'''

#basic_score...kajori was for asserting freqs always added to 1000
#uncommented to get tests to pass 4/24/15
    
#score_quintet takes a rooted quintet tree Q as input-is ladderized output of get_rooted_quintet
#def score_quintet(Q,l, treelist):
    #rooted_shape_str = Q.as_newick_string()
    #rooted_shape = ''
    #y = ['(', ')', ',']
    #for i in range(len(rooted_shape_str)):
        #if rooted_shape_str[i] in y:
            #rooted_shape = rooted_shape + rooted_shape_str[i]
    #shapes = ['(((,),),(,))', '((((,),),),)', '(((,),(,)),)']
    #U = get_dist(l,treelist)
    #[u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15] = U[0][1]
    #scorefuncs= [inv51, inv52, inv53]
    #if rooted_shape in shapes:
        #score_function = scorefuncs[shapes.index(rooted_shape)]
        #score = score_function(u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15)
    #else:
        #print 'error: quintet topology ' + Q.as_newick_string() + '  not in list'
        #score = 0
    #return score

#S is unrooted Astral tree, l is list of five taxa labels from S, i is the index of edge we want to root at in set 2n-3
def get_rooted_quintet(S,l,i):
    T = dendropy.Tree(S)
    edgelist = [e for e in T.postorder_edge_iter()]
    root_edge = edgelist[i]
    T.reroot_at_edge(root_edge)
    T.retain_taxa_with_labels(l)
    T.ladderize(ascending=False)
    return T

#This function scores an edge based on all five-element subsets of the leaf set of the species tree.  This will be replaced or 
#modified to score an edge based on 2n-3 five-elements subsets where n is the numbers of species.  i indexes in the edges in 
#find _best_edge_by_total_quintet_score(S,treelist):
def total_quintet_score(S,i,treelist):
    #T needs to be S rooted at edge i!!!!
    T = dendropy.Tree(S)
    T.deroot() 
    edgelist = [e for e in T.postorder_edge_iter()]
    #root_edge = edgelist[i]
    T.reroot_at_edge(edgelist[i])
    L = [n.taxon.label for n in T.leaf_nodes()]
    #print L
    Q1 = list(itertools.combinations(L,5))
    #print Q1 
    Q = [list(Q1[k]) for k in range(len(Q1))]
    print Q
    #Qtreeslist =[get_rooted_quintet(T,l,i) for l in Q]
    scorelist = []
    print 'calculating score for ' + str(i+1) + 'th edge'
    #for j in range(len(Qtreeslist)):

    for j in range(len(Q)):
        H,U,str = basic_score_quintet(T, Q[j],treelist)
        scorelist.append(basic_score_quintet(T, Q[j],treelist))
    totalscore = sum(scorelist)
    #print scorelist
    return totalscore,U,str

    
############### KAJORI BEGIN   ###############   
#This function scores an edge based on all five-element subsets of the leaf set of the species tree.  This will be replaced or 
#modified to score an edge based on 2n-3 five-elements subsets where n is the numbers of species.  i indexes in the edges in 
#find _best_edge_by_total_quintet_score(S,treelist):
# total_quintet_score_kajori is a dfs based approach to find the quintet
def total_quintet_score_kajori(S,i,treelist):
    #print 'Inside  total_quintet_score'
    #T needs to be S rooted at edge i!!!!
    T = dendropy.Tree(S)
    node_label = {}
    #print 'node_labels'
    #print node_label
    nodelist = [n for n in T.postorder_node_iter()]
    for n in nodelist:
    	node_label[n.oid]=[]
    	if (n.is_leaf()): 
    		node_label[n.oid].append(n.taxon.label)
    	else :
    		for child in n.child_nodes():
    			node_label[n.oid].append(max(node_label[child.oid]))
		#print 'node_label[node_id]= ',node_label[n.oid]
		#print ' %%%%%%%%%%%%%%%%%  \n'
			
    #print ' Next loop after for loop'
    edgelist = [e for e in T.postorder_edge_iter()]
    e=edgelist[i]
    quintet=[]
    visited_node=[]
    q = Queue.Queue()
    q.put(e.tail_node)
    q.put(e.head_node)
    while (q.empty()==False and len(quintet)<5):
        node_id=q.get()
        if (node_id is not None):
            #print 'node_id',node_id,'type(node_id)',type(node_id),'node_id.oid',type(node_id.oid)
            #for each node label
            for label in node_label[node_id.oid]:
                if label not in quintet and  len(quintet)<5:
                    quintet.append(label)
            #for each node 
            if node_id.oid not in visited_node:
                visited_node.append(node_id.oid)
                for n in node_id.get_adjacent_nodes():
                    if n.oid not in visited_node:
                        q.put(n)
                        
    # quintet
    T.reroot_at_edge(e)
   
   
    H,U = basic_score_quintet(T,quintet,treelist)
    #print ' score = ',H, 'quintet = ',quintet, ' U =',U
    assert sum(U)==1000
    return (H,U,quintet)



# this function takes as input postion of an edge in postorder iterations and returns taxons on its both sides
#unit test written for it
def nearest_quintet_kajori(S,i):
    T= dendropy.Tree(S)
    T.deroot()
    T.encode_splits()
    T.update_splits()
    DS = T.split_edges
    edgelist = [e for e in T.postorder_edge_iter()]
    INVDS = {v: k for k, v in DS.items()}
    
    split_hash_bitmask=INVDS[edgelist[i]]
    no_leaf_nodes=T.leaf_nodes()
    #print ' split_hash_bitmask =', split_hash_bitmask
    taxon_1=taxon_with_split_edge(T,split_hash_bitmask) #input is a hash_bitamsk input_no, finds tha taxon set associated with it
    #print taxon_1,len(taxon_1)
    
    #boundary conditions check
    if (len(taxon_1)==len(no_leaf_nodes)):
        taxon_2=[]
    elif (len(taxon_1)==0 ):
         taxon_2=[node_id.taxon.label for node_id in no_leaf_nodes]
    else:
        T.prune_taxa_with_labels(taxon_1)
        no_leaf_nodes=T.leaf_nodes()
        taxon_2=[node_id.taxon.label for node_id in no_leaf_nodes]
        
    #print 'taxon_1',taxon_1, ' taxon_2',taxon_2
    return taxon_1,taxon_2
    
# this function sort_list_distance_kajori takes as input a list of taxons 
#and sorts them on the basis of their distance from the root
#unit test written for it
def sort_list_distance_kajori(T,taxon):
    node_list=[n for n in T.level_order_node_iter()]
    quintet=[]
    for node in node_list:
        if (node.is_leaf() and node.taxon.label in taxon and len(quintet)<5):
            quintet.append(node.taxon.label)
    return quintet

# this function finds the quinters closest to the root
# S = TreeList ,i = index of edge under investigation  dir= tail/head
#total_quintet_score_distance_kajori along with nearest_quintet_kajori endsures that the edge is there in the induced subgraph 
def total_quintet_score_distance_kajori(S,i,treelist,user_function):
    T= dendropy.Tree(S)
    T.deroot()
    taxon_1,taxon_2=nearest_quintet_kajori(copy.deepcopy(S),i)
    
    edgelist_postorder = [e for e in T.postorder_edge_iter()]
    pos=0
    
    '''
    for e in edgelist_postorder:
        
        e.length=1
        if( i==19):
            print ' pos = ',pos ,' edge =', e.oid, '  ',s
            pos=pos+1
            
    '''
    
    #HACK 2 #TALK TO RUTH
    if (i!=len(edgelist_postorder)-1): T.reroot_at_edge(edgelist_postorder[i])
    #print 'taxon_1,taxon_2',taxon_1,taxon_2
    
    quintet_1=sort_list_distance_kajori(T,taxon_1)
    quintet_2=sort_list_distance_kajori(T,taxon_2)
    
    #print 'quintet_1,quintet_2',quintet_1,quintet_2
    #HACK 3
    if (i==len(edgelist_postorder)-1): 
        if (quintet_2==[]) : 
            quintet=quintet_1
        else: quintet= quintet_2
        T= dendropy.Tree(S)     
        H,U,error_str = basic_score_quintet(T,quintet,treelist,user_function)
        #assert sum(U)==1000
        return (H,U,quintet)       
    
    
    quintet=[]
    (flag_1,flag_2,pos_1,pos_2)=(0,0,0,0)
    while (pos_1+pos_2 < 5 and pos_1 < len(quintet_1)  and pos_2 < len(quintet_2) ):
        #print 'quintet_1[pos_1] ', quintet_1[pos_1], T.find_node_with_taxon_label(quintet_1[pos_1]).taxon.label
        #print ' .distance_from_root() ', T.find_node_with_taxon_label(quintet_1[pos_1]).distance_from_root()
        dist_1=T.find_node_with_taxon_label(quintet_1[pos_1]).distance_from_root()
        dist_2=T.find_node_with_taxon_label(quintet_2[pos_2]).distance_from_root()
        if(dist_1<dist_2):
            quintet.append(quintet_1[pos_1])
            pos_1=pos_1+1
            flag_1=1
        else:
            quintet.append(quintet_2[pos_2])
            pos_2=pos_2+1
            flag_2=1
    while (pos_1+pos_2 < 5 and pos_1 < len(quintet_1)):
        quintet.append(quintet_1[pos_1])
        pos_1=pos_1+1
        flag_1=1
    while (pos_1+pos_2 < 5 and pos_2 < len(quintet_2)):
        quintet.append(quintet_2[pos_2])
        pos_2=pos_2+1
        flag_2=1
    
    #print 'pos_1,pos_2,len(quintet)', pos_1,pos_2,len(quintet)
    if(flag_1==0):
        quintet[4]=quintet_1[0]
    elif (flag_2==0):
        quintet[4]=quintet_2[0]
    #print 'final quintet',quintet,'len',len(quintet)
    
    
    T= dendropy.Tree(S)     
    T.reroot_at_edge(edgelist_postorder[i])
    H,U,error_str = basic_score_quintet(T,quintet,treelist,user_function)
    #assert sum(U)==1000
    return (H,U,quintet)            
    
    #find the edge id associated with the outlier , label indicates the taxon label of the outlier node
def find_edge_associated_with_outlier_kajori(S, label):
    #print 'find_edge_associated_with_outlier'
    T = dendropy.Tree(S)
    node_oid=T.find_node_with_taxon_label(label)
    edge=node_oid.incident_edges()
    edgelist = [e.oid for e in T.postorder_edge_iter()]
    pos_edge=edgelist.index(edge[0].oid) 
    return pos_edge
    # Removes nodes that belong to the input_bitmask from S
   
   #input is a bitamsk input_no, it removes the edge associated with it
def remove_split_bitmask_kajori(T,input_no):
    T.encode_splits()
    T.update_splits()
    DS = T.split_edges
    if input_no not in DS.keys():
        print 'Split_BitMask NOT PRESENT IN INPUT TREE',DS.keys()
        return
    e=DS[input_no]  #returns the edge associated with the hash bitmask input_no
    #print 'type (e)',type(e)
    node=T.mrca(split_bitmask=e.split_bitmask)
    T.prune_subtree(node)
    #T.print_plot()
    return T
    
    #input is a hash_bitamsk input_no, finds tha taxon set associated with it
def taxon_with_split_edge(S,input_no):
    T= dendropy.Tree(S)
    T.deroot()
    T.encode_splits()
    T.update_splits()
    DS = T.split_edges
    if input_no not in DS.keys():
        print 'Split_BitMask NOT PRESENT IN INPUT TREE'
        return
    e=DS[input_no]  #returns the edge associated with the hash bitmask input_no
    #print 'type (e)',type(e)
    node_id=T.mrca(split_bitmask=e.split_bitmask)
    #Queue Initialization
    taxon_set=[]
    q = Queue.Queue()
    q.put(node_id)
    
    while (q.empty()==False):
        node_id=q.get()
        if (node_id.is_leaf()):
            taxon_set.append(node_id.taxon.label)
        else:
            for n in node_id.child_nodes():
                q.put(n)
    #print taxon_set
    return taxon_set

        
#given a quintet finds the edges induced by the quintet
#outputs the edge indexes in the post_order iteration
def find_induced_edge_indexes(S,quintet):
    T= dendropy.Tree(S)
    T.encode_splits()
    T.update_splits()
    DS = T.split_edges
    ES = [e for e in T.postorder_edge_iter()]
    INVDS = {v: k for k, v in DS.items()}
    
    ESBITS = [INVDS[ES[i]] for i in range(len(ES))] #dict that gives bitmasks in same order as edge list should be
    edges_included=[]
    taxon_set=[n.taxon.label for n in T.leaf_nodes()]
    pos=0
    for split_hash_bitmask in ESBITS:
        taxon_1=taxon_with_split_edge(T,split_hash_bitmask)
        #print 'taxon_1 ',taxon_1, ' pos ',pos
        #ensures that the elements from the quintet are on both sides on the edge
        #print 'taxon_1 =',taxon_1,' pos ',pos
        if ((len(set(taxon_1) & set(quintet))>0) and ( len((set(taxon_set)-set(taxon_1)) & set(quintet))>0)):
            edges_included.append(pos)
            #print 'included '
        #else:
            #print 'taxon_1 =',taxon_1
        #corner cases
        #THIS HACK WORKS ONLY FOR THE 10-taXON DATASET. DISCUSS WITH RUTH ABOUT THIS
        elif (len(set(taxon_set)-set(taxon_1)) ==0) and ( '0' in quintet) :
            edges_included.append(pos)
            #print 'included '
        pos=pos+1
    #print 'edges_included =',edges_included
    return edges_included
            
#given a quintet returns the scores of the induced edges in the postorder_edge_iter
#if the edge index in the postorder_edge_iter is not induced the score is reported as -99
def edge_score_on_quintet(S,quintet,treelist,user_function):
    T= dendropy.Tree(S)
    ES = [e for e in T.postorder_edge_iter()]
    edge_list=find_induced_edge_indexes(copy.deepcopy(S),quintet)
    #print 'len induced_edge ', len(edge_list), ' list ', edge_list
    #print 'before scoring'
    #print(S.as_ascii_plot())
    score=[]
    violations=''
   
    for index in range(0,len(ES)):
        T_copy=copy.deepcopy(T)
        if (index in edge_list):
            #print 'index=',index, type(ES[index])
            T_copy.reroot_at_edge(ES[index])
            H,U,error_str= basic_score_quintet(T_copy,quintet,treelist,user_function)
            score.append(H)
           
            
            if bool(error_str): print '\n For edge index '+ str(index)+ ' \n  '+error_str +' \n U= ',U
        else:
            score.append(-99)
            #assert sum(U)==1000
     
    #print violations
    return score,U

#it takes as input a tree and prints the tree with the edge scores 
#if score=True, output= edge_index_in_post_order_iteration/score_of _that_edge
#else output= edge_index_in_post_order_iteration
#since it was not possible to print the edge labels in ascii,
# they are printed as node labels 
# the node labels show the labels of the edges from its parent node
def my_print_tree(S_original,score,file_name):
    
    S=copy.deepcopy(S_original) #so that we do not change the labels of the original tree
    node_list=[n for n in S.postorder_node_iter()]
    #print 'len', len(score), len(node_list)
    for pos in range(0,len(node_list)-1): #I'm stopping at the second to last edge because I think in the postorder traversal the final node is redundant due to Dendropy's 'seed node' structure
        #print pos,score[pos],type(node_list[pos])
       
        if( score[pos]==-99 and node_list[pos].is_leaf()): 
            node_list[pos].taxon.label=node_list[pos].taxon.label+'/'+str(pos)
        elif( score[pos]!=-99 and node_list[pos].is_leaf()): 
            node_list[pos].taxon.label=node_list[pos].taxon.label+'/'+str(pos)+'/'+str(score[pos])
        elif( score[pos]==-99 and node_list[pos].is_internal()): 
            node_list[pos].label=str(pos)
        else:
            node_list[pos].label=str(pos)+'/'+str(score[pos])
    S.print_plot(show_internal_node_labels=True)
    S.write_to_path(file_name,'newick')
    return S

#it takes as input a tree, the quintet, 
#the scores associated with those edges
# Output : it removes those edges wich are not part of the induced edge subset.
def format_tree(S,quintet,score,filename):
    S.retain_taxa_with_labels(quintet,delete_outdegree_one=False)
    #print score
    while -99 in score: score.remove(-99)
    #print score
    node_list=[n for n in S.postorder_node_iter()]
    e=[e for e in S.postorder_edge_iter()]
    #print ' len = ',len(node_list),len(e),len(score)
    
    
    for pos in range(len(node_list)-1): # the last edge is the edge added by dendropy. So I am not using it
        if (node_list[pos].is_leaf()):node_list[pos].taxon.label=node_list[pos].taxon.label+'/'+str(score[pos])
        else:node_list[pos].label=str(score[pos])
    S.print_plot(show_internal_node_labels=True)
    S.write_to_path(filename,'newick')
    return S
    
#takes a tree as input and prints the newick string formatting the branch lengths
def print_newick_string(S):
    T_temp=copy.deepcopy(S.as_newick_string())
    temp=re.sub(':[^\)^,]+', '', T_temp)
    print temp 

#takes a 5-taxon tree (T_5) as input and finds whether the edges beside the root are getting the best edge 
#if yes return 1 if no other edge other than the root edges get the min score
#if yes return 2 if no other edge other than the root edges get the min score
#return 0 if both the root edges > min score
#return 1 both root edges have same score and only the root edge gets the min score
#return 2 both root edges have same score but there is some other node with the min score
#return 3 both root edges do not have same score but only the root edge gets the min score
#return 4 both root edges do not  have same score and there is some other node with the min score
#this function will be helpful for the discussion section
#Pay attention to the returned values
#for 0 - return root score+ min score
#rest root score+ number of edges other than root edges with  in score
def scores_edges_root(T_5,score):
    node_list=[ n for n in T_5.level_order_node_iter()]
    #check if the root is the node_list[0] - check passed because node_list[0] has bel None
    #for n in  node_list:
    #    if (n.is_leaf()): print n.taxon.label
    #    else:print n.label
    root=node_list[0]
    root_edges=[ e for e in root.get_incident_edges()]
    edge_list=[e for e in T_5.postorder_edge_iter()]
    #print ' scores_edges_root len(score) = ',len(score),' len(edge_list) ',len(edge_list)
    my_dict={}
    for index in range(len(score)):
        my_dict[edge_list[index]]=score[index]
    #print ' dictionary ',  my_dict.keys() , '\n  ***' ,root_edges
    print 'scores of root edges'
    root_score=[]
    for e in set(root_edges) & set(my_dict.keys()):
        root_score.append(my_dict[e])
    if ( min(root_score) > min(score)):
        print 'root is not properly scored'
        return 0,'root score = '+ str(root_score[0]), min(score)
    if (root_score[0]==root_score[1]):
        print ' root gets min score and both root edges get same score'
        if (score.count(min(score))==2): 
            print ' No edge other than root gets the best score '
            return 1,'root score = '+ str(root_score[0]),0
        else:  return 2,'root score = '+ str(root_score[0]), score.count(min(score)) - 2
    else:
        print ' root gets min score'
        if (score.count(min(score))==1): 
            print ' No edge other than root gets the best score '
            return 3, 'root score = '+ str(root_score[0]),0
        else:  return 4, 'root score = '+ str(root_score[0]),score.count(min(score))-1
        
        
        
    
############### KAJORI END ###############   
    

def find_best_edge_by_total_quintet_score(S,treelist):
    T = dendropy.Tree(S) 
    T.deroot()
    edgelist = [e for e in T.postorder_edge_iter()]
    numedges = len(edgelist)
    Scores = [total_quintet_score(S,j,treelist) for j in range(numedges)] 
    best_edge = Scores.index(min(Scores))
    T.reroot_at_edge(edgelist[best_edge])
    return T

#quartetsfile is a string name of file, file must be in working directory
class QuartetsInfo:
    def __init__(self,quartetsfile):
        self.filename = quartetsfile
        #self.quartet_dict = {}
        #with open(quartetsfile) as f:
            #for line in f:
                #(key,val) = line.split()
                #self.quartet_dict[key] = int(val)
    def quartet_dict(self):
        d = {}
        with open(self.filename) as f:
            for line in f:
                (key,val) = line.split()
                d[key] = int(val)
        return d

#L is a list of four labels that are strings-alphabetical order required
    def get_freqs(self,L):
        s = self.quartet_dict()
        #q = ['(('+L[0]+','+L[1]+'),('+L[2]+','+L[3]+'));','(('+L[0]+','+L[2]+'),('+L[1]+','+L[3]+'));','(('+L[0]+','+L[3]+'),('+L[1]+','+L[2]+'));']
        #u = [0,0,0]
        #for i in range(3):
        #    if q[i] in s:
        #        u[i] = s[q[i]]
        u = [0,0,0]
        q0 = ['(('+L[0]+','+L[1]+'),('+L[2]+','+L[3]+'));','(('+L[0]+','+L[1]+'),('+L[3]+','+L[2]+'));','(('+L[1]+','+L[0]+'),('+L[2]+','+L[3]+'));','(('+L[1]+','+L[0]+'),('+L[3]+','+L[2]+'));', '(('+L[2]+','+L[3]+'),('+L[0]+','+L[1]+'));','(('+L[2]+','+L[3]+'),('+L[1]+','+L[0]+'));','(('+L[3]+','+L[2]+'),('+L[0]+','+L[1]+'));','(('+L[3]+','+L[2]+'),('+L[1]+','+L[0]+'));']
        q1 = ['(('+L[0]+','+L[2]+'),('+L[1]+','+L[3]+'));','(('+L[0]+','+L[2]+'),('+L[3]+','+L[1]+'));','(('+L[2]+','+L[0]+'),('+L[1]+','+L[3]+'));','(('+L[2]+','+L[0]+'),('+L[3]+','+L[1]+'));','(('+L[1]+','+L[3]+'),('+L[0]+','+L[2]+'));','(('+L[1]+','+L[3]+'),('+L[2]+','+L[0]+'));','(('+L[3]+','+L[1]+'),('+L[0]+','+L[2]+'));','(('+L[3]+','+L[1]+'),('+L[2]+','+L[0]+'));']
        q2 = ['(('+L[0]+','+L[3]+'),('+L[1]+','+L[2]+'));', '(('+L[0]+','+L[3]+'),('+L[2]+','+L[1]+'));','(('+L[3]+','+L[0]+'),('+L[1]+','+L[2]+'));', '(('+L[3]+','+L[0]+'),('+L[2]+','+L[1]+'));','(('+L[1]+','+L[2]+'),('+L[0]+','+L[3]+'));','(('+L[1]+','+L[2]+'),('+L[3]+','+L[0]+'));','(('+L[2]+','+L[1]+'),('+L[0]+','+L[3]+'));','(('+L[2]+','+L[1]+'),('+L[3]+','+L[0]+'));']
        for i in range(8):
            if q0[i] in s:
                u[0] = s[q0[i]]
        for i in range(8):
            if q1[i] in s:
                u[1] = s[q1[i]]
        for i in range(8):
            if q2[i] in s:
                u[2] = s[q2[i]]
        return u
        
#L same as get_freqs, i,j are indices of L specified as cherry of interest
    def quartet_score(self,L,i,j):
        a = [i,j]
        a.sort()
        f1, f2, f3 = self.get_freqs(L)
        if a in [[0,1], [2,3]]:
            u1, u2, u3 = f1, f2, f3
        elif a in [[0,2], [1,3]]:
            u1, u2, u3 = f2, f1, f3
        elif a in [[0,3], [1,2]]:
            u1, u2, u3 = f3, f1, f2
        else:
            print 'problem with quartet score input indices'
            print L, i, j
        a12 = min(u1-u2,0)
        a13 = min(u1-u3,0)
        score = abs(u2-u3) + (-1)*a12 + (-1)*a13
        return score

    def quartet_labels_dict(self):
        g = self.quartet_dict()
        h = {}
        gkeys = g.keys()
        for i in range(len(gkeys)):
            keycopy = gkeys[i]
            labellist = keycopy.split(',')
            for j in range(4):
                labellist[j] = labellist[j].translate(None,string.punctuation)        
            h[gkeys[i]] = labellist
        return h

#D1,D2  are two taxa labels as strings in alphabetical order
    def score_double(self,D1,D2):
        g = self.quartet_dict()
        h = self.quartet_labels_dict()
        #print g
        #print h
        score = 0
        for k in g:
            if '(' + D1 + ',' + D2 + ')' in k:
                #print k
                newscore = self.quartet_score(h[k], h[k].index(D1), h[k].index(D2))
                #print newscore
                score = score + newscore
            if '(' + D2 + ',' + D1 + ')' in k:
                #print k
                newscore = self.quartet_score(h[k], h[k].index(D1), h[k].index(D2))
                #print newscore
                score = score + newscore
        return score

    #def penalty_score(self, subset1, subset2):
        
    def score_tree(self,treequartetsfile):
        treelist = []
        with open(treequartetsfile) as f:
            for line in f:
                (keygen,val) = line.split()
                labellist = keygen.split(',')
                for j in range(4):
                    labellist[j] = labellist[j].translate(None,string.punctuation)        
                treelist.append(labellist)
        #print treelist
        score = 0
        for k in range(len(treelist)):
                newscore = self.quartet_score(treelist[k],0,1)
                #print newscore
                score = score + newscore
        return score

class QuartetsInfoInequalitiesOnly:
    def __init__(self,quartetsfile):
        self.filename = quartetsfile
        #self.quartet_dict = {}
        #with open(quartetsfile) as f:
            #for line in f:
                #(key,val) = line.split()
                #self.quartet_dict[key] = int(val)
    def quartet_dict(self):
        d = {}
        with open(self.filename) as f:
            for line in f:
                (key,val) = line.split()
                d[key] = int(val)
        return d

#L is a list of four labels that are strings
    def get_freqs(self,L):
        s = self.quartet_dict()
        #q = ['(('+L[0]+','+L[1]+'),('+L[2]+','+L[3]+'));','(('+L[0]+','+L[2]+'),('+L[1]+','+L[3]+'));','(('+L[0]+','+L[3]+'),('+L[1]+','+L[2]+'));']
        #u = [0,0,0]
        #for i in range(3):
        #    if q[i] in s:
        #        u[i] = s[q[i]]
        u = [0,0,0]
        q0 = ['(('+L[0]+','+L[1]+'),('+L[2]+','+L[3]+'));','(('+L[0]+','+L[1]+'),('+L[3]+','+L[2]+'));','(('+L[1]+','+L[0]+'),('+L[2]+','+L[3]+'));','(('+L[1]+','+L[0]+'),('+L[3]+','+L[2]+'));', '(('+L[2]+','+L[3]+'),('+L[0]+','+L[1]+'));','(('+L[2]+','+L[3]+'),('+L[1]+','+L[0]+'));','(('+L[3]+','+L[2]+'),('+L[0]+','+L[1]+'));','(('+L[3]+','+L[2]+'),('+L[1]+','+L[0]+'));']
        q1 = ['(('+L[0]+','+L[2]+'),('+L[1]+','+L[3]+'));','(('+L[0]+','+L[2]+'),('+L[3]+','+L[1]+'));','(('+L[2]+','+L[0]+'),('+L[1]+','+L[3]+'));','(('+L[2]+','+L[0]+'),('+L[3]+','+L[1]+'));','(('+L[1]+','+L[3]+'),('+L[0]+','+L[2]+'));','(('+L[1]+','+L[3]+'),('+L[2]+','+L[0]+'));','(('+L[3]+','+L[1]+'),('+L[0]+','+L[2]+'));','(('+L[3]+','+L[1]+'),('+L[2]+','+L[0]+'));']
        q2 = ['(('+L[0]+','+L[3]+'),('+L[1]+','+L[2]+'));', '(('+L[0]+','+L[3]+'),('+L[2]+','+L[1]+'));','(('+L[3]+','+L[0]+'),('+L[1]+','+L[2]+'));', '(('+L[3]+','+L[0]+'),('+L[2]+','+L[1]+'));','(('+L[1]+','+L[2]+'),('+L[0]+','+L[3]+'));','(('+L[1]+','+L[2]+'),('+L[3]+','+L[0]+'));','(('+L[2]+','+L[1]+'),('+L[0]+','+L[3]+'));','(('+L[2]+','+L[1]+'),('+L[3]+','+L[0]+'));']
        for i in range(8):
            if q0[i] in s:
                u[0] = s[q0[i]]
        for i in range(8):
            if q1[i] in s:
                u[1] = s[q1[i]]
        for i in range(8):
            if q2[i] in s:
                u[2] = s[q2[i]]
        return u
        
        
#L same as get_freqs, i,j are indices of L specified as cherry of interest
    def quartet_score(self,L,i,j):
        a = [i,j]
        a.sort()
        f1, f2, f3 = self.get_freqs(L)
        if a in [[0,1], [2,3]]:
            u1, u2, u3 = f1, f2, f3
        elif a in [[0,2], [1,3]]:
            u1, u2, u3 = f2, f1, f3
        elif a in [[0,3], [1,2]]:
            u1, u2, u3 = f3, f1, f2
        else:
            print 'problem with quartet score input indices'
            print L, i, j
        a12 = min(u1-u2,0)
        a13 = min(u1-u3,0)
        score = (-1)*a12 + (-1)*a13
        return score

    def quartet_labels_dict(self):
        g = self.quartet_dict()
        h = {}
        gkeys = g.keys()
        for i in range(len(gkeys)):
            keycopy = gkeys[i]
            labellist = keycopy.split(',')
            for j in range(4):
                labellist[j] = labellist[j].translate(None,string.punctuation)        
            h[gkeys[i]] = labellist
        return h

#D1,D2  are two taxa labels as strings in alphabetical order
    def score_double(self,D1,D2):
        g = self.quartet_dict()
        h = self.quartet_labels_dict()
        #print g
        #print h
        score = 0
        for k in g:
            if '(' + D1 + ',' + D2 + ')' in k:
                #print k
                newscore = self.quartet_score(h[k], h[k].index(D1), h[k].index(D2))
                #print newscore
                score = score + newscore
        return score

    #def penalty_score(self, subset1, subset2):
        
    def score_tree(self,treequartetsfile):
        treelist = []
        with open(treequartetsfile) as f:
            for line in f:
                (keygen,val) = line.split()
                labellist = keygen.split(',')
                for j in range(4):
                    labellist[j] = labellist[j].translate(None,string.punctuation)        
                treelist.append(labellist)
        #print treelist
        score = 0
        for k in range(len(treelist)):
                newscore = self.quartet_score(treelist[k],0,1)
                #print newscore
                score = score + newscore
        return score


class SubsetPenalties:
    def __init__(self,labels,setlist,quartetsfile):
        self.quartetsinfo = QuartetsInfo(quartetsfile)
        self.setlist = setlist
        self.matrix_maker = matrixmaker.MatrixMaker(labels,setlist)
        self.matrix = self.matrix_maker.matrix()
        self.labels = self.matrix_maker.labels

    #def pairs(self):
        #clist = copy.copy(self.setlist)
        #clabels = copy.copy(self.labels)
        #p = []
        #while len(clist) > 0:
            #one = clist.pop(0)
            #two = list(set(clabels) - set(one))
            #if two in clist:
                #p.append([one, two])
                #clist.remove(two)
            #else:
                #print "missing pair from input setlist"
        #return p
        
        
    def add_quartets(self,set1,set2):
        stuntlabels = []
        for i in range(len(self.labels)):
            stuntlabels.append(self.labels[i])
        compl = set1 + set2
        #print stuntlabels, compl
        for s in compl:
            stuntlabels.remove(s)
            #print stuntlabels
        #for t in set2:
            #stuntlabels.remove(t)
            #print stuntlabels
        A = itertools.combinations(stuntlabels,2)
        AA = list(A)
        SminusApairs = [list(AA[j]) for j in range(len(AA))]
        A1A2 = []
        for k in range(len(set1)):
            for l in range(len(set2)):
                A1A2.append([set1[k],set2[l]])
        AQ = []
        for m in range(len(A1A2)):
            for n in range(len(SminusApairs)):
                tmp = A1A2[m] + SminusApairs[n]
                tmp.sort()
                AQ.append([tmp, tmp.index(A1A2[m][0]), tmp.index(A1A2[m][1])])
        #print AQ
        return AQ

    def subtract_quartets(self,set1,set2):
        a1 = itertools.combinations(set1,2)
        a2 = itertools.combinations(set2,2)
        aa1 = list(a1)
        aa2 = list(a2)
        A1 = [list(aa1[i]) for i in range(len(aa1))]
        A2 = [list(aa2[j]) for j in range(len(aa2))]
        SQ = []
        for m in range(len(A1)):
            for n in range(len(A2)):
                tmp = A1[m] + A2[n]
                tmp.sort()
                SQ.append([tmp,tmp.index(A1[m][0]), tmp.index(A1[m][1])])
        #print SQ
        return SQ

#need to enter score1, score2 for set1, set2 now: later can make this flexible? sets of size 1 and 2 are 0 and fixed. 
    def penalty_score(self,set1,set2,score1,score2):
        g = self.quartetsinfo.quartet_dict()
        h = self.quartetsinfo.quartet_labels_dict()
        #labs = self.labels
        #print labs, set1, set2
        AQ = self.add_quartets(set1,set2)
        #print AQ
        SQ = self.subtract_quartets(set1,set2)
        #print SQ
        score = score1 + score2
        for i in range(len(AQ)):
            score = score+self.quartetsinfo.quartet_score(AQ[i][0], AQ[i][1], AQ[i][2])
            #print score
        for j in range(len(SQ)): 
            score = score-self.quartetsinfo.quartet_score(SQ[j][0], SQ[j][1], SQ[j][2])
        return score

#problem is finding min NONZERO entry here. 
    #def scored_matrix(self):
        #SM = np.copy(self.matrix)
        #for j in range(len(self.setlist)):
            #SM[j,j] = np.inf
        #for i in range(len(self.setlist)):
            #if len(self.setlist[i]) == 1:
                #SM[i,i] = 0
            #if len(self.setlist[i]) == 2:
                #m = self.setlist.index([self.setlist[i][0]])
                #n = self.setlist.index([self.setlist[i][1]])
                #SM[i,m] = self.penalty_score([self.setlist[i][0]], [self.setlist[i][1]], 0, 0)
                #SM[i,n] = self.penalty_score([self.setlist[i][0]], [self.setlist[i][1]], 0, 0)
            #if len(self.setlist[i]) > 2:
                #clabels = copy.copy(self.setlist[i])
                #clist = copy.copy(self.setlist)
                #for s in clist:
                    #if SM[i,self.setlist.index(s)] == np.inf:
                        #clist.remove(s)
                #while len(clist) > 0:
                        #one = clist.pop(0)
                        #two = list(set(clabels) - set(one))
                        #one.sort()
                        #two.sort()
                        #if two in clist:
                            #clist.remove(two)
                            #onedex = self.setlist.index(one)
                            #twodex = self.setlist.index(two)
                            #SM[i,onedex] = self.penalty_score(one, two, min(SM[onedex]), min(SM[twodex]))
                            #SM[i,twodex] = self.penalty_score(one, two, min(SM[onedex]), min(SM[twodex]))
                        #else: 
                            #SM[i,onedex] = np.inf
                            #SM[i,twodex] = np.inf
        #print SM
        #return SM

    def scored_matrix(self):
        SM = np.zeros((len(self.setlist),len(self.setlist)))
        SM.fill(np.inf)
        for i in range(len(self.setlist)):
            if len(self.setlist[i]) == 1:
                SM[i,i] = 0
            if len(self.setlist[i]) == 2:
                m = self.setlist.index([self.setlist[i][0]])
                n = self.setlist.index([self.setlist[i][1]])
                SM[i,m] = self.penalty_score([self.setlist[i][0]], [self.setlist[i][1]], 0, 0)
                SM[i,n] = self.penalty_score([self.setlist[i][0]], [self.setlist[i][1]], 0, 0)
            if len(self.setlist[i]) > 2:
                clabels = copy.copy(self.setlist[i])
                clist = copy.copy(self.setlist)
                while len(clist) > 0:
                        one = clist.pop(0)
                        onedex = self.setlist.index(one)
                        if self.matrix[i,onedex] == 1:
                            two = list(set(clabels) - set(one))
                            one.sort()
                            two.sort()
                            if two in clist:
                                clist.remove(two)
                                twodex = self.setlist.index(two)
                                SM[i,onedex] = self.penalty_score(one, two, min(SM[onedex]), min(SM[twodex]))
                                SM[i,twodex] = self.penalty_score(one, two, min(SM[onedex]), min(SM[twodex]))
                        #else: 
        return SM

    #def clades(self):
        #SM = self.scored_matrix()
        #L = []
        #temprows =[self.setlist.index(self.setlist[-1]), self.setlist.index(self.setlist[-1])]
        #while temprow > 0:
                #onedex = SM[temprow].argmin()
                #one = self.setlist[onedex]
                #two = list(set(self.setlist[temprow])-set(one))
                #two.sort()
                #twodex = self.setlist.index(two)



    def clades(self):
        SM = self.scored_matrix()
        #L = []
        L = [self.setlist[-1]]
        tempsets = [self.setlist[-1]]
        while len(tempsets) > 0:
            tmp = tempsets.pop(0)
            if len(tmp) == 1:
                L.append(tmp)
            else:
                rowdex = self.setlist.index(tmp)
                onedex = SM[rowdex].argmin()
                one = self.setlist[onedex]
                two = list(set(self.setlist[rowdex]) - set(one))
                two.sort()
                L.append(one)
                L.append(two)
                if len(one) > 1:
                    tempsets.append(one)
                if len(two) > 1:
                    tempsets.append(two)
        return L

    #def tree(self):
        #L = self.clades()
        #start = L.pop(0)
        #s = ''
        #for i in range(len(start)-1):
            #s = s + start[i]+ ','
        #s = s + start[-1]
        #s = '(' + s + ')'
        #while len(L) > 0:
            #tmp = L.pop(0)
            #if len(tmp)> 1:
                #spot = s.index(tmp[0])
                #t = '('
                #for j in range(len(tmp)-1):
                    #t = t + tmp[j] +','
                #t = t + tmp[-1] + ')'
                #for k in range(len(tmp)):
                    #s.replace(tmp[k] + ',','')
                    #s.replace(tmp[k],'')
                #s = s[:spot] + t + s[spot:]
        #return s

    def tree(self):
        L = self.clades()
        start = L.pop(0)
        s = ''
        for i in range(len(start)):
            s = s + start[i]+ ','
        #s = s + start[-1]
        s = '(' + s + ')'
        while len(L) > 0:
            tmp = L.pop(0)
            if len(tmp) > 1:
                first = tmp[0]
                #find the first taxon, to use first letter as string index place
                spot = s.index(first[0])
                t = '('
                for j in range(len(tmp)):
                    t = t + tmp[j] + ','
                t = t + '),'
                for k in range(len(tmp)):
                    s = s.replace(tmp[k] + ',', '')
                    s = s.replace(tmp[k], '')
                s = s[:spot] + t + s[spot:] 
        s = s.replace(',)',')')
        s = s + ';'
        return s

class SubsetPenaltiesInequalitiesOnly:
    def __init__(self,labels,setlist,quartetsfile):
        self.quartetsinfo = QuartetsInfoInequalitiesOnly(quartetsfile)
        self.setlist = setlist
        self.matrix_maker = matrixmaker.MatrixMaker(labels,setlist)
        self.matrix = self.matrix_maker.matrix()
        self.labels = self.matrix_maker.labels

    #def pairs(self):
        #clist = copy.copy(self.setlist)
        #clabels = copy.copy(self.labels)
        #p = []
        #while len(clist) > 0:
            #one = clist.pop(0)
            #two = list(set(clabels) - set(one))
            #if two in clist:
                #p.append([one, two])
                #clist.remove(two)
            #else:
                #print "missing pair from input setlist"
        #return p
        
        
    def add_quartets(self,set1,set2):
        stuntlabels = []
        for i in range(len(self.labels)):
            stuntlabels.append(self.labels[i])
        compl = set1 + set2
        #print stuntlabels, compl
        for s in compl:
            stuntlabels.remove(s)
            #print stuntlabels
        #for t in set2:
            #stuntlabels.remove(t)
            #print stuntlabels
        A = itertools.combinations(stuntlabels,2)
        AA = list(A)
        SminusApairs = [list(AA[j]) for j in range(len(AA))]
        A1A2 = []
        for k in range(len(set1)):
            for l in range(len(set2)):
                A1A2.append([set1[k],set2[l]])
        AQ = []
        for m in range(len(A1A2)):
            for n in range(len(SminusApairs)):
                tmp = A1A2[m] + SminusApairs[n]
                tmp.sort()
                AQ.append([tmp, tmp.index(A1A2[m][0]), tmp.index(A1A2[m][1])])
        #print AQ
        return AQ

    def subtract_quartets(self,set1,set2):
        a1 = itertools.combinations(set1,2)
        a2 = itertools.combinations(set2,2)
        aa1 = list(a1)
        aa2 = list(a2)
        A1 = [list(aa1[i]) for i in range(len(aa1))]
        A2 = [list(aa2[j]) for j in range(len(aa2))]
        SQ = []
        for m in range(len(A1)):
            for n in range(len(A2)):
                tmp = A1[m] + A2[n]
                tmp.sort()
                SQ.append([tmp,tmp.index(A1[m][0]), tmp.index(A1[m][1])])
        #print SQ
        return SQ

#need to enter score1, score2 for set1, set2 now: later can make this flexible? sets of size 1 and 2 are 0 and fixed. 
    def penalty_score(self,set1,set2,score1,score2):
        g = self.quartetsinfo.quartet_dict()
        h = self.quartetsinfo.quartet_labels_dict()
        #labs = self.labels
        #print labs, set1, set2
        AQ = self.add_quartets(set1,set2)
        #print AQ
        SQ = self.subtract_quartets(set1,set2)
        #print SQ
        score = score1 + score2
        for i in range(len(AQ)):
            score = score+self.quartetsinfo.quartet_score(AQ[i][0], AQ[i][1], AQ[i][2])
            #print score
        for j in range(len(SQ)): 
            score = score-self.quartetsinfo.quartet_score(SQ[j][0], SQ[j][1], SQ[j][2])
        return score

#problem is finding min NONZERO entry here. 
    #def scored_matrix(self):
        #SM = np.copy(self.matrix)
        #for j in range(len(self.setlist)):
            #SM[j,j] = np.inf
        #for i in range(len(self.setlist)):
            #if len(self.setlist[i]) == 1:
                #SM[i,i] = 0
            #if len(self.setlist[i]) == 2:
                #m = self.setlist.index([self.setlist[i][0]])
                #n = self.setlist.index([self.setlist[i][1]])
                #SM[i,m] = self.penalty_score([self.setlist[i][0]], [self.setlist[i][1]], 0, 0)
                #SM[i,n] = self.penalty_score([self.setlist[i][0]], [self.setlist[i][1]], 0, 0)
            #if len(self.setlist[i]) > 2:
                #clabels = copy.copy(self.setlist[i])
                #clist = copy.copy(self.setlist)
                #for s in clist:
                    #if SM[i,self.setlist.index(s)] == np.inf:
                        #clist.remove(s)
                #while len(clist) > 0:
                        #one = clist.pop(0)
                        #two = list(set(clabels) - set(one))
                        #one.sort()
                        #two.sort()
                        #if two in clist:
                            #clist.remove(two)
                            #onedex = self.setlist.index(one)
                            #twodex = self.setlist.index(two)
                            #SM[i,onedex] = self.penalty_score(one, two, min(SM[onedex]), min(SM[twodex]))
                            #SM[i,twodex] = self.penalty_score(one, two, min(SM[onedex]), min(SM[twodex]))
                        #else: 
                            #SM[i,onedex] = np.inf
                            #SM[i,twodex] = np.inf
        #print SM
        #return SM

    def scored_matrix(self):
        SM = np.zeros((len(self.setlist),len(self.setlist)))
        SM.fill(np.inf)
        for i in range(len(self.setlist)):
            if len(self.setlist[i]) == 1:
                SM[i,i] = 0
            if len(self.setlist[i]) == 2:
                m = self.setlist.index([self.setlist[i][0]])
                n = self.setlist.index([self.setlist[i][1]])
                SM[i,m] = self.penalty_score([self.setlist[i][0]], [self.setlist[i][1]], 0, 0)
                SM[i,n] = self.penalty_score([self.setlist[i][0]], [self.setlist[i][1]], 0, 0)
            if len(self.setlist[i]) > 2:
                clabels = copy.copy(self.setlist[i])
                clist = copy.copy(self.setlist)
                while len(clist) > 0:
                        one = clist.pop(0)
                        onedex = self.setlist.index(one)
                        if self.matrix[i,onedex] == 1:
                            two = list(set(clabels) - set(one))
                            one.sort()
                            two.sort()
                            if two in clist:
                                clist.remove(two)
                                twodex = self.setlist.index(two)
                                SM[i,onedex] = self.penalty_score(one, two, min(SM[onedex]), min(SM[twodex]))
                                SM[i,twodex] = self.penalty_score(one, two, min(SM[onedex]), min(SM[twodex]))
                        #else: 
        return SM

    #def clades(self):
        #SM = self.scored_matrix()
        #L = []
        #temprows =[self.setlist.index(self.setlist[-1]), self.setlist.index(self.setlist[-1])]
        #while temprow > 0:
                #onedex = SM[temprow].argmin()
                #one = self.setlist[onedex]
                #two = list(set(self.setlist[temprow])-set(one))
                #two.sort()
                #twodex = self.setlist.index(two)



    def clades(self):
        SM = self.scored_matrix()
        #L = []
        L = [self.setlist[-1]]
        tempsets = [self.setlist[-1]]
        while len(tempsets) > 0:
            tmp = tempsets.pop(0)
            if len(tmp) == 1:
                L.append(tmp)
            else:
                rowdex = self.setlist.index(tmp)
                onedex = SM[rowdex].argmin()
                one = self.setlist[onedex]
                two = list(set(self.setlist[rowdex]) - set(one))
                two.sort()
                L.append(one)
                L.append(two)
                if len(one) > 1:
                    tempsets.append(one)
                if len(two) > 1:
                    tempsets.append(two)
        return L

    #def tree(self):
        #L = self.clades()
        #start = L.pop(0)
        #s = ''
        #for i in range(len(start)-1):
            #s = s + start[i]+ ','
        #s = s + start[-1]
        #s = '(' + s + ')'
        #while len(L) > 0:
            #tmp = L.pop(0)
            #if len(tmp)> 1:
                #spot = s.index(tmp[0])
                #t = '('
                #for j in range(len(tmp)-1):
                    #t = t + tmp[j] +','
                #t = t + tmp[-1] + ')'
                #for k in range(len(tmp)):
                    #s.replace(tmp[k] + ',','')
                    #s.replace(tmp[k],'')
                #s = s[:spot] + t + s[spot:]
        #return s

    def tree(self):
        L = self.clades()
        start = L.pop(0)
        s = ''
        for i in range(len(start)):
            s = s + start[i]+ ','
        #s = s + start[-1]
        s = '(' + s + ')'
        while len(L) > 0:
            tmp = L.pop(0)
            if len(tmp) > 1:
                first = tmp[0]
                #find the first taxon, to use first letter as string index place
                spot = s.index(first[0])
                t = '('
                for j in range(len(tmp)):
                    t = t + tmp[j] + ','
                t = t + '),'
                for k in range(len(tmp)):
                    s = s.replace(tmp[k] + ',', '')
                    s = s.replace(tmp[k], '')
                s = s[:spot] + t + s[spot:] 
        s = s.replace(',)',')')
        s = s + ';'
        return s

#dendropy_clades_tree takes list of Taxa and list of clades, returns a tree
def dendropy_clades_tree(Taxa,cladelist):
    trees = ''
    for i in range(len(cladelist)):
        #print i
        compl = list(set(Taxa)-set(cladelist[i]))
        #print compl
        if len(compl) > 0:
            half1 = '(('
            for k in range(len(cladelist[i])-1):
                half1 = half1 + cladelist[i][k] + ','
            half1 = half1 + cladelist[i][-1] +'),'
            half2 = '('
            for j in range(len(compl)-1):
                half2 = half2 + compl[j] + ','
            half2  = half2 + compl[-1] + '));'
            tree = half1+half2
            trees = trees+tree
    tmplist = dendropy.TreeList.get_from_string(trees,'newick')
    finaltree = tmplist.consensus(min_freq=0)
    return finaltree

        


#powerset makes the list of all subsets of a label set
def powerset(labellist):
    P = []
    for j in range(1,len(labellist)+1):
        tmplist = list(itertools.combinations(labellist,j))
        tmplist = [list(tmplist[i]) for i in range(len(tmplist))]
        for k in range(len(tmplist)):
            tmplist[k].sort()
        P = P + tmplist
    return P

class QuintetsMaker:
    def __init__(self,labels):
        self.labels = labels

#here l is list of five labels
    def quintets(self,l):
        t1 = '((' +l[0]+ ',' +l[1]+'),' +l[2]+ ',(' +l[3]+ ',' +l[4]+ '));\n'
        t2 = '((' +l[0]+ ',' +l[1]+'),' +l[3]+ ',(' +l[2]+ ',' +l[4]+ '));\n'
        t3 = '((' +l[0]+ ',' +l[1]+'),' +l[4]+ ',(' +l[2]+ ',' +l[3]+ '));\n'
        t4 = '((' +l[0]+ ',' +l[2]+'),' +l[1]+ ',(' +l[3]+ ',' +l[4]+ '));\n'
        t5 = '((' +l[0]+ ',' +l[2]+'),' +l[3]+ ',(' +l[1]+ ',' +l[4]+ '));\n'
        t6 = '((' +l[0]+ ',' +l[2]+'),' +l[4]+ ',(' +l[1]+ ',' +l[3]+ '));\n'
        t7 = '((' +l[0]+ ',' +l[3]+'),' +l[1]+ ',(' +l[2]+ ',' +l[4]+ '));\n'
        t8 = '((' +l[0]+ ',' +l[3]+'),' +l[2]+ ',(' +l[1]+ ',' +l[4]+ '));\n'
        t9 = '((' +l[0]+ ',' +l[3]+'),' +l[4]+ ',(' +l[1]+ ',' +l[2]+ '));\n'
        t10 = '((' +l[0]+ ',' +l[4]+'),' +l[1]+ ',(' +l[2]+ ',' +l[3]+ '));\n'
        t11 = '((' +l[0]+ ',' +l[4]+'),' +l[2]+ ',(' +l[1]+ ',' +l[3]+ '));\n'
        t12 = '((' +l[0]+ ',' +l[4]+'),' +l[3]+ ',(' +l[1]+ ',' +l[2]+ '));\n'
        t13 = '((' +l[1]+ ',' +l[2]+'),' +l[0]+ ',(' +l[3]+ ',' +l[4]+ '));\n'
        t14 = '((' +l[1]+ ',' +l[3]+'),' +l[0]+ ',(' +l[2]+ ',' +l[4]+ '));\n'
        t15 = '((' +l[1]+ ',' +l[4]+'),' +l[0]+ ',(' +l[2]+ ',' +l[3]+ '));\n'
        q = [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15]
        return q
    
    def fives(self):
        s = list(itertools.combinations(self.labels,5))
        return s
    
    def allquintets(self):
        f = self.fives()
        All = []
        for j in range(len(f)):
            All = All + self.quintets(f[j])
        return All
#give filename as a string
    def makefile(self,filename):
        f = open(filename,'w')
        for k in range(len(self.allquintets())):
            f.write(self.allquintets()[k])
        f.close()

class FileClades:
    def __init__(self,filename):
        self.filename = filename

    def clades(self):
        C = []
        f = open(self.filename, 'r')
        for line in f:
            if '{' in line:
                C.append(line)
        f.close()
        for i in range(len(C)):
            C[i] = C[i].replace('{','')
        for j in range(len(C)):
            C[j] = C[j].replace('}','')
        for k in range(len(C)):
            C[k] = C[k].replace('\n','')
        for n in range(len(C)):
            C[n] = C[n].replace(' ','')
        D = []
        for a in range(len(C)):
            if len(C[a]) > 0:
                c = C[a].split(',')
                c.sort()
                D.append(c)
        D.reverse()
        Taxa = []
        for i in range(len(D)):
            if len(D[i]) == 1:
                Taxa = Taxa + D[i]
        Taxa.sort()
        #Taxa.sort(key=float)
        if (Taxa in D) == False:
            D.append(Taxa)
        return D



#labels here is a list of just five labels
#def quintets(labels):
    

#this build tree takes an output from a clades function from SubsetPenalties.clades
    

#This class was not working well, replace with score_matrix function in class SubsetPenalties 
#matrix input called m  here is an instance of matrixmaker.MatrixMaker(labels,setlist)
#class ScoredMatrix:
    #def __init__(self,labels,setlist,quartetsfile):
        #self.start_matrix = matrixmaker.MatrixMaker(labels,setlist)
        #self.quartetsfile = quartetsfile
        #self.labels = self.start_matrix.labels
        #self.setlist = self.start_matrix.setlist
        #self.quartetsinfo = QuartetsInfo(quartetsfile)
        #self.matrixy = self.start_matrix.matrix()
        #self.quartet_dict = {}

    #def scored_matrix(self):
        #K  = self.matrixy
        #L = self.setlist
        #for i in range(len(L)): 
            #if len(L[i]) == 2:
               #for j in range(i-1):
                    #if K[i, j] == 1:
                        #K[i,j] = self.quartetsinfo.score_double(L[i][0], L[i][1])
        #return K            
        
