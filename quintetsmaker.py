import itertools

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


