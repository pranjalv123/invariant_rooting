import numpy as np


class MatrixMaker:
    def __init__(self, labels, setlist):
        self.labels = labels
        self.setlist = setlist
        self.binlabels = {}
        for i in range(len(labels)):
            self.binlabels[labels[i]] = 2**i

    def binset(self, setlistentry):
        b = 0
        for j in range(len(setlistentry)):
            b = b + self.binlabels[setlistentry[j]]
        return b
       
    def inclusion(self, set1, set2):
        mask1 = self.binset(set1)
        mask2 = self.binset(set2)
        if mask1 == (mask1 & mask2):
            return 1
        else:
            return np.inf

    def matrix(self):
        #M = np.zeros((len(self.setlist), len(self.setlist)))
        K = []
        for i in range(len(self.setlist)):
            row = [self.inclusion(self.setlist[j], self.setlist[i]) for j in range(len(self.setlist))]
            K.append(row)
        M = np.array(K)
        return M
   
        
