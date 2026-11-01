class MultiSet:
    
    def __init__(self, multiset={}):
        self.multiset=multiset
    
    def __add__(self,other):
        ms=self.multiset.copy()
        for key in other.multiset:
            if key in ms:
                ms[key] += other.multiset[key]
            else:
                ms[key] = other.multiset[key]
        return MultiSet(ms)
        
    def __mul__(self, n):
        ms=self.multiset.copy()
        for key in ms:
            ms[key] *= n
        return MultiSet(ms)
    
    def __str__(self):
        return str(self.multiset)  