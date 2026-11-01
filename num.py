# elementary number theory

class ModInt:
    
    def __init__(self,n, modulus):
        self.modulus = modulus
        self.rep=n%modulus

    def int(n, modulus):
        return ModInt(n, modulus)

    def inverse(self):
        if gcd(self.rep, self.modulus) != 1:
            return ValueError("Division by non-unit")
        for i in range(self.modulus):
            if i*self.rep % self.modulus == 1:
                return ModInt(i, self.modulus)

    def is_zero(self):
        return self.rep == 0

    def is_one(self):
        return self.rep == 1

    def __add__(self, other):
        return ModInt(self.rep+other.rep, self.modulus)
        
    def __mul__(self,other):
        return ModInt(self.rep * other.rep, self.modulus)
    
    def __sub__(self,other):
        return ModInt(self.rep - other.rep, self.modulus)
        
    def __truediv__(self,other):
        return self * other.inverse()
        
    def __eq__(self,other):
        return self.rep==other.rep

    def __neg__(self):
        return ModInt(-self.rep, self.modulus)
        
    def __hash__(self):
        return hash((self.rep, self.modulus))

    def __str__(self):
        return str(self.rep)

def gcd(n, m):
    if m == 0:
        return n
    r = n % m
    return gcd(m, r)