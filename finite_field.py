from poly import *
from num import *

class FiniteField:

    def __init__(self, char, prim_poly, val=None):
        self.char = char
        self.prim_poly = prim_poly
        self.val = val
        if self.val is None:
            self.val = Polynomial(ModInt, char) #initialise at 0

    def primitive_element(char, prim_poly):
        return FiniteField(char, prim_poly, Polynomial(ModInt, char, coeff=[ModInt(0, char), ModInt(1, char)]))

    def from_p_ary(n, char, prim_poly):
        coeff=[]
        while n > 0:
            coeff.append(ModInt(n, char))
            n -= (n % char)
            n = n//char
        return FiniteField(char, prim_poly, val = Polynomial(ModInt, char, coeff=coeff))

    def same_int(self, n): # creates element n in the same finite field as self
        return FiniteField.int(n, self.char, self.prim_poly)

    def int(n, char, prim_poly):
        return FiniteField(char, prim_poly, Polynomial.int(n, ModInt, char))

    def is_zero(self):
        return (self == FiniteField.int(0, self.char, self.prim_poly))

    def is_one(self):
        return (self == FiniteField.int(1, self.char, self.prim_poly))

    def inverse(self):
        return FiniteField(self.char, self.prim_poly, val=Polynomial.ext_gcd(self.val, self.prim_poly)[0] % self.prim_poly)

    def __hash__(self):
        return hash((self.char, self.prim_poly, self.val))

    def __add__(self, other):
        return FiniteField(self.char, self.prim_poly, (self.val+other.val) % self.prim_poly)

    def __sub__(self, other):
        return FiniteField(self.char, self.prim_poly, (self.val-other.val) % self.prim_poly)

    def __mul__(self, other):
        return FiniteField(self.char, self.prim_poly, (self.val * other.val) % self.prim_poly)

    def __truediv__(self, other):
        return self * other.inverse()

    def __str__(self):
        s=""
        val=self.val
        for i in range(val.deg+1):
            if not val.coeff[i].is_zero():
                s = s+(str(val.coeff[i]) if not val.coeff[i].is_one() else "")+"\u03B1^"+str(i)
                if i<val.deg:
                    s+="+"
        return s

    def __eq__(self, other):
        return (self.val == other.val)

    def __neg__(self):
        return FiniteField(self.char, self.prim_poly, -self.val)

