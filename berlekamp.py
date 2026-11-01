from linalg import *
from poly import *
from num import *
from multiset import *
from finite_field import *
#import sys
#sys.setrecursionlimit(10) for debugging

p=5 #char of field
prim_poly=Polynomial(ModInt, p, coeff=[ModInt(0,p), ModInt(1,p)])
m=prim_poly.deg #deg of field over F_p
q=p**m #size of field

def cyclotomic(n): #recursively calculates cyclotomic polynomial
    product = Polynomial.int(1, FiniteField, p, prim_poly)
    for d in range(1, n):
        if n % d == 0:
            product *= cyclotomic(d)
    xn_minus_one = Polynomial(FiniteField, p, prim_poly)
    xn_minus_one.set_value(0, FiniteField.int(-1, p, prim_poly))
    xn_minus_one.set_value(n, FiniteField.int(1, p, prim_poly))
    return xn_minus_one//product

# computes the matrix representation of the Frobenius r |-> r^q on the ring R=F[X]/(f), where F=F_q
def frobenius_matrix(f):
    rows=[]
    d=f.deg
    for i in range(d):
        xip=Polynomial(FiniteField, p, prim_poly)
        xip.set_value(i*q,FiniteField.int(1, p, prim_poly)) # the polynomial X^{iq}
        r=xip % f
        rows.append(r.coeff + [FiniteField.int(0, p, prim_poly) for _ in range(d-r.deg-1)]) # need to add (d-r.deg-1) zeros to get a vector of length d
    return Matrix(FiniteField, p, prim_poly, rows=rows).transpose()

def no_of_factors(f):
    if f.deg == 0:
        return 0
    h=Polynomial.gcd(f,f.diff())
    if h.deg == f.deg:
        new_coeff=[]
        for i in range(f.deg//p+1):
            new_coeff.append(f.coeff[p*i])
        return p*no_of_factors(Polynomial(FiniteField, p, prim_poly, coeff=new_coeff))
    if h.deg > 0:
        return no_of_factors(h)+no_of_factors(Polynomial.long_div(f,h)[0])
    return f.deg - (frobenius_matrix(f)-Matrix.id_mat(f.deg, FiniteField, p, prim_poly)).rank() #nullity of F-I, where F is the matrix of the Frobenius

# returns multiset of irreducible factors
def berlekamp(f):
    if f.deg == 0:
        return MultiSet()
    f_d = f.diff()
    if f_d.is_zero():
        new_coeff=[]
        for i in range(f.deg//p+1):
            new_coeff.append(f.coeff[p*i])
        return berlekamp(Polynomial(FiniteField, p, prim_poly, coeff=new_coeff))*p
    h=Polynomial.gcd(f,f_d)
    if h.deg > 0:
        return berlekamp(h.norm()) + berlekamp((f//h).norm())
    mat=frobenius_matrix(f)-Matrix.id_mat(f.deg, FiniteField, p, prim_poly)
    mat.gauss()
    fact_poly=Matrix.null_space_basis(mat, f.deg)
    if len(fact_poly)==1: # f is irreducible
        return MultiSet({f:1})
    nontriv=False #found a non-trivial factorization of f
    i=0
    S=[]
    while i<len(fact_poly) and not nontriv:
        h_i=Polynomial(FiniteField, p, prim_poly, coeff=fact_poly[i])
        coeff_vec=[ModInt(0, p) for _ in range(m)]
        done = False
        while not done: # loop through elements of the finite field lexicographically
            a = FiniteField(p, prim_poly, Polynomial(ModInt, p, coeff=coeff_vec.copy())) # the element a = f(alpha) where alpha is a primitive element, and the coefficients of f are given by coeff_vec
            a_poly=Polynomial(FiniteField, p, prim_poly, coeff=[a]) # the constant polynomial a
            g_ia=Polynomial.gcd(f,h_i-a_poly).norm()
            if g_ia.deg>0 and g_ia.deg<f.deg:
                nontriv=True # found a non-trivial factorisation!
                S.append(g_ia)
            done = True
            # increment lexicographically
            for i in range(m-1, -1, -1):
                if coeff_vec[i] != ModInt(-1, p):
                    coeff_vec[i] += ModInt(1, p)
                    done = False
                    break
                else:
                    coeff_vec[i]=ModInt(0,p)
        i+=1
    factors=MultiSet()
    for g in S:
        factors += berlekamp(g)
    return factors


f=Polynomial(FiniteField, p, prim_poly)
f.set_value(0, FiniteField.int(1, p, prim_poly))
f.set_value(4, FiniteField.int(-10, p, prim_poly))
f.set_value(8, FiniteField.int(1, p, prim_poly))
print("POLYNOMIAL:", f)
factors=berlekamp(f)
print("NUMBER OF FACTORS:", len(factors.multiset.items()))
for key,value in factors.multiset.items():
   print("FACTOR: ", key, "  MULTIPLICITY: ", value)