class Polynomial:
    
    def __init__(self, coeff_ring, *coeff_ring_init_args, coeff = None):
        self.coeff_ring = coeff_ring
        self.coeff_ring_init_args = coeff_ring_init_args
        self.coeff=coeff
        if coeff is None:
            self.coeff = [self.coeff_int(0)]
        self.deg=None
        self.lc=None
        self.update() # sets the degree, leading coefficient, and truncates leading zeros

    def coeff_int(self, n): # creates the integer n in the same coeff ring as self
        return self.coeff_ring.int(n, *self.coeff_ring_init_args)

    def poly_int(self, n): #creates the constant polynomial n over the same coeff ring as self
        return Polynomial.int(n, self.coeff_ring, *self.coeff_ring_init_args)

    def int(n, coeff_ring, *coeff_ring_init_args): # creates the integer n as a polynomial
        return Polynomial(coeff_ring, *coeff_ring_init_args, coeff=[coeff_ring.int(n, *coeff_ring_init_args)])

    def is_zero(self):
        return (self == self.poly_int(0))

    def is_one(self):
        return (self == self.poly_int(1))

    def __add__(self,other):
        new_coeff=[]
        for i in range(max(self.deg, other.deg)+1):
            if i>self.deg:
                new_coeff.append(other.coeff[i])
            elif i>other.deg:
                new_coeff.append(self.coeff[i])
            else:
                new_coeff.append(self.coeff[i]+other.coeff[i])
        return Polynomial(self.coeff_ring, *self.coeff_ring_init_args, coeff=new_coeff)
        
    def __sub__(self,other):
        new_coeff=[]
        for i in range(max(self.deg, other.deg)+1):
            if i>self.deg:
                new_coeff.append(-other.coeff[i])
            elif i>other.deg:
                new_coeff.append(self.coeff[i])
            else:
                new_coeff.append(self.coeff[i]-other.coeff[i])
        return Polynomial(self.coeff_ring, *self.coeff_ring_init_args, coeff=new_coeff)
    
    def __mul__(self,other):
        new_coeff=[]
        for i in range(self.deg+other.deg+1):
            new_coeff.append(self.coeff_int(0))
            for j in range(i+1):
                if j<len(self.coeff) and i-j<len(other.coeff):
                    new_coeff[i] = new_coeff[i] + (self.coeff[j]*other.coeff[i-j])
        return Polynomial(self.coeff_ring, *self.coeff_ring_init_args, coeff=new_coeff)

    def __floordiv__(self, other):
        return Polynomial.long_div(self, other)[0]
        
    def __eq__(self, other):
        if self.deg != other.deg:
            return False
        for i in range(self.deg+1):
            if self.coeff[i] != other.coeff[i]:
                return False
        return True
        
    def __hash__(self):
        return hash(tuple(self.coeff))

    def __mod__(self, other):
        return Polynomial.long_div(self, other)[1]

    def __neg__(self):
        return self.poly_int(0)-self

    def __str__(self):
        s=""
        for i in range(self.deg+1):
            if not self.coeff[i].is_zero():
                s = s+("("+str(self.coeff[i])+")" if not self.coeff[i].is_one() else "")+"X^"+str(i)
                if i<self.deg:
                    s+="+"
        return s
        
    def set_value(self, pos, val):
        if pos > self.deg:
            self.coeff += [self.coeff_int(0) for _ in range(pos-self.deg-1)]
            self.coeff.append(val)
        else:
            self.coeff[pos]=val
        self.update()
    
    def diff(self): # formal differentiation
        new_coeff=[]
        for i in range(1,self.deg+1):
            new_coeff.append(self.coeff_int(i)*self.coeff[i])
        return Polynomial(self.coeff_ring, *self.coeff_ring_init_args, coeff=new_coeff)
        
    def norm(self): # returns polynomial divided by leading coefficient
        new_coeff=[]
        c=self.coeff_int(1)/self.coeff[self.deg]
        for i in range(self.deg+1):
            new_coeff.append(c*self.coeff[i])
        return Polynomial(self.coeff_ring, *self.coeff_ring_init_args, coeff=new_coeff)

    def gcd(f,g):
        if g.is_zero():
            return f
        ld=Polynomial.long_div(f,g)
        return Polynomial.gcd(g,ld[1]).norm()

    def ext_gcd(f,g): # returns (a, b, d) such that d=gcd(f,g) and af+bg=d
        if g.is_zero(): # returns (1,0,f)
            return (Polynomial(f.coeff_ring, *f.coeff_ring_init_args, coeff=[f.lc.inverse()]), f.poly_int(0), f.norm())
        q,r=Polynomial.long_div(f,g) # f=qg+r
        a,b,d=Polynomial.ext_gcd(g,r) # ag+br=d
        # bf+(a-bq)g=d
        return (b, a-b*q, d)

    def long_div(f, g, q=None):
        if q is None:
            q = f.poly_int(0)
        if f.deg<g.deg or f.is_zero():
            return (q, f) # (q, f) = (quotient, remainder)
        h=f.poly_int(0)
        h.set_value(f.deg-g.deg, f.lc/g.lc)
        return Polynomial.long_div(f-h*g, g, q+h)

    def update(self):
        self.deg = 0
        self.lc = self.coeff_int(0)
        for i in range(len(self.coeff)-1, -1, -1):
            if self.coeff[i] != self.coeff_int(0):
                self.deg = i
                self.lc = self.coeff[i]
                break
            elif i != 0:
                self.coeff.pop()
