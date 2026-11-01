from functools import cmp_to_key

# MATRIX ARITHMETIC
class Matrix:
    
    def __init__(self, coeff_ring, *coeff_ring_args, rows=[]):
        self.coeff_ring=coeff_ring
        self.coeff_ring_args=coeff_ring_args
        self.rows=rows
    
    def __add__(self, other):
        if self.calc_dim() != other.calc_dim():
            raise ValueError("Cannot add/subtract matrices of different dimensions")
        rows=[]
        dim=self.calc_dim()
        for i in range(dim[0]):
            rows.append([])
            for j in range(dim[1]):
                rows[i].append(self.rows[i][j]+other.rows[i][j])
        return Matrix(self.coeff_ring, *self.coeff_ring_args, rows=rows)
    
    def __sub__(self, other):
        if self.calc_dim() != other.calc_dim():
            raise ValueError("Cannot add/subtract matrices of different dimensions")
        rows=[]
        dim=self.calc_dim()
        for i in range(dim[0]):
            rows.append([])
            for j in range(dim[1]):
                rows[i].append(self.rows[i][j]-other.rows[i][j])
        return Matrix(self.coeff_ring, *self.coeff_ring_args, rows=rows)
        
    def __mul__(self, other):
        r1,c1=self.calc_dim()
        r2,c2=other.calc_dim()
        if c1 != r2:
            raise ValueError("Wrong dimensions for multiplication!")
        rows=[]
        for i in range(r1):
            rows.append([]) # add the i-th row
            for j in range(c2):
                rows[i].append(self.coeff_ring.int(0, *self.coeff_ring_args)) # initialise (i,j)-th entry to be 0
                for k in range(c1):
                    rows[i][j] += self.rows[i][k] * other.rows[k][j]
        return Matrix(self.coeff_ring, *self.coeff_ring_args, rows=rows)

    def copy(self):
        return Matrix(self.coeff_ring, *self.coeff_ring_args, rows=self.rows)
    
    # gaussian elimination, computes rref
    def gauss(self):
        self.rows=sorted(self.rows, key=cmp_to_key(compare)) #sort the rows so that the matrix is upper triangular
        r, c=self.calc_dim() # r = no. of rows, c = no. of columns
        i=0
        j=0
        done=True
        while i<r and j<c:
            if self.rows[i][j] != self.coeff_ring.int(0, *self.coeff_ring_args): 
                # we are at the leading entry on the i-th row
                le=self.rows[i][j]
                if le != self.coeff_ring.int(1, *self.coeff_ring_args):
                    for x in range(j,c): #normalise so the leading entry is 1
                        self.rows[i][x] *= self.coeff_ring.int(1, *self.coeff_ring_args)/le
                k=i+1
                while k<r and self.rows[k][j] != self.coeff_ring.int(0, *self.coeff_ring_args):
                    le=self.rows[k][j]
                    for x in range(j,c):
                        self.rows[k][x] -= le*self.rows[i][x]
                    k+=1
                if k>i+1:
                    done=False #something was changed; we weren't in REF yet
                i=k
            else:
                j+=1 # we are not at the leading entry on the i-th row, move right
        if done:
            return
        self.gauss()

    def rank(self):
        M=self.copy()
        M.gauss()
        r,c=M.calc_dim()
        i=0
        j=0
        rank=0
        while i<r and j<c:
            if M.rows[i][j] != self.coeff_ring.int(0, *self.coeff_ring_args):
                rank+=1
                i+=1
            else:
                j+=1
        return rank
    
    def transpose(self):
        rows=[]
        r,c=self.calc_dim()
        for j in range(c):
            rows.append([])
            for i in range(r):
                rows[j].append(self.rows[i][j])
        return Matrix(self.coeff_ring, *self.coeff_ring_args, rows=rows)

    def calc_dim(self):
        return (len(self.rows), len(self.rows[0]))
        
    def __str__(self):
        dim=self.calc_dim()
        s=""
        for i in range(dim[0]):
            for j in range(dim[1]):
                s+=str(self.rows[i][j])
                if j<dim[1]-1:
                    s+=","
            s+="\n"
        return s

    # returns n*n identity matrix
    def id_mat(n, coeff_ring, *coeff_ring_args):
        rows=[]
        for i in range(n):
            rows.append([])
            for j in range(n):
                rows[i].append(coeff_ring.int(0, *coeff_ring_args))
            rows[i][i]=coeff_ring.int(1, *coeff_ring_args)
        return Matrix(coeff_ring, *coeff_ring_args, rows=rows)

    # given matrix in RREF, finds null space
    def null_space_basis(matrix, u): # u is the number of unknowns

        if len(matrix.rows) == 0: # everything is a solution to an empty system of equations
            id=[]
            for i in range(u):
                id.append([])
                for j in range(u):
                    id[i].append([matrix.coeff_ring.int(0, *matrix.coeff_ring_args)])
                id[i][i]=matrix.coeff_ring.int(1, *matrix.coeff_ring_args)
            return id
        
        f=0 # the leading zeros in the first row are free parameters; f counts the number of leading zeros
        while f<len(matrix.rows[0]) and matrix.rows[0][f] == matrix.coeff_ring.int(0, *matrix.coeff_ring_args):
            f+=1
        new_matrix=matrix.copy()
        # truncate first row and first f+1 columns
        new_matrix.rows=new_matrix.rows[1:]
        for i in range(len(new_matrix.rows)):
            new_matrix.rows[i]=new_matrix.rows[i][f+1:]
        partial_sol=Matrix.null_space_basis(new_matrix, u-f-1)
        sol=[]
        for i in range(f+len(partial_sol)): # f + len(partial_sol) is the total number of free parameters
            sol.append([])
            for j in range(u):
                if i<f:
                    if j == i:
                        sol[i].append(matrix.coeff_ring.int(1, *matrix.coeff_ring_args))
                    else:
                        sol[i].append(matrix.coeff_ring.int(0, *matrix.coeff_ring_args))
                else:
                    if j <= f:
                        sol[i].append(matrix.coeff_ring.int(0, *matrix.coeff_ring_args))
                    else:
                        sol[i][f] -= matrix.rows[0][j]*partial_sol[i-f][j-f-1]
                        sol[i].append(partial_sol[i-f][j-f-1])
        return sol

def compare(x,y): #comparison function for sorting rows during gaussian elimination
    # if first non-zero entry of x occurs earlier than in y, then x<y
    for i in range(len(x)):
        if not x[i].is_zero():
            return -1
        if not y[i].is_zero():
            return 1
    return 0