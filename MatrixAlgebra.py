'''
some efforts to implement Matrix Algebra:
in particular, PLU decomposition:

v.0.0.1 working
v.0.0.2 (to add/rewrite) :
    + list vs. generator
    + function for eye martix
    + column vector as columns not as rows
    + classes vs. procedures
    + error handling
'''
def t(m):
    """transpose the matrix
    """
    return zip(*m)
    
def mprint(m):
    """print the matrix
    """
    for row in m:
        print(*("{:5.2f}".format(a) for a in row))

def mult(a,b):
    """multiply two matrices
    """
    return [[sum(i*j for i,j in zip(row, col)) for col in zip(*b)] for row in a]
    
def lu(a):
    """recursive PLU decomposition
    """
    if len(a)==1:
        return [[1]], [[1]], a
    else :
        q = mpivot(a)
        b = mult (q,a)
        b11,r,c,c1 = mblocks(b)
        k = len(c1)
        b0 = [[float(i==j)/b11[0][0] for j in range(k)] for i in range(k)]
        b1 = list(mult(t(c),r))
        b1 = mult(b0,b1)
        b1 = sub(c1,b1)
        #mprint(b1)
       
        return lumap(*lu(b1),b11, r,c,q, k)


def lumap( p1,l1,u1, b11, r, c, q, k):
    """'formating' recursive output
    """
    k = len(p1)
    p = mbuild([[1]], [[0]*k],[[0]*k], p1)
    p = mult(t(q),p)
    
    p2 = [[float(i==j)/b11[0][0] for j in range(k)] for i in range(k)]
    p2 = list(mult(p2, list(t(c))))
    l = mbuild([[1]], [[0]*k], list(t(mult(t(p),p2))), l1)
    
    u = mbuild(b11, r, [[0]*k], u1)
    #print(k)
    #mprint(l)
    
    return p, l, u
        
def mblocks(m):
    """decompose matrix on blocks:
        |  m1 row |
    m = |
        | col  m2 |
    """
    m1 = [[m[0][0]]]
    row = [m[0][1:]]
    col =[list(t(m))[0][1:]]  
    m2 = list(t(list(t(m[1:]))[1:]))
    #mprint(m1)
    #mprint(row)
    #mprint(col)
    #mprint(m2)
    return m1, row, col, m2

def mbuild(m1, row, col, m2):
    """build matrix from blocks
    """
    mx = col
    mx += list(t(m2))
    m = [m1[0]+row[0]]
    m += list(t(mx))
    #mprint(m)
    return m
    
    
def mpivot(m):
    """pivoting matrix for m, used in Doolittle's method."""
    n = len(m)                                                                                                                                                                                         
    eye = [[float(i==j) for i in range(n)] for j in range(n)]                                                                                                                                                                                                                                                                                                                                                          
    #mprint(eye)
    #print("-----")
    for j in range(n):
        row = max(range(j, n), key=lambda i: abs(m[i][j]))
        if j != row:
            
            # Swap the rows                                                                                                                                                                                                                            
            eye[j], eye[row] = eye[row], eye[j]

    return eye

def sub(a,b):
    """substructing matrices a -b 
    """
    n = len(a)
    return [[a[i][j]- b[i][j] for j in range(n)] for i in range(n)]  
    
    
#a = [[2,4,1],
#     [4,-10,2],
#     [1,2,4]]
a = [[1,3,5],
     [2,4,7],
     [1,1,0]]
     

# some debugging below    
#b = [[1,0,0],
#     [-4/7,1,0],
#     [-1/7,0,1]]
     

#d = mult(b,a)    
#p = mpivot(a)

#c = mult(p,a)
#d = mult(b,c)
#mprint(c)
#mprint(d)
#lu(a, p=[], l=[], u=[])
#m1, row, col, m2 = mblocks(a)
#print(m1)
#mprint(t(t(row)))
#print(col)
#print(m2)

#c= mbuild(m1, row, col, m2)
#print(c)
for m in lu(a):
    mprint(m)
    print("----")
p,l,u = lu(a)
mprint(mult(l,u))
