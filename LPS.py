"""
The main function, LPS(p,q), generates the Lubotzky, Phillips and Sarnak
Ramanujan graph for given odd prime inputs p,q. This implementation follows a
blueprint suggested by Randy Elzinga at
https://mast.queensu.ca/~ctardif/lps/LPSSup.pdf
"""

import sys
import csv
import time
import itertools

import numpy as np
from sympy.abc import a, b, c, d
# from sympy.solvers.diophantine import diop_general_sum_of_squares
from sympy.solvers.diophantine.diophantine import diop_general_sum_of_squares


"""
Basic modular arithmetic/number theoretic functions
"""

#Finds inverse of integer a mod m.
def modinv(a, m):
    while a < 0:
        a += m
    g, x, y = iterative_egcd(a, m)
    if g != 1:
        print(repr(a) + "^{-1} mod", m, "DNE")
        return None
    else:
        return x % m

#Used to find mod inverses
def iterative_egcd(a, b):
    x,y, u,v = 0,1, 1,0
    while a != 0:
        q,r = b//a,b%a; m,n = x-u*q,y-v*q # use x//y for floor "floor division"
        b,a, x,y, u,v = a,r, u,v, m,n
    return b, x, y

#Used to error check that the input p and q are prime
def isPrime(a):
    if a==1:
        return False
    else:
        return all(a % i for i in range(2, a))

# Taken from:
#https://codereview.stackexchange.com/questions/43210/
#tonelli-shanks-algorithm-implementation-of-prime-modular-square-root
#Determines if a is quadratic residual modulo odd prime p.
def legendre_symbol(a, p):
    a=int(a)
    p=int(p)
    ls = pow(a, (p - 1)//2, p)
    if ls == p - 1:
        return -1
    return ls

# Taken from:
#https://codereview.stackexchange.com/questions/43210/
#tonelli-shanks-algorithm-implementation-of-prime-modular-square-root
#Solves the equation x^2=a mod p and returns a list of x solutions. Uses the
#Tonelli-Shanks algorithm.
def prime_mod_sqrt(a, p):
    #for Python 3 compatibility with the pow function.
    a=int(a)
    p=int(p)
    a %= p

    # Simple case
    if a == 0:
        return [0]
    if p == 2:
        return [a]

    # Check solution existence on odd prime
    if legendre_symbol(a, p) != 1:
        return []

    # Simple case
    if p % 4 == 3:
        x = pow(a, (p + 1)//4, p)
        return [x, p-x]

    # Factor p-1 on the form q * 2^s (with Q odd)
    q, s = p - 1, 0
    while q % 2 == 0:
        s += 1
        q //= 2
    # Select a z which is a quadratic non resudue modulo p
    z = 1
    while legendre_symbol(z, p) != -1:
        z += 1
    c = pow(z, q, p)

    # Search for a solution
    x = pow(a, (q + 1)//2, p)
    t = pow(a, q, p)
    m = s
    while t != 1:
        # Find the lowest i such that t^(2^i) = 1
        i, e = 0, 2
        for i in range(1, m):
            if pow(t, e, p) == 1:
                break
            e *= 2

        # Update next value to iterate
        b = pow(c, 2**(m - i - 1), p)
        x = (x * b) % p
        t = (t * b * b) % p
        c = (b * b) % p
        m = i

    return [x, p-x]

"""
Generating set of the Cayley graph.
Elements are 2X2 matrices

[ x11 x12 ]
[ x21 x22 ]

stored as a vector/tuple [x11, x12, x21, x22].
"""
#Finds canonical solutions to diophantine equation a^2+b^2+c^2+d^2=0 mod p
#Sympy's solver only returns non-negative solutions, as an unordered set.
def get_canonical(p):
    sol = [sorted(list(s), key=lambda e: e%2 == 0) for s in
           diop_general_sum_of_squares(a**2 + b**2 + c**2 + d**2 - p,
                                       limit=8*(p+1))] #limit is max num of sol
    return sol

#Given canonical solution to diophantine equation, finds all solutions, up to
#signs and order.
def sol_expand(x):
    l=[]
    for i in range(0,4):
        if x[i] != 0:
            l.append([x[i],-x[i]])
        else:
            l.append([0])
    solSigns=itertools.product(l[0],l[1],l[2],l[3])
    sols=[]
    for x in solSigns:
        sols.extend(list(set(list(itertools.permutations(x)))))
    return np.array(list(set(sols)))

#Finds restricted set of solutions to diophantine equation that will be used
#to create generating set of the LPS graph.
def get_sol_generators(p):
    total_sol = np.concatenate(list(map(sol_expand, get_canonical(p))))
    print("total_sols")
    for sol in total_sol:
        print(sol)
    print("*" * 75)
    if p % 4 == 1:
        bool_sol = np.all([total_sol[:,0] > 0, total_sol[:,0] % 2 == 1], axis=0)
    elif p % 4 == 3:
        cond_1 = np.all([total_sol[:,0] > 0, total_sol[:,0] % 2 == 0], axis=0)
        cond_2 = np.all([total_sol[:,0] == 0, total_sol[:,1] > 0], axis=0)
        bool_sol = np.any([cond_1, cond_2], axis=0)
    else:
        print("Not a valid prime number.")
    return total_sol[bool_sol]

#Helper function to find solutions to x^2+y^2+1=0 mod p, used in creating
#generating set.
def solve_mod(q):
    for (x,y) in itertools.product(range(0,q),range(0,q)):
        if (x**2+y**2+1)%q==0:
            return x,y

#Returns a reduced set of generators for the LPS graph. If both s and s^-1 are
#in the generating set, only s is kept.
def get_generators(p,q):
    S=[]
    # chuck out the inverses
    gen_sols=get_sol_generators(p)
    print("gen_sols before reduction:")
    print(gen_sols)
    print("*" * 75)
    for i in range(gen_sols.shape[0]):
        idx_rows=None
        for j in range(i+1, gen_sols.shape[0]):
            if np.array_equal(np.concatenate([gen_sols[i, :1], -1*gen_sols[i, 1:]], axis=None), gen_sols[j, :]):
                idx_rows=j
                break
        if idx_rows!=None:
            gen_sols = np.delete(gen_sols, idx_rows, 0)
    x,y=solve_mod(q)
    print("gen_sols after reduction")
    print(gen_sols)
    for s in gen_sols:
        gen=[(s[0]+s[1]*x+s[3]*y)%q,
            (-s[1]*y+s[2]+s[3]*x)%q,
            (-s[1]*y-s[2]+s[3]*x)%q,
            (s[0]-s[1]*x-s[3]*y)%q]
        print("s: ", s)
        print("gen: ", gen)
        S.append(tuple(gen))
    return S

"""
Generate the group elements of PSL or PGL.
Elements are 2X2 matrices

[ x11 x12 ]
[ x21 x22 ]

stored as a vector/tuple [x11, x12, x21, x22].
"""

#Helper function returning first nonzero element in 2X2 matrix, recorded as a
#vector [x11,x12,x21,x22].
def first_nonzero(x,q):
    if x[0]%q !=0:
        val=x[0]%q
    else:
        val=x[2]%q
    return val

#Returns matrix in standard form for PGL. Each element is a representative of
#the equivalence class in which the first nonzero element in the first
#column is 1.
def sMat_PGL(x,q):
    if x[0]%q !=0:
        x=[(modinv(x[0],q)*xi)%q for xi in x]
    else:
        x=[(modinv(x[2],q)*xi)%q for xi in x]
    return x

#Returns matrix in standard form for PSL.
#The determinant of each element X is a square of some sigma mod q. The
#representative is either sigma^(-1)*X or -sigma^(-1)*X chosen so that first
#nonzero element of the first column is in {1,...,q-1/2}
def sMat_PSL(x,q):
    det=(x[0]*x[3]-x[1]*x[2])%q
    sigma=modinv(prime_mod_sqrt(det,q)[0],q)
    x_pos=[(sigma*xi)%q for xi in x]
    x_neg=[(-1*sigma*xi)%q for xi in x]
    if first_nonzero(x_pos,q) in set(range(1,int((q-1)/2)+1)):
        return x_pos
    else:
        return x_neg

#Returns elements of PGL(2,q).
def PGL(q):
    C=[]
    for x in itertools.product(range(0,q),range(0,q),range(0,q),range(0,q)):
        if (x[0]*x[3]-x[1]*x[2])%q != 0 :
            C.append(tuple(sMat_PGL(x,q)))
    return set(C)

#Returns elements of PSL(2,q). First finds element of PGL, then constructs
#PSL as the subgroup.
def PSL(q):
    C=PGL(q)
    P=[]
    for c in C:
        det=(c[0]*c[3]-c[1]*c[2])%q
        if prime_mod_sqrt(det,q): #if determinant is quadratic residue
            P.append(tuple(sMat_PSL(c,q))) #append matrix in standard form
    return P

"""
Generate the edge set of LPS(p,q).
For speed, each group element (represented in standard form) is mapped to a
unique integer between 1 and the group order.
"""

#Maps representative element of PGL to int between 1 and group order.
def mat2int_PGL(s,q):
    if s[0]%q != 0:
        x=s[2]
        a=(x+1)
        b=(s[1]*(x+1)-s[3])%q
        c=(s[3]-s[1]*x)%q
    else:
        a=0
        b=s[3]
        c=s[1]
    if c == 0:
        print('WHOA, MISTAKE')
    return (c-1)*(q**2+q)+b*(q+1)+a+1

# Map integer between 1 and group order to element of PGL
def int2mat_PGL(k,q):
    a = (k-1)%(q+1)
    kprime = (k - a - 1)/(q+1)
    b = kprime%q
    c = (kprime - b)/q+1
    if a == 0:
        return sMat_PGL([0,c, 1,b ],q)
    else:
        return sMat_PGL([1,b+c, a-1 , b*(a-1)+c*a ], q)
        

#Maps a representative element of PSL to int between 1 and group order.
def mat2int_PSL(s,q):
    if s[0]%q != 0:
        a=s[0]
        b=s[2]
        c=s[1]
    else:
        a=s[2]
        b=q
        c=s[3]
    return c*(q**2-1)/2 +b*(q-1)/2 + a

# Maps an integer between 1 and the group order to the associated element of PSL
def int2mat_PSL(k,q):
    a = k%( (q-1)/2)
    if a == 0:
        a = (q-1)/2
    kprime = (k-a)/( (q-1)/2)
    b = kprime % (q+1)
    c = (kprime - b)/(q+1)
    if b == q:
        return sMat_PSL([0,q-modinv(a,q),a,c],q)
    else:
        return sMat_PSL([a,c,b, (1+b*c)*modinv(a,q) ],q)

    
#OBSOLETE! Slow method for generating PGL edge set. Searching all vertices to
#find index of endpoint of an edge is too expensive.
def create_edges_OLD(p,q):
    E=[]
    S=get_generators(p,q)
    V=list(PGL(q))
    t = time.time()
    for s in S:
        srcInd=0
        for v in V:
            vs_prod=[v[0]*s[0]+v[1]*s[2], v[0]*s[1]+v[1]*s[3],
                    v[2]*s[0]+v[3]*s[2], v[2]*s[1]+v[3]*s[3]]
            vs_prod=sMat_PGL(vs_prod,q)
            tarInd=np.where((np.array(list(V)) == vs_prod).all(axis=1))
            tarInd=tarInd[0][0]
            E.append([srcInd, tarInd])
            srcInd=srcInd+1
    elapsed = time.time() - t
    print(repr(int(elapsed))+ ' seconds to compute')
    return E

#Main function, generates LPS graph for inputs p and q.
def LPS(p,q):
    if not isPrime(p) or not isPrime(q) or p==2 or q==2 or p==q:
        print('Error: p and q must be distinct, odd primes')
        return None
    E=[]
    gens={}
    nodes={}
    S=get_generators(p,q)
    if legendre_symbol(p,q) == -1: #PGL
        if p+1 > ((q**3 - q) - 1):
            print('Error: the degree cannot exceed the group order')
            return None
        print('Group is PGL, ' + repr(int(p+1))+ '-regular bipartite graph on '+
            repr(int((q**3-q)))+ ' vertices')
        t = time.time()
        V=list(PGL(q))
        elapsed = time.time() - t
        print('Vertices generated in '+ repr(int(elapsed)) + ' seconds')
        t = time.time()
        for s in S:
            s=sMat_PGL(s,q)
            gens[mat2int_PGL(s,q)]=s
            for v in V:
                nodes[mat2int_PGL(v,q)]=v
                vs_prod=[(v[0]*s[0]+v[1]*s[2])%q, (v[0]*s[1]+v[1]*s[3])%q,
                        (v[2]*s[0]+v[3]*s[2])%q, (v[2]*s[1]+v[3]*s[3])%q]
                vs_prod=sMat_PGL(vs_prod,q)
                tarInd=mat2int_PGL(vs_prod,q)
                srcInd=mat2int_PGL(v,q)
                E.append([srcInd, tarInd])
        elapsed = time.time() - t
        print('Edges generated in '+ repr(int(elapsed))+ ' seconds')
    elif legendre_symbol(p,q)==1: #PSL
        if p+1 > ((q**3 - q)/2 - 1):
            print('Error: the vertex degree cannot exceed the group order')
            return None
        print('Group is PSL, ' + repr(int(p+1))+ '-regular graph on '+
            repr(int((q**3-q)/2))+ ' vertices')
        t = time.time()
        V=list(PSL(q))
        elapsed = time.time() - t
        print('Vertices generated in '+ repr(int(elapsed)) + ' seconds')
        t = time.time()
        for s in S:
            s=sMat_PSL(s,q)
            gens[mat2int_PSL(s,q)]=s
            for v in V:
                nodes[mat2int_PSL(v,q)]=v
                vs_prod=[(v[0]*s[0]+v[1]*s[2])%q, (v[0]*s[1]+v[1]*s[3])%q,
                        (v[2]*s[0]+v[3]*s[2])%q, (v[2]*s[1]+v[3]*s[3])%q]
                vs_prod=sMat_PSL(vs_prod,q)
                tarInd=mat2int_PSL(vs_prod,q)
                srcInd=mat2int_PSL(v,q)
                E.append([int(srcInd), int(tarInd)])
        elapsed = time.time() - t
        print('Edges generated in '+repr(int(elapsed))+ ' seconds')
    return E, gens, nodes


if __name__ == '__main__':
    # p = int(sys.argv[-2])
    # q = int(sys.argv[-1])
    p = 23
    q = 5
    E = LPS(p,q)
    if E == None:
        sys.exit()
    ## save edge list
    with open("LPS_{}_{}.csv".format(p, q), "w") as f:
        writer = csv.writer(f)
        for e in E:
            writer.writerow(e)
        f.close()
