# Use at own risk!

import json

# Initialize global polynomial ring with variable t
R.<t> = QQ['t']

# Dictionary of g_polynomials g([F,P]) of all faces F of a lattice polytope P (initialized with backend='normaliz'):
def g_dictionary(P):
    poset = P.face_lattice()
    g_dict = dict((F.as_polyhedron(),-1) for F in poset)
    g_dict[P] = 1
    ls_pos = poset.sorted(list(poset))
    ls_pos.remove(poset.top())
    for k in range(0,len(ls_pos)):
        F = ls_pos[len(ls_pos)-1-k]
        interval_pos = poset.subposet(poset.closed_interval(F,poset.top()))
        ls_interval_pos = list(interval_pos)
        ls_interval_pos.remove(interval_pos.bottom())
        h = 0
        for G in ls_interval_pos:
                h += (t-1)^(interval_pos.rank(G)-1) * g_dict[G.as_polyhedron()]
        g_help = (1-t)*h
        e = ceil(1.0*interval_pos.rank()/2) - 1
        coeff = g_help.coefficients(sparse=False)
        g = 0
        for n in range(0,min(e+1,len(coeff))):
            g += t^n * coeff[n]
        g_dict[F.as_polyhedron()] = g
    return g_dict

# The local h*-polynomial (also called l*-polynomial) of a lattice polytope P (initialized with backend='normaliz'):
def l_star_polynomial(P):
    g_dict = g_dictionary(P)
    lstar = (-1)^(P.dim() + 1) * g_dict[P.face_lattice().bottom().as_polyhedron()]
    poset = P.face_lattice()
    lst_pos = list(poset)
    lst_pos.remove(poset.bottom())
    for F in lst_pos:
        FPolytope = F.as_polyhedron()
        hstarF = FPolytope.ehrhart_series(t).numerator()
        lstar += hstarF * g_dict[FPolytope] * (-1)^(P.dim() - FPolytope.dim())
    return lstar

# For a lattice simplex Delta the g-polynomial g([F,Delta]) always equals 1
def l_star_polynomial_simplex(Delta):
    lstar = (-1)^(Delta.dim() + 1)
    poset = Delta.face_lattice()
    lst_pos = list(poset)
    lst_pos.remove(poset.bottom())
    for F in lst_pos:
        FPolytope = F.as_polyhedron()
        hstarF = FPolytope.ehrhart_series(t).numerator()
        lstar += hstarF * (-1)^(Delta.dim() - FPolytope.dim())
    return lstar

def fat_or_flat(P):
    Polyt = Polyhedron(P,backend='normaliz')
    lstar = l_star_polynomial(Polyt)
    if lstar != 0:
        hstar = Polyt.ehrhart_series().numerator()
        if hstar.degree() != lstar.degree():
            return False
    return True

def search_not_fat_or_flat(ls_polytopes):
    ls_non_fat_or_flat = []
    for P in ls_polytopes:
        if not fat_or_flat(P):
            ls_non_fat_or_flat.append(P)
    return ls_non_fat_or_flat

def search_thin(ls_polytopes):
    ls_thin_polytopes = []
    for Polyt in ls_polytopes:
        P = Polyhedron(Polyt,backend='normaliz')
        if l_star_polynomial(P) == 0:
            ls_thin_polytopes.append(P)
    return ls_thin_polytopes

def search_thin_simplex(ls_simplices):
    ls_thin_simplices = []
    for Polyt in ls_simplices:
        Delta = Polyhedron(Polyt,backend='normaliz')
        lstar = l_star_polynomial_simplex(Delta)
        print(lstar)
        if lstar == 0:
            ls_thin_simplices.append(Polyt)
    return ls_thin_simplices

def search_thin_simplex_homogenized(ls_simplices):
    ls_thin_simplices = []
    for Polyt in ls_simplices:
        Delta = Polyhedron(Polyt,backend='normaliz')
        lstar = l_star_polynomial_simplex(Delta)
        print(lstar)
        if lstar == 0:
            ls_thin_simplices.append(homogenize(Polyt))
    return ls_thin_simplices

def search_thin_4_polytopes(ls_polytopes):
    ls_thin_polytopes = []
    i = 1
    for Polyt in ls_polytopes:
        P = Polyhedron(Polyt,backend='normaliz')
        nrVert = len(P.vertices_list())
        if nrVert == 5:
            if l_star_polynomial_simplex(P) == 0:
                ls_thin_polytopes.append(P)
                print("Thin simplex!"+str(i))
        else:
            if l_star_polynomial(P) == 0:
                ls_thin_polytopes.append(P)
                print("Thin polytope!"+str(i))
        i = i+1
    return ls_thin_polytopes

def fat_or_flat_simplex(Delta):
    Polyt = Polyhedron(Delta,backend='normaliz')
    lstar = l_star_polynomial_simplex(Polyt)
    if lstar != 0:
        hstar = Polyt.ehrhart_series().numerator()
        if hstar.degree() != lstar.degree():
            return False
    return True

def search_not_fat_or_flat_simplex(ls_simplices):
    ls_non_fat_or_flat = []
    for Delta in ls_simplices:
        if not fat_or_flat_simplex(Delta):
            ls_non_fat_or_flat.append(Delta)
    return ls_non_fat_or_flat

def thin_but_not_trivially_thin(ls_polytopes):
    ls_not_trivially_thin = []
    for P in ls_polytopes:
        Polyt = Polyhedron(P,backend='normaliz')
        hstar = Polyt.ehrhart_series().numerator()
        if 2*hstar.degree() > Polyt.dim():
            if l_star_polynomial(Polyt) == 0:
                ls_not_trivially_thin.append(P)
    return ls_not_trivially_thin

def thin_but_not_trivially_thin_simplex(ls_simplices):
    ls_not_trivially_thin = []
    for Delta in ls_simplices:
        Polyt = Polyhedron(Delta,backend='normaliz')
        hstar = Polyt.ehrhart_series().numerator()
        if 2*hstar.degree() > Polyt.dim():
            if l_star_polynomial_simplex(Polyt) == 0:
                ls_not_trivially_thin.append(Delta)
    return ls_not_trivially_thin

# With respect to the classification by Francisco Santos of empty lattice 4-simplices:
def make_empty_4_simplex(Vol, b):
    i = 0
    inv = 0
    for i in range(0,5):
        if gcd(Vol,b[i]) == 1:
            inv = b[i].inverse_mod(Vol)
            break
    b2 = [((-inv*x) % Vol) for x in b]
    b2[i] = -1
    # b2Sum must be multiple of Vol:
    b2Sum = sum(b2)
    b2.pop(i)
    b2[0] = b2[0] + Vol - b2Sum
    Delta = Polyhedron([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],b2],backend='normaliz')
    return Delta

# For processing a polytope in Polymake:
def homogenize(P):
    return [([1]+pt) for pt in P]

# Only for full-dimensional P!
def isSpanning(P):
    Polyt = Polyhedron(P,backend='normaliz')
    Polyt = Polyt.translation([-n for n in P[0]])
    A = matrix(Polyt.integral_points())
    A = A.transpose()
    elemDiv = A.elementary_divisors()
    if elemDiv[len(elemDiv)-1] == 1:
        return true
    else:
        return false

# P is given as list of points whose convex hull is the polytope for which we want to check if it is a join.
def isJoin(P):
    Polyt = Polyhedron(P)
    d = Polyt.dim()
    vertices = Polyt.vertices_list()
    if len(vertices) == (d+1):
        return True
    vertP = Set(tuple(v) for v in vertices)
    for i in range(0,floor((d-1)/2)+1):
        for F in Polyt.faces(i):
            vertF = Set(tuple(v) for v in ((F.as_polyhedron()).vertices_list()))
            Q = Polyhedron([list(v) for v in list(vertP.difference(vertF))])
            if Q.dim() == (d-i-1):
                return True
    return False

# The cayleyJoinFactors function below needs the input to define a *full-dimensional* polyhedron!
# P is assumed to be given as list of points whose convex hull is the polytope.
# The function returns the list of all pairs (as lists of length 2) of faces of P such that P is 
# the Cayley join of the two faces. Of course, the output list can be empty. The dimension of the first factor 
# of each pair cannot exceed the dimension of the second, by construction.
def cayleyJoinFactors(P):
    result = []
    Polyt = Polyhedron(P,backend='normaliz')
    d = Polyt.dim()
    vertices = Polyt.vertices_list()
    vertP = Set(tuple(v) for v in vertices)
    for i in range(0,floor((d-1)/2)+1):
        for F in Polyt.faces(i):
            vertF = Set(tuple(v) for v in ((F.as_polyhedron()).vertices_list()))
            Q = Polyhedron([list(v) for v in list(vertP.difference(vertF))],backend='normaliz')
            if Q.dim() == (d-i-1):
                FPolyt = Polyhedron([list(v) for v in list(vertF)],backend='normaliz')
                vF = FPolyt.vertices_list()[0]
                vQ = Q.vertices_list()[0]
                FPolyt2 = FPolyt.translation([-n for n in vF])
                Q2 = Q.translation([-n for n in vQ])
                lsF = FPolyt2.vertices_list()
                lsF.remove([0 for i in range(0,len(vertices[0]))])
                lsQ = Q2.vertices_list()
                lsQ.remove([0 for i in range(0,len(vertices[0]))])
                ls = lsF + lsQ
                normal = ((matrix(ls).transpose()).kernel()).basis()[0]
                if abs(normal.dot_product(vector(ZZ,vF)) - normal.dot_product(vector(ZZ,vQ))) == 1:
                    result.append([FPolyt.vertices_list(),Q.vertices_list()])
    return result

# Only for full-dimensional polytopes!
def isCayleyJoin(P):
    Polyt = Polyhedron(P,backend='normaliz')
    d = Polyt.dim()
    vertices = Polyt.vertices_list()
    vertP = Set(tuple(v) for v in vertices)
    for i in range(0,floor((d-1)/2)+1):
        for F in Polyt.faces(i):
            vertF = Set(tuple(v) for v in ((F.as_polyhedron()).vertices_list()))
            Q = Polyhedron([list(v) for v in list(vertP.difference(vertF))],backend='normaliz')
            if Q.dim() == (d-i-1):
                FPolyt = Polyhedron([list(v) for v in list(vertF)],backend='normaliz')
                vF = FPolyt.vertices_list()[0]
                vQ = Q.vertices_list()[0]
                FPolyt2 = FPolyt.translation([-n for n in vF])
                Q2 = Q.translation([-n for n in vQ])
                lsF = FPolyt2.vertices_list()
                lsF.remove([0 for i in range(0,len(vertices[0]))])
                lsQ = Q2.vertices_list()
                lsQ.remove([0 for i in range(0,len(vertices[0]))])
                ls = lsF + lsQ
                normal = ((matrix(ls).transpose()).kernel()).basis()[0]
                if abs(normal.dot_product(vector(ZZ,vF)) - normal.dot_product(vector(ZZ,vQ))) == 1:
                    return True
    return False
