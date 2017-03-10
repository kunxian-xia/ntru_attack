from sage.modules.free_module_integer import IntegerLattice
from sage.stats.distributions.discrete_gaussian_polynomial import DiscreteGaussianDistributionPolynomialSampler
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
from sage.rings.arith import next_prime,euler_phi
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.integer_ring import IntegerRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.functional import log
from sage.functions.other import sqrt
from almost_inverse import inverse

precision = 20

def trace(t, G):
    return sum( [tau(t) for tau in G] )
def norm(t, G):
    return prod( [tau(t) for tau in G] ) 

# elem is in O_K
def mod_q(elem, q):
    def half(a):
        if a < q/2:
	    return a
 	else:
	    return a-q

    z = elem.parent().gen()
    n = len(elem.list())
    
    return sum([half(mod(coerce(Integer, elem[i]), q).lift()) * z**i for i in range(n)])

def NTRU(h, K, q):
    basis = K.integral_basis()
    H = h.matrix()

    return IntegerLattice(block_matrix( [[1, H], [0, q]]), lll_reduce=True)

def NTRU_subfield(hprime, q, nprime, r):
    z = hprime.parent().gen()

    mat = []
    for i in range(nprime):
        coordinate = (hprime * z**(r*i)).vector().list()
        mat.append( [coordinate[r*j] for j in range(nprime)] ) 

    Hprime = matrix(mat)
    
    return IntegerLattice(block_matrix([ [1, Hprime], [0, q]]), lll_reduce=True)

def attack(m, q, r = 4, sigma = 3.0, subfield_only=False):
    K = CyclotomicField(m, 'z')
    z = K.gen()
    OK = K.ring_of_integers()
    G = K.galois_group()

    n = euler_phi(m)
    mprime = m / r 
    nprime = euler_phi(mprime)
    Gprime = [tau for tau in G if tau(z**r) == z**r]

    R = PolynomialRing(IntegerRing(),'a')
    a = R.gen()
    phim = a**n + 1  
    D = DiscreteGaussianDistributionIntegerSampler(sigma)
    
    print "sampling f,g"
    while True:
        f = sum([D()*z**i for i in range(n)])
        fx = sum([f[i]*a**i for i in range(n)])

        res = inverse(fx, phim, q)
        if res[0]:
            f_inv = sum([res[1][i]*z**i for i in range(n)])
            print "f_inv * f = %s (mod %d)" %((f*f_inv).mod(q), q)
            break
  
    g = sum([D()*z**i for i in range(n)])
    print "done sampling f, g"

    #h = [g*f^{-1)]_q
    h = (g*f_inv).mod(q)

    lognorm_f = log(f.vector().norm(), 2)
    lognorm_g = log(g.vector().norm(), 2)

    print "f*h - g = %s" %  (f*h-g).mod(q)
    print "log q = ", log(q, 2).n(precision)
    print "log |f| = %s, log |g| = %s" %( lognorm_f.n(precision), 
                                         lognorm_g.n(precision) )
    print "log |(f,g)| = ", log(sqrt(f.vector().norm()**2 + g.vector().norm()**2), 2).n(precision)
    
    print "begin computing N(f), N(g), N(h), Tr(h), fbar"
    fprime = norm(f, Gprime)
    gprime = norm(g, Gprime)
    hprime = norm(h, Gprime).mod(q)
    htr = trace(h, Gprime)
    fbar = prod([tau(f) for tau in Gprime[1:] ])
    print "end computing N(f), N(g), N(h), Tr(h), fbar"

    lognorm_fp = log(fprime.vector().norm(), 2)
    lognorm_gp = log(gprime.vector().norm(), 2)

    print "%d * log |f| - log |f'| = %s" %(r, r * lognorm_f.n(precision) - lognorm_fp.n(precision))
    print "log |(f', g')| = ", log(sqrt(fprime.vector().norm()**2 + gprime.vector().norm()**2), 2).n(precision)
    print "log |N(f), Tr(g fbar)| = ", log( sqrt(fprime.vector().norm()**2 + 
                                                trace(g*fbar, Gprime).vector().norm()**2), 2).n(precision)
    
    #(fprime, gprime) lies in the lattice \Lambda_hprime^q
    print "f'*h' - g' = %s " % (hprime*fprime - gprime).mod(q)
    print "N(f) Tr(h) - Tr(g fbar) = %s" % (htr*fprime - trace(g*fbar, Gprime)).mod(q)

    if not subfield_only:
        ntru_full = NTRU(h, K, q)
        full_sv = ntru_full.shortest_vector()
    
        print "log |v| = %s" % log(full_sv.norm(), 2).n(precision)

    ntru_subfield = NTRU_subfield(hprime, q, nprime, r)
    ntru_trace_subfield = NTRU_subfield(htr, q, nprime, r)

    print "begin computing Shortest Vector of subfield lattice"
    norm_sv = ntru_subfield.shortest_vector() 
    tr_sv = ntru_trace_subfield.shortest_vector()
    print "end computing Shortest Vector of subfield lattice"

    norm_xp = sum([coerce(Integer, norm_sv[i])*z**(r*i) for i in range(nprime)])
    tr_xp = sum([coerce(Integer, tr_sv[i])*z**(r*i) for i in range(nprime)])

    #test if xprime belongs to <fprime>
    mat = []
    for i in range(nprime):
        coordinate = (fprime * z**(r*i)).vector().list()
        mat.append( [coordinate[r*j] for j in range(nprime)] ) 
    FL = IntegerLattice(mat)
    print norm_sv[:nprime] in FL
    print tr_sv[:nprime] in FL

    norm_x = norm_xp 
    norm_y = mod_q(norm_x*h, q)

    tr_x = tr_xp
    tr_y = mod_q(tr_x *h, q)
    
    print "Norm map: log |(x,y)| = ", log( sqrt(norm_x.vector().norm()**2 + norm_y.vector().norm()**2), 2).n(precision)
    print "Trace map: log |(x,y)| = ", log( sqrt(tr_x.vector().norm()**2 + tr_y.vector().norm()**2), 2).n(precision)

# f,g is a small polynomial: each coefficient is either 0,1,-1

for _ in range(4):
    attack(m=128, q = next_prime(2**10), r = 4, sigma=2, subfield_only=True)
