from sage.modules.free_module_integer import IntegerLattice
from sage.stats.distributions.discrete_gaussian_polynomial import DiscreteGaussianDistributionPolynomialSampler
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler

from almost_inverse import inverse

precision = 20

# elem is in O_K
def mod_q(elem, q):
    def half(a):
        if a < q/2:
	    return a
 	else:
	    return a-q

    z = elem.parent().gen()
    n = len(elem.list())
    
    return sum([half(mod(coerce(Integer, elem[i]), q).lift()) * z^i for i in range(n)])

def NTRU(h, K, q):
    basis = K.integral_basis()
    H = h.matrix()

    return IntegerLattice(block_matrix( [[1, H], [0, q]]), lll_reduce=True)

def NTRU_subfield(hprime, q, nprime, r):
    z = hprime.parent().gen()

    mat = []
    for i in range(nprime):
        coordinate = (hprime * z^(r*i)).vector().list()
        mat.append( [coordinate[r*j] for j in range(nprime)] ) 

    Hprime = matrix(mat)
    
    return IntegerLattice(block_matrix([ [1, Hprime], [0, q]]), lll_reduce=True)

def attack(m, q, r = 4, sigma = 3.0, subfield_only=False):
    K.<z> = CyclotomicField(m)
    OK = K.ring_of_integers()

    n = euler_phi(m)
    mprime = m / r 
    nprime = euler_phi(mprime)

    phim = a^n + 1    
    R.<a> = ZZ['a']
    D = DiscreteGaussianDistributionIntegerSampler(sigma)
    
    while True:
        f = sum([D()*z^i for i in range(n)])
        fx = sum([f[i]*a^i for i in range(n)])

        res = inverse(fx, phim, q)
        if res[0]:
            f_inv = sum([res[1][i]*z^i for i in range(n)])
            print "f_inv * f = %s (mod %d)" %((f*f_inv).mod(q), q)
            break
        
    g = sum([D()*z^i for i in range(n)])

    #h = [g*f^{-1)]_q
    h = (g*f_inv).mod(q)

    lognorm_f = log_b(f.vector().norm(), 2)
    lognorm_g = log_b(g.vector().norm(), 2)

    print "f*h - g = %s" %  (f*h-g).mod(q)
    print "log q = ", log_b(q, 2).n(precision)
    print "log |f| = %s, log |g| = %s", %( lognorm_f.n(precision), 
                                            lognorm_g.n(precision) )
    print "log |(f,g)| = ", log_b(sqrt(f.vector().norm()^2 + g.vector().norm()^2), 2).n(precision)
    
    fprime = prod([tau(f) for tau in K.galois_group() if tau(z^r) == z^r])
    gprime = prod([tau(g) for tau in K.galois_group() if tau(z^r) == z^r])
    hprime = prod([tau(h) for tau in K.galois_group() if tau(z^r) == z^r])

    lognorm_fp = log_b(fprime.vector().norm(), 2)
    lognorm_gp = log_b(gprime.vector().norm(), 2)

    print "%d * log |f| - log |f'| = %s" %(r, r * lognorm_f.n(precision) - lognorm_fp.n(precision))
    print "log |(f', g')| = ", log_b(sqrt(fprime.vector().norm()^2 + gprime.vector().norm()^2), 2).n(precision)

    #(fprime, gprime) lies in the lattice \Lambda_hprime^q
    print "f'*h' - g' = %s " % (hprime*fprime - gprime).mod(q)

    if not subfield_only:
        ntru_full = NTRU(h, K, q)
        full_sv = ntru_full.shortest_vector()
    
        print "log |v| = %s" % log_b(full_sv.norm(), 2).n(precision)

    ntru_subfield = NTRU_subfield(hprime, q, nprime, r)
    
    sub_sv = ntru_subfield.shortest_vector()    
    xprime = sum([coerce(Integer, sub_sv[i])*z^(r*i) for i in range(nprime)])

    #test if xprime belongs to <fprime>
    mat = []
    for i in range(nprime):
        coordinate = (fprime * z^(r*i)).vector().list()
        mat.append( [coordinate[r*j] for j in range(nprime)] ) 
    FL = IntegerLattice(mat)
    print sub_sv[:nprime] in FL

    x = xprime 
    y = mod_q(x*h, q)
    
    print "log |(x,y)| = ", log_b( sqrt(x.vector().norm()^2 + y.vector().norm()^2), 2).n(precision)
    return (x, y, f,g, h)

# f,g is a small polynomial: each coefficient is either 0,1,-1

for _ in range(2):
    attack(m=256, q = next_prime(2^100), sigma=4, subfield_only=True)
for _ in range(2):
    attack(m=512, q = next_prime(2^100), sigma=4, subfield_only=True)