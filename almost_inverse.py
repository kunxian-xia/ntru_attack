from sage.rings.finite_rings.integer_mod import Mod

# f is assumed to be integral polynomial in x
def mod_q(f, q):
    n = f.degree()
    x = f.parent().gen()
    
    return sum([Mod(f[i], q)*(x**i) for i in range(n+1)])

# a \in Z[x]
# m \in Z[x]
# find b s.t. a*b = 1 (mod p)
# p is prime 
def inverse_p(a, m, p):
    a = mod_q(a, p)
    m = mod_q(m, p)

    return almost_inverse(a, m)
    
def inverse(a, m, q):
    if q.is_prime_power():
        (p, r) = q.perfect_power()
        res = inverse_p(a, m, p)
        if res[0]:
            b = res[1]
            x = a.parent().gen()
            
            b = sum([b[i].lift()*(x**i) for i in range(b.degree()+1)])
            s = p
            while (s < q):
                s = s**2
                b = (b*(2-a*b)).mod(m)
                b = sum([ Mod(b[i], s).lift()*(x**i) for i in range(b.degree()+1)])
            res = (res[0], b)
        return res
    else:
        raise ValueError
    
#a, m are polynomials in Z_p[x]
# m = x^N + 1, N = 2^k or 
# m = x^N - 1, N is prime
def almost_inverse(a, m):
    f, g = a, m
    b, c = 1, 0
    N, k = m.degree(), 0
    x = m.parent().gen()
    
    while True:

        if f == 0:
            return (False, )
        while f[0] == 0:
            f = f >> 1
            c = c*x
            k = k+1
        if f.degree() == 0:
            t = Mod(-k, N).lift()
            s = (k+t)/N
            b = (x**t)*(f[0]**(-1))*b
            if Mod(s, 2).lift() == 1:
                b = -b
            res= (True, b.mod(m))

            return res
        if f.degree() < g.degree():
            f,g = g, f
            b,c = c, b
        u = f[0]*(g[0]**(-1))
        f = f- u*g
        b = b- u*c
    
