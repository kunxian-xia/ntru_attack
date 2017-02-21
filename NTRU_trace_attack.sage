"""
map h <- NTRU(O_K, q) to tr_{K/L}(h)

"""
def trace(t, G):
    return sum( [tau(t) for tau in G] )
def norm(t, G):
    return prod( [tau(t) for tau in G] ) 

def mod_q(t, q, z):
    n = len(t.list())
    return sum( [mod(t[i], q).lift()*z^i for i in range(n)] )

def trace_attack(m, q, r):
    K.<z> = CyclotomicField(m)
    R = K.ring_of_integers()
    n = euler_phi(m)
    G = K.galois_group()

    #subfield 
    Gprime = [tau for tau in G if tau(z^r) == z^r]

    #f, g
    f = R.random_element()
    g = R.random_element()

    finv = f.inverse_mod(R.ideal(q))

    #h
    h = mod_q(g*finv, q, z)
    
    #compute trace of h, 
    htr = trace(h, Gprime)

    fno = norm(f, Gprime)
    finv_no = norm(finv, Gprime)
    fbar = prod([tau(f) for tau in Gprime[1:] ])

    #print fbar*g
    print mod_q(trace(f, Gprime), r, z)
    print mod_q(trace(g, Gprime), r, z)
    #print mod_q(trace(fbar*g, Gprime), r, z)

    #assert fno*finv_no = 1 (mod q)
    print mod_q(fno*finv_no, q, z)
    
    print mod_q(htr*fno - trace(g*fbar, Gprime), q, z)
    