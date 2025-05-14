''' by JDU '''

# ----------------------------------------------------------------------------------------------

from sage.all import *      # See Cryptography/README.md 

# ----------------------------------------------------------------------------------------------

def pollards_rho_DLP(h, g):
    """ Returns the discrete logarithm of h to base g (mod p). """
    """ h and g must be elements of Zp* (with prime p), and be given in the form of mod(p)-objects. """
    """ g must be a generator of the group Zp*. """

    if not h or not g:
        raise ValueError("Inputs (h) and (g) must be non-zero!")
        
    if not hasattr(h, "modulus") or not hasattr(g, "modulus"):
        raise TypeError("Inputs (h) and (g) must be mod-objects!")

    if not h.modulus() == g.modulus():
        raise ValueError("Inputs (h) and (g) must have the same modulus!")

    p = g.modulus()
    n = g.multiplicative_order()
    
    if n != p-1:
        raise ValueError("g must be a generator of Zp*, with prime p!")

    def f(x, a, b):
        r = int(x) % 3
        if r == 1:
            return (h*x, (a+1) %n, b)
        elif r == 2:
            return (x**2, (2*a) %n, (2*b) %n)
        elif r == 0:
            return (g*x, a, (b+1) %n)

    xi, ai, bi = x2i, a2i, b2i = 1, 0, 0

    while True:
        xi, ai, bi = f(xi, ai, bi)
        x2i, a2i, b2i = f(*f(x2i, a2i, b2i))       
        if xi == x2i:
            s = (ai - a2i) %n
            t = (b2i - bi) %n
            d = gcd(s, n)
            t_d = t//d 
            n_d = n//d 
            inv = inverse_mod(s//d, n_d)
            for k in range(d):
                x = (t_d*inv + k*n_d) %n
                if g**x == h:
                    return x

# ----------------------------------------------------------------------------------------------
