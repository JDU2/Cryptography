''' DSS(DSA) Sinature creation and verification '''
''' by JDU '''

# ----------------------------------------------------------------------------------------------

from sage.all import *      # See Cryptography/README.md 

# ----------------------------------------------------------------------------------------------

def DSS_sign(m, prvkey):
    """ """

    (p, q, g, X) = prvkey
    
    if not isinstance(m, (int, Integer)) or not all(isinstance(i, (int, Integer)) for i in prvkey):
        raise TypeError("Inputs (m, prvkey) must consist of integers!") 

    k = ZZ.random_element(2, q-1)
    r = power_mod(g, k, p) % q
    s = inverse_mod(k, q) * (m + X*r) % q
    sig = (r, s)
    
    return sig

# ----------------------------------------------------------------------------------------------

def DSS_verify(m, sig, pubkey):
    """ """
    
    (p, q, g, Y) = pubkey
    (r, s) = sig

    if not isinstance(m, (int, Integer)) or not all(isinstance(i, (int, Integer)) for i in sig) or not all(isinstance(i, (int, Integer)) for i in pubkey):
        raise TypeError("Inputs (m, sig, pubkey) must consist of integers!")
    if not (0 < r < q): 
        raise ValueError("Signature component (r) must be within valid bounds! (0 < r < q)")
    if not (0 < s < q):
        raise ValueError("Signature component (s) must be within valid bounds! (0 < s < q)")

    t = inverse_mod(s, q) * m % q
    u = inverse_mod(s, q) * r % q
    v = power_mod(g, t, p) * power_mod(Y, u, p) % p % q
    
    return v == r

# ----------------------------------------------------------------------------------------------

