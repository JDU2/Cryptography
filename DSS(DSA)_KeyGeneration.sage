''' DSS(DSA) key generation '''
''' by JDU '''

# ----------------------------------------------------------------------------------------------

from sage.all import *      # See Cryptography/README.md 

# ----------------------------------------------------------------------------------------------

def DSS_DSA_generate_key_pair(c):
    """ Returns a random DSS(DSA) key pair based on the input value (c), a menu from 1-4, which determines """
    """ the security standard of the key pair, by specifying the bit size of the prime modulus p and prime q. """

    if not isinstance(c, (int, Integer)):
        raise TypeError("Input (c) must be an integer!")   
    if not 1 <= c <= 4:
        raise ValueError("Input (c) must be within 1 <= c <= 4")

    if c == 1:
        bs_p=1024
        bs_q=160
    if c == 2:       
        bs_p=2048
        bs_q=224
    if c == 3:    
        bs_p=2048
        bs_q=256
    if c == 4:    
        bs_p=3072
        bs_q=256

    # q must have exactly bs_q bits: 2^(bs_q -1) < q < 2^(bs_q)
    q = random_prime((2**bs_q), False, 2**(bs_q-1)+1) # upper bound (1st arg!) not inclusive
    
    while True:
        k = ZZ.random_element(2**(bs_p - bs_q), 2**(bs_p - (bs_q-1))) # upper bound (2nd arg) not inclusive
        p = k * q + 1 # plus 1 makes it an odd number
        if (p-1) % q == 0:
            if p.nbits() == bs_p:
                if is_pseudoprime(p): # flag=0 (default): performs the Baillie-PSW test, a strong combination of Miller-Rabin and Lucas test
                    break
                    
    while True:
        h = ZZ.random_element(2, p-1) # upper bound (2nd arg) not inclusive
        if power_mod(h, (p-1)//q, p) != 1:
            g = power_mod(h, (p-1)//q, p)
            break
            
    X = ZZ.random_element(2, q-1) # upper bound (2nd arg) not inclusive
    Y = power_mod(g, X, p)
    
    return (p, q, g, X), (p, q, g, Y) # (private key), (public key)

# ----------------------------------------------------------------------------------------------
