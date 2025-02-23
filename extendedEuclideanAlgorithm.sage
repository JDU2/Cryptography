''' Extended Euclidean Algorithm (XGCD) '''
''' by JDU '''

# -------------------------------------------------------------------------------------------------------------------

from sage.all import *    # See Cryptography/README.md

# -------------------------------------------------------------------------------------------------------------------

def xgcd(a, b):
    """ Returns the greatest common divisor of two integers (a) and (b), and the Bézout coefficients a' and b' """
    """ Bézout's identity formula:  a' * a + b' * b = gcd(a, b) """
    """                             where:  a' is the modular inverse of a (mod b) """
    """                             and:    b' is the modular inverse of b (mod a) """

    if not isinstance(a, (int, Integer)):
        raise TypeError("Input (a) must be an integer!")

    if not isinstance(b, (int, Integer)):
        raise TypeError("Input (b) must be an integer!")
    
    if b == 0:
        # End point of the call chain
        # Start of recursive upwards propagation from here
        return abs(a), sgn(a), 0
        
    q = a // b  # Quotient
    r = a % b   # Remainder
    
    # Recursive downwards propagation until b == 0
    # Then recursive variable assignment from upwards propagation
    d, b_star, r_star = xgcd(b, r)
    
    # Recursive updating of values
    # (Renaning of return values is optional)
    gcd_ab = d
    a_prime = r_star
    b_prime = b_star - q * r_star
    
    # Recursive upwards propagation until final return
    return gcd_ab, a_prime, b_prime

# -------------------------------------------------------------------------------------------------------------------
