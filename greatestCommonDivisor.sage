''' 2 Implementations of the Euclidean Algorithm '''
''' by JDU '''

# ----------------------------------------------------------------------------------------------

from sage.all import *    # See Cryptography/README.md

# ----------------------------------------------------------------------------------------------

def gcd_rec(a, b):
    """ Returns the greatest common divisor (GCD) of two integers (a) and (b). """
    """ This is the recursive version of the euclidean algorithm. """

    if not isinstance(a, (int, Integer)):
        raise TypeError("Input (a) must be an integer!")

    if not isinstance(b, (int, Integer)):
        raise TypeError("Input (b) must be an integer!")

    if b == 0:
        return abs(a)
    
    return gcd_rec(b, a % b)

# ----------------------------------------------------------------------------------------------

def gcd_itr(a, b):
    """ Returns the greatest common divisor (GCD) of two integers (a) and (b). """
    """ This is the iterative version of the euclidean algorithm. """

    if not isinstance(a, (int, Integer)):
        raise TypeError("Input (a) must be an integer!")

    if not isinstance(b, (int, Integer)):
        raise TypeError("Input (b) must be an integer!")

    while(b):
        a, b = b, a % b

    return abs(a)

# ----------------------------------------------------------------------------------------------
