''' Primality Test Functions '''
''' by JDU '''

# ------------------------------------------------------------------------------------------------------------

from sage.all import *    # See Cryptography/README.md

# ------------------------------------------------------------------------------------------------------------

def primalityTest_1(n, itr):
    """ Tests whether an integer (n) is a prime number, for (itr) iterations. """
    """ Characteristics: A simple approach without the reliance on prime number tools. """
    """ Note: If the given integer for iterations (itr) is not sufficient, it prints out the certainty level of """
    """       primality (in %) and the max amount of iterations required to determine guaranteed primality for (n). """
    """ Warning: This method gets very slow for (itr) above 10^7. """

    # Test prime factors 2 and 3
    if (not n % 2) or (not n % 3): 
        return False

    itrLimit = itr      
    sr = floor(sqrt(n))
    max_itr = ((sr-1)//6)+1

    # Test all other potential prime factors
    for i in range(7, sr, 6):
        if (not n % i) or (not n % (i-2)): 
            return False
        itr -= 1
        if not itr:
            certainty = (1+itrLimit*6)*100/float(sr)
            print(f"Result: Input (n) is possibly prime! ... Certainty: {certainty:.18f} % ... Increase (itr) for more certainty!")
            print(f"Note: This method requires {max_itr} iterations for 100% certainty, if (n) is prime.\n")
            return
 
    return True

# ------------------------------------------------------------------------------------------------------------
