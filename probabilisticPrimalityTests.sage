''' Probabilistic primality tests - Fermat's and Miller-Rabin '''
''' by JDU '''

# ----------------------------------------------------------------------------------------------------

from sage.all import *      # See Cryptography/README.md 

# ----------------------------------------------------------------------------------------------------

def fermatsPrimalityTest(n, itr = 50):
    """ Tests whether an integer (n) is probably prime, for (itr) iterations. """
    """ The probabilitiy of this function to return a false "true" is approximately 1/(2^itr), with disregard to carmichael numbers."""
    """ Note: If (n) is a so-called carmichael number, this test will always return "true" despite (n) being composite. """
    """       Therefore you should use fermat's primality test in combination with a deterministic primality test, that tests the """
    """       range of potential carmichael factors from 3 to n^(1/3), or choose a more robust probabilistic test (Miller-Rabin). """

    if not isinstance(n, (int, Integer)) or not isinstance(itr, (int, Integer)):
        raise TypeError("Inputs (n, itr) must be integers!")

    if itr < 1: raise ValueError("Input (itr) must be >= 1")

    if n < 2:
        return False
    
    if n == 2:
        # Out of test range
        return True
    
    for _ in range(itr):
        a = 0
        # Assure coprimeness
        while gcd(a, n) != 1:
            a = ZZ.random_element(2, n)

        # By fermat's little theorem ...
        if power_mod(a, n-1, n) != 1:
            # n is definitely composite
            return False  
    
    # n is probably prime
    return True  

# ----------------------------------------------------------------------------------------------------

def millerRabinPrimalityTest(n, itr = 50):
    """ Tests whether an integer (n) is probably prime, for (itr) iterations. """
    """ The probabilitiy of this function to return a false "true" is approximately 1/(4^itr)."""

    if not isinstance(n, (int, Integer)) or not isinstance(itr, (int, Integer)):
        raise TypeError("Inputs (n, itr) must be integers!")

    if itr < 1: raise ValueError("Input (itr) must be >= 1")

    if n < 2:
        return False

    # Write n-1 as 2^r * d
    r, d = 0, n - 1

    while d % 2 == 0:
        r += 1
        d //= 2

    # Witness loop
    for _ in range(itr):
        
        a = ZZ.random_element(2, n-1)
        x = power_mod(a, d, n)

        if x == 1 or x == n-1:
            # Probably prime in this round
            continue

        for _ in range(r-1):
            # Square x up to r-1 times
            x = power_mod(x, 2, n)
            if x == n-1:
                # Probably prime in this round
                break
        else:
            # Definitely composite
            return False

    # Probably prime    
    return True

# ----------------------------------------------------------------------------------------------------
