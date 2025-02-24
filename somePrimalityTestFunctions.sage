''' 2 Primality Test Functions '''
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

    if not isinstance(n, (int, Integer)) or not isinstance(itr, (int, Integer)):
        raise TypeError("Inputs (n, itr) must be integers!")

    if n < 2: raise ValueError("Input (n) must be >= 2")

    if itr < 1: raise ValueError("Input (itr) must be >= 1")

    if n <= 3: return True

    # Test prime factors 2 and 3    
    if (n % 2 == 0) or (n % 3 == 0): 
        return False
     
    itrLimit = itr      
    sr = floor(sqrt(n))
    max_itr = ((sr-1)//6)+1
    
    # Adjustments for least possible iterations
    if (sr % 6) != 5: max_itr -= 1  
    if (sr % 6) != 1: sr += 1

    # Test all other potential prime factors
    for i in range(5, sr, 6):
        if (n % i == 0) or (n % (i+2) == 0):
            return False
        itr -= 1
        if itr == 0:
            if itrLimit == max_itr:
                # Special case
                break
            certainty = itrLimit*100/float(max_itr)
            print(f"Result: Input (n) is possibly prime! ... Certainty: {certainty:.18f} % ... Increase (itr) for more certainty!\n")
            print(f"Note: This method requires {max_itr} iterations for 100% certainty, if (n) is prime.\n")
            return None
 
    return True

# ------------------------------------------------------------------------------------------------------------

def primalityTest_2(n, itr):
    """ Tests whether an integer (n) is a prime number, for (itr) iterations. """
    """ Characteristics: This is an alternative approach which is at least in theory more efficient, as it only tests  """
    """                  real prime numbers as potential factors of (n), in contrast to primalityTest_1(), which also  """
    """                  tests certain composites besides primes. However, the reliance on prime counting tools makes  """
    """                  it overall slower than primalityTest_1(), despite using less iterations. """
    """ Note: If the given integer for iterations (itr) is not sufficient, it prints out the certainty level of """
    """       primality (in %) and the max amount of iterations required to determine guaranteed primality for (n). """
    """ Warning: This method gets very slow for (itr) above 10^6. """

    if not isinstance(n, (int, Integer)) or not isinstance(itr, (int, Integer)):
        raise TypeError("Inputs (n, itr) must be integers!")

    if n < 2: raise ValueError("Input (n) must be >= 2")

    if itr < 1: raise ValueError("Input (itr) must be >= 1")

    itrLimit = itr      
    sr = floor(sqrt(n))

    if sr >= 2:
        # Max number of primes to test
        max_nop = prime_pi(previous_prime(sr+1))
    else:
        max_nop = 1

    # Test all potential prime factors
    for i in primes(sr+1):
        if n % i == 0: 
            return False
        itr -= 1
        if itr == 0:
            if itrLimit == max_nop:
                # Special case
                break
            certainty = (itrLimit-itr)*100/float(max_nop)
            print(f"Result: Input (n) is possibly prime! ... Certainty: {certainty:.18f} % ... Increase (itr) for more certainty!")
            print(f"Note: This method requires {max_nop} iterations for 100% certainty, if (n) is prime.\n")
            return None
 
    return True

# ------------------------------------------------------------------------------------------------------------
