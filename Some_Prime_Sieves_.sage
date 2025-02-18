''' 3 Prime Sieves - For a better understanding of the distribution of primes '''
''' by JDU '''

# ----------------------------------------------------------------------------------------------

from sage.all import *      # For full sagemath functionality, you might need this import in 
                            # case your sagemath kernel did not properly initialize, which can 
                            # happen when working with jupyter notebook files within vs code

# ----------------------------------------------------------------------------------------------

def prime_sieve_1(n):
    """ Returns a list of all prime numbers from 2 up to n (inclusive). """
    """ Characteristics: A faster alternative to the optimized sieve of eratosthenes. """
    """ Warning: Your execution environment might run out of memory and crash for inputs above 10^6. """
    
    if n < 2: return []

    # Initialize the sieve
    n1 = n+1
    sieve = [True] * n1            # One extra element (n+1) such that index 0 does represent integer 0
    sieve[0] = sieve[1] = False    # 0 and 1 are not prime
    sr = floor(sqrt(n))
    
    # Mark non-primes divisible by 2 
    # These are the elements of residue classes [0],[2],[4] mod 6, starting from 4
    for x in range(4, n1, 2):
        sieve[x] = False
    
    # Mark non-primes divisible by 3 (but not by 2)
    # These are the elements of residue class [3] mod 6, starting from 9
    for x in range(9, n1, 6):
        sieve[x] = False

    # Mark non-primes of residue class [5] mod 6, starting from 5
    for i in range(5, sr+1, 6):
        if sieve[i]:
            for x in range(i**2, n1, 2*i):  
                sieve[x] = False

    # Mark non-primes of residue class [1] mod 6, starting from 7
    for j in range(7, sr+1, 6):
        if sieve[j]:
            for x in range(j**2, n1, 2*j):
                sieve[x] = False

    return [p for p in range(n1) if sieve[p]]

# ----------------------------------------------------------------------------------------------

def prime_sieve_2(n, s1=True, s2=True):
    """ Returns two complementary lists for prime numbers from 5 up to n (inclusive). """
    """ Characteristics: Faster than sieve 1 as it does not consider primes 2, 3, and their multiples. """
    """ s1, s2: Relates to the two sectors in which primes >= 5 can appear; more details in the code. """
    """         Default values are "True". Set one of them to "False" to turn one sector off. """
    """ Warning: Your execution environment might run out of memory and crash for inputs above 10^6. """
    
    if n < 5 or (not s1 and not s2): return []

    # Initialize the sieve
    n1 = n+1
    sieve = [True] * n1     # One extra element (n+1) such that index 0 does represent integer 0
    sr = floor(sqrt(n))

    # Any prime number >= 5 can only be in one of the following sectors
    sector1 = range(5, n1, 6)    # elements of residue class [5] mod 6, starting from 5
    sector2 = range(7, n1, 6)    # elements of residue class [1] mod 6, starting from 7
    
    if s1:
        # Mark non-primes of residue class [5] mod 6, starting from 5
        for i in range(5, sr+1, 6):
            if sieve[i]:
                for x in range(i**2, n1, 2*i):  
                    sieve[x] = False
    
    if not s2: 
        return [p for p in sector1 if sieve[p]]    # Returns only primes from sector 1

    if s2:
        # Mark non-primes of residue class [1] mod 6, starting from 7
        for j in range(7, sr+1, 6):
            if sieve[j]:
                for x in range(j**2, n1, 2*j):
                    sieve[x] = False

    if not s1: 
        return [p for p in sector2 if sieve[p]]    # Returns only primes from sector 2

    # Returns the primes from both sectors
    return sorted([p for p in sector1 if sieve[p]] + [p for p in sector2 if sieve[p]])

# ----------------------------------------------------------------------------------------------

def prime_sieve_3(n):
    """ Returns a list of all prime numbers from 2 up to n (inclusive). """
    """ Characteristics: Slightly slower than sieve 1, """ 
    """                  but operates without setting any sieve element more than once to "false", """
    """                  by using a pattern that alternates the distance to the next non-primes. """
    """ Warning: Your execution environment might run out of memory and crash for inputs above 10^6. """
    
    if n < 2: return []

    # Initialize the sieve
    n1 = n+1
    sieve = [True] * n1            # One extra element (n+1) such that index 0 does represent integer 0
    sieve[0] = sieve[1] = False    # 0 and 1 are not prime
    sr = floor(sqrt(n))
    
    # Mark non-primes divisible by 2 
    # These are the elements of residue classes [0],[2],[4] mod 6, starting from 4
    for x in range(4, n1, 2):
        sieve[x] = False
    
    # Mark non-primes divisible by 3 (but not by 2)
    # These are the elements of residue class [3] mod 6, starting from 9
    for x in range(9, n1, 6):
        sieve[x] = False

    # Mark non-primes of residue class [5] mod 6, starting from 5
    for i in range(5, sr+1, 6):
        if sieve[i]:
            x = i**2
            small = 2*i
            big = 4*i
            while x < n1:   # Alternate distance between 2*i and 4*i
                sieve[x] = False 
                x += small   # Start with 2*i
                if x < n1:
                    sieve[x] = False
                x += big

    # Mark non-primes of residue class [1] mod 6, starting from 7
    for j in range(7, sr+1, 6):
        if sieve[j]:
            x = j**2
            small = 2*j
            big = 4*j
            while x < n1:   # Alternate distance between 4*j and 2*j 
                sieve[x] = False 
                x += big   # Start with 4*j
                if x < n1:
                    sieve[x] = False
                x += small

    return [p for p in range(n1) if sieve[p]]

# ----------------------------------------------------------------------------------------------
