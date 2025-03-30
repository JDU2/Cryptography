''' These are the sieves from https://github.com/JDU2/Cryptography/blob/main/somePrimeSieves.sage, '''
''' that are now using bitsets (from sagemath) instead of (bit-)lists, which use less memory '''
''' by JDU '''

# ---------------------------------------------------------------------------------------------------------

from sage.all import *            # See Cryptography/README.md

# ---------------------------------------------------------------------------------------------------------

def prime_sieve_1(n):
    """ Returns a list of all prime numbers from 2 up to n (inclusive). """
    """ Characteristics: A faster alternative to the optimized sieve of eratosthenes. """
    """ Warning: Your execution environment might run out of memory and crash for inputs above 10^9. """

    if not isinstance(n, (int, Integer)):
        raise TypeError("Input (n) must be an integer!")
    
    if n < 2: return []

    # Initialize the sieve
    n1 = n+1
    sieve = Bitset(n1 * '1')              # One extra bit (n+1) such that index 0 does represent integer 0
    sieve.discard(0); sieve.discard(1)    # 0 and 1 are not prime
    sr = floor(sqrt(n))
    
    # Mark non-primes divisible by 2 
    # These are the elements of residue classes [0],[2],[4] mod 6, starting from 4
    sieve.difference_update(Bitset(range(4, n1, 2)))
    
    # Mark non-primes divisible by 3 (but not by 2)
    # These are the elements of residue class [3] mod 6, starting from 9
    sieve.difference_update(Bitset(range(9, n1, 6)))

    # Mark non-primes of residue class [5] mod 6, starting from 5
    for i in range(5, sr+1, 6):
        if i in sieve:
            sieve.difference_update(Bitset(range(i**2, n1, 2*i)))

    # Mark non-primes of residue class [1] mod 6, starting from 7
    for j in range(7, sr+1, 6):
        if j in sieve:
            sieve.difference_update(Bitset(range(j**2, n1, 2*j)))
    
    return list(sieve)

# ----------------------------------------------------------------------------------------------

def prime_sieve_2(n, s1=True, s2=True):
    """ Returns a list of all prime numbers from 5 up to n (inclusive), which is composed out of two subsets (residue class [1] and [5] mod 6). """
    """ Characteristics: Faster than sieve 1 as it does not consider primes 2, 3, and their multiples. """
    """ s1, s2: Relates to the two sectors in which primes >= 5 can appear (mod 6). """
    """         Set one of them to "false" if you want only primes from one sector. """
    """ Warning: Your execution environment might run out of memory and crash for inputs above 10^9. """

    if not isinstance(n, (int, Integer)):
        raise TypeError("Input (n) must be an integer!")

    if not isinstance(s1, bool) or not isinstance(s2, bool):
        raise TypeError("Inputs (s1, s2) must be booleans!")

    if n < 5 or (not s1 and not s2): return []

    # Initialize the sieve
    n1 = n+1
    sieve = Bitset(n1 * '1')    # One extra bit (n+1) such that index 0 does represent integer 0
    sr = floor(sqrt(n))

    # Any prime number >= 5 can only be in one of the following sectors
    sector1 = range(5, n1, 6)     # Elements of residue class [5] mod 6, starting from 5
    sector2 = range(7, n1, 6)     # Elements of residue class [1] mod 6, starting from 7
    
    if s1:
        # Mark non-primes of residue class [5] mod 6, starting from 5
        for i in range(5, sr+1, 6):
            if i in sieve:
                sieve.difference_update(Bitset(range(i**2, n1, 2*i)))
    
    if not s2: 
        # Returns only primes from sector 1
        return [p for p in sector1 if p in sieve]

    if s2:
        # Mark non-primes of residue class [1] mod 6, starting from 7
        for j in range(7, sr+1, 6):
            if j in sieve:
                sieve.difference_update(Bitset(range(j**2, n1, 2*j)))

    if not s1:
        # Returns only primes from sector 2
        return [p for p in sector2 if p in sieve]    

    # Returns the primes from both sectors
    return sorted([p for p in sector1 if p in sieve] + [p for p in sector2 if p in sieve])

# ----------------------------------------------------------------------------------------------

def prime_sieve_3(n):
    """ Returns a list of all prime numbers from 2 up to n (inclusive). """
    """ Characteristics: Slightly slower than sieve 1, """ 
    """                  but operates without setting any sieve element more than once to "false", """
    """                  by using a pattern that alternates the distance to the next non-primes. """
    """ Warning: Your execution environment might run out of memory and crash for inputs above 10^9. """

    if not isinstance(n, (int, Integer)):
        raise TypeError("Input (n) must be an integer!")

    if n < 2: return []

    # Initialize the sieve
    n1 = n+1
    sieve = Bitset(n1 * '1')              # One extra bit (n+1) such that index 0 does represent integer 0
    sieve.discard(0); sieve.discard(1)    # 0 and 1 are not prime
    sr = floor(sqrt(n))
    
    # Mark non-primes divisible by 2 
    # These are the elements of residue classes [0],[2],[4] mod 6, starting from 4
    sieve.difference_update(Bitset(range(4, n1, 2)))
    
    # Mark non-primes divisible by 3 (but not by 2)
    # These are the elements of residue class [3] mod 6, starting from 9
    sieve.difference_update(Bitset(range(9, n1, 6)))

    # Mark non-primes of residue class [5] mod 6, starting from 5
    for i in range(5, sr+1, 6):
        if i in sieve:
            x = i**2
            small = 2*i
            big = 4*i
            while x < n1:    # Alternate distance between 2*i and 4*i
                sieve.discard(x) 
                x += small   # Start with 2*i
                if x < n1:
                    sieve.discard(x)
                x += big

    # Mark non-primes of residue class [1] mod 6, starting from 7
    for j in range(7, sr+1, 6):
        if j in sieve:
            x = j**2
            small = 2*j
            big = 4*j
            while x < n1:   # Alternate distance between 4*j and 2*j 
                sieve.discard(x) 
                x += big    # Start with 4*j
                if x < n1:
                    sieve.discard(x)
                x += small

    return list(sieve)

# ----------------------------------------------------------------------------------------------
