'''Four prime sieves'''
'''by JDU'''

# --------------------------------------------------------------------------------------------

def prime_sieve_1_hsmw(n):
    """ Returns a list of all prime numbers from 2 to n (inclusive)"""
    """ Characteristics: A faster alternative to the sieve of eratosthenes """
    """ Warning: Your system might run out of memory and crash for inputs above 10^6 """
    
    if n < 2: return []

    # Initialize the sieve
    n1 = n+1
    sieve = [True] * n1
    sieve[0] = sieve[1] = False  # 0 and 1 are not prime
    sr = ceil(sqrt(n))
    
    # Mark non-primes divisible by 2 and 3
    for x in range(2, 4):
        for i in range(x*x, n1, x):
            sieve[i] = False
    
    # Mark non-primes starting from 5
    for i in range(5, sr):
        if sieve[i]:
            for k in range(i*i, n1, 2*i):  
                sieve[k] = False

    return [x for x in range(n1) if sieve[x]]

# --------------------------------------------------------------------------------------------

def prime_sieve_2_hsmw(n):
    """ Returns a list of all prime numbers from 2 to n (inclusive)"""
    """ Characteristics: A faster alternative to the sieve of eratosthenes """
    """ This is a slight variation to "sieve 1" with almost identical runtime """
    """ Warning: Your system might run out of memory and crash for inputs above 10^6 """
    
    if n < 2: return []

    # Initialize the sieve
    n1 = n+1
    sieve = [True] * n1
    sieve[0] = sieve[1] = False
    sr = ceil(sqrt(n))

    # Mark non-primes divisible by 2 and 3
    for x in range(2, 4):
        for i in range(x*x, n1, x):
            sieve[i] = False

    # Mark non-primes of residue class [5] mod 6, starting from 5
    for i in range(5, sr, 6):
        if sieve[i]:
            for k in range(i*i, n1, 2*i):  
                sieve[k] = False

    # Mark non-primes of residue class [1] mod 6, starting from 7
    for j in range(7, sr, 6):
        if sieve[j]:
            for k in range(j*j, n1, 2*j):
                sieve[k] = False

    return [x for x in range(n1) if sieve[x]]

# --------------------------------------------------------------------------------------------

def prime_sieve_3_hsmw(n,s1=True,s2=True):
    """ Returns two complementary lists for prime numbers from 5 to n (inclusive) """
    """ Characteristics: Even faster than the other ones (sieve 1&2) due to the output format """
    """ s1,s2: relates to the two sectors in which primes >= 5 can appear, more details in the code """
    """        default values are "True", set one of them to "False" to turn one sector off """
    """ Warning: Your system might run out of memory and crash for inputs above 10^6 """
    
    if n < 5: return

    # Initialize the sieve
    n1 = n+1
    sieve = [True] * (n1)
    sr = ceil(sqrt(n))

    # the following for-loop is not necessary for this sieve, but it does improve the total runtime
    # Mark non-primes divisible by 2 and 3
    for x in range(2, 4):
        for i in range(x*x, n1, x):
            sieve[i] = False
    
    # Mark non-primes starting from 5
    for i in range(5, sr):
        if sieve[i]:
            for k in range(i*i, n1, 2*i):  
                sieve[k] = False

    # any prime >= 5 can only be in one of the following sectors
    sector1 = range(5, n1, 6)  # elements of residue class [5] mod 6, starting from 5
    sector2 = range(7, n1, 6)  # elements of residue class [1] mod 6, starting from 7

    if not s1: 
        return [x for x in sector2 if sieve[x]]
    if not s2: 
        return [x for x in sector1 if sieve[x]]
        
    return [x for x in sector1 if sieve[x]], [x for x in sector2 if sieve[x]]

# --------------------------------------------------------------------------------------------

def prime_sieve_4_hsmw(n):
    """ Returns a list of all prime numbers from 2 to n (inclusive) """
    """ Characteristics: Multiple times slower than the other ones (sieve 1-3), """ 
    """                  but operates without setting any sieve element twice to "False" """
    """ Warning: Your system might run out of memory and crash for inputs above 10^6 """
    
    if n < 2: return []

    # Initialize the sieve
    n1 = n+1
    sieve = [True] * n1
    sieve[0] = sieve[1] = False  # 0 and 1 are not prime
    sr = ceil(sqrt(n))

    # Mark non-primes divisible by 2 and 3
    for x in range(2, 4):
        for i in range(x*x, n1, x):
            sieve[i] = False
       
    for i in range(5, sr, 6):
        if sieve[i]:
            val = i * i
            toggle = True # Start with 2*i
            while val < n1:
                sieve[val] = False
                if toggle: 
                    val += 2 * i
                else: 
                    val += 4 * i
                toggle = not toggle  # Alternate between 2*i and 4*i

    for j in range(7, sr, 6):
        if sieve[j]:
            val = j * j
            toggle = False  # Start with 4*j
            while val < n1:
                sieve[val] = False
                if toggle:
                    val += 2 * j
                else:
                    val += 4 * j
                toggle = not toggle  # Alternate between 4*j and 2*j 

    return [x for x in range(n1) if sieve[x]]

# --------------------------------------------------------------------------------------------