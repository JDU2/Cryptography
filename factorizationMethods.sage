''' Factorization methods - Based on Trial division, Fermat's method, and Dixon's method '''
''' by JDU '''

# ----------------------------------------------------------------------------------------------------

from sage.all import *    # See Cryptography/README.md

# ------------------------------------------------------------------------------------------------------------

def trialDivisionFactorization(n, itr):
    """ Searches for a non-trivial factor of a composite integer (n) by trial division, for up to (itr) iterations. """
    """ Characteristics: A simple optimized approach without the reliance on prime counting tools. Uses a step size of 6. """
    """                  If a factor is found it returns the factorization based on that factor, otherwise "none". """
    """ Note: If the given number of iterations (itr) is not sufficient, it prints out the the max amount of iterations """
    """       required to find a factor of (n) for the theoretical worst case, where (n) would be a perfect square of a """
    """       prime number. Additionally, it prints out the coverage level of worst case iterations (in %). """
    """ Warning: This method becomes slow for (itr) above 10^8 (required for prime factors > 6 x 10^8). """

    if not isinstance(n, (int, Integer)) or not isinstance(itr, (int, Integer)):
        raise TypeError("Inputs (n, itr) must be integers!")

    if n < 4:       raise ValueError("Input (n) must be >= 4")
    if is_prime(n): raise ValueError("Input (n) must be a composite number!")
    if itr < 1:     raise ValueError("Input (itr) must be >= 1")

    # Test prime factors 2 and 3    
    if (n % 2 == 0): 
        return 2, n//2
    if (n % 3 == 0): 
        return 3, n//3
        
    itrLimit = itr      
    sr = floor(sqrt(n))
    wc_itr = ((sr-1)//6)+1
    
    # Adjustments for special properties
    if (sr % 6) != 5: wc_itr -= 1  
    if (sr % 6) != 1: sr += 1

    # Test all other potential prime factors
    for i in range(5, sr, 6):
        if (n % i == 0):
            return i, n//i
        if (n % (i+2) == 0):
            return i+2, n//(i+2)
        itr -= 1
        if itr == 0:
            wc_coverage = itrLimit*100/float(wc_itr)
            print(f"Result: No prime factors found within {itrLimit} iterations! \n")
            print(f"This method requires {wc_itr} iterations for the worst case (=> if (n) is a perfect square of a prime number). \n")
            print(f"Coverage of worst case iteration amount: {wc_coverage:.18f} % \n")
            return None

# ----------------------------------------------------------------------------------------------------

def fermatsFactorization(n, itr = 2*10**6):
    """ Returns two distinct non-trivial factors (a) and (b) of a positive odd integer (n), """
    """ after performing at most (itr) iterations, otherwise returns "false". """
    """ Fermat's formula: n = ab = (x+y)(x-y) = x²-y², where n is odd and a > b > 0. """       
    """ Characteristics: Optimized step size by sorting out prime factor 3. """
    """                  If you prefer a simpler implementation, then just sort out prime factor 2 """
    """                  (standard) and use a step size of 2, which does not rely on mod 6 analysis. """
    """ Note: Fermat's method works best if the difference of (a) and (b) is expected to be small. """
    """       If no factors were found within (itr) iterations, it prints the maximum number """
    """       of iterations required for a definite result (=> factorization or primeness). """
    """ Warning: This method becomes slow for performed (itr) above 2 x 10^6. """

    if not isinstance(n, (int, Integer)) or not isinstance(itr, (int, Integer)):
        raise TypeError("Inputs (n, itr) must be integers!") 
    if n < 0:
        raise ValueError("Input (n) must be positive!")
    if n % 2 == 0:
        raise ValueError("Input (n) must be odd!")
    if itr < 1: 
        raise ValueError("Input (itr) must be >= 1")

    # Sorting out prime factor 3
    if n % 3 == 0: 
        if n == 3: 
            print(f"Input ({n}) is a prime number!")
            return False
        return n//3, 3

    # Sorting out squares (a must be greater than b!)
    if is_square(n): 
        print(f"Input ({n}) is a square value!")
        return False

    # Initial values for x and y
    x = floor(sqrt(n))+1
    y = 0

    # Worst case: (b) = 5 (=> smallest prime factor after 2 and 3)
    max_y = ((n/5) - 5)/2       
    max_x = sqrt(n + max_y**2)

    # Now that we sorted out prime factors 2 and 3 ...
    # (n) can only be in residue classes [1],[5] mod 6
    if n % 6 == 1:
        # (x) can only be in residue classes [1],[2],[4],[5] mod 6
        while True:
            # Adjust start value for (x)
            if x%6 == 1 or x%6 == 4:
                x_1st_step, x_2nd_step = 1, 2
                break
            elif x%6 == 2 or x%6 == 5:
                x_1st_step, x_2nd_step = 2, 1
                break
            x += 1

    # Remaining case
    elif n % 6 == 5:
        # (x) can only be in residue classes [0],[3] mod 6
        while True:
            # Adjust start value for (x)
            if x%6 == 0 or x%6 == 3:
                x_step = 3
                break
            x += 1

    # Worst case iterations for compositeness (or primeness) of (n)
    max_itr = floor((max_x - x)/3)+1

    # Now search for a perfect square of (y) by stepping up (x) (=> x²-n = y²)
    if n % 6 == 1:
        for _ in range(itr):
            if is_square(x**2 - n):
                y = sqrt(x**2 - n)
                break
            x += x_1st_step
            if is_square(x**2 - n):
                y = sqrt(x**2 - n)
                break
            x += x_2nd_step

    elif n % 6 == 5:
        for _ in range(itr):
            if is_square(x**2 - n):
                y = sqrt(x**2 - n)
                break
            x += x_step

    # Failure: No integer square value of (y) has been found 
    if not y:
        if itr < max_itr:
            print(f"No non-trivial factors (a),(b) found within {itr} iterations!")
            print(f"{max_itr} iterations required (at worst case) for a definite result.")
        else: 
            print(f"Input ({n}) is a prime number!")
        return False

    # Factorize (n)
    a = (x+y)
    b = (x-y)

    # Success
    return a, b
# ----------------------------------------------------------------------------------------------------

def dixonsFactorization(n, B, B_fn = False):
    """ Returns two non-trivial factors of integer (n), based on a random search with a given smoothness bound (B or B_fn), otherwise returns "None". """
    """ From Fermat:     n = (x+y)(x-y) = x²-y² """    
    """ Dixons's Idea:   n = (x+y)(x-y) = x²-y² = 0 (mod n) """
    """                  Therefore, x² = y² (mod n) """
    """                  Now, find random values x² (mod n) that are B-smooth and collect the exponent vectors (mod 2) of their prime factors over the whole factor base determined by (B). """
    """                  Solve the system of linear equations to find combinations of exponent vectors (by vector addition) that result in a zero vector (mod 2). """
    """                  Then, take each of these linearly dependent combinations of exponent vectors and construct a congruence of squares => x² = y² (mod n), """
    """                  where x is constructed by multiplying their corresponding x values together (mod n) and y² is constructed by multiplying """
    """                  their corresponding y² component values (=> x² mod n) together, then taken the square root of the product (mod n) to get y. """
    """                  Then, check whether gcd(x-y,n) or gcd(x+y,n) results in a non-trivial factor of (n). """
    """ Characteristics: Implementation based on the suggestions from the book "applied cryptanalysis" by Stamp and Low. """
    """                  This implementation uses (-1) as an additional entry in the factor base and searches for modular """
    """                  numbers between (-n/2) and (n/2), which allows finding more B-smooth numbers within the smoothness bound. """
    """                  (B) can be chosen manually or from a menu of functions that are dependent on (n). """
    """                  To selcect such function choose (B_fn) between [1,2,3] and set (B) to zero or "False". """
    """ Note: The success in finding a non-trivial factor of a composite integer (n) depends on """
    """       the amount of relations (=> B-smooth numbers x² mod n) and on the size of (B). """
    """       This implementation collects len(factor_base)+1 relations, which is the minimum amount """
    """       that guarantees to yield at least one linear dependency (=> zero vector mod 2). """
    """       However, in some cases all the linear dependencies are based on x and y values, """
    """       where x+y = n and/or x = y, and both of these conditions can (but don't """
    """       always have to) result in finding a trivial factor of (n), which returns "None". """

    if not isinstance(n, (int, Integer)) or not isinstance(B, (int, Integer)):
        raise TypeError("Inputs (n, B) must be integers!")  
    if n < 4: 
        raise ValueError("Input (n) must be >= 4")   
    if not B:
        # two good lower bounds for B
        if B_fn == 1: 
            B = isqrt(sqrt(n)) # B = n^(1/4), good for n < 10^11
        elif B_fn == 2: 
            B = int(log(n)**2) # B = ln(n)², good for n > 10^12
        # good balanced upper bound for B
        elif B_fn == 3: 
            B = ceil(exp(sqrt(log(n)*log(log(n))/2))) # B = e^((ln(n)*ln(ln(n))/2)^(1/2))
        else: 
            raise ValueError("Input (B_fn) must be in [1,2,3]")       
    if B < 2: 
        raise ValueError("Input (B) must be >= 2")
    
    factor_base = [-1]+[p for p in prime_range(B+1)]
    x_components_candidates, y2_components_candidates, relations = [], [], []
    sr = floor(sqrt(n))+1

    while len(relations) < len(factor_base)+1: # this guarantees to yield at least one linear dependency

        while True:  
            # find numbers x such that (x^2 % n) is B-smooth (=> relation)
            x = randint(sr, n-1)
            if x in x_components_candidates:
                continue
            x2_mod_n = power_mod(x, 2, n)
            if not x2_mod_n:
                # 0 can't get factored!
                # but we can directly check the relation
                # note: y = sqrt(x2_mod_n)%n = 0
                d = gcd(x, n)
                if 1 < d < n:
                    return d, n//d # factorization succeeded
            if x2_mod_n > n/2:
                # reduce absolute value by modular shift into negative space
                x2_mod_n -= n
                        
            # start trial division over factor base (excluding -1)
            exponents = [0] * len(factor_base)
            a = x2_mod_n
            for i, p in enumerate(factor_base[1::],1):
                while a % p == 0:
                    a //= p
                    exponents[i] += 1
                if abs(a) == 1:
                    # new B-smooth relation found!
                    exponent_vector_mod_2 = [e % 2 for e in exponents]
                    break
            if abs(a) > 1:
                # relation is not B-smooth
                continue
            if x2_mod_n < 0:
                exponent_vector_mod_2[0] = 1
            elif sum(exponent_vector_mod_2) == 0: 
                # check single relation
                y = isqrt(x2_mod_n) % n
                # find non-trivial factors of n
                d = gcd(x - y, n)
                if 1 < d < n:
                    return d, n//d # factorization succeeded
                d = gcd(x + y, n)
                if 1 < d < n:
                    return d, n//d # factorization succeeded
            # exit inner while-loop
            break
        
        x_components_candidates.append(x)
        y2_components_candidates.append(x2_mod_n)
        # let exponent vectors mod 2 represent the relations
        relations.append(exponent_vector_mod_2)

    M = Matrix(GF(2), len(factor_base), len(relations))
    for i, exp_v in enumerate(relations):
        M.set_column(i, exp_v)

    # solve system of linear equations
    null_space = M.right_kernel()
    for v in null_space.basis():
        # construct x for congruence of squares x²=y² (mod n)
        x = prod(x_components_candidates[i] for i in range(len(v)) if v[i] == 1) % n
        # construct y² for congruence of squares x²=y² (mod n)
        y2 = prod(y2_components_candidates[i] for i in range(len(v)) if v[i] == 1)
        # get y
        y = isqrt(y2) % n
        # find non-trivial factors of n
        d = gcd(x - y, n)
        if 1 < d < n:
            return d, n//d # factorization succeeded
        d = gcd(x + y, n)
        if 1 < d < n:
            return d, n//d # factorization succeeded
            
    return None # factorization failed

# ----------------------------------------------------------------------------------------------------
