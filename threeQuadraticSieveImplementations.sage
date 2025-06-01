''' 3 implementations of the Quadratic Sieve Factoring Algorithm, labeled 1,2,3 '''
''' by JDU '''

# ------------------------------------------------------------------------------------------------------------

from sage.all import *    # See Cryptography/README.md

# ------------------------------------------------------------------------------------------------------------

def QS_Factorization_1(n, I_sf, B, B_fn = False):
    """ Attempts to factor an integer (n) based on the Quadratic Sieve method with smoothness bound (B or B_fn). """
    """ The endpoints of the sieving interval [-I, I] are scaled by the factor (I_sf) applied on the smoothness bound. """
    """ Returns two non-trivial factors of n, or "None" if unsuccessful. """
    """ Characteristics: Uses "fast log base 2 approximations" for sieving. """
    """                  The smoothness bound can be chosen manually (B) or from a menu of functions of (n). """
    """                  To selcect such function choose (B_fn) between [1,2,3] and set (B) to zero or "False". """
    """ Recommendations: For efficiency reasons, start with a reasonably small choice for the smoothness bound. """
    """                  Start with a small scaling factor (around 0.5), then slowly increase it (up until 2.0  """
    """                  or 3.0) before further increasing the smoothness bound and resetting the scaling factor. """

    if not isinstance(n, (int, Integer)) or not isinstance(I_sf, (float, int, Integer)) or not isinstance(B, (int, Integer)):
        raise TypeError("Inputs (n, B) must be integers, (I_sf) must be float or integer") 
    if n < 4: 
        raise ValueError("Input (n) must be >= 4")
    if is_square(n): # important check; avoids infinite loops
        s = isqrt(n)
    if n % 2 == 0:
        return 2, n//2
    if I_sf <= 0:
        raise ValueError("Input (I_sf) must be > 0")
    if not B:
        # two good lower bounds for B
        if B_fn == 1: 
            B = isqrt(sqrt(n)) # B = n^(1/4), good for n <= 10^11
        elif B_fn == 2: 
            B = int(log(n)**2) # B = ln(n)², good for n >= 10^12
        # good balanced upper bound for B
        elif B_fn == 3: 
            B = ceil(exp(sqrt(log(n)*log(log(n))/2))) # B = e^((ln(n)*ln(ln(n))/2)^(1/2))
        else: 
            raise ValueError("Input (B_fn) must be in [1,2,3]")       
    if B < 2: 
        raise ValueError("Input (B) must be >= 2")

    sr = floor(sqrt(n))
    I = int(B*I_sf) # scaling of interval endpoints

    # the quadratic polynomial is given as Q(z)=(sr+z)²-n with z in interval [-I,I].
    # since we want to construct congruences of squares x²=y²(mod n), for which we will construct x and y²,
    # compute the "x components candidates" (= (sr+z) values) for the whole interval [-I,I]
    x_cc = [sr + z for z in range(-I, I+1)]
    # compute the "y² components candidates" (= Q(z) values) for the whole interval [-I,I]
    y2_cc = [x**2 -n for x in x_cc]
    
    factor_base = []
    sieve = [0] * len(y2_cc)

    def fast_prime_power(n, p):
        if p == 2:
            return 1 << (n & -n).bit_length() - 1
        p_e = 1
        while n % p == 0:
            n //= p
            p_e *= p
        return p_e

    # start sieving
    for p in prime_range(B+1):
        # Find the square roots of n modulo p, if they exist
        if p == 2:
            roots = [1] # Because n is odd
        else:
            ls = legendre_symbol(n, p)
            if ls == 0:
                # non-trivial factor of n found
                return p, n//p
            elif ls == 1:
                roots = mod(n, p).sqrt(all=true)
                factor_base.append(p)
            else:
                # no modular square roots exist (ls == -1)
                continue
        for r in roots:
            # Determine the start index for sieving
            start_index = (r - sr + I) %p
            for i in range(start_index, len(y2_cc), p):
                # add approximated logarithmic weight
                sieve[i] += fast_prime_power(y2_cc[i], p).bit_length()
    
    # calculate log(base 2) approximations (=> thresholds)
    y2_log = [abs(i).bit_length()-1 for i in y2_cc]
    
    x_cc_sieved = []
    y2_cc_sieved = []

    # evaluate the sieve elements
    for x, y2, s, l in zip(x_cc, y2_cc, sieve, y2_log):
        # compare the accumulated approximated logarithmic weights against the thresholds
        if s > l:
            # store the sieved candidates
            x_cc_sieved.append(x)
            y2_cc_sieved.append(y2)

    factor_base = [-1]+factor_base
    relations = []
    
    for s in y2_cc_sieved:
        f_dict = dict(factor(s)) # dictionary maps prime factors to their exponents
        if s < 0:
            f_dict[-1] = 1 # needs to be set manually!
        exponent_vector_mod_2 = [f_dict.get(p,0)%2 for p in factor_base]
        # collect new relation (represented by exponent vector mod 2)
        relations.append(exponent_vector_mod_2)
    
    M = Matrix(GF(2), zip(*relations)) # zip() transposes columns into rows

    # solve the linear algebra
    null_space = M.right_kernel() # solves M * v = 0 (mod 2) 
    for v in null_space.basis(): # .basis() converts the null space into a list of vectors v
        x = prod(x_cc_sieved[i] for i in range(len(v)) if v[i] == 1)%n
        y2 = prod(y2_cc_sieved[i] for i in range(len(v)) if v[i] == 1)
        y = isqrt(y2)%n
        d = gcd(n, x-y)
        if 1 < d < n:
            return d, n//d  # Factorization succeeded
        d = gcd(n, x+y)
        if 1 < d < n:
            return d, n//d  # Factorization succeeded
            
    return None  # Factorization failed

# ------------------------------------------------------------------------------------------------------------

def QS_Factorization_2(n, I_sf, B, B_fn = False):
    """ Attempts to factor an integer (n) based on the Quadratic Sieve method with smoothness bound (B or B_fn). """
    """ The endpoints of the sieving interval [-I, I] are scaled by the factor (I_sf) applied on the smoothness bound. """
    """ Returns two non-trivial factors of n, or "None" if unsuccessful. """
    """ Characteristics: Uses "reproduction by multiplication" for sieving. """
    """                  The smoothness bound can be chosen manually (B) or from a menu of functions of (n). """
    """                  To selcect such function choose (B_fn) between [1,2,3] and set (B) to zero or "False". """
    """ Recommendations: For efficiency reasons, start with a reasonably small choice for the smoothness bound. """
    """                  Start with a small scaling factor (around 0.5), then slowly increase it (up until 2.0  """
    """                  or 3.0) before further increasing the smoothness bound and resetting the scaling factor. """

    if not isinstance(n, (int, Integer)) or not isinstance(I_sf, (float, int, Integer)) or not isinstance(B, (int, Integer)):
        raise TypeError("Inputs (n, B) must be integers, (I_sf) must be float or integer") 
    if n < 4: 
        raise ValueError("Input (n) must be >= 4")
    if is_square(n): # important check; avoids infinite loops
        s = isqrt(n)
    if n % 2 == 0:
        return 2, n//2
    if I_sf <= 0:
        raise ValueError("Input (I_sf) must be > 0")
    if not B:
        # two good lower bounds for B
        if B_fn == 1: 
            B = isqrt(sqrt(n)) # B = n^(1/4), good for n <= 10^11
        elif B_fn == 2: 
            B = int(log(n)**2) # B = ln(n)², good for n >= 10^12
        # good balanced upper bound for B
        elif B_fn == 3: 
            B = ceil(exp(sqrt(log(n)*log(log(n))/2))) # B = e^((ln(n)*ln(ln(n))/2)^(1/2))
        else: 
            raise ValueError("Input (B_fn) must be in [1,2,3]")       
    if B < 2: 
        raise ValueError("Input (B) must be >= 2")

    sr = floor(sqrt(n))
    I = int(B*I_sf) # scaling of interval endpoints

    # the quadratic polynomial is given as Q(z)=(sr+z)²-n with z in interval [-I,I].
    # since we want to construct congruences of squares x²=y²(mod n), for which we will construct x and y²,
    # compute the "x components candidates" (= (sr+z) values) for the whole interval [-I,I]
    x_cc = [sr + z for z in range(-I, I+1)]
    # compute the "y² components candidates" (= Q(z) values) for the whole interval [-I,I]
    y2_cc = [x**2 -n for x in x_cc]
    
    factor_base = []
    sieve = [1] * len(y2_cc)

    def fast_prime_power(n, p):
        if p == 2:
            return 1 << (n & -n).bit_length() - 1
        p_e = 1
        while n % p == 0:
            n //= p
            p_e *= p
        return p_e
    
    # start sieving
    for p in prime_range(B+1):
        # Find the square roots of n modulo p, if they exist
        if p == 2:
            roots = [1] # Because n is odd
        else:
            ls = legendre_symbol(n, p)
            if ls == 0:
                # non-trivial factor of n found
                return p, n//p
            elif ls == 1:
                roots = mod(n, p).sqrt(all=true)
                factor_base.append(p)
            else:
                # no modular square roots exist (ls == -1)
                continue
        for r in roots:
            # Determine the start index for sieving
            start_index = (r - sr + I) %p
            for i in range(start_index, len(y2_cc), p):
                # multiply by prime power
                sieve[i] *= fast_prime_power(y2_cc[i], p)
  
    x_cc_sieved = []
    y2_cc_sieved = []

    # evaluate the sieve elements
    for x, y2, s in zip(x_cc, y2_cc, sieve):
        # compare the recreated B-smooth Q(z) values against the real Q(z) values
        if s == abs(y2):
            # store the sieved candidates
            x_cc_sieved.append(x)
            y2_cc_sieved.append(y2)
 
    factor_base = [-1]+factor_base
    relations = []
    
    for s in y2_cc_sieved:
        f_dict = dict(factor(s)) # dictionary maps prime factors to their exponents
        if s < 0:
            f_dict[-1] = 1 # needs to be set manually!
        exponent_vector_mod_2 = [f_dict.get(p,0)%2 for p in factor_base]
        # collect new relation (represented by exponent vector mod 2)
        relations.append(exponent_vector_mod_2)
    
    M = Matrix(GF(2), zip(*relations)) # zip() transposes columns into rows

    # solve the linear algebra
    null_space = M.right_kernel() # solves M * v = 0 (mod 2) 
    for v in null_space.basis(): # .basis() converts the null space into a list of vectors v
        x = prod(x_cc_sieved[i] for i in range(len(v)) if v[i] == 1)%n
        y2 = prod(y2_cc_sieved[i] for i in range(len(v)) if v[i] == 1)
        y = isqrt(y2)%n
        d = gcd(n, x-y)
        if 1 < d < n:
            return d, n//d  # Factorization succeeded
        d = gcd(n, x+y)
        if 1 < d < n:
            return d, n//d  # Factorization succeeded
            
    return None  # Factorization failed

# ------------------------------------------------------------------------------------------------------------

def QS_Factorization_3(n, I_sf, B, B_fn = False):
    """ Attempts to factor an integer (n) based on the Quadratic Sieve method with smoothness bound (B or B_fn). """
    """ The endpoints of the sieving interval [-I, I] are scaled by the factor (I_sf) applied on the smoothness bound. """
    """ Returns two non-trivial factors of n, or "None" if unsuccessful. """
    """ Characteristics: Uses "reduction by division" for sieving. """
    """                  The smoothness bound can be chosen manually (B) or from a menu of functions of (n). """
    """                  To selcect such function choose (B_fn) between [1,2,3] and set (B) to zero or "False". """
    """ Recommendations: For efficiency reasons, start with a reasonably small choice for the smoothness bound. """
    """                  Start with a small scaling factor (around 0.5), then slowly increase it (up until 2.0  """
    """                  or 3.0) before further increasing the smoothness bound and resetting the scaling factor. """

    if not isinstance(n, (int, Integer)) or not isinstance(I_sf, (float, int, Integer)) or not isinstance(B, (int, Integer)):
        raise TypeError("Inputs (n, B) must be integers, (I_sf) must be float or integer") 
    if n < 4: 
        raise ValueError("Input (n) must be >= 4")
    if is_square(n): # important check; avoids infinite loops
        s = isqrt(n)
        return s, s        
    if n % 2 == 0:
        return 2, n//2
    if I_sf <= 0:
        raise ValueError("Input (I_sf) must be > 0")
    if not B:
        # two good lower bounds for B
        if B_fn == 1: 
            B = isqrt(sqrt(n)) # B = n^(1/4), good for n <= 10^11
        elif B_fn == 2: 
            B = int(log(n)**2) # B = ln(n)², good for n >= 10^12
        # good balanced upper bound for B
        elif B_fn == 3: 
            B = ceil(exp(sqrt(log(n)*log(log(n))/2))) # B = e^((ln(n)*ln(ln(n))/2)^(1/2))
        else: 
            raise ValueError("Input (B_fn) must be in [1,2,3]")       
    if B < 2: 
        raise ValueError("Input (B) must be >= 2")

    sr = floor(sqrt(n))
    I = int(B*I_sf) # scaling of interval endpoints

    # the quadratic polynomial is given as Q(z)=(sr+z)²-n with z in interval [-I,I].
    # since we want to construct congruences of squares x²=y²(mod n), for which we will construct x and y²,
    # compute the "x components candidates" (= (sr+z) values) for the whole interval [-I,I]
    x_cc = [sr + z for z in range(-I, I+1)]
    # compute the "y² components candidates" (= Q(z) values) for the whole interval [-I,I]
    y2_cc = [x**2 -n for x in x_cc]
    
    factor_base = []
    sieve = y2_cc[:] # copy list
    
    # start sieving
    for p in prime_range(B+1):
        # Find the square roots of n modulo p, if they exist
        if p == 2:
            roots = [1] # Because n is odd
        else:
            ls = legendre_symbol(n, p)
            if ls == 0:
                # non-trivial factor of n found
                return p, n//p
            elif ls == 1:
                roots = mod(n, p).sqrt(all=true)
                factor_base.append(p)
            else:
                # no modular square roots exist (ls == -1)
                continue
        for r in roots:
            # Determine the start index for sieving
            start_index = (r - sr + I) %p
            for i in range(start_index, len(y2_cc), p):
                # divide by prime
                while sieve[i] % p == 0:
                    sieve[i] //= p
    
    x_cc_sieved = []
    y2_cc_sieved = []
    
    # evaluate the sieve elements
    for x, y2, s in zip(x_cc, y2_cc, sieve):
        if abs(s) == 1:
            # store the sieved candidates
            x_cc_sieved.append(x)
            y2_cc_sieved.append(y2)

    factor_base = [-1]+factor_base
    relations = []
    
    for s in y2_cc_sieved:
        f_dict = dict(factor(s)) # dictionary maps prime factors to their exponents
        if s < 0:
            f_dict[-1] = 1 # needs to be set manually!
        exponent_vector_mod_2 = [f_dict.get(p,0)%2 for p in factor_base]
        # collect new relation (represented by exponent vector mod 2)
        relations.append(exponent_vector_mod_2)
    
    M = Matrix(GF(2), zip(*relations)) # zip() transposes columns into rows

    # solve the linear algebra
    null_space = M.right_kernel() # solves M * v = 0 (mod 2) 
    for v in null_space.basis(): # .basis() converts the null space into a list of vectors v
        x = prod(x_cc_sieved[i] for i in range(len(v)) if v[i] == 1)%n
        y2 = prod(y2_cc_sieved[i] for i in range(len(v)) if v[i] == 1)
        y = isqrt(y2)%n
        d = gcd(n, x-y)
        if 1 < d < n:
            return d, n//d  # Factorization succeeded
        d = gcd(n, x+y)
        if 1 < d < n:
            return d, n//d  # Factorization succeeded
            
    return None  # Factorization failed

# ------------------------------------------------------------------------------------------------------------
