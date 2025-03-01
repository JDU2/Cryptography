''' Factorization methods - Fermat's and Pollard's '''
''' by JDU '''

# ----------------------------------------------------------------------------------------------------

from sage.all import *      # See Cryptography/README.md 

# ----------------------------------------------------------------------------------------------------

def fermatsFactorization(n, itr = 10**6):
    """ Returns two non-trivial factors (a) and (b) of integer (n) after a maximum of (itr) iterations or "false". """
    """ Fermat's formula:     n = ab = (t+s)(t-s) = t²-s² """
    """ Characteristics: Optimized step size by sorting out prime factors 2 and 3. (If you prefer a simpler variant of this code, """
    """                  then just sort out prime factor 2 and use a step size of 2 instead, which does not rely on mod 6 analysis) """
    """ Note: Non-trivial factors do exclude 1 and itself (n). """
    """       If (n) is composite, then (b) is its smallest possible prime factor. """
    """       Prints the maximum number of iterations required for a definite """
    """       result, if no factors were found within (itr) iterations. """

    if not isinstance(n, (int, Integer)) or not isinstance(itr, (int, Integer)):
        raise TypeError("Inputs (n, itr) must be integers!")

    if itr < 1: raise ValueError("Input (itr) must be >= 1")

    if n <= 3:
        if n == 2 or n == 3: print("Input (n) is a prime number!")
        return False

    if n % 2 == 0: return n//2, 2
    if n % 3 == 0: return n//3, 3

    if is_square(n): 
        print("Input (n) is a square value!")
        return False

    # Initial values for t and s
    t = floor(sqrt(n))+1
    s = 0

    # Worst case: (b) = 5 (smallest prime factor after 2 and 3)
    max_s = ((n/5) - 5)/2       
    max_t = sqrt(n + max_s**2)
    
    # (n) can only be in residue classes [1],[5] mod 6
    if n % 6 == 1:
        # (t) can only be in residue classes [1],[2],[4],[5] mod 6
        while True:
            # Adjust start value for (t)
            if t%6 == 1 or t%6 == 4:
                t_1st_step, t_2nd_step = 1, 2
                break
            elif t%6 == 2 or t%6 == 5:
                t_1st_step, t_2nd_step = 2, 1
                break
            t += 1

    elif n % 6 == 5:
        # (t) can only be in residue classes [0],[3] mod 6
        while True:
            # Adjust start value for (t)
            if t%6 == 0 or t%6 == 3:
                t_step = 3
                break
            t += 1

    # Worst case iterations for compositeness
    max_itr = floor((max_t - t)/3)+1

    # Limit iterations for case (n) is prime
    if itr > max_itr:
        itr = max_itr

    # (n) can only be in residue classes [1],[5] mod 6
    if n % 6 == 1:
        for _ in range(itr):
            if is_square(t**2 - n):
                s = sqrt(t**2 - n)
                break
            t += t_1st_step
            if is_square(t**2 - n):
                s = sqrt(t**2 - n)
                break
            t += t_2nd_step

    elif n % 6 == 5:
        for _ in range(itr):
            if is_square(t**2 - n):
                s = sqrt(t**2 - n)
                break
            t += t_step

    # No integer square value of (s) was found 
    if not s:
        if itr < max_itr:
            print(f"No non-trivial factors (a),(b) found within {itr} iterations!")
            print(f"{max_itr} iterations required (at worst case) for a definite result.")
        else: 
            print("Input (n) is a prime number!")
        return False
    
    a = (t+s)
    b = (t-s)
    
    return a, b

# ----------------------------------------------------------------------------------------------------
