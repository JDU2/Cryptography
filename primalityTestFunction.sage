''' Primality Test Function '''
''' by JDU '''

# ------------------------------------------------------------------------------------------------------------

def primalityTest(n, itr):
    """ Tests whether an integer (n) is a prime number. """
    """ If the given number of iterations (itr) is not sufficient to determine primality, it prints """
    """ out the certainty level of primality and the maximum number of iterations for the worst case. """

    # Test prime factors 2 and 3
    if (not n%2) or (not n%3): 
        return False

    itrLimit = itr      
    sr = floor(sqrt(n))
    wc_itr = ((sr-1)//6)+1   # Worst case iterations

    # Test all other potential prime factors
    for i in range(7, sr+1, 6):
        if not n%i or not n%(i-2): 
            return False
        itr-=1
        if not itr:
            print(f"Result: (n) is possibly prime! ... Certainty: {(1+itrLimit*6)*100/float(sr):.15f} % ... Increase (itr) for more certainty!")
            print(f"If (n) is a prime then in the worst case it would take up to {wc_itr} iterations for 100% certainty.\n")
            return
        
    return True

# ------------------------------------------------------------------------------------------------------------
