

def long_range_empirical_nontrivial_prob(bits, num_N=100, trials_per_N=10**5):
    """
    Generates 'num_N' random distinct odd semiprimes N, with primes up to 'bits' bits. For each N, performs 'trials_per_N'
    trials to find x,y such that x^2 ≡ y^2 (mod N), for random unique x in range [1, N-1].
    Returns the average probability for x ≠ ±y (mod N), and total congruences of squares tested.
    """

    semiprimes = set()
    while len(semiprimes) < num_N:
        p = random_prime(2**bits+1, lbound=3) # upper bound is not inclusive, therefore +1
        q = random_prime(2**bits+1, lbound=3)
        N = p * q
        if N not in semiprimes:
            semiprimes.add(N)

    total_cong_of_squares_tested = 0
    total_nontrivial = 0

    for N in semiprimes:
        cong_of_squares_tested = 0
        nontrivial = 0
        x_tested = set()

        for _ in range(trials_per_N):
            x = randint(1, N-1)
            if x in x_tested:
                continue
            y2 = power_mod(x,2,N)
            if is_square(y2):
                cong_of_squares_tested += 1
                y = isqrt(y2)
                if (x%N != y%N) and (x%N != (-y)%N):
                    nontrivial += 1
            x_tested.add(x)
            
        total_cong_of_squares_tested += cong_of_squares_tested
        total_nontrivial += nontrivial

    if total_cong_of_squares_tested == 0:
        raise ValueError("No congruences of squares found!")

    avg_prob = total_nontrivial / total_cong_of_squares_tested
    return avg_prob, total_cong_of_squares_tested

