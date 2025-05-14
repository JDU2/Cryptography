def pohlig_hellman(h, g, *pf):
    """ Returns the discrete logarithm of h to base g (mod p). """
    """ h and g must be elements of Zp* (with prime p), and be given in the form of mod(p)-objects. """
    """ g must be a generator of the group Zp*, """
    """ and pf is the prime factorization [(p1,e1),...,(p?,e?)] of the order of the group Zp*. """

    if not h or not g or not pf:
        raise ValueError("Inputs (h),(g),(pf) must be non-zero!")
        
    if not hasattr(h, "modulus") or not hasattr(g, "modulus"):
        raise TypeError("Inputs (h) and (g) must be mod-objects!")

    if not h.modulus() == g.modulus():
        raise ValueError("Inputs (h) and (g) must have the same modulus!")

    p = g.modulus()
    n = g.multiplicative_order()
    
    if n != p-1:
        raise ValueError("g must be a generator of Zp*, with prime p!")
    
    congr_x = []

    for (p_i, e_i) in pf:
        pi_ei = p_i**e_i
        g_i = g**(n // p_i)
        h_k = h
        x_i = 0
        for k in range(e_i):
            exp = n // p_i**(k+1)
            h_k = h_k**exp
            a_k = discrete_log(h_k, g_i, ord=p_i) # uses BSGS as default
            x_i += a_k * p_i**k
            h_k = h * g**(-x_i)
        congr_x.append((x_i, pi_ei))
    
    x = crt([x_i for (x_i, pi_ei) in congr_x], [pi_ei for (x_i, pi_ei) in congr_x])
    return x
