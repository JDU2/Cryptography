# ------------------------------------------------------------------------------------------------------------

from sage.all import *    # See Cryptography/README.md

# ------------------------------------------------------------------------------------------------------------

def rsa_NFS_to_symmetric(rsa_bits):
    c = (64/9)**(1/3)
    x = rsa_bits
    # computational effort in bits (for the Number Field Sieve on RSA), where x is the RSA modulus in bits
    ce_bits = ceil(exp(c*(x*log(2))**(1/3)*log(x*log(2))**(2/3))).bit_length()
    # Note: If the computational effort in bits for an exhaustive search, on a symmetric cryptosystem, is 'ce_bits', 
    # then it can find a symmetric key up to a bit-size of 'ce_bits' 
    return ce_bits

def symmetric_to_rsa_NFS(sym_bits):
    c = (64/9)**(1/3)
    # f(x) := (work factor for NFS in bits) - (symmetric key in bits)
    f = lambda x: ceil(exp(c*(x*log(2))**(1/3)*log(x*log(2))**(2/3))).bit_length() - sym_bits
    min_x = 512 # lower bound for RSA bits
    max_x = 32768 # upper bound for RSA bits
    if f(min_x) > 0:
        raise ValueError("sym_bits too small for RSA ≥512 bits")
    if f(max_x) < 0:
        raise ValueError("sym_bits too large for RSA ≤32768 bits")
    # attempts to find x for f(x)=0, within given bounds for x, where x is the RSA modulus in bits
    return ceil(find_root(f, min_x, max_x))

# ------------------------------------------------------------------------------------------------------------

def rsa_QS_to_symmetric(rsa_bits):
    x = rsa_bits
    # computational effort in bits (for the Quadratic Sieve on RSA), where x is the RSA modulus in bits
    ce_bits = ceil(exp(sqrt(x*log(2)*log(x*log(2))))).bit_length()
    # Note: If the computational effort in bits for an exhaustive search, on a symmetric cryptosystem, is 'ce_bits', 
    # then it can find a symmetric key up to a bit-size of 'ce_bits' 
    return ce_bits

# ------------------------------------------------------------------------------------------------------------

def compare_ce_QS_vs_NFS():
    for x in range(8,512+1,8):
        print(f"{x}-bit RSA modulus:  Computational effort in bits:  QS: {rsa_QS_to_symmetric(x)}, NFS: {rsa_NFS_to_symmetric(x)}")
    return

# ------------------------------------------------------------------------------------------------------------
