''' Message Encoding And Decoding Functions '''
''' by JDU '''

# ----------------------------------------------------------------------------------------------------

from sage.all import *    # See Cryptography/README.md

# ----------------------------------------------------------------------------------------------------

myDefaultChars = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz(.,;:!?)[<+-*/=>]@| "    # (last character is a blank space)

# ----------------------------------------------------------------------------------------------------

def str2num(s, supported = myDefaultChars):
    """ Encodes a string (s) to a positive integer if its characters are supported. """
    
    if not isinstance(s, str):
        raise TypeError("Input (s) must be a string!")

    if not isinstance(supported, str):
        raise TypeError("Input (supported) must be a string!")

    listDigits, skipped = [], []

    for c in s:
        if c in supported:
            # each character (c) in (s) is represented as a digit value ranging from 0 to len(supported)
            listDigits.append(supported.find(c))
        elif c not in skipped:
            # only list each skipped character once
            skipped.append(c)   

    if skipped:
        print(f"Warning: The characters [{', '.join(skipped)}] weren't supported! Consider adding them to myDefaultChars!")
        return None
    
    # reverses the order of the list elements. now the least significant digit becomes the FIRST entry
    listDigits.reverse()

    # applies the weighting of each digit (according to our given base) and adds them up 
    return ZZ(listDigits, len(supported))   

# ----------------------------------------------------------------------------------------------------


def str2num_Alt(s, supported = myDefaultChars):
    """ Encodes a string (s) to a positive integer if its characters are supported. """
    """ This is an alternative implementation of str2num(). """

    if not isinstance(s, str):
        raise TypeError("Input (s) must be a string!")

    if not isinstance(supported, str):
        raise TypeError("Input (supported) must be a string!")

    n, base, weight, skipped = 0, len(supported), 1, []

    # reverses the order of the string elements. now it starts with the last character
    s = s[::-1]

    for c in s:
        if c in supported:
            n += supported.find(c) * weight
            weight *= base
        elif c not in skipped:
            # only list each skipped character once
            skipped.append(c)

    if skipped:
        print(f"Warning: The characters [{', '.join(skipped)}] weren't supported! Consider adding them to myDefaultChars!")
        return None
    
    return n

# ----------------------------------------------------------------------------------------------------

def num2str(n, supported = myDefaultChars):
    """ Decodes an integer (n) back to a string. """
    """ Note: To ensure compatibility with your str2num() function, they """
    """       both must use the exact same input value for "supported". """
    
    if not isinstance(n, (int, Integer)):
        raise TypeError("Input (n) must be an integer!")
    
    if not isinstance(supported, str):
        raise TypeError("Input (supported) must be a string!")

    # reads in the digits of (n) according to our given base
    listDigits = n.digits(len(supported))

    # reverses the order of the list elements. now the least significant digit becomes the LAST entry
    listDigits.reverse()

    return "".join([supported[d] for d in listDigits])

# ----------------------------------------------------------------------------------------------------

def num2str_Alt(n, supported = myDefaultChars):

    """ Decodes an integer (n) back to a string. """
    """ This is an alternative implementation of num2str(). """
    """ Note: To ensure compatibility with your str2num() function, they """
    """       both must use the exact same input value for "supported". """

    if not isinstance(n, (int, Integer)):
        raise TypeError("Input (n) must be an integer!")
    
    if not isinstance(supported, str):
        raise TypeError("Input (supported) must be a string!")

    s = ""
    base = len(supported)

    while n > 0:
        s += supported[n % base]
        n //= base 

    # reverses the order of the string elements
    return s[::-1]

# ----------------------------------------------------------------------------------------------------
