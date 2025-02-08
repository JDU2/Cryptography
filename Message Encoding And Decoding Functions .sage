'''Some Message Encoding And Decoding Functions'''
'''by J-Wells'''

# ----------------------------------------------------------------------------------------------------

myDefaultChars = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz(.,;:!?)[<+-*/=>]@| "   # (last character is a blank space)

# ----------------------------------------------------------------------------------------------------

def str2num(myStr, supported = myDefaultChars):
    """ converts any string to a positive integer """
    
    listDigits, skipped = [], []
    
    for c in myStr:
        if c in supported:
            # each character in myStr is represented as a digit value ranging from 0 to len(supported)
            listDigits.append(supported.find(c))
        elif c not in skipped:
            skipped.append(c)
            
    if len(listDigits) < len(myStr):
        print(f"Warning: The characters [{', '.join(skipped)}] weren't supported! Consider adding them to myDefaultChars!")
        return None
        
    # reverses the order of the list elements such that we have the least significant digit as the first entry
    listDigits.reverse()
    
    # applies the weighting of each digit (according to our given base) and adds them up 
    return ZZ(listDigits, len(supported))   
