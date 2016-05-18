import sys

def reverse_and_complement( s ):
    """Takes a DNA string, takes the complement and returns it reversed"""
    
                
    r = ''
    for c in s:
        if c == 'a':
            r = 't' + r
        elif c == 'A':
            r = 'T' + r
        elif c == 'c':
            r = 'g' + r
        elif c == 'C':
            r = 'G' + r
        elif c == 'g':
            r = 'c' + r
        elif c == 'G':
            r = 'C' + r
        elif c == 't':
            r = 'a' + r
        elif c == 'T':
            r = 'A' + r
        elif c == 'n':
            r = 'n' + r
        elif c == 'N':
            r = 'N' + r
        elif c == 'x':
            r = 'x' + r
        elif c == 'X':
            r = 'X' + r
        # ambiguity codes, see
        # http://www.in-silico.com/s_restriction/Nucleotide_ambiguity_code.html:
        elif c == 's':
            r = 's' + r
        elif c == 'S':
            r = 'S' + r
        elif c == 'm':
            r = 'k' + r
        elif c == 'M':
            r = 'K' + r
        elif c == 'k':
            r = 'm' + r
        elif c == 'K':
            r = 'M' + r
        elif c == 'w':
            r = 'w' + r
        elif c == 'W':
            r = 'W' + r

        elif c == 'r':
            r = 'y' + r
        elif c == 'R':
            r = 'Y' + r
            
        elif c == 'y':
            r = 'r' + r
        elif c == 'Y':
            r = 'R' + r

        elif c == 'v':
            r = 'b' + r
        elif c == 'V':
            r = 'B' + r

        elif c == 'h':
            r = 'd' + r
        elif c == 'H':
            r = 'D' + r

        elif c == 'd':
            r = 'h' + r
        elif c == 'D':
            r = 'H' + r

        elif c == 'b':
            r = 'v' + r
        elif c == 'B':
            r = 'V'+ r


        else:
            sys.exit("reversecomplement: %s is not DNA"%s)
    
                
    return r


