def reverseSeq(seq,q=True):
    """
    gets nucleutide sequence, returns it's complementary sequance, accepts only ACTG
    by defualt function asserts behaviour, change q to False to avoid program termination because of error
    :param seq: String
    :param q: boolean, True to quit over illegal input
    :return: comp String
    """
    compNucs = {
        'A' : 'T',
        'T' : 'A',
        'G' : 'C',
        'C' : 'G',
        '-' : '-'
    }
    rslt = ""
    try:
        for loc in range(len(seq)):
            rslt+= compNucs[seq[loc]]
    except KeyError as e:
        print("seqTool Error: Illegal char:", seq[loc], "in sequence:",seq ,e)
        if q:
            quit(1)
    except TypeError as e:
        print("seqTool Error: Wrong type was given to reverse function:", type(seq), e)
        if q:
            quit(2)
    return rslt

def legalSeq(seq,legalChars):
    """
    gets sequence, and an iterable collection of chars to check if the sequence is made up only from them
    :param seq: String
    :param legalChars: String
    :return: Boolean
    """
    for chr in seq:
        if chr not in legalChars:
            return False
    return True
