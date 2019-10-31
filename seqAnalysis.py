#modules:
from Scripts.casDBFile import * #in production add : Scripts.
from Scripts.seqType import *
from Scripts.seqTools import *
from Scripts.Main import *
#libraries:
from Bio.pairwise2 import *

def parseSeq(WT,VAR):
    try:
        loc = -1 #variation location
        WT = str.upper(WT)
        VAR = str.upper(VAR)
        start = "-"
        if not(legalSeq(WT+VAR, "-AGTCagtc")):
            return "Error: Illegal sign, input must only contain: a,g,t,c,A,G,T,C,-"
        if WT == VAR:
            return "Error: Nothing to calculate, sequences match"
        elif  len(WT) < len(VAR):
            return "Error: Reference allele must be equal or longer than variant"
        elif len(WT) > len(VAR):
            VAR = align.globalxx(WT, VAR)[0]
            VAR = VAR[1]
            print("fVAR: "+str(VAR))
        for i in range(len(WT)):
            if start == "-" and VAR[i] != "-":  # not at position fix by allignment
                start = VAR[i]  # passed position fix
            if WT[i] != VAR[i] and start != "-":
                if loc != -1:
                    return "Error: maximum of 1 variations are allowed between sequences (multiple differences)"
                loc = i  # where to cut sequence
        try:
            seq3 = WT[0:loc]
            refallele = WT[loc]
            variations = [VAR[loc]]
            seq5 = WT [loc+1:loc+21]
        except Exception as e:
            return "Error: Please provide at least 20 BP from the variation up and down stream"
        if len(seq3) < 20:
            return "Error: Please provide at least 20 BP from the variation up and down stream"

        if (100 > len(seq3 + refallele + seq5) > 40):  # sequence size limitations
            return seqType("NA", "NA", "NA", seq3, seq5, refallele, variations, "NA", "NA", "NA", "NA","NA")
        else:
            return "Error: length must be 40-100 BP"
        return seq
    except Exception as e:
        return "Error: Unexpected script error"

def analyseSeq(seq):
    """
    gets sequence and searches all of the PAMs on it compared to it'so analyze, only seqTypes are handled
        returns dictionary
    :param seq:
    :return:
    """
    if isinstance(seq,str): #got parsing error, passing it on, nothing to process here
        return seq
    currentWT = seq.seq3[-20:] + seq.wt + seq.seq5[:20] #cut irrelevant areas
    rslt = {} #{"CAS Type":"Locations   and comp or not"}
    for variant in seq.varSeqList:  # Go over all variants in seq>variants
        for CAS in casDB:  # Go over all PAMseqs
            try:
                currentVar = seq.seq3 + variant + seq.seq5
                compwt = reverseSeq(currentWT)[::-1]
                compVar = reverseSeq(currentVar)[::-1]
                matchPositions = match(CAS, currentVar, len(seq.seq3), rtrnLoc=True)
                if matchPositions != []:  # if there are any locations that match on variant
                    for matchIndex in matchPositions:  # go over each of these locations
                        if (not (match(CAS, currentWT, matchIndex,specific=True))):  # check if the WT does not match in these specific locations
                            if rslt.get(CAS) == None:  # update DB - new entry
                                rslt[CAS] = getCasName(CAS) + ": " + CAS + " on index: " + str(matchIndex)
                            else:
                                rslt[CAS] += ", " + str(matchIndex)
                matchPositions = match(CAS, compVar, len(seq.seq5), rtrnLoc=True)  # again for complementary
                if matchPositions != []:
                    for matchIndex in matchPositions:
                        if (not (match(CAS, compwt, matchIndex, specific=True))):
                            if rslt.get(CAS) == None:  # update DB - new entry
                                rslt[CAS] = getCasName(CAS) + ": "+ CAS + " on index: comp. " + str(matchIndex)
                            else:
                                rslt[CAS] += ", comp. " + str(matchIndex)
            except Exception as e:
                print("failed to match due to an unexpected error:", e, "\n Statistics unreliable")
    return rslt

def parseWebRslts(rslts):
    if isinstance(rslts,str):
        return rslts
    if rslts == {}:
        return "No results found"
    final = ""
    for rslt in rslts:
        print("rslt in rslts:",rslts[rslt])
        final += str(rslts[rslt])+".\n"
    return final

def process_sequence(WT,VAR):
    return parseWebRslts((analyseSeq(parseSeq(WT,VAR))))

def tests():
    print("<<<<<<<<<<Starting tests sequence:>>>>>>>>>>>")
    WT = "AAAAAAAAAAAAAATTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCC" #G to C
    VAR ="AAAAAAAAAAAAAATTTTTTTTTTTTTTTTCGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCC"
    test(WT,VAR)
    print("G to C - PASS")
    WT = "AAAAAAAAAAAAAATTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCC" #indle in VAR
    VAR ="AAAAAAAAAAAAAATTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCC"
    test(WT, VAR)
    print("indle in var - PASS")
    WT = "AAAAAAAAAAAAAATTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCC" #indle in WT
    VAR ="AAAAAAAAAAAAAATTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCC"
    test(WT, VAR)
    print("indle in WT - FAIL")
    WT = "AAAAAAAAAAAAAATTTTTTTTTTTTTTTTGGGGGGGGGGGGggGGGGCCCCCcCCCCCCCCCCCCC" #small letters
    VAR ="AaaAAAAAAAAAAATTTTTTTTTTTTTTTTCGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCC"
    test(WT, VAR)
    print("small letters - PASS")
    WT = "AAAAAAAAAAAAAA???TTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCC" #undefined sign
    VAR ="AAAAAAAAAAAAAATTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCC"
    test(WT, VAR)
    print("illegal sign - FAIL")
    WT = "AAAAAAAAAAAAAATTTTTTTTTTTTTTTTGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCC"  # multiple indles in WT
    VAR ="AAAAAAAAAAAAAATTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCC"
    test(WT, VAR)
    print("multi indle WT - FAIL")
    WT = "AAAAAAAAAAAAAATTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCC"  # multiple indles in VAR
    VAR ="AAAAAAAAAAAAAATTTTTTTTTTTTTTTTGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCC"
    test(WT, VAR)
    print("multi indle VAR - depends where, FAIL")
    WT = "AAAAAAAAAAAAAATTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCC"  # multiple variations
    VAR ="AAAAAATAAAAAAATTTTTTTATTTTTTTTGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCC"
    test(WT, VAR)
    print("multi variation - FAIL")
    WT = "AAAAAAAAAAAAAATTTTTTTTTTTTTTTT"  # Too Short
    VAR ="AAAAAAAAAAAAAGTTTTTTTTTTTTTTTT"
    test(WT, VAR)
    print("Too short - FAIL")
    WT = "AAAAAAAAAAAAAATTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCC" #short var
    VAR ="AAAATTTTTTTTTTTTTTTTCGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCC"
    test(WT, VAR)
    print("short var at start - PASS")

def test(WT,VAR):
    print("-------------")
    print("with: ")
    print("WT:  ", WT)
    print("VAR: ", VAR)
    print(process_sequence(WT, VAR))