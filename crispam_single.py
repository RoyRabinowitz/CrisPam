#AUTHORS: Shiri Almog shirialmog1@gmail.com, Roy Darnell
from Scripts.casDBFile_new import casDB, CASletter, CASnames ,casDic
from Scripts.seqType import *
from Scripts.seqTools import *
from Bio.Seq import Seq
from Scripts.BEseqType import *
#Libraries
import xml.etree.ElementTree as ET
import csv
import time

badParse = 0
AVERAGERSLTS = 0
totalSequances = 0

def match(CAS,sequence,seq3len,rtrnLoc=False,specific=False):
    """
    gets full Sequence and CAS's PAM, returns Match if the PAM matches the Variant somewhere
    if rtrnLoc = True, returns location of the match instead of True, or False if none found
    :param CAS: String
    :param Variant: String
    :param specific: Boolean - check only the specific location given
    :return: if rtrnLoc -> list with int locations; else Boolean
    """
    window = len(CAS)
    movement = window
    rslt = []
    start = seq3len - window + 1
    if specific == True: #don't move
        movement = 1
        start = seq3len
    for i in range(movement): #window size options
        match = True #flag
        for j in range(window): #compare single chars in window
            try:
                if sequence[start+i+j] not in CASletter[CAS[j]]:
                    match = False
                    # print("false match")  # debug
                    break
            except Exception as e:
                match = False
                print("Corrupt Data, skipping variant",e)
        if match:
            if rtrnLoc:
                rslt.append(start + i)
            else:
                return True
    if rtrnLoc:
        return rslt
    else:
        return False

def matchOnDB(snp):
    """
    :param DB: gets seqType DB
    :return: rsltsDB - {key-String wtSequence:[seqType,[CAS1,CAS2...CASn],...]}
    """
    rslt={}
    currentWT = snp.seq5+snp.wt+snp.seq3
    tempSeq = duplicate(snp)
    variant=snp.varSeqList
    for CAS in casDB: #Go over all PAMseqs
        print (CAS)
        try:
           currentVar = snp.seq5+variant+snp.seq3
           compwt = reverseSeq(currentWT)[::-1] #on complementary we work with the reversed antisense, so the PAM remains the same, but the window starting point changes to seq5 length
           compVar = reverseSeq(currentVar)[::-1]
           matchPositions = match(CAS, currentVar, len(snp.seq5), rtrnLoc=True)
           print (matchPositions)
           rev = False

           if matchPositions != []: #if there are any locations that match on variant
               for matchIndex in matchPositions: #go over each of these locations
                   casname = getCas(CAS)
                   # if (not (match(CAS,currentWT,matchIndex,specific=True))):
                   for name in casname:
                       valid = True
                       for cass in casDic[name]:
                           if (match(cass, currentWT, matchIndex,specific=True)):
                               print (2)
                               loc=match(CAS, currentWT, matchIndex, rtrnLoc=True,specific=True)
                               valid = False
                               break

                           if (match(CAS, currentWT, matchIndex, specific=True)):
                               loc = match(CAS, currentWT, matchIndex, rtrnLoc=True, specific=True)
                               valid = False
                               break

                           if name == "SpCas9" and (match("NAG", currentWT, matchIndex,specific=True)):
                               valid = False
                               break


                           if valid==True: #check if the wt does not match in these specific locations
                               orig_seq, RC = get_printables(snp, variant,CAS, matchIndex)
                               if name in rslt:
                                   rslt[name].update({matchIndex: [orig_seq,RC, CAS,rev]})
                               else:
                                   rslt[name] = {matchIndex:[orig_seq,RC, CAS,rev]}

           rev = True

           matchPositions = match(CAS, compVar, len(snp.seq5), rtrnLoc=True) #again for complementary
           if matchPositions != []:
               for matchIndex in matchPositions:
                   casname = getCas(CAS)
                   # if (not (match(CAS, compwt, matchIndex,specific=True))):
                   print (casname)
                   print (matchIndex)
                   print (compVar)
                   print (compwt)
                   for name in casname:
                       valid = True
                       for cass in casDic[name]:
                           print (cass)
                           if (match(cass, compwt, matchIndex,specific=True)):
                               print (2)
                               loc=match(cass, compwt, matchIndex,rtrnLoc=True,specific=True)
                               print(5)
                               print (loc)
                               valid = False
                               break

                           if (match(CAS, compwt, matchIndex, specific=True)):
                               print (3)
                               loc = match(CAS, compwt, matchIndex, rtrnLoc=True, specific=True)
                               print(loc)
                               valid = False
                               break

                           if name == "SpCas9" and (match("NAG", compwt, matchIndex,specific=True)):
                               valid = False
                               break

                           if valid==True:
                                print (4)
                                orig_seq,RC=get_printables(snp, variant, CAS, matchIndex)
                                if name in rslt:
                                    rslt[name].update({matchIndex:[orig_seq,RC,CAS, rev]})
                                else:
                                    rslt[name] = {matchIndex:[orig_seq,RC, CAS , rev]}



        except Exception as e:
            print("failed to match due to an unexpected error:",e,"\n Statistics unreliable")
    return rslt

def get_printables(snp, variant, CAS, loc):
    #original
    # start=min(len(snp.seq5),30)*(-1)
    # end=min(len(snp.seq3),30)
    start=min(len(snp.seq3),30)*(-1)
    end=min(len(snp.seq3),30)
    if loc+len(CAS)<(len(snp.seq3)+1):
        orig_seq = snp.seq5[start:loc] + "<class style='color:#8A2BE2'>" +snp.seq5[loc:]+"<b>" + variant + "</b></class>" + snp.seq3[:end]
        corrected=snp.seq5[start:loc] + "<class style='color:#8A2BE2'>" +snp.seq5[loc:] +"<b>"+snp.wt+ "</b></class>" + snp.seq3[:end]
        rev_seq3 = Seq(snp.seq3).reverse_complement()
        rev_seq5 = Seq(snp.seq5).reverse_complement()
        RC = rev_seq3[:loc] + "<class style='color:#8A2BE2'>" + rev_seq3[loc:] + "<b>" + Seq(
            variant).reverse_complement() + "</b></class>" + rev_seq5
        # RC = Seq(snp.seq3[:end]).reverse_complement() + "<class style='color:blue'><b>" + Seq(
        #     variant).reverse_complement() + "</b>" + Seq(snp.seq5[loc:]).reverse_complement() + "</class>" + Seq(
        #     snp.seq5[start:loc]).reverse_complement()

        RC_corr = Seq(snp.seq3[:end]).reverse_complement() + "<class style='color:blue'><b>" + Seq(
            snp.wt).reverse_complement() + "</b>" + Seq(snp.seq5[loc:]).reverse_complement() + "</class>" + Seq(
            snp.seq5[start:loc]).reverse_complement()
    else:
        add=loc+len(CAS)-(len(snp.seq3)+1)
        orig_seq=snp.seq5[start:loc] + "<class style='color:#8A2BE2'>" +snp.seq5[loc:] +"<b>" +variant+"</b>"+ snp.seq3[:add]+"</class>" + snp.seq3[add:end]
        corrected=snp.seq5[start:loc] + "<class style='color:#8A2BE2'>" +snp.seq5[loc:] +"<b>" +snp.wt+"</b>"+ snp.seq3[:add]+"</class>" + snp.seq3[add:end]
        rev_seq3=Seq(snp.seq3).reverse_complement()
        rev_seq5=Seq(snp.seq5).reverse_complement()
        RC=rev_seq3[:loc] + "<class style='color:#8A2BE2'>" +rev_seq3[loc:]+"<b>"+ Seq(variant).reverse_complement()+"</b>"+ rev_seq5[:add]+"</class>"+ rev_seq5[add:]
        # RC = Seq(snp.seq3[add:end]).reverse_complement() + "<class style='color:blue'>" + Seq(snp.seq3[:add]).reverse_complement()+ "<b>"+Seq(
        #     variant).reverse_complement() + "</b>" + Seq(snp.seq5[loc:]).reverse_complement() + "</class>" + Seq(
        #     snp.seq5[start:loc]).reverse_complement()

        RC_corr = Seq(snp.seq3[add:end]).reverse_complement() + "<class style='color:blue'>" + Seq(
            snp.seq3[:add]).reverse_complement() + "<b>" + Seq(
            snp.wt).reverse_complement() + "</b>" + Seq(snp.seq5[loc:]).reverse_complement() + "</class>" + Seq(
            snp.seq5[start:loc]).reverse_complement()
    return orig_seq,RC


def gRNA(seq,PAM,loc=True):
    """
    Gets a sequance, PAM and Variant location, returns calculated guideRNA
    :param loc: Variant location, if none passed, assume seq is seqType
    :param seq: seqType, or String and loc provided
    :param PAM: String
    :return:
    """
    if loc: #seqType
        loc = match(PAM, seq, len(seq.seq3), loc=True)
        return reverseSeq(seq[loc:loc + len(PAM)])
    loc = match(PAM,seq,loc,loc=True)
    return reverseSeq(seq[loc:loc+len(PAM)])

def getCasName(Cas):
    if CASnames[Cas] != None:
        return CASnames[Cas]
    else:
        return "Na"

def getCas(Cas):
    li=[]
    for cass in casDic:
        if Cas in casDic[cass]:
            li.append(cass)
    return li

def delete_chars(field):
    field_new=field.replace(" ","")
    for char in field_new:
        if char!='A' and char!='G' and char!='T' and char!='C':
            field_new=field_new.replace(char,"")
    return field_new

def main(upSeq, downSeq, mutation, wt,personalPAM,start,end,fromto,stream):
    upSeq_ = delete_chars(upSeq.upper())
    downSeq_ = delete_chars(downSeq.upper())
    wt_ = delete_chars(wt.upper())
    mutation_ = delete_chars(mutation.upper())
    snp=seqType(1234, 1234, 1234, upSeq_, downSeq_, wt_, mutation_, 1, 1, 1, 1, 1,1)
    if personalPAM != "":
        casDB.append(personalPAM)
        casDic.update({'User customized BE': personalPAM})
        CASnames.update({personalPAM: 'User customized BE'})
    original_sequence=upSeq_+ "<b>"+ wt_+"</b>"+downSeq_
    mutation_sequence = upSeq_ +"<b>"+ mutation_ +"</b>"+ downSeq_
    rslt = matchOnDB(snp) #results DB is a dict, by definition duplicates are not possible\allowed in a set
    print (rslt)
    return rslt,original_sequence,mutation_sequence

if __name__ == "__main__":
    main()