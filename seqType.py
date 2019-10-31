class seqType:
    def __init__(self,snpID,gene,name,seq3,seq5,wt,variants,chromo,repl
                 ,orientation,mutatype,clinical,clinVar="",ranks=[],rslts=[],
                 residue="",readingFrame="",aaPosition=""):
        self.snpID = snpID
        self.geneID = gene
        self.geneName= name
        self.seq3 = seq3
        self.seq5 = seq5
        self.wt = wt
        self.varSeqList = variants #list of strings (possible variants to refallele)
        self.chromo = chromo
        self.baseRepl = repl
        self.clinVar = clinVar
        self.rank = ranks
        self.rslts = rslts #[("VAR1 DNA",[CAS1,...,CASn]),...,("VARn DNA",[CAS1,...,CASn])]
        self.orientation = orientation
        self.mutatype = mutatype
        self.clinical = clinical
        self.aaPosition = aaPosition
        self.readingFrame = readingFrame
        self.residue = residue

    def __repr__(self):
        return "Id's:"+self.snpID+" , "+self.geneID+\
               " \nname:"+self.geneName+\
               " \nWTseq:"+self.seq3[-20:]+"<"+self.wt+">"+self.seq5[:20]+\
               " \nvariants: "+str(self.varSeqList)+\
               " \nchromo: "+self.chromo+\
               " \nBase replacement:"+self.baseRepl+\
               " \nResults: "+str(self.rslts)

    def __eq__(self,other):
        if other == None:
            return False
        if not (isinstance(other,seqType)):
            raise Exception()
        if other.snpID == self.snpID:
            if other.geneID == self.geneID and other.geneName == self.geneName: ###finish adding values for comparison!!!!<<<<>>>>
                return True
        return False

    def __hash__(self):
        return hash(self.seq3[-20:]+self.wt+self.seq5[:20])

    def rankCas(self):
        """
        gets .cas field, uses outside ranking system script
        :return: list of ranked cases?
        """
        #to be implemented
        print("unavailable")

    # def __repr__(self): #for debug
    #     return (self.seq3[-20:]+self.wt+self.seq5[:20])

def duplicate(seq):
    """
    creates a deep copy of a seqType object
    :param seq: SeqType
    :return: seqType
    """
    rslts = []
    for rslt in seq.rslts:
        rslts.append(rslt)
    return seqType(
        seq.snpID,
        seq.geneID,
        seq.geneName,
        seq.seq3,
        seq.seq5,
        seq.wt,
        seq.varSeqList,
        seq.chromo,
        seq.baseRepl,
        seq.orientation,
        seq.mutatype,
        seq.clinical,
        clinVar=seq.clinVar,
        ranks=seq.rank,
        rslts=rslts,
        residue= seq.residue,
        readingFrame = seq.readingFrame,
        aaPosition=seq.aaPosition)