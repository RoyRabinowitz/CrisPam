#AUTHORS: Shiri Almog shirialmog1@gmail
casDB = ["NGG","NGCG","NGAG","NGAN","NNGRRT","NNNRRT","NNNNGATT","NNRGAAW","NAAAAC","NNNVRYM","TTTV","NANG","NNGTGA","TTN","NRRH","NRCH","TATV","TYCV","TGTV","VTTV","TTV","TTTT","TTCN","NGGN","NG", "NRTH","NGGNG"]
CASletter = {
                "N":["A","C","T","G"],
                "A":["A"],
                "C":["C"],
                "T":["T"],
                "G":["G"],
                "U":["U"],
                "W": ["A","T"],
                "S": ["C","G"],
                "M": ["A","C"],
                "K": ["G","T"],
                "R": ["A","G"],
                "Y": ["C","T"],
                "B": ["C","G","T"],
                "D": ["A","G","T"],
                "H": ["A","C","T"],
                "V": ["A","C","G"],

}

casDic= {
    "SpCas9":["NGG"],
    "VRER-SpCas9":["NGCG"],
    "EQR-SpCas9":["NGAG"],
    "VQR-SpCas9":["NGAN"],
    "SaCas9":["NNGRRT"],
    "KKH-SaCas9": ["NNNRRT"],
    "NmCas9":["NNNNGATT" ],
    "St1Cas9":["NNRGAAW"],
    "TdCas9":["NAAAAC"],
    "CjCas9":["NNNVRYM"],
    "SpCas9-(From-Strep.-past.)":["NNGTGA"],
    "LbCas12a&AsCas12a":["TTTV"],
    "FnCas12a":["TTV"],
    "xCas9":["NG"] ,
    "AsCas12a-RVR":["TATV"],
    "AsCas12a-RR":["TYCV"],
    "enAsCas12a":["TGTV","VTTV","TTTT","TTCN"],
    "St3Cas9":["NGGNG"],
    "SpCas9-NG":["NG","NANG"],
    "SpCas9-NRRH":["NRRH","NGGN"],
    "SpCas9-NRTH":["NRTH", "NGGN"],
    "SpCas9-NRCH": ["NRCH", "NGGN"],
    "Cas12e-(CasX)":["TTCN"]
    }

CASnames = {
    "NGG":"SpCas9",
    "NGCG":"VRER-SpCas9",
    "NGAG":"EQR-SpCas9",
    "NGAN":"VQR-SpCas9",
    "NNGRRT":"SaCas9",
    "NNNRRT":"KKH-SaCas9",
    "NNNNGATT":"NmCas9",
    "NNRGAAW":"St1Cas9",
    "NAAAAC":"TdCas9",
    "NNNVRYM":"CjCas9",
    "NNGTGA":"SpCas9-(From-Strep)",
    "TTTV":"LbCas12a&AsCas12a",
    "TTV":"FnCas12a",
    "NG":["xCas9","SpCas9-NG"],
    "TATV":"AsCas12a-RVR",
    "TYCV":"AsCas12a-RR",
    "TGTV":"enAsCas12a",
    "VTTV":"enAsCas12a",
    "TTTT": "enAsCas12a",
    "TTCN": ["enAsCas12a","Cas12e-(CasX)"],
    "NANG": "SpCas9-NG",
    "NRRH": "SpCas9-NRRH",
    "NGGN":["SpCas9-NRRH","SpCas9-NRTH", "SpCas9-NRCH"],
    "NRTH": "SpCas9-NRTH",
    "NRCH": "SpCas9-NRCH",
    "NGGNG": "St3Cas9",






}