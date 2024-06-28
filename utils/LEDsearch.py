import re
import pandas as pd

# Amino acid pair in letter3
aa_mapping = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"
}

def AA_code_3to1(ori_aa3):
    one_letter = {v: k for k, v in aa_mapping.items()}
    return one_letter.get(ori_aa3.upper())

def initial_pre():
    global df_lfc
    df_lfc = pd.read_csv("./data/lfc_data.csv")
    df_LED = df_lfc[df_lfc["lfc"] < -1.8]
    
    df = pd.DataFrame()
    df["LEDipep1"] = df_LED.apply(lambda row: AA_code_3to1(row["aa1"]) + AA_code_3to1(row["aa2"]), axis=1)
    df["position"] = ""  # Initialize position as a string type
    return df

def LEDipepSearch(aa_seq): 
    df = initial_pre()
    
    LED1_list = df["LEDipep1"].tolist()
    df = df.set_index("LEDipep1")
    LEDipeps_loc = []
    for L in LED1_list:
        pos_rec = []
        for match in re.finditer(L, aa_seq):
            pos = match.start()
            pos_rec.append(f"{pos+1}_{pos+2}")
            LEDipeps_loc.append(f"{L[0]}{pos+1}_{L[1]}{pos+2}")
            #print(f"Find {L} seq in position {pos+1}, the seq is {aa_seq[pos-1:pos+3]}")
        df.at[L, "position"] = ",".join(pos_rec)
    #print(LEDipeps_loc)    
    return df, LEDipeps_loc

def lfcScore(dipep="AA"):
    df = pd.read_csv("./data/dipep_score.csv")
    df = df.set_index("dipep")
    return df.loc[dipep, "lfc"]
    
# Example usage
#aa_seq = "MELLILKANAITTILTAVTFCFASGQNITEEFYQSTCSAVSKGYLGALRTGWYTSVITIELSNIKENKCNGTDAKVKLIKQELDKYKNAVTDLQLLMQSTPATGSGSAICSGVAVCKVLHLEGEVNKIKSALLSTNKAVVSLSNGVSVLTFKVLDLKNYIDKQLLPILNKQSCSIPNIETVIEFQQKNNRLLEITREFSVNAGVTTPVSTYMLTNSELLSLINDMPITNDQKKLMSNNVQIVRQQSYSIMCIIKEEVLAYVVQLPLYGVIDTPCWKLHTSPLCTTNTKEGSNICLTRTDRGWYCDNAGSVSFFPQAETCKVQSNRVFCDTMNSRTLPSEVNLCNVDIFNPKYDCKIMTSKTDVSSSVITSLGAIVSCYGKTKCTASNKNRGIIKTFSNGCDYVSNKGVDTVSVGNTLYCVNKQEGQSLYVKGEPIINFYDPLVFPSDEFDASISQVNEKINQSLAFIRKSDELLHNVNAGKSTTNIMITTIIIVIIVILLSLIAVGLLLY"
#df, LEDipeps_loc = LEDipepSearch(aa_seq)
#print(df)
#print(LEDipeps_loc)

