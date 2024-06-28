from Bio.SeqUtils import seq1, seq3
import pandas as pd
import math, os

# Calculate the lfc score of a dipeptide
def lfcScore(dipep="AA"):
    global df
    df = pd.read_csv("./data/dipep_score.csv")
    df = df.set_index("dipep")
    # lfc call from df by df["lfc"].loc["AA"]
    return df["lfc"].loc[dipep]

# Input peptide_sequence (len>1), output overall_score
def seq_lfcScore(seq):
    scores = []
    sum = 0
    for i in range(0, len(seq) - 1):
        aa = seq[i:i+2]
        sum += 2 ** lfcScore(aa)
        scores.append(lfcScore(aa))
    #print([seq[i:i+2] for i in range(0, len(seq) - 1)])
    #print(scores)
    all_score = math.log2(sum / len(seq))
    #print(f"{seq} Seq_Score: {format(all_score, '.3f')}\n")
    return all_score

# Example usage:
#-seq_lfcScore("CIIK")+seq_lfcScore("CIIK")


# input subsequence and location, output the conservative substitution string    
def ConservScore(R_ori,R_mut):
    conscore = con_ref[R_ori].loc[R_mut]
    return conscore





def Batch_subScore(site_opt,tar_seq):
    global con_ref
    con_ref = pd.read_csv("./data/BLOSUM62.csv")
    con_ref = con_ref.set_index('res')
    df = pd.DataFrame(columns=["strategy","Seq_ori","Seq_mut","ΔSeqScore"])
    
    resi1 = int(site_opt.split("_")[0][1:])
    seq_ori = tar_seq[resi1-2:resi1+2]
    Rseq_ori = [f"{seq_ori}_1_{site_opt.split('_')[0]}",f"{seq_ori}_2_{site_opt.split('_')[1]}"]
        
    #generate mutational plan of all aa sustitution
    for Rseq in Rseq_ori:
        aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'] # all aa for mutation
        seq_ori = Rseq.split("_")[0]
        R_pos = int(Rseq.split("_")[1])
        R_ori = Rseq.split("_")[2][0]
        aas.remove(R_ori)
        for aa in aas:
            R_mut = aa
            seq_mut = seq_ori[0:R_pos]+aa+seq_ori[R_pos+1:]
            aasubscore_change = -seq_lfcScore(seq_ori)+seq_lfcScore(seq_mut)
            ResiPos_mut = Rseq.split("_")[2]+aa
            df.loc[len(df)] = [ResiPos_mut,seq_ori,seq_mut,aasubscore_change]
    
    df["Conservative_Score"] = df["strategy"].apply(lambda x: ConservScore(x[0],x[-1]))
    return df
    
    
    
    
    
    
    
    
