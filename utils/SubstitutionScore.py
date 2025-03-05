from Bio.SeqUtils import seq1, seq3
import pandas as pd
import math
import os

# Load data once to optimize performance
df_lfc = pd.read_csv("./data/dipep_score.csv", dtype={"dipep": str}).set_index("dipep")

# Ensure index is uppercase (if needed) and explicitly handle "NA" as a string
df_lfc.index = df_lfc.index.str.upper().fillna("NA")  

# Ensure index is uppercase and remove NA values
df_lfc.index = df_lfc.index.str.upper()
df_lfc = df_lfc.dropna()  # Remove rows with NA values


# Calculate the LFC score of a dipeptide
def lfcScore(dipep="AA"):
    dipep = dipep.upper()  # Ensure uppercase
    if dipep in df_lfc.index:
        return df_lfc["lfc"].loc[dipep]
    else:
        print(f"Warning: Dipeptide {dipep} not found in dataset. Returning default score 0.")
        return 0  # Default value to prevent KeyError


# Compute overall LFC score for a given sequence (length > 1)
def seq_lfcScore(seq):
    if len(seq) < 2:
        print("Warning: Sequence too short for LFC scoring.")
        return 0

    sum_scores = 0
    for i in range(len(seq) - 1):
        dipep = seq[i:i+2]
        score = lfcScore(dipep)
        sum_scores += 2 ** score

    all_score = math.log2(sum_scores / (len(seq) - 1))
    return all_score


# Get conservative substitution score from BLOSUM62
def ConservScore(R_ori, R_mut):
    try:
        con_ref = pd.read_csv("./data/BLOSUM62.csv").set_index('res')
        if R_ori in con_ref.index and R_mut in con_ref.columns:
            return con_ref.loc[R_ori, R_mut]
        else:
            print(f"Warning: No substitution score for {R_ori} -> {R_mut}. Returning 0.")
            return 0
    except Exception as e:
        print(f"Error reading BLOSUM62.csv: {e}")
        return 0


# Compute substitution score for a batch of mutations
def Batch_subScore(site_opt, tar_seq):
    df = pd.DataFrame(columns=["strategy", "Seq_ori", "Seq_mut", "ΔSeqScore"])

    try:
        resi1 = int(site_opt.split("_")[0][1:])
    except (IndexError, ValueError) as e:
        print(f"Error: Invalid site_opt format '{site_opt}'. Expected format like 'H208_R209'.")
        return df  # Return empty DataFrame

    seq_ori = tar_seq[resi1-2:resi1+2]
    Rseq_ori = [f"{seq_ori}_1_{site_opt.split('_')[0]}", f"{seq_ori}_2_{site_opt.split('_')[1]}"]

    all_aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    for Rseq in Rseq_ori:
        try:
            seq_ori = Rseq.split("_")[0]
            R_pos = int(Rseq.split("_")[1])
            R_ori = Rseq.split("_")[2][0]

            if R_ori not in all_aa:
                print(f"Warning: Unexpected amino acid '{R_ori}' in {Rseq}. Skipping.")
                continue

            aas = [aa for aa in all_aa if aa != R_ori]  # Remove original residue

            for aa in aas:
                seq_mut = seq_ori[:R_pos] + aa + seq_ori[R_pos+1:]
                aasubscore_change = -seq_lfcScore(seq_ori) + seq_lfcScore(seq_mut)
                ResiPos_mut = Rseq.split("_")[2] + aa
                df.loc[len(df)] = [ResiPos_mut, seq_ori, seq_mut, aasubscore_change]

        except Exception as e:
            print(f"Error processing mutation batch for {Rseq}: {e}")
            continue

    df["Conservative_Score"] = df["strategy"].apply(lambda x: ConservScore(x[0], x[-1]))
    return df

