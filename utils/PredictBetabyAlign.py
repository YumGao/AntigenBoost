from Bio.PDB import DSSP, PDBParser, PDBIO, Select
from Bio.SeqUtils import seq1
from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio.pairwise2 import format_alignment



class AtomSelect(Select):
    def accept_atom(self, atom):
        # Define criteria to select atoms to keep
        if atom.get_parent().get_id()[0] == " ":
            # Keep atoms not associated with water or small molecules
            return True
        else:
            # Exclude atoms associated with water or small molecules
            return False

def clean_pdb(input_pdb, output_pdb):
    # Create PDB parser
    parser = PDBParser()

    # Parse the PDB file
    structure = parser.get_structure("pdb_structure", input_pdb)

    # Define a new structure for selected atoms
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, select=AtomSelect())
    print("Removed water and small molecules.")

# Example usage:
#input_pdb = "4WE8.pdb"
#output_pdb = "4WE8_clean.pdb"
#clean_pdb(input_pdb, output_pdb)

def obtain_full_sequence_from_PDB(pdb_file):
    # Open the PDB file
    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    # Initialize sequence string
    sequence = ''

    # Iterate through lines to find the sequence in the header
    for line in lines:
        if line.startswith('SEQRES'):
            # Extract residue names from SEQRES records
            residues = line.split()[4:]
            # Convert 3-letter code to 1-letter code
            sequence += ''.join([seq1(residue) for residue in residues])
    return sequence


def obtain_residues_from_PDB(pdb_file, chain_id="A"):
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('X', pdb_file)
    
    # Get the first model in the structure
    model = structure[0]
    chain = model[chain_id]
    residues = chain.get_residues()
    
    # Create a list to store the sequence and residue IDs
    sequence = []
    residue_ids = []

    for residue in residues:
        if residue.get_id()[0] == ' ':
            sequence.append(seq1(residue.get_resname()))
            residue_ids.append(residue.get_id()[1])
    
    # Join the sequence list into a string
    sequence_str = ''.join(sequence)

    return sequence_str, residue_ids

#example
#pdb_file = "4WE8.pdb"
#sequence, residue_ids = obtain_residues_from_PDB(pdb_file, chain_id="A")
#print("Sequence:", sequence)
#print("Residue IDs:", residue_ids)
    

def print_beta_strand_residues(pdb_file, chain_id="A"):
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('X', pdb_file)

    # Get the first model in the structure
    model = structure[0]

    # Run DSSP on the model
    dssp = DSSP(model, pdb_file)

    record_res_loc = []
    # Iterate over residues in DSSP output
    for key in dssp.keys():
        residue_info = dssp[key]
        # Check if the residue is part of a beta strand and matches the chain ID
        full_chain_id = structure[0][key[0]].get_full_id()[2]
        if residue_info[2] == 'E' and full_chain_id == chain_id:
            res_name = residue_info[1]
            res_id = key[1][1]
            record_res_loc.append((res_name, res_id))

    return record_res_loc
#return format:[('T', 29), ('E', 30), ('E', 31), ('F', 32), ('Y', 33), ('S', 38), ('A', 39), ('V', 40)]

    
def align_sequences(seq1, seq2, gap_open_penalty=-4, gap_extension_penalty=-0.5):
    # Load the BLOSUM62 matrix for scoring
    blosum62 = substitution_matrices.load("BLOSUM62")

    # Perform the alignment using globalds (global alignment with defined scoring matrix)
    alignments = pairwise2.align.globalds(
        seq1, seq2, blosum62, 
        gap_open_penalty, gap_extension_penalty
    )

    # Get the best alignment (highest score)
    best_alignment = alignments[0]
    aligned_seq1, aligned_seq2, score, start, end = best_alignment

    # Calculate identity and similarity scores
    identity = 0
    similarity = 0
    total_matches = 0  # Total matching positions (ignoring gaps)
    for a, b in zip(aligned_seq1, aligned_seq2):
        if a != '-' and b != '-':  # Ignore positions with gaps
            total_matches += 1
            if a == b:
                identity += 1
            if blosum62[a][b] > 0:  # Consider as a similarity if the BLOSUM62 score is positive
                similarity += 1

    # Calculate the percentages
    identity_percentage = (identity / total_matches) * 100 if total_matches > 0 else 0
    similarity_percentage = (similarity / total_matches) * 100 if total_matches > 0 else 0

    return best_alignment, identity_percentage, similarity_percentage
    

    


def save_alignment_to_fasta(alignment, output_file='align.fasta'):
    aligned_seq1, aligned_seq2, score, start, end = alignment
    with open(output_file, 'w') as f:
        f.write(f">aligned_seq1\n{aligned_seq1}\n")
        f.write(f">aligned_seq2\n{aligned_seq2}\n")
    

def find_beta_strand_residues(ref_seq, tar_seq, ref_resid, ref_beta_id, aligned_ref_seq, aligned_tar_seq):
    # Initialize variables to keep track of positions
    ref_index = 0
    tar_index = 0

    # List to store beta strand residues in target sequence
    tar_beta_id = []

    # Iterate through aligned sequences
    for ref_res, tar_res in zip(aligned_ref_seq, aligned_tar_seq):
        # Check if current ref_res is not a gap and ref_index is within bounds
        if ref_res != '-' and ref_index < len(ref_resid):
            # Check if this residue is on a beta strand in the reference sequence
            if (ref_res, ref_resid[ref_index]) in ref_beta_id:
                # If tar_res is not a gap, add to the target beta strand residues
                if tar_res != '-':
                    tar_beta_id.append((tar_res, tar_index + 1))
            # Move to the next reference residue
            ref_index += 1
        # Move to the next target residue if it is not a gap
        if tar_res != '-':
            tar_index += 1

    return tar_beta_id



def predict_beta_tarseq(tar_seq,ref_pdb,chain_id):
    
    
#obtain the reference residue id of residues located on betastrand
    ref_beta_id = print_beta_strand_residues(ref_pdb,chain_id)
    print(f"Found {len(ref_beta_id)} residues located on beta strand.\n")#,ref_beta_id)

# remove small molecules
    cleaned_pdb = ref_pdb.replace(".pdb","_cleaned.pdb")
    clean_pdb(ref_pdb, cleaned_pdb)
    print("Save cleaned reference pdb structure into:", cleaned_pdb)

# obtain reference sequence data and residue id
    ref_resid = obtain_residues_from_PDB(cleaned_pdb,chain_id)
    ref_seq_full = obtain_full_sequence_from_PDB(ref_pdb)
    ref_seq, ref_resid = obtain_residues_from_PDB(cleaned_pdb,chain_id)
    print(f"ref_seq:\n{ref_seq}")
    #print(f"ref_residue_id:\n{ref_resid}")



#sequence alignment for ref_seq and tar_seq
    best_alignment,identity_score, similarity_score =  align_sequences(tar_seq, ref_seq_full)
    print("Pair-wise alignment with BLOSUM62:")
    print("Identity:",identity_score)
    print("Similarity:",similarity_score)
    best_alignment,identity_score, similarity_score =  align_sequences(tar_seq, ref_seq)
    aligned_tar_seq, aligned_ref_seq, score, start, end = best_alignment
    
# Print the alignment
    #print("Pair-wise alignment results\n",format_alignment(*best_alignment))
    save_alignment_to_fasta(best_alignment)


# obtain the target res_id on betastrand.
    tar_beta_id = find_beta_strand_residues(ref_seq, tar_seq, ref_resid, ref_beta_id, aligned_ref_seq, aligned_tar_seq)
    print("Predicted Beta Strand Residues for the target sequence:", tar_beta_id)
    return tar_beta_id



