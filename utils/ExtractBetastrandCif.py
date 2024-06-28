from Bio.PDB import PDBParser
from Bio.PDB import MMCIFParser, DSSP

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

    



def get_beta_strand_residues(cif_file, chain_id="A"):
    # Parse the CIF file
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('AlphaFold_model', cif_file)
    
    # Get the first model
    model = structure[0]
    
    # Run DSSP on the model
    dssp = DSSP(model, cif_file)
    
    # Collect residues that are part of beta strands
    beta_strand_residues = []
    for key in dssp.keys():
        chain, res_id = key
        if chain == chain_id and dssp[key][2] == 'E':  # 'E' denotes a residue in a beta strand
            beta_strand_residues.append((dssp[key][1], res_id[1]))
    
    return beta_strand_residues

# Example usage with your local CIF file
#cif_file = "RSV_AF3_model_0.cif"  # Replace with the actual path to your CIF file
#chain_id = "A"  # Replace with the actual chain ID if it's different
#beta_strand_residues = get_beta_strand_residues(cif_file, chain_id)
#print("Beta Strand Residues:\n", beta_strand_residues)










