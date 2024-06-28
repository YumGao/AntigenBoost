import requests, os
from utils.LEDsearch import LEDipepSearch
from utils.PredictBetabyAlign import predict_beta_tarseq
from utils.SubstitutionScore import Batch_subScore 
from utils.Plot import Dot_plot
from utils.ExtractBetastrandCif import get_beta_strand_residues
import matplotlib.pyplot as plt
import argparse

       

def download_pdb(pdb_id,ref_pdb):
    
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    print(response)
    if response.status_code == 200:
        with open(ref_pdb, 'wb') as f:
            f.write(response.content)
        print(f"PDB file {ref_pdb} downloaded successfully.")
    else:
        print(f"Failed to download PDB file {ref_pdb} ")

#download_pdb("8WSQ")



def select_beta_strand_sites(LEDipeps_loc, tar_beta_id):
    # Convert tar_beta_id to a set for faster lookup
    beta_set = set(tar_beta_id)

    # Initialize a list to store the selected sites
    selected_sites = []

    # Iterate through each site in LEDipeps_loc
    for site in LEDipeps_loc:
        # Split the site into two residues
        residue1, residue2 = site.split('_')
        res1_name, res1_id = residue1[0], int(residue1[1:])
        res2_name, res2_id = residue2[0], int(residue2[1:])
        
        # Check if either residue is in the beta strand set
        if (res1_name, res1_id) in beta_set or (res2_name, res2_id) in beta_set:
            selected_sites.append(site)
    
    return selected_sites

def main():
# Create the parser
    parser = argparse.ArgumentParser(description='Process some files.')

# Add the arguments
    parser.add_argument('-s', '--target_sequence', required=True, help='Target sequence for optimization')
    parser.add_argument('-o', '--output', required=True, help='Output folder path')

    parser.add_argument('-chid', '--chain_id', required=True, default='A',help='Chain id to used of the input structure') # set defaut as "A"



    parser.add_argument('-pred', '--pre_ref_struc', help='Predicted structure of the target sequence')
    parser.add_argument('-id', '--pdb_id', help='The PDB id of homology reference structure')
    global args
    args = parser.parse_args()

    tar_seq = args.target_sequence #Input the target sequence for optimization
    output_folder_path = args.output # output folder

    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)
        print(f"The output folder was created at {output_folder_path}")
    else:
        print(f"The output folder already exists at {output_folder_path}")

    chain_id = args.chain_id

# use predicted structure as input
    if args.pre_ref_struc:
    # The user provided a local path to the predicted structure
    ## use cif file as input
        cif_file = args.pre_ref_struc
        if not os.path.isfile(cif_file):
            parser.error("The file specified with -pred does not exist.")
    # Process the local structure file
        print(f"Processing local structure from: {cif_file}")
        tar_beta_id = get_beta_strand_residues(cif_file, chain_id)
    
    
    
## use pdb_id as input (homology structure)    
    elif args.pdb_id:    
        pdb_id = args.pdb_id.upper() #input pdb id to fetch ref structure
        ref_pdb = f"{pdb_id}.pdb"
        download_pdb(pdb_id, ref_pdb)
        tar_beta_id = predict_beta_tarseq(tar_seq,ref_pdb,chain_id)
    
    else:
    # Neither argument was provided, which should not happen due to required=True
        parser.error("Please provide either a local structure path with -pred or a PDB ID with -id.")

#Run AntigenBoost to search for LEDipep sites
    df, LEDipeps_loc = LEDipepSearch(tar_seq)
#print("Found LEDipeps:\n",df)
    selected_sites = select_beta_strand_sites(LEDipeps_loc, tar_beta_id)


# Print the result
    print("Identified LEDipep sites on the target sequence:", selected_sites)

    plt.switch_backend('Agg')
    print("Let's do optimization!!")

    for i in range(len(selected_sites)):
        site_opt = selected_sites[i]
        df_strategy = Batch_subScore(site_opt,tar_seq)
        print("---------------------------",site_opt,"---------------------------")
        print(df_strategy)
        csv_file = os.path.join(output_folder_path, f'MutationalScore_{site_opt}.csv')
        df_strategy.to_csv(csv_file, index=False)


    # Use the Agg backend for headless environments
        csv_file = os.path.join(output_folder_path, f'{site_opt}.png')     
        Dot_plot(df_strategy, filename=csv_file)

    print("Done!")



if __name__ == "__main__":
    main()

	





