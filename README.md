## Overview

AntigenBoost performs a sequence optimization process using several utility functions. It downloads PDB files, predicts beta strands, identifies LEDipep sites, and generates mutational scores and plots.

## Requirements

- Python 3.x
- Required Python libraries:
  - `biopython`
  - `pandas`
  - `matplotlib`
  - `requests`

Install the required libraries using the following command:

```bash
pip install -r requirements.txt
```


## Usage
### Running the Script
To run the script, use the following command:
```bash
python main.py -s <target_sequence> -o <output_folder> -chid <chain_id> [-pred <predicted_structure>] [-id <pdb_id>]
```

### Arguments
-s, --target_sequence (required): The target sequence for optimization.
-o, --output (required): The output folder path where results will be saved.
-chid, --chain_id (required): The chain ID to be used from the input structure. Default is 'A'.
-pred, --pre_ref_struc (optional): Path to the predicted structure of the target sequence (cif file supported).
-id, --pdb_id (optional): The PDB ID of the homology reference structure.

### Example
```bash
python main.py -s "MKTAYIAKQRQISFVKSHFSRQDILDLWQYHIYEKYK" -o "./output" -chid "A" -id "1A2B"
```
This example will:
1. Download the PDB file for the given PDB ID (1A2B).
2. Predict beta strands from the homology reference structure.
3. Search for LEDipep sites on the target sequence.
4. Identify LEDipep sites on the target sequence.
5. Generate and save mutational score CSV files and corresponding plots.


```bash
python main.py -s "MKTAYIAKQRQISFVKSHFSRQDILDLWQYHIYEKYK" -o "./output" -chid "A" -pred "predict_structure.cif"
```
This example will:
1. Use the local predicted structure file (predict_structure.cif) as input.
2. Extract beta strand residues from the predicted structure.
3. Search for LEDipep sites on the target sequence.
4. Identify LEDipep sites on the target sequence.
5. Generate and save mutational score CSV files and corresponding plots.

### Output
The script generates:

CSV files containing mutational scores for each identified LEDipep site.
PNG files with corresponding plots of mutational scores for each identified LEDipep site.
The output files are saved in the specified output folder.


