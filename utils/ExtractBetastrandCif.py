from pathlib import Path

from Bio.PDB import DSSP, MMCIFParser, PDBParser


def _load_structure(structure_file):
    structure_path = Path(structure_file)
    suffix = structure_path.suffix.lower()

    if suffix == ".pdb":
        parser = PDBParser(QUIET=True)
    elif suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError(
            f"Unsupported structure format: {suffix}. "
            "Supported formats: .pdb, .cif, .mmcif"
        )

    return parser.get_structure("structure_model", str(structure_path))


def _get_beta_strand_residues(structure_file, chain_id="A"):
    structure = _load_structure(structure_file)
    model = structure[0]
    dssp = DSSP(model, str(structure_file))

    beta_strand_residues = []
    for key in dssp.keys():
        residue_info = dssp[key]
        if residue_info[2] == "E" and key[0] == chain_id:
            beta_strand_residues.append((residue_info[1], key[1][1]))

    return beta_strand_residues


def print_beta_strand_residues(structure_file, chain_id="A"):
    return _get_beta_strand_residues(structure_file, chain_id)


def get_beta_strand_residues(structure_file, chain_id="A"):
    return _get_beta_strand_residues(structure_file, chain_id)
