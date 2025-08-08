import argparse
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem

def clean_pdb_chain(input_pdb, output_pdb, chain_id="A"):
    """
    Save only ATOM records from the specified chain from input_pdb to output_pdb.
    """
    with open(input_pdb, "r") as fin, open(output_pdb, "w") as fout:
        for line in fin:
            if line.startswith("ATOM") and line[21].strip() == chain_id:
                fout.write(line)

def prepare_protein_pdbqt(pdb_file, pdbqt_file):
    """
    Prepare a protein PDBQT file from a cleaned PDB file using Meeko and RDKit.
    """
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    if mol is None:
        raise RuntimeError(f"Failed to read protein from {pdb_file} with RDKit")
    mol = Chem.AddHs(mol)
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    writer = PDBQTWriterLegacy()
    result = writer.write_string(setups[0])
    pdbqt_str = result[0] if isinstance(result, tuple) else result
    with open(pdbqt_file, "w") as f:
        f.write(pdbqt_str)
    print(f"PDBQT file written to {pdbqt_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Pre-clean PDB and prepare protein PDBQT using Meeko and RDKit.")
    parser.add_argument('--pdb', required=True, help='Input protein PDB file (may contain water/ligands/multiple chains)')
    parser.add_argument('--clean_pdb', default='clean_protein.pdb', help='Output cleaned PDB file (only one chain)')
    parser.add_argument('--pdbqt', required=True, help='Output protein PDBQT file')
    parser.add_argument('--chain', default='A', help='Chain identifier to keep (default: A)')
    args = parser.parse_args()

    # Step 1: Clean PDB (keep only one chain, only ATOM records)
    clean_pdb_chain(args.pdb, args.clean_pdb, chain_id=args.chain)
    print(f"Cleaned PDB saved as {args.clean_pdb} (chain {args.chain})")

    # Step 2: Prepare PDBQT using Meeko
    prepare_protein_pdbqt(args.clean_pdb, args.pdbqt)




'''
import subprocess
import os
import tempfile
import shutil

def prepare_protein(input_pdb, output_pdb):
    """
    Prepares a protein PDB file for docking:
    1. Fix missing atoms/residues with pdbfixer
    2. Add hydrogens with reduce
    3. Remove waters with OpenBabel
    """
    tmp_dir = tempfile.mkdtemp()
    try:
        fixed_pdb = os.path.join(tmp_dir, 'fixed.pdb')
        fixed_h_pdb = os.path.join(tmp_dir, 'fixed_h.pdb')
        clean_pdb = os.path.join(tmp_dir, 'clean.pdb')

        # Step 1: Fix missing atoms and residues with pdbfixer
        subprocess.run([
            'pdbfixer', input_pdb, 
            '--output=' + fixed_pdb, 
            '--add-atoms=heavy', '--add-residues'
        ], check=True)

        # Step 2: Add hydrogens and protonate with reduce
        with open(fixed_h_pdb, 'w') as out_f:
            subprocess.run(['reduce', fixed_pdb], stdout=out_f, check=True)

        # Step 3: Remove waters with OpenBabel
        subprocess.run([
            'obabel', fixed_h_pdb, '-O', clean_pdb, '-xr'
        ], check=True)

        # Copy to final output
        shutil.copyfile(clean_pdb, output_pdb)
        print(f"Protein preparation completed: {output_pdb}")

    finally:
        shutil.rmtree(tmp_dir)

def pdb_to_pdbqt_receptor(pdb_path, pdbqt_path):
    """
    Converts protein from PDB to PDBQT using OpenBabel.
    """
    cmd = ['obabel', pdb_path, '-O', pdbqt_path, '-xh', '--partialcharge', 'gasteiger']
    subprocess.run(cmd, check=True)

if __name__ == '__main__':
    # Input and output paths
    input_pdb = '/mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/data/protein/1ere.pdb'
    prepped_pdb = '/mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/data/protein/1ere_prep_new.pdb'
    pdbqt_path = '/mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/data/protein/1ere_prep.pdbqt'

    # Step 1-3: Prepare protein (fix, protonate, remove waters)
    prepare_protein(input_pdb, prepped_pdb)

    # Step 4: Convert to PDBQT
    pdb_to_pdbqt_receptor(prepped_pdb, pdbqt_path)
    print(f'Converted to PDBQT: {pdbqt_path}')
'''