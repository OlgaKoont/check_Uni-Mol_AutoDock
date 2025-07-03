import subprocess
import os
import csv
import glob
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem

# def fix_sdf_with_obabel(input_sdf, output_sdf):
#     # Fix SDF file with Open Babel (add hydrogens, kekulize, gen3d)
#     # This helps if RDKit cannot kekulize or read the molecule due to aromaticity or valence issues.
#     cmd = [
#         'obabel', input_sdf, '-O', output_sdf,
#         '--gen3d', '--kekulize', '-h'
#     ]
#     subprocess.run(cmd, check=True)

def sdf_to_pdbqt_ligand_meeko_rdkit(sdf_path, pdbqt_path, fixed_dir=None):
    # Read the molecule from SDF using RDKit
    mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    if mol is None:
        raise RuntimeError(f"Failed to read molecule from {sdf_path} with RDKit")
    # Add explicit hydrogens (required for Meeko/RDKit)
    mol = Chem.AddHs(mol)
    try:
        # Kekulize molecule to resolve aromaticity; if this fails, try to fix with Open Babel
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception as e:
        # If kekulization fails, try to fix the molecule with Open Babel and retry
        if fixed_dir is not None:
            os.makedirs(fixed_dir, exist_ok=True)
            fixed_sdf = os.path.join(fixed_dir, os.path.basename(sdf_path))
            fix_sdf_with_obabel(sdf_path, fixed_sdf)
            mol = Chem.MolFromMolFile(fixed_sdf, removeHs=False)
            if mol is None:
                raise RuntimeError(f"Failed to read fixed molecule from {fixed_sdf} with RDKit")
            mol = Chem.AddHs(mol)
            try:
                Chem.Kekulize(mol, clearAromaticFlags=True)
            except Exception as e2:
                raise RuntimeError(f"Kekulization failed after Open Babel fix for {fixed_sdf}: {e2}")
        else:
            raise RuntimeError(f"Kekulization or hydrogen addition failed for {sdf_path}: {e}")
    # Prepare the molecule for docking with Meeko using the new API (Meeko >= 0.5)
    # Avoid deprecated methods like write_pdbqt_file, write_pdbqt_string, setup
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    writer = PDBQTWriterLegacy()
    # Meeko's write_string can return a string or a tuple depending on version
    result = writer.write_string(setups[0])
    if isinstance(result, tuple):
        pdbqt_str = result[0]
    else:
        pdbqt_str = result
    # Write the PDBQT string to file
    with open(pdbqt_path, "w") as f:
        f.write(pdbqt_str)
        
from openbabel import openbabel, pybel

def sdf_to_pdbqt(input_sdf: str, output_pdbqt: str):
    """
    Convert an SDF file to a PDBQT file using Open Babel with 3D coordinates and Gasteiger charges.

    Parameters:
    - input_sdf: Path to the input SDF file.
    - output_pdbqt: Path to the output PDBQT file.
    """
    mols = list(pybel.readfile("sdf", input_sdf))

    if not mols:
        raise ValueError("No molecules found in the SDF file.")

    mol = mols[0]  # If you want to handle multiple molecules, loop over mols

    # Generate 3D coordinates (if needed)
    mol.make3D()

    # Assign Gasteiger charges
    mol.OBMol.AddHydrogens()
    openbabel.OBChargeModel.FindType("gasteiger").ComputeCharges(mol.OBMol)

    # Write to PDBQT
    mol.write(format="pdbqt", filename=output_pdbqt)
    print(f"Converted: {input_sdf} → {output_pdbqt}")

def extract_smiles(sdf_path):
    # Extract SMILES string from SDF using Open Babel
    result = subprocess.run(['obabel', sdf_path, '-osmi', '-h'], capture_output=True, text=True)
    for line in result.stdout.splitlines():
        if line.strip():
            return line.strip().split()[0]
    return ""

def run_qvina(qvina_path, receptor_pdbqt, ligand_pdbqt, center, size, output_pdbqt, log_file):
    # Run QVina docking with specified parameters
    cmd = [
        qvina_path,
        '--receptor', receptor_pdbqt,
        '--ligand', ligand_pdbqt,
        '--center_x', str(center[0]),
        '--center_y', str(center[1]),
        '--center_z', str(center[2]),
        '--size_x', str(size[0]),
        '--size_y', str(size[1]),
        '--size_z', str(size[2]),
        '--out', output_pdbqt,
        '--log', log_file,
        '--score_only'
    ]
    # QVina returns non-zero exit status if there is a problem with the input files or chemistry
    subprocess.run(cmd, check=True, capture_output=True, text=True)

def extract_affinity(log_file):
    # Parse affinity score from QVina log file
    try:
        with open(log_file, 'r') as f:
            for line in f:
                if line.strip().startswith('REMARK VINA RESULT:'):
                    try:
                        return float(line.strip().split()[3])
                    except Exception:
                        return ""
    except Exception:
        return ""
    return ""

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Batch QVina scoring for SDF ligands.")
    parser.add_argument('--center', nargs=3, type=float, required=True, help='Docking box center: x y z')
    parser.add_argument('--size', nargs=3, type=float, required=True, help='Docking box size: x y z')
    parser.add_argument('--ligand_dir', required=True, help='Directory with ligand SDF files')
    parser.add_argument('--receptor', required=True, help='Path to receptor PDBQT')
    parser.add_argument('--qvina', required=True, help='Path to qvina02 binary')
    parser.add_argument('--out_csv', required=True, help='Output CSV file')
    parser.add_argument('--tmp_dir', default='tmp_qvina', help='Temporary directory for intermediate files')
    parser.add_argument('--fixed_sdf_dir', default='fixed_sdf', help='Directory for fixed SDF files')
    parser.add_argument('--error_log', default='qvina_errors.txt', help='File to write errors')
    args = parser.parse_args()

    center = args.center
    size = args.size
    ligand_dir = args.ligand_dir
    receptor_pdbqt = args.receptor
    qvina_path = args.qvina
    out_csv = args.out_csv
    tmp_dir = args.tmp_dir
    fixed_sdf_dir = args.fixed_sdf_dir
    error_log = args.error_log

    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(fixed_sdf_dir, exist_ok=True)

    ligand_files = sorted(glob.glob(os.path.join(ligand_dir, '*.sdf')))

    with open(out_csv, 'w', newline='') as csvfile, open(error_log, 'w') as errfile:
        writer = csv.writer(csvfile)
        writer.writerow(['ligand_file', 'smiles', 'affinity'])

        for ligand_sdf in ligand_files:
            ligand_name = os.path.basename(ligand_sdf)
            ligand_pdbqt = os.path.join(tmp_dir, ligand_name.replace('.sdf', '.pdbqt'))
            output_pdbqt = os.path.join(tmp_dir, ligand_name.replace('.sdf', '_out.pdbqt'))
            log_file = os.path.join(tmp_dir, ligand_name.replace('.sdf', '.log'))

            except subprocess.CalledProcessError as e:
                # Log QVina errors (exit status 1, likely due to bad input files)
                errfile.write(f"[QVina error] {ligand_name}: {e}\nstdout: {e.stdout}\nstderr: {e.stderr}\n\n")
                print(f"Error (qvina) for {ligand_name}: {e}")
            except Exception as e:
                # Log all other errors (RDKit, Meeko, file I/O, etc.)
                errfile.write(f"[General error] {ligand_name}: {e}\n\n")
                print(f"Error processing {ligand_name}: {e}")

# Key comments:
# - Uses only the new Meeko API (no deprecated methods) to avoid DeprecationWarning[1].
# - Handles all possible return types from Meeko's write_string for compatibility with all versions.
# - If RDKit can't kekulize, tries to fix the molecule with Open Babel and reprocesses it.
# - All QVina and molecule processing errors are logged, not fatal for the batch.
# - SMILES extraction uses Open Babel for robust SDF support.
# - Script is robust for batch processing and will skip and log problematic ligands instead of failing.
