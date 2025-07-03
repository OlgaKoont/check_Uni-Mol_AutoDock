import os
import argparse
import csv
import multiprocessing
from functools import partial
from rdkit import Chem
from vina import Vina
import shutil
import uuid
import subprocess

def prepare_rigid_receptor_pdbqt(input_pdb, output_pdbqt, chain_to_keep=None):
    filtered_pdb = output_pdbqt.replace('.pdbqt', '_filtered.pdb')
    with open(input_pdb, 'r') as fin, open(filtered_pdb, 'w') as fout:
        for line in fin:
            if line.startswith('ATOM'):
                if chain_to_keep is None or line[21] == chain_to_keep:
                    fout.write(line)
    cmd = [
        'obabel', filtered_pdb,
        '-O', output_pdbqt,
        '-xr',
        '--addpolarh',
        '--partialcharge', 'gasteiger'
    ]
    subprocess.run(cmd, check=True)
    if os.path.exists(filtered_pdb):
        os.remove(filtered_pdb)
    if not os.path.exists(output_pdbqt):
        raise RuntimeError(f"Failed to create receptor PDBQT: {output_pdbqt}")
    print(f"Rigid receptor prepared without ROOT tags: {output_pdbqt}")
    return output_pdbqt

def convert_sdf_to_pdbqt(sdf_path, output_dir):
    base_name = os.path.splitext(os.path.basename(sdf_path))[0]
    pdbqt_path = os.path.join(output_dir, base_name + '.pdbqt')
    mol2_path = os.path.join(output_dir, base_name + '.mol2')
    cmd1 = f"obabel {sdf_path} -O {mol2_path} 2>/dev/null"
    os.system(cmd1)
    cmd2 = f"obabel {mol2_path} -O {pdbqt_path} --partialcharge gasteiger --AddPolarH 2>/dev/null"
    os.system(cmd2)
    if os.path.exists(mol2_path):
        os.remove(mol2_path)
    if os.path.exists(pdbqt_path):
        return pdbqt_path
    else:
        raise RuntimeError(f"Failed to convert {sdf_path} to pdbqt.")

def process_ligand(sdf_path, receptor_pdbqt, center, box_size, output_csv=None):
    try:
        mol = Chem.MolFromMolFile(sdf_path)
        smiles = Chem.MolToSmiles(mol) if mol else ""
        tmp_dir = f"tmp_pdbqt_{uuid.uuid4()}"
        os.makedirs(tmp_dir, exist_ok=True)
        ligand_pdbqt = convert_sdf_to_pdbqt(sdf_path, tmp_dir)
        v = Vina(sf_name='vina')
        v.set_receptor(receptor_pdbqt)
        v.set_ligand_from_file(ligand_pdbqt)
        v.compute_vina_maps(center=center, box_size=box_size)
        affinity = v.score()[0]
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        result = {'file': sdf_path, 'smiles': smiles, 'affinity': affinity}
        if output_csv:
            with open(output_csv, 'a', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=['file', 'smiles', 'affinity'])
                writer.writerow(result)
        return result
    except Exception as e:
        print(f"Error processing {sdf_path}: {e}")
        if 'tmp_dir' in locals() and os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        return {'file': sdf_path, 'smiles': '', 'affinity': None}

def main():
    parser = argparse.ArgumentParser(description='Rigid receptor docking scoring with AutoDock Vina')
    parser.add_argument('-d', '--sdf_dir', required=True, help='Directory with ligand SDF files')
    parser.add_argument('-p', '--protein_pdb', required=True, help='Protein PDB file for receptor preparation')
    parser.add_argument('-o', '--output', required=True, help='CSV output file')
    parser.add_argument('-j', '--jobs', type=int, default=2, help='Number of parallel jobs')
    parser.add_argument('--center', nargs=3, type=float, required=True, metavar=('X', 'Y', 'Z'), help='Docking box center')
    parser.add_argument('--box_size', nargs=3, type=float, default=[40,40,40], metavar=('X','Y','Z'), help='Docking box size')
    parser.add_argument('--chain', type=str, default=None, help='Chain to keep in receptor (default: all chains)')
    args = parser.parse_args()

    receptor_pdbqt = os.path.splitext(args.protein_pdb)[0] + '_prepared.pdbqt'
    if not os.path.exists(receptor_pdbqt):
        print(f"Preparing receptor PDBQT from {args.protein_pdb} ...")
        prepare_rigid_receptor_pdbqt(args.protein_pdb, receptor_pdbqt, chain_to_keep=args.chain)
    else:
        print(f"Using existing receptor PDBQT: {receptor_pdbqt}")

    sdf_files = [os.path.join(args.sdf_dir, f) for f in os.listdir(args.sdf_dir) if f.endswith('.sdf')]

    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['file', 'smiles', 'affinity'])
        writer.writeheader()

    worker = partial(process_ligand,
                     receptor_pdbqt=receptor_pdbqt,
                     center=args.center,
                     box_size=args.box_size,
                     output_csv=args.output)

    with multiprocessing.Pool(processes=args.jobs) as pool:
        results = pool.map(worker, sdf_files)

    print(f"\nProcessed ligands: {len(results)}")
    print(f"Results saved to: {args.output}")

if __name__ == '__main__':
    main()
