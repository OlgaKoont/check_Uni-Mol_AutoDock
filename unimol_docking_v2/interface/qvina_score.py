import subprocess
import os

def sdf_to_pdbqt_ligand(sdf_path, pdbqt_path):
    """
    Converts ligand from SDF to PDBQT using OpenBabel.
    """
    cmd = ['obabel', sdf_path, '-O', pdbqt_path, '-xh', '--gen3d', '--partialcharge', 'gasteiger']
    subprocess.run(cmd, check=True)

def run_qvina(qvina_path, receptor_pdbqt, ligand_pdbqt, center, size, output_pdbqt, log_file):
    """
    Runs qvina02 to calculate binding energy for an existing ligand pose.
    """
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
    subprocess.run(cmd, check=True)

if __name__ == '__main__':
    # === CONFIGURABLE PARAMETERS ===
    # Directory with input files (protein and ligand)
    input_dir = '/mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/interface/predict_sdf_1tqn'
    pdb_dir = '/mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/data/protein'
    # Directory for results
    output_dir = '/mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/interface/predict_sdf_1tqn/qvina'
    # Path to QVina binary
    qvina_path = '/mnt/c/Users/466259/qvina/bin/qvina02'

    # File names (can be changed)
    ligand_sdf = '4wnv_predict_95.sdf'
    receptor_pdbqt = '1tqn.pdbqt'
    ligand_pdbqt = 'ligand.pdbqt'
    output_pdbqt = 'ligand_out.pdbqt'
    log_file = 'qvina.log'

    # Docking box parameters (adjust for your system)
    center = [-17.114, -24.426, -11.868]
    size = [91.528, 146.97, 121.868]

    # === PATH PREPARATION ===
    ligand_sdf_path = os.path.join(input_dir, ligand_sdf)
    receptor_pdbqt_path = os.path.join(pdb_dir, receptor_pdbqt)
    ligand_pdbqt_path = os.path.join(output_dir, ligand_pdbqt)
    output_pdbqt_path = os.path.join(output_dir, output_pdbqt)
    log_file_path = os.path.join(output_dir, log_file)

    os.makedirs(output_dir, exist_ok=True)

    # === FILE CONVERSION ===
    sdf_to_pdbqt_ligand(ligand_sdf_path, ligand_pdbqt_path)

    # === RUN QVina ===
    run_qvina(qvina_path, receptor_pdbqt_path, ligand_pdbqt_path, center, size, output_pdbqt_path, log_file_path)

    # === OUTPUT BINDING ENERGY ===
    with open(log_file_path, 'r') as f:
        for line in f:
            if line.strip().startswith('REMARK VINA RESULT:'):
                print('Binding energy:', line.strip())
