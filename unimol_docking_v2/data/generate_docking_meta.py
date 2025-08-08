#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
import os
import argparse

def calculated_docking_grid_sdf(work_path, json_path, pdbid, ligid, add_size=10):
    input_path = os.path.join(work_path, f'{ligid}.sdf')
    os.makedirs(json_path, exist_ok=True)
    output_grid = os.path.join(json_path, pdbid + '.json')
    mol = Chem.MolFromMolFile(str(input_path), sanitize=False)
    coords = mol.GetConformer(0).GetPositions().astype(np.float32)
    min_xyz = [min(coord[i] for coord in coords) for i in range(3)]
    max_xyz = [max(coord[i] for coord in coords) for i in range(3)]
    center = np.mean(coords, axis=0)
    size = [abs(max_xyz[i] - min_xyz[i]) for i in range(3)]
    center_x, center_y, center_z = center
    size_x, size_y, size_z = size
    size_x += add_size
    size_y += add_size
    size_z += add_size
    grid_info = {
        "center_x": float(center_x),
        "center_y": float(center_y),
        "center_z": float(center_z),
        "size_x": float(size_x),
        "size_y": float(size_y),
        "size_z": float(size_z)
    }
    with open(output_grid, 'w') as f:
        json.dump(grid_info, f, indent=4)
    print(f"{pdbid}: Center=({center_x:.6f}, {center_y:.6f}, {center_z:.6f}), Size=({size_x:.6f}, {size_y:.6f}, {size_z:.6f})")

def create_meta_info_from_ligands(ligand_path, output_csv, default_pdbid=None):
    """
    ??????? meta_info.csv ? ????????? pdb_code ? lig_code ?? ?????? ?????? ? ????? ligand_path.
    ?????????, ??? ????? ????? ??????: <pdbid>_<ligid>.sdf
    ???? default_pdbid ??????, ?? ????? ??????????? ??? ???? ???????.
    """
    records = []
    for filename in os.listdir(ligand_path):
        if filename.endswith(".sdf"):
            name = filename[:-4]  # ?????? ??????????
            if default_pdbid is not None:
                pdbid = default_pdbid
                ligid = name
            else:
                parts = name.split('_', 1)
                if len(parts) == 2:
                    pdbid, ligid = parts
                else:
                    # ???? ?????? ?? ?????????, ????? ?????????? ??? ?????? pdbid ??? unknown
                    pdbid = "unknown"
                    ligid = name
            records.append({"pdb_code": pdbid, "lig_code": ligid})
    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    print(f"Meta info CSV created: {output_csv}")

def main():
    parser = argparse.ArgumentParser(description="Generate docking grid JSON and meta CSV for Uni-Mol Docking")
    parser.add_argument('--protein_path', type=str, required=True, help='???? ? ????? ? ??????? (PDB)')
    parser.add_argument('--ligand_path', type=str, required=True, help='???? ? ????? ? ????????? (SDF)')
    parser.add_argument('--meta_info_file', type=str, required=True, help='CSV ? ????????? pdb_code, lig_code. ???? ????? ???, ?? ????? ??????')
    parser.add_argument('--add_size', type=float, default=10, help='????? ? ??????? ??????-?????')
    parser.add_argument('--output_grid_path', type=str, default=None, help='????? ??? JSON ??????-?????')
    parser.add_argument('--output_meta_csv', type=str, default='posebuster_428_one2one.csv', help='???????? CSV ???? ??? ?????')
    parser.add_argument('--default_pdbid', type=str, default=None, help='???? ??? ??????? ????????? ? ?????? ?????, ??????? pdbid')
    args = parser.parse_args()

    # ???? meta_info_file ?? ??????????, ??????? ???
    if not os.path.isfile(args.meta_info_file):
        print(f"???? {args.meta_info_file} ?? ??????, ??????? ?? ?????? ?????? ? {args.ligand_path}")
        create_meta_info_from_ligands(args.ligand_path, args.meta_info_file, default_pdbid=args.default_pdbid)

    df = pd.read_csv(args.meta_info_file)
    pdb_ids = list(df['pdb_code'].values)
    lig_ids = list(df['lig_code'].values)

    grid_path = args.output_grid_path or f'posebuster428_grid{int(args.add_size)}'
    os.makedirs(grid_path, exist_ok=True)

    # ????????? JSON ??????-?????
    for pdbid, ligid in tqdm(zip(pdb_ids, lig_ids), total=len(pdb_ids)):
        calculated_docking_grid_sdf(args.ligand_path, grid_path, pdbid, ligid, add_size=args.add_size)

    # ???????? CSV ??? ????????? ???????
    df_out = pd.DataFrame(columns=['input_protein', 'input_ligand', 'input_docking_grid', 'output_ligand_name'])
    predict_name_suffix = 'predict'
    for i, (pdbid, ligid) in enumerate(zip(pdb_ids, lig_ids)):
        single_protein_path = os.path.abspath(os.path.join(args.protein_path, pdbid + '.pdb'))
        single_ligand_path = os.path.abspath(os.path.join(args.ligand_path, f'{pdbid}_{ligid}.sdf'))
        single_grid_path = os.path.abspath(os.path.join(grid_path, pdbid + '.json'))
        predict_name = f'{pdbid}_{predict_name_suffix}'
        df_out.loc[i] = [single_protein_path, single_ligand_path, single_grid_path, predict_name]

    print(df_out.info())
    print(df_out.head(3))
    df_out.to_csv(args.output_meta_csv, index=False)
    print(f"Meta CSV saved to {args.output_meta_csv}")

if __name__ == "__main__":
    main()
