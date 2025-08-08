#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import pandas as pd

protein_path = "/mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/data/protein/1ere_prep.pdb"
ligand_dir = "/mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/data/ligand/ligands_1ere_ic50"
grid_path = "/mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/data/docking_grids/1ere.json"
output_dir = "predict_sdf"

ligand_files = sorted([f for f in os.listdir(ligand_dir) if f.endswith(".sdf")])

rows = []
for i, lig_file in enumerate(ligand_files, 1):
    ligand_path = os.path.join(ligand_dir, lig_file)
    output_name = f"4wnv_predict_{i}"
    rows.append({
        "input_protein": protein_path,
        "input_ligand": ligand_path,
        "input_docking_grid": grid_path,
        "output_ligand_name": output_name
    })

df = pd.DataFrame(rows)
df.to_csv("batch_one2many_1ere_ic50.csv", index=False)
