#!/bin/bash

python3 plus_qvina_as.py \
  --center -17.114 -24.426 -11.868 \
  --size 91.528 146.97 121.868 \
  --ligand_dir /mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/interface/predict_sdf_1tqn \
  --receptor /mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/data/protein/1tqn.pdbqt \
  --qvina /mnt/c/Users/466259/qvina/bin/qvina02 \
  --out_csv /mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/interface/predict_sdf_1tqn/results_1tqn.csv \
  --tmp_dir /mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/interface/predict_sdf_1tqn/tmp_1tqn \
  --error_log /mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/interface/predict_sdf_1tqn/1tqn_qvina_errors.txt
