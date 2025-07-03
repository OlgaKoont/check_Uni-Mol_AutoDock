#!/bin/bash

python scoring.py \
  -d /mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/interface/predict_sdf_4wnv_ki_07 \
  -p /mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/data/protein/4wnv_meeko.pdb \
  -o results_4wnv_ki_07.csv \
  --center  6.0 23.95 -3.86 --box_size 20 20 20 -j 4

#python scoring.py \
#  -d /mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/interface/predict_sdf_4k7a_ic50 \
#  -p /mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/data/protein/4k7a_prep.pdb \#
#  -o results_4k7a_ic50.csv \
#  --center  -17.204 7.265 -12.319 --box_size 40 40 40 -j 4

#python scoring.py \
#  -d /mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/interface/predict_sdf_1tqn \
#  -p /mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/data/protein/1tqn_prep.pdb \
#  -o results_1tqn.csv \
#  --center  -17.114 -24.426 -11.868 --box_size 40 40 40 -j 4

#python3 plus_qvina_as.py \
#  --center 34.569 47.309 90.104 \
#  --size 167.99 141.678 292.256 \
#  --ligand_dir /mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/interface/predict_sdf_1ere \
#  --receptor /mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/data/protein/1ere.pdbqt \
#  --qvina /mnt/c/Users/466259/qvina/bin/qvina02 \
#  --out_csv /mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/interface/predict_sdf_1ere/results_1ere.csv \
#  --tmp_dir /mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/interface/predict_sdf_1ere/tmp_1ere \
#  --error_log /mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/interface/predict_sdf_1ere/1ere_qvina_errors.txt
