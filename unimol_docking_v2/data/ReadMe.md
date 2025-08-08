1. Для подготовки директории с .sdf файлами лигандов используйте split_sdf.py (пример запуска можно найти в файле)
2. Для подготовки batch_one2many используйте generate_docking_meta.py
   
Пример запуска:

python generate_docking_meta.py   --protein_path "/mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/data/protein"   --ligand_path "/mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/data/ligand/ligands_4wnv_ic50"   --meta_info_file "/mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/data/meta_info_4wnv_ic50.csv"   --default_pdbid 4wnv   --add_size 10   --output_grid_path "/mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/data/docking_grids"   --output_meta_csv "/mnt/c/Users/466259/Uni-Mol/unimol_docking_v2/data/posebuster_428_one2one_4wnv_ic50.csv"
