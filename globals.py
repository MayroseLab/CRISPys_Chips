# The short name of the organism from plaza.
# organism_abbreviation = "ath"

import os
import pickle
import pandas as pd

# get the path of the scripts directory
CODE_PATH = os.path.dirname(os.path.realpath(__file__))

# for crispritz you need to give a path to folder where there are files with fatsat of each chromosome
# you can make those files using splitfasta https://pypi.org/project/split-fasta/
genome_by_chr_path = "/groups/itay_mayrose/udiland/crispys_chips_arabidopsis/ath_split_files"

pam_file_path = "/groups/itay_mayrose/udiland/crispys_chips_arabidopsis/pam"

gff_file = "/groups/itay_mayrose/udiland/crispys_chips_arabidopsis/annotation.all_transcripts.all_features.ath.gff3.gz"

# To create the gff db file for gffutils in create_family_folders, run the following:
# db = gffutils.create_db(path_to_gff,path_to_output_database, merge_strategy="create_unique")

# -----------------------------------------------------------
# Binary file paths
# -----------------------------------------------------------

goldoff_model_path = "/groups/itay_mayrose/udiland/crispys_off_target/regression_union_log_max.xgb"

crispritz_script_path = "/groups/itay_mayrose/udiland/miniconda3/envs/crunch/bin/crispritz.py"

# Models of scoring functions loadings
moff_loaded_model = None
moff_mtx1 = None
moff_mtx2 = None

# -----------------------------------------------------------
# Off-target file parameters
# -----------------------------------------------------------

# any off-target with score above the threshold will cause the candidate to be filtered
feature_score_dict = {"gene": 0.15, "Any": 0.15,  "exon": 0.15}

# regions that will be ignored
ignore_regions = {"pseudogene", "pseudogenic_transcript"}

# -------------------------------------------------------------------------------
# I used the 'make_family_genes_dict_from_plaza_tbl' function from 'create_families_folders.py' to save pickle file of family_name:gene_in_family dictionary
# this file will be used to ignore off-targets that belong to family genes

familiy_genes_dict_path = "/groups/itay_mayrose/udiland/crispys_chips_arabidopsis/fam_genes_dict.p"