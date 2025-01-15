#!/user/bin/python3

from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
from pathlib import Path
import argparse
import os
import json
import scipy
import gzip
import subprocess

parent = "/space/scratch/gemma-single-cell-data-ensembl-id/mus_musculus/"


study_names = os.listdir(parent)
for study_name in study_names:
    study_path = os.path.join(parent, study_name)
    outdir = f"/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/processed_data/mus_musculus/{study_name}"
    # Run the Python script
    meta_path = f"/space/gemmaData/metadata/{study_name}/ct_final.tsv"
    try:
        meta_df=pd.read_csv(meta_path, sep="\t")
    except Exception as e:
        print(f"Error reading {meta_path}: {e}")
        continue
    command = f"python gemma_test_preprop.py --study_path {study_path} --study_name {study_name} --meta_path {meta_path} --outdir {outdir}"
    subprocess.run(command, shell=True)
