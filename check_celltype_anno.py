#!/bin/python

import anndata as ad
import pandas as pd
import os
import sys

# Check for command-line arguments
if len(sys.argv) < 2:
    print("Usage: python script_name.py <topdir>")
    sys.exit(1)

# Get the top directory from command-line arguments
topdir = sys.argv[1]

# Search for files matching *.h5ad
h5ad_files = [
    os.path.join(root, file)
    for root, _, files in os.walk(topdir)
    for file in files if file.endswith(".h5ad")
]

failed_files = []

for filepath in h5ad_files:
    try:
        adata = ad.read_h5ad(filepath)  # Try to open the file
        metaname = os.path.basename(filepath).replace(".h5ad", "_obs.tsv")
        dirname = os.path.dirname(filepath)
        meta = adata.obs
        meta.to_csv(os.path.join(dirname, metaname), sep="\t")
    except Exception as e:  # Catch errors
        failed_files.append(filepath)

log_file = os.path.join(topdir, "failed_files.log")  # File to store paths of corrupted files

with open(log_file, "w") as f:
    f.write("\n".join(failed_files))

print(f"Processing complete. Corrupted file paths have been written to {log_file}")
