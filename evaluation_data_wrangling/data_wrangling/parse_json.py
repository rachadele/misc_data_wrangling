#!/bin/python

import json
import os
import pandas as pd

parent = "/space/scratch/gemma-single-cell-data-ensembl-id/"
#json_file = "/space/scratch/gemma-single-cell-data-ensembl-id/GSE198014/GSE198014_meta.json"

dirs = [os.path.join(parent,d) for d in os.listdir(parent) \
        if os.path.isdir(os.path.join(parent, d))]



def process_sample(sample):
    sample_dict = {}
    
    # Extract basic details using .get()
    sample_id = sample.get("id")
    accession = sample.get("accession", {}).get("accession")
    organism = sample.get("arrayDesign", {}).get("taxon", {}).get("scientificName", "").replace(" ", "_").lower()
    
    # Add basic details to the dictionary
    sample_dict["sample_id"] = sample_id
    sample_dict["accession"] = accession
    sample_dict["organism"] = organism
    
    # Extract characteristics
    for characteristic in sample.get("sample", {}).get("characteristics", []):
        category = characteristic.get('category', '')
        value = characteristic.get('value', '')
        if category:
            sample_dict[category] = value
    
    # Extract additional details if available
    bio_source = sample.get("sample", {}).get("BioSource")
    gene = sample.get("Gene")
    #ncbi_id = sample.get("accession", {}).get("externalDatabase", {}).get("id")
    
    if bio_source:
        sample_dict["BioSource"] = bio_source
    if gene:
        sample_dict["Gene"] = gene
  #  if ncbi_id:
      #  sample_dict["NCBI_id"] = ncbi_id
    for item in sample["sample"]["factorValues"]:
        characteristics = item.get("characteristics")
        if characteristics:
            for characteristic in characteristics:
                if characteristic.get("category") == "developmental stage":
                    sample_dict["developmental_stage"] = characteristic.get("value")
        
    # Convert to DataFrame
    return pd.DataFrame([sample_dict])

    
json_files = []
for dir in dirs:
    for file in os.listdir(dir):
        if file.endswith(".json"):
            json_files.append(os.path.join(parent,dir,file))
    
study_meta = {}
            
# Process each JSON file
for json_file in json_files:
    sample_meta = pd.DataFrame()
    with open(json_file) as f:
        print(f"Processing file: {json_file}")
        data = json.load(f)
        for sample in data["data"]:
            # Process each sample and collect its metadata
            sample_data = process_sample(sample)
            # Append to the DataFrame
            sample_meta = pd.concat([sample_meta, pd.DataFrame(sample_data)], ignore_index=True)
        # Add the DataFrame to the study metadata dictionary
        study_meta[json_file.split("/")[-1].split("_")[0]] = sample_meta
        
# Save the metadata to a tsv file
for study in study_meta:
    study_meta[study].to_csv(f"{study}_metadata.tsv", sep="\t", index=False)
    print(f"Saved metadata for study {study} to {study}_metadata.tsv")