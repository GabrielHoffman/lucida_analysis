#!/usr/bin/env python3

import sys
import argparse
import anndata as ad
from pathlib import Path

def main():

  # Define arguments
  #-------------------

  parser = argparse.ArgumentParser(
      description="Concat h5ad files from the same set of samples"
  )

  parser.add_argument("--files",
      type=Path, 
      required=True,
      help="Input .h5ad file")

  parser.add_argument("--output",      
      type=Path, 
      required=True,
      help="Output .h5ad file")

  # parse args
  args = parser.parse_args()

  # read list of h5a files
  files = open(args.files, "r").read().split('\n')
  
  # remove empty strings
  files = [x for x in files if x.strip()]

  print(files)

  # 1. Read each file into a list of AnnData objects
  adatas = [ad.read_h5ad(path) for path in files]

  # 2. Concatenate them along the observations (cells) axis
  adataComb = ad.concat(adatas, 
    axis = 0, 
    join = "outer", 
    index_unique = "-",)

  adataComb.write_h5ad( args.output )



if __name__ == "__main__":
  main()
