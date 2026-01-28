#!/usr/bin/env python
from Bio import SeqIO
import glob
import os
import sys

gene_id = sys.argv[1]
ffn_files = glob.glob(sys.argv[2])
output_fasta = sys.argv[3]

all_gene_ids = set()
with open(gene_id) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if parts:
            gid = parts[0]
            if gid:
                all_gene_ids.add(gid)

with open(output_fasta, "w") as out_f:
  for ffn_file in ffn_files:
    for record in SeqIO.parse(ffn_file, "fasta"):
      header_id = record.id.split()[0]
      if header_id in all_gene_ids:
        out_f.write(f">{record.id}\n{record.seq}\n")
