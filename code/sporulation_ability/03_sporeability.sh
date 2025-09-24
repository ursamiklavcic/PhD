#!/bin/bash

#SBATCH --job-name=blast
#SBATCH --nodes=1
#SBATCH --nodelist=heracpu02.nlzoh.si
#SBATCH --cpus-per-task=32 
#SBATCH --mem=50G
#SBATCH --time=14-12:00:00

# Species that were found by MetaPhlAn4. 

## Download genomes 
# Input file with GTDB accessions in the third column
genomes_dir="/home/nlzoh.si/ursmik1/projects/thesis/data/sporulation_ability/genomes"

mkdir -p "$genomes_dir"

# Read accession numbers from the third column 
cut -f3 /home/nlzoh.si/ursmik1/projects/thesis/data/sporulation_ability/mpa_accession_manual.tsv | while read -r accession; do
    # Clean up whitespace
    accession=$(echo "$accession" | tr -d '[:space:]')

    if [[ "$accession" == GCA_* || "$accession" == GCF_* ]]; then
        echo "Downloading genome for accession: $accession"
        curl -s -OJX GET \
          "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${accession}/download?include_annotation_type=GENOME_FASTA&filename=${accession}.zip" \
          -H "Accept: application/zip"

       # Move to output directory if download succeeded
        if [ -f "${accession}.zip" ]; then
            mv "${accession}.zip" "$genomes_dir/"
            echo "Downloaded: ${accession}"
        fi
    else
        echo "Invalid accession format: $accession"
    fi
done

# Extract .fna genome files from .zip files 
mkdir -p "${genomes_dir}/extracted_fna"  # ensure dir exists

for zip in "${genomes_dir}"/*.zip; do
    echo "Processing $zip..."
    temp_dir=$(mktemp -d)
    unzip -q "$zip" -d "$temp_dir"
    find "$temp_dir" -type f -name "*.fna" | while read -r fna_file; do
        cp "$fna_file" "${genomes_dir}/extracted_fna/"
        echo "Extracted: $(basename "$fna_file")"
    done
    rm -rf "$temp_dir"
done


# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate blast

# Set paths
genes_dir="/home/nlzoh.si/ursmik1/projects/spore_evo_dynamics/data/sporulation_signature/gene_aa_sequences"
genomes_dir="/home/nlzoh.si/ursmik1/projects/thesis/data/sporulation_ability/genomes/extracted_fna"
output="/home/nlzoh.si/ursmik1/projects/thesis/data/sporulation_ability"

mkdir -p "${output}/blast_tmp"

# Output file with header
output_file="${output}/mpa_blast_results.tsv"
echo -e "gene_name\tgenome_id\tevalue\tidentity" > "$output_file"

# Loop over each gene FASTA file
for gene_file in "$genes_dir"/*.fasta; do
    gene_name=$(basename "$gene_file" .fasta)

    # Loop over each genome FASTA file
    for genome_file in "$genomes_dir"/*.fna; do
        genome_id=$(basename "$genome_file" .fna)

        # Make BLAST database if missing
        if [ ! -f "${genome_file}.nhr" ]; then
            makeblastdb -in "$genome_file" -dbtype nucl
        fi

        # Temporary BLAST output
        tmp_out="${output}/blast_tmp/${gene_name}_${genome_id}.tmp"

        # Run BLAST
        tblastn -query "$gene_file" -db "$genome_file" -outfmt "6 qseqid sseqid evalue pident" > "$tmp_out"

        # Filter and format results: e-value < 1e-5 & identity > 30
        awk -v gene="$gene_name" -v genome="$genome_id" \
            '$3 < 1e-5 && $4 > 30 {print gene, genome, $3, $4}' OFS="\t" "$tmp_out" >> "$output_file"
    done
done

# Clean up
rm -r "${output}/blast_tmp"