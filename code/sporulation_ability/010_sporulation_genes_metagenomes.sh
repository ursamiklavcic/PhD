#!/bin/bash
#SBATCH --job-name=mapping
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=1-00:00:00
#SBATCH --output=map.out

# ------------------------------- User input -----------------------------------------
genes_dir="/home/nlzoh.si/ursmik1/projects/spore_evo_dynamics/data/sporulation_signature/gene_nc_sequences"
FASTQ_DIR="/home/storage/raw_data/BGI/F23A910000041_BGI_metagenomic/raw_fq"
index_dir="/home/nlzoh.si/ursmik1/projects/thesis/data/sporulation_ability/metagenomes/index"
data_dir="/home/nlzoh.si/ursmik1/projects/thesis/data/sporulation_ability/metagenomes/data"


set -e
eval "$(conda shell.bash hook)"
source /home/nlzoh.si/larbez1/miniconda3/etc/profile.d/conda.sh 
conda activate /home/nlzoh.si/larbez1/miniconda3/envs/qgene_env

# Create output directories 
mkdir -p "$index_dir"
mkdir -p "$data_dir"

# Combine all sporulation genes into one file!
FASTA="${index_dir}/nc_sporulation_genes.fasta" 
cat "$genes_dir"/*.fasta > "$FASTA"

# ==============================================================================
# CREATE GTF FILE
# ==============================================================================

GTF="${index_dir}/sporulation_genes.gtf"

echo "Creating GTF annotation..."

awk '
BEGIN{OFS="\t"}

/^>/{
    if(seqid!=""){
        print seqid,"custom","CDS",1,seqlen,".","+",".","gene_id \""seqid"\""
    }

    seqid=substr($1,2)
    seqlen=0
    next}
{
    seqlen += length($0)}
END{
    print seqid,"custom","CDS",1,seqlen,".","+",".","gene_id \""seqid"\""}
' "$FASTA" > "$GTF"

# ==============================================================================
# GENE LENGTHS
# ==============================================================================

GL="${index_dir}/gene_lengths.txt"

echo "Calculating gene lengths..."

awk '
BEGIN{FS="\t"; OFS="\t"}
{
    split($9,a,"\"")
    gene=a[2]
    len=$5-$4+1
    print gene,len
}
' "$GTF" > "$GL"

# ==============================================================================
# BUILD BOWTIE2 INDEX
# ==============================================================================

INDEX_PREFIX="${index_dir}/sporulation_db"

if [ ! -f "${INDEX_PREFIX}.1.bt2" ]; then

    echo "Building Bowtie2 index..."

    bowtie2-build "$FASTA" "$INDEX_PREFIX"
fi

# Mapping
for FILE in $FASTQ_DIR/*1.fq; do
  while [ "$(jobs -p | wc -l)" -ge "$SLURM_NTASKS" ]; do
    wait -n
  done
  SAMPLE=$(basename $FILE .1.fq)
  srun --ntasks=1 bowtie2 --very-sensitive -p $SLURM_CPUS_PER_TASK -x "$INDEX_PREFIX" -1 $FILE -2 $FASTQ_DIR/${SAMPLE}.2.fq -S $data_dir/${SAMPLE}.sam &
done
wait

# Sort
for FILE in $FASTQ_DIR/*1.fq; do
  while [ "$(jobs -p | wc -l)" -ge "$SLURM_NTASKS" ]; do
    wait -n
  done
  SAMPLE=$(basename $FILE .1.fq)
  srun --ntasks=1 samtools sort --threads $SLURM_CPUS_PER_TASK -o $data_dir/${SAMPLE}.sorted.bam -O bam $data_dir/${SAMPLE}.sam &
done
wait

# Filter low quality mapping
for FILE in "$FASTQ_DIR"/*1.fq; do
    while [ "$(jobs -p | wc -l)" -ge "$SLURM_NTASKS" ]; do
        wait -n
    done
    SAMPLE=$(basename "$FILE" .1.fq)
    srun --ntasks=1  samtools view -b -q 20 "$data_dir/${SAMPLE}.sorted.bam" > "$data_dir/${SAMPLE}.filtered.bam" &
done
wait

# Count reads
for FILE in $FASTQ_DIR/*1.fq; do
  while [ "$(jobs -p | wc -l)" -ge "$SLURM_NTASKS" ]; do
    wait -n
  done
  SAMPLE=$(basename $FILE .1.fq)
  srun --ntasks=1 htseq-count -f bam -r pos -t CDS -i gene_id "$data_dir/${SAMPLE}.sorted.bam" "$GTF" > "$data_dir/${SAMPLE}.count"  &
done
wait


# Combine counts
files=("$data_dir"/*.count)
combined_counts="${data_dir}/combined_counts.txt"
first_file="${files[0]}"

cut -f1,2 "$first_file" > "$combined_counts"

for file in "${files[@]:1}"; do
    cut -f2 "$file" > temp_col
    paste "$combined_counts" temp_col > temp_combined
    mv temp_combined "$combined_counts"
done

header="Gene_ID\t$(basename "$first_file" .count)"

for file in "${files[@]:1}"; do
    header="$header\t$(basename "$file" .count)"
done

sed -i "1i$header" "$combined_counts"

rm -f temp_col
rm -f "$data_dir"/*.sam

conda deactivate

