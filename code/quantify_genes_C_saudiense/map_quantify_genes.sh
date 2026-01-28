#!/bin/bash
#SBATCH --job-name=mapping
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=1-00:00:00
#SBATCH --output=map.out

# ------------------------------- User input -----------------------------------------
GENE_ID="/home/storage/finished_projects/LaraB/rezultati/Quantify_genes_Ursa/SCRIPTS/id.txt"
FFN="/home/nlzoh.si/ursmik1/projects/C_saudiense/004_prokka/*/*.ffn"
FASTQ_DIR="/home/storage/finished_projects/LaraB/rezultati/Quantify_genes_Ursa/FASTQ"
OUTPUT="/home/storage/finished_projects/LaraB/rezultati/Quantify_genes_Ursa/RESULTS"
# ------------------------------- Python scripts -------------------------------------
find_fasta="/home/storage/finished_projects/LaraB/rezultati/Quantify_genes_Ursa/SCRIPTS/find_fasta.py"
tpm_table="/home/storage/finished_projects/LaraB/rezultati/Quantify_genes_Ursa/SCRIPTS/tpm_table.py"
# ------------------------------------------------------------------------------------

set -e
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /home/nlzoh.si/larbez1/miniconda3/envs/qgene_env

# --------------------- Generate fasta, gtf, fasta index, gene_lengths ----------------
pomozna_mapa=$OUTPUT/"index"
podatki=$OUTPUT/"data"
mkdir -p $pomozna_mapa $podatki
FASTA=$pomozna_mapa/"all_filtered_genes.fasta"
GTF=$pomozna_mapa/"all_filtered_genes.gtf"
GL=$pomozna_mapa/"gene_lengths"

# Find fasta genes
if ! ls $FASTA >/dev/null 2>&1; then
  sed -i 's/[\r]//g' $GENE_ID
  $find_fasta $GENE_ID "$FFN" $FASTA
fi

# Create gtf from fasta
if ! ls $GTF >/dev/null 2>&1; then
  awk 'BEGIN{OFS="\t"} /^>/{header=$0; gsub(/^>/,"",header); split(header,a," "); gene_id=a[1]; seqid=gene_id; getline seq; seqlen=length(seq); print seqid,"PROKKA","CDS",1,seqlen,".","+",".","gene_id \""gene_id"\""}' $FASTA > $GTF
fi

# Determine gene lengths
if ! ls $GL >/dev/null 2>&1; then
  cut -f4,5,9 $GTF | sed 's/gene_id //g' | gawk '{print $3, $2-$1+1}' | tr ' ' '\t' > $GL
fi

# Index fasta
if ! ls $pomozna_mapa/*.bt2 >/dev/null 2>&1; then
  fasta_name=$(basename $FASTA .fasta)
  bowtie2-build $FASTA $pomozna_mapa/$fasta_name
fi
# ------------------------------------------------------------------------------------
# ------------------------------- Mapping --------------------------------------------
# Map
for FILE in $FASTQ_DIR/*1.fq; do
  while [ "$(jobs -p | wc -l)" -ge "$SLURM_NTASKS" ]; do
    wait -n
  done
  SAMPLE=$(basename $FILE .1.fq)
  srun --ntasks=1 bowtie2 -p $SLURM_CPUS_PER_TASK -x $pomozna_mapa/$fasta_name -1 $FILE -2 $FASTQ_DIR/${SAMPLE}.2.fq -S $podatki/${SAMPLE}.map.sam &
done
wait

# Sort
for FILE in $FASTQ_DIR/*1.fq; do
  while [ "$(jobs -p | wc -l)" -ge "$SLURM_NTASKS" ]; do
    wait -n
  done
  SAMPLE=$(basename $FILE .1.fq)
  srun --ntasks=1 samtools sort --threads $SLURM_CPUS_PER_TASK -o $podatki/${SAMPLE}.map.sorted.bam -O bam $podatki/${SAMPLE}.map.sam &
done
wait

# Add read groups for picard to work
for FILE in $FASTQ_DIR/*1.fq; do
  while [ "$(jobs -p | wc -l)" -ge "$SLURM_NTASKS" ]; do
    wait -n
  done
  SAMPLE=$(basename $FILE .1.fq)
  srun --ntasks=1 picard AddOrReplaceReadGroups -I $podatki/${SAMPLE}.map.sorted.bam -O $podatki/${SAMPLE}.map.sorted.RG.bam -RGID $SAMPLE -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM $SAMPLE &
done
wait

# Mark duplicates
for FILE in $FASTQ_DIR/*1.fq; do
  while [ "$(jobs -p | wc -l)" -ge "$SLURM_NTASKS" ]; do
    wait -n
  done
  SAMPLE=$(basename $FILE .1.fq)
  srun --ntasks=1 picard MarkDuplicates -I $podatki/${SAMPLE}.map.sorted.RG.bam -O $podatki/${SAMPLE}.map.markdup.bam -M $podatki/${SAMPLE}.map.markdup.metrics -AS true --VALIDATION_STRINGENCY LENIENT -MAX_FILE_HANDLES 1000 --REMOVE_DUPLICATES true &
done
wait

# Count
for FILE in $FASTQ_DIR/*1.fq; do
  while [ "$(jobs -p | wc -l)" -ge "$SLURM_NTASKS" ]; do
    wait -n
  done
  SAMPLE=$(basename $FILE .1.fq)
  srun --ntasks=1 htseq-count -r pos -t CDS -f bam $podatki/${SAMPLE}.map.markdup.bam $GTF > $podatki/${SAMPLE}.count &
done
wait

# TPM
for FILE in $FASTQ_DIR/*1.fq; do
  SAMPLE=$(basename $FILE .1.fq)    
  $tpm_table -n $SAMPLE -c $podatki/${SAMPLE}.count -i <(echo -e "$SAMPLE\t100") -l $GL > $podatki/${SAMPLE}.tpm
done
wait
# ------------------------------------------------------------------------------------
# --------------------------- Combine TPM results ------------------------------------

OUT=$OUTPUT/"combine_tpm.txt"
FS="\t"

> $OUT
first=true
for file in $podatki/*.tpm; do
    if $first; then
        awk -v FS="$FS" -v OFS="$FS" '{print $1, $2}' "$file" > $OUTPUT/tmp_main
        first=false
    else
        awk -v FS="$FS" -v OFS="$FS" '{print $2}' "$file" > $OUTPUT/tmp_col
        paste -d "$FS" $OUTPUT/tmp_main $OUTPUT/tmp_col > $OUTPUT/tmp_new
        mv $OUTPUT/tmp_new $OUTPUT/tmp_main
    fi
done
# Add gene names
awk -v FS="$FS" -v OFS="$FS" 'NR==FNR{map[$1]=$2;next}{print map[$1],$0}' $GENE_ID $OUTPUT/tmp_main > $OUT
rm $podatki/*.sam $podatki/*.bam $OUTPUT/tmp_main $OUTPUT/tmp_col

# Combine counts
files=($podatki/*.count)
output=$OUTPUT/"combined_counts.txt"
first_file="${files[0]}"
paste_header=$(head -n1 "$first_file")
cut -f1,2 "$first_file" > "$output"

for file in "${files[@]:1}"; do
    cut -f2 "$file" > temp_col
    paste "$output" temp_col > temp_combined
    mv temp_combined "$output"
done
header="Gene_ID\t$(basename "$first_file" .count)"
for file in "${files[@]:1}"; do
    header="$header\t$(basename "$file" .count)"
done
sed -i "1i$header" "$output"
rm -f temp_col

conda deactivate
