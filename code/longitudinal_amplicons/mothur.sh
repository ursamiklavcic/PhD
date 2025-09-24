#!/bin/bash


# Script is in majority prepared as recommended in Mothur MiSeq SOP, adjusted for V3V4 region of 16S rRNA.
# All important outputs have prefix final.* , the rest can be discarded after analysis
# Mothur writes large .txt files during the analysis, therefore it is more efficient to run the script on /fast/

# To run this script submit and sbatch job with this script at the end. Example:
# sbatch --time=3-12:00 --mem=200GB --cpus-per-task=84 --job-name=mothur --nodelist=heracpu02.nlzoh.si mothur.sh

eval "$(conda shell.bash hook)"
source activate mothur

mothur

# Set the number of processors, as mothur will use all you've got ! 
set.current(processors=84)

# Stability file, 1st column name of sample, 2nd column forward read, 3rd column reverse read
make.file(inputdir=/volumes/homehpc/storage/finished_projects/ursa/longitudinal_amplicons2_2023/raw_reads, outputdir=/volumes/homehpc/storage/finished_projects/ursa/longitudinal_amplicons2_2023/mothur, type=gz, prefix=stability)
# Makes contigs
make.contigs(file=stability.files) 
set.current(inputdir=/volumes/homehpc/storage/finished_projects/ursa/longitudinal_amplicons2_2023/mothur)

# Sequence data and group identity for each sequence
summary.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table)

# remove sequences that are not good enough
# maxambig = max number of ambiguous base = 0
# minlength = min length of the sequence = 2,5% percentil
# maxlength = 97,5% percentil
# maxhomop = max length of homopolymer = 8
screen.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table, maxambig=0, minlength=439, maxlength=465, maxhomop=8)

# Only align unique sequences 
unique.seqs(fasta=stability.trim.contigs.good.fasta, count=stability.contigs.good.count_table)

# Aligning sequences to the reference database
align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=/volumes/homehpc/storage/DB/silva/silva.V3V4.pcr.align)
summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table)
screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, start=6428, end=23440)

# Remove gaps in the sequences 
filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.)

# to remove any redundancy
unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table)

# Pre-cluster; this is basicly making ASVs
pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=4)

# Remove chimeras
chimera.vsearch(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)

#Classify 
classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.count_table, reference=/volumes/homehpc/storage/DB/RDP/trainset9_032012.pds.fasta, taxonomy=/volumes/homehpc/storage/DB/RDP/trainset9_032012.pds.tax)

# Remove all that is not Bacteria 
remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)

# Rename files into final, so you can find them later! 
rename.file(fasta=current, count=current, taxonomy=current, prefix=final)

# OTU
cluster.split(fasta=final.fasta, count=final.count_table, taxonomy=final.taxonomy, taxlevel=4, cutoff=0.03)
make.shared(list=final.opti_mcc.list, count=final.count_table, label=0.03)
classify.otu(list=final.opti_mcc.list, count=final.count_table, taxonomy=final.taxonomy, label=0.03)

# ASV 
#make.shared(count=final.count_table)
#classify.otu(list=final.asv.list, count=final.count_table, taxonomy=final.taxonomy, label=ASV)
#get.oturep(column=stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.dist, list=final.asv.list, fasta=final.fasta, count=final.count_table) 
