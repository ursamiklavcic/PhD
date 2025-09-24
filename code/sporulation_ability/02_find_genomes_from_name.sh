# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate datasets

while IFS=$'\t' read -r clade_name species_name || [ -n "$species_name" ]; do
  echo "Processing: $species_name" >&2
  accessions=$(datasets summary genome taxon "$species_name" --reference --limit 1 < /dev/null | jq -r '
    if (.reports | length) == 0 then
      "NA"
    else
      [.reports[].accession] | join("\t")
    end
  ')
  echo -e "${clade_name}\t${species_name}\t${accessions}"
done < ~/projects/thesis/data/sporulation_ability/spore_determine.tsv > mpa_accession.tsv


# For GTDB 
#while read -r clade_name || [ -n "$clade_name" ]; do
#  echo "Processing: $clade_name" >&2
#  accessions=$(datasets summary genome taxon "$clade_name" --reference --limit 1 < /dev/null | jq -r '
#    if (.reports | length) == 0 then
#      "NA"
#    else
#      [.reports[].accession] | join("\t")
#    end
#  ')
#  echo -e "${clade_name}\t${accessions}"
#done < ~/projects/thesis/data/sporulation_ability/spore_determine_gtdb.tsv > gtdb_accessions.tsv

