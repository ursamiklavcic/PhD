# Base directory
base_dir="/home/storage/finished_projects/UrsaM/spore_evo_dynamics/PStrain"

# Output file
output="/home/nlzoh.si/ursmik1/projects/thesis/data/PStrain/PStrain_all.txt"

# Write header once
echo -e "Species\trel_abund_species\tStrain\tfreq_strain\trel_abund_strain\tsample" > "$output"

# Loop through all subdirectories inside base_dir
for sample_dir in "$base_dir"/*/; do
  # Remove trailing slash and extract just the folder name (sample ID)
  sample_id=$(basename "$sample_dir")
  
  # Path to the file
  file_path="$sample_dir/result/strain_RA.txt"
  
  # Check if the file exists
  if [[ -f "$file_path" ]]; then
    # Append file contents to output adding the sample_id as last column
    awk -v sid="$sample_id" '{
      # Add sample ID as last column for each line
      print $0 "\t" sid
    }' "$file_path" >> "$output"
  else
    echo "Warning: File not found $file_path"
  fi
done
sed '1! s/ \+/\t/g' "$output" >  /home/nlzoh.si/ursmik1/projects/thesis/data/PStrain/PStrain_all.tsv
