# load the files in an array
input_files=(homabay_CSP_full_protein.fas Kilifi_CSP_full_protein.fas Kisumu_CSP_full_protein.fas)

# Create a lop for analyzing the individual files
for input_file in "${input_files[@]}"
do
    
    # Extracting the basename without the  extension
    file_name=$(basename "${input_file}" .fas)
    
    
    # Creating the output files
    th2r_file="${file_name}_th2r.tsv"
    th3r_file="${file_name}_th3r.tsv"
    
    echo "Extracting the th2r and th3r sequences in ${file_name}"
    
    seqkit fx2tab "${input_file}" | cut -f 2| cut -c310-327 | sort | uniq -c  > "${th2r_file}"
    seqkit fx2tab "${input_file}" | cut -f 2| cut -c352-363 | sort | uniq -c  > "${th3r_file}"
    
done
