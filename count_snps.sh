#! /bin/bash 

# Prepare reference file
seqkit fx2tab reference_sequences.fasta | cut -f 2 > ref_seq.tsv

# Calculate the reference length
ref_length=$(awk '{print length}' ref_seq.tsv)

# Input files
input_files=(homabay_CSP_full_protein.fas Kilifi_CSP_full_protein.fas Kisumu_CSP_full_protein.fas)

# Create a for loop for preparing files 
for file in "${input_files[@]}"
do
    # Get the basename of the input file
    base_name=$(basename "${file}" .fas)
    
    # Define the name of the output file 
    output_file="${base_name}_processed.tsv"

    echo "Pre-processing ${file}"

    seqkit fx2tab "${file}" | cut -f 2 | tr -d '*' > "${output_file}"  
    
    echo "Output saved as ${output_file}"

    # Use awk to calculate line lengths in the output file
    awk '{print length}' "${output_file}" | while read -r line_length
    do
        if [ "$ref_length" -ne "$line_length" ]; then
            echo "Files have different line lengths. Skipping ${file}"
            continue 2  # Continue to the next iteration of the outer loop
        fi
    done

    # Concatenate the content of ref_seq.tsv at the beginning of the output file
    cat ref_seq.tsv "${output_file}" > "${output_file}.temp"
    mv "${output_file}.temp" "${output_file}"

    echo "Joined file saved as ${output_file}"
done 

