#! /bin/bash

# Select nucleotide files of interest

protein_file=(homabay_CSP_full_protein.fas Kilifi_CSP_full_protein.fas Kisumu_CSP_full_protein.fas)


# Prepare nucleotides for population genetic analysis. In this case extract CSP C-terminal region from position 909-1140

for i in "${protein_file[@]}"
do
    
    # Create basename of individual samples"
    base_name=$(basename ${i} .fas )
    
    # create output files
    output_files="${base_name}_csp_crr.fasta"
    temp_seq="${base_name}_seq_crr.tsv"
    temp_header="${base_name}_header_crr.tsv"
    
    
    # Extract the c-terminal region of CSP
    seqkit fx2tab ${i} |  cut -f 2 | cut -c105-272 > "${temp_seq}"
    
    # Extract the gene id from the file
    seqkit fx2tab ${i} |  cut -f 1 > "${temp_header}"
    
    # paste the gene_id data to the c-terminal region of CSP
    paste -d "\t" ${temp_header} ${temp_seq} | seqkit tab2fx  > ${output_files}
done

