#! /bin/bash

# Select nucleotide files of interest

fasta_file=(kisumu_CSP_full_nucleotide.fas kilifi_CSP_full_nucleotide.fas Homabay_full_nucleotide_csp.fas)


# Prepare nucleotides for population genetic analysis. In this case extract CSP C-terminal region from position 909-1140

for i in "${fasta_file[@]}"
do
    
    # Create basename of individual samples"
    base_name=$(basename ${i} .fas )
    
    # create output files
    output_files="${base_name}_csp_c_terminal.fasta"
    temp_seq="${base_name}_seq.tsv"
    temp_header="${base_name}_header.tsv"
    
    
    # Extract the c-terminal region of CSP
    seqkit fx2tab ${i} | tail +2 | cut -f 2 | cut -c909-1140 > "${temp_seq}"
    
    # Extract the gene id from the file
    seqkit fx2tab ${i} | tail +2 | cut -f 1 > "${temp_header}"
    
    # paste the gene_id data to the c-terminal region of CSP
    paste -d "\t" ${temp_header} ${temp_seq} | seqkit tab2fx  > ${output_files}
done

