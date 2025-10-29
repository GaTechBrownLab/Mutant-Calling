process clustalo {
    
    publishDir "${params.outdir}/mutant_calling_output/${gene_ID}/alignments/clustalo", mode: 'copy'

    input:
        tuple val(gene_ID), path(input_file), path(ref)

    output:
        tuple val(gene_ID), path("*_aln.fasta")
        
    script:
    """
    new_file=\$(head -1 "$input_file")
    new_file_clean=\${new_file//>}.faa
    mv "${input_file}" "\${new_file_clean}"

    cat "$ref" <(echo) "\${new_file_clean}" > "\${new_file_clean%.*}_PAO1.faa"
    sed -i 's/ //g' "\${new_file_clean%.*}_PAO1.faa"

    clustalo -i \${new_file_clean%.*}_PAO1.faa -o \${new_file_clean%.*}_aln.fasta

    """
}