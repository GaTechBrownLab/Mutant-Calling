#!/bin/bash

# Inputs
input_CCS="flat_df_annotated_Adult_CF.csv"
subset_to_pull="20"
data="data"
outdir="Results"
multifasta_input=""
input="full"

mkdir -p $outdir

# Subset specific genes
if [ $input = "subset" ]; then

    echo "Subsetting genes"

    # Pull just top genes
    python reformat_top_genes.py \
        "$input_CCS" \
        "$subset_to_pull" \
        "$data" \
        "$outdir"

    # Pull top genes from multifasta
    seqkit grep -f ${outdir}/Top_20_genes.txt ${data}/shap_genes_of_interest_blocked_model_prot.fasta > ${outdir}/shap_genes_of_interest_blocked_model_prot_sub.fasta

    multifasta_input="${outdir}/shap_genes_of_interest_blocked_model_prot_sub.fasta"
fi

# Run blast
# PAO1 genome retrieved from pseudomonas db
makeblastdb -in ${data}/Pseudomonas_aeruginosa_PAO1_107_cleaned.faa -parse_seqids -out ${outdir}/PAO1_prot/PAO1_prot -dbtype prot

blastp \
    -query ${multifasta_input} \
    -task blastp \
    -db "${outdir}/PAO1_prot/PAO1_prot" \
    -outfmt "6 qseqid sseqid sacc pident nident qlen length evalue slen qstart qend sstart send" \
    -evalue 0.01 \
    -num_threads 24 \
    >> "${outdir}/PAO1_blast_comparison_new.txt"

# Reformat blast
# PAO1 annotation retrieved from pseudomonas db

python reformat_blast_new.py \
    "${outdir}/PAO1_blast_comparison_new.txt" \
    "$data" \
    "$outdir"

# Pull PAO1 genes that match to locus tags
# PAO1 nucleotide and protein files retrieved from pseudomonas db
mkdir -p ${outdir}/PAO1_nucleotide

while read -r id; do
    seqkit grep -n -p "$id" \
        "${data}/Pseudomonas_aeruginosa_PAO1_107_cleaned.ffn" \
        -o "${outdir}/PAO1_nucleotide/${id}_nucl.fna"
done < "${outdir}/PA_locus_tags.txt"

mkdir -p ${outdir}/PAO1_protein

while read -r id; do
    seqkit grep -n -p "$id" \
        "${data}/Pseudomonas_aeruginosa_PAO1_107_cleaned.faa" \
        -o "${outdir}/PAO1_protein/${id}_prot.faa"
done < "${outdir}/PA_locus_tags.txt"