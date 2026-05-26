# Mutant-Calling

This Nextflow pipeline uses identified genes or regions of interest in a reference genome to identify the status of the gene or region in a set of target genome. It then aligns to the reference, and identifies amino acid substitutions, deletions, and insertions to determine the potential functionality of the gene. If enabled, the function of the gene can be assigned and the ancestral state can be determined.

## Dependencies

**Nextflow**  
Install nextflow following the instructions at https://www.nextflow.io/docs/latest/getstarted.html.

**Anaconda**
This pipeline is enabled with conda. Install conda at https://anaconda.org/.

All dependencies must be added to path.

## Installation
**Via github:**  
```bash 
git clone git@github.com:GaTechBrownLab/Mutant-Calling.git
```

**Via nextflow:** 
```bash 
nextflow pull GaTechBrownLab/Mutant-Calling
```

## Usage

**Basic usage:**  
```bash
nextflow run main.nf -with-conda --input_genes ./data/genes/*.fna --input_prots ./data/prot/*.faa
    --input_muts ./data/mutation_patterns/*.txt --host_genomes ./data/hosts/*.fna
    --outdir /results --pres_abs true --graphing false
```

## Pipeline overview


```mermaid
graph TD;
    id1[/"Input genome fasta files (.fna)"/] ==> id1_1[Make blast database];
    id3[/"Input reference genes (.fna)"/] ==> id2;
    subgraph nuclgraph [" "]
        id1_1 ==> id2["Run **blastn**"];
        id2 ==> id6[Identify complete, incomplete, and missing genes];
        id6 ==> id7[For **complete** genes, pull the fasta sequence with **bedtools**];
        id7 ==> id8[Translate nucleotide to protein sequence using **Transeq**]; 
    end
    id4[/"Input reference proteins (.faa)"/] ==> id5[Run **blastp**];
    subgraph protgraph [" "]
        id8 ==> id5;
        id5 ==> id9[Identify WT and mutant proteins];
        id9 ==> id10[Cluster mutant proteins with **cd-hit**];
        id10 ==> id11[Identify a reference sequence per cluster, and align against reference genes with **clustalo**];
        id11 ==> id12[Identify amino acid substitutions, insertions, and deletions with **biopython**]
    end
    id3 ==> id11;
    id12 ==> id13[\"**Final output:** Assigned mutations and functionality to all genes within reference genome set"\]
    id14[/"Input mutation patterns associated with a non-functional protein per gene"/] ==> id13

	classDef nuclsteps fill:#235e8d,font-size:25px,stroke:#000000,color:#FFFFFF;
	class id1_1,id2,id6,id7,id8 nuclsteps;

    classDef protsteps fill:#28885d,font-size:25px,stroke:#000000,color:#FFFFFF;
	class id5,id9,id10,id11,id12 protsteps;

    classDef inputs fill:#132157,font-size:25px,stroke:#000000,color:#FFFFFF;
    class id1,id3,id4,id14 inputs;

    classDef output fill:#06402B,font-size:25px,stroke:#000000,color:#FFFFFF;
    class id13 output;

    classDef nuclgraphcol fill:#91bfe4,font-size:25px,stroke:#000000,color:#FFFFFF;
    class nuclgraph nuclgraphcol;

    classDef protgraphcol fill:#96e0bf,font-size:25px,stroke:#000000,color:#FFFFFF;
    class protgraph protgraphcol;

    linkStyle default stroke:black,stroke-width:2px;

```
