# Mutant-Calling

This Nextflow pipeline uses identified genes or regions of interest in a reference genome to identify the status of the gene or region in a set of target genome. It then aligns to the reference, and identifies amino acid substitutions, deletions, and insertions to determine the potential functionality of the gene. If enabled, the function of the gene can be assigned and the ancestral state can be determined.

## Dependencies

**Nextflow**  
Install nextflow following the instructions at https://www.nextflow.io/docs/latest/getstarted.html.

**Micromamba**
This pipeline is enabled with micromamba. Install micromamba at https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html.

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
**Inputs**
Input gene fna files must be in the output {gene}_nucl.fna, and the header must be in >{gene} format. Similarily, all input proteins must be in {gene}_prot.faa, with the header >{gene}. Host genomes must have a .fna extension. Mutation patterns must be in the [Original AA][Position][New AA] format. If the second AA is variable, write it as [Original AA][Position].

All example test inputs are stored in data, and the test host genomes are stored in a zip file that must be unzipped prior to running.

**Basic usage:**  
```bash
nextflow run main.nf -with-conda --input_genes "data/genes/*.fna" --input_prots "data/prot/*.faa" --input_muts "data/mutation_patterns/*.txt" --host_genomes "data/hosts/*.fna" --pres_abs true --graphing false
```

**Outputs**
Outputs with information on the status of each gene (functional, non-functional, or alternative function) are in the final_outputs folder. All other intermediate outputs can be found in the mutant_calling_output folder.


## Pipeline overview


```mermaid
graph TD;
    id1[/"`Input genome fasta files (.fna)`"/] ==> id1_1["`Make blast database`"];
    id3[/"`Input reference genes (.fna)`"/] ==> id2;
    subgraph nuclgraph [" "]
        id1_1 ==> id2["`Run **BLASTn**`"];
        id2 ==> id6["`Identify complete, incomplete, and missing genes`"];
        id6 ==> id7["`For **complete** genes, extract the fasta sequence (**bedtools**)`"];
        id7 ==> id8["`Translate nucleotide sequence (**transeq**)`"]; 
    end
    id4[/"`Input reference proteins (.faa)`"/] ==> id5["`Run **BLASTp**`"];
    subgraph protgraph [" "]
        id8 ==> id5;
        id5 ==> id9["`Identify WT and mutant proteins`"];
        id9 ==> id10["`Cluster mutant proteins (**cd-hit**)`"];
        id10 ==> id11["`Isolate the reference sequence per cluster, and align against reference genes (**clustalo**)`"];
        id11 ==> id12["`Identify amino acid substitutions, insertions, and deletions (**biopython**)`"]
    end
    id4 ==> id11;
    id12 ==> id13[\"`**Final output:** Assigned mutations and functionality to all genes within reference genome set`"\]
    id14[/"`Input mutation patterns associated with a non-functional protein`"/] ==> id13

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
