This experiment was performed to compare sensitivities of alignment-independent and dependent methods for calling differentiall expressed genes. Alignment methods used here are:

1. Star (alignment dependent)
2. Kallisto (alignment independent)
3. Salmon (alignment independent) 

Shortly, alignment independent methods work by assigning the  most likely transcript for each read while alignment dependent methods align the reads to a reference genome. 
While alignment independent methods are significantly faster than the alignment dependent methods, alignment dependent methods are necessary for discovering novel features.
For the purposes of this experiment, I focused on both gene level quantification for calculating loss of information through different workflows.

## Data Simulation:

Experimental data was simulated using Polyester R package which uses negative binomial distribution and an input fasta file for transcripts to simulate the dataset. Data was only 
simulated for transcripts on human chromosome 22 and ~10% of the transcripts were set to be differentially expressed between 2 groups where each group contained 8 replicates.
Input data was quantified by summing created reads per gene.

## Evaluation:

Performance of different feature mapping methods were evaluated by comparing percentage of genes that were correctly aligned and identified as differentially expressed.
Absolute relative difference (ARD) and Pearson's correletion coefficient was also caulcuated between the input and output samples. 

## Results

| Method | Mean ARD | R | % of correctly Aligned Genes  | % of correctly identified DEGS|
|--------|----------|---|-------------------------------|-------------------------------|
| First b| Second ce| Th| fourth cell                   |                               |
| Second | foo      | **| baz                           |                               | 
| Third l| quux     | ba| bar                           |                               |
|--------+----------+---+-------------------------------+-------------------------------|

{: .custom-class #custom-id}









