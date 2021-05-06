# Software for Prox-seq

Prox-seq is a single-cell sequencing assay developed by Tay Lab at the University of Chicago. Prox-seq can be used to obtain gene expression profile, and protein complex information in single cells. This is achieved by leveraging proximity ligation assay (PLA) with single-cell RNA-sequencing techniques.

There are two modes of operation for Prox-seq, one is droplet-based (eg, Drop-seq) and one is plate-based (eg, Smart-seq2). Prox-seq can readily be used with any poly-A-based single-cell RNA sequencing techniques.

#### Prox-seq read alignment

A Java program is used to align PLA sequencing reads to a reference list of antibody barcode, and returns a matrix of PLA product counts.

Examples of how to use the program to align sequencing data from Drop-seq and droplet-based pipelines are provided.

#### Prox-seq data analysis

Besides mRNA data, Prox-seq also provides counts of PLA products. The PLA count data can be directly analyzed with protein count data analysis pipelines, such as the Seurat package in RStudio. Furthermore, PLA count data can be used to calculate protein abundance (similar to CITE-seq and REAP-seq), and protein complex abundance. The ProxseqFunctions.py file contains the functions used for these kinds of calculations.

Please refer to [this example](https://github.com/tay-lab/Prox-seq/blob/master/PLA_data_analysis_example.ipynb) of how to analyze Prox-seq count data.
