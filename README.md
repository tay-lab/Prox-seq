# Software for Prox-seq

Prox-seq is a single-cell sequencing assay developed by Tay Lab at the University of Chicago. Prox-seq can be used to obtain gene expression profile and protein proximity information from single cells. There are two modes of operation for Prox-seq, one is droplet-based (eg, Drop-seq) Drop-seq and one is plate-based (eg, Smart-seq2). Prox-seq can readily be used with any poly-A-based single-cell RNA sequencing techniques.

#### PLA alignment

A Java program is used to align PLA sequencing reads to a reference list of antibody barcode, and returns a matrix of PLA product counts.

Examples of how to use the program to align sequencing data from Drop-seq and Smart-seq2 pipelines are provided.

#### PLA product analysis

The ProxseqFunctions.py file contains functions used to estimate protein complex abundance.

Please refer to [this example](https://github.com/tay-lab/Prox-seq/blob/master/PLA_data_analysis_example.ipynb) of how to predict the protein complex abundance from the Prox-seq count data.
