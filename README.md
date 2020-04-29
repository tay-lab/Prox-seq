# Software for Prox-seq

Prox-seq is a single-cell sequencing assay developed by Tay Lab at the University of Chicago. Prox-seq can be used to obtain gene expression profile and protein proximity information from single cells. There are two modes of operation for Prox-seq, one based on Drop-seq and one based on Smart-seq2. However, Prox-seq can be used with any poly-A-based single-cell RNA sequencing techniques.

#### PLA alignment

A Java program is used to align PLA sequencing reads to a reference list of antibody barcode, and returns a matrix of PLA product counts.

Examples of how to use the program to align sequencing data from Drop-seq and Smart-seq2 pipelines are provided.

#### PLA product analysis

The ProxseqFunctions.py file contains functions used to estimate protein complex abundance.

An example of how to use the functions to simulate PLA data, and estimate the protein complex abundance from it is provided.
