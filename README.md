# Software for Prox-seq

Prox-seq is a single-cell sequencing assay developed by Tay Lab at the University of Chicago. Prox-seq can be used to obtain gene expression profile, and protein complex information in single cells. This is achieved by leveraging proximity ligation assay (PLA) with single-cell RNA-sequencing techniques.

There are two modes of operation for Prox-seq, one is droplet-based (eg, Drop-seq) and one is plate-based (eg, Smart-seq2). Prox-seq can readily be used with any poly-A-based single-cell RNA sequencing techniques.

## Prox-seq read alignment

A Java program is used to align PLA sequencing reads to a reference list of antibody barcode, and returns a matrix of PLA product counts.

Please refer to [the manual](https://github.com/tay-lab/Prox-seq/blob/master/alignment_manual.pdf) on how to use the alignment program.

## Prox-seq data analysis

Besides mRNA data, Prox-seq also provides counts of PLA products. The PLA count data can be directly analyzed with protein count data analysis pipelines, such as the Seurat package in RStudio. Furthermore, PLA count data can be used to calculate protein abundance (similar to CITE-seq and REAP-seq), and protein complex abundance. The [ProxseqFunctions.py](https://github.com/tay-lab/Prox-seq/blob/master/ProxseqFunctions.py) file contains the functions used for these kinds of calculations. Details about these functions are available below.

Please refer to [this example](https://github.com/tay-lab/Prox-seq/blob/master/PLA_data_analysis_example.ipynb) of how to analyze Prox-seq count data.

## Details ProxseqFunctions.py

### Standard analysis

**calculateProteinAbundance(data, sep=':')**

Calculate protein abundance by summing the counts of each target, regardless of the target of probe A and B.

    Parameters
    ----------
    data : pandas data frame
        Columns are cell barcodes, rows are PLA products
    sep : string, optional
        The separator format for PLA product.
        Default is ':'.

    Returns
    -------
    Returns a pandas data frame.
    Each row is the total abundance of proteins 1, 2,...

**calculateProbeAbundance(data, sep=':')**

Calculate probe abundance by summing the counts of probes A and B of each protein target.

    Parameters
    ----------
    data : pandas data frame
        Columns are cell barcodes, rows are PLA products
    sep : string, optional
        The separator format for PLA product.
        Default is ':'.

    Returns
    -------
    Returns a pandas data frame.
    Each row is the total abundance of probe A1, A2,... and B1, B2,...

**calculateExpected(data, PLA_list=None, sep=':')**

Calculate the expected random count of a PLA product, if no protein interactions exist in the data.

    Parameters
    ----------
    data : pandas data frame
        Input digital PLA count matrix.
    PLA_list: list, optional
        List of PLA products for which expected count is calculated.
        If None (the default), calculate expected count for all PLA products.

    Returns
    -------
    A data frame of expected count (rows = PLA, columns = single cells).

### Protein complex prediction

**estimateComplexes(data, non_complex=[], mean_cutoff=1, p_cutoff=0.05, p_adjust=True, sym_weight=0.25, df_guess=None, nIter=200, tol=5, sep=':')**

Estimate complex abundance by iteratively solving a system of quadratic equations. The system of equations is set up based on the expected random count of each PLA product.

    Parameters
    ----------
    data : pandas data frame
        Input digital PLA expression matrix (PLA products x single cells).
    non_complex : list
        List of PLA products or proteins that do no form protein complexes.
        Example: X:Y means X:Y does not form a complex, while X means X does
        not form complexes with any other proteins.
        Default is [].
    non_express : list
        List of protein targets that do not form complexes (eg, isotype antibodies).
        Default is [].
    mean_cutoff : float
        PLA products whose estimated complex abundance at each iteration fails
        the 1-sided t-test sample mean > mean_cutoff is kept as 0.
        Default is 1.
    p_cutoff : float
        The alpha level to decide if the 1-sided t-test is sinificant.
        Default is 0.05.
    p_adjust : boolean
        Whether to perform FDR correction for the one-sided t-test.
        Default is True.
    sym_weight : float (0 <= sym_weight <= 1).
        The weight factor used to enforce symmetry condition. 0 means no enforcement.
        Default is 0.25.
    df_guess : pandas data frame
        First guesses of true complex abundance (must be the same shape as data).
        If None (the default), use 0 as the first guess.
    nIter : int
        Max number of iterations to perform.
    tol : float
        If the change in solution between current and last iteration is below
        this value, convergence is reached.
        Default is 5.
    sep : string, optional
        The separator convention in the names of PLA complexes.
        Default is ':'.

    Returns
    -------
    A data frame with the same shape as df, containing estimated complex abundance
