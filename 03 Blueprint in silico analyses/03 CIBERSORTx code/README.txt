Order of Operations of CIBERSORTx code:

(1) cibersort_pure_matrix.R

Odd CT Out Comparisons Select FC folder:

(2) cibersort_DESeq2.R

(3) cibersort_mix_matrix.R, cibersort_pure_matrix_partial.R (can be run simultaneously)

Folders:

Odd CT Out Comparisons Select FC: DESeq2 code compares each cell type with the collective two other cell types (1 v 2 and 3, 2 v 1 and 3, 3 v 1 and 2), restricts the results to genes that have a log2 fold change > 0 (clusters upregulated in one cell type compared to the other two cell types), and picks the 'best' clusters from these results by choosing the clusters with the lowest p-values from each comparison.