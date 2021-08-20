Supporting Information for Estimating Cell Type Composition using Isoform Expression One Gene at a Time by Heiling, Wilson, Rashid, Sun, and Ibrahim

Contained in this .zip folder are the code files used to perform the analyses given in the paper, as well as several relevant materials.

00 _prepare_gene_anno:

This folder contains the code needed to create the non-overlapping exon-level gene model information required for the running of the IsoDeconvMM method. The given code steps create the .bed file and the knownIsoforms list object, which are both needed as input into the IsoDeconvMM() function.

01 _Step_1_get_counts:

This code was used to create the Blueprint pure reference exon set count files from the raw Blueprint data.

02 IsoDeconvMM code used in paper:

This folder contains the code for the IsoDeconvMM method. These code files were sourced when running the analyses in the paper. The code to run this method can also be downloaded as a R package from the GitHub repo https://github.com/hheiling/IsoDeconvMM

03 Blueprint in silico analyses:

This folder contains the following code and materials:

Subfolder '01 data': Code examines the Blueprint meta info and identifies the Blueprint pure reference sample files of interest to use in the Blueprint in silico analyses

Subfolder 'Blueprint_Materials': Contains several materials needed to run the analyses.

Subfolder '02 IsoDeconvMM code': Contains the code that split the pure reference samples into appropriate subfolders (pure reference vs mixture creation), created the mixture files, ran isoDetector, and ran the Blueprint in silico analyses.

Subfolder '03 CIBERSORTx code': Contains the code that created the input materials for the CIBERSORTx procedure. This folder contains its own README for more instructions.

The .Rmd file "Blueprint_IsoDetect_Results.Rmd" analyzes the isoDetect results (isoDetect function of the isoform package) and chooses the best transcript clusters to use in the analyses.

The remaining code analyses the results and creates the graphics seen in the paper. 

04 Simulation analyses:

The Rmarkdown documents contain the instructions for how to create the data for the Simulations section of the paper. The folder 'geneModel_code' contains additional code files needed for the simulation of the data, and the folder 'Human_Materials' contains relevant materials needed for the simulation of the data.

Subfolders 'sim1_code', 'sim2_code', and 'sim3_code' contain additional code needed to create the mixture files from the simulated pure reference files and the code needed to run IsoDeconvMM on the simulated data.

05 Additional Web Appendix B analyses:

This folder contains additional code and materials needed for a couple of additional analyses presented in the Web Appendix B. All code files are numbered in the order in which they needed to be run.


