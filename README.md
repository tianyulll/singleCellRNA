# singleCellRNA

I used these scripts to analyze a large-scale single-cell RNA seq dataset.
For specific project reason, it was run on both raw and filtered counts.
Shell scripts used to accommadate memory requirements for raw data analysis.

Workflow:
- Count with cellranger
- Preprocess with normalization, PCA, selecting highly variable genes
- Integrate 16 datasets based on highly variable genes using rPCA
- Dimensional reduction
- Identify cells with LTR readthrough
- Compare distance in PC space between identified cells
