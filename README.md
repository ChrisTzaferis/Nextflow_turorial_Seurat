# Nextflow turorial Seurat
A nextflow tutorial for scRNA-seq analysis using SeuratV5 (nf-core hackathon March 2025)

In the nextflow.config file you can define the following parameters:
1) <b>input_folder:</b> a path to a directory containing the output of cellranger. It should contain the following files <b><i>barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz or barcodes.tsv, genes.tsv, and matrix.mtx</b></i>.
2) <b>output_name:</b> the name that will be used to store the processed Seurat object.
3) <b>mt_pattern:</b> a pattern for mitochondrial genes. Use "^MT-" for human data or "^mt-" for mouse data.

After successful execution, the processed RDS object is stored in a subfolder within the work directory created during the run.

In order to run the analysis open a terminal and use the following command: <b>nextflow run main.nf -profile conda -resume</b>
