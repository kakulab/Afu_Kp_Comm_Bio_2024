# Aspergillus fumigatus and Klebsiella pneumoniae interaction
Dual-RNAseq Data Analysis workflow for Aspergillus fumigatus and Klebsiella pneumoniae Co-culture

This repository contains a pipeline for the primary analysis of Illumina short-read RNAseq data obtained from the co-culture of Aspergillus fumigatus and Klebsiella pneumoniae. It includes the retrieval of genomic data required for analysis from the NCBI, quality control of the raw data, trimming and mapping of reads, removal of reads mapping to ribosomal RNAs and counting the number of reads mapping to genes. Only the RNAseq raw data in .bam or .fastq format (compressed or uncompressed) have to be provided by the user.

## Tools required for analysis:
samtools (http://www.htslib.org/)

bedtools (https://bedtools.readthedocs.io/en/latest/)

FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

MultiQC (https://multiqc.info/)

cutadapt (https://cutadapt.readthedocs.io/en/stable/)

NextGenMap (https://github.com/Cibiv/NextGenMap/wiki)

DeepTools (https://deeptools.readthedocs.io/en/develop/)

HTSeq (https://htseq.readthedocs.io/en/release_0.11.1/#)

All the above-mentioned tools have to be included in yout PATH environment.

# Usage:
Clone the repository by typing: "https://github.com/kakulab/Afu_Kp_Comm_Bio_2024.git" and copy the raw data into the main directory.

If necessary, modify the adapter sequence for read trimming in the config_file.txt. By default, the Illumina TruSeq adapter sequence is used. For other adapter sequences, refer to the following resource: https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences-1000000002694-11.pdf

Run the analysis script by typing the following command: bash analysis_script.sh

After the pipeline has finished change into the diff_expr_analysis directory and use the edgeR_script.R script as a basis for differential expression analysis in R.

