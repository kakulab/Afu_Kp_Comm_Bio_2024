Test your NGS analysis pipeline first on a small input file.
To create a small file out of one of your samples do this:

samtools view -bs 0.1 name_of_your_file.bam > small_file.bam       # randomly subsamples 10% of the reads from the input bam file

For the analysis you have to put your raw bam file(s), a genome file (.fasta) and a .gff file into a folder, which represents your working directory.
Then adjust the working directory in your analysis script. (WKDIR=path_to_your_working_directory)

To execute your analysis script you have to copy it into your home directory first.
Then execute your script using this command:

~/name_of_your_script.sh


If the analysis script cannot be executed, make it executable with the following command and try again.

chmod 777 ~/name_of_your_script.sh

Always check some aligned samples in the IGV and compare to the counts at the end.
To look at the aligned and sorted bam files in the IGV you have to index them with the following command:

samtools index name_of_file.bam.sam.bam.sorted.bam

This creates an index file (.bai), which has to be kept in the same folder as the bam file.


For the edgeR analysis you have to use the 2 scripts; edgeR_script1.R and edgeR_script2.R
The first one you don't have to modify if you have 3 replicates and if you want to filter out genes with less than 1 counts per million (see script). 
In the second you have to change the sample and group names where indicated.
You also need your count files and a Target.txt file containing sample info for the two conditions you want to compare. Just modify the existing Targets.txt file.
