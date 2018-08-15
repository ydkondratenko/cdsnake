# cdsnake
Snakemake pipeline using CD-HIT-OTU-MiSeq utilities to extract OTUs from paired-end reads

To use pipeline you need to edit Snakefile according to your starting conditions:

1.	Place your paired-end 16S reads in a folder. Reads should be unzipped prior to pipeline usage (should have an extension “.fastq”). Forward and reverse reads should have “_R1_” and “_R2_” in filenames. Alternatively, filenames for forward and reverse reads can end with “_1” and “_2”, prior to “.fastq” extension.
2.	Specify the name of this folder in INPUT_FOLDER parameter.
3.	Install CD-HIT-OTU-Miseq as described here. Follow steps 1-4 of Installation section.
4.	Specify the path to cd-hit folder in CDHIT parameter.
5.	Specify the name of your 16S database in DATABASE parameter. Note, that database should be formatted as described in 4 step of Installation of CD-HIT-OTU-Miseq (this reformatting places sequences with more specific annotations are at the beginning of the file). You can download reformatted versions of Greengenes and Silva databases here.
6.	Specify a path to trimmomatic.jar file in TRIMMOMATIC parameter.
7.	You can specify p and q parameters to better fit your data. p and q defines length of parts of R1 and R2 reads with a decent quality. For example, you can have the quality of first 220 nucleotides of your R1 reads of about 30 with a sudden drop of quality after 220 nucleotide. In this case you should set value of p parameter equal to 220. You can examine quality of your reads using graphic reports of fastqc or use other method of quality accession. 
8.	Parameters THREADS and MINLEN specify settings for trimmomatic, which is used for data filtering. Trimmomatic will try to use THREADS number of threads and will discard cropped reads if they become shorter than 100 nucleotides.

To run pipeline, you need to install Snakemake. Then place your edited Snakefile in input folder containing reads and cd into this folder. To run pipeline with default Snakemake settings, use command snakemake. To modify general parameters of pipeline (such as number of cores to use) use Snakemake parameters as described here.
Results of clustering for pooled samples will be stored in folder “pooled” in input folder. OTU table for all samples will be saved in file “OTU.txt” in “pooled” folder.
