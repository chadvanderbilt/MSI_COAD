# MSI_COAD

This repository is the code base to support the submission to Microbiome journal related to intratumoral microbiome in Colorectal Carcinoma, titled "The Microbiome of Microsatellite Instable Colorectal Carcinoma".


The environment required to generate the taxonomy files must have the following bioinformatics packages in the $PATH installed with conda or otherwise:
1) bwa
2) samtools
3) blastn

If using conda the specific versions used for this paper can be installed using the base_environment.yaml included in repository with command:

```bash
conda env create -f base_environment.yaml
```


Also  clone KronaTools repository <https://github.com/marbl/Krona/tree/master/KronaTools> and follow documentation for setting up.  Ensure the that the path to KronaTools is available for below. 

Download the nucleotide database from NCBI for blast queries. 


```bash
rsync -av --progress rsync://ftp.ncbi.nlm.nih.gov/blast/db/nt* /path/to/nt/folder/
```

This will download the tar files from the ncbi website.  The files then need to be decompressed which will generate a summary file titled nt. 

```bash
cd /path/to/nt/folder/
for file in *.gz
do
tar -zxvpf "$file"
rm "$file"
done
```

microbiome_generate_subscript__LSF.sh is a shell script for running microbiome analysis on paired end fastq files.  The script generates LSF submission scripts for running each sample in parallel.  

The first section of the file requires modification of the paths to the KronoTools directory, the downloaded nt database, and the human reference genome. 

The script will generate the needed directories for intermediate files and files associated with running each individual job. 

In line 26, the path to the fastq.gz files must be specified with appropriate wild cards to access all directories with fastq files. Ensure the that the fastq files selected are only the forward reads (typically R1 is in the file name). 

To run execute 
```bash
cd ./MSI_COAD
bash microbiome_generate_subscript__LSF.sh /home/user/microbiome_directory
```

This script will generate and submit jobs to LSF.  Upon completion of jobs, the /home/user/microbiome_directory/tax_out will contain microbiome taxonomy files named in same fashion as the fastq.gz input files. The intermediate blast files are deleted to save memory.  If user would like to keep the blast files, line 74 can be commented out.  

The awk.sh script consolidates the taxonomy files to a single file where each row contains the file nam and can be executed as follows:

```bash
cd ./MSI_COAD
bash awk.sh /home/user/microbiome_directory/tax_out
```

This will generate a consolidated taxonomy file in /home/user/microbiome_directory/. 

The final script ran takes the consolidated taxonomy file in and generate databases of species and genus for all of the cases with the readcount of each genus and species respectively tallied.   The script will generate these files as outputs and will also perform two class comparison alpha diverisity and enrichment analysis based on a manifest file that provides the case id (column  named 'DMP_ASSAY_ID') and additional binary column. The example provided has one collumn will be called MSI_H.  The manifest file can contain any number of classes and enrichment analysis and alpha diversity files will be generated for each comparison.

All R scripts were developed with R version 4.2.0.  The shebang points to the R script which should be modified to match the users environment. To ensure all package dependencies are met in users environment, please run R_packages_install.R before executint the main script. 
