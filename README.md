# MSI_COAD

This repository is the code base to support the submission to Microbiome related to intratumoral microbiome in Colorectal Carcinoma.


The environment required to generate the taxonomy files must have the following bioinformatics packages in the $PATH installed with conda or otherwise:
1) bwa
2) samtools
3) blastn

If using conda the specific versions used for this paper can be installed using the base_environment.yaml included in repository with command:

```bash
conda env create -f base_environment.yaml
```


Also  clone KronaTools repository <https://github.com/marbl/Krona/tree/master/KronaTools> and follow documentation for setting up.  Ensure the that the path to KronaTools is available for below. 

Download the nucleotidate database from ncbi by running the following command

```bash
rsync -av --progress rsync://ftp.ncbi.nlm.nih.gov/blast/db/nt* /path/to/nt/folder/
```

This will download the tar files from the ncbi website.  The then need to be decompressed. 

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


