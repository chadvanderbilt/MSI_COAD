#!/bin/bash

# provide base path when executing script as the first argument
BASE_PATH=$1

echo $BASE_PATH
#Set R environment variables
export LD_LIBRARY_PATH=/data/iacobuzc/vanderbc/anaconda3/lib:$LD_LIBRARY_PATH
export R_LIBS_USER='/lila/data/iacobuzc/vanderbc/anaconda3/lib/R/library'

# Provde path to to executable files
export BIN=/data/vanderbilt/silent_validation/scripts



export GENE_LIST=/data/vanderbilt/silent_validation/manifest/cancerGeneList.tsv
# Script will genererate the following directories:
mkdir -p $BASE_PATH/fastq
mkdir -p $BASE_PATH/BAM
mkdir -p $BASE_PATH/fasta
mkdir -p $BASE_PATH/sub_scripts
mkdir -p $BASE_PATH/stout
mkdir -p $BASE_PATH/blast
mkdir -p $BASE_PATH/tax_out


cd /data/vanderbilt/

# Provide the path with wild cards to the fastq files

for file in "$BASE_PATH"/Project_*/*/*/*_R1_001.fastq.gz; do
   echo $file
   IFS="/" read -ra path_parts <<< "$file"
   
   # Extracting relevant components from the path_parts array

   fastq="${path_parts[-1]}"
   
   # Creating the 'case' variable based on the extracted components
   case="test_${fastq}"
   
   echo $case
   rm "$BASE_PATH"/sub_scripts/"$case".sh
   echo "#!/bin/bash" >> "$BASE_PATH"/sub_scripts/"$case".sh
   echo "#BSUB -J RNA"$case"" >> "$BASE_PATH"/sub_scripts/"$case".sh
   echo "#BSUB -oo "$BASE_PATH"/stout/"$case".st_out" >> "$BASE_PATH"/sub_scripts/"$case".sh
   echo "#BSUB -W 165:00" >> "$BASE_PATH"/sub_scripts/"$case".sh
   echo "#BSUB -R \"rusage[mem=16GB]\"" >> "$BASE_PATH"/sub_scripts/"$case".sh
   echo "#BSUB -q cpuqueue" >> "$BASE_PATH"/sub_scripts/"$case".sh
   echo "#BSUB -n 8" >> "$BASE_PATH"/sub_scripts/"$case".sh
   echo "#BSUB -L /bin/bash" >> "$BASE_PATH"/sub_scripts/"$case".sh
   echo "Ref=/data/vanderbilt/new/mutect_ref/Homo_sapiens_assembly19.fasta" >> "$BASE_PATH"/sub_scripts/"$case".sh
   echo "cd /data/vanderbilt/" >> "$BASE_PATH"/sub_scripts/"$case".sh
   echo "bwa mem \
         -M \
         -t 8 \
         -v 3 \
         \$Ref \
         $file \
         ${file/_R1/_R2} \
         | samtools sort -O BAM \
         | samtools view -f 4 -o "$BASE_PATH"/BAM/"$case".unmapped.bam;" >> "$BASE_PATH"/sub_scripts/"$case".sh
   echo "samtools fasta \
   "$BASE_PATH"/BAM/"$case".unmapped.bam > \
   "$BASE_PATH"/fasta/"$case".fasta;" >> "$BASE_PATH"/sub_scripts/"$case".sh
   echo "blastn \
   -evalue 1e-10 -word_size 28 -db \
   /data/iacobuzc/vanderbc/Krona/KronaTools/ntBlast/ntBlast/nt \
   -query \
   "$BASE_PATH"/fasta/"$case".fasta \
   -outfmt 7 -perc_identity 90 -num_threads 4 > \
   "$BASE_PATH"/blast/"$case".fasta.bl;" >> "$BASE_PATH"/sub_scripts/"$case".sh

   echo "/lila/data/vanderbilt/Krona/KronaTools/scripts/ClassifyBLAST.pl \
   -s "$BASE_PATH"/blast/"$case".fasta.bl \
   -o "$BASE_PATH"/tax_out/"$case".fasta.bl.tax;" >> /"$BASE_PATH"/sub_scripts/"$case".sh
   echo "rm "$BASE_PATH"/blast/"$case".fasta.bl" >> /"$BASE_PATH"/sub_scripts/"$case".sh
   chmod +x "$BASE_PATH"/sub_scripts/"$case".sh
   bsub < "$BASE_PATH"/sub_scripts/"$case".sh

done;