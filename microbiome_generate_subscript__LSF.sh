#!/bin/bash

# provide base path when executing script as the first argument
BASE_PATH=$1

echo $BASE_PATH

# Provde path to to executable KronaTools files
export KT_BIN=/path/to/Krona/KronaTools
export REF_GEN=/path/to/human/reference/genome.fasta
export NT_DB=/path/to/nt/folder/


# Script will genererate the following directories:
mkdir -p $BASE_PATH/fastq
mkdir -p $BASE_PATH/BAM
mkdir -p $BASE_PATH/fasta
mkdir -p $BASE_PATH/sub_scripts
mkdir -p $BASE_PATH/stout
mkdir -p $BASE_PATH/blast
mkdir -p $BASE_PATH/tax_out


cd $BASE_PATH

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
   echo "Ref="$REF_GEN"" >> "$BASE_PATH"/sub_scripts/"$case".sh
   echo "cd "$BASE_PATH"" >> "$BASE_PATH"/sub_scripts/"$case".sh
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
   "$NT_DB"/nt \
   -query \
   "$BASE_PATH"/fasta/"$case".fasta \
   -outfmt 7 -perc_identity 90 -num_threads 4 > \
   "$BASE_PATH"/blast/"$case".fasta.bl;" >> "$BASE_PATH"/sub_scripts/"$case".sh

   echo ""$KT_BIN"/scripts/ClassifyBLAST.pl \
   -s "$BASE_PATH"/blast/"$case".fasta.bl \
   -o "$BASE_PATH"/tax_out/"$case".fasta.bl.tax;" >> /"$BASE_PATH"/sub_scripts/"$case".sh
   echo "rm "$BASE_PATH"/blast/"$case".fasta.bl" >> /"$BASE_PATH"/sub_scripts/"$case".sh
   chmod +x "$BASE_PATH"/sub_scripts/"$case".sh
   bsub < "$BASE_PATH"/sub_scripts/"$case".sh

done;
