#!/usr/bin/env bash
path=$1

directory=$(dirname "$path")

echo $directory

rm $directory/awk_out_final_full.txt
rm $directory/awk_tmp2.txt
rm $directory/awk_tmp3.txt

awk '{print $0, FILENAME}' \
   $path/*.tax >> \
   $directory/awk_tmp2.txt

cat $directory/awk_tmp2.txt | grep -v "taxID" > $directory/awk_tmp3.txt

while IFS= read -r row; do
   read var1 var2 var3 var5 < <(echo $row);
   file=$(basename "$var5")
   file_rm="${file%.fasta.bl.tax}"
   echo $file_rm
   echo $var1 $var2 $file_rm >> $directory/awk_out_final_full.txt
done < $directory/awk_tmp3.txt
