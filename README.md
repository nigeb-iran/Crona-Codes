# Crona-Codes 
This Tutorial is for extracting nucleotide variation from SARS-COV2 genome. The sequence file must be  download from gisaid.org database.

Required OS. and file:

1-Linux Ubuntu

2-bwa

3-samtools

4-Python >= 3.0

Stages:
1- First downlod sequence file from gisaid.org and rename the file to "sequence.fasta".

2- Run the single.py in directory. This script convert a multi-line file to sinle-line fasta file.

3- Download metadata file from gisaid.org and rename it to "metadata".Run the following commands in Terminal to generate a new metadata file.

awk –f “\t” ‘{print $4}’  metadata > id 

awk –f “\t” ‘{print $5}’ metadata >year

awk –f “\t” ‘{print $6}’ metadata >country

paste –d “/” id year country > metadata2

4- Run the convert-id-seq.py script. By running this script, the sequences of each country will be extracted. 

5- Run the following command to write the name of the countries to a file. 

ls . > country

Additional words or characters other than the names of the countries must be removed from the file


7- Run the "final-mutation.sh" script in directory
