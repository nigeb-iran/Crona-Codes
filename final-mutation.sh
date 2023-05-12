#!/bin/bash
if [ -e country-mutation ]; then
rm -r country-mutation
fi
if [ ! -e hasrak ]; then
mkdir hasrak
fi
mkdir country-mutation
for f in $(cat country);
do 
mkdir country-mutation/$f
if [ ! -e hasrak/per-country ]; then
mkdir hasrak/per-country
fi
if [ ! -e hasrak/per-country/"$f" ]; then
mkdir hasrak/per-country/"$f"
fi
awk -F '\t' '{ if ($8 !~ "frameshift" && $8 !~ "deletion" && $6 !~ "intergeniregion" && $6 !="") {print $0}}' $f/"$f".fasta.vcf.ann.final > country-mutation/"$f".fasta.vcf.ann.final

awk -F "\t" 'length($4) >= 1' country-mutation/"$f".fasta.vcf.ann.final > country-mutation/"$f"_all.mutation_sorted
sed -i 's/\t\t/\t/g' country-mutation/"$f"_all.mutation_sorted
 awk -F '\t' 'BEGIN { OFS = "\t" } { if( $2 <= 805 ) $3 = "NSP1"; \
                                     if( $2 >  805 && $2 <= 2719  ) $3 = "NSP2" ; \
                                     if( $2 >  2719 && $2 <= 8554  ) $3 = "NSP3" ; \
                                     if( $2 >  8554 && $2 <= 10054  ) $3 = "NSP4" ; \
                                     if( $2 >  10054 && $2 <= 10972  ) $3 = "3Clikeproteinase" ; \
                                     if( $2 >  10972 && $2 <= 11842  ) $3 = "nsp6" ; \
                                     if( $2 >  11842 && $2 <= 12091  ) $3 = "nsp7" ; \
                                     if( $2 >  12091 && $2 <= 12685  ) $3 = "nsp8" ; \
                                     if( $2 >  12685 && $2 <= 13024  ) $3 = "nsp9" ; \
                                     if( $2 >  13024 && $2 <= 13441  ) $3 = "nsp10" ; \
                                     if( $2 >  13441 && $2 <= 16237  ) $3 = "Rdrp" ; \
                                     if( $2 >  16237 && $2 <= 18039  ) $3 = "helicase" ; \
                                     if( $2 >  18039 && $2 <= 19620  ) $3 = "exonuclease" ; \
                                     if( $2 >  19620 && $2 <= 20658  ) $3 = "endoRNAse" ; \
                                     if( $2 >  20658 && $2 <= 21552  ) $3 = "methyltransferase" ; \
                                     if( $2 >  21563 && $2 <= 25384  ) $3 = "Spike" ; \
                                     if( $2 >  25393 && $2 <= 26220  ) $3 = "ORF3aprotein" ; \
                                     if( $2 >  26245 && $2 <= 26472  ) $3 = "Envelopprotein" ; \
                                     if( $2 >  26523 && $2 <= 27191  ) $3 = "MemberaneGlycoprotein" ; \
                                     if( $2 >  27202 && $2 <= 27387  ) $3 = "ORF6protein" ; \
                                     if( $2 >  27394 && $2 <= 27759  ) $3 = "ORF7aprotein" ; \
                                     if( $2 >  27894 && $2 <= 28259  ) $3 = "ORF8aprotein" ; \
                                     if( $2 >  28274 && $2 <= 29533  ) $3 = "Nucleocapsid" ; \
                                     if( $2 >  29558 && $2 <= 29674  ) $3 = "ORF10protein" ; \
$1="'$f'" ;print $0  }' country-mutation/"$f"_all.mutation_sorted > country-mutation/"$f"_all.mutation_sorted_lable.l
awk -F "\t" 'BEGIN{OFS ="\t"} length($4) == 1 {print $0}' country-mutation/"$f"_all.mutation_sorted_lable.l > country-mutation/"$f"_all.mutation_sorted_lable.e1 

#awk 'BEGIN {FS=OFS="\t"}{gsub("-","_",$1)}1' country-mutation/"$f"_all.mutation_sorted_lable.e1 > country-mutation/"$f"_all.mutation_sorted_lable.e100 && mv country-mutation/"$f"_all.mutation_sorted_lable.e100 country-mutation/"$f"_all.mutation_sorted_lable.e1

sed -i 's/,/\t/g' country-mutation/"$f"_all.mutation_sorted_lable.e1
awk -F "\t" '{print $1 "\t" $3 "\t" "\t" $5 "\t" $6 "\t" $2 "-"$7 "\t" $8"\t"$9 }' country-mutation/"$f"_all.mutation_sorted_lable.e1 > country-mutation/"$f"_all.mutation_sorted_lable.e12 

awk -F "\t" 'BEGIN{OFS ="\t"} length($4) == 3 {print $0}' country-mutation/"$f"_all.mutation_sorted_lable.l > country-mutation/"$f"_all.mutation_sorted_lable.e3

awk 'BEGIN {FS=OFS="\t"}{gsub("-","_",$1)}1' country-mutation/"$f"_all.mutation_sorted_lable.e3 > country-mutation/"$f"_all.mutation_sorted_lable.e100 && mv country-mutation/"$f"_all.mutation_sorted_lable.e100 country-mutation/"$f"_all.mutation_sorted_lable.e3

sed -i 's/,/\t/g;s/-/\t/g' country-mutation/"$f"_all.mutation_sorted_lable.e3
awk -F "\t" '{print $1 "\t" $3 "\t" $6 "\t" $7 "\t" $2 "-" $9 "\t" $10 "\t" $11 "\n"  $1 "\t" $3 "\t" $6 "\t" $8 "\t" $2 "-" $12 "\t" $13 "\t" $14}' country-mutation/"$f"_all.mutation_sorted_lable.e3 > country-mutation/"$f"_all.mutation_sorted_lable.e32

awk -F "\t" 'BEGIN{OFS ="\t"} length($4) == 5 {print $0}' country-mutation/"$f"_all.mutation_sorted_lable.l > country-mutation/"$f"_all.mutation_sorted_lable.e5

awk 'BEGIN {FS=OFS="\t"}{gsub("-","_",$1)}1' country-mutation/"$f"_all.mutation_sorted_lable.e5 > country-mutation/"$f"_all.mutation_sorted_lable.e100 && mv country-mutation/"$f"_all.mutation_sorted_lable.e100 country-mutation/"$f"_all.mutation_sorted_lable.e5

sed -i 's/,/\t/g;s/-/\t/g' country-mutation/"$f"_all.mutation_sorted_lable.e5
awk -F "\t" '{print $1 "\t" $3 "\t" $7 "\t" $8"\t" $2 "-" $11 "\t" $12 "\t" $13 "\n"  $1 "\t" $3 "\t" $7 "\t" $9 "\t" $2 "-" $14 "\t" $15 "\t" $16 "\n"  $1 "\t" $3 "\t" $7 "\t" $10 "\t" $2 "-" $17 "\t" $18 "\t" $19 }' country-mutation/"$f"_all.mutation_sorted_lable.e5 > country-mutation/"$f"_all.mutation_sorted_lable.e52
paste -d  "\n" country-mutation/"$f"_all.mutation_sorted_lable.e12  country-mutation/"$f"_all.mutation_sorted_lable.e32  country-mutation/"$f"_all.mutation_sorted_lable.e52 > country-mutation/"$f"_all.mutation_sorted_lable
/bin/rm country-mutation/*.e*
sort -t$'\t' -n -k 4 country-mutation/"$f"_all.mutation_sorted_lable > country-mutation/"$f"_all.mutation_sorted_lable.e
mv country-mutation/"$f"_all.mutation_sorted_lable.e country-mutation/"$f"_all.mutation_sorted_lable
sed -i '/^$/d' country-mutation/"$f"_all.mutation_sorted_lable
mv country-mutation/"$f"_all.mutation_sorted_lable country-mutation/"$f"/"$f"_all.mutation_sorted_lable
done
