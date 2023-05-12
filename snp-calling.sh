#!/bin/bash
#please correct country name eg. Usa.fasta 
##clean directory

#/bin/rm -f ./*/*.final


#/bin/rm -f ./*/*.html


#/bin/rm -f ./*/*.bcf


#/bin/rm -f ./*/*.vcf


#/bin/rm -f ./*/*.bai


#/bin/rm -f ./*/*.ann


#/bin/rm -f ./*/*.mut


#/bin/rm -f ./*/*.I16


#/bin/rm -f ./*/*.header


#/bin/rm -f ./*/*.csv


#/bin/rm -f ./*/*.tsv


#/bin/rm -f ./*/*.fai


#/bin/rm -f ./*/*.bam


#/bin/rm -f ./*/*.txt


#/bin/rm -f ./*/*.aa

#/bin/rm -f ./*/*.sam

#/bin/rm -f ./*/*.eff
#/bin/rm -f ./*/*.fastaa*

#bwa index hcov-refrence-genome/hcov-refrence.fa
wc -l country > number-country
b=$(awk '{print $1}' number-country)
	for ((i=1;i<=$b;i++))                                                                                 
         do
           id=$(sed -n "$i"p country)
           mkdir $id
           #cp "$id".fas $id/$id.fasta
           sed -i 's/--/NN/g;s/C-/CN/g;s/T-/TN/g;s/A-/AN/g;s/G-/GN/g;s/-C/NC/g;s/-T/NT/g;s/-A/NA/g;s/-G/NG/g' $id.fas
           #seqkit fx2tab    $id.fas > $id.fas-single
           #mv $id.fas-single $id.fas
           #sed -i 's/^/>/g' $id.fas
          # sed -i 's/\t/\n/g' $id.fas
           #sed -i '/^$/d' $id.fas
           split -l 50000 $id.fas $id/$id.fasta
	   #/bin/rm $id/$id.fasta
           for g in $(ls $id/$id.fasta*);
           do 
           bwa mem -t 40 hcov-refrence-genome/hcov-refrence.fa $g > $g-1.sam
           samtools view -S -b $g-1.sam > $g-1.bam
           samtools sort $g-1.bam -o $g.sorted.bam
           samtools faidx $g
           samtools index $g.sorted.bam
           done
           /bin/rm $id/$id.fasta.sorted.bam
           bcftools mpileup -B --max-depth 500000 -F 0.00002 -f  hcov-refrence-genome/hcov-refrence.fa $id/$id.fasta*.sorted.bam -a FORMAT/DP4,INFO/AD > $id/$id.fasta.bcf
           
           bcftools view $id/$id.fasta.bcf > $id/$id.fasta.vcf
          java -jar ~/snpEff/snpEff.jar -c ~/snpEff/snpEff.config hcov-refrence $id/$id.fasta.vcf -stats $id/$id.fasta.html > $id/$id.fasta.vcf.ann
          #java -Xmx1G -jar snpEff.jar build -gff3 hcov-refrence
          ###java -jar ~/snpEff/SnpSift.jar extractFields $id/$id.fasta.vcf.ann  "ANN[*].AA"  > $id/$id.fasta.mut
          #sed -i  's/^p.*fs\t//g;s/^p.*fs//g;s/^p.//g;s/\tp./-/g;s/\t//g' $id/$id.fasta.mut
          #java -jar ~/snpEff/SnpSift.jar extractFields $id/$id.fasta.vcf.ann  "ANN[*].EFFECT" > $id/$id.fasta.eff
          java -jar ~/snpEff/SnpSift.jar extractFields $id/"$id".fasta.vcf.ann CHROM POS REF ALT AD "EFF[0].CODON" "ANN[0].AA" "ANN[0].EFFECT" "EFF[1].CODON" "ANN[1].AA" "ANN[1].EFFECT" "EFF[2].CODON" "ANN[2].AA" "ANN[2].EFFECT" "EFF[4].CODON" "ANN[4].AA" "ANN[4].EFFECT"|sed 's/downstream_gene_variant/\t/g;s/upstream_gene_variant/\t/g;s/frameshift_variant&missense_variant/\t/g;s/frameshift_variant&start_lost//g;s/c..[1-9]*.delTins...//g;s/c..[1-9]*delCins...//g;s/c..[1-9]*delGins...//g;s/c..[1-9]*delAins...//g;s/c.-[1-9]*.>.//g;s/c.\*[1-9]*.delTins<\*>//g;s/c.\*[1-9]*.delCins<\*>//g;s/c.\*[1-9]*.delGins<\*>//g;s/c.\*[1-9]*.delAins<\*>//g;s/c.\*[1-9]*....//g;s/,<\*>//g;s/c.[0-9]*del.ins//g;s/p.*fs//g;s/\t*\t/\t/g;s/c.[0-9]*//g;s/,0\t/\t/g;s/<\*>/\t/g;s/\t\t\t/\t\t/g;s/\tp./\t/g;s/frameshift_variant&stop_lost&missense_variant//g;s/frameshift_variant&stop_lost&spli_region_variant//g;s/frameshift_variant&stop_lost//g;s/frameshift_variant&spli_region_variant&synonymous_variant//g;s/start_lost//g;s/stop_lost&spli_region_variant//g' > $id/$id.fasta.vcf.ann.final
         
          sed -i "1d" $id/$id.fasta.vcf.ann.final
	  #sed -i 's/p.*.fs//g;s/\tframeshift_variant&missense_variant//g;s/^\t//g;s/p.//g;s/_variant//g' $id/$id.fasta.mut
          #awk '{print $2"\t"$3}' $id/$id.fasta.mut > $id/$id.fasta.mut.aa
          ###sed -i "1d" $id/$id.fasta.mut
          ###sed -i "1d" $id/$id.fasta.eff
          #bcftools query  -f '%I16\n' $id/$id.fasta.bcf | awk -F "," '{print $1}' > $id/$id.fasta.I16
          ###/usr/local/bcf/bin/bcftools query  -f '%CHROM	%POS	%ID	%REF	%ALT	%AD\n' $id/$id.fasta.bcf > $id/$id.fasta.header
	  ###sed -i 's/,<\*>//g;s/<\*>//g;s/,0$//g' $id/$id.fasta.header
          ###paste $id/$id.fasta.header $id/$id.fasta.mut  $id/$id.fasta.eff > $id/$id.fasta.csv
          ###sed -i "s/NC_045512.2/$id/g" $id/$id.fasta.csv
          #paste $id/$id.fasta.tsv $id/$id.fasta.mut.aa > $id/$id.fasta.csv
          #sed 's/p.Met1fs//g;s/p.Ter122delins??????//g' $id/$id.fasta.csv > $id/$id.fasta.csv3
          #awk 'BEGIN{FS=OFS="\t"}{if ($9 != "") { $8=$8"-"$9;$9=""}  print $0}' $id/$id.fasta.csv3 > $id/$id.fasta.csv2
          #awk -F '\t' '{ if( $8 != "") print $0}' $id/$id.fasta.csv2 > $id/$id.fasta.csv-final
          #sed -i 's/,<\*>//g;s/<\*>//g' $id/$id.fasta.csv-final
          #sed -i 's/\tp./\t/g' $id/$id.fasta.csv-final
          #sed -i 's/,<\*>//g;s/<\*>//g' $id/$id.fasta.csv2     
          #sed -i 's/\tp./\t/g' $id/$id.fasta.csv2                   
          #mv -f $id/$id.fasta.csv2  $id/$id.fasta.csv
          #cp $id/$id.fasta.csv ../web/$f/
          #cp $id/$id.fasta.html ../web/$f/
          #seqkit sort -n $id/$id.fasta > $id/$id.fasta.sort
          #sed -i '/^>/! {s/\-/N/g}' $id/$id.fasta.sort
          #cdhit -i $id/$id.fasta.sort -o $id/$id.fasta.cdhit -c 1 -G 1 -s 0.999 -M 0
          #cat $id/$id.fasta.cdhit >> genes/$f/all-cdhit
          #/bin/cp genes/$f/all-cdhit-metadata ../nextstrain/$f
          #/bin/cp genes/$f/all-cdhit ../nextstrain/$f
          #rm -i -f $id/$id.fasta.csv3
          #rm -i -f $id/$id.fasta.csv2
        done

