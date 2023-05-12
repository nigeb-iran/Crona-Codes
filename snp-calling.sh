#!/bin/bash
#please correct country name eg. Usa.fasta 
##clean directory
bwa index hcov-refrence-genome/hcov-refrence.fa
wc -l country > number-country
b=$(awk '{print $1}' number-country)
	for ((i=1;i<=$b;i++))                                                                                 
         do
           id=$(sed -n "$i"p country)
           mkdir $id
           sed -i 's/--/NN/g;s/C-/CN/g;s/T-/TN/g;s/A-/AN/g;s/G-/GN/g;s/-C/NC/g;s/-T/NT/g;s/-A/NA/g;s/-G/NG/g' $id.fas
                      split -l 50000 $id.fas $id/$id.fasta
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
          java -jar ~/snpEff/SnpSift.jar extractFields $id/"$id".fasta.vcf.ann CHROM POS REF ALT AD "EFF[0].CODON" "ANN[0].AA" "ANN[0].EFFECT" "EFF[1].CODON" "ANN[1].AA" "ANN[1].EFFECT" "EFF[2].CODON" "ANN[2].AA" "ANN[2].EFFECT" "EFF[4].CODON" "ANN[4].AA" "ANN[4].EFFECT"|sed 's/downstream_gene_variant/\t/g;s/upstream_gene_variant/\t/g;s/frameshift_variant&missense_variant/\t/g;s/frameshift_variant&start_lost//g;s/c..[1-9]*.delTins...//g;s/c..[1-9]*delCins...//g;s/c..[1-9]*delGins...//g;s/c..[1-9]*delAins...//g;s/c.-[1-9]*.>.//g;s/c.\*[1-9]*.delTins<\*>//g;s/c.\*[1-9]*.delCins<\*>//g;s/c.\*[1-9]*.delGins<\*>//g;s/c.\*[1-9]*.delAins<\*>//g;s/c.\*[1-9]*....//g;s/,<\*>//g;s/c.[0-9]*del.ins//g;s/p.*fs//g;s/\t*\t/\t/g;s/c.[0-9]*//g;s/,0\t/\t/g;s/<\*>/\t/g;s/\t\t\t/\t\t/g;s/\tp./\t/g;s/frameshift_variant&stop_lost&missense_variant//g;s/frameshift_variant&stop_lost&spli_region_variant//g;s/frameshift_variant&stop_lost//g;s/frameshift_variant&spli_region_variant&synonymous_variant//g;s/start_lost//g;s/stop_lost&spli_region_variant//g' > $id/$id.fasta.vcf.ann.final
         
          sed -i "1d" $id/$id.fasta.vcf.ann.final
	      done
