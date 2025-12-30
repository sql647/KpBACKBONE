VCF2Matrix.pl ../out_vcf/pop1_5000.vcf 0 > pop1_5000.matrix
snp_tongyong.pl pop1_5000.matrix > pop1_5000.fa
FastTree -nt pop1_5000.fa > pop1_5000.FastTree
