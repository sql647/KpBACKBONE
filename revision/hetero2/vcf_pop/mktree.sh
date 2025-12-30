VCF2Matrix.pl ../out_vcf/pop1_10000.vcf 0 > pop1_10000.matrix
snp_tongyong.pl pop1_10000.matrix > pop1_10000.fa
FastTree -nt pop1_10000.fa > pop1_10000.FastTree
