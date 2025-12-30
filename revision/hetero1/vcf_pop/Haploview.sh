VCF2Matrix.pl pop1maf002.recode.vcf 0 > pop1maf002.matrix
perl /home/public/PublicSoftware/TRSNPMatrix2MajorAMinorT.pl pop1maf002.matrix > pop1maf002.matrix.MajorAMinorT
perl /home/public/PublicSoftware/snp2info.pl pop1maf002.matrix.MajorAMinorT pop1maf002.matrix.MajorAMinorT.info
perl /home/public/PublicSoftware/snp2ped.pl pop1maf002.matrix.MajorAMinorT pop1maf002.matrix.MajorAMinorT.ped
java -jar /home/public/PublicSoftware/Haploview.jar -memory 102400 -n -info pop1maf002.matrix.MajorAMinorT.info -pedfile pop1maf002.matrix.MajorAMinorT.ped -maxdistance 30 -dprime -minGeno 0.6 -minMAF 0 -hwcutoff 0  -log pop1maf002.matrix.MajorAMinorT.log ; gzip pop1maf002.matrix.MajorAMinorT.ped.LD
LD_decay.StatRR_avg_v2.pl pop1maf002.matrix.MajorAMinorT.ped.LD.gz pop1maf002.matrix.MajorAMinorT.statRR.avg.matrix
LD_decay.Stat_StepAvg.pl pop1maf002.matrix.MajorAMinorT.statRR.avg.matrix 10 > pop1maf002.matrix.MajorAMinorT.statRR.avg.matrix.step10
StatAll.pl pop1maf002.matrix.MajorAMinorT.statRR.avg.matrix.step10 > pop1maf002.matrix.MajorAMinorT.statRR.avg.matrix.step10.All.stat


