## Quick start ##
```
git clone https://github.com/zhangrengang/SLRhunter
cd SLRhunter

# install
pip3 install .

chmod +x bin/* /src/*
export PATH=$PATH:`pwd`/bin/:`pwd`/src/

```
## Dependencies ##
  - kmc
  - bwa
  - pandepth
  - gemma

## kmer-based GWAS ##

Input files:

`gwas.design`:
```
m1      male    2_CleanData/CleanData/m1_filter_R1.fastq.gz     2_CleanData/CleanData/m1_filter_R2.fastq.gz
m2      male    2_CleanData/CleanData/m2_filter_R1.fastq.gz     2_CleanData/CleanData/m2_filter_R2.fastq.gz
m3      male    2_CleanData/CleanData/m3_filter_R1.fastq.gz     2_CleanData/CleanData/m3_filter_R2.fastq.gz
...
f1      female  2_CleanData/CleanData/f1_filter_R1.fastq.gz     2_CleanData/CleanData/f1_filter_R2.fastq.gz
f2      female  2_CleanData/CleanData/f2_filter_R1.fastq.gz     2_CleanData/CleanData/f2_filter_R2.fastq.gz
f3      female  2_CleanData/CleanData/f3_filter_R1.fastq.gz     2_CleanData/CleanData/f3_filter_R2.fastq.gz
...
```
`phenotype.txt`:
```
cat gwas.design| cut -f2 | sed 's/female/0/;s/male/1/' > phenotype.txt
```

Reference genome: `ref/ref.fa`.


Count kmers from reads and run GWAS
```
sdrhunter.py -s gwas.design -lower_count 5 -p 60 -g ref/ref.fa

bedcov=sdr-results/k31.kmer.mat
pheno=phenotype.txt
N=`head $pheno -n1 | awk '{print NF}'`
ns=`seq 1 $N`
for n in $ns
do
	gemma -g $bedcov.geno$g -a $bedcov.snps -p $pheno -lm 1 -n $n -o $pre 
done

```
There are also some other models for gemma.


## coverage-based GWAS ##
Get depth coverage from BAM files:
```
tmpdir=tmp
mkdir $tmpdir -p
bin=200
samtools faidx ref/ref.fa
python ~/src/fai2bed.py ref/ref.fa.fai $bin > ref/ref.fa.bed
bed=ref/ref.fa.bed

for q in 1 0
do

	for sm in `cut -f 1 $design`
	do
        BAM=3_MappingResult/$sm.bam
        pandepth  -i $BAM -t 4 -o $tmpdir/$sm.$q -q $q -b $bed && \
        zcat $tmpdir/$sm.$q.bed.stat.gz | sed '1d' | cut -f 7 > $tmpdir/$sm.$q.cov
    done

	covs=`cat $design | awk -v dd=$tmpdir -v q=$q '{print dd"/"$1"."q".cov"}' | tr "\n" " "`
	paste ref/ref.fa.bed $covs > all_samples.q$q.bedcov
done

```
The bin size can be increased to reduce memory costs, or decreased to refine resolution.

Then run GWAS:

```
pheno=phenotype.txt
N=`head $pheno -n1 | awk '{print NF}'`
ns=`seq 1 $N`

for q in 0 1
do
	bedcov=all_samples.q$q.bedcov
	Gemma.py cov2gemma $bedcov fold=1 > $bedcov.geno1 &
	Gemma.py cov2gemma $bedcov fold=2 > $bedcov.geno2 &
	Gemma.py cov2gemma0 $bedcov > $bedcov.geno0 &
done
wait

for q in 0 1
do
	bedcov=all_samples.q$q.bedcov

	# linear model
	for g in 0 1 2
	do
		for n in $ns
		do
			pre=cov.q$q.geno$g.lm.n$n
			gemma -g $bedcov.geno$g -p $pheno -a $bedcov.snps -lm 1 -n $n -o $pre 

		done
	done

done

```
