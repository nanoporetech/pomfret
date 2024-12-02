Pomfret: methylation-assisted phase block joiner. 

Pomfret exploits information from base modification in variant-homogeneous regions
to try to determine the read phasing and close phase gaps. 
Pomfret require no post-processing steps from a regular 
variant phasing pipeline. The minimum inputs are a sorted and index 
BAM containing read alignments and a phased VCF. The output is 
a GTF file containing meth-phased phase blocks 
and a meth-phased VCF, where variants are kept intact and their 
phase block IDs may be updated. 
See below for other input combinations. 

# Getting started

```
# Install pomfret (requires gcc, zlib and htslib)
git clone https://git.oxfordnanolabs.local/xfeng/pomfret/
cd pomfret && make
# if htslib is at a custom location, make sure it is 
# in the LD_LIBRARY_PATH and compile with 
# DIR_HTSLIB=/path/bin/htslib make

# Run on test data
pomfret methphase -o out -c 60 --vcf example/variants.vcf.gz example/phased.bam 2>log

# Run pomfret
# Join phase blocks:
pomfret methphase -t4 -o output --vcf phased.vcf reads_haplotagged.bam 2>log # produces a vcf updated methylation-joined phase block IDs, and a bam where read haplotags are updated.
pomfret methphase -t4 --vcf phased.vcf --bam-is-untagged --dont-write-bam --gtf blocks.gtf reads.bam 2>log # haplotag the input bam based on the vcf, then try to join phase blocks specified by the gtf (overrides vcf phase blocks). Outputs an updated vcf.
# Other utilities:
pomfret varhaptag -o output phased.vcf reads.bam 2>log # produces a haplotagged bam and a tsv mapping reads to the haplotags
pomfret methstat -o out.tsv --vcf phased.vcf reads.bam 2>log  # lists heterozygous-looking methylation sites
# -h or --help prints options.
# Input files assumed to be persistent and static throughout the run. 
# Bam must be sorted and indexed. Vcf/gtf/tsv files assumed to be sorted, they can be plain or gz'd.
```

# TBD

- methphase: detect read coverage and make -c optional.
- methphase & varhaptag: paralleled input haplotagging.
- methphase & varhaptag: let paralleled methods also have single thread version.
