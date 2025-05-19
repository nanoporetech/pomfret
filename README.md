Pomfret: methylation-assisted phase block joiner. 

Pomfret exploits information from base modification in variant-homogeneous regions
to try to determine the read phasing and close phase gaps. 
Pomfret require no post-processing steps from a regular 
variant phasing pipeline. The input is a sorted and index 
BAM containing read alignments and a phased VCF. The output is 
a GTF file containing meth-phased phase blocks 
and a meth-phased VCF, where variants are kept intact while their 
phase block IDs may be updated. 

# Getting started

```
# Install pomfret (requires gcc, zlib and htslib)
git clone https://git.oxfordnanolabs.local/xfeng/pomfret/
cd pomfret && make
# if htslib is at a custom location, make sure it is 
# in the LD_LIBRARY_PATH and compile with 
# DIR_HTSLIB=/path/bin/htslib make

# Run on test data
pomfret methphase -o out -c 60 --write-bam --vcf example/variants.vcf.gz example/phased.bam 2>log

# Run pomfret
# Join phase blocks:
pomfret methphase -t4 -o prefix --vcf phased.vcf tagged.bam 2>log  # produces prefix.mp.gtf and prefix.mp.vcf 

# Use GTF to specify phase blocks instead
pomfret methphase -t4 -o prefix --gtf phased.gtf tagged.bam 2>log  # no vcf output

# Override phase blocks in VCF with GTF, also infer
# read haplotags of the input bam using the VCF as a preprocessing step:
pomfret methphase -t4 -o prefix --vcf phased.vcf --bam-is-untagged --gtf blocks.gtf untagged.bam 2>log

# Report success, error and no-join of meth-phasing in a subset of the phased regions: 
pomfret report -t4 -o prefix --vcf phased.vcf tagged.bam 2>log
# or:
pomfret report -t4 -o prefix --vcf phased.vcf --chunk-size 50000 --chunk-stride 100000 --bam-is-untagged untagged.bam 2>log

# -h or --help prints options.
# Input files assumed to be persistent and static throughout the run. 
# Bam must be sorted and indexed. Vcf/gtf/tsv files assumed to be sorted, they can be plain or gz'd.
```

# Evaluations
