Pomfret: methylation-assisted phase block joiner. 

Pomfret exploits information from base modification in variant-homogeneous regions
to try to determine the read phasing and close phase gaps. 
Pomfret require no post-processing steps from a regular 
variant phasing pipeline. The input is a sorted and index 
BAM containing read alignments (MD tag required) and a phased VCF. 
The output is 
a methylation-phased VCF where variants are kept intact while their 
phase block IDs may be updated and
a GTF file containing the phase blocks.

# Getting started

Pomfret builds on linux mackOS. 
Building requires gcc, zlib and [HTSlib](https://github.com/samtools/htslib). 
If HTSlib is at a custom location, 
please have it in the environment variable `LD_LIBRARY_PATH` and 
set `DIR_HTSLIB` to its location.

```
# build pomfret
git clone https://git.oxfordnanolabs.local/xfeng/pomfret/
cd pomfret && make
# with custom HTSlib installation, do instead: cd pomfret && DIR_HTSLIB=/path/to/htslib make

# Run on test data
pomfret methphase -o out -c 60 --write-bam --vcf example/variants.vcf.gz example/phased.bam 2>log

# Run pomfret
# A basic run:
pomfret methphase -t4 -o prefix --vcf phased.vcf tagged.bam 2>log  # produces prefix.mp.gtf and prefix.mp.vcf 

# Use GTF to specify phase blocks instead:
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
# Bam must be sorted and indexed. VCF, GTF and tsv files assumed to be sorted, they can be plain or gzipped.
```

# Evaluations

Basecalled and aligned HG002 readsets was obtained from [Jan 2025 public data release](https://epi2me.nanoporetech.com/giab-2025.01/).
There were two runs, and each was basecalled by either the SUP or the HAC model.
For each basecalling setup, 
here me also merged the BAM files with `samtools merge -@32`
and subsampled them to different depths to be used in evaluations
by `samtools view -bh -s INT.FLOAT -@24` where the provided parameter is
calculated as the target depth divded by 90.
The depths were verified with `samtools depth` afterwards and were within expectation.

We performed variant calling with [wf-human-variation workflow[(https://github.com/epi2me-labs/wf-human-variation):
`nextflow run epi2me-labs/wf-human-variation -r v2.4.1 --out_dir wfhuman.out --bam input.bam --override_basecaller_cfg dna_r10.4.1_e8.2_400bps_MODEL@v5.0.0 --sex XY --ref hg38.fna --bed chr_all_all.bed --sample_name hg002 --threads 32 --snp -profile singularity`
where the MODEL was either "sup"  or "hac", according to the input.
We ran whatshap~\cite{martin2016whatshap} (version 2.4.dev2+gbe88057) with
``whatshap phase -o phased.vcf
--ignore-read-groups
--reference=hg38.fna hg002.wf_snp.vcf.gz input.bam''.
GTF file containing phase blocks was generated with
``whatshap stats -{}-gtf=phased.vcf.gtf phased.vcf''.
VCF is compressed and indexed.

For variant-based phasing, we used whatshap (version 2.4.dev2+gbe88057) or hapcut2 (v1.3.4, commit 0eb707). 
The commands for whatshap were
`whatshap phase --only-snvs -o phased.vcf
--ignore-read-groups
--reference=hg38.fna hg002.wf_snp.vcf.gz input.bam`,
which should match wf-human-variation workflow's
`--phase` option (see [code here](https://github.com/epi2me-labs/wf-human-variation/blob/036d1ac08ce5abb751463e385eb9b611ac43096d/modules/local/wf-human-snp.nf#L507))
.
The commands for hapcut2 were
`extractHAIRS --bam input.bam --vcf hg002.wf\_snp.vcf --ref hg38.fa --ont 1 --out out_fragment''
then 
`HAPCUT2 --fragments out_fragment --VCF hg002.wf_snp.vcf --output out_final`
.
Reference genome was the GCA\_000001405.15 no-alt analysis set.

For evaluation, we used `whatshap compare --names truth,phased  --ploidy 2 --ignore-sample-name --only-snvs --tsv-pairwise compare.tsv ref.vcf in.vcf` 
and `whatshap stats --tsv stats.tsv in.vcf`.
The reference VCF was GIAB's HG002 T2TQ100-V1.1 VCF.
