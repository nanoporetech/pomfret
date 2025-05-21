# Pomfret

**Pomfret** is a methylation-assisted phase block joiner that leverages **5-methylcytosine (5mC)** signals to improve read phasing and resolve phase block discontinuities.

Haplotype-specific 5mC methylation is common in **humans and other vertebrates**, and can reveal long-range haplotype structure even in variant-sparse regions.
Pomfret uses this information to bridge phase blocks and produce more complete phasing results.

---

## Overview

Pomfret integrates directly into a standard variant phasing workflow.
No additional post-processing is required.

**Inputs:**
- A **sorted and indexed BAM** file with alignments (MD tag required)
- A **phased VCF** file

**Outputs:**
- A **VCF** file with updated phase block IDs (variant calls are unchanged)
- A **GTF** file defining the methylation-informed phase blocks

---

## Getting started

### Build Instructions

Pomfret builds on **Linux** and **macOS**. Dependencies:

- ``gcc``
- ``zlib``
- [HTSlib](https://github.com/samtools/htslib)

```bash
# Clone and build Pomfret
git clone https://git.oxfordnanolabs.local/xfeng/pomfret/
cd pomfret
make
```

## Quick Start

Run Pomfret on the provided example data:

```bash
pomfret \
    methphase \
	-o out \
	-c 60 \
	--write-bam \
	--vcf example/variants.vcf.gz \
	example/phased.bam
```

Standard usage:
```bash
pomfret \
    methphase \
	-o prefix \
	--vcf phased.vcf \
	tagged.bam
```

---

## Additional Options

* To define phase blocks with a GTF file instead of a VCF, use the ``--gtf`` option.
* For untagged BAMs (e.g. from Dorado), use ``--bam-is-untagged``.
  Pomfret will infer read haplotypes from the basecalled sequence and the VCF.

---

## Reporting Phase Join Quality

Pomfret also provides a ``report`` command to assess methylation phasing quality:

```
# Report success, error and no-join of meth-phasing in a subset of the phased regions:
pomfret \
    report \
	-o prefix \
	--vcf phased.vcf \
	tagged.bam
```

The report can be focued to specific regions using ``--chunk-size 50000 --chunk-stride 100000``.

---

## Notes

* Add the ``-t`` option to specify the number of threads to use.
* All input files must remain unchanged during execution.
* BAM files must be sorted and indexed.
* VCF, GTF, and TSV files must be sorted, and may be plain or gzipped.
* Use ``-h`` or ``--help`` for command options.

For custom HTSlib install: ``LD_LIBRARY_PATH=/path/to/htsliblib DIR_HTSLIB=/path/to/htslib make``

---

# Evaluations

Basecalled and aligned HG002 readsets was obtained from [Jan 2025 public data release](https://epi2me.nanoporetech.com/giab-2025.01/).
There were two runs, and each was basecalled by either the SUP or the HAC model.
For each basecaller, BAM files were merged with `samtools merge -@32`.
BAMs were then sampled to different depths to be used in evaluations by `samtools view -bh -s INT.FLOAT -@24`.

We performed variant calling with [wf-human-variation workflow](https://github.com/epi2me-labs/wf-human-variation):
`nextflow run epi2me-labs/wf-human-variation -r v2.4.1 --out_dir wfhuman.out --bam input.bam --override_basecaller_cfg dna_r10.4.1_e8.2_400bps_MODEL@v5.0.0 --sex XY --ref hg38.fna --bed chr_all_all.bed --sample_name hg002 --threads 32 --snp -profile singularity`
where the MODEL was either ``sup``  or ``hac``, according to the input.

whatshap~\cite{martin2016whatshap} (version 2.4.dev2+gbe88057) was run as follows:

```
whatshap \
    phase \
	-o phased.vcf \
    --ignore-read-groups \
    --reference=hg38.fna \
	hg002.wf_snp.vcf.gz \
	input.bam
```

GTF file containing phase blocks was generated with ``whatshap stats -{}-gtf=phased.vcf.gtf phased.vcf``.
VCF was compressed and indexed.

For variant-based phasing, we used whatshap (version 2.4.dev2+gbe88057) or hapcut2 (v1.3.4, commit 0eb707). 
Whatshap was run as follows:

```
whatshap \
    phase \
	--only-snvs \
	-o phased.vcf \
	--ignore-read-groups \
	--reference=hg38.fna \
	hg002.wf_snp.vcf.gz \
	input.bam
```
This matches [wf-human-variation workflow's ``--phase`` option](https://github.com/epi2me-labs/wf-human-variation/blob/036d1ac08ce5abb751463e385eb9b611ac43096d/modules/local/wf-human-snp.nf#L507).

Hapcut was run as follows:
```
extractHAIRS \
	--bam input.bam \
	--vcf hg002.wf_snp.vcf \
	--ref hg38.fa \
	--ont 1 \
	--out out_fragment
HAPCUT2 \
	--fragments out_fragment \
	--VCF hg002.wf_snp.vcf \
	--output out_final
```

Reference genome was the GCA_000001405.15 no-alt analysis set.
Pomfret was ran as follows:
```
pomfret \
	methphase \
	-o output \
	--vcf phased.vcf.gz \
	-u \
	-t 32 \
	input.bam
```

This command assumes the input BAM is not haplotagged and performs tagging as a preprocessing step.
``pomfret methphase -t32`` runs took 20-30 minutes with peak memory around 2.5GiB. 

For evaluation, ``whatshap`` was run as follows:
```
whatshap \
	compare \
	--names truth,phased \
	--ploidy 2 \
	--ignore-sample-name \
	--only-snvs \
	--tsv-pairwise compare.tsv \
	ref.vcf \
	in.vcf
whatshap \
	stats \
	--tsv stats.tsv \
	in.vcf
```

The reference VCF was GIAB's HG002 T2TQ100-V1.1 VCF.

|dataset    | phasing  | phase block N50| switch error rate|
|:----      |:------   | :------             | :------             |
|combined-hac-30x | whatshap | 0.69Mb | 0.128%|
|combined-hac-30x | pomfret-whatshap | 1.52Mb | 0.127%|
|combined-hac-60x | whatshap | 1.00Mb | 0.121%|
|combined-hac-60x | pomfret-whatshap | 2.74Mb | 0.121%|
|combined-sup-30x | whatshap | 0.77Mb | 0.103%|
|combined-sup-30x | pomfret-whatshap | 1.45Mb | 0.104%|
|combined-sup-60x | whatshap | 1.00Mb | 0.115%|
|combined-sup-60x | pomfret-whatshap | 2.19Mb | 0.116%|
