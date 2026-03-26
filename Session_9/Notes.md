# Session 9 — Variant Calling Fundamentals

This session introduces the core logic behind variant calling: how we go from aligned reads (SAM/BAM) to identifying SNPs and small indels. The goal is conceptual understanding, not running full pipelines or using heavy tools.

## What Variant Calling Actually Is
Variant calling attempts to answer:
“Given all the reads covering this genomic position, what is the most likely true base (or genotype)?”

Because sequencing data contains:
- errors
- mismatches
- mapping uncertainty
- PCR artefacts
- low‑coverage regions

…variant calling is probabilistic, not deterministic.

## Pileups: The Foundation of Variant Calling
A pileup is simply all reads stacked at a given reference position.
Example:
Reference: A
Reads:     A A A A G A A

Here, most reads support A, but one read supports G.
Variant callers examine:
- base counts
- base qualities
- mapping qualities
- strand orientation
- depth
- local alignment context

…and decide whether the evidence supports a variant.

## Types of Variants
SNPs
Single‑base substitutions
Example: A → G
Indels
Insertions or deletions relative to the reference
Example:
- Insertion: A → AT
- Deletion: AT → A

Complex variants
Small combinations of SNPs + indels

## Why Variant Calling Is Hard
Variant callers must distinguish true variants from:
- sequencing errors
- low‑quality bases
- misaligned reads
- repeats
- strand bias
- low coverage
- PCR duplicates

This is why modern callers use Bayesian models, likelihoods, and heuristics.

## VCF: The Output Format
Variant callers produce a VCF (Variant Call Format) file.
Key fields:
CHROM - Chromosome
POS - Position
REF - Reference Base
ALT - Alternate Base(s)
QUAL - Confidence Score
INFO - Extra annotations (Depth, Allele frequency, filters)

Example:
chr1    105    .    A    G    47    PASS    DP=15;AF=0.2

## Deliverable for This Session
A small Python script that:
- Reads a SAM file
- Builds a simple pileup (counts A/C/G/T at each position)
- Prints positions where a non‑reference base exceeds a threshold

This is a toy variant caller, but it teaches the core logic behind real tools like:
- bcftools
- FreeBayes
- GATK HaplotypeCaller

Please find in SAM_Parser.py in the Session_9 folder

## What I Learned
- How pileups represent evidence for or against a variant
- Why variant calling is fundamentally probabilistic
- How sequencing and mapping errors influence calls
- The structure and purpose of VCF files
- How to implement a minimal variant‑calling heuristic in Python
