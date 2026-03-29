# Session 10 — Variant Calling & VCF Output
This session extends the pileup engine from Session 9 into a minimal but functional variant caller. The goal is to take base‑level pileup information and convert it into a list of high‑confidence variants, then export them in VCF (Variant Call Format) — the standard format used across genomics.

# Overview

Variant calling in this session follows a simple, transparent logic:
- Identify the most common base at each genomic position.
- Compare it to the reference base.
- Compute:
- ALT count
- Depth (DP)
- Allele fraction (AF)
- Apply filtering rules to decide whether a position is a variant.
- Export all passing variants to a valid VCF v4.2 file.
  
This produces a lightweight but realistic variant‑calling workflow.

## Variant Calling Rules

A position is considered a variant if:
- most_common_base != reference_base
- depth ≥ 2
- alt_count ≥ 1
- alt_fraction ≥ 0.2
  
These thresholds mimic the simplest form of variant filtering used in real pipelines.

## Variant Table

Before writing to VCF, variants are summarised in a human‑readable table:

CHROM - Chromosome name

POS - 1-based genomic position

REF - Reference allele

ALT - Alternate allele

DEPTH - Total read depth

REF_COUNT - Number of reads supporting REF

ALT_COUNT - Number of reads supporting ALT

ALT_FRACTION - ALT_COUNT/DEPTH

This table is printed with fixed‑width formatting for readability.

## VCF Output

VCF (Variant Call Format) is the standard for storing variant calls.

A minimal valid VCF contains:

Header


##fileformat=VCFv4.2

#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO



Data rows

Each variant is written as:

CHROM  POS  .  REF  ALT  .  PASS  DP=<depth>;AF=<allele_fraction>

- . is used for fields not implemented (ID, QUAL).
- PASS indicates the variant met all filtering criteria.
- INFO stores depth (DP) and allele fraction (AF).

Spacing/alignment does not affect machine parsing — VCF is tab‑delimited.

## File Naming

VCF files typically use:
- .vcf — plain text
- .vcf.gz — bgzip‑compressed
- .bcf — binary VCF (bcftools)
  
This project outputs standard .vcf files.

## What This Session Achieves

By the end of Session 10, the pipeline can:
- Parse SAM alignments
- Build a pileup
- Identify variants
- Filter low‑confidence sites
- Produce a clean variant table
- Export a valid VCF file
  
This completes the core logic of a simple variant caller.
