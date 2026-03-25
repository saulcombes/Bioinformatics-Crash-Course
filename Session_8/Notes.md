# Session 8 — Alignment & SAM File Fundamentals

This session introduces the core ideas behind sequence alignment and the structure of SAM files. The goal is to build intuition for how aligners work and how to interpret their output, without running real tools yet.

## What Alignment Is

Alignment is the process of taking sequencing reads and determining where they came from in the reference genome.
An aligner tries to answer:
“Given this read, where in the genome is it most likely to originate?”

Because genomes contain repeats, errors, and variation, this is not trivial. Modern aligners use fast heuristics to approximate the best match.

## How Modern Aligners Work (High‑Level)

Most aligners follow the same conceptual steps:
1. Indexing the Reference
The genome is pre‑processed into a searchable data structure (FM‑index, minimizers, k‑mers).
This allows fast lookup of short sequences.
2. Seeding
Short exact matches (“seeds”) between the read and the reference are found.
These seeds act as anchors.
3. Extension
The aligner extends outward from each seed, allowing mismatches, insertions, and deletions.
4. Scoring
Each possible alignment is scored based on:
- matches
- mismatches
- gaps
- penalties
The best‑scoring alignment is chosen.
5. Output
The result is written to a SAM (or BAM) file.

## What a SAM File Contains

A SAM file is a tab‑delimited text file with two parts:
1. Header (starts with @)
Describes:
- reference genome
- chromosome lengths
- software used
- read groups
Example:
@SQ SN:chr1 LN:248956422
@PG ID:bwa PN:bwa VN:0.7.17


2. Alignment Section (one line per read)
Each read has 11 mandatory fields:
|  |  | 
|  |  | 
|  |  | 
|  |  | 
|  |  | 
|  |  | 
|  |  | 
|  |  | 
|  |  | 
|  |  | 
|  |  | 
|  |  | 



## Key Concepts

CIGAR Strings
Describe how the read aligns:
- M = match/mismatch
- I = insertion
- D = deletion
- S = soft‑clipped
- H = hard‑clipped
Example:
98M2S → 98 aligned bases, 2 soft‑clipped.
  
MAPQ (Mapping Quality)
A measure of confidence (0–60):
- 60 = unique, confident alignment
- 0 = ambiguous or unmapped
  
FLAG
A bitwise code describing read status:
- 0 = mapped, forward
- 16 = mapped, reverse
- 4 = unmapped

## Deliverable for This Session

A small Python script that:
- reads a SAM file
- skips headers
- extracts key fields (QNAME, RNAME, POS, MAPQ, CIGAR)
- prints them for each read
This builds practical SAM literacy and prepares for later sessions involving BAM files, variant calling, and QC.

## What I Learned

- How aligners find approximate matches efficiently
- Why alignment is probabilistic, not exact
- How to interpret the essential fields in a SAM file
- How CIGAR strings encode mismatches, indels, and clipping
- How MAPQ reflects alignment confidence
- How to parse SAM files programmaticall
