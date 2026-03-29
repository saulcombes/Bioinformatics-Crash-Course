# Session 11 — CIGAR Parsing & Alignment Walking

## Overview
Session 11 focuses on interpreting SAM alignments by converting each read’s CIGAR string into a base-by-base description of how the read aligns to the reference genome. This alignment-walking engine is the foundation for all downstream variant calling.

## Key Concepts

### 1. CIGAR Operations
CIGAR strings describe how a read aligns to the reference. Each operation may consume read bases, reference bases, or both.

M  Match or mismatch (consumes read and reference)
I  Insertion relative to reference (consumes read only)
D  Deletion relative to reference (consumes reference only)
S  Soft-clip (consumes read only, not aligned)
*  Unmapped read

Soft-clipped bases (S) are present in the read sequence but do not align to the reference.

### 2. Read Length Rule
A valid SAM alignment must satisfy:

Sum of read-consuming CIGAR operations = length of SEQ

Read-consuming operations are M, I, S, =, X.

### 3. Reverse-Strand Reads
Reads with FLAG 16 are aligned to the reverse strand.

Before walking the CIGAR:
- reverse-complement the read sequence
- apply the CIGAR left-to-right as written

## CIGAR Parsing
CIGAR strings are parsed into (length, operation) pairs.

Example:
"10M1I10M" becomes [('10','M'), ('1','I'), ('10','M')]

Unmapped reads (CIGAR = * or FLAG 4) are skipped.

## Alignment Walking
The walker processes each CIGAR operation and emits alignment events:

(ref_pos, base) for aligned bases
("INS", base) for insertions
(ref_pos, "DEL") for deletions
soft-clips produce no events

Reference positions are 1-based to match SAM conventions.

Example for 10M1I10M:
[(100,'A'), (101,'C'), ..., ('INS','T'), (110,'G'), ...]

## Event Lists
Each read produces a list describing exactly how it aligns to the reference.

Examples:

read001 -> 50M -> events at positions 100–149

read002 -> 40M2S -> events at positions 100–139

read003 -> unmapped -> []

read004 -> 10M1I10M -> events plus an insertion

read005 -> 5S20M -> events at positions 100–119

These event lists feed directly into Session 12.

## Achievements in Session 11
- Correct CIGAR parsing
- Correct handling of M, I, D, S
- Correct handling of unmapped reads
- Correct handling of reverse-strand reads
- Accurate reference coordinate tracking
- Correct generation of alignment event lists
- Validation against a test SAM file

This completes the alignment-walking engine.
