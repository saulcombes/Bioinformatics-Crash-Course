# Session 12 — Variant Calling, Filtering & VCF Output

## Overview

Session 12 brings together everything from the alignment and CIGAR‑parsing stages to produce real variant calls.  
By the end of this session, the pipeline can:

- build a pileup from aligned reads  
- apply variant‑calling rules  
- detect SNPs and insertions  
- apply basic filtering  
- output a valid VCF file  

This is the first session where the pipeline behaves like a miniature real‑world variant caller.

---

## 12.1 — Pileup Formation

Reads are walked base‑by‑base using their SAM POS and CIGAR operations.  
For each `(chrom, pos)`:

- bases from all reads are collected  
- depth is computed implicitly  
- insertions and deletions are recorded separately  

Reference sequences are loaded from FASTA and stored as **strings**, indexed with `pos - 1` because:

- FASTA strings → 0‑based  
- SAM POS → 1‑based  

This correction is essential for accurate REF base lookup.

---

## 12.2 — Variant Calling Logic

Variants are called using simple but biologically meaningful rules.

### SNPs

A SNP is called when:

- depth ≥ 2 (or 3)  
- ALT count ≥ 1  
- ALT fraction ≥ 0.2  
- ALT count > REF count  

This automatically filters out:

- ties  
- noise  
- low‑support mismatches  

### Insertions
Insertions are represented as:

REF = reference_base ALT = reference_base + inserted_sequence

This is already **left‑normalized**, matching standard VCF conventions.

---

## 12.3 — Filtering (Implicit + Explicit)

Even before writing the VCF, the variant list already applies several filters:

- minimum depth  
- minimum ALT evidence  
- minimum allele fraction  
- REF vs ALT comparison  
- strand consistency (implicit)  
- left‑normalisation for indels  

These rules mirror the behaviour of real variant callers like FreeBayes and GATK.

---

## 12.4 — VCF Output
A minimal but valid VCF is produced with:

##fileformat=VCFv4.2 #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO

Each variant written as:

chrom   pos   .   REF   ALT   .   PASS   .

This VCF is fully standards‑compliant and can be loaded into IGV, bcftools, or any downstream tool.

---

## Summary
By the end of Session 12, the pipeline can:

- interpret SAM alignments  
- walk CIGAR strings  
- build a pileup  
- call SNPs and insertions  
- apply biologically meaningful filters  
- output a valid VCF  

This completes the core of a functional variant‑calling engine.

