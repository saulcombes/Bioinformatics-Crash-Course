# Session 13 — Base Quality Scores, Evidence Weighting, and Low‑Quality Filtering

## Overview
This session introduces base quality scores (Phred QUAL) into the variant‑calling pipeline.
By incorporating QUAL, the pileup becomes more reliable, and low‑confidence sequencing noise is prevented from influencing variant calls.

---

## 1. Base Quality Scores (Phred)

Each base in a read has a quality score encoded as an ASCII character in the SAM QUAL field.

Conversion:
Q = ord(qual_char) - 33

Interpretation:
- Q20 → ~1% error probability
- Q30 → ~0.1% error probability
- Q37 → ~0.02% error probability

We use Q20 as the default threshold for filtering.

---

## 2. Event Structure

walk_cigar() produces a list of events per read.
Each event is one of the following:

Match / mismatch:
(ref_pos, base, qual_char)

Insertion:
("INS", base, qual_char)

Deletion:
(ref_pos, "DEL")

These events are the raw material for building the pileup.

---

## 3. Building the Pileup

The pileup groups all high‑quality observations by genomic position.

Each stored observation is a dictionary:
{
    "base": base,
    "qual": q_score
}

Low‑quality bases (Q < 20) are skipped, but the read continues to be processed.
This prevents positional drift and ensures correct genomic alignment.

---

## 4. Insertions

Insertions are collected in a second pass.

Consecutive "INS" events are grouped into a single inserted sequence and anchored to the previous reference position.

Example:
ref:  A C T G
read: A C T G + "AA"

Stored as:
pos 103 → insertion "AA"

---

## 5. Deletions

Deletions are counted per reference position.

Example:
ref:  A C T G
read: A - - G

Stored as:
pos 101 → deletion length 2

---

## 6. Why Quality Filtering Matters

Without quality filtering:
- a single low‑quality miscall can create false SNPs
- insertions/deletions can be triggered by noise
- pileup becomes unreliable

With Q20 filtering:
- only confident bases contribute
- high‑quality reads still support the position even if another read is low‑quality
- evidence is weighted correctly across reads

This mirrors how real variant callers behave.

