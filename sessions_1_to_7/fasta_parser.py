
def count_bases_from_fasta(filepath):
   with open(filepath, "r") as f:
      sequence = ""
      for line in f:
         if line.startswith(">"):
            continue
         sequence += line.strip()

   counts = {"A":0,"T":0,"G":0,"C":0}
   for base in sequence:
      if base in counts:
         counts[base] += 1

   return sequence,counts
   
seq,counts = count_bases_from_fasta("test_sequence.fasta")

print("Sequence length:", len(seq))
print(counts)

