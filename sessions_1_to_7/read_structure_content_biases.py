from itertools import zip_longest
from collections import Counter
import pandas as pd
import numpy as np

def read_length_distribution(sequences):
    seqlen = []
    for line in sequences:
        seqlen.append(len(line))

    return seqlen

def n_content_per_cycle(sequences):
    padded = list(zip_longest(*sequences, fillvalue="N"))
    n_percentages = []

    for cycle in padded:
        total = len(cycle)
        n_count = cycle.count("N")
        n_percent = (n_count / total) * 100
        n_percentages.append(n_percent)

    return n_percentages

def overrepresented_sequences(seqs, lengths =[ 12,10,8 ], min_fraction = 0.01):
    seqend = []
    for L in lengths:
        suffixes = [s[-L:]for s in seqs if len(s) >= L]
        counts = Counter(suffixes)

        total=len(seqs)
        threshold = total*min_fraction

        for seq,count in counts.items():
            if count >= threshold:
                percentage = (count/total)*100
                seqend.append((seq,percentage))

    print("")
    df = pd.DataFrame(seqend,columns=["Sequence", "Percentage Incidence"])
    df.index = df.index +1
    df = df.reset_index().rename(columns={"index": "Rank"})
    return df

def k_mer_content(sequences):
    kmerseqs = []
    for line in sequences:
        end = len(line) - 4
        x=0
        y=5
        kmerseqs.append(line[x:y])
        while True:
          x += 1
          y += 1 
          if x == end:
              break
          kmerseqs.append(line[x:y])
    
    counts = Counter(kmerseqs)
    kres = [{kmer:kcount} for kmer,kcount in counts.items()]
    kres1=kres
    percents= []
    flags = []
    total = len(kmerseqs)
    for kmer,kcount in counts.items():
        if kcount == 1:
            kres.remove({kmer:kcount})
        else:
            percent = (kcount/total)*100
            percents.append(percent)
            if percent >= 0.1:
                flags.append("*")
            else:
                flags.append(".")

    kresman = [(list(d.keys())[0], list(d.values())[0]) for d in kres]
    df = pd.DataFrame(kresman,columns=["K-mer","Count"])
    df["Percentage Incidence(%)"]=percents
    df["Flags"]=flags

    return df,counts

##ACGTACGTACGT

from collections import Counter

from collections import Counter

def adapter_candidates(seqs,
                       scan_tail=30,
                       adapter_lengths=(8, 10, 12),
                       min_support=50):
    """
    Very simple adapter inference:
    - look only in the last `scan_tail` bases
    - for each length L, count all k-mers of length L in that region
    - return the most frequent k-mer(s) per L above a support threshold
    """

    adapters = []

    for L in adapter_lengths:
        counts = Counter()

        for seq in seqs:
            tail = seq[-scan_tail:]
            if len(tail) < L:
                continue
      
            for i in range(len(tail) - L + 1):
                kmer = tail[i:i+L]
                counts[kmer] += 1

        if not counts:
            continue

        kmer, support = counts.most_common(1)[0]
        if support >= min_support:
            adapters.append(kmer)

    return adapters

def sequence_trimming(seqs, quals, adapter_candidates):
    adapter_candidates = sorted(adapter_candidates, key=len, reverse=True)

    trimmed_seqs = []
    trimmed_quals = []
    lengths_removed = []

    for seq, qual in zip(seqs, quals):
        removed = 0

        for adapt in adapter_candidates:
            if seq.endswith(adapt):
                removed = len(adapt)
                break  

        if removed > 0:
            trimmed_seqs.append(seq[:-removed])
            trimmed_quals.append(qual[:-removed])
        else:
            trimmed_seqs.append(seq)
            trimmed_quals.append(qual)

        lengths_removed.append(removed)

    return trimmed_seqs, trimmed_quals, lengths_removed
        



