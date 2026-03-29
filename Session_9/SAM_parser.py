file = "test.SAM"
reference = "SAMreference.fa"

def SAM_parser(filepath):
    q = []


    with open(filepath,"r") as f:
        for line in f:
            line = line.strip()
            fields = line.strip().split("\t")
            if line.startswith("@"):
                continue
            if not line:
                continue
            if len(fields) < 11:
                continue
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]
            seq = fields[9]

            q.append((qname, flag,rname,pos,mapq,cigar,seq))

    return q

from collections import defaultdict

def pileup_formation(q, referenceog):
    reference = {}
    current_chr = None
    seq_chunks = []

    with open(referenceog, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_chr:
                    reference[current_chr] = "".join(seq_chunks)
                current_chr = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)

        # Save last chromosome
        if current_chr:
            reference[current_chr] = "".join(seq_chunks)

    # Pileup structure
    pileup = defaultdict(lambda: defaultdict(int))

    # Process reads
    for qname, flag, rname, pos, mapq, cigar, seq in q:
        if rname == "*":
            continue  # unmapped

        pos = int(pos)  # SAM is 1-based

        for i, base in enumerate(seq):
            genome_pos = pos + i
            pileup[(rname, genome_pos)][base] += 1

    return pileup,reference
    
def print_pileup_table(pileup,reference):
    print("CHROM\tPOS\tRef\tA\tC\tG\tT\tDEPTH\tFLAG")
    for (chrom, pos), counts in sorted(pileup.items()):
        seq_len = len(reference[chrom])
        if pos - 1 >= seq_len:
            continue
        ref = reference[chrom][pos-1]
        A = counts.get("A", 0)
        C = counts.get("C", 0)
        G = counts.get("G", 0)
        T = counts.get("T", 0)
        depth = A + C + G + T
        most_prom = max(counts, key=counts.get)
        flag="."
        if most_prom != ref:
            flag = "*"
        print(f"{chrom}\t{pos}\t{ref}\t{A}\t{C}\t{G}\t{T}\t{depth}\t{flag}")


q = SAM_parser(file)
for qname,flag,rname,pos,mapq,cigar,seq in q:
    print(f"Read: {qname}")
    print(f"    Chromosome: {rname}")
    print(f"    Position: {pos}")
    print(f"    MAPQ: {mapq}")
    print(f"    CIGAR: {cigar}")
    print(f"    Sequence: {seq}")

pileup,reference = pileup_formation(q,reference)
print_pileup_table(pileup,reference)
