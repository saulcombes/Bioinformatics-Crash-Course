from collections import Counter
file = "Basics/test.SAM"
reference = "Basics/SAMreference.fa"

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
            qual = fields[10]

            q.append((qname, flag,rname,pos,mapq,cigar,seq,qual))

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
    for qname, flag, rname, pos, mapq, cigar, seq,qual in q:
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

def variant_caller(pileup, reference):
    # Variant rules = depth >= 2 or 3, Alt fraction >= 0.2, Alt count >= 1
    for (chrom,pos), counts in sorted(pileup.items()):
        seq_length = len(reference[chrom])
        if pos-1 >= seq_length:
            continue
        ref = reference[chrom][pos-1]
        A = counts.get("A", 0)
        C = counts.get("C", 0)
        G = counts.get("G", 0)
        T = counts.get("T", 0)
        depth = A + C + G + T
        most_prom = max(counts, key=counts.get)
        if most_prom == ref:
            pileup.pop((chrom,pos))
            continue
        alt_count = counts.get(most_prom)
        alt_fraction = alt_count/depth
        if depth >= 2 and alt_count >= 1 and alt_fraction >= 0.2:
            continue
        else:
            pileup.pop((chrom,pos))

    return pileup

def print_variant_table(pileup,reference):
    print(f"{'CHROM':<6} {'POS':<6} {'REF':<4} {'ALT':<4} {'DEPTH':<6} {'REF_COUNT':<10} {'ALT_COUNT':<10} {'ALT_FRACTION':<12}")

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
        ref_count = counts.get(ref, 0)
        alt_count = counts.get(most_prom, 0 )
        alt_fraction = alt_count/depth
        print(f"{chrom:<6} {pos:<6} {ref:<4} {most_prom:<4} {depth:<6} {ref_count:<10} {alt_count:<10} {alt_fraction:<12.2f}")

def VCF_writer(vpileup, reference):
    filename = str(input("Please enter filename: "))
    with open(f"{filename}.vcf", "w") as f:
        f.write("##fileformat=VCFv4.2\n")
    with open(f"{filename}.vcf","a") as f:
        f.write(f"{'#CHROM':<6}\t{'POS':<6}\t{'ID':<4}\t{'REF':<4}\t{'ALT':<4}\t{'QUAL':<5}\t{'FILTER':<10}\t{'INFO':<8}\n")
        for (chrom, pos), counts in sorted(vpileup.items()):
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
            ref_count = counts.get(ref, 0)
            alt_count = counts.get(most_prom, 0 )
            alt_fraction = alt_count/depth
            f.write(f"{chrom:<6}\t{pos:<6}\t{'.':<4}\t{ref:<4}\t{most_prom:<4}\t{'.':<5}\t{'PASS':<10}\t{f'DP={depth};AF={alt_fraction}':<8}\n")

def cigar_interpret(cigar,reference):
        cigarT = []
        num = ""
        for char in cigar:
            if char == "*":
                continue
            if char.isdigit():
                num += char
            else:
                cigarT.append((num, char))
                num = ""
        return cigarT

def walk_cigar(cigar_list,seq,qual,ref,start_pos):
    ref_i = start_pos
    read_i = 0
    events = []

    for length,op in cigar_list:
            length = int(length)
            if op == "M":
                for _ in range(length):
                    events.append((ref_i,seq[read_i],qual[read_i]))
                    ref_i += 1
                    read_i += 1
                
            elif op == "I":
                for _ in range(length):
                    events.append(("INS",seq[read_i],qual[read_i]))
                    read_i += 1
            
            elif op == "D":
                for _ in range(length):
                    events.append((ref_i, "DEL"))
                    ref_i += 1
            
            elif op == "S":
                read_i += length

    return events

def variant_pileup(all_events):
    pileup = {}
    insertions = {}
    deletions = {}

    for chrom, events in all_events:

        pileup.setdefault(chrom, {})
        insertions.setdefault(chrom, {})
        deletions.setdefault(chrom, {})

        for ev in events:
            pos, base, qual = ev
            qual = ord(qual)-33
            if qual <= 19:
                continue
            if pos == "INS":
                continue
            ref_pos = pos
            if base not in ("DEL",):
                pileup[chrom].setdefault(ref_pos, []).append({
                    "base": base,
                    "qual": qual
                })
            if base == "DEL":
                deletions[chrom][ref_pos] = deletions[chrom].get(ref_pos, 0) + 1

        i = 0
        while i < len(events):
            pos, base, qual = events[i]
            if pos == "INS":
                ins_seq = []
                j = i
                while j < len(events) and events[j][0] == "INS":
                    ins_seq.append(events[j][1])
                    j += 1
                if i > 0:
                    prev_ref_pos = events[i - 1][0]
                    insertions[chrom].setdefault(prev_ref_pos, []).append("".join(ins_seq))

                i = j
            else:
                i += 1

    return pileup, insertions, deletions

def call_variants_complex(pileup,insertions,deletions,reference):
    variants = []

    for chrom in pileup:
        for pos in sorted(pileup[chrom].keys()):
            ref_base = reference[chrom][pos-1]
            base_entries = pileup[chrom][pos]
            bases = [entry["base"] for entry in base_entries]
            counts = Counter(bases)
            ref_count = counts.get(ref_base,0)
            alt_counts = {b:c for b,c in counts.items() if b != ref_base}
            if not alt_counts:
                pass
            else:
                alt_base,alt_count = max(alt_counts.items(),key=lambda x: x[1])
                if alt_count> ref_count:
                    variants.append((chrom,pos,ref_base,alt_base))
            if chrom in insertions and pos in insertions[chrom]:
                for ins in insertions[chrom][pos]:
                    ref=ref_base
                    alt=ref_base+ins
                    variants.append((chrom,pos,ref,alt))
            if chrom in deletions and pos in deletions[chrom]:
                del_len = deletions[chrom][pos]
                deleted_bases = ''.join(reference[chrom][pos+i] for i in range(1,del_len+1))
                ref=ref_base+deleted_bases
                alt = ref_base
                variants.append((chrom,pos,ref,alt))
        return variants

def VCF_writer_complex(vpileup, reference):
    filename = str(input("Please enter filename: "))
    with open(f"{filename}.vcf", "w") as f:
        f.write("##fileformat=VCFv4.2\n")
    with open(f"{filename}.vcf","a") as f:
        f.write(f"{'#CHROM':<6}\t{'POS':<6}\t{'ID':<4}\t{'REF':<4}\t{'ALT':<4}\t{'QUAL':<5}\t{'FILTER':<10}\t{'INFO':<8}\n")
        for (chrom,pos,ref,alt) in sorted(vpileup):
            seq_len = len(reference[chrom])
            if pos - 1 >= seq_len:
                continue
            ref = reference[chrom][pos-1]
            f.write(f"{chrom:<6}\t{pos:<6}\t{'.':<4}\t{ref:<4}\t{alt:<4}\t{'.':<5}\t{'PASS':<10}\t{'.':<8}\n")

q = SAM_parser(file)
for qname,flag,rname,pos,mapq,cigar,seq,qual in q:
    print(f"Read: {qname}")
    print(f"    Chromosome: {rname}")
    print(f"    Position: {pos}")
    print(f"    MAPQ: {mapq}")
    print(f"    CIGAR: {cigar}")
    print(f"    Sequence: {seq}")
    print(f"    Qual: {qual}")
    print(f"    Calc: {len(seq)} / {len(qual)}")


pileup,reference = pileup_formation(q,reference)
print_pileup_table(pileup,reference)
vpileup = variant_caller(pileup,reference)
print("")
print_variant_table(vpileup,reference)
all_events = []
for qname,flag,rname,pos,mapq,cigar,seq,qual in q:
        read = cigar_interpret(cigar,reference)
        if flag ==4:
            all_events.append((rname,[]))
            continue
        if flag == 16:
            comp = str.maketrans("ACGTacgt","TGCAtgca")
            seq = seq.translate(comp)[::-1]
            qual = qual[::-1]
        events = walk_cigar(read,seq,qual,reference,pos)
        all_events.append((rname,events))



pileup,insertions,deletions = variant_pileup(all_events)
variants = call_variants_complex(pileup,insertions,deletions,reference)
VCF_writer_complex(variants, reference)
