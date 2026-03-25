file = "Basics/test.SAM"

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

            q.append((qname, flag,rname,pos,mapq,cigar))

    return q

q = SAM_parser(file)
for qname,flag,rname,pos,mapq,cigar in q:
    print(f"Read: {qname}")
    print(f"    Chromosome: {rname}")
    print(f"    Position: {pos}")
    print(f"    MAPQ: {mapq}")
    print(f"    CIGAR: {cigar}")

