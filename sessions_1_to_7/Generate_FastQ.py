import random

outfile = "testdata_10000.fastq"

N = 10000

adapters = [
    "ACGTACGTACGT",    
    "TTGACCGT",   
    "GATCTAGCTA"       
]

bases = ["A", "C", "G", "T"]

def random_seq(length):
    return "".join(random.choice(bases) for _ in range(length))

def random_qual(length):
    return "I" * length

with open(outfile, "w") as f:
    for i in range(1, N + 1):

        insert_len = random.randint(50, 150)
        insert = random_seq(insert_len)

        if random.random() < 0.4:
            adapter = random.choice(adapters)
            seq = insert + adapter
        else:
            seq = insert

        qual = random_qual(len(seq))

        f.write(f"@read{i}\n")
        f.write(seq + "\n")
        f.write("+\n")
        f.write(qual + "\n")