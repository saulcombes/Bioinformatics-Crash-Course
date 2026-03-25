
def parse_multi_fasta(filepath):
    sequences = {}
    current_header = None
    current_seq = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:]
                current_seq= []
            else:
                current_seq.append(line)

            if current_header:
                sequences[current_header] = "".join(current_seq)

    return sequences

def count_bases_multi(sequences_dict):
    results={}
    for header, seq in sequences_dict.items():
        counts = {
            "A": seq.count("A"),
            "T": seq.count("T"),
            "G": seq.count("G"),
            "C": seq.count("C")
        }
        results[header] = counts
    return results
