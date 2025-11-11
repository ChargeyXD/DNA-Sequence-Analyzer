# analysis.py
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction  # Import gc_fraction instead of GC

def analyze_sequence(dna_seq):
    seq = Seq(dna_seq.upper())

    # Base frequency
    base_counts = {base: seq.count(base) for base in "ATCG"}

    # GC Content
    gc_content = gc_fraction(seq) * 100  # gc_fraction returns a fractionajskda
    print("GC Content:", gc_content)

    # Reverse Complement
    rev_comp = seq.reverse_complement()

    # Call mutation detection
    # mutations = detect_mutations(seq)

    # Return or print results
    return {
        "base_counts": base_counts,
        "gc_content": gc_content,
        "reverse_complement": str(rev_comp),
        # "mutations": mutations
    }
