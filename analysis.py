# analysis.py
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from mutations import detect_mutations  # Import mutation detection

def analyze_sequence(dna_seq):
    seq = Seq(dna_seq.upper())

    # Base frequency
    base_counts = {base: seq.count(base) for base in "ATCG"}

    # GC Content
    gc_content = GC(seq)

    # Reverse Complement
    rev_comp = seq.reverse_complement()

    # Call mutation detection (you can pass seq or string as needed)
    mutations = detect_mutations(seq)

    # Here you could modify this function to return all these results
    # For now, just print or you can update app.py to display
    print("Base Counts:", base_counts)
    print("GC Content:", gc_content)
    print("Reverse Complement:", rev_comp)
    print("Mutations:", mutations)
