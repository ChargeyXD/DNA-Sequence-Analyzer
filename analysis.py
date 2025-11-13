# analysis.py
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from mutations import detect_mutations

def analyze_sequence(dna_seq, compare_seq=None):
    """
    Analyze a DNA sequence for various properties.
    
    Parameters:
    - dna_seq: Primary DNA sequence string
    - compare_seq: Optional second sequence for mutation detection
    
    Returns:
    - Dictionary with analysis results
    """
    # Clean and convert to uppercase
    seq = Seq(dna_seq.upper().replace(" ", "").replace("\n", ""))
    
    # Base frequency count
    base_counts = {
        'A': seq.count('A'),
        'T': seq.count('T'),
        'C': seq.count('C'),
        'G': seq.count('G')
    }
    
    # GC Content (as percentage)
    gc_content = gc_fraction(seq) * 100.0
    
    # Reverse Complement
    rev_comp = str(seq.reverse_complement())
    
    # Build results dictionary
    results = {
        'sequence': str(seq),
        'base_counts': base_counts,
        'gc_content': gc_content,
        'reverse_complement': rev_comp
    }
    
    # Optional mutation detection
    if compare_seq:
        compare_seq_clean = compare_seq.upper().replace(" ", "").replace("\n", "")
        mutations = detect_mutations(str(seq), compare_seq_clean)
        results['mutations'] = mutations
    
    return results
