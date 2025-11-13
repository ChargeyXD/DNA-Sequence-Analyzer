# mutations.py

def detect_mutations(original_seq, mutated_seq):
    """
    Compare two DNA sequences and detect point mutations.
    
    Parameters:
    - original_seq: Original DNA sequence string
    - mutated_seq: Mutated/altered DNA sequence string
    
    Returns:
    - Dictionary with mutation information
    """
    # Handle length mismatch
    if len(original_seq) != len(mutated_seq):
        return {
            'has_mutations': True,
            'mutation_count': 'Length mismatch',
            'positions': [],
            'error': f"Sequences have different lengths: {len(original_seq)} vs {len(mutated_seq)}"
        }
    
    # Find all positions where sequences differ
    mutations = []
    for i, (base1, base2) in enumerate(zip(original_seq, mutated_seq)):
        if base1 != base2:
            mutations.append({
                'position': i + 1,  # 1-indexed for readability
                'original': base1,
                'mutated': base2
            })
    
    return {
        'has_mutations': len(mutations) > 0,
        'mutation_count': len(mutations),
        'positions': mutations
    }

def calculate_mutation_rate(original_seq, mutated_seq):
    """
    Calculate the mutation rate between two sequences.
    
    Returns:
    - Float: Percentage of bases that differ
    """
    if len(original_seq) != len(mutated_seq):
        return None
    
    differences = sum(1 for a, b in zip(original_seq, mutated_seq) if a != b)
    return (differences / len(original_seq)) * 100.0
