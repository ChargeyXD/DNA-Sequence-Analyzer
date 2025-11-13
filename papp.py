# app.py
import streamlit as st
from analysis import analyze_sequence

def user_interface():
    st.title("🧬 DNA Sequence Analysis Tool")
    st.markdown("Analyze DNA sequences for base composition, GC content, reverse complement, and mutations.")
    
    # Input section
    st.subheader("Input DNA Sequence")
    seq_input = st.text_area(
        "Enter DNA Sequence (A, T, C, G only)", 
        height=150,
        placeholder="Example: ATCGATCGATCG"
    )
    
    # Optional: mutation detection
    st.subheader("Mutation Detection (Optional)")
    compare_seq = st.text_area(
        "Enter a second sequence to compare for mutations",
        height=100,
        placeholder="Leave empty to skip mutation detection"
    )
    
    # Analyze button
    if st.button("Analyze Sequence"):
        if seq_input:
            results = analyze_sequence(seq_input, compare_seq if compare_seq else None)
            display_results(results)
        else:
            st.warning("⚠️ Please enter a DNA sequence to analyze.")

def display_results(results):
    st.success("✅ Analysis Complete!")
    
    # Display input sequence
    st.subheader("Input Sequence")
    st.code(results['sequence'], language=None)
    
    # Base frequency
    st.subheader("Base Frequency")
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Adenine (A)", results['base_counts']['A'])
    col2.metric("Thymine (T)", results['base_counts']['T'])
    col3.metric("Cytosine (C)", results['base_counts']['C'])
    col4.metric("Guanine (G)", results['base_counts']['G'])
    
    # GC Content
    st.subheader("GC Content")
    st.metric("GC Percentage", f"{results['gc_content']:.2f}%")
    st.progress(results['gc_content'] / 100)
    
    # Reverse Complement
    st.subheader("Reverse Complement")
    st.code(results['reverse_complement'], language=None)
    
    # Mutations (if available)
    if 'mutations' in results and results['mutations']:
        st.subheader("Mutation Analysis")
        if results['mutations']['has_mutations']:
            st.warning(f"🧬 {results['mutations']['mutation_count']} mutation(s) detected!")
            for mutation in results['mutations']['positions']:
                st.write(f"Position {mutation['position']}: {mutation['original']} → {mutation['mutated']}")
        else:
            st.info("No mutations detected. Sequences are identical.")

if __name__ == "__main__":
    user_interface()
