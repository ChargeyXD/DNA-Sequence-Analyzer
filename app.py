# app.py
import streamlit as st
from analysis import analyze_sequence

def user_interface():
    st.title("DNA Sequence Analysis Tool")
    seq_input = st.text_area("Enter DNA Sequence", height=150)
    if st.button("Analyze"):
        if seq_input:
            analyze_sequence(seq_input)
        else:
            st.warning("Please enter a DNA sequence.")

if __name__ == "__main__":
    user_interface()
