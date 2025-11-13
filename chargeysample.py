# app.py
import streamlit as st
from analysis import analyze_sequence

def user_interface():
    st.image(image="design\\dna")
    st.title("DNA Sequence Analysis Tool",anchor='')
    st.markdown("Analyze DNA sequences for base composition, GC content, reverse complement, and mutations.")


if __name__ == "__main__":
    user_interface()