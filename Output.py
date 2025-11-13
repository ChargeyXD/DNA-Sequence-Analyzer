# pages/Output.py
"""
Output page for DNA Sequence Analysis Tool (Streamlit multipage)
Reads results from st.session_state['dna_results'] and displays them.
"""

import streamlit as st
import io
import json

st.set_page_config(page_title="DNA Analysis - Output", page_icon="🧬", layout="wide")

# Re-declare helper functions used by the display (kept identical to input page functions)
def format_mutation_table(mutations: dict) -> str:
    if not mutations:
        return "No mutation information."
    if not mutations.get("has_mutations"):
        return "No mutations detected. Sequences are identical (within compared region)."
    lines = ["| Position | Original | Mutated |", "|---:|:---:|:---:|"]
    for m in mutations.get("positions", []):
        lines.append(f"| {m['position']} | {m['original']} | {m['mutated']} |")
    if mutations.get("lengths") and mutations["lengths"]["seq"] != mutations["lengths"]["compare_seq"]:
        lines.append(
            f"\n**Note:** sequences differ in length: input={mutations['lengths']['seq']} vs compare={mutations['lengths']['compare_seq']}"
        )
    return "\n".join(lines)

def results_to_json(results: dict) -> str:
    return json.dumps(results, indent=2)

# Page header
st.title("📤 Output — DNA Analysis Results")
st.markdown("This page displays the results of the most recent analysis (from the Input page).")

dna_saved = st.session_state.get('dna_results')

if not dna_saved:
    st.warning("No analysis results found. Go to the Input page, enter a sequence and click Analyze.")
    st.markdown("[← Back to Input page](/)", unsafe_allow_html=True)
    st.stop()

results = dna_saved.get("results", {})
show_plot = dna_saved.get("show_plot", True)

st.success("Results loaded from session state.")

# Display results (kept same structure as display_results)
st.subheader("Input Sequence")
st.code(results.get("sequence", ""), language=None)

st.subheader("Base Frequency")
base_counts = results.get("base_counts", {"A": 0, "T": 0, "C": 0, "G": 0})
colA, colT, colC, colG = st.columns(4)
colA.metric("Adenine (A)", base_counts.get("A", 0))
colT.metric("Thymine (T)", base_counts.get("T", 0))
colC.metric("Cytosine (C)", base_counts.get("C", 0))
colG.metric("Guanine (G)", base_counts.get("G", 0))

st.subheader("GC Content")
gc = results.get("gc_content", 0.0)
st.metric("GC Percentage", f"{gc:.2f}%")
st.progress(min(max(gc / 100.0, 0.0), 1.0))

# Streamlit native bar chart (same as in the input page)
if show_plot:
    try:
        import pandas as pd
        labels = ["A", "T", "C", "G"]
        counts = [base_counts.get(lbl, 0) for lbl in labels]
        df = pd.DataFrame({"Base": labels, "Count": counts}).set_index("Base")
        st.subheader("Base Composition Chart")
        st.bar_chart(df)
    except Exception as e:
        st.warning(f"Could not render chart: {e}")

st.subheader("Reverse Complement")
st.code(results.get("reverse_complement", ""), language=None)

st.subheader("Mutation Analysis")
if "mutations" in results:
    mutations = results["mutations"]
    if mutations.get("has_mutations"):
        st.warning(f"🧬 {mutations.get('mutation_count', 0)} mutation(s) detected!")
        st.markdown(format_mutation_table(mutations))
    else:
        st.info("No mutations detected. Sequences are identical (within compared region).")
else:
    st.info("No mutation comparison was provided in the analysis.")

# Download results button
st.markdown("---")
json_bytes = json.dumps(results, indent=2).encode("utf-8")
st.download_button("Download results (JSON)", data=io.BytesIO(json_bytes), file_name="dna_analysis_results.json", mime="application/json")

st.markdown("---")
st.markdown("[← Back to Input page](/)", unsafe_allow_html=True)
