# app.py
"""
DNA Sequence Analysis Tool (Streamlit)

This file keeps ALL your original logic exactly the same,
with ONLY ONE CHANGE:
➡ Replaced the matplotlib bar chart with Streamlit's st.bar_chart()
"""

from typing import Dict, List, Optional
import streamlit as st
import json
import io

# Try to import user's analysis.analyze_sequence if available, otherwise use fallback
try:
    from analysis import analyze_sequence  # type: ignore
except Exception:
    def analyze_sequence(seq: str, compare_seq: Optional[str] = None) -> Dict:
        """Fallback analysis: computes base counts, GC content, reverse complement,
        and simple mutation detection (position-by-position) if compare_seq provided."""
        seq = seq.upper().replace("\n", "").replace("\r", "").strip()
        base_counts = {b: seq.count(b) for b in ("A", "T", "C", "G")}
        total_bases = sum(base_counts.values()) or 1
        gc_content = (base_counts["G"] + base_counts["C"]) / total_bases * 100
        complement_map = str.maketrans("ATCG", "TAGC")
        reverse_complement = seq.translate(complement_map)[::-1]

        results = {
            "sequence": seq,
            "base_counts": base_counts,
            "gc_content": gc_content,
            "reverse_complement": reverse_complement,
        }

        if compare_seq:
            cmp = compare_seq.upper().replace("\n", "").replace("\r", "").strip()
            min_len = min(len(seq), len(cmp))
            positions = []

            for i in range(min_len):
                if seq[i] != cmp[i]:
                    positions.append({"position": i + 1, "original": seq[i], "mutated": cmp[i]})

            has_mutations = len(positions) > 0 or len(seq) != len(cmp)

            results["mutations"] = {
                "has_mutations": has_mutations,
                "mutation_count": len(positions) + (1 if len(seq) != len(cmp) else 0),
                "positions": positions,
                "lengths": {"seq": len(seq), "compare_seq": len(cmp)}
            }

        return results


# ---------- Helper functions ----------
def validate_and_clean_sequence(seq: str, ignore_non_atcg: bool = False) -> str:
    """Uppercase, strip whitespace, optionally remove non A/T/C/G characters."""
    if seq is None:
        return ""
    s = seq.upper().replace(" ", "").replace("\n", "").replace("\r", "").strip()

    if ignore_non_atcg:
        cleaned = "".join(ch for ch in s if ch in ("A", "T", "C", "G"))
        return cleaned

    invalid = [ch for ch in s if ch not in ("A", "T", "C", "G")]
    if invalid:
        raise ValueError(f"Invalid characters found: {sorted(set(invalid))}")
    return s


def format_mutation_table(mutations: Dict) -> str:
    if not mutations:
        return "No mutation information."
    if not mutations.get("has_mutations"):
        return "No mutations detected. Sequences are identical."

    lines = ["| Position | Original | Mutated |", "|---:|:---:|:---:|"]
    for m in mutations.get("positions", []):
        lines.append(f"| {m['position']} | {m['original']} | {m['mutated']} |")

    if mutations.get("lengths") and mutations["lengths"]["seq"] != mutations["lengths"]["compare_seq"]:
        lines.append(
            f"\n**Note:** sequences differ in length: input={mutations['lengths']['seq']} vs compare={mutations['lengths']['compare_seq']}"
        )

    return "\n".join(lines)


def results_to_json(results: Dict) -> str:
    return json.dumps(results, indent=2)


# ---------- Streamlit UI ----------
st.set_page_config(page_title="DNA Sequence Analysis Tool", page_icon="🧬", layout="wide")


def user_interface():

    st.title("🧬 DNA Sequence Analysis Tool")
    st.markdown("Analyze DNA sequences for base composition, GC content, reverse complement, and mutations.")

    # --- Sidebar controls ---
    st.sidebar.header("Settings & Input")
    ignore_non_atcg = st.sidebar.checkbox("Ignore non-A/T/C/G characters (strip them)", value=True)
    show_plot = st.sidebar.checkbox("Show base composition plot", value=True)
    display_progress_as_percent = st.sidebar.checkbox("Show GC progress bar", value=True)

    st.sidebar.markdown("---")
    st.sidebar.subheader("Quick samples")

    if st.sidebar.button("Insert sample: short"):
        st.session_state["seq_input"] = "ATGCGTACGAT"

    if st.sidebar.button("Insert sample: longer"):
        st.session_state["seq_input"] = "ATGCGTACGATGCTAGCTAGCTGACTGATCGATCGATCGTTAGC"

    st.sidebar.markdown("---")
    st.sidebar.info("You can paste a sequence, upload a FASTA/plain text file, or compare with a second sequence.")

    # ---------- Main Input Area ----------
    col1, col2 = st.columns([2, 1])

    with col1:

        st.subheader("Input DNA Sequence")

        seq_input = st.text_area(
            "Enter DNA Sequence (A, T, C, G only)",
            key="seq_input",
            height=200,
            placeholder="Example: ATCGATCGATCG"
        )

        # File upload
        uploaded = st.file_uploader(
            "Or upload a FASTA / text file with sequence",
            type=["txt", "fasta", "fa"]
        )

        if uploaded is not None:
            try:
                raw = uploaded.getvalue().decode("utf-8")
            except Exception:
                raw = uploaded.getvalue().decode("latin-1")

            lines = [l.strip() for l in raw.splitlines() if not l.startswith(">")]
            from_file_seq = "".join(lines)

            st.info("Uploaded file processed and appended.")

            seq_input = (seq_input or "") + ("\n" if seq_input else "") + from_file_seq
            st.session_state["seq_input"] = seq_input

        st.subheader("Mutation Detection (Optional)")

        compare_seq = st.text_area(
            "Enter a second sequence to compare for mutations",
            key="compare_seq_input",
            height=120,
            placeholder="Example: ATCGATCGATCG"
        )

        analyze_clicked = st.button("Analyze Sequence")

    with col2:

        st.subheader("Actions")

        if st.button("Clear inputs"):
            st.session_state["seq_input"] = ""
            st.session_state["compare_seq_input"] = ""
            st.experimental_rerun()

        st.markdown("---")
        st.subheader("Output options")

        download_results_placeholder = st.empty()

    # ---------- When Analyze is clicked ----------
    if analyze_clicked:

        try:
            cleaned_seq = validate_and_clean_sequence(seq_input or "", ignore_non_atcg=ignore_non_atcg)
        except ValueError as e:
            st.error(f"Input validation failed: {e}")
            st.stop()

        cleaned_compare = None

        if compare_seq:
            try:
                cleaned_compare = validate_and_clean_sequence(compare_seq or "", ignore_non_atcg=ignore_non_atcg)
            except ValueError as e:
                st.error(f"Compare sequence validation failed: {e}")
                st.stop()

        if not cleaned_seq:
            st.warning("⚠️ Please enter a valid DNA sequence.")
        else:
            results = analyze_sequence(cleaned_seq, cleaned_compare if cleaned_compare else None)

            display_results(
                results,
                show_plot=show_plot,
                progress_bar=display_progress_as_percent,
                download_container=download_results_placeholder
            )


def display_results(results: Dict, show_plot: bool = True, progress_bar: bool = True, download_container=None):

    st.success("✅ Analysis Complete!")

    # --------------------------------------------------
    # INPUT SEQUENCE
    # --------------------------------------------------
    st.subheader("Input Sequence")
    st.code(results.get("sequence", ""))


    # --------------------------------------------------
    # BASE FREQUENCY
    # --------------------------------------------------
    st.subheader("Base Frequency")

    base_counts = results.get("base_counts", {"A": 0, "T": 0, "C": 0, "G": 0})

    colA, colT, colC, colG = st.columns(4)
    colA.metric("Adenine (A)", base_counts.get("A", 0))
    colT.metric("Thymine (T)", base_counts.get("T", 0))
    colC.metric("Cytosine (C)", base_counts.get("C", 0))
    colG.metric("Guanine (G)", base_counts.get("G", 0))


    # --------------------------------------------------
    # GC CONTENT
    # --------------------------------------------------
    st.subheader("GC Content")

    gc = results.get("gc_content", 0.0)
    st.metric("GC Percentage", f"{gc:.2f}%")

    if progress_bar:
        st.progress(min(max(gc / 100.0, 0.0), 1.0))


    # --------------------------------------------------
    # STREAMLIT NATIVE BAR CHART (REPLACED MATPLOTLIB)
    # --------------------------------------------------
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


    # --------------------------------------------------
    # REVERSE COMPLEMENT
    # --------------------------------------------------
    st.subheader("Reverse Complement")
    st.code(results.get("reverse_complement", ""))


    # --------------------------------------------------
    # MUTATIONS
    # --------------------------------------------------
    if "mutations" in results:

        st.subheader("Mutation Analysis")
        mutations = results["mutations"]

        if mutations.get("has_mutations"):
            st.warning(f"🧬 {mutations.get('mutation_count', 0)} mutation(s) detected!")
            st.markdown(format_mutation_table(mutations))
        else:
            st.info("No mutations detected. Sequences are identical.")


    # --------------------------------------------------
    # DOWNLOAD RESULTS
    # --------------------------------------------------
    if download_container is not None:
        json_str = results_to_json(results)
        buf = io.BytesIO(json_str.encode("utf-8"))

        download_container.download_button(
            label="Download results (JSON)",
            data=buf,
            file_name="dna_analysis_results.json",
            mime="application/json"
        )


    # --------------------------------------------------
    # UTILITIES
    # --------------------------------------------------
    st.markdown("---")
    st.subheader("Utilities")
    st.write("Copy sequence:")
    st.code(results.get("sequence", ""))


if __name__ == "__main__":
    user_interface()
