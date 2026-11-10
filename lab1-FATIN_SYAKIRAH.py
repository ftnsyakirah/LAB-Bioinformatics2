import streamlit as st
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd

# PAGE DESIGN
st.set_page_config(
    page_title="Protein Sequence Analyzer",
    page_icon="🧬",
    layout="wide"
)

# CUSTOM CSS

st.markdown("""
<style>

.stApp {
    background-color: #ff1493;
}

/* Make page use full width */
.block-container {
    max-width: 100% !important;
    padding-top: 2rem;
    padding-left: 3rem;
    padding-right: 3rem;
}

/* Main text */
h1, h2, h3, p, label {
    color: white !important;
}

/* Input box */
.stTextInput input {
    background-color: white !important;
    color: black !important;
}

/* Text Area (Protein Sequence) */
.stTextArea textarea {
    background-color: white !important;
    color: black !important;
}

/* Retrieve Button */
.stButton > button {
    background-color: #28a745;
    color: white !important;
    font-weight: bold;
    font-size: 18px;
    border-radius: 10px;
    width: 220px;
    height: 50px;
    border: none;
}

.stButton > button:hover {
    background-color: #218838;
}

/* Metrics */
[data-testid="stMetricValue"] {
    color: white !important;
}

[data-testid="stMetricLabel"] {
    color: white !important;
}

</style>
""", unsafe_allow_html=True)

# FUNCTIONS retrieve_data(id)

def retrieve_data(id):

    Entrez.email = "ftnsyakirah@gmail.com"

    handle = Entrez.efetch(
        db = "protein",
        id = id,
        rettype = "fasta",
        retmode = "text"
    )

    record = SeqIO.read(handle, "fasta")

    return record

from Bio.SeqUtils.ProtParam import ProteinAnalysis

# FUNCTIONS get_basic_analysis(sequence)

def get_basic_analysis(sequence):

    seq_analysis = ProteinAnalysis(str(sequence))

    results = {
        "length": len(sequence),
        "composition": seq_analysis.count_amino_acids(),
        "molecular_weight": seq_analysis.molecular_weight(),
        "isoelectric_point": seq_analysis.isoelectric_point()
    }

    return results

# MAIN PAGE

st.title("🧬 Protein Sequence Analyzer 💖")

st.write(
    "Retrieve a protein sequence using a UniProt ID and perform basic protein analysis."
)

protein_id = st.text_input(
    "Enter UniProt ID",
    value="P04637"
)

if st.button("🔍 Retrieve Protein"):

    try:

        record = retrieve_data(protein_id)
        analysis = get_basic_analysis(record.seq)
        st.success("Protein retrieved successfully! ✅")
        st.subheader("📋 Protein Information")
        st.write("**Name:**", record.name)
        st.write("**Description:**", record.description)
        st.write("**Sequence Length:**", analysis["length"])
        st.subheader("🧬 Protein Sequence")

        st.text_area(
            "Sequence",
            str(record.seq),
            height=200
        )

        st.subheader("📊 Basic Analysis")
        col1, col2, col3 = st.columns(3)

        col1.metric(
            "Length",
            analysis["length"]
        )

        col2.metric(
            "Molecular Weight",
            f"{analysis['molecular_weight']:.2f}"
        )

        col3.metric(
            "Isoelectric Point",
            f"{analysis['isoelectric_point']:.2f}"
        )

        composition_df = pd.DataFrame(
            analysis["composition"].items(),
            columns=["Amino Acid", "Count"]
        )

        st.subheader("🧪 Amino Acid Composition")

        st.dataframe(
            composition_df,
            use_container_width=True
        )

    except Exception:
        st.error("❌ Invalid UniProt ID or retrieval failed.")