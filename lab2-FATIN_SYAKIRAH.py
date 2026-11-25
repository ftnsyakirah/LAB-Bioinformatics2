import streamlit as st
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# NOTE: Biopython is used ONLY by Lab 1 (retrieval + basic analysis).
# The alignment algorithms in Lab 2 are written from scratch (no library).

# =====================================================================
# PAGE DESIGN
# =====================================================================
st.set_page_config(
    page_title="Protein Sequence Analyzer",
    page_icon="🧬",
    layout="wide"
)

st.markdown("""
<style>
.stApp { background-color: #ff1493; }

.block-container {
    max-width: 100% !important;
    padding-top: 2rem;
    padding-left: 3rem;
    padding-right: 3rem;
}

h1, h2, h3, p, label { color: white !important; }

.stTextInput input, .stNumberInput input {
    background-color: white !important;
    color: black !important;
}

.stTextArea textarea {
    background-color: white !important;
    color: black !important;
}

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
.stButton > button:hover { background-color: #218838; }

[data-testid="stMetricValue"] { color: white !important; }
[data-testid="stMetricLabel"] { color: white !important; }

/* Tabs */
button[data-baseweb="tab"] p { color: white !important; font-size: 18px; }
</style>
""", unsafe_allow_html=True)


# =====================================================================
# LAB 1 FUNCTIONS
# =====================================================================
def retrieve_data(id):
    Entrez.email = "ftnsyakirah@gmail.com"
    handle = Entrez.efetch(
        db="protein",
        id=id,
        rettype="fasta",
        retmode="text"
    )
    record = SeqIO.read(handle, "fasta")
    return record


def get_basic_analysis(sequence):
    seq_analysis = ProteinAnalysis(str(sequence))
    results = {
        "length": len(sequence),
        "composition": seq_analysis.count_amino_acids(),
        "molecular_weight": seq_analysis.molecular_weight(),
        "isoelectric_point": seq_analysis.isoelectric_point()
    }
    return results


# =====================================================================
# LAB 2 FUNCTIONS (alignment implemented from scratch)
# =====================================================================
def clean_sequence(text):
    """Remove FASTA header lines, whitespace and digits; return uppercase."""
    lines = [l for l in text.strip().splitlines() if not l.startswith(">")]
    seq = "".join(lines)
    return "".join(ch for ch in seq if ch.isalpha()).upper()


def align(seq1, seq2, match, mismatch, gap, mode="global"):
    """
    Pairwise alignment with a linear gap penalty.
      mode = "global" -> Needleman-Wunsch
      mode = "local"  -> Smith-Waterman

    Returns: score, aligned1, aligned2, match_line, score_matrix, path (list of (i, j))
    """
    n, m = len(seq1), len(seq2)
    local = (mode == "local")

    H = [[0] * (m + 1) for _ in range(n + 1)]       # score matrix
    P = [[""] * (m + 1) for _ in range(n + 1)]      # pointers: D, U, L, S(top)

    # Initialisation
    if not local:
        for i in range(1, n + 1):
            H[i][0] = i * gap
            P[i][0] = "U"
        for j in range(1, m + 1):
            H[0][j] = j * gap
            P[0][j] = "L"

    # Fill matrix
    best_score, best_pos = 0, (0, 0)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diag = H[i - 1][j - 1] + s
            up = H[i - 1][j] + gap
            left = H[i][j - 1] + gap

            best = max(diag, up, left)
            if local and best <= 0:
                H[i][j] = 0
                P[i][j] = "S"          # stop
            else:
                H[i][j] = best
                # tie-break priority: diagonal > up > left
                if best == diag:
                    P[i][j] = "D"
                elif best == up:
                    P[i][j] = "U"
                else:
                    P[i][j] = "L"

            if local and H[i][j] > best_score:
                best_score, best_pos = H[i][j], (i, j)

    # Starting cell of the traceback
    if local:
        i, j = best_pos
        score = best_score
    else:
        i, j = n, m
        score = H[n][m]

    # Traceback
    a1, a2, mid = [], [], []
    path = []
    while True:
        path.append((i, j))
        if local:
            if H[i][j] == 0:
                break
        else:
            if i == 0 and j == 0:
                break

        d = P[i][j]
        if d == "D":
            a1.append(seq1[i - 1])
            a2.append(seq2[j - 1])
            mid.append("|" if seq1[i - 1] == seq2[j - 1] else " ")
            i, j = i - 1, j - 1
        elif d == "U":
            a1.append(seq1[i - 1])
            a2.append("-")
            mid.append(" ")
            i -= 1
        elif d == "L":
            a1.append("-")
            a2.append(seq2[j - 1])
            mid.append(" ")
            j -= 1
        else:
            break

    a1 = "".join(reversed(a1))
    a2 = "".join(reversed(a2))
    mid = "".join(reversed(mid))
    path.reverse()

    return score, a1, a2, mid, H, path


def matrix_to_styled_df(seq1, seq2, H, path):
    """DataFrame of the score matrix with the backtrace path highlighted."""
    # unique labels are required by pandas Styler
    rows = ["-"] + [f"{c} ({i})" for i, c in enumerate(seq1, 1)]
    cols = ["-"] + [f"{c} ({j})" for j, c in enumerate(seq2, 1)]
    df = pd.DataFrame(H, index=rows, columns=cols)

    path_set = set(path)
    end_cell = path[-1] if path else None

    def highlight(_):
        styles = pd.DataFrame("", index=df.index, columns=df.columns)
        for (i, j) in path_set:
            styles.iloc[i, j] = "background-color: #ffe14d; color: black; font-weight: bold;"
        if end_cell:
            styles.iloc[end_cell[0], end_cell[1]] = (
                "background-color: #ff8c00; color: black; font-weight: bold;"
            )
        return styles

    return df.style.apply(highlight, axis=None)


def show_alignment(a1, a2, mid, width=60):
    """Print the alignment in blocks so long sequences don't overflow."""
    blocks = []
    for k in range(0, len(a1), width):
        blocks.append(
            f"Seq1: {a1[k:k + width]}\n"
            f"      {mid[k:k + width]}\n"
            f"Seq2: {a2[k:k + width]}\n"
        )
    st.code("\n".join(blocks), language=None)


def sequence_input(label, key, default):
    """Let the user type a sequence OR fetch it using the Lab 1 retrieval."""
    mode = st.radio(
        f"{label} source",
        ["Type / paste sequence", "Retrieve by UniProt ID (Lab 1)"],
        key=f"{key}_mode",
        horizontal=True
    )

    if mode == "Type / paste sequence":
        return clean_sequence(
            st.text_area(label, value=default, height=130, key=f"{key}_text")
        )

    pid = st.text_input("UniProt ID", value="P04637", key=f"{key}_id")
    if st.button("Retrieve", key=f"{key}_btn"):
        try:
            record = retrieve_data(pid)
            st.session_state[f"{key}_seq"] = str(record.seq)
            st.session_state[f"{key}_desc"] = record.description
        except Exception:
            st.error("❌ Invalid UniProt ID or retrieval failed.")

    seq = st.session_state.get(f"{key}_seq", "")
    if seq:
        st.caption(f"{st.session_state.get(f'{key}_desc', '')} — length {len(seq)}")
    return seq


# =====================================================================
# PAGES
# =====================================================================
def lab1_page():
    st.write(
        "Retrieve a protein sequence using a UniProt ID and perform basic protein analysis."
    )

    protein_id = st.text_input("Enter UniProt ID", value="P04637")

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
            st.text_area("Sequence", str(record.seq), height=200)

            st.subheader("📊 Basic Analysis")
            col1, col2, col3 = st.columns(3)
            col1.metric("Length", analysis["length"])
            col2.metric("Molecular Weight", f"{analysis['molecular_weight']:.2f}")
            col3.metric("Isoelectric Point", f"{analysis['isoelectric_point']:.2f}")

            composition_df = pd.DataFrame(
                analysis["composition"].items(),
                columns=["Amino Acid", "Count"]
            )
            st.subheader("🧪 Amino Acid Composition")
            st.dataframe(composition_df, use_container_width=True)

        except Exception:
            st.error("❌ Invalid UniProt ID or retrieval failed.")


def lab2_page():
    st.write(
        "Align two sequences using **Needleman-Wunsch (global)** or "
        "**Smith-Waterman (local)**, implemented from scratch."
    )

    # --- Sequences ---
    c1, c2 = st.columns(2)
    with c1:
        seq1 = sequence_input("Sequence 1", "s1", "GATTACA")
    with c2:
        seq2 = sequence_input("Sequence 2", "s2", "GCATGCU")

    # --- Settings ---
    st.subheader("⚙️ Alignment Settings")
    s1c, s2c, s3c, s4c = st.columns(4)
    mode_label = s1c.selectbox(
        "Alignment type", ["Global (Needleman-Wunsch)", "Local (Smith-Waterman)"]
    )
    match = s2c.number_input("Match score", value=1, step=1)
    mismatch = s3c.number_input("Mismatch score", value=-1, step=1)
    gap = s4c.number_input("Gap score (penalty)", value=-2, step=1)

    mode = "global" if mode_label.startswith("Global") else "local"

    if st.button("▶ Align"):
        if not seq1 or not seq2:
            st.error("❌ Please provide both sequences.")
            return

        score, a1, a2, mid, H, path = align(seq1, seq2, match, mismatch, gap, mode)

        # --- Score + stats ---
        st.subheader("🏆 Alignment Result")
        aligned_len = len(a1)
        matches = mid.count("|")
        gaps = a1.count("-") + a2.count("-")
        identity = (matches / aligned_len * 100) if aligned_len else 0.0

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Alignment Score", int(score))
        m2.metric("Aligned Length", aligned_len)
        m3.metric("Identity", f"{identity:.1f}%")
        m4.metric("Gaps", gaps)

        # --- Aligned sequences ---
        st.subheader("🧬 Aligned Sequences")
        if aligned_len == 0:
            st.info("No positive-scoring local alignment found.")
        else:
            show_alignment(a1, a2, mid)

        # --- Matrix + backtrace ---
        st.subheader("🧮 Score Matrix with Backtrace Path")
        if len(seq1) > 100 or len(seq2) > 100:
            st.warning(
                "Sequences are long (>100), so the matrix display is skipped "
                "to keep the app responsive. The score and alignment above are still correct."
            )
        else:
            st.caption(
                "🟨 yellow = backtrace path   |   🟧 orange = cell where the traceback ends"
            )
            st.dataframe(matrix_to_styled_df(seq1, seq2, H, path))


# =====================================================================
# MAIN PAGE
# =====================================================================
st.title("🧬 Protein Sequence Analyzer 💖")

tab1, tab2 = st.tabs(["📋 Lab 1: Protein Retrieval", "🔬 Lab 2: Pairwise Alignment"])

with tab1:
    lab1_page()

with tab2:
    lab2_page()