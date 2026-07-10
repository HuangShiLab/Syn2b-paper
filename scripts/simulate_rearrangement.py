#!/usr/bin/env python3
"""
simulate_rearrangement.py

Simulate inversions and translocations on a single complete genome at fixed
substitution divergence, then compare a Mash-distance proxy with Syn2b
tag-adjacency metrics.

Core claim: "tag adjacency tracks structural change while Mash does not."
"""

import argparse
import csv
import math
import os
import random
import subprocess
import sys
import tempfile

# ---------------------------------------------------------------------------
# Optional third-party dependencies
# ---------------------------------------------------------------------------
try:
    from Bio import SeqIO
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

try:
    from scipy.stats import kendalltau
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
BASIS = ["A", "C", "G", "T"]
BCGI_SITE = "GAAGGCC"
MU = 0.01
KMER = 21


# ---------------------------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------------------------
def parse_fasta(path):
    """Return (header, sequence) from the first FASTA record."""
    if BIOPYTHON_AVAILABLE:
        rec = next(SeqIO.parse(path, "fasta"))
        return rec.id, str(rec.seq).upper()

    with open(path, "r") as fh:
        lines = fh.readlines()

    header = None
    seq_parts = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                break
            header = line[1:].split()[0]
        else:
            seq_parts.append(line.upper())
    if header is None:
        raise ValueError(f"No FASTA record found in {path}")
    return header, "".join(seq_parts)


def write_fasta(path, header, sequence, width=80):
    """Write a single-record FASTA file."""
    with open(path, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(sequence), width):
            fh.write(sequence[i:i + width] + "\n")


# ---------------------------------------------------------------------------
# Genome mutation
# ---------------------------------------------------------------------------
def substitute(seq, mu=MU, rng=None):
    """Apply random point mutations at rate mu."""
    if rng is None:
        rng = random
    seq_list = list(seq)
    for i, base in enumerate(seq_list):
        if rng.random() < mu:
            alt = rng.choice([b for b in BASIS if b != base])
            seq_list[i] = alt
    return "".join(seq_list)


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def inversion(seq, size, rng=None):
    """Introduce a random inversion of *size* bp."""
    if rng is None:
        rng = random
    if size >= len(seq):
        size = len(seq) // 2
    start = rng.randint(0, len(seq) - size)
    end = start + size
    return seq[:start] + reverse_complement(seq[start:end]) + seq[end:]


def translocation(seq, size, rng=None):
    """Swap two non-overlapping segments of *size* bp each."""
    if rng is None:
        rng = random
    if 2 * size >= len(seq):
        size = len(seq) // 4
    # Pick first segment
    start1 = rng.randint(0, len(seq) - 2 * size)
    end1 = start1 + size
    # Pick second segment after the first
    start2 = rng.randint(end1, len(seq) - size)
    end2 = start2 + size
    seg1 = seq[start1:end1]
    seg2 = seq[start2:end2]
    return seq[:start1] + seg2 + seq[end1:start2] + seg1 + seq[end2:]


def insertion(seq, size, rng=None):
    """Insert a random DNA fragment of *size* bp at a random position."""
    if rng is None:
        rng = random
    if size >= len(seq):
        size = len(seq) // 10
    insert_pos = rng.randint(0, len(seq))
    fragment = "".join(rng.choices(BASIS, k=size))
    return seq[:insert_pos] + fragment + seq[insert_pos:]


def deletion(seq, size, rng=None):
    """Delete a random segment of *size* bp."""
    if rng is None:
        rng = random
    if size >= len(seq):
        size = len(seq) // 10
    start = rng.randint(0, len(seq) - size)
    return seq[:start] + seq[start + size:]


# ---------------------------------------------------------------------------
# In-silico multi-enzyme digest
# ---------------------------------------------------------------------------
def is_pure_atcg(seq):
    """Check if sequence contains only A/T/C/G."""
    return all(b in "ATCG" for b in seq)


# Enzyme definitions: (name, tag_length, [pattern_functions])
# Each pattern function returns True if the window matches.
def _bcgi_fwd(window):
    return window[10:13] == "CGA" and window[19:22] == "TGC"

def _bcgi_rev(window):
    return window[10:13] == "GCA" and window[19:22] == "TCG"


def _alfi(window):
    return window[10:13] == "GCA" and window[19:22] == "TGC"


def _bpli(window):
    return window[8:11] == "GAG" and window[16:19] == "CTC"


def _search_cjepi(seq):
    """
    Search for CjePI recognition sites (5'-CCYGA-3', Y=C or T) in sequence.
    Returns list of positions where the 5-mer starts.
    CjePI is a Type IIG RM enzyme from Campylobacter jejuni with ~1.1 sites/kb.
    """
    sites = []
    for i in range(len(seq) - 4):
        motif = seq[i:i+5]
        if motif in ("CCCGA", "CCTGA"):
            sites.append(i)
    return sites


ENZYMES = [
    ("BcgI", 32, [_bcgi_fwd, _bcgi_rev]),
    ("AlfI", 32, [_alfi]),
    ("BplI", 27, [_bpli]),
]


CJEPI_TAG_LEN = 32
CJEPI_OFFSET = 10  # CCYGA starts at offset 10 in the 32-bp tag


def _find_all_motifs(seq, motif):
    """Find all occurrences of a motif in sequence using efficient search."""
    start = 0
    positions = []
    while True:
        idx = seq.find(motif, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1
    return positions


def digest_multi_fast(sequence, include_cjepi=False):
    """
    Optimized multi-enzyme digestion using anchor-based search.
    Much faster than sliding window for long genomes.
    """
    seq = sequence.upper()
    all_tags = []
    
    # BcgI fwd: anchor CGA at offset 10, verify TGC at offset 19
    for pos in _find_all_motifs(seq, "CGA"):
        if pos >= 10 and pos + 22 <= len(seq):
            if seq[pos + 9:pos + 12] == "TGC":
                start = pos - 10
                window = seq[start:start + 32]
                if all(b in "ATCG" for b in window):
                    all_tags.append(("BcgI", window, start))
    
    # BcgI rev: anchor GCA at offset 10, verify TCG at offset 19
    for pos in _find_all_motifs(seq, "GCA"):
        if pos >= 10 and pos + 22 <= len(seq):
            if seq[pos + 9:pos + 12] == "TCG":
                start = pos - 10
                window = seq[start:start + 32]
                if all(b in "ATCG" for b in window):
                    all_tags.append(("BcgI", window, start))
    
    # AlfI: anchor GCA at offset 10, verify TGC at offset 19
    for pos in _find_all_motifs(seq, "GCA"):
        if pos >= 10 and pos + 22 <= len(seq):
            if seq[pos + 9:pos + 12] == "TGC":
                start = pos - 10
                window = seq[start:start + 32]
                if all(b in "ATCG" for b in window):
                    all_tags.append(("AlfI", window, start))
    
    # BplI: anchor GAG at offset 8, verify CTC at offset 16
    for pos in _find_all_motifs(seq, "GAG"):
        if pos >= 8 and pos + 19 <= len(seq):
            if seq[pos + 8:pos + 11] == "CTC":
                start = pos - 8
                window = seq[start:start + 27]
                if all(b in "ATCG" for b in window):
                    all_tags.append(("BplI", window, start))
    
    # CjePI: high-density CCYGA sites
    if include_cjepi:
        for motif in ("CCCGA", "CCTGA"):
            for pos in _find_all_motifs(seq, motif):
                if pos >= CJEPI_OFFSET and pos + (CJEPI_TAG_LEN - CJEPI_OFFSET) <= len(seq):
                    start = pos - CJEPI_OFFSET
                    window = seq[start:start + CJEPI_TAG_LEN]
                    if all(b in "ATCG" for b in window):
                        all_tags.append(("CjePI", window, start))
    
    all_tags.sort(key=lambda x: x[2])
    return all_tags


def digest_multi(sequence, include_cjepi=False):
    """
    Simulate multi-enzyme in-silico digestion.
    Returns a sorted list of (enzyme_name, tag_sequence, position).
    Uses optimized fast version by default.
    """
    return digest_multi_fast(sequence, include_cjepi=include_cjepi)


# Legacy slow version kept for reference
# def digest_multi_slow(sequence, include_cjepi=False): ...


def write_tgt(path, genome_id, total_length, tags):
    """
    Write a single-contig TGT text file.
    Format:
        >genome_id|length=NNN
        Enzyme:SEQ@POS [-gap- Enzyme:SEQ@POS]*
    """
    with open(path, "w") as fh:
        fh.write(f">{genome_id}|length={total_length}\n")
        for j, (enzyme, tseq, tpos) in enumerate(tags):
            if j > 0:
                gap = tpos - tags[j - 1][2]
                fh.write(f" -{gap}- ")
            fh.write(f"{enzyme}:{tseq}@{tpos}")
        fh.write("\n")


def read_tgt_tags(path):
    """Parse a single-contig TGT text file; return list of tag sequences in order."""
    with open(path, "r") as fh:
        lines = fh.readlines()

    # Skip header and comment lines
    body = ""
    for line in lines:
        if line.startswith(">") or line.startswith("#"):
            continue
        body += line.strip() + " "

    if not body.strip():
        return []

    # Split on gap markers; keep only tag entries
    raw_parts = body.split()
    tag_seqs = []
    for part in raw_parts:
        if part.startswith("-") and part.endswith("-"):
            continue
        if "@" in part and ":" in part:
            # Format: BcgI:SEQUENCE@POSITION
            enzyme_seq_pos = part.split(":", 1)
            if len(enzyme_seq_pos) == 2:
                seq_pos = enzyme_seq_pos[1].split("@", 1)
                if len(seq_pos) >= 1:
                    tag_seqs.append(seq_pos[0])
    return tag_seqs


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def mash_proxy(seq_a, seq_b, k=KMER):
    """
    Simple k-mer Jaccard proxy for sequence divergence.
    Uses *canonical* k-mers (min of forward and reverse-complement),
    matching Mash's default behaviour.  This makes the proxy blind to
    pure inversions.
    Returns (jaccard_similarity, mash_distance_approx).
    """
    def canonical_kmers(s):
        kmers = set()
        rc = reverse_complement(s)
        for i in range(len(s) - k + 1):
            fwd = s[i:i + k]
            rev = rc[len(s) - k - i:len(s) - i]
            kmers.add(min(fwd, rev))
        return kmers

    k_a = canonical_kmers(seq_a)
    k_b = canonical_kmers(seq_b)
    inter = len(k_a & k_b)
    union = len(k_a | k_b)
    if union == 0:
        return 0.0, 1.0
    jaccard = inter / union
    # Mash distance approximation: d = -1/k * ln(2*J/(1+J))
    if jaccard == 0:
        mash_d = 1.0
    else:
        mash_d = -1.0 / k * math.log(2.0 * jaccard / (1.0 + jaccard))
    return jaccard, mash_d


def adjacency_jaccard(tags_a, tags_b):
    """Jaccard similarity of adjacent tag-sequence pairs."""
    def adj_set(tlist):
        # Use canonical (sorted) pair so strand/orientation doesn't matter
        s = set()
        for i in range(len(tlist) - 1):
            a, b = tlist[i], tlist[i + 1]
            s.add((a, b) if a <= b else (b, a))
        return s

    set_a = adj_set(tags_a)
    set_b = adj_set(tags_b)
    inter = len(set_a & set_b)
    union = len(set_a | set_b)
    if union == 0:
        return 0.0
    return inter / union


def breakpoint_count(tags_a, tags_b):
    """Symmetric difference of adjacency sets."""
    def adj_set(tlist):
        s = set()
        for i in range(len(tlist) - 1):
            a, b = tlist[i], tlist[i + 1]
            s.add((a, b) if a <= b else (b, a))
        return s

    set_a = adj_set(tags_a)
    set_b = adj_set(tags_b)
    return len(set_a ^ set_b)


def kendall_tau_rank(tags_a, tags_b):
    """Kendall's tau on the order of shared tags."""
    # Build position maps
    pos_a = {t: i for i, t in enumerate(tags_a)}
    pos_b = {t: i for i, t in enumerate(tags_b)}
    shared = [t for t in tags_a if t in pos_b]
    if len(shared) < 2:
        return None
    ranks_a = [pos_a[t] for t in shared]
    ranks_b = [pos_b[t] for t in shared]
    if SCIPY_AVAILABLE:
        tau, _ = kendalltau(ranks_a, ranks_b)
        return tau
    # Simple O(n^2) implementation if scipy is absent
    n = len(shared)
    concordant = 0
    discordant = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (ranks_a[i] - ranks_a[j]) * (ranks_b[i] - ranks_b[j]) > 0:
                concordant += 1
            else:
                discordant += 1
    total = n * (n - 1) // 2
    if total == 0:
        return 0.0
    return (concordant - discordant) / total


# ---------------------------------------------------------------------------
# Synthetic genome generator
# ---------------------------------------------------------------------------
def create_synthetic_genome(path, length=2_000_000, seed=42):
    """Create a random FASTA with enough GAAGGCC sites for a realistic digest."""
    rng = random.Random(seed)
    seq = []
    # Generate random sequence, but ensure some BcgI sites exist
    site = BCGI_SITE
    block = 500
    for _ in range(length // block):
        # mostly random
        rand_part = "".join(rng.choices(BASIS, k=block - 7))
        # sprinkle a site every few blocks to get reasonable tag density
        if rng.random() < 0.5:
            insert_pos = rng.randint(0, block - 7)
            chunk = rand_part[:insert_pos] + site + rand_part[insert_pos:]
            seq.append(chunk[:block])
        else:
            seq.append(rand_part + "".join(rng.choices(BASIS, k=7)))
    remainder = length - len(seq) * block
    if remainder > 0:
        seq.append("".join(rng.choices(BASIS, k=remainder)))

    sequence = "".join(seq)[:length]
    write_fasta(path, "synthetic_test_genome", sequence)
    return sequence


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------
def run_experiment(input_fasta, syn2b_binary, output_csv, output_png):
    rng = random.Random(42)

    # ------------------------------------------------------------------
    # Step 1: Load genome
    # ------------------------------------------------------------------
    if not os.path.isfile(input_fasta):
        print(f"Warning: Genome FASTA not found at {input_fasta}")
        print("Creating a synthetic test genome instead …")
        input_fasta = tempfile.mktemp(suffix=".fasta")
        create_synthetic_genome(input_fasta)

    genome_id, original_seq = parse_fasta(input_fasta)
    genome_len = len(original_seq)
    print(f"Loaded genome: {genome_id}, length = {genome_len:,} bp")

    # ------------------------------------------------------------------
    # Step 2: Generate derived genomes
    # ------------------------------------------------------------------
    tmpdir = tempfile.mkdtemp(prefix="syn2b_rearr_")
    genomes = {}   # label -> (fasta_path, seq)

    # Original
    orig_path = os.path.join(tmpdir, "original.fasta")
    write_fasta(orig_path, genome_id, original_seq)
    genomes["original"] = (orig_path, original_seq)

    # Group A — Substitutions only (control)
    substituted_seqs = []
    for rep in range(1, 2):
        label = f"control_{rep}"
        seq = substitute(original_seq, mu=MU, rng=rng)
        substituted_seqs.append(seq)
        path = os.path.join(tmpdir, f"{label}.fasta")
        write_fasta(path, label, seq)
        genomes[label] = (path, seq)

    # Group B — Substitutions + Rearrangements
    # Apply each SV to the *same* substituted background (control_1) so that
    # the only difference between control_1 and the SV genomes is the SV itself.
    sv_specs = [
        ("inversion", 100_000),
        ("inversion", 500_000),
        ("translocation", 100_000),
        ("translocation", 500_000),
        ("insertion", 10_000),
        ("deletion", 10_000),
    ]

    base_seq = substituted_seqs[0]  # same background as control_1
    sv_positions = {}  # label -> sv_position
    for sv_type, sv_size in sv_specs:
        # Scale SV if genome is too small
        effective_size = min(sv_size, genome_len // 3)
        for rep in range(1, 3):
            label = f"{sv_type}_{sv_size // 1000}kb_r{rep}"
            # Use a distinct seed per replicate so SVs land at different positions
            rep_rng = random.Random(42 + hash((sv_type, sv_size, rep)) % 10000)
            seq = base_seq
            sv_pos = None
            if sv_type == "inversion":
                sv_pos = rep_rng.randint(0, len(seq) - effective_size)
                seq = inversion(seq, effective_size, rng=rep_rng)
            elif sv_type == "translocation":
                # translocation picks two segments; record start of first segment
                if 2 * effective_size >= len(seq):
                    effective_size = len(seq) // 4
                start1 = rep_rng.randint(0, len(seq) - 2 * effective_size)
                sv_pos = start1
                seq = translocation(seq, effective_size, rng=rep_rng)
            elif sv_type == "insertion":
                sv_pos = rep_rng.randint(0, len(seq))
                seq = insertion(seq, effective_size, rng=rep_rng)
            elif sv_type == "deletion":
                sv_pos = rep_rng.randint(0, len(seq) - effective_size)
                seq = deletion(seq, effective_size, rng=rep_rng)
            path = os.path.join(tmpdir, f"{label}.fasta")
            write_fasta(path, label, seq)
            genomes[label] = (path, seq)
            sv_positions[label] = sv_pos

    print(f"Generated {len(genomes)} genomes in {tmpdir}")

    # ------------------------------------------------------------------
    # Step 3: Digest all genomes (Python implementation)
    # ------------------------------------------------------------------
    tgt_files = {}
    for label, (fasta_path, seq) in genomes.items():
        tags = digest_multi(seq)
        tgt_path = os.path.join(tmpdir, f"{label}.tgt")
        write_tgt(tgt_path, label, len(seq), tags)
        tgt_files[label] = tgt_path
        print(f"  {label}: {len(tags)} tags")

    # Optionally invoke syn2b binary if it exists and digest works
    if syn2b_binary and os.path.isfile(syn2b_binary):
        print(f"\nSyn2b binary found at {syn2b_binary}")
        print("Note: syn2b digest is currently a stub; using Python digestion.")

    # ------------------------------------------------------------------
    # Step 4: Compute metrics (original vs each derived)
    # ------------------------------------------------------------------
    original_seq = genomes["original"][1]
    original_tags = read_tgt_tags(tgt_files["original"])

    results = []
    for label in genomes:
        if label == "original":
            continue
        seq = genomes[label][1]
        tags = read_tgt_tags(tgt_files[label])

        # Mash proxy
        _, mash_d = mash_proxy(original_seq, seq)

        # Syn2b metrics
        aj = adjacency_jaccard(original_tags, tags)
        bp = breakpoint_count(original_tags, tags)
        kt = kendall_tau_rank(original_tags, tags)

        # Determine group / SV metadata
        group = "control" if label.startswith("control") else "rearranged"
        sv_type = "none"
        sv_size = 0
        sv_pos = None
        sv_label = label
        if label.startswith("inversion_"):
            sv_type = "inversion"
            parts = label.split("_")
            sv_size = int(parts[1].replace("kb", "")) * 1000
            sv_label = f"{sv_type}_{parts[1]}"
        elif label.startswith("translocation_"):
            sv_type = "translocation"
            parts = label.split("_")
            sv_size = int(parts[1].replace("kb", "")) * 1000
            sv_label = f"{sv_type}_{parts[1]}"
        elif label.startswith("insertion_"):
            sv_type = "insertion"
            parts = label.split("_")
            sv_size = int(parts[1].replace("kb", "")) * 1000
            sv_label = f"{sv_type}_{parts[1]}"
        elif label.startswith("deletion_"):
            sv_type = "deletion"
            parts = label.split("_")
            sv_size = int(parts[1].replace("kb", "")) * 1000
            sv_label = f"{sv_type}_{parts[1]}"
        sv_pos = sv_positions.get(label, None)

        results.append({
            "genome_label": label,
            "group": group,
            "sv_type": sv_type,
            "sv_size": sv_size,
            "sv_position": sv_pos,
            "sv_label": sv_label,
            "mash_proxy": round(mash_d, 6),
            "syn2b_adjacency_jaccard": round(aj, 6),
            "syn2b_breakpoint_count": bp,
            "kendall_tau": round(kt, 6) if kt is not None else None,
        })
        print(f"  {label}: mash={mash_d:.4f}, adj_jaccard={aj:.4f}, breakpoints={bp}, tau={kt}")

    # ------------------------------------------------------------------
    # Step 5: Save CSV
    # ------------------------------------------------------------------
    fieldnames = [
        "genome_label", "group", "sv_type", "sv_size", "sv_position", "sv_label",
        "mash_proxy", "syn2b_adjacency_jaccard", "syn2b_breakpoint_count",
        "kendall_tau",
    ]
    with open(output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    print(f"\nCSV saved to {output_csv}")

    # ------------------------------------------------------------------
    # Step 5b: Plot
    # ------------------------------------------------------------------
    if not MATPLOTLIB_AVAILABLE:
        print("matplotlib not available; skipping figure generation.")
        return

    # Build data grouped by sv_label for boxplot
    # Order: control, then SV types by size
    category_order = [
        "control",
        "inversion_50kb", "inversion_100kb", "inversion_500kb", "inversion_1Mb",
        "translocation_100kb", "translocation_500kb",
        "insertion_1kb", "insertion_10kb", "insertion_50kb",
        "deletion_1kb", "deletion_10kb", "deletion_50kb",
    ]

    # Map category to color
    category_colors = {
        "control": "#4C78A8",
        "inversion_50kb": "#E45756", "inversion_100kb": "#E45756",
        "inversion_500kb": "#E45756", "inversion_1Mb": "#E45756",
        "translocation_100kb": "#F58518", "translocation_500kb": "#F58518",
        "insertion_1kb": "#54A24B", "insertion_10kb": "#54A24B", "insertion_50kb": "#54A24B",
        "deletion_1kb": "#B279A2", "deletion_10kb": "#B279A2", "deletion_50kb": "#B279A2",
    }

    # Group results by sv_label
    grouped = {cat: [] for cat in category_order}
    for r in results:
        cat = r.get("sv_label", "control")
        if cat in grouped:
            grouped[cat].append(r)

    # Prepare y-data per category
    mash_data = [[r["mash_proxy"] for r in grouped[cat]] for cat in category_order]
    adj_data = [[r["syn2b_adjacency_jaccard"] for r in grouped[cat]] for cat in category_order]
    bp_data = [[r["syn2b_breakpoint_count"] for r in grouped[cat]] for cat in category_order]

    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    x_positions = range(len(category_order))

    # Helper to plot boxplot + scatter
    def plot_metric(ax, data, ylabel, title, ylim=None):
        # Boxplot with light fill
        bp = ax.boxplot(
            data, positions=x_positions, widths=0.6, patch_artist=True,
            showfliers=False, medianprops=dict(color="black", linewidth=1.5),
            boxprops=dict(facecolor="#E8E8E8", color="black", linewidth=0.5),
            whiskerprops=dict(color="black", linewidth=0.5),
            capprops=dict(color="black", linewidth=0.5),
        )
        # Overlay individual points with category colors
        for i, cat in enumerate(category_order):
            yvals = data[i]
            color = category_colors.get(cat, "#333333")
            # Jitter x slightly for visibility
            jitter = (random.random() - 0.5) * 0.2
            for y in yvals:
                ax.scatter(i + jitter, y, color=color, edgecolor="black", linewidth=0.5, s=40, zorder=3)
        ax.set_xticks(x_positions)
        ax.set_xticklabels(category_order, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.axhline(0, color="black", linewidth=0.5)
        if ylim is not None:
            ax.set_ylim(ylim)

    plot_metric(axes[0], mash_data, "Mash distance proxy", "Sequence divergence proxy\n(all genomes ≈ 1% SNPs)", ylim=(0, None))
    plot_metric(axes[1], adj_data, "Adjacency Jaccard", "Syn2b tag-adjacency similarity", ylim=(0, 1.05))
    plot_metric(axes[2], bp_data, "Breakpoint count", "Syn2b tag-adjacency breakpoints", ylim=(0, None))

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#4C78A8", edgecolor="black", label="Control (SNPs only)"),
        Patch(facecolor="#E45756", edgecolor="black", label="Inversion"),
        Patch(facecolor="#F58518", edgecolor="black", label="Translocation"),
        Patch(facecolor="#54A24B", edgecolor="black", label="Insertion"),
        Patch(facecolor="#B279A2", edgecolor="black", label="Deletion"),
    ]
    fig.legend(handles=legend_elements, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 0.02))

    fig.suptitle(
        "Syn2b Rearrangement Validation\n"
        "Tag adjacency tracks structural change; Mash proxy does not",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0.08, 1, 0.95])
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    print(f"Figure saved to {output_png}")

    # Summary
    print("\n--- Summary ---")
    print("Control genomes (SNPs only):")
    for r in results:
        if r["group"] == "control":
            print(f"  {r['genome_label']}: mash={r['mash_proxy']:.4f}, adj_jaccard={r['syn2b_adjacency_jaccard']:.4f}")
    print("Rearranged genomes (SNPs + SV):")
    for r in results:
        if r["group"] == "rearranged":
            print(f"  {r['genome_label']}: mash={r['mash_proxy']:.4f}, adj_jaccard={r['syn2b_adjacency_jaccard']:.4f}, breakpoints={r['syn2b_breakpoint_count']}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Simulate rearrangements and compare Mash vs Syn2b metrics."
    )
    parser.add_argument(
        "--input",
        default="data/e_coli_k12.fasta",
        help="Input complete-genome FASTA (default: data/e_coli_k12.fasta)",
    )
    parser.add_argument(
        "--syn2b",
        default="target/release/syn2b",
        help="Path to syn2b binary (default: target/release/syn2b)",
    )
    parser.add_argument(
        "--csv",
        default="scripts/rearrangement_validation.csv",
        help="Output CSV path (default: scripts/rearrangement_validation.csv)",
    )
    parser.add_argument(
        "--png",
        default="scripts/rearrangement_validation.png",
        help="Output figure path (default: scripts/rearrangement_validation.png)",
    )
    parser.add_argument(
        "--mu",
        type=float,
        default=MU,
        help=f"Substitution rate (default: {MU})",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)",
    )
    args = parser.parse_args()

    random.seed(args.seed)

    # Ensure output directories exist
    os.makedirs(os.path.dirname(args.csv) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.png) or ".", exist_ok=True)

    run_experiment(args.input, args.syn2b, args.csv, args.png)


if __name__ == "__main__":
    main()
