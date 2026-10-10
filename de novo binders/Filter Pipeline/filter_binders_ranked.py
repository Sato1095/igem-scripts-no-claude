#!/usr/bin/env python3
"""
filter_binders.py — Filter and rank de novo binder designs from an
RFdiffusion → ProteinMPNN → AlphaFold pipeline.

Reads a FASTA whose headers carry per-design metrics, e.g.:

    >design:1 n:2|mpnn:1.288|plddt:0.910|i_ptm:0.806|i_pae:4.017|rmsd:0.757
    <target_seq>/<peptide_seq>

Applies a metric filter, de-duplicates, scores by a composite z-score,
flags backbone self-consistency and (optionally) a target anchor motif,
and writes a ranked CSV shortlist.

Stdlib only — no external dependencies.

Examples
--------
    # Standard tier, default output
    python filter_binders.py design_fasta.txt

    # Stringent tier, custom output, MDM2 p53-anchor motif scan on
    python filter_binders.py design_fasta.txt -o top.csv \
        --tier stringent --anchor-motif

    # Fully custom thresholds
    python filter_binders.py design_fasta.txt \
        --iptm 0.75 --ipae 8 --rmsd 1.8 --plddt 0.82
"""

import argparse
import csv
import re
import statistics
import sys
from collections import defaultdict

# --------------------------------------------------------------------------
# Metric direction: +1 = higher is better, -1 = lower is better.
# Add/remove keys here if your headers carry different fields.
# --------------------------------------------------------------------------
METRIC_DIRECTION = {"iptm": +1, "ipae": -1, "rmsd": -1, "plddt": +1}

# Preset filter tiers: (iptm_min, ipae_max, rmsd_max, plddt_min)
TIERS = {
    "stringent": (0.80, 5.0, 1.5, 0.85),
    "standard":  (0.70, 10.0, 2.0, 0.80),
    "lenient":   (0.60, 10.0, 2.5, 0.70),
}

HEADER_RE = re.compile(
    r">design:(?P<design>\d+)\s+n:(?P<n>\d+)\|"
    r"mpnn:(?P<mpnn>[\d.]+)\|"
    r"plddt:(?P<plddt>[\d.]+)\|"
    r"i_ptm:(?P<iptm>[\d.]+)\|"
    r"i_pae:(?P<ipae>[\d.]+)\|"
    r"rmsd:(?P<rmsd>[\d.]+)"
)


def _normalise_header(name):
    key = re.sub(r"[^a-z0-9]+", "_", name.strip().lower()).strip("_")
    return {
        "id": "design", "name": "design", "design_id": "design",
        "sequence": "pep", "peptide": "pep", "inhibitor": "pep",
        "inhibitor_sequence": "pep", "seq": "pep",
        "i_ptm": "iptm", "iptm": "iptm", "i_pae": "ipae",
        "ipae": "ipae", "interface_pae": "ipae",
        "interface_rmsd": "rmsd", "rmsd": "rmsd", "plddt": "plddt",
        "mpnn": "mpnn"
    }.get(key, key)


def parse_table(path):
    """Parse a TXT/TSV/CSV table containing inhibitor metrics."""
    with open(path, newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
            reader = csv.DictReader(fh, dialect=dialect)
        except csv.Error:
            fh.seek(0)
            reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            return []
        fmap = {_normalise_header(x): x for x in reader.fieldnames}
        required = {"iptm", "ipae", "rmsd", "plddt"}
        missing = required - set(fmap)
        if missing:
            raise ValueError("missing required metric columns: " + ", ".join(sorted(missing)))
        entries = []
        for rownum, row in enumerate(reader, 2):
            try:
                raw_design = row.get(fmap.get("design", ""), "") or rownum - 1
                entries.append({
                    "design": int(float(raw_design)),
                    "n": rownum - 1,
                    "mpnn": float(row.get(fmap.get("mpnn", ""), 0) or 0),
                    "plddt": float(row[fmap["plddt"]]),
                    "iptm": float(row[fmap["iptm"]]),
                    "ipae": float(row[fmap["ipae"]]),
                    "rmsd": float(row[fmap["rmsd"]]),
                    "pep": (row.get(fmap.get("pep", ""), "") or "").strip(),
                })
            except (ValueError, TypeError) as exc:
                print(f"WARNING: skipping row {rownum}: {exc}", file=sys.stderr)
        return entries


def parse_input(path):
    """Parse FASTA-style input or a tabular inhibitor-metrics file."""
    with open(path) as fh:
        first = next((ln.strip() for ln in fh if ln.strip()), "")
    return parse_fasta(path) if first.startswith(">") else (parse_table(path), None)


def parse_fasta(path):
    """Return a list of design dicts parsed from the metric-carrying FASTA."""
    entries, target = [], None
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]

    if len(lines) % 2 != 0:
        print(f"WARNING: odd line count ({len(lines)}); last record may be "
              "truncated.", file=sys.stderr)

    for i in range(0, len(lines) - 1, 2):
        m = HEADER_RE.match(lines[i])
        if not m:
            print(f"WARNING: skipping unparseable header: {lines[i][:60]}",
                  file=sys.stderr)
            continue
        seq = lines[i + 1]
        if "/" in seq:
            tgt, pep = seq.split("/", 1)
            target = target or tgt
        else:
            pep = seq  # single-chain output; no target given
        g = m.groupdict()
        entries.append({
            "design": int(g["design"]), "n": int(g["n"]),
            "mpnn": float(g["mpnn"]), "plddt": float(g["plddt"]),
            "iptm": float(g["iptm"]), "ipae": float(g["ipae"]),
            "rmsd": float(g["rmsd"]), "pep": pep,
        })
    return entries, target


def deduplicate(entries):
    """Collapse identical (design, peptide) pairs (MPNN sometimes repeats)."""
    seen, unique, dropped = set(), [], 0
    for e in entries:
        key = (e["design"], e["pep"])
        if key in seen:
            dropped += 1
            continue
        seen.add(key)
        unique.append(e)
    return unique, dropped


def passes(e, iptm_min, ipae_max, rmsd_max, plddt_min):
    return (e["iptm"] >= iptm_min and e["ipae"] <= ipae_max
            and e["rmsd"] <= rmsd_max and e["plddt"] >= plddt_min)


def add_composite_scores(entries, weights=None):
    """Z-score inhibitor metrics and combine them into one composite score."""
    weights = weights or {k: 1.0 for k in METRIC_DIRECTION}
    stats = {}
    for k in METRIC_DIRECTION:
        vals = [e[k] for e in entries]
        mu = statistics.mean(vals)
        sd = statistics.pstdev(vals) or 1.0
        stats[k] = (mu, sd)
    for e in entries:
        e["composite"] = sum(
            direction * weights.get(k, 1.0) * (e[k] - stats[k][0]) / stats[k][1]
            for k, direction in METRIC_DIRECTION.items()
        )


def anchor_flags(pep, aromatic="FWY", hydrophobic="LIMFV"):
    """Flag the p53-style i, i+4, i+7 helical anchor register and Trp presence.

    p53 engages MDM2 via Phe19 / Trp23 / Leu26 — one helical face, spacing
    i, i+4, i+7. Generalizes to any helix-in-cleft binder; tune the residue
    classes for a different target.
    """
    has_trp = "W" in pep
    has_register = any(
        pep[i] in (aromatic + hydrophobic)
        and pep[i + 4] in aromatic
        and pep[i + 7] in hydrophobic
        for i in range(len(pep) - 7)
    )
    return has_trp, has_register


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Filter and rank binder designs from a metric-carrying FASTA.")
    ap.add_argument("fasta", nargs="?", help="Input inhibitor metrics TXT/CSV/TSV or FASTA.")
    ap.add_argument("--select-file", action="store_true", help="Open a file picker to select the inhibitor metrics file.")
    ap.add_argument("-o", "--output", default="ranked_candidates.csv",
                    help="Output CSV path (default: ranked_candidates.csv).")
    ap.add_argument("--tier", choices=TIERS, default="standard",
                    help="Preset filter tier (default: standard). "
                         "Any explicit threshold flag overrides the preset.")
    ap.add_argument("--iptm", type=float, help="Min interface pTM.")
    ap.add_argument("--ipae", type=float, help="Max interface PAE.")
    ap.add_argument("--rmsd", type=float, help="Max design-vs-AF RMSD (Å).")
    ap.add_argument("--plddt", type=float, help="Min pLDDT.")
    ap.add_argument("--top", type=int, default=0,
                    help="Keep only the top N after ranking (0 = keep all).")
    ap.add_argument("--no-dedup", action="store_true",
                    help="Do not collapse duplicate (design, sequence) pairs.")
    ap.add_argument("--anchor-motif", action="store_true",
                    help="Score the p53-style i,i+4,i+7 anchor register.")
    ap.add_argument("--weights", default="1,1,1,1",
                    help="Composite weights: iptm,ipae,rmsd,plddt.")
    args = ap.parse_args(argv)

    if args.select_file:
        try:
            import tkinter as tk
            from tkinter import filedialog
            root = tk.Tk(); root.withdraw()
            selected = filedialog.askopenfilename(
                title="Select inhibitor metrics file",
                filetypes=[("Text/CSV files", "*.txt *.tsv *.csv"), ("All files", "*.*")]
            )
            root.destroy()
        except Exception as exc:
            sys.exit(f"Could not open file picker: {exc}")
        if not selected:
            sys.exit("No input file selected.")
        args.fasta = selected

    if not args.fasta:
        ap.error("provide an input file or use --select-file")

    try:
        weight_values = [float(x.strip()) for x in args.weights.split(",")]
        if len(weight_values) != 4 or any(x < 0 for x in weight_values):
            raise ValueError
    except ValueError:
        ap.error("--weights must contain four non-negative numbers: iptm,ipae,rmsd,plddt")
    weights = dict(zip(METRIC_DIRECTION, weight_values))

    # Resolve thresholds: start from tier preset, apply any explicit overrides.
    iptm_min, ipae_max, rmsd_max, plddt_min = TIERS[args.tier]
    if args.iptm is not None:  iptm_min = args.iptm
    if args.ipae is not None:  ipae_max = args.ipae
    if args.rmsd is not None:  rmsd_max = args.rmsd
    if args.plddt is not None: plddt_min = args.plddt

    try:
        entries, target = parse_input(args.fasta)
    except (OSError, ValueError) as exc:
        sys.exit(f"Could not parse input file: {exc}")
    if not entries:
        sys.exit("No parseable designs found — check the FASTA format.")
    total = len(entries)

    dropped = 0
    if not args.no_dedup:
        entries, dropped = deduplicate(entries)

    # Backbone self-consistency: how many surviving seqs per backbone pass.
    per_backbone_total = defaultdict(int)
    per_backbone_pass = defaultdict(int)
    for e in entries:
        per_backbone_total[e["design"]] += 1
        if passes(e, iptm_min, ipae_max, rmsd_max, plddt_min):
            per_backbone_pass[e["design"]] += 1

    kept = [e for e in entries
            if passes(e, iptm_min, ipae_max, rmsd_max, plddt_min)]
    if not kept:
        sys.exit("No designs passed the filter — try a looser tier.")

    add_composite_scores(kept, weights=weights)
    kept.sort(key=lambda e: -e["composite"])
    if args.top > 0:
        kept = kept[:args.top]

    # Write CSV.
    fields = ["rank", "design", "n", "peptide", "length",
              "iptm", "i_pae", "rmsd", "plddt", "mpnn", "composite",
              "backbone_pass_count"]
    if args.anchor_motif:
        fields += ["has_Trp", "anchor_register"]
    with open(args.output, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for rank, e in enumerate(kept, 1):
            row = [rank, e["design"], e["n"], e["pep"], len(e["pep"]),
                   f"{e['iptm']:.3f}", f"{e['ipae']:.2f}", f"{e['rmsd']:.2f}",
                   f"{e['plddt']:.3f}", f"{e['mpnn']:.2f}", f"{e['composite']:.3f}",
                   per_backbone_pass[e["design"]]]
            if args.anchor_motif:
                trp, reg = anchor_flags(e["pep"])
                row += ["Y" if trp else "N", "Y" if reg else "N"]
            w.writerow(row)

    # Console summary.
    uniq_bb = len({e["design"] for e in kept})
    robust = sorted(d for d, c in per_backbone_pass.items()
                    if c == per_backbone_total[d] and c > 1)
    print(f"Parsed {total} designs"
          + (f" ({dropped} duplicates removed)" if dropped else ""))
    print(f"Target length: {len(target) if target else 'n/a'} | "
          f"Tier: {args.tier} "
          f"(iptm>={iptm_min}, ipae<={ipae_max}, rmsd<={rmsd_max}, "
          f"plddt>={plddt_min})")
    print(f"Composite weights (iptm, ipae, rmsd, plddt): "
          f"{weights['iptm']}, {weights['ipae']}, {weights['rmsd']}, {weights['plddt']}")
    print(f"Passed filter: {len(kept)} sequences across {uniq_bb} backbones")
    print(f"Fully-consistent backbones (all surviving MPNN seqs pass): "
          f"{len(robust)}")
    if robust:
        print("  " + ", ".join(map(str, robust)))
    print(f"Wrote ranked shortlist → {args.output}")


if __name__ == "__main__":
    main()
