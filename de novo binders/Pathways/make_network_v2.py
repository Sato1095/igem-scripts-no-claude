#!/usr/bin/env python3
"""
make_network_v2.py -- derive network_v2 from the OmniPath-built network by
applying the documented literature corrections found during rule construction.

Input : final/edges.csv           (output of build_network.py, 31 nodes/117 edges)
Output: network_v2_edges.csv      (29 nodes / 99 edges)
        network_v2_nodes.csv
        network_v2_changelog.csv  (every change, with its reason)

Every correction below was argued from biology or from database artefact, NOT
from the proteomics. This script is the record of that step; before it existed
the transformation was done ad hoc and Table 1 was not regenerable.

  python3 make_network_v2.py --in final/edges.csv --out .
"""

import argparse
import csv
import os
from collections import defaultdict

# ---------------------------------------------------------------------------
# Corrections. Each entry: (source, target) -> (action, reason)
# action is "sign:<new>" or "drop"
# ---------------------------------------------------------------------------
CORRECTIONS = {
 ("AKT1", "FOXO3"): ("sign:inhibits",
   "AKT phosphorylates FOXO3 to exclude it from the nucleus; recorded as "
   "activation, which would invert the entire PI3K arm"),
 ("MAPK3", "MAP2K1"): ("sign:inhibits",
   "ERK-to-MEK feedback is inhibitory (T292); mirror of the MAPK1 edge already "
   "corrected in decisions_final.csv"),
 ("MAPK1", "MAPK3"): ("drop",
   "ERK1 and ERK2 do not regulate one another; mirror of MAPK3->MAPK1, "
   "already dropped"),
 ("RAF1", "MAPK1"): ("drop", "RAF does not phosphorylate ERK; bypasses MEK"),
 ("RAF1", "MAPK3"): ("drop", "RAF does not phosphorylate ERK; bypasses MEK"),
 ("TP53", "KRAS"): ("drop", "KRAS is a clamped perturbation input"),
 ("AKT1", "MDM4"): ("drop", "MDM4 is a clamped constant input"),
 ("MDM2", "MDM4"): ("drop", "MDM4 is a clamped constant input"),
 ("APAF1", "CYCS"): ("drop",
   "reverse record; cytochrome c is released from mitochondria, not produced "
   "by APAF1. Retaining it makes the apoptosome self-sustaining"),
 ("BAX", "BAK1"): ("drop",
   "BAX and BAK are parallel effectors; neither activates the other"),
 ("TP53", "CASP3"): ("drop",
   "transcriptional bypass: lets p53 fire the executioner caspase without "
   "mitochondrial permeabilisation, collapsing the decision being modelled"),
 ("TP53", "CASP7"): ("drop", "transcriptional bypass of the intrinsic route"),
 ("TP53", "CASP9"): ("drop", "transcriptional bypass of apoptosome assembly"),
 ("FOXO3", "CASP3"): ("drop", "transcriptional bypass of the intrinsic route"),
 ("APAF1", "CASP3"): ("drop", "indirect; the route runs through CASP9"),
}

DROP_NODES = {
 "PARP1": "PARP1's only remaining in-edges are transcriptional (TP53, MAPK1); "
          "--terminal-sinks removed CASP3->PARP1, so PARP1 cannot function as "
          "a cleavage readout. Apoptosis is read from CASP3/CASP7 instead",
}

MERGE = {"PIK3CA": "PI3K", "PIK3R1": "PI3K"}
MERGE_REASON = ("p85 and p110 are recorded with mutual activation, which locks "
                "permanently ON once either fires; merged into a single PI3K "
                "node")

MANUAL_ADD = [{
 "source": "KRAS", "target": "PI3K", "sign": "activates",
 "sign_basis": "manual_decision", "mechanism": "binding to p110alpha RBD",
 "mechanism_source": "MANUAL", "n_resources": "1", "n_references": "1",
 "curation_effort": "2", "resources": "SIGNOR", "pmids": "",
 "provenance": "MANUAL_ADD",
 "rationale": "RAS binds the p110alpha RBD. Single-resource in OmniPath and "
              "therefore excluded by the evidence threshold, leaving the whole "
              "PI3K/AKT arm disconnected from the perturbed node. Restored "
              "deliberately; without it KRAS degradation cannot reach AKT, "
              "BAD, FOXO3 or GSK3B and the model cannot address its question",
}]

FIXED_INPUTS = {"KRAS": "input (perturbation)", "MDM2": "input (perturbation)",
                "MDM4": "input (constant)"}
READOUTS = {"CASP3": "readout (apoptosis)", "CASP7": "readout (apoptosis)",
            "CDKN1A": "readout (arrest)"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True,
                    help="final/edges.csv from build_network.py")
    ap.add_argument("--out", default=".", help="output directory")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.inp, newline="", encoding="utf-8-sig")))
    cols = list(rows[0].keys())
    for extra in ("mechanism_flag", "mechanism_ptm", "mechanism_ptm_resources"):
        if extra not in cols:
            cols.append(extra)

    kept, log = {}, []
    for r in rows:
        r = {c: r.get(c, "") for c in cols}
        s, t = r["source"].strip().upper(), r["target"].strip().upper()

        if s in DROP_NODES or t in DROP_NODES:
            node = s if s in DROP_NODES else t
            log.append([s, t, r["sign"], "node removed", DROP_NODES[node]])
            continue

        if (s, t) in CORRECTIONS:
            action, reason = CORRECTIONS[(s, t)]
            if action == "drop":
                log.append([s, t, r["sign"], "dropped", reason])
                continue
            new = action.split(":", 1)[1]
            log.append([s, t, r["sign"], f"sign -> {new}", reason])
            r["sign"] = new
            r["sign_basis"] = "manual_decision"
            r["provenance"] = "OmniPath+MANUAL_DECISION"
            r["rationale"] = reason

        s2, t2 = MERGE.get(s, s), MERGE.get(t, t)
        if s2 == t2:
            log.append([s, t, r["sign"], "removed (internal to merged node)",
                        MERGE_REASON])
            continue
        if (s2, t2) != (s, t):
            log.append([s, t, r["sign"], f"remapped -> {s2}->{t2}", MERGE_REASON])
        r["source"], r["target"] = s2, t2

        key = (s2, t2)
        if key in kept:
            a, b = int(kept[key]["curation_effort"] or 0), int(r["curation_effort"] or 0)
            if b > a:
                kept[key] = r
            log.append([s, t, r["sign"], f"merged duplicate of {s2}->{t2}",
                        "better-supported record retained"])
        else:
            kept[key] = r

    for add in MANUAL_ADD:
        row = {c: add.get(c, "") for c in cols}
        kept[(add["source"], add["target"])] = row
        log.append([add["source"], add["target"], add["sign"], "ADDED",
                    add["rationale"]])

    edges = sorted(kept.values(), key=lambda r: (r["source"], r["target"]))
    os.makedirs(args.out, exist_ok=True)

    p = os.path.join(args.out, "network_v2_edges.csv")
    with open(p, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader(); w.writerows(edges)

    A, I = defaultdict(int), defaultdict(int)
    for r in edges:
        (A if r["sign"] == "activates" else I)[r["target"]] += 1
    nodes = sorted({r["source"] for r in edges} | {r["target"] for r in edges})
    with open(os.path.join(args.out, "network_v2_nodes.csv"), "w",
              newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "role", "in_activators", "in_inhibitors", "out_degree"])
        for n in nodes:
            role = FIXED_INPUTS.get(n) or READOUTS.get(n) or "intermediate"
            w.writerow([n, role, A[n], I[n],
                        sum(1 for r in edges if r["source"] == n)])

    with open(os.path.join(args.out, "network_v2_changelog.csv"), "w",
              newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["source", "target", "original_sign", "action", "reason"])
        w.writerows(log)

    print(f"  in : {len(rows)} edges")
    print(f"  out: {len(edges)} edges, {len(nodes)} nodes")
    print(f"  changelog: {len(log)} entries -> network_v2_changelog.csv")


if __name__ == "__main__":
    main()
