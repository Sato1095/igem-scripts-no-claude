#!/usr/bin/env python3
"""
run_screen.py -- enumerate rule variants and score them against the
pre-registered criteria in screen_criteria_v2.md.

Method: for each of the 2^5 = 32 variants, and each of the 4 conditions, the
attractor is computed by EXACT value percolation (clamp the perturbed inputs,
then repeatedly fix any node whose function has become constant over every
assignment of its still-free inputs; iterate to a fixpoint). Where percolation
determines every node, the resulting state is provably the unique attractor,
reachable from every initial state. No sampling is involved.

Criteria are read from this file's CRITERIA dict, which transcribes
screen_criteria_v2.md sections C1 and C2. C3 was withdrawn before running.
The both_degraded condition contributes no criteria; it is the prediction.

Every variant is reported, pass or fail.
"""

import argparse
import itertools
import json
import os
import re
import sys

IDENT = re.compile(r"\b(?!and\b|or\b|not\b|True\b|False\b)[A-Za-z_][A-Za-z0-9_]*\b")


# ---- the five varied rules (form A = default convention, form B = alternative)
VARIANTS = {
 "R1": {"A": {"BCL2L1": "not (BAD or BCL2L11 or BBC3 or PMAIP1 or TP53)"},
        "B": {"BCL2L1": "not (BAD and (BBC3 or PMAIP1))"}},
 "R2": {"A": {"BBC3": "TP53 or FOXO3"},
        "B": {"BBC3": "TP53"}},
 "R3": {"A": {"CDKN1A": "(TP53 or FOXO3) and not (GSK3B or MDM2 or MAPK1)"},
        "B": {"CDKN1A": "(TP53 or FOXO3) and not MDM2"}},
 "R4": {"A": {"BAX":  "(TP53 or BCL2L11 or GSK3B or BBC3) and not (BCL2L1 or MCL1 or AKT1)",
              "BAK1": "not (BCL2L1 or MCL1)"},
        "B": {"BAX":  "(TP53 or BCL2L11 or GSK3B or BBC3) and not ((BCL2L1 and MCL1) or AKT1)",
              "BAK1": "not (BCL2L1 and MCL1)"}},
 "R5": {"A": {"TP53": "not (MDM2 or MDM4)"},
        "B": {"TP53": "not (MDM2 and MDM4)"}},
}

CONDITIONS = {
 "control":        {"KRAS": True,  "MDM2": True,  "MDM4": True},
 "MDM2_degraded":  {"KRAS": True,  "MDM2": False, "MDM4": True},
 "KRAS_degraded":  {"KRAS": False, "MDM2": True,  "MDM4": True},
 "both_degraded":  {"KRAS": False, "MDM2": False, "MDM4": True},
}

# ---- pre-registered criteria (screen_criteria_v2.md, C1 and C2) -------------
CRITERIA = {
 "control": {"MAPK1": 1, "MAPK3": 1, "AKT1": 1, "TP53": 0,
             "MCL1": 1, "BCL2L1": 1, "CASP3": 0, "CASP7": 0},
 "MDM2_degraded": {"TP53": 1, "CDKN1A": 1, "BBC3": 1, "CASP3": 0, "CASP7": 0,
                   "MAPK1": 1, "MAPK3": 1, "BCL2L1": 1},
}
REPORT_NODES = ["MAPK1", "AKT1", "TP53", "CDKN1A", "BBC3", "PMAIP1", "BAD",
                "BCL2L11", "MCL1", "BCL2L1", "BAX", "BAK1", "CYCS", "XIAP",
                "CASP9", "CASP3", "CASP7"]


def parse(path):
    r = {}
    for ln in open(path, encoding="utf-8"):
        ln = ln.split("#")[0].strip()
        if ln:
            t, e = ln.split("=", 1)
            r[t.strip()] = e.strip()
    return r


def percolate(rules, clamp, max_free=18):
    code = {n: compile(e, "<r>", "eval") for n, e in rules.items()}
    fixed = dict(clamp)
    while True:
        changed = False
        for n, expr in rules.items():
            if n in fixed:
                continue
            ins = IDENT.findall(expr)
            free = sorted({v for v in ins if v not in fixed})
            if len(free) > max_free:
                continue
            vals = set()
            for combo in itertools.product([False, True], repeat=len(free)):
                env = dict(fixed); env.update(dict(zip(free, combo)))
                for v in set(ins) | {n}:
                    env.setdefault(v, False)
                vals.add(bool(eval(code[n], {"__builtins__": {}}, env)))
                if len(vals) > 1:
                    break
            if len(vals) == 1:
                fixed[n] = vals.pop(); changed = True
        if not changed:
            return fixed


def write_report(rows, path, rules_path):
    """Regenerate screen_report.md. Previously done by an ad-hoc script, which
    made the report un-reproducible; it is part of the tool now."""
    P = [r for r in rows if r["pass"]]
    nodes = ["TP53", "CDKN1A", "BCL2L1", "MCL1", "BAX", "CYCS", "CASP3", "CASP7"]
    L = ["# Rule-variant screen — results\n",
         "Criteria: `screen_criteria_v2.md` (fixed before running).",
         "Method: exact value percolation, 32 variants x 4 conditions. No sampling.",
         f"Base rules: `{os.path.basename(rules_path)}`. Script: `run_screen.py`.\n",
         f"**{len(P)} of {len(rows)} variants satisfy all criteria.**\n",
         "## 1. Full grid\n", "| variant | pass | failed criteria |", "|---|---|---|"]
    for r in rows:
        L.append(f"| `{r['variant']}` | {'**PASS**' if r['pass'] else 'fail'} | "
                 f"{'; '.join(r['failures']) if r['failures'] else '—'} |")
    L.append("\n## 2. Phenotype across all variants\n")
    L.append("D = apoptosis (CASP3 = CASP7 = 1)\n")
    L.append("| variant | pass | control | MDM2 | KRAS | BOTH |")
    L.append("|---|---|---|---|---|---|")
    def death(a): return a["CASP3"] == 1 and a["CASP7"] == 1
    for r in rows:
        A = r["attractors"]
        cells = " | ".join("D" if death(A[c]) else "."
                           for c in ["control", "MDM2_degraded",
                                     "KRAS_degraded", "both_degraded"])
        L.append(f"| `{r['variant']}` | {'yes' if r['pass'] else ''} | {cells} |")
    n_both = sum(1 for r in rows if death(r["attractors"]["both_degraded"]))
    n_mdm2 = sum(1 for r in rows if death(r["attractors"]["MDM2_degraded"]))
    n_kras = sum(1 for r in rows if death(r["attractors"]["KRAS_degraded"]))
    n_syn = sum(1 for r in rows if death(r["attractors"]["both_degraded"])
                and not death(r["attractors"]["MDM2_degraded"])
                and not death(r["attractors"]["KRAS_degraded"]))
    L.append(f"\nApoptosis reached in: control {sum(1 for r in rows if death(r['attractors']['control']))}/{len(rows)}, "
             f"MDM2 alone {n_mdm2}/{len(rows)}, KRAS alone {n_kras}/{len(rows)}, "
             f"combination {n_both}/{len(rows)}.")
    L.append(f"\n**True synergy** (combination kills, neither single agent does): "
             f"{n_syn}/{len(rows)}.")
    L.append("\n## 3. Attractors of the passing variants\n")
    for r in P:
        L.append(f"### `{r['variant']}`\n")
        L.append("| condition | " + " | ".join(nodes) + " |")
        L.append("|---|" + "---|" * len(nodes))
        for c, a in r["attractors"].items():
            L.append(f"| {c} | " + " | ".join(str(a[n]) for n in nodes) + " |")
        L.append("")
    open(path, "w", encoding="utf-8").write("\n".join(L) + "\n")
    print(f"  wrote {path}")


def main():
    ap = argparse.ArgumentParser(description="32-variant rule screen")
    ap.add_argument("--rules", required=True, help="base rules file")
    ap.add_argument("--out-json", default="screen_results.json")
    ap.add_argument("--out-md", default="screen_report.md")
    args = ap.parse_args()
    base = parse(args.rules)
    ids = list(VARIANTS)
    rows = []
    for combo in itertools.product("AB", repeat=len(ids)):
        label = "".join(f"{i}{c}" for i, c in zip(ids, combo))
        rules = dict(base)
        for i, c in zip(ids, combo):
            rules.update(VARIANTS[i][c])
        rec = {"variant": label, "forms": dict(zip(ids, combo)),
               "attractors": {}, "failures": [], "undetermined": []}
        for cond, clamp in CONDITIONS.items():
            fx = percolate(rules, clamp)
            free = [n for n in rules if n not in fx]
            if free:
                rec["undetermined"].append(f"{cond}: {','.join(free)}")
            rec["attractors"][cond] = {n: int(fx.get(n, -1)) for n in REPORT_NODES}
            for node, want in CRITERIA.get(cond, {}).items():
                got = fx.get(node)
                if got is None or int(got) != want:
                    rec["failures"].append(f"{cond}.{node}="
                                           f"{'?' if got is None else int(got)}"
                                           f"(want {want})")
        rec["pass"] = not rec["failures"] and not rec["undetermined"]
        rows.append(rec)
    json.dump(rows, open(args.out_json, "w"), indent=1)
    write_report(rows, args.out_md, args.rules)

    npass = sum(r["pass"] for r in rows)
    print(f"{len(rows)} variants, {npass} pass\n")
    hdr = f"{'variant':<22}{'pass':<6}{'n_fail':<8}first failures"
    print(hdr); print("-" * 100)
    for r in rows:
        print(f"{r['variant']:<22}{'YES' if r['pass'] else 'no':<6}"
              f"{len(r['failures']):<8}{'; '.join(r['failures'][:3])}")
    return rows


if __name__ == "__main__":
    main()
