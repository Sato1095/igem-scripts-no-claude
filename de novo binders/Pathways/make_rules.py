#!/usr/bin/env python3
"""
make_rules.py -- draft, validate and export the Boolean rules for a network
built by build_network.py.

The edge list cannot tell you whether two activators combine with AND or OR.
That is the modeller's call, so this tool drafts a starting point, marks the
nodes where the default is probably wrong, and then CHECKS your hand-edited
rules back against the network.

  1. draft     python3 make_rules.py --edges network_v2_edges.csv \\
                   --nodes network_v2_nodes.csv \\
                   --fixed KRAS MDM2 MDM4 --out rules_draft.txt
  2. edit      fix the TODO nodes by hand -> rules_v3_unfitted.txt
  3. check     python3 make_rules.py --check rules_v3_unfitted.txt \\
                   --edges network_v2_edges.csv --nodes network_v2_nodes.csv \\
                   --fixed KRAS MDM2 MDM4
  4. export    ... --check rules_v3_unfitted.txt --export-maboss model_v3.bnd
               ... --check rules_v3_unfitted.txt --export-boolnet model_v3.bnet

Nothing here is hard-coded to a particular network: --edges and --nodes are
required arguments. The example paths are just the current artifacts.

The check verifies, for every input of every rule, that the rule's actual
monotone behaviour matches the sign of the edge: an activator must be able to
turn the node on, an inhibitor must be able to turn it off. It also reports
edges your rules silently ignore, and relations your rules invent.

Rule syntax: TARGET = expression, using and / or / not and parentheses.
Lines starting with # are comments.
"""

import argparse
import csv
import itertools
import os
import random
import re
import sys
from collections import defaultdict

RULE_RE = re.compile(r"^\s*([A-Za-z0-9_]+)\s*=\s*(.+?)\s*$")
IDENT_RE = re.compile(r"\b(?!and\b|or\b|not\b|True\b|False\b)[A-Za-z_][A-Za-z0-9_]*\b")
# a mechanism naming a physical complex is the usual reason a rule needs AND
COMPLEX_HINT = re.compile(r"complex|binding|apoptosome|pore|dimer", re.I)


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def sanitize(name):
    return name.strip().upper().replace("-", "_")


def load_network(edges_path, nodes_path):
    edges = []
    for r in read_csv(edges_path):
        s, t = sanitize(r["source"]), sanitize(r["target"])
        sign = (r.get("sign") or "").strip().lower()
        if sign not in ("activates", "inhibits"):
            print(f"  skipping unsigned edge {s}->{t} (sign={sign!r})")
            continue
        edges.append({"source": s, "target": t, "sign": sign,
                      "mechanism": r.get("mechanism", ""),
                      "effort": int(float(r.get("curation_effort") or 0))})
    nodes = {}
    if nodes_path:
        for r in read_csv(nodes_path):
            nodes[sanitize(r["id"])] = r.get("role", "")
    for e in edges:
        nodes.setdefault(e["source"], "")
        nodes.setdefault(e["target"], "")
    return edges, nodes


def draft(edges, nodes, out_path, fixed):
    acts, inhs = defaultdict(list), defaultdict(list)
    for e in edges:
        (acts if e["sign"] == "activates" else inhs)[e["target"]].append(e)

    lines = [
        "# Boolean rules drafted by make_rules.py -- EDIT BEFORE USE.",
        "# The default for every node is  (any activator) AND NOT (any inhibitor).",
        "# That is wrong wherever inputs must co-occur (complex formation,",
        "# apoptosome assembly, cooperative pore formation). Nodes where that is",
        "# likely are marked TODO below.",
        "#",
        "# Syntax: TARGET = expression, with and / or / not and parentheses.",
        "",
    ]
    todo = []
    for n in sorted(nodes):
        a = sorted(acts[n], key=lambda e: -e["effort"])
        i = sorted(inhs[n], key=lambda e: -e["effort"])
        role = nodes.get(n, "")

        if n in fixed:
            lines.append(f"# {n}: fixed input (declared with --fixed)")
            lines.append(f"{n} = {n}")
            lines.append("")
            continue
        if not a and not i:
            lines.append(f"# {n}: no regulators in the network -> input node"
                         + (f"  [{role}]" if role else ""))
            lines.append(f"{n} = {n}")
            lines.append("")
            continue

        lines.append(f"# {n}"
                     + (f"  [{role}]" if role else "")
                     + f"  ({len(a)} activator(s), {len(i)} inhibitor(s))")
        for e in a + i:
            mark = "+" if e["sign"] == "activates" else "-"
            mech = e["mechanism"] or "mechanism not reported"
            lines.append(f"#   {mark} {e['source']:<10} effort={e['effort']:<5} {mech}")

        needs_and = (len(a) > 1 and any(COMPLEX_HINT.search(e["mechanism"] or "")
                                        for e in a))
        if needs_and or len(a) + len(i) >= 5:
            reason = ("activators include a complex/binding mechanism -- do they "
                      "have to co-occur?" if needs_and else
                      f"{len(a) + len(i)} inputs -- check the combination logic")
            lines.append(f"#   TODO: {reason}")
            todo.append(n)

        expr_a = " or ".join(e["source"] for e in a) if a else ""
        expr_i = " or ".join(e["source"] for e in i) if i else ""
        if a and i:
            la = f"({expr_a})" if len(a) > 1 else expr_a
            li = f"({expr_i})" if len(i) > 1 else expr_i
            rule = f"{la} and not {li}"
        elif a:
            rule = expr_a
        else:
            li = f"({expr_i})" if len(i) > 1 else expr_i
            rule = f"not {li}"
        lines.append(f"{n} = {rule}")
        lines.append("")

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  wrote {out_path}  ({len(nodes)} rules)")
    if todo:
        print(f"  {len(todo)} node(s) marked TODO: {', '.join(todo)}")
    return todo


def parse_rules(path):
    rules = {}
    order = []
    for ln, raw in enumerate(open(path, encoding="utf-8"), 1):
        line = raw.split("#")[0].strip()
        if not line:
            continue
        m = RULE_RE.match(line)
        if not m:
            sys.exit(f"ERROR: {path}:{ln}: cannot parse rule: {raw.strip()!r}")
        tgt, expr = sanitize(m.group(1)), m.group(2)
        if tgt in rules:
            sys.exit(f"ERROR: {path}:{ln}: duplicate rule for {tgt}")
        rules[tgt] = expr
        order.append(tgt)
    return rules, order


def inputs_of(expr):
    return sorted(set(IDENT_RE.findall(expr)))


def monotonicity(expr, inputs, trials=400):
    """Classify each input's effect by evaluating the rule both ways."""
    code = compile(expr, "<rule>", "eval")
    others = [x for x in inputs]
    result = {}
    for target in inputs:
        pos = neg = False
        rest = [x for x in others if x != target]
        if len(rest) <= 12:
            combos = itertools.product([False, True], repeat=len(rest))
        else:
            combos = ([random.choice([False, True]) for _ in rest]
                      for _ in range(trials))
        for combo in combos:
            env = dict(zip(rest, combo))
            env[target] = False
            lo = bool(eval(code, {"__builtins__": {}}, env))
            env[target] = True
            hi = bool(eval(code, {"__builtins__": {}}, env))
            if hi and not lo:
                pos = True
            if lo and not hi:
                neg = True
            if pos and neg:
                break
        result[target] = ("mixed" if pos and neg else
                          "activates" if pos else
                          "inhibits" if neg else "inert")
    return result


def check(rules, order, edges, nodes, fixed):
    problems, warnings = [], []
    edge_sign = {(e["source"], e["target"]): e["sign"] for e in edges}

    for n in sorted(nodes):
        if n not in rules:
            problems.append(f"no rule for node {n}")
    for tgt, expr in rules.items():
        if tgt not in nodes:
            warnings.append(f"rule for {tgt}, which is not a network node")
        for src in inputs_of(expr):
            if src not in nodes and src != tgt:
                problems.append(f"{tgt}: rule references unknown node {src}")

    used = set()
    for tgt, expr in rules.items():
        ins = inputs_of(expr)
        try:
            mono = monotonicity(expr, ins)
        except SyntaxError as exc:
            problems.append(f"{tgt}: cannot evaluate rule ({exc})")
            continue
        for src in ins:
            if src == tgt:
                continue          # self-reference: an input node holding state
            used.add((src, tgt))
            declared = edge_sign.get((src, tgt))
            if declared is None:
                problems.append(
                    f"{tgt}: rule uses {src}, but no such edge exists")
            elif mono[src] == "inert":
                warnings.append(
                    f"{tgt}: {src} appears in the rule but never changes the "
                    f"outcome (masked by other terms)")
            elif mono[src] == "mixed":
                warnings.append(
                    f"{tgt}: {src} acts as activator in some states and "
                    f"inhibitor in others (edge says {declared})")
            elif mono[src] != declared:
                problems.append(
                    f"{tgt}: {src} is {declared} in the network but behaves as "
                    f"{mono[src]} in the rule")

    for (s, t), sign in edge_sign.items():
        if (s, t) not in used and t not in fixed:
            warnings.append(f"edge {s} -{'>' if sign == 'activates' else '|'} "
                            f"{t} is not used by any rule")
    return problems, warnings


def export_maboss(rules, order, path):
    with open(path, "w", encoding="utf-8") as fh:
        for n in order:
            expr = rules[n].replace(" and ", " AND ").replace(" or ", " OR ")
            expr = re.sub(r"\bnot\b", "NOT", expr)
            fh.write(f"Node {n} {{\n  logic = {expr};\n"
                     f"  rate_up = @logic ? 1 : 0;\n"
                     f"  rate_down = @logic ? 0 : 1;\n}}\n\n")
    cfg = os.path.splitext(path)[0] + ".cfg"
    with open(cfg, "w", encoding="utf-8") as fh:
        for n in order:
            fh.write(f"{n}.istate = 0;\n")
        fh.write("\nthread_count = 4;\nmax_time = 50;\ntime_tick = 0.5;\n"
                 "sample_count = 5000;\n")
    print(f"  wrote {path} and {cfg}")


def export_boolnet(rules, order, path):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("targets, factors\n")
        for n in order:
            expr = rules[n].replace(" and ", " & ").replace(" or ", " | ")
            expr = re.sub(r"\bnot\b", "!", expr)
            fh.write(f"{n}, {expr}\n")
    print(f"  wrote {path}")


def main():
    ap = argparse.ArgumentParser(description="Draft and validate Boolean rules")
    ap.add_argument("--edges", required=True)
    ap.add_argument("--nodes")
    ap.add_argument("--out", help="draft mode: rules file to write")
    ap.add_argument("--check", help="check mode: rules file to validate")
    ap.add_argument("--fixed", nargs="*", default=[],
                    help="nodes held as constant inputs (drug targets, "
                         "context constants)")
    ap.add_argument("--export-maboss")
    ap.add_argument("--export-boolnet")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    random.seed(args.seed)

    fixed = {sanitize(g) for g in args.fixed}
    edges, nodes = load_network(args.edges, args.nodes)
    print(f"  network: {len(nodes)} nodes, {len(edges)} signed edges")

    if args.out:
        draft(edges, nodes, args.out, fixed)
        return

    if not args.check:
        sys.exit("give either --out (draft) or --check (validate)")

    rules, order = parse_rules(args.check)
    problems, warnings = check(rules, order, edges, nodes, fixed)

    for w in warnings:
        print(f"  WARN  {w}")
    for p in problems:
        print(f"  FAIL  {p}")
    print(f"\n  {len(rules)} rules, {len(problems)} problem(s), "
          f"{len(warnings)} warning(s)")

    if problems:
        sys.exit(1)
    if args.export_maboss:
        export_maboss(rules, order, args.export_maboss)
    if args.export_boolnet:
        export_boolnet(rules, order, args.export_boolnet)


if __name__ == "__main__":
    main()
