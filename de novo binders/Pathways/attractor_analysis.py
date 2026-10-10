#!/usr/bin/env python3
"""
attractor_analysis.py -- exact and sampled attractor analysis for a Boolean
rule file, with the method documented in the output.

Two independent procedures are run per condition:

  (1) VALUE PERCOLATION (exact, deterministic)
      Clamp the perturbed inputs, then repeatedly fix any node whose function
      has become constant over every assignment of its still-free inputs.
      Iterate to a fixpoint. The result is the smallest TRAP SPACE containing
      the clamped subspace: a subcube closed under the update, which no
      trajectory can leave. If percolation fixes every node, that single state
      is the unique attractor and it is reachable from EVERY initial state --
      a proof, not a sample.

  (2) EXHAUSTIVE TERMINAL-SCC SEARCH (exact, if the residual space is small)
      Over the nodes percolation left free, build the full asynchronous state
      transition graph and find its terminal strongly connected components
      (Tarjan). Terminal SCCs are exactly the attractors. This distinguishes
      fixed points from cyclic/complex attractors.

  (3) RANDOM-START SAMPLING (inexact, for comparison)
      Random walks from uniformly random initial states, reporting which
      attractor each converges to. Reported only to show agreement with the
      exact result; it proves nothing on its own.

Usage:
  python3 attractor_analysis.py --rules rules.txt \
      --clamp control:KRAS=1,MDM2=1,MDM4=1 \
      --clamp MDM2off:KRAS=1,MDM2=0,MDM4=1 \
      [--samples 120] [--out report.md]
"""

import argparse
import itertools
import random
import re
import sys
from collections import defaultdict

IDENT = re.compile(r"\b(?!and\b|or\b|not\b|True\b|False\b)[A-Za-z_][A-Za-z0-9_]*\b")


def parse_rules(path):
    rules = {}
    for ln in open(path, encoding="utf-8"):
        ln = ln.split("#")[0].strip()
        if not ln:
            continue
        t, e = ln.split("=", 1)
        rules[t.strip()] = e.strip()
    return rules


def inputs_of(expr):
    return sorted(set(IDENT.findall(expr)))


def percolate(rules, code, clamp, max_free=18):
    """Smallest trap space containing the clamped subcube. Exact."""
    fixed = dict(clamp)
    rounds = 0
    while True:
        rounds += 1
        changed = False
        for n, expr in rules.items():
            if n in fixed:
                continue
            ins = inputs_of(expr)
            free = [v for v in ins if v not in fixed]
            if len(free) > max_free:
                continue
            vals = set()
            for combo in itertools.product([False, True], repeat=len(free)):
                env = dict(fixed)
                env.update(dict(zip(free, combo)))
                env.setdefault(n, False)
                for v in ins:
                    env.setdefault(v, False)
                vals.add(bool(eval(code[n], {"__builtins__": {}}, env)))
                if len(vals) > 1:
                    break
            if len(vals) == 1:
                fixed[n] = vals.pop()
                changed = True
        if not changed:
            return fixed, rounds


def terminal_sccs(rules, code, fixed, free):
    """Exhaustive Tarjan over the asynchronous STG of the residual subspace."""
    k = len(free)
    idx = {n: i for i, n in enumerate(free)}

    def succ(state):
        env = dict(fixed)
        env.update({n: bool(state >> idx[n] & 1) for n in free})
        out = set()
        for n in free:
            v = bool(eval(code[n], {"__builtins__": {}}, env))
            if v != env[n]:
                out.add(state ^ (1 << idx[n]))
        return out

    N = 1 << k
    index = {}
    low = {}
    onstack = {}
    stack = []
    comps = []
    for root in range(N):                      # iterative Tarjan
        if root in index:
            continue
        work = [(root, iter(succ(root)))]
        index[root] = low[root] = len(index)
        stack.append(root); onstack[root] = True
        while work:
            v, it = work[-1]
            advanced = False
            for w in it:
                if w not in index:
                    index[w] = low[w] = len(index)
                    stack.append(w); onstack[w] = True
                    work.append((w, iter(succ(w))))
                    advanced = True
                    break
                elif onstack.get(w):
                    low[v] = min(low[v], index[w])
            if advanced:
                continue
            work.pop()
            if work:
                low[work[-1][0]] = min(low[work[-1][0]], low[v])
            if low[v] == index[v]:
                comp = []
                while True:
                    w = stack.pop(); onstack[w] = False; comp.append(w)
                    if w == v:
                        break
                comps.append(comp)
    term = []
    for comp in comps:
        s = set(comp)
        if all(succ(x) <= s for x in comp):
            term.append(sorted(comp))
    return term


def sample_walks(rules, code, clamp, nodes, n_samples, steps, seed=0):
    rng = random.Random(seed)
    free = [n for n in nodes if n not in clamp]
    hits = defaultdict(int)
    for _ in range(n_samples):
        st = {n: rng.choice([False, True]) for n in nodes}
        st.update(clamp)
        for _ in range(steps):
            n = rng.choice(free)
            st[n] = bool(eval(code[n], {"__builtins__": {}}, st))
        hits[tuple(st[n] for n in nodes)] += 1
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rules", required=True)
    ap.add_argument("--clamp", action="append", required=True,
                    help="NAME:NODE=0/1,NODE=0/1")
    ap.add_argument("--samples", type=int, default=120)
    ap.add_argument("--steps", type=int, default=6000)
    ap.add_argument("--max-exhaustive", type=int, default=22,
                    help="max residual free nodes for the exact SCC search")
    ap.add_argument("--out")
    args = ap.parse_args()

    rules = parse_rules(args.rules)
    code = {n: compile(e, "<r>", "eval") for n, e in rules.items()}
    nodes = sorted(rules)
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    say(f"# Attractor analysis\n")
    say(f"Rules: `{args.rules}`  ({len(nodes)} nodes, "
        f"state space 2^{len(nodes)} = {2**len(nodes):,})\n")

    for spec in args.clamp:
        name, assigns = spec.split(":", 1)
        clamp = {}
        for a in assigns.split(","):
            k, v = a.split("=")
            clamp[k.strip()] = bool(int(v))
        say(f"\n## Condition: {name}")
        say(f"Clamped: {', '.join(f'{k}={int(v)}' for k, v in clamp.items())}  "
            f"-> reduced network of {len(nodes)-len(clamp)} free nodes "
            f"({2**(len(nodes)-len(clamp)):,} states)\n")

        fixed, rounds = percolate(rules, code, clamp)
        free = [n for n in nodes if n not in fixed]
        say(f"**Percolation** converged in {rounds} rounds; "
            f"fixed {len(fixed)}/{len(nodes)} nodes, {len(free)} left free.")
        if not free:
            on = sorted(n for n, v in fixed.items() if v)
            say(f"\nEvery node is determined, so the clamped subspace percolates "
                f"to a single state. That state is a fixed point, it is the "
                f"ONLY attractor, and it is reachable from every initial "
                f"state. This is exact.\n")
            say(f"ON: {', '.join(on) if on else '(none)'}")
            say(f"OFF: {', '.join(n for n in nodes if not fixed[n])}")
        else:
            say(f"Free after percolation: {', '.join(free)}\n")
            if len(free) <= args.max_exhaustive:
                term = terminal_sccs(rules, code, fixed, free)
                say(f"**Exhaustive terminal-SCC search** over "
                    f"{2**len(free):,} residual states: {len(term)} "
                    f"attractor(s).")
                for i, comp in enumerate(term, 1):
                    kind = "fixed point" if len(comp) == 1 else \
                           f"cyclic/complex attractor ({len(comp)} states)"
                    say(f"\n- attractor {i}: {kind}")
                    if len(comp) == 1:
                        st = dict(fixed)
                        st.update({n: bool(comp[0] >> j & 1)
                                   for j, n in enumerate(free)})
                        say(f"  ON: {', '.join(sorted(n for n in nodes if st[n]))}")
            else:
                say(f"Residual space too large for exhaustive search "
                    f"(> 2^{args.max_exhaustive}).")

        hits = sample_walks(rules, code, clamp, nodes, args.samples,
                            args.steps)
        say(f"\n**Random-start sampling** ({args.samples} uniformly random "
            f"initial states, {args.steps} asynchronous updates each): "
            f"{len(hits)} distinct end state(s).")
        if len(hits) == 1 and not free:
            end = dict(zip(nodes, list(hits)[0]))
            agree = all(end[n] == fixed[n] for n in nodes)
            say(f"Sampling agrees with the exact result: "
                f"{'yes' if agree else 'NO -- investigate'}.")
        say(f"\nNote: sampling covers {args.samples} of "
            f"{2**(len(nodes)-len(clamp)):,} states "
            f"({args.samples/2**(len(nodes)-len(clamp)):.2e} of the space). "
            f"It is corroboration, not proof; the percolation result above is "
            f"the exact one.")

    if args.out:
        open(args.out, "w", encoding="utf-8").write("\n".join(L) + "\n")
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
