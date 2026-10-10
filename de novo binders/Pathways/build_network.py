#!/usr/bin/env python3
"""
build_network.py -- build a signed, directed, provenance-tracked network for a
KRAS / MDM2 -> apoptosis Boolean model, entirely from OmniPath.

NO BIOLOGY IS HARD-CODED IN THIS FILE.

The only inputs that encode a scientific choice arrive from outside the code:

  --anchors        the perturbed proteins (your experiment)
  --phenotype-re   a regex matched against annotation records (your endpoint
                   definition, e.g. apoptosis AND cell-cycle arrest)
  --decisions      a reviewable CSV of human sign calls, each requiring a PMID,
                   a rationale and a name
  --expression     a cell-line expression table (e.g. DepMap for A549)

Pipeline
--------
  1. fetch      raw TSVs from the OmniPath web service, cached + checksummed
  2. qualify    directionality and sign resolved from consensus fields;
                evidence thresholds applied; every rejection given a reason
  3. mechanism  joined from the enzyme-substrate endpoint (modification+residue)
  4. expand     undirected k-hop ball around the anchors (hub-capped)
  5. phenotype  annotations queried for the candidates only; matches become sinks
  6. prune      keep nodes on directed anchor->sink paths within a cutoff
  7. context    optional expression filter for the cell line
  8. review     unsigned/conflicting on-path edges -> review_queue.csv
  9. report     nodes, edges, exclusions, and a manifest.json recording every
                URL, checksum, parameter and count for the run

Outputs
-------
  nodes.csv, edges.csv, excluded_nodes.csv, excluded_edges.csv,
  review_queue.csv, manifest.json, raw/ (untouched API responses + .meta.json)

Offline use: pass --interactions-tsv / --enzsub-tsv / --annotations-tsv to work
from files you exported yourself; the manifest records their checksums instead.
"""

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time
from collections import defaultdict, deque
from urllib.parse import urlencode
from urllib.request import urlopen

SCRIPT_VERSION = "2.0"
BASE = "https://omnipathdb.org"
TRUEISH = {"1", "true", "t", "yes", "y"}

csv.field_size_limit(10 ** 7)


# --------------------------------------------------------------------------
# fetching (raw responses are cached verbatim and checksummed)
# --------------------------------------------------------------------------

def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _looks_like_error(raw):
    """OmniPath returns HTTP 200 with a plain-text error banner on bad params."""
    head = raw[:400].decode("utf-8", "replace")
    first = head.splitlines()[0] if head.splitlines() else ""
    if "not entirely good" in head.lower() or head.lstrip().startswith("<"):
        return True
    return "\t" not in first          # a real TSV header is tab-separated


def fetch(endpoint, params, cache_dir, refresh=False, timeout=300):
    """GET an OmniPath endpoint, cache the raw bytes, write a .meta.json."""
    os.makedirs(cache_dir, exist_ok=True)
    qs = urlencode(sorted(params.items()))
    tag = hashlib.sha256(f"{endpoint}?{qs}".encode()).hexdigest()[:12]
    path = os.path.join(cache_dir, f"{endpoint}_{tag}.tsv")
    if os.path.exists(path) and not refresh:
        return path
    url = f"{BASE}/{endpoint}?{qs}"
    print(f"  GET {url}", flush=True)
    with urlopen(url, timeout=timeout) as resp:
        raw = resp.read()
    if _looks_like_error(raw):
        msg = raw[:1500].decode("utf-8", "replace")
        sys.exit(f"\nERROR: the OmniPath service rejected this request:\n"
                 f"  {url}\n\n--- server said ---\n{msg}\n-------------------\n"
                 f"Fix the offending parameter above and re-run. Nothing was "
                 f"cached, so no stale error file is left behind.")
    with open(path, "wb") as fh:
        fh.write(raw)
    with open(path + ".meta.json", "w") as fh:
        json.dump({"url": url, "endpoint": endpoint, "params": params,
                   "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                time.gmtime()),
                   "bytes": len(raw), "sha256": _sha256(path)}, fh, indent=2)
    return path


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --------------------------------------------------------------------------
# qualification: direction, sign, evidence
# --------------------------------------------------------------------------

def _b(v):
    return str(v).strip().lower() in TRUEISH


def _i(v):
    try:
        return int(float(str(v).strip()))
    except (TypeError, ValueError):
        return 0


def _pmids(ref):
    """OmniPath references look like 'SIGNOR:12345678;SPIKE:23456789'."""
    out = []
    for token in (ref or "").replace(",", ";").split(";"):
        p = token.strip().split(":")[-1].strip()
        if p.isdigit() and p not in out:
            out.append(p)
    return out


def _resources(src):
    out = []
    for token in (src or "").replace(",", ";").split(";"):
        # sub-resources arrive as 'SIGNOR_CollecTRI'; keep the primary name
        r = token.strip().split("_")[0]
        if r and r not in out:
            out.append(r)
    return out


def resolve_sign(row):
    """Sign from OmniPath's own fields. Consensus wins; conflicts are kept."""
    cs, ci = _b(row.get("consensus_stimulation")), _b(row.get("consensus_inhibition"))
    s, i = _b(row.get("is_stimulation")), _b(row.get("is_inhibition"))
    if cs and not ci:
        return "activates", "consensus"
    if ci and not cs:
        return "inhibits", "consensus"
    if s and not i:
        return "activates", "is_stimulation"
    if i and not s:
        return "inhibits", "is_inhibition"
    if s and i:
        return "AMBIGUOUS", "resource_conflict"
    return "UNSIGNED", "no_sign_reported"


def qualify(rows, min_resources, min_references, keep_undirected=False,
            min_effort=0):
    """Collapse to one record per (source, target); attach evidence + reasons."""
    kept, rejected = {}, []
    for r in rows:
        s = (r.get("source_genesymbol") or "").strip().upper()
        t = (r.get("target_genesymbol") or "").strip().upper()
        if not s or not t:
            continue
        if s == t:
            rejected.append((s, t, "self_loop"))
            continue
        if not _b(r.get("is_directed")) and not keep_undirected:
            rejected.append((s, t, "not_directed"))
            continue

        refs = _pmids(r.get("references"))
        res = _resources(r.get("sources"))
        n_ref = _i(r.get("n_references")) or len(refs)
        n_res = _i(r.get("n_resources")) or len(res)
        if n_res < min_resources:
            rejected.append((s, t, f"resources {n_res} < {min_resources}"))
            continue
        if n_ref < min_references:
            rejected.append((s, t, f"references {n_ref} < {min_references}"))
            continue
        ce = _i(r.get("curation_effort"))
        if min_effort and ce < min_effort:
            rejected.append((s, t, f"curation_effort {ce} < {min_effort}"))
            continue

        sign, basis = resolve_sign(r)
        rec = {"source": s, "target": t, "sign": sign, "sign_basis": basis,
               "type": (r.get("type") or "").strip(),
               "resources": res, "pmids": refs,
               "n_resources": n_res, "n_references": n_ref,
               "curation_effort": _i(r.get("curation_effort")),
               "provenance": "OmniPath"}

        key = (s, t)
        if key not in kept:
            kept[key] = rec
            continue
        # same pair reported by more than one dataset: union the evidence
        m = kept[key]
        m["resources"] = list(dict.fromkeys(m["resources"] + res))
        m["pmids"] = list(dict.fromkeys(m["pmids"] + refs))
        m["n_resources"] = max(m["n_resources"], n_res)
        m["n_references"] = max(m["n_references"], n_ref)
        m["curation_effort"] = max(m["curation_effort"], rec["curation_effort"])
        if m["type"] != rec["type"] and rec["type"]:
            m["type"] = ";".join(dict.fromkeys(
                [x for x in (m["type"], rec["type"]) if x]))
        if m["sign"] != sign:
            real = {m["sign"], sign} - {"UNSIGNED"}
            if len(real) == 1:
                m["sign"], m["sign_basis"] = real.pop(), "single_dataset"
            else:
                m["sign"], m["sign_basis"] = "AMBIGUOUS", "dataset_conflict"
    return kept, rejected


def attach_mechanism(edges, enzsub_rows, exclude_re=None):
    """Mechanism comes from OmniPath, never from guessing.

    Post-translational detail is joined from the enzyme-substrate endpoint and
    kept in its OWN column with its own resources, so a questionable site can
    be traced back. It never overwrites the transcriptional label: an edge that
    is transcriptional in the interactions table stays labelled that way, and
    carrying both is itself a flag worth reviewing.
    """
    mech = defaultdict(set)
    ptm_ref = defaultdict(set)
    ptm_res = defaultdict(set)
    for r in enzsub_rows or []:
        e = (r.get("enzyme_genesymbol") or "").strip().upper()
        s_ = (r.get("substrate_genesymbol") or "").strip().upper()
        if not e or not s_:
            continue
        res_list = _resources(r.get("sources"))
        if exclude_re:
            kept_res = [x for x in res_list if not exclude_re.search(x)]
            raw = (r.get("sources") or "")
            if exclude_re.search(raw) and not kept_res:
                continue          # record rests only on excluded resources
            res_list = kept_res or res_list
        mod = (r.get("modification") or "").strip()
        res_t = (r.get("residue_type") or "").strip()
        off = (r.get("residue_offset") or "").strip()
        label = mod + (f"@{res_t}{off}" if res_t and off else "")
        if label:
            mech[(e, s_)].add(label)
        ptm_ref[(e, s_)].update(_pmids(r.get("references")))
        ptm_res[(e, s_)].update(res_list)

    for key, m in edges.items():
        is_txn = "transcriptional" in m["type"]
        ptm = ";".join(sorted(mech.get(key, ())))
        m["mechanism_ptm"] = ptm
        m["mechanism_ptm_resources"] = ";".join(sorted(ptm_res.get(key, ())))
        if is_txn and ptm:
            m["mechanism"] = "transcriptional_regulation;+PTM"
            m["mechanism_source"] = "OmniPath:interactions.type+enzsub"
            m["mechanism_flag"] = "PTM_RECORD_ON_TRANSCRIPTIONAL_EDGE"
        elif is_txn:
            m["mechanism"] = "transcriptional_regulation"
            m["mechanism_source"] = "OmniPath:interactions.type"
            m["mechanism_flag"] = ""
        elif ptm:
            m["mechanism"] = ptm
            m["mechanism_source"] = "OmniPath:enzsub"
            m["mechanism_flag"] = ""
        else:
            m["mechanism"] = ""
            m["mechanism_source"] = "not_reported"
            m["mechanism_flag"] = ""
        extra = [x for x in sorted(ptm_ref.get(key, ())) if x not in m["pmids"]]
        m["pmids"] = m["pmids"] + extra


# --------------------------------------------------------------------------
# topology
# --------------------------------------------------------------------------

def build_adjacency(edges):
    fwd, rev, und = defaultdict(set), defaultdict(set), defaultdict(set)
    for (s, t) in edges:
        fwd[s].add(t)
        rev[t].add(s)
        und[s].add(t)
        und[t].add(s)
    return fwd, rev, und


def bfs_ball(adj, seeds, hops, hub_degree):
    """Undirected k-hop ball; traversal stops at promiscuous hubs.

    Also returns the frontier -- the nodes one step beyond the ball -- so the
    boundary of the search is written to excluded_nodes.csv instead of being
    invisible.
    """
    dist = {s: 0 for s in seeds}
    q = deque(dist)
    hubs, frontier = set(), set()
    while q:
        n = q.popleft()
        if dist[n] >= hops:
            frontier.update(m for m in adj.get(n, ()) if m not in dist)
            continue
        if hub_degree and n not in seeds and len(adj.get(n, ())) > hub_degree:
            hubs.add(n)
            frontier.update(m for m in adj.get(n, ()) if m not in dist)
            continue
        for m in adj.get(n, ()):
            if m not in dist:
                dist[m] = dist[n] + 1
                q.append(m)
    return dist, hubs, frontier - set(dist)


def bfs_dist(adj, seeds):
    dist = {s: 0 for s in seeds}
    q = deque(dist)
    while q:
        n = q.popleft()
        for m in adj.get(n, ()):
            if m not in dist:
                dist[m] = dist[n] + 1
                q.append(m)
    return dist


# --------------------------------------------------------------------------
# phenotype sinks, discovered from annotations of the candidate set
# --------------------------------------------------------------------------

def uniprot_to_symbol(interaction_rows):
    """Derive an accession -> symbol map from the interactions table itself,
    so no extra identifier query is needed."""
    m = {}
    for r in interaction_rows:
        for acc_col, sym_col in (("source", "source_genesymbol"),
                                 ("target", "target_genesymbol")):
            acc = (r.get(acc_col) or "").strip()
            sym = (r.get(sym_col) or "").strip().upper()
            if acc and sym:
                m[acc] = sym
    return m


def phenotype_sinks(candidates, regex, resources, cache_dir, refresh,
                    annotations_tsv=None, batch=100, dump=0, acc2sym=None,
                    field="value", list_only=False):
    """Query annotations for the CANDIDATES ONLY, then regex-match the records.

    The match runs over every field of each annotation row, so this does not
    depend on the exact schema of any one annotation resource. Complex
    identifiers (FOS_JUN, PIK3CA_PIK3R1 ...) are split into their members
    before querying, because the endpoint takes proteins, not complexes.
    """
    pat = re.compile(regex, re.I)

    # complexes arrive from the TF datasets as MEMBER1_MEMBER2; query members
    queryable, complex_members = set(), {}
    for g in candidates:
        if "_" in g:
            members = [m for m in g.split("_") if m]
            complex_members[g] = members
            queryable.update(members)
        else:
            queryable.add(g)

    rows = []
    if annotations_tsv:
        rows = read_tsv(annotations_tsv)
    else:
        genes = sorted(queryable)
        for i in range(0, len(genes), batch):
            chunk = genes[i:i + batch]
            # NB: the annotations endpoint does NOT accept `genesymbols`;
            # sending it makes the service reject the whole request.
            p = {"proteins": ",".join(chunk), "format": "tsv"}
            if resources:
                p["resources"] = resources
            path = fetch("annotations", p, cache_dir, refresh)
            got = read_tsv(path)
            if not got:
                print(f"    WARNING: batch {i // batch + 1} returned 0 rows "
                      f"({os.path.basename(path)})")
            rows.extend(got)

    print(f"    annotation rows fetched: {len(rows)}")
    if rows:
        print(f"    annotation columns: {', '.join(list(rows[0].keys())[:12])}")

    # map member -> owning complex so a complex inherits its members' labels
    owner = defaultdict(set)
    for cx, members in complex_members.items():
        for m in members:
            owner[m].add(cx)

    hits = defaultdict(set)
    matched_labels = defaultdict(int)
    term_genes = defaultdict(set)
    seen_values, in_scope = set(), 0
    for r in rows:
        g = (r.get("genesymbol") or "").strip().upper()
        if not g:
            acc = (r.get("uniprot") or r.get("record_id") or "").strip()
            g = (acc2sym or {}).get(acc, "").upper()
        targets = ({g} & candidates) | owner.get(g, set())
        if not targets:
            continue
        in_scope += 1
        label = (r.get("value") or r.get("label") or "").strip()
        if len(seen_values) < 4000 and label:
            seen_values.add(label)
        # match one named column when it exists (lets you anchor a regex
        # precisely); fall back to the whole row for unknown schemas
        subject = (r.get(field) or "").strip() if field in r else \
            " ".join(str(v) for v in r.values() if v)
        if pat.search(subject):
            src = (r.get("source") or r.get("resource") or "?").strip()
            matched_labels[f"{src}:{label}"] += 1
            term_genes[label].update(targets)
            for t in targets:
                hits[t].add(f"{src}:{label}"[:120])

    if list_only:
        print(f"\n  terms matching /{regex}/ among your {len(candidates)} "
              f"candidates, rarest first\n  (rare = specific = a better model "
              f"output; common = whole-pathway membership)\n")
        for lab, genes in sorted(term_genes.items(), key=lambda x: len(x[1])):
            gl = ", ".join(sorted(genes)[:14])
            more = f" (+{len(genes) - 14} more)" if len(genes) > 14 else ""
            print(f"  [{len(genes):>3}]  {lab}\n         {gl}{more}")
        print(f"\n  {len(term_genes)} terms shown. Pick the ones that name an "
              f"OUTCOME, then re-run with\n  --phenotype-re "
              f"'^(TERM_A|TERM_B)$'\n")
        sys.exit(0)

    print(f"    rows about candidate genes: {in_scope}")
    print(f"    genes matching /{regex}/ in field '{field}': {len(hits)}")
    if matched_labels:
        print("    annotation terms that produced sinks (term: gene rows):")
        for lab, n in sorted(matched_labels.items(),
                             key=lambda x: -x[1])[:20]:
            print(f"      {n:>6}  {lab}")
    frac = len(hits) / max(len(candidates), 1)
    if frac > 0.25:
        print(f"\n    WARNING: {len(hits)} of {len(candidates)} candidate "
              f"nodes ({frac:.0%}) were classed as phenotype sinks.\n"
              f"    A sink is a model OUTPUT. At this breadth the path prune "
              f"does nothing and the\n    network will be too dense to write "
              f"Boolean rules for. Narrow --phenotype-re to\n    specific "
              f"terms from the list above (anchor it, e.g. "
              f"'^GOBP_EXECUTION_PHASE_OF_APOPTOSIS$').\n")
    if not hits and seen_values:
        sample = sorted(seen_values)[:dump or 25]
        print("    no match -- sample of the annotation vocabulary actually "
              "returned:")
        for v in sample:
            print(f"      {v}")
    return hits


# --------------------------------------------------------------------------
# human decisions (the replacement for hard-coded bridges)
# --------------------------------------------------------------------------

REQUIRED_DECISION_COLS = ["source", "target", "sign", "pmid", "rationale"]


def apply_decisions(edges, path, allow_add):
    """Apply reviewed sign calls. Every row needs a PMID, rationale and name."""
    applied, refused = [], []
    if not path:
        return applied, refused
    for r in read_tsv(path) if path.endswith(".tsv") else _read_csv(path):
        sign_peek = (r.get("sign") or "").strip().lower()
        needed = [c for c in REQUIRED_DECISION_COLS
                  if not (c == "pmid" and sign_peek == "drop")]
        missing = [c for c in needed if not (r.get(c) or "").strip()]
        if missing:
            refused.append((r.get("source"), r.get("target"),
                            "missing " + ",".join(missing)))
            continue
        s, t = r["source"].strip().upper(), r["target"].strip().upper()
        sign = r["sign"].strip().lower()
        if sign not in ("activates", "inhibits", "drop"):
            refused.append((s, t, f"bad sign '{sign}'"))
            continue
        key = (s, t)
        status = (r.get("status") or "").strip().upper()
        if key not in edges:
            if not (allow_add and status == "ADD"):
                refused.append((s, t, "absent from OmniPath; needs status=ADD "
                                      "and --allow-manual-edges"))
                continue
            edges[key] = {"source": s, "target": t, "type": "",
                          "resources": [], "pmids": [], "n_resources": 0,
                          "n_references": 0, "curation_effort": 0,
                          "mechanism": (r.get("mechanism") or "").strip(),
                          "mechanism_source": "MANUAL", "provenance": "MANUAL_ADD"}
        m = edges[key]
        if sign == "drop":
            m["sign"], m["sign_basis"] = "DROP", "manual_decision"
        else:
            m["sign"] = sign
            m["sign_basis"] = "manual_decision"
        # 'from_edge' = the call rests on the references OmniPath already
        # carries for this interaction; nothing new is asserted
        pm = (r.get("pmid") or "").strip()
        if pm.lower() not in ("from_edge", ""):
            m["pmids"] = list(dict.fromkeys(m.get("pmids", []) + _pmids(pm)))
        elif pm.lower() == "from_edge" and not m.get("pmids"):
            refused.append((s, t, "pmid=from_edge but the edge has no "
                                  "references of its own"))
            continue
        m["rationale"] = r["rationale"].strip()
        if m["provenance"] != "MANUAL_ADD":
            m["provenance"] = "OmniPath+MANUAL_DECISION"
        applied.append((s, t, sign, m["provenance"]))
    return applied, refused


def _read_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def suggest_nodes(edges, scope, d_f, d_b, min_links, limit):
    """Rank nodes OUTSIDE the declared scope by how well they connect it.

    These are proposals for review, never auto-included: a node earns a row by
    touching several scope members and sitting on a short anchor->sink route.
    """
    links = defaultdict(set)
    effort = defaultdict(int)
    resources = defaultdict(set)
    degree = defaultdict(set)
    for (a, b), m in edges.items():
        degree[a].add(b)
        degree[b].add(a)
        for n, other in ((a, b), (b, a)):
            if n not in scope and other in scope:
                links[n].add(other)
                effort[n] = max(effort[n], m["curation_effort"])
                resources[n].update(m["resources"])
    rows = []
    for n, partners in links.items():
        if len(partners) < min_links:
            continue
        f, b = d_f.get(n), d_b.get(n)
        route = (f + b) if (f is not None and b is not None) else ""
        tot = max(len(degree[n]), 1)
        spec = len(partners) / tot
        # specificity demotes promiscuous hubs: a TF wired to 900 genes that
        # happens to touch 20 of yours is not evidence your scope is missing it
        score = len(partners) * spec
        rows.append((round(score, 2), len(partners), tot, round(spec, 3), n,
                     effort[n], len(resources[n]), route,
                     ";".join(sorted(partners)[:12])))
    rows.sort(key=lambda r: -r[0])
    return [[n, k, tot, spec, score, route, eff, nres, part]
            for score, k, tot, spec, n, eff, nres, route, part in rows[:limit]]


def main():
    ap = argparse.ArgumentParser(
        description="Build a provenance-tracked signed network from OmniPath")
    ap.add_argument("--anchors", nargs="+", required=True,
                    help="perturbed proteins, e.g. KRAS MDM2")
    ap.add_argument("-d", "--dir", default="out")
    ap.add_argument("--cache", default="out/raw")
    ap.add_argument("--refresh", action="store_true")

    ap.add_argument("--datasets", default="omnipath,kinaseextra,pathwayextra,"
                                          "collectri")
    ap.add_argument("--organism", default="9606")
    ap.add_argument("--min-resources", type=int, default=2)
    ap.add_argument("--min-references", type=int, default=1)
    ap.add_argument("--keep-undirected", action="store_true")
    ap.add_argument("--min-curation-effort", type=int, default=0,
                    help="drop interactions below this curation effort; the "
                         "thin tail of an aggregated interactome is mostly "
                         "single-report noise")

    ap.add_argument("--hops", type=int, default=2,
                    help="undirected radius around the anchors")
    ap.add_argument("--hub-degree", type=int, default=250,
                    help="stop traversal at nodes with more neighbours (0=off)")
    ap.add_argument("--cutoff", type=int, default=6,
                    help="max d(anchor->node) + d(node->sink)")

    ap.add_argument("--phenotype-re",
                    default=r"apopto|cell[_ ]cycle[_ ]arrest|p53",
                    help="regex defining the model's output phenotypes")
    ap.add_argument("--annotation-resources",
                    default="MSigDB,CancerSEA,SIGNOR")
    ap.add_argument("--exclude-ptm-resources",
                    default=r"ProtMapper|Sparser|RLIMS|REACH|TRIPS|"
                            r"BEL-Large-Corpus|Medscan|EVEX",
                    help="drop enzyme-substrate records that rest only on "
                         "these resources (default: text-mining readers, which "
                         "routinely swap enzyme and substrate). '' disables.")
    ap.add_argument("--terminal-sinks", action="store_true",
                    help="drop out-edges of sink nodes: a node declared as a "
                         "model READOUT should not regulate the network "
                         "(OmniPath carries its full biology regardless)")
    ap.add_argument("--nodes", metavar="FILE",
                    help="SCOPE MODE: a declared node list (one symbol per "
                         "line). Fetches every qualified edge among them and "
                         "ranks candidate additions instead of expanding.")
    ap.add_argument("--suggest", type=int, default=40,
                    help="scope mode: how many candidate additions to rank")
    ap.add_argument("--suggest-min-links", type=int, default=3,
                    help="scope mode: minimum edges to the declared set")
    ap.add_argument("--list-terms", metavar="REGEX",
                    help="explore mode: list every annotation term matching "
                         "REGEX with the candidate genes it covers, then exit")
    ap.add_argument("--phenotype-field", default="value",
                    help="annotation column the regex is matched against "
                         "(default: value; use '*' to match the whole row)")
    ap.add_argument("--dump-vocabulary", type=int, default=25,
                    help="on zero matches, print this many observed labels")
    ap.add_argument("--sinks", nargs="*", default=[],
                    help="override the discovered sinks (not recommended)")

    ap.add_argument("--expression", help="CSV with columns gene,value")
    ap.add_argument("--expression-min", type=float, default=1.0)

    ap.add_argument("--decisions", help="reviewed sign calls (CSV)")
    ap.add_argument("--allow-manual-edges", action="store_true",
                    help="permit decisions.csv rows with status=ADD")

    ap.add_argument("--interactions-tsv")
    ap.add_argument("--enzsub-tsv")
    ap.add_argument("--annotations-tsv")
    args = ap.parse_args()

    os.makedirs(args.dir, exist_ok=True)
    anchors = {a.upper() for a in args.anchors}
    manifest = {"script_version": SCRIPT_VERSION,
                "run_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "argv": sys.argv[1:], "parameters": vars(args), "inputs": [],
                "counts": {}}

    def note_input(path):
        manifest["inputs"].append({"path": os.path.abspath(path),
                                   "sha256": _sha256(path)})

    # 1. fetch ---------------------------------------------------------------
    fields = ["sources", "references", "curation_effort", "n_references",
              "n_resources", "type"]
    if args.interactions_tsv:
        ia_path = args.interactions_tsv
    else:
        ia_path = fetch("interactions",
                        {"datasets": args.datasets, "genesymbols": "yes",
                         "organisms": args.organism, "format": "tsv",
                         "fields": ",".join(fields)},
                        args.cache, args.refresh)
    note_input(ia_path)
    ia_rows = read_tsv(ia_path)

    if args.enzsub_tsv:
        es_path = args.enzsub_tsv
    else:
        es_path = fetch("enzsub",
                        {"genesymbols": "yes", "organisms": args.organism,
                         "format": "tsv",
                         "fields": "sources,references,curation_effort"},
                        args.cache, args.refresh)
    note_input(es_path)
    es_rows = read_tsv(es_path)

    # 2-3. qualify + mechanism ----------------------------------------------
    edges, rejected = qualify(ia_rows, args.min_resources,
                              args.min_references, args.keep_undirected,
                              args.min_curation_effort)
    excl = re.compile(args.exclude_ptm_resources, re.I) \
        if args.exclude_ptm_resources else None
    attach_mechanism(edges, es_rows, excl)
    manifest["counts"]["interactions_fetched"] = len(ia_rows)
    manifest["counts"]["edges_qualified"] = len(edges)
    manifest["counts"]["edges_rejected"] = len(rejected)

    missing_anchors = [a for a in anchors
                       if not any(a in k for k in edges)]
    if missing_anchors:
        sys.exit(f"ERROR: anchors absent from the qualified network: "
                 f"{', '.join(missing_anchors)} -- loosen the evidence "
                 f"thresholds or widen --datasets")

    # 4. expand --------------------------------------------------------------
    fwd, rev, und = build_adjacency(edges)
    in_network_all = {g for k in edges for g in k}

    # ---- SCOPE MODE: annotate a declared node set, suggest additions --------
    if args.nodes:
        note_input(args.nodes)
        scope = {l.strip().upper() for l in open(args.nodes, encoding="utf-8")
                 if l.strip() and not l.startswith("#")}
        absent = sorted(scope - in_network_all)
        scope &= in_network_all
        sinks_d = {g.upper() for g in args.sinks} & scope
        final = {k: v for k, v in edges.items()
                 if k[0] in scope and k[1] in scope}
        if args.terminal_sinks and sinks_d:
            cut = sorted(f"{a}->{b}" for (a, b) in final if a in sinks_d)
            final = {k: v for k, v in final.items() if k[0] not in sinks_d}
            manifest["terminal_sink_edges_removed"] = cut
            print(f"  terminal sinks: removed {len(cut)} out-edge(s) from "
                  f"readout node(s) {', '.join(sorted(sinks_d))}")
        d_f = bfs_dist(fwd, anchors & in_network_all)
        d_b = bfs_dist(rev, sinks_d or (anchors & in_network_all))

        applied, refused = apply_decisions(final, args.decisions,
                                           args.allow_manual_edges)
        final = {k: v for k, v in final.items() if v["sign"] != "DROP"}
        review = {k: v for k, v in final.items()
                  if v["sign"] in ("AMBIGUOUS", "UNSIGNED")}
        sugg = suggest_nodes(edges, scope, d_f, d_b, args.suggest_min_links,
                             args.suggest)

        os.makedirs(args.dir, exist_ok=True)

        def w2(name, header, rows_):
            pth = os.path.join(args.dir, name)
            with open(pth, "w", newline="", encoding="utf-8") as fh:
                wr = csv.writer(fh)
                wr.writerow(header)
                wr.writerows(rows_)
            print(f"  {pth}  ({len(rows_)} rows)")

        hdr = ["source", "target", "sign", "sign_basis", "mechanism",
               "mechanism_source", "mechanism_flag", "mechanism_ptm",
               "mechanism_ptm_resources", "n_resources", "n_references",
               "curation_effort", "resources", "pmids", "type", "provenance",
               "rationale"]

        def erow(m):
            return [m["source"], m["target"], m["sign"], m["sign_basis"],
                    m.get("mechanism", ""), m.get("mechanism_source", ""),
                    m.get("mechanism_flag", ""), m.get("mechanism_ptm", ""),
                    m.get("mechanism_ptm_resources", ""),
                    m["n_resources"], m["n_references"], m["curation_effort"],
                    ";".join(m["resources"]), ";".join(m["pmids"]),
                    m.get("type", ""), m["provenance"], m.get("rationale", "")]

        w2("nodes.csv", ["id", "role", "d_from_anchor", "d_to_sink",
                         "in_degree", "out_degree"],
           [[n, "anchor" if n in anchors else
             ("sink" if n in sinks_d else "intermediate"),
             d_f.get(n, ""), d_b.get(n, ""),
             sum(1 for k in final if k[1] == n),
             sum(1 for k in final if k[0] == n)] for n in sorted(scope)])
        w2("edges.csv", hdr, [erow(m) for _, m in sorted(final.items())])
        w2("review_queue.csv", hdr + ["decision", "pmid", "rationale_in"],
           [erow(m) + ["", "", ""] for _, m in sorted(review.items())])
        w2("suggestions.csv", ["gene", "links_to_model", "total_degree",
                               "specificity", "score", "route_length",
                               "max_curation_effort", "n_resources",
                               "connects_to"], sugg)
        w2("scope_not_in_omnipath.csv", ["gene", "reason"],
           [[g, "no qualified edge at the current evidence thresholds"]
            for g in absent])

        manifest["mode"] = "scope"
        manifest["counts"].update(
            {"scope_nodes": len(scope), "scope_absent": len(absent),
             "edges_final": len(final), "edges_needing_review": len(review),
             "suggestions": len(sugg),
             "manual_decisions_applied": len(applied),
             "manual_decisions_refused": len(refused)})
        manifest["manual_decisions_refused"] = refused
        with open(os.path.join(args.dir, "manifest.json"), "w") as fh:
            json.dump(manifest, fh, indent=2)
        dens = len(final) / max(len(scope), 1)
        print(f"\n  scope: {len(scope)} nodes, {len(final)} edges "
              f"({dens:.1f} edges/node), {len(review)} needing review")
        if absent:
            print(f"  {len(absent)} scope gene(s) had no qualified edge: "
                  f"{', '.join(absent[:10])}")
        print(f"  {len(sugg)} candidate additions ranked in suggestions.csv")
        return

    # When sinks are declared, seed the search from BOTH ends: the anchor->sink
    # route is then covered at half the radius, keeping the candidate set small.
    declared = {g.upper() for g in args.sinks}
    in_network = {g for k in edges for g in k}
    missing_sinks = sorted(declared - in_network)
    if missing_sinks:
        print(f"  WARNING: declared sink(s) absent from the qualified network "
              f"and ignored: {', '.join(missing_sinks)}\n"
              f"           (lower --min-resources / --min-references, or check "
              f"the symbol)")
    declared &= in_network
    seeds = anchors | declared
    ball, hubs, frontier = bfs_ball(und, seeds, args.hops, args.hub_degree)
    candidates = set(ball)
    manifest["declared_sinks_used"] = sorted(declared)
    manifest["declared_sinks_missing"] = missing_sinks
    manifest["counts"]["candidates_in_ball"] = len(candidates)
    manifest["hubs_not_traversed"] = sorted(hubs)
    manifest["counts"]["frontier_nodes"] = len(frontier)

    # 5. phenotype sinks -----------------------------------------------------
    if declared:
        sinks = set(declared)
        sink_evidence = {s: {"manual:--sinks"} for s in sinks}
    else:
        sink_evidence = phenotype_sinks(candidates,
                                        args.list_terms or args.phenotype_re,
                                        args.annotation_resources, args.cache,
                                        args.refresh, args.annotations_tsv,
                                        dump=args.dump_vocabulary,
                                        acc2sym=uniprot_to_symbol(ia_rows),
                                        field=args.phenotype_field,
                                        list_only=bool(args.list_terms))
        sinks = set(sink_evidence)
    sinks -= anchors
    if not sinks:
        sys.exit("ERROR: no phenotype sinks found -- widen --phenotype-re, "
                 "--annotation-resources or --hops")
    manifest["counts"]["sinks_discovered"] = len(sinks)

    # 6. prune ---------------------------------------------------------------
    sub = {k: v for k, v in edges.items()
           if k[0] in candidates and k[1] in candidates}
    sfwd, srev, _ = build_adjacency(sub)
    d_f = bfs_dist(sfwd, anchors)
    d_b = bfs_dist(srev, sinks)

    kept_nodes = {}
    dropped_nodes = {n: f"beyond --hops {args.hops} from an anchor"
                     for n in frontier}
    for n in candidates:
        f, b = d_f.get(n), d_b.get(n)
        if n in anchors:
            kept_nodes[n] = ("anchor", 0, b if b is not None else -1)
        elif n in sinks and f is not None:
            kept_nodes[n] = ("sink", f, 0)
        elif f is None and b is None:
            dropped_nodes[n] = "off any anchor->sink path (unreached both ways)"
        elif f is None:
            dropped_nodes[n] = "not reachable from an anchor"
        elif b is None:
            dropped_nodes[n] = "cannot reach a phenotype sink"
        elif f + b > args.cutoff:
            dropped_nodes[n] = f"path too long ({f}+{b}={f + b} > {args.cutoff})"
        else:
            kept_nodes[n] = ("intermediate", f, b)

    # 7. expression context --------------------------------------------------
    if args.expression:
        note_input(args.expression)
        expr = {}
        for r in _read_csv(args.expression):
            g = (r.get("gene") or "").strip().upper()
            try:
                expr[g] = float(r.get("value"))
            except (TypeError, ValueError):
                continue
        for n in list(kept_nodes):
            if n in anchors:
                continue
            if n in expr and expr[n] < args.expression_min:
                dropped_nodes[n] = (f"not expressed in the cell line "
                                    f"({expr[n]} < {args.expression_min})")
                del kept_nodes[n]
            elif n not in expr:
                kept_nodes[n] = kept_nodes[n] + ("no_expression_record",)

    keep = set(kept_nodes)
    final = {k: v for k, v in sub.items() if k[0] in keep and k[1] in keep}

    # 8. decisions + review queue -------------------------------------------
    applied, refused = apply_decisions(final, args.decisions,
                                       args.allow_manual_edges)
    final = {k: v for k, v in final.items() if v["sign"] != "DROP"}
    review = {k: v for k, v in final.items()
              if v["sign"] in ("AMBIGUOUS", "UNSIGNED")}
    manifest["counts"]["edges_final"] = len(final)
    manifest["counts"]["edges_needing_review"] = len(review)
    manifest["counts"]["manual_decisions_applied"] = len(applied)
    manifest["counts"]["manual_decisions_refused"] = len(refused)
    manifest["counts"]["manual_added_edges"] = sum(
        1 for v in final.values() if v["provenance"] == "MANUAL_ADD")
    manifest["manual_decisions_refused"] = refused

    # 9. report --------------------------------------------------------------
    def w(name, header, rows):
        p = os.path.join(args.dir, name)
        with open(p, "w", newline="", encoding="utf-8") as fh:
            wr = csv.writer(fh)
            wr.writerow(header)
            wr.writerows(rows)
        print(f"  {p}  ({len(rows)} rows)")

    w("nodes.csv", ["id", "role", "d_from_anchor", "d_to_sink",
                    "phenotype_evidence", "flags"],
      [[n, v[0], v[1], v[2], ";".join(sorted(sink_evidence.get(n, []))),
        ";".join(v[3:])] for n, v in sorted(kept_nodes.items())])

    edge_header = ["source", "target", "sign", "sign_basis", "mechanism",
                   "mechanism_source", "mechanism_flag", "mechanism_ptm",
                   "mechanism_ptm_resources", "n_resources", "n_references",
                   "curation_effort", "resources", "pmids", "type",
                   "provenance", "rationale"]

    def edge_row(m):
        return [m["source"], m["target"], m["sign"], m["sign_basis"],
                m.get("mechanism", ""), m.get("mechanism_source", ""),
                m.get("mechanism_flag", ""), m.get("mechanism_ptm", ""),
                m.get("mechanism_ptm_resources", ""),
                m["n_resources"], m["n_references"], m["curation_effort"],
                ";".join(m["resources"]), ";".join(m["pmids"]),
                m.get("type", ""), m["provenance"], m.get("rationale", "")]

    w("edges.csv", edge_header,
      [edge_row(m) for _, m in sorted(final.items())])
    w("review_queue.csv", edge_header + ["decision", "pmid", "rationale_in"],
      [edge_row(m) + ["", "", ""] for _, m in sorted(review.items())])
    w("excluded_nodes.csv", ["gene", "reason"],
      sorted(dropped_nodes.items()))
    w("excluded_edges.csv", ["source", "target", "reason"],
      sorted(set(rejected))[:100000])

    with open(os.path.join(args.dir, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"\n  nodes kept: {len(keep)}   edges: {len(final)}   "
          f"needing review: {len(review)}")
    if review:
        print("  -> resolve review_queue.csv into a decisions CSV "
              "(source,target,sign,pmid,rationale) and re-run "
              "with --decisions")


if __name__ == "__main__":
    main()
