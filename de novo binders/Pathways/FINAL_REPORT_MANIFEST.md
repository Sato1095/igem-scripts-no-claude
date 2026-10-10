# Final report — file manifest

Project: Boolean model of MDM2 and KRAS degradation in NSCLC (A549).
All checksums are SHA-256 of the files as generated.

---

## TIER 1 — Main figures and tables

| file | role | sha256 (first 12) |
|---|---|---|
| `network_v2_edges.csv` | **Table 1.** The model network: 29 nodes, 99 signed directed edges, each with sign, sign basis, mechanism, resources, PMIDs, curation effort and provenance. | 38fd2a69f60a |
| `network_v2_nodes.csv` | **Table 2.** Node list with role (input / intermediate / readout) and in/out degree. | 6f743d1572dd |
| `rules_v3_unfitted.txt` | **Table 3.** Boolean rules. Every departure from the default convention is marked OVERRIDE with its reason inline. | 2d96124f1e8d |
| `attractor_report.md` | **Table 4 / Figure 1.** Exact attractor per condition by value percolation. All 29 nodes determined in all four conditions; unique attractor provably reachable from all 2^26 states. | b50a4a60e7ab |
| `screen_report.md` | **Table 5 / Figure 2.** Rule-variant screen: full 32-row grid, which rule forms are determined, failure modes, predictions of passing variants, discriminating experiment. | d73b8a2d402e |

## TIER 2 — Supplementary methods and data

| file | role | sha256 (first 12) |
|---|---|---|
| `model_scope_v1.txt` | **Supp. Table S1.** Declared node scope from literature, with module comments. The stated boundary of the model. | de1bc39496de |
| `decisions_final.csv` | **Supp. Table S2.** 25 manual edge adjudications: 5 sign calls, 20 drops, each with rationale and reference basis. | 6f60eab683d9 |
| `screen_criteria_v2.md` | **Supp. Methods.** Pre-registered pass criteria, fixed before any variant was run, including the documented v1→v2 changes (C3 withdrawal, XIAP removal, R4 redefinition). | cc829977c479 |
| `rules_justification.md` | **Supp. Methods.** Three-tier justification of every rule: mechanical, structural override, or fitted. Contains the explicit statement of which rules were data-influenced. | 1b47438ea861 |
| `validation_both_arms_A549.csv` | **Supp. Table S3.** Both proteomics arms mapped onto model nodes, with log2 fold change, p-values and directional calls. | 2548abe823c9 |
| `screen_results.json` | **Supp. Data.** Machine-readable full screen output: every variant, every condition, full attractor, failed criteria. | 975e5d480d13 |
| `network_v2_changelog.csv` | **Supp. Table S4.** All 28 changes from the 117-edge network to network_v2, each with original sign, action and reason. | 1cdfd5266397 |

## TIER 3 — Code (repository / Zenodo)

| file | role | sha256 (first 12) |
|---|---|---|
| `build_network.py` | **Step 1.** Network construction from OmniPath with provenance tracking, evidence thresholds and a run manifest. Produces `final/edges.csv` (31 nodes / 117 edges). | 23a636da85ba |
| `make_network_v2.py` | **Step 2.** Applies the documented literature corrections to produce network_v2 (29 nodes / 99 edges) plus `network_v2_changelog.csv`. Verified to reproduce `network_v2_edges.csv` exactly. | 509203151fde |
| `attractor_analysis.py` | **Step 3.** Exact attractor analysis: value percolation + exhaustive terminal-SCC search + random-start sampling for corroboration. | 576a1ece9a13 |
| `run_screen.py` | **Step 4.** The 32-variant rule screen executing the pre-registered criteria; also regenerates `screen_report.md`. | 1de0d3e810be |
| `make_rules.py` | **Step 3.** Rule drafting, sign-consistency validation, and MaBoSS/BoolNet export. | 695bdcdfd308 |

### Pipeline, end to end

Run `run_pipeline.ps1` (Windows) or `run_pipeline.sh` (POSIX) to execute all
of it in one command. The individual steps are:

```
build_network.py   --anchors KRAS MDM2 --sinks CASP3 CASP7 PARP1 CDKN1A \
                   --nodes model_scope_v1.txt --min-resources 2 \
                   --min-references 1 --min-curation-effort 5 \
                   --terminal-sinks --decisions decisions_final.csv -d final
make_network_v2.py --in final/edges.csv --out .
make_rules.py      --edges network_v2_edges.csv --nodes network_v2_nodes.csv \
                   --fixed KRAS MDM2 MDM4 --out rules_draft.txt
                   # hand-edited to rules_v3_unfitted.txt; overrides documented
                   # inline and in rules_justification.md
attractor_analysis.py --rules rules_v3_unfitted.txt --clamp ... --out attractor_report.md
run_screen.py      --rules rules_v3_unfitted.txt \
                   --out-json screen_results.json --out-md screen_report.md
```

Only the rules step involves manual work; every other artifact regenerates
from its inputs.

**Also required in the repository, not generated here:**
- `run01/raw/` — the archived OmniPath snapshot (raw API responses with URL, timestamp and checksum per file). Without this the network cannot be regenerated, because OmniPath updates continuously.
- `final/manifest.json` — the run manifest for the 31-node build.

## NOT for the report

| file | why |
|---|---|
| `rules_v2.txt` | The hand-fitted rule set. Superseded by the screen, which does the same job by a stated procedure. Keep privately as a record. |
| `screen_criteria_v1.md` | Superseded by v2. Keep only if you want to show the criteria revision history. |
| `rules_default_OR_OR.txt`, `rules_exceptions.md` | Intermediate working files. |
| `model_v2.bnd/.bnet/.cfg`, `cond_*.cfg` | Exports of the fitted rule set; regenerate from `rules_v3_unfitted.txt` plus the chosen variant before use. |
| `validation_nutlin_A549.csv` | Subsumed by `validation_both_arms_A549.csv`. |

---

## Claims the files support

1. A provenance-tracked signed network for KRAS/MDM2 signalling and intrinsic
   apoptosis, with all manual curation documented.
2. Exact attractor analysis: one fixed point per condition, no cyclic
   attractors, outcome independent of initial state.
3. MDM2 degradation alone never produces apoptosis — 0 of 32 rule variants.
4. Combined degradation produces apoptosis in 28 of 32 variants.
5. True synergy (combination kills, neither single agent does) in 4 of 32,
   requiring both BCL-XL as a coincidence detector and p53-dependent PUMA.
6. The discriminating experiment: PUMA induction after KRAS degradation in
   p53-intact A549, or FOXO3 knockdown in the same background.

## Claims the files do NOT support

- That the model is validated. Consistency with one drug arm, using rules
  selected after inspecting model–data disagreements, is not validation.
- That KRAS degradation alone is ineffective. The dataset used to suggest this
  did not significantly reduce KRAS (0.91x, p=0.49) and was withdrawn from the
  criteria for that reason.
- Anything about the combination in cells. No dataset in this project perturbs
  both targets.

## Held back / outstanding

- Bortezomib and seven other drugs in the A549 proteomics remain unused — the
  only genuine out-of-sample test available.
- The in-silico single-node perturbation screen (29 nodes x 2 states x 4
  conditions) has not been run.
- Nutlin inhibits MDM2–p53 binding; the model clamps MDM2 to 0. These are not
  the same perturbation, and the clamped model cannot reproduce the observed
  rise in MDM2 protein.
