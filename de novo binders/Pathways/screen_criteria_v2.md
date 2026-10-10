# Pre-registered criteria for the rule-variant screen -- v2

**Version 2.** Changes from v1, all made BEFORE any variant was run:
(a) C3 (KRAS-degraded criteria) demoted from pass criteria to a consistency
note, because the AdPROM dataset did not significantly reduce KRAS;
(b) the XIAP criterion removed from C2 as unsatisfiable by construction;
(c) R4 redefined to vary the guard logic at BAX and BAK1 together, since BAK1
reaches CYCS independently and varying BAX alone would have had no effect.

**Status: fixed before any variant is run. Nothing below may be
changed after seeing results. If a criterion turns out to be wrong, it is
amended in a new version of this file with the reason and the date, and the
screen is re-run in full.**

Model: network_v2 (29 nodes, 99 edges). Variants differ only in the five rules
listed in section 3.

---

## 1. What this screen is and is not

It is **model selection**: of the rule combinations that are defensible on
biological grounds, which are consistent with the observed single-agent
responses in A549?

It is **not validation**. The criteria are built from the same proteomics that
motivated inspecting these rules in the first place, so consistency here is a
constraint, not evidence of correctness.

Two acknowledged contaminations, to be stated in any write-up:

- The five varied rules were chosen after seeing which rules produced
  model-data disagreements. The screen therefore tests *which alternatives are
  consistent*, not whether the model was correct a priori.
- `TP53` (R5) was previously set to the AND form with the Nutlin data in view.
  It is returned to the screen as a free variant so the choice is made by the
  stated procedure rather than by hand.

## 2. Data used, and data held back

**Used as pass criteria:** the Nutlin arm only
(validation_both_arms_A549.csv).

**Used as a consistency note, not a criterion:** the AdPROM RAS-degradation
arm — see C3 below for why.

**Held back for out-of-sample testing:** Bortezomib and the remaining seven
drugs. These must not be consulted until the screen is complete and a variant
set has been chosen.

**Not available at all:** the combination arm. No dataset here perturbs both
targets, so the combination prediction cannot be tested by any means currently
in hand.

## 3. Varied rules (2 forms each -> 32 variants)

| id | node | form A | form B |
|----|------|--------|--------|
| R1 | BCL2L1 | `not (BAD or BCL2L11 or BBC3 or PMAIP1 or TP53)` | `not (BAD and (BBC3 or PMAIP1))` |
| R2 | BBC3 | `TP53 or FOXO3` | `TP53` |
| R3 | CDKN1A | `(TP53 or FOXO3) and not (GSK3B or MDM2 or MAPK1)` | `(TP53 or FOXO3) and not MDM2` |
| R4 | BAX **and** BAK1 guards | BAX: `... and not (BCL2L1 or MCL1 or AKT1)`; BAK1: `not (BCL2L1 or MCL1)` — both guards must be absent | BAX: `... and not ((BCL2L1 and MCL1) or AKT1)`; BAK1: `not (BCL2L1 and MCL1)` — either guard alone suffices |
| R5 | TP53 | `not (MDM2 or MDM4)` | `not (MDM2 and MDM4)` |

Form A is in every case the default produced by the "one activator / one
inhibitor suffices" convention. Form B is the stated alternative.

Note on R1: BAD is the canonical high-affinity BCL-XL antagonist, so form A is
arguably the better-supported binding biology. If form B is required to satisfy
the criteria, that tension must be reported rather than resolved silently.

All other 24 rules are held fixed at their rules_v3 values.

## 4. Conditions

| condition | KRAS | MDM2 | MDM4 |
|---|---|---|---|
| control | 1 | 1 | 1 |
| MDM2_degraded | 1 | 0 | 1 |
| KRAS_degraded | 0 | 1 | 1 |
| both_degraded | 0 | 0 | 1 |

`both_degraded` is run and reported but contributes NO pass criteria. It is the
prediction.

## 5. Pass criteria

All criteria are directional (node ON or OFF). No thresholds are applied to
fold changes; a Boolean node cannot be compared to a magnitude. Observed values
are given as mean log2 fold change for reference only.

### C1 — control: a viable proliferating A549 cell
Structural requirement, not from data. A rule set that cannot produce a
growing cell is rejected before any perturbation is considered.

| node | required | reason |
|---|---|---|
| MAPK1, MAPK3 | 1 | constitutive KRAS G12S signalling |
| AKT1 | 1 | PI3K arm active |
| TP53 | 0 | held off by MDM2/MDM4 |
| MCL1, BCL2L1 | 1 | guards intact |
| CASP3, CASP7 | 0 | no basal apoptosis |

### C2 — MDM2 degraded: p53 response without death
From the Nutlin arm.

| node | required | observed log2FC | call |
|---|---|---|---|
| TP53 | 1 | +2.24 | UP |
| CDKN1A | 1 | +2.13 | UP |
| BBC3 | 1 | +0.76 | UP |
| CASP3 | 0 | +0.04 | unchanged |
| CASP7 | 0 | +0.08 | unchanged |
| MAPK1, MAPK3 | 1 | -0.11, -0.29 | unchanged (MAPK arm not engaged) |
| BCL2L1 | 1 | +0.78 | UP — BCL-XL is NOT lost under MDM2 inhibition |

The BCL2L1 criterion is expected to be the most discriminating, since it
directly constrains R1. It is included because the observation is unambiguous
(+0.78 log2, consistent replicates) and because a guard that survives MDM2
inhibition is the only way the model can distinguish single agent from
combination.

### C3 — WITHDRAWN as a pass criterion

The AdPROM dataset did **not** significantly reduce KRAS (0.91x, p=0.49),
while HRAS fell to 0.46x (p=0.014). The observed absence of downstream change
is therefore consistent with two incompatible readings: that KRAS degradation
does nothing in A549, or that KRAS was never meaningfully degraded. These
cannot be distinguished in this data, and RAS isoforms share most tryptic
peptides, making isoform-level quantification unreliable in either direction.

A pass criterion built on an ambiguous null would force the model to be
INSENSITIVE to KRAS loss, potentially encoding an experimental artifact as a
structural feature — the opposite of what a model built to predict KRAS
degradation should do.

**Reclassified as a consistency note.** Because KRAS was not significantly
reduced, these cells were effectively still in the control condition, and the
proteome did not move. This is weak corroboration that the control attractor
is stable. It is not a test of the KRAS arm.

**Consequence:** the KRAS_degraded and both_degraded outcomes are now both
PURE PREDICTIONS, untested by any data in hand.

## 6. Scoring

Each variant is scored by exact value percolation under each condition. A
variant PASSES if it satisfies every criterion in C1, C2 and C3.

Reported for every one of the 32 variants, pass or fail:
- the full attractor per condition
- which individual criteria failed
- the predicted `both_degraded` outcome

The complete grid is published. Reporting only the passing variants is not
acceptable.

## 7. Interpretation rules, fixed in advance

- **Exactly one variant passes** — report it as the sole rule set consistent
  with the constraints, and state that consistency is not validation.
- **Several pass** — report the set. State which criterion would discriminate
  them and what experiment would provide it. Do not choose among them on
  aesthetic grounds.
- **None pass** — do not relax criteria to obtain a winner. Report the failure
  and the candidate explanations: Boolean abstraction conflating protein
  abundance with activity; missing nodes (NF-kB arm, BCL2); MDM2 *inhibition*
  (nutlin) treated as equivalent to *degradation*; the HRAS/KRAS ambiguity in
  C3; single timepoint versus steady-state attractor.

## 8. Known limitations of the screen itself

- With C3 withdrawn, roughly 13 binary constraints remain and many are
  correlated (CASP3, CASP7 and CASP9 move together; TP53 drives CDKN1A and
  BBC3). Effective constraint count is nearer 4. The screen is weakly
  determining and some varied rules are expected to come out UNDETERMINED —
  for example R2 (BBC3), since FOXO3 is off in the MDM2-degraded condition and
  both forms then behave identically. Undetermined rules are reported as such,
  not resolved by preference.
- Nutlin inhibits MDM2-p53 binding; the model clamps MDM2 to 0. These are not
  the same perturbation. MDM2 protein rose 2.84 log2 under nutlin precisely
  because p53 transcribes it — a behaviour the clamped model cannot show.
- PMAIP1 (NOXA), BCL2L11 (BIM) and MDM4 were not detected in either dataset
  and are therefore unconstrained by any criterion.
- A Boolean attractor is a steady state; the proteomics is a single timepoint.
  Agreement in direction is the most that can be claimed.

---

**Approved by:** ____________________  **Date:** ____________
