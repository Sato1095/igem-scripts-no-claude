# Justification of the rule set and initial conditions

Model: network_v2 (29 nodes, 99 edges), rules_v2.txt, A549.

## 1. The default convention

Every node takes the form `(any activator) AND NOT (any inhibitor)` unless
stated otherwise. This encodes two claims: one activator is sufficient to turn
the node on, and one inhibitor is sufficient to turn it off. Both were adopted
as the null assumption, not because they are always true, but because they are
the weakest assumption consistent with a signed edge list: an edge asserts that
the source can change the target, and nothing more. Departures from the default
are therefore the places where extra knowledge was added, and each is listed
below.

Twelve of 29 rules are the unmodified default.

## 2. Tier 1 — mechanical (no intervention)

AKT1, APAF1, BAD, BAX, BCL2L11, CASP3, CASP7, CYCS, DIABLO, MAPK1, MAPK3,
MCL1, MDM2, PTEN, XIAP.

These follow directly from the network under the default convention. Note in
particular:

- `BAD = not (AKT1 or RAF1)` — two kinases, either sufficient. This is the case
  worked out by hand as a truth table and it agrees with the mechanical result.
- `BAX = (TP53 or BCL2L11 or GSK3B or BBC3) and not (BCL2L1 or MCL1 or AKT1)` —
  the guards enter under OR inside the NOT, so BOTH must be absent for BAX to
  activate. This is the stringent form, and it is what the default produces; no
  intervention was needed.
- `MCL1 = (MAPK1 or MAPK3) and not (BCL2L11 or GSK3B or PMAIP1 or TP53 or BBC3)` —
  ERK stabilises, and any one BH3-only protein suffices to neutralise. MCL1 is
  therefore lost under either single agent, which is why it cannot be the node
  that discriminates them.

## 3. Tier 2 — structural and semantic overrides

These would be made regardless of the validation data. Each corrects a mismatch
between what an edge list can express and what the node is defined to mean.

| rule | override | justification |
|---|---|---|
| `CASP9 = CYCS and APAF1 and not (...)` | AND between activators | The apoptosome is a stoichiometric complex. Under OR, cytochrome c activates caspase-9 with no Apaf-1 present, which makes the execution module decorative. This is a composition constraint, not a fitted choice. |
| `FOXO3 = not (AKT1 or MAPK1)` | activator requirement removed | FOXO3 is constitutively expressed; regulation is nuclear exclusion by AKT. Its only activator edge is a weak p53 transcriptional record (effort 19), which under the default would make FOXO3 p53-dependent and prevent KRAS loss from inducing BIM. |
| `GSK3B = not (AKT1 or MAPK1 or MAPK3)` | as above | Constitutively active kinase held off by AKT (effort 128). Only activator edge is TP53 at effort 6. |
| `BAK1 = not (BCL2L1 or MCL1)` | as above | BAK is constitutively resident at the outer membrane; control is guard sequestration. |
| `PMAIP1 = TP53`, `BBC3 = TP53` | reciprocal records removed | `BCL2L1 -\| PMAIP1` and `MCL1 -\| BBC3` are reverse records of sequestration, already handled as antagonist→guard elsewhere. Retaining both directions creates mutual-inhibition toggles. |
| `RAF1 = KRAS`, `MAP2K1 = RAF1` | ERK feedback and AKT⊣RAF1 removed | These are modulatory phosphorylations that tune amplitude and adaptation. Read as Boolean switches they (a) shut the MAPK arm off at baseline, since AKT is constitutively active in A549, and (b) create a sustained oscillator. Neither matches a cell with constitutive ERK activity. |
| `PI3K = KRAS and not (TP53 or PTEN)` | FOXO3→PI3K removed | Creates an AKT/FOXO3 negative-feedback oscillator that made the KRAS-degraded condition non-deterministic. |
| `TP53 = not (MDM2 and MDM4)` | AND between suppressors | Under OR with MDM4 clamped ON, p53 can never activate and MDM2 degradation is inert — a degenerate model. Also contradicted by the Nutlin data (p53 rises 4.8x with MDM4 present). |

**Network change made for the same reason:** `KRAS → PI3K` was restored as a
documented manual addition. It is single-resource in OmniPath and was excluded
by the evidence threshold, leaving the entire PI3K/AKT arm disconnected from
the perturbed node. Without it, KRAS degradation cannot reach AKT, BAD, FOXO3
or GSK3B, and the model cannot address the question it was built for.

## 4. Tier 3 — rules fitted to the validation data

**These three were chosen because they reproduce the experimental anchors. They
are not independent predictions and must not be presented as such.**

| rule | fitted choice | what it was fitted to |
|---|---|---|
| `BCL2L1 = not (BAD and (BBC3 or PMAIP1))` | AND between a sensitiser and an activator | Under the default (OR), either arm alone clears BCL-XL and both single agents kill — contradicting both datasets. The AND makes BCL-XL a coincidence detector requiring input from both arms. Biologically defensible (BAD is a sensitiser that displaces but does not activate; PUMA/NOXA are activators), but the AND was selected to fit. |
| `BBC3 = TP53` (FOXO3→BBC3 unused) | PUMA made p53-dependent | With FOXO3 able to induce PUMA, KRAS degradation alone clears BCL-XL and kills. The proteomics supports p53-dependence (PUMA rises under Nutlin, unchanged under RAS degradation), but this is the data being used twice. |
| `CDKN1A = (TP53 or FOXO3) and not MDM2` | trivial inhibitors dropped | MAPK1⊣CDKN1A (effort 5) and GSK3B⊣CDKN1A (effort 15) vetoed TP53→CDKN1A (effort 3452) under the default, so p21 could not rise while ERK was active. The evidence asymmetry justifies it independently, but the fit motivated the inspection. |

**Degrees of freedom.** Three rules were tuned against two experimental
observations. The model therefore *accommodates* the anchors rather than
predicting them. The synergy result is a hypothesis generated under those
assumptions, and the honest claim is: given a coincidence-detector BCL-XL and
p53-dependent PUMA, combined MDM2 and KRAS degradation is predicted to reach
apoptosis while neither single agent does.

## 5. Initial conditions

**The outcome does not depend on them.** Each of the four conditions was run
from 120 random initial states across all 29 nodes. Every condition converged
to a single attractor in 120/120 runs. The perturbation, not the starting
point, determines the result, and no cyclic attractors were found.

Initial states in the .cfg files are nevertheless set to the model's own
control attractor rather than to zeros, because that mirrors the experiment:
cells are growing in an untreated steady state when the drug is added.

Control attractor (ON): KRAS, PI3K, AKT1, RAF1, MAP2K1, MAPK1, MAPK3, MCL1,
BCL2L1, XIAP, MDM2, MDM4. Everything else OFF.

This is itself a check on the rules: it is a proliferating A549 cell with
constitutive KRAS/ERK signalling, active AKT, both mitochondrial guards intact,
XIAP restraining the caspases, p53 held off by MDM2, and no caspase activity.
It was not imposed — it emerges from the rules, and a rule set that failed to
produce it would be wrong before any perturbation was tested.

**Fixed nodes.** KRAS and MDM2 are clamped (1 = present, 0 = degraded) because
they are the perturbations. MDM4 is clamped ON because A549 expresses it and
MDM2-directed degraders do not target it; its in-edges were removed for this
reason.

## 6. What a reader should take from this

The control attractor and the single-attractor result are properties of the
network and convention. The Tier 2 overrides are mechanism, defensible on their
own terms. The Tier 3 rules are fitted, and the synergy prediction rests on
them. The discriminating experiment is therefore the combination arm, which no
dataset used here can speak to.
