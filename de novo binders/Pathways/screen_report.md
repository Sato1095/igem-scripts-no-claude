# Rule-variant screen — results

Criteria: `screen_criteria_v2.md` (fixed before running).
Method: exact value percolation, 32 variants x 4 conditions. No sampling.
Base rules: `rules_v3_unfitted.txt`. Script: `run_screen.py`.

**4 of 32 variants satisfy all criteria.**

## 1. Full grid

| variant | pass | failed criteria |
|---|---|---|
| `R1AR2AR3AR4AR5A` | fail | MDM2_degraded.TP53=0(want 1); MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BBC3=0(want 1) |
| `R1AR2AR3AR4AR5B` | fail | MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BCL2L1=0(want 1) |
| `R1AR2AR3AR4BR5A` | fail | MDM2_degraded.TP53=0(want 1); MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BBC3=0(want 1) |
| `R1AR2AR3AR4BR5B` | fail | MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BCL2L1=0(want 1) |
| `R1AR2AR3BR4AR5A` | fail | MDM2_degraded.TP53=0(want 1); MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BBC3=0(want 1) |
| `R1AR2AR3BR4AR5B` | fail | MDM2_degraded.BCL2L1=0(want 1) |
| `R1AR2AR3BR4BR5A` | fail | MDM2_degraded.TP53=0(want 1); MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BBC3=0(want 1) |
| `R1AR2AR3BR4BR5B` | fail | MDM2_degraded.BCL2L1=0(want 1) |
| `R1AR2BR3AR4AR5A` | fail | MDM2_degraded.TP53=0(want 1); MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BBC3=0(want 1) |
| `R1AR2BR3AR4AR5B` | fail | MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BCL2L1=0(want 1) |
| `R1AR2BR3AR4BR5A` | fail | MDM2_degraded.TP53=0(want 1); MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BBC3=0(want 1) |
| `R1AR2BR3AR4BR5B` | fail | MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BCL2L1=0(want 1) |
| `R1AR2BR3BR4AR5A` | fail | MDM2_degraded.TP53=0(want 1); MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BBC3=0(want 1) |
| `R1AR2BR3BR4AR5B` | fail | MDM2_degraded.BCL2L1=0(want 1) |
| `R1AR2BR3BR4BR5A` | fail | MDM2_degraded.TP53=0(want 1); MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BBC3=0(want 1) |
| `R1AR2BR3BR4BR5B` | fail | MDM2_degraded.BCL2L1=0(want 1) |
| `R1BR2AR3AR4AR5A` | fail | MDM2_degraded.TP53=0(want 1); MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BBC3=0(want 1) |
| `R1BR2AR3AR4AR5B` | fail | MDM2_degraded.CDKN1A=0(want 1) |
| `R1BR2AR3AR4BR5A` | fail | MDM2_degraded.TP53=0(want 1); MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BBC3=0(want 1) |
| `R1BR2AR3AR4BR5B` | fail | MDM2_degraded.CDKN1A=0(want 1) |
| `R1BR2AR3BR4AR5A` | fail | MDM2_degraded.TP53=0(want 1); MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BBC3=0(want 1) |
| `R1BR2AR3BR4AR5B` | **PASS** | — |
| `R1BR2AR3BR4BR5A` | fail | MDM2_degraded.TP53=0(want 1); MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BBC3=0(want 1) |
| `R1BR2AR3BR4BR5B` | **PASS** | — |
| `R1BR2BR3AR4AR5A` | fail | MDM2_degraded.TP53=0(want 1); MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BBC3=0(want 1) |
| `R1BR2BR3AR4AR5B` | fail | MDM2_degraded.CDKN1A=0(want 1) |
| `R1BR2BR3AR4BR5A` | fail | MDM2_degraded.TP53=0(want 1); MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BBC3=0(want 1) |
| `R1BR2BR3AR4BR5B` | fail | MDM2_degraded.CDKN1A=0(want 1) |
| `R1BR2BR3BR4AR5A` | fail | MDM2_degraded.TP53=0(want 1); MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BBC3=0(want 1) |
| `R1BR2BR3BR4AR5B` | **PASS** | — |
| `R1BR2BR3BR4BR5A` | fail | MDM2_degraded.TP53=0(want 1); MDM2_degraded.CDKN1A=0(want 1); MDM2_degraded.BBC3=0(want 1) |
| `R1BR2BR3BR4BR5B` | **PASS** | — |

## 2. Phenotype across all variants

D = apoptosis (CASP3 = CASP7 = 1)

| variant | pass | control | MDM2 | KRAS | BOTH |
|---|---|---|---|---|---|
| `R1AR2AR3AR4AR5A` |  | . | . | D | D |
| `R1AR2AR3AR4AR5B` |  | . | . | D | D |
| `R1AR2AR3AR4BR5A` |  | . | . | D | D |
| `R1AR2AR3AR4BR5B` |  | . | . | D | D |
| `R1AR2AR3BR4AR5A` |  | . | . | D | D |
| `R1AR2AR3BR4AR5B` |  | . | . | D | D |
| `R1AR2AR3BR4BR5A` |  | . | . | D | D |
| `R1AR2AR3BR4BR5B` |  | . | . | D | D |
| `R1AR2BR3AR4AR5A` |  | . | . | D | D |
| `R1AR2BR3AR4AR5B` |  | . | . | D | D |
| `R1AR2BR3AR4BR5A` |  | . | . | D | D |
| `R1AR2BR3AR4BR5B` |  | . | . | D | D |
| `R1AR2BR3BR4AR5A` |  | . | . | D | D |
| `R1AR2BR3BR4AR5B` |  | . | . | D | D |
| `R1AR2BR3BR4BR5A` |  | . | . | D | D |
| `R1AR2BR3BR4BR5B` |  | . | . | D | D |
| `R1BR2AR3AR4AR5A` |  | . | . | D | D |
| `R1BR2AR3AR4AR5B` |  | . | . | D | D |
| `R1BR2AR3AR4BR5A` |  | . | . | D | D |
| `R1BR2AR3AR4BR5B` |  | . | . | D | D |
| `R1BR2AR3BR4AR5A` |  | . | . | D | D |
| `R1BR2AR3BR4AR5B` | yes | . | . | D | D |
| `R1BR2AR3BR4BR5A` |  | . | . | D | D |
| `R1BR2AR3BR4BR5B` | yes | . | . | D | D |
| `R1BR2BR3AR4AR5A` |  | . | . | . | . |
| `R1BR2BR3AR4AR5B` |  | . | . | . | D |
| `R1BR2BR3AR4BR5A` |  | . | . | . | . |
| `R1BR2BR3AR4BR5B` |  | . | . | . | D |
| `R1BR2BR3BR4AR5A` |  | . | . | . | . |
| `R1BR2BR3BR4AR5B` | yes | . | . | . | D |
| `R1BR2BR3BR4BR5A` |  | . | . | . | . |
| `R1BR2BR3BR4BR5B` | yes | . | . | . | D |

Apoptosis reached in: control 0/32, MDM2 alone 0/32, KRAS alone 24/32, combination 28/32.

**True synergy** (combination kills, neither single agent does): 4/32.

## 3. Attractors of the passing variants

### `R1BR2AR3BR4AR5B`

| condition | TP53 | CDKN1A | BCL2L1 | MCL1 | BAX | CYCS | CASP3 | CASP7 |
|---|---|---|---|---|---|---|---|---|
| control | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 |
| MDM2_degraded | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| KRAS_degraded | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 1 |
| both_degraded | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 1 |

### `R1BR2AR3BR4BR5B`

| condition | TP53 | CDKN1A | BCL2L1 | MCL1 | BAX | CYCS | CASP3 | CASP7 |
|---|---|---|---|---|---|---|---|---|
| control | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 |
| MDM2_degraded | 1 | 1 | 1 | 0 | 1 | 1 | 0 | 0 |
| KRAS_degraded | 0 | 0 | 0 | 0 | 1 | 1 | 1 | 1 |
| both_degraded | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 1 |

### `R1BR2BR3BR4AR5B`

| condition | TP53 | CDKN1A | BCL2L1 | MCL1 | BAX | CYCS | CASP3 | CASP7 |
|---|---|---|---|---|---|---|---|---|
| control | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 |
| MDM2_degraded | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| KRAS_degraded | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 0 |
| both_degraded | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 1 |

### `R1BR2BR3BR4BR5B`

| condition | TP53 | CDKN1A | BCL2L1 | MCL1 | BAX | CYCS | CASP3 | CASP7 |
|---|---|---|---|---|---|---|---|---|
| control | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 0 |
| MDM2_degraded | 1 | 1 | 1 | 0 | 1 | 1 | 0 | 0 |
| KRAS_degraded | 0 | 0 | 1 | 0 | 1 | 1 | 0 | 0 |
| both_degraded | 1 | 1 | 0 | 0 | 1 | 1 | 1 | 1 |

