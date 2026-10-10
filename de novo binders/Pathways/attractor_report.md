# Attractor analysis

Rules: `rules_v3_unfitted.txt`  (29 nodes, state space 2^29 = 536,870,912)


## Condition: control
Clamped: KRAS=1, MDM2=1, MDM4=1  -> reduced network of 26 free nodes (67,108,864 states)

**Percolation** converged in 6 rounds; fixed 29/29 nodes, 0 left free.

Every node is determined, so the clamped subspace percolates to a single state. That state is a fixed point, it is the ONLY attractor, and it is reachable from every initial state. This is exact.

ON: AKT1, BCL2L1, KRAS, MAP2K1, MAPK1, MAPK3, MCL1, MDM2, MDM4, PI3K, RAF1, XIAP
OFF: APAF1, BAD, BAK1, BAX, BBC3, BCL2L11, CASP3, CASP7, CASP9, CDKN1A, CYCS, DIABLO, FOXO3, GSK3B, PMAIP1, PTEN, TP53

**Random-start sampling** (120 uniformly random initial states, 6000 asynchronous updates each): 1 distinct end state(s).
Sampling agrees with the exact result: yes.

Note: sampling covers 120 of 67,108,864 states (1.79e-06 of the space). It is corroboration, not proof; the percolation result above is the exact one.

## Condition: MDM2_degraded
Clamped: KRAS=1, MDM2=0, MDM4=1  -> reduced network of 26 free nodes (67,108,864 states)

**Percolation** converged in 5 rounds; fixed 29/29 nodes, 0 left free.

Every node is determined, so the clamped subspace percolates to a single state. That state is a fixed point, it is the ONLY attractor, and it is reachable from every initial state. This is exact.

ON: APAF1, BAK1, BAX, BBC3, CYCS, DIABLO, KRAS, MAP2K1, MAPK1, MAPK3, MDM4, PMAIP1, PTEN, RAF1, TP53
OFF: AKT1, BAD, BCL2L1, BCL2L11, CASP3, CASP7, CASP9, CDKN1A, FOXO3, GSK3B, MCL1, MDM2, PI3K, XIAP

**Random-start sampling** (120 uniformly random initial states, 6000 asynchronous updates each): 1 distinct end state(s).
Sampling agrees with the exact result: yes.

Note: sampling covers 120 of 67,108,864 states (1.79e-06 of the space). It is corroboration, not proof; the percolation result above is the exact one.

## Condition: KRAS_degraded
Clamped: KRAS=0, MDM2=1, MDM4=1  -> reduced network of 26 free nodes (67,108,864 states)

**Percolation** converged in 6 rounds; fixed 29/29 nodes, 0 left free.

Every node is determined, so the clamped subspace percolates to a single state. That state is a fixed point, it is the ONLY attractor, and it is reachable from every initial state. This is exact.

ON: APAF1, BAD, BAK1, BAX, BBC3, BCL2L11, CASP3, CASP7, CASP9, CYCS, DIABLO, FOXO3, GSK3B, MDM2, MDM4
OFF: AKT1, BCL2L1, CDKN1A, KRAS, MAP2K1, MAPK1, MAPK3, MCL1, PI3K, PMAIP1, PTEN, RAF1, TP53, XIAP

**Random-start sampling** (120 uniformly random initial states, 6000 asynchronous updates each): 1 distinct end state(s).
Sampling agrees with the exact result: yes.

Note: sampling covers 120 of 67,108,864 states (1.79e-06 of the space). It is corroboration, not proof; the percolation result above is the exact one.

## Condition: both_degraded
Clamped: KRAS=0, MDM2=0, MDM4=1  -> reduced network of 26 free nodes (67,108,864 states)

**Percolation** converged in 6 rounds; fixed 29/29 nodes, 0 left free.

Every node is determined, so the clamped subspace percolates to a single state. That state is a fixed point, it is the ONLY attractor, and it is reachable from every initial state. This is exact.

ON: APAF1, BAD, BAK1, BAX, BBC3, BCL2L11, CASP3, CASP7, CASP9, CYCS, DIABLO, FOXO3, GSK3B, MDM4, PMAIP1, TP53
OFF: AKT1, BCL2L1, CDKN1A, KRAS, MAP2K1, MAPK1, MAPK3, MCL1, MDM2, PI3K, PTEN, RAF1, XIAP

**Random-start sampling** (120 uniformly random initial states, 6000 asynchronous updates each): 1 distinct end state(s).
Sampling agrees with the exact result: yes.

Note: sampling covers 120 of 67,108,864 states (1.79e-06 of the space). It is corroboration, not proof; the percolation result above is the exact one.
