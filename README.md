# Unitary Contextual Equivalence & the Bell-State Ricochet Identity

Zero-dependency numerical verification of a lossless representation-learning mechanism derived from post-selected quantum teleportation (P-CTC / Bell-state projection).

## Mathematical Formulation

In standard classical architectures, learning invariance across two distinct inputs (e.g., two syntactically distinct clauses or feature configurations that are equivalent under one context) is typically achieved by **rank-deficient contraction**—multiplying the distinguishing degrees of freedom by zero. Once erased, that distinction cannot be recovered when the evaluation context shifts.

In a unitary Hilbert-space parameterization ($U^\dagger U = I$ on $\mathcal{H} = (\mathbb{C}^2)^{\otimes n}$), information destruction is forbidden. Using the **Ricochet Property** of the maximally entangled state $|\Phi^+\rangle = \frac{1}{\sqrt{d}}\sum_{i=0}^{d-1} |i\rangle \otimes |i\rangle$:

$$(A \otimes I)|\Phi^+\rangle = (I \otimes A^T)|\Phi^+\rangle$$

applying $U_1$ to register 1 and $U_2^*$ to register 2 followed by Bell-state post-selection yields the normalized **Hilbert–Schmidt inner product**:

$$\langle \Phi^+ | (U_1 \otimes U_2^*) | \Phi^+ \rangle = \frac{1}{d}\,\mathrm{Tr}(U_1 U_2^\dagger)$$

which equals $1$ if and only if $U_1 = U_2$ (up to global phase). When equivalence is enforced only on a measured register while leaving ancilla/phase registers unmeasured ($\mathrm{Tr}_{-i}[\hat{O}_i \rho]$), the unitary operator rotates non-invariant structure into orthogonal subspaces rather than annihilating it—preserving full state resolution for context-dependent downstream observables.

## Running Locally

Zero third-party dependencies (standard library Python 3 only):

```bash
python3 bell_ricochet.py
python3 -m unittest test_bell_ricochet.py
```
