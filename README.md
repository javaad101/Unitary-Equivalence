# Plain English

We want to be able to say that a hot dog is a sandwich in some ways and that it is not a sandwich in others. We are forcing a machine to learn that two things can be equal for some purposes and non-equivalent for others. The ability to force equivalences and non-equivalences upon a machine is useful because it allows us to make a known representation of fact in language that is comprehensible to us as humans, and to force the machine to somehow internalize that fact, without necessarily having to express a view on how exactly it wishes to do so—it allows the machine to represent facts to itself in whatever manner is useful to the machine. The ability to do this allows us to convey statements in the nature of logic without having to make trillion-token models.

*Note: The idea is presented in a manner that I believe to be representative of substantially all of the logic that I am trying to capture. I've omitted more involved discussions for things that are not presently relevant here for economy of expression and time. The balance of this document was generated to summarize what is intended to be conveyed hereby with greater detail without having to figure out how to do so myself. I would've written you a shorter letter, but I didn't have the time.*

---

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
