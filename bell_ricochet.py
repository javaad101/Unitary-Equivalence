"""
Bell-State Ricochet & Unitary Contextual Equivalence (Zero-Dependency Reference)
Author: Javaad Ali

Self-contained complex Hilbert-space implementation demonstrating:
  1. The Bell-state Ricochet identity: (A x I)|Phi^+> = (I x A^T)|Phi^+>
  2. Post-selected Bell overlap = normalized Hilbert-Schmidt inner product:
         <Phi^+| (U1 x U2*) |Phi^+> = (1/d) Tr(U1 U2^dagger)
  3. Why unitary operators in U(2^n) preserve syntactic/contextual degrees of
     freedom by rotating non-invariant structure into orthogonal subspaces,
     unlike rank-deficient scalar projections (`count = sum(bits)`).
"""

from __future__ import annotations
import cmath
import math
import random

Matrix = list[list[complex]]
Vector = list[complex]


def zeros_vec(n: int) -> Vector:
    return [0j] * n


def eye(n: int) -> Matrix:
    return [[1.0 + 0j if i == j else 0j for j in range(n)] for i in range(n)]


def transpose(A: Matrix) -> Matrix:
    return [[A[i][j] for i in range(len(A))] for j in range(len(A[0]))]


def conj_matrix(A: Matrix) -> Matrix:
    return [[val.conjugate() for val in row] for row in A]


def dagger(A: Matrix) -> Matrix:
    return [[A[i][j].conjugate() for i in range(len(A))] for j in range(len(A[0]))]


def matmul(A: Matrix, B: Matrix) -> Matrix:
    m, k, n = len(A), len(B), len(B[0])
    out = [[0j] * n for _ in range(m)]
    for i in range(m):
        row_a = A[i]
        out_row = out[i]
        for p in range(k):
            a_ip = row_a[p]
            if a_ip != 0j:
                row_b = B[p]
                for j in range(n):
                    out_row[j] += a_ip * row_b[j]
    return out


def matvec(A: Matrix, v: Vector) -> Vector:
    return [sum(A[i][j] * v[j] for j in range(len(v))) for i in range(len(A))]


def vdot(u: Vector, v: Vector) -> complex:
    return sum(u[i].conjugate() * v[i] for i in range(len(u)))


def l2_norm_diff(u: Vector, v: Vector) -> float:
    return math.sqrt(sum(abs(u[i] - v[i]) ** 2 for i in range(len(u))))


def kron(A: Matrix, B: Matrix) -> Matrix:
    ra, ca = len(A), len(A[0])
    rb, cb = len(B), len(B[0])
    out = [[0j] * (ca * cb) for _ in range(ra * rb)]
    for i in range(ra):
        for j in range(ca):
            a_ij = A[i][j]
            if a_ij == 0j:
                continue
            base_r, base_c = i * rb, j * cb
            for p in range(rb):
                out_row = out[base_r + p]
                row_b = B[p]
                for q in range(cb):
                    out_row[base_c + q] = a_ij * row_b[q]
    return out


def trace(A: Matrix) -> complex:
    return sum(A[i][i] for i in range(len(A)))


def random_unitary(d: int, seed: int = 42) -> Matrix:
    """Samples a unitary matrix in U(d) via complex Gram-Schmidt QR orthonormalization."""
    rng = random.Random(seed)
    cols: list[Vector] = []
    for _ in range(d):
        v = [complex(rng.gauss(0.0, 1.0), rng.gauss(0.0, 1.0)) for _ in range(d)]
        for basis in cols:
            proj = vdot(basis, v)
            v = [v[i] - proj * basis[i] for i in range(d)]
        norm = math.sqrt(sum(abs(x) ** 2 for x in v))
        cols.append([x / norm for x in v])
    # Convert column list to row-major matrix
    return [[cols[j][i] for j in range(d)] for i in range(d)]


def bell_state(d: int) -> Vector:
    """Constructs maximally entangled generalized Bell state |Phi^+> in C^(d^2)."""
    psi = zeros_vec(d * d)
    amp = 1.0 / math.sqrt(d)
    for i in range(d):
        psi[i * d + i] = complex(amp, 0.0)
    return psi


def ricochet_lhs(A: Matrix) -> Vector:
    """Evaluates (A x I)|Phi^+>."""
    d = len(A)
    return matvec(kron(A, eye(d)), bell_state(d))


def ricochet_rhs(A: Matrix) -> Vector:
    """Evaluates (I x A^T)|Phi^+>."""
    d = len(A)
    return matvec(kron(eye(d), transpose(A)), bell_state(d))


def postselected_bell_overlap(U1: Matrix, U2: Matrix) -> complex:
    """Evaluates <Phi^+| (U1 x U2*) |Phi^+>."""
    d = len(U1)
    phi = bell_state(d)
    op = kron(U1, conj_matrix(U2))
    return vdot(phi, matvec(op, phi))


def hilbert_schmidt_overlap(U1: Matrix, U2: Matrix) -> complex:
    """Analytical Hilbert-Schmidt inner product (1/d) Tr(U1 U2^dagger)."""
    d = len(U1)
    return trace(matmul(U1, dagger(U2))) / d


def encode_binary_features_to_state(bits: list[int]) -> Vector:
    """Maps an n-bit discrete specification losslessly into C^(2^n)."""
    dim = 1 << len(bits)
    idx = 0
    for b in bits:
        idx = (idx << 1) | (b & 1)
    vec = zeros_vec(dim)
    vec[idx] = 1.0 + 0j
    return vec


def apply_single_qubit_ry(state: Vector, n_qubits: int, target_qubit: int, theta: float) -> Vector:
    """Applies single-qubit R_y(theta) unitary in O(2^n) directly in C^(2^n)."""
    c = math.cos(theta / 2.0)
    s = math.sin(theta / 2.0)
    bit_pos = n_qubits - 1 - target_qubit
    step = 1 << bit_pos
    out = list(state)
    for i in range(0, len(state), step * 2):
        for j in range(i, i + step):
            a0 = state[j]
            a1 = state[j + step]
            out[j] = c * a0 - s * a1
            out[j + step] = s * a0 + c * a1
    return out


def run_contextual_equivalence_demo() -> dict[str, float]:
    """
    8-qubit (256-dimensional Hilbert space) experiment contrasting:
      - Lossy scalar contraction (`count = sum(bits)`), erasing bit-identity
      - Lossless Hilbert-space encoding under unitary transformations and
        context-conditioned observables.
    """
    spec_a = [1, 1, 1, 1, 1, 1, 1, 0]  # 7 active features (q7 = 0)
    spec_b = [1, 1, 1, 1, 1, 1, 0, 1]  # 7 active features (q7 = 1)

    psi_a = encode_binary_features_to_state(spec_a)
    psi_b = encode_binary_features_to_state(spec_b)

    # Context 1 observable: identity on q7 (invariant to q6 <-> q7 swap)
    yield_ctx1_a = vdot(psi_a, psi_a).real
    yield_ctx1_b = vdot(psi_b, psi_b).real

    # Context 2 observable: projector P_{q7=0} (annihilates states with q7 = 1)
    proj_a = [psi_a[idx] if (idx & 1) == 0 else 0j for idx in range(256)]
    proj_b = [psi_b[idx] if (idx & 1) == 0 else 0j for idx in range(256)]
    yield_ctx2_a = vdot(psi_a, proj_a).real
    yield_ctx2_b = vdot(psi_b, proj_b).real

    # Unitarity check: rotating any register preserves exact L2 state separation
    rot_a = apply_single_qubit_ry(psi_a, 8, target_qubit=7, theta=math.pi / 3.0)
    rot_b = apply_single_qubit_ry(psi_b, 8, target_qubit=7, theta=math.pi / 3.0)

    return {
        "scalar_projection_a": float(sum(spec_a)),
        "scalar_projection_b": float(sum(spec_b)),
        "ctx1_yield_a": yield_ctx1_a,
        "ctx1_yield_b": yield_ctx1_b,
        "ctx2_yield_a": yield_ctx2_a,
        "ctx2_yield_b": yield_ctx2_b,
        "l2_distance_before_unitary": l2_norm_diff(psi_a, psi_b),
        "l2_distance_after_unitary": l2_norm_diff(rot_a, rot_b),
    }


if __name__ == "__main__":
    U = random_unitary(8, seed=42)
    ricochet_err = l2_norm_diff(ricochet_lhs(U), ricochet_rhs(U))
    overlap = postselected_bell_overlap(U, U)
    hs_diff = abs(overlap - hilbert_schmidt_overlap(U, U))
    demo = run_contextual_equivalence_demo()

    print(f"Ricochet identity ||(U x I)|Phi+> - (I x U^T)|Phi+>||_2 = {ricochet_err:.2e}")
    print(f"Post-selected Bell overlap <Phi+|(U x U*)|Phi+>         = {overlap.real:.12f}")
    print(f"Analytical (1/d) Tr(U U^dagger) discrepancy             = {hs_diff:.2e}")
    print("8-qubit (256-dim) contextual equivalence:", demo)
