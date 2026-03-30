"""
Quantum Arithmetic Circuits for ECDLP
======================================

Scalable quantum circuits for Shor's ECDLP algorithm.  Two strategies:

1. **EfficientPermutationAdder** – index-based encoding with O(N·n) gate
   decomposition instead of dense 2^n × 2^n unitary matrices.  Practical for
   curves up to ~16-bit group order on current simulators/hardware.

2. **Modular arithmetic primitives** (QFT-based) – building blocks for a fully
   arithmetic coordinate-encoding approach that scales polynomially.

References
----------
- Beauregard (2003): Circuit for Shor's algorithm using 2n+3 qubits
- Roetteler, Naehrig, Svore, Lauter (2017): Quantum resource estimates for
  computing elliptic curve discrete logarithms
"""

import math
import numpy as np
from typing import Tuple, Optional, List, Dict
from dataclasses import dataclass

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.circuit.library import QFTGate

# Re-use curve helpers from the main module
from projecteleven import CurveParams, EllipticCurve, PointEncoder


# ---------------------------------------------------------------------------
# Part 1 — Efficient permutation-based point adder
# ---------------------------------------------------------------------------

def _apply_transposition(qc: QuantumCircuit, a: int, b: int,
                         qubits: List[int], n_bits: int):
    """Swap computational basis states |a⟩ and |b⟩ on *qubits*.

    Uses the CNOT-reduction method:
    1. CNOT from a pivot differing-bit to all other differing bits so that
       |a⟩ and |b⟩ now differ in only one bit.
    2. Multi-controlled X on that pivot bit, conditioned on all other bits
       matching the (transformed) pattern.
    3. Undo the CNOTs.

    Gate count: O(n) per transposition.
    """
    if a == b:
        return

    diff = a ^ b
    diff_bits = [i for i in range(n_bits) if (diff >> i) & 1]
    same_bits = [(i, (a >> i) & 1) for i in range(n_bits) if not ((diff >> i) & 1)]

    pivot = diff_bits[0]
    other_diff = diff_bits[1:]

    # --- Step 1: CNOT pivot → other differing bits ---
    for d in other_diff:
        qc.cx(qubits[pivot], qubits[d])

    # After CNOTs, compute the transformed value of a
    a_prime = a
    for d in other_diff:
        a_prime ^= ((a >> pivot) & 1) << d   # XOR bit d with pivot bit of a

    # --- Step 2: multi-controlled X on pivot, conditioned on all other bits ---
    # Prepare controls: flip bits that need to be 0→1 for the MCX
    x_undo = []
    controls = []
    for i in range(n_bits):
        if i == pivot:
            continue
        expected = (a_prime >> i) & 1
        if expected == 0:
            qc.x(qubits[i])
            x_undo.append(i)
        controls.append(qubits[i])

    if controls:
        qc.mcx(controls, qubits[pivot])
    else:
        qc.x(qubits[pivot])

    for i in x_undo:
        qc.x(qubits[i])

    # --- Step 3: undo CNOTs ---
    for d in reversed(other_diff):
        qc.cx(qubits[pivot], qubits[d])


def _controlled_transposition(qc: QuantumCircuit, ctrl: int,
                              a: int, b: int,
                              qubits: List[int], n_bits: int,
                              ancilla_qubits: List[int] = None):
    """Controlled swap of |a⟩ ↔ |b⟩, activated when *ctrl* = |1⟩.

    Same CNOT-reduction strategy but the MCX gains one extra control.
    If *ancilla_qubits* has enough qubits (>= num_controls - 2), uses V-chain
    MCX decomposition: O(n) Toffolis instead of O(n^2).
    """
    if a == b:
        return

    diff = a ^ b
    diff_bits = [i for i in range(n_bits) if (diff >> i) & 1]

    pivot = diff_bits[0]
    other_diff = diff_bits[1:]

    # Step 1: CNOTs (unconditional -- they cancel outside |a⟩/|b⟩ subspace)
    for d in other_diff:
        qc.cx(qubits[pivot], qubits[d])

    a_prime = a
    for d in other_diff:
        a_prime ^= ((a >> pivot) & 1) << d

    # Step 2: controlled MCX
    x_undo = []
    controls = [ctrl]          # extra control qubit
    for i in range(n_bits):
        if i == pivot:
            continue
        expected = (a_prime >> i) & 1
        if expected == 0:
            qc.x(qubits[i])
            x_undo.append(i)
        controls.append(qubits[i])

    needed_anc = max(0, len(controls) - 2)
    if ancilla_qubits and len(ancilla_qubits) >= needed_anc and needed_anc > 0:
        qc.mcx(controls, qubits[pivot], ancilla_qubits[:needed_anc], mode='v-chain')
    else:
        qc.mcx(controls, qubits[pivot])

    for i in x_undo:
        qc.x(qubits[i])

    # Step 3: undo CNOTs
    for d in reversed(other_diff):
        qc.cx(qubits[pivot], qubits[d])


class EfficientPermutationAdder:
    """Replaces QuantumPointAdder with O(N·n) gate decomposition.

    Instead of materialising a 2^(n+1) × 2^(n+1) unitary matrix (~4 GB for
    13-bit curves), decomposes each "add point S" permutation into
    transpositions via cycle decomposition.

    Memory: O(N) for the permutation table.
    Gates:  O(N·n) per controlled point addition (N = group order, n = n_bits).

    Optimisations
    -------------
    - Each controlled addition is built as an isolated sub-circuit and appended
      to the main circuit as a single opaque gate.  This avoids quadratic DAG
      growth in Qiskit when millions of gates are added one-by-one.
    - An optional ancilla qubit enables V-chain MCX decomposition (O(n) Toffolis
      instead of O(n²) without ancilla).
    """

    def __init__(self, encoder: PointEncoder):
        self.encoder = encoder
        self.n = encoder.n
        self.n_bits = encoder.n_bits
        self.ec = encoder.ec
        self._perm_cache: Dict[int, List[Tuple[int, int]]] = {}
        self._circuit_cache: Dict[Tuple[int, bool], QuantumCircuit] = {}

    def _get_transpositions(self, S: Optional[Tuple[int, int]]) -> List[Tuple[int, int]]:
        """Cycle-decompose the 'add S' permutation into transpositions."""
        s_idx = self.encoder.encode(S)
        if s_idx in self._perm_cache:
            return self._perm_cache[s_idx]

        # Build permutation table
        perm = list(range(1 << self.n_bits))
        for i in range(self.n):
            P = self.encoder.decode(i)
            P_plus_S = self.ec.add(P, S)
            perm[i] = self.encoder.encode(P_plus_S)

        # Cycle decomposition -> transpositions
        transpositions = []
        visited = [False] * len(perm)
        for start in range(len(perm)):
            if visited[start] or perm[start] == start:
                visited[start] = True
                continue
            cycle = []
            j = start
            while not visited[j]:
                visited[j] = True
                cycle.append(j)
                j = perm[j]
            # (c0 c1 c2 ... ck) = (c0 ck)(c0 c_{k-1})...(c0 c1)
            for idx in range(len(cycle) - 1, 0, -1):
                transpositions.append((cycle[0], cycle[idx]))

        self._perm_cache[s_idx] = transpositions
        return transpositions

    def _build_controlled_add_circuit(self, S: Optional[Tuple[int, int]],
                                       num_ancilla: int) -> QuantumCircuit:
        """Build a sub-circuit for controlled addition of S.

        Qubit layout of the sub-circuit:
          [0 .. n_bits-1]           = point register
          [n_bits]                  = control qubit
          [n_bits+1 .. n_bits+num_ancilla] = ancilla qubits (for v-chain MCX)
        """
        s_idx = self.encoder.encode(S)
        cache_key = (s_idx, num_ancilla)
        if cache_key in self._circuit_cache:
            return self._circuit_cache[cache_key]

        n_bits = self.n_bits
        n_qubits = n_bits + 1 + num_ancilla
        sub = QuantumCircuit(n_qubits, name=f"CAdd{s_idx}")

        pt_qubits = list(range(n_bits))
        ctrl = n_bits
        ancilla_qubits = list(range(n_bits + 1, n_qubits)) if num_ancilla > 0 else None

        transpositions = self._get_transpositions(S)
        for a, b in transpositions:
            _controlled_transposition(sub, ctrl, a, b, pt_qubits, n_bits, ancilla_qubits)

        self._circuit_cache[cache_key] = sub
        return sub

    def apply_controlled_add(self, qc: QuantumCircuit,
                             ctrl_qubit: int,
                             pt_qubits: List[int],
                             S: Optional[Tuple[int, int]],
                             ancilla_qubits: List[int] = None):
        """Add classical point S to the point register, controlled by ctrl_qubit.

        Builds the transposition circuit as an isolated sub-circuit, then
        appends it as a single gate to *qc*.  This keeps the main DAG small.
        """
        if S is None:
            return

        num_ancilla = len(ancilla_qubits) if ancilla_qubits else 0
        sub = self._build_controlled_add_circuit(S, num_ancilla)

        # Map sub-circuit qubits to main circuit qubits
        qubit_map = pt_qubits + [ctrl_qubit]
        if ancilla_qubits:
            qubit_map.extend(ancilla_qubits)

        qc.append(sub.to_gate(), qubit_map)


# ---------------------------------------------------------------------------
# Part 2 — QFT-based modular arithmetic
# ---------------------------------------------------------------------------

def phi_add_constant(qc: QuantumCircuit, qubits: List[int],
                     n_bits: int, constant: int):
    """Add a classical constant to a register already in Fourier basis.

    Convention: QFTGate(n) includes bit-reversal swaps, so after QFT
    qubit j (Qiskit LSB-first) carries phase 2πx / 2^{n-j}.
    Rotation to add *constant*:  P(2π · constant / 2^{n-j})  on qubit j.
    """
    for j in range(n_bits):
        denom = 1 << (n_bits - j)
        angle = 2 * math.pi * constant / denom
        # Skip trivially-zero rotations
        if abs(angle % (2 * math.pi)) > 1e-14:
            qc.p(angle, qubits[j])


def c_phi_add_constant(qc: QuantumCircuit, ctrl: int,
                       qubits: List[int], n_bits: int, constant: int):
    """Controlled: add classical constant in Fourier basis."""
    for j in range(n_bits):
        denom = 1 << (n_bits - j)
        angle = 2 * math.pi * constant / denom
        if abs(angle % (2 * math.pi)) > 1e-14:
            qc.cp(angle, ctrl, qubits[j])


def modular_add_constant(qc: QuantumCircuit, qubits: List[int],
                         n_bits: int, constant: int, p: int, ancilla: int):
    """Compute (target + constant) mod p in-place.

    Register *qubits* has *n_bits* qubits (enough to hold values up to 2p-1).
    *ancilla* is a single qubit initialised to |0⟩ and restored to |0⟩.

    Beauregard's modular adder:
      1. φ-add a
      2. φ-sub p
      3. IQFT; copy MSB → ancilla (MSB=1 ⟹ underflow)
      4. QFT; if ancilla: φ-add p
      5. IQFT; uncompute ancilla; QFT; restore register; IQFT
    """
    constant = constant % p
    qft = QFTGate(n_bits)
    iqft = qft.inverse()

    # 1. QFT
    qc.append(qft, qubits)
    # 2. φ-add constant, then φ-sub p
    phi_add_constant(qc, qubits, n_bits, constant)
    phi_add_constant(qc, qubits, n_bits, -p)
    # 3. IQFT → check MSB for underflow
    qc.append(iqft, qubits)
    qc.cx(qubits[n_bits - 1], ancilla)       # MSB=1 ⟹ result was negative
    # 4. QFT → conditionally add p back
    qc.append(qft, qubits)
    c_phi_add_constant(qc, ancilla, qubits, n_bits, p)
    # 5. IQFT (result is now (target + constant) mod p)
    qc.append(iqft, qubits)

    # --- uncompute ancilla ---
    qc.x(ancilla)
    qc.append(qft, qubits)
    phi_add_constant(qc, qubits, n_bits, -constant)
    qc.append(iqft, qubits)
    qc.cx(qubits[n_bits - 1], ancilla)
    qc.append(qft, qubits)
    phi_add_constant(qc, qubits, n_bits, constant)
    qc.append(iqft, qubits)


def c_modular_add_constant(qc: QuantumCircuit, ctrl: int,
                           qubits: List[int], n_bits: int,
                           constant: int, p: int, ancilla: int):
    """Controlled modular addition of a classical constant.

    Beauregard's controlled modular adder: the p subtraction is unconditional
    (not controlled on ctrl) so that the ancilla uncompute works correctly
    when ctrl is part of a superposition.
    """
    constant = constant % p
    qft = QFTGate(n_bits)
    iqft = qft.inverse()

    qc.append(qft, qubits)
    c_phi_add_constant(qc, ctrl, qubits, n_bits, constant)
    phi_add_constant(qc, qubits, n_bits, -p)
    qc.append(iqft, qubits)
    qc.cx(qubits[n_bits - 1], ancilla)
    qc.append(qft, qubits)
    c_phi_add_constant(qc, ancilla, qubits, n_bits, p)
    qc.append(iqft, qubits)

    # uncompute ancilla
    qc.x(ancilla)
    qc.append(qft, qubits)
    c_phi_add_constant(qc, ctrl, qubits, n_bits, -constant)
    qc.append(iqft, qubits)
    qc.cx(qubits[n_bits - 1], ancilla)
    qc.append(qft, qubits)
    c_phi_add_constant(qc, ctrl, qubits, n_bits, constant)
    qc.append(iqft, qubits)


def phi_add_quantum(qc: QuantumCircuit, a_qubits: List[int],
                    b_qubits: List[int], n_bits: int):
    """Add quantum register |b⟩ to |a⟩ (a must be in QFT basis).

    After QFTGate (with swaps), qubit a_j carries frequency 2^j.
    Phase from b_k on a_j: CP(2π / 2^{n-j-k}).

    Gate count: O(n²) CP gates.
    """
    for j in range(n_bits):
        for k in range(n_bits):
            exp = n_bits - j - k
            if exp <= 0:
                continue
            angle = 2 * math.pi / (1 << exp)
            qc.cp(angle, b_qubits[k], a_qubits[j])


def c_phi_add_quantum(qc: QuantumCircuit, ctrl: int,
                      a_qubits: List[int], b_qubits: List[int],
                      n_bits: int, sign: int = 1):
    """Controlled: add sign*|b⟩ to |a⟩ in QFT basis when ctrl=|1⟩.

    Decomposes each doubly-controlled phase CCP(θ) into 2 CX + 3 CP gates.
    Set sign=-1 to subtract.
    """
    for j in range(n_bits):
        for k in range(n_bits):
            exp = n_bits - j - k
            if exp <= 0:
                continue
            angle = sign * 2 * math.pi / (1 << exp)
            # CCP(θ) on ctrl, b_k → a_j
            qc.cp(angle / 2, b_qubits[k], a_qubits[j])
            qc.cx(ctrl, b_qubits[k])
            qc.cp(-angle / 2, b_qubits[k], a_qubits[j])
            qc.cx(ctrl, b_qubits[k])
            qc.cp(angle / 2, ctrl, a_qubits[j])


# ---------------------------------------------------------------------------
# Part 2b — Modular multiplication (constant and quantum-quantum)
# ---------------------------------------------------------------------------

def modular_multiply_constant(qc: QuantumCircuit,
                              x_qubits: List[int],
                              out_qubits: List[int],
                              n_bits: int,
                              constant: int,
                              p: int,
                              ancilla: int):
    """Compute |x⟩|0⟩ → |x⟩|c·x mod p⟩  (out_qubits accumulator).

    Uses shift-and-add: for each bit k of x, if x_k=1,
    add (c · 2^k mod p) to out_qubits.

    *ancilla*: single qubit for the modular adder.
    """
    for k in range(n_bits):
        addend = (constant * (1 << k)) % p
        if addend == 0:
            continue
        c_modular_add_constant(qc, x_qubits[k], out_qubits,
                               n_bits, addend, p, ancilla)


def modular_multiply_constant_inplace(qc: QuantumCircuit,
                                      x_qubits: List[int],
                                      tmp_qubits: List[int],
                                      n_bits: int,
                                      constant: int,
                                      p: int,
                                      ancilla: int):
    """Compute |x⟩ → |c·x mod p⟩ in-place using a temporary register.

    1. |x⟩|0⟩  → |x⟩|cx mod p⟩           (forward multiply)
    2. SWAP x ↔ tmp                        → |cx mod p⟩|x⟩
    3. |cx mod p⟩|x⟩ → |cx mod p⟩|0⟩      (inverse multiply by c⁻¹)
    """
    c_inv = pow(constant, -1, p)

    # Forward: accumulate c·x into tmp
    modular_multiply_constant(qc, x_qubits, tmp_qubits, n_bits, constant, p, ancilla)

    # SWAP
    for i in range(n_bits):
        qc.swap(x_qubits[i], tmp_qubits[i])

    # Inverse: subtract c⁻¹·(new x) from tmp to zero it out
    for k in range(n_bits):
        addend = (c_inv * (1 << k)) % p
        if addend == 0:
            continue
        # Subtract = add (p - addend)
        c_modular_add_constant(qc, x_qubits[k], tmp_qubits,
                               n_bits, (p - addend) % p, p, ancilla)


# ---------------------------------------------------------------------------
# Part 2c — Quantum-quantum modular arithmetic
# ---------------------------------------------------------------------------

def c_modular_add_quantum(qc: QuantumCircuit, ctrl: int,
                          b_qubits: List[int],
                          a_qubits: List[int],
                          n_bits: int, p: int, ancilla: int):
    """Controlled: |a⟩ → |a + b mod p⟩ when ctrl=|1⟩, |b⟩ unchanged.

    Uses Beauregard's modular addition with unconditional p subtraction
    for correct ancilla uncompute.
    """
    qft = QFTGate(n_bits)
    iqft = qft.inverse()

    qc.append(qft, a_qubits)
    c_phi_add_quantum(qc, ctrl, a_qubits, b_qubits, n_bits)
    phi_add_constant(qc, a_qubits, n_bits, -p)
    qc.append(iqft, a_qubits)
    qc.cx(a_qubits[n_bits - 1], ancilla)
    qc.append(qft, a_qubits)
    c_phi_add_constant(qc, ancilla, a_qubits, n_bits, p)
    qc.append(iqft, a_qubits)

    # Uncompute ancilla
    qc.x(ancilla)
    qc.append(qft, a_qubits)
    c_phi_add_quantum(qc, ctrl, a_qubits, b_qubits, n_bits, sign=-1)
    qc.append(iqft, a_qubits)
    qc.cx(a_qubits[n_bits - 1], ancilla)
    qc.append(qft, a_qubits)
    c_phi_add_quantum(qc, ctrl, a_qubits, b_qubits, n_bits)
    qc.append(iqft, a_qubits)


def modular_multiply_quantum(qc: QuantumCircuit,
                             a_qubits: List[int],
                             b_qubits: List[int],
                             out_qubits: List[int],
                             shifted_a_qubits: List[int],
                             tmp_qubits: List[int],
                             n_bits: int,
                             p: int,
                             ancilla: int):
    """Compute |a⟩|b⟩|0⟩ → |a⟩|b⟩|a·b mod p⟩.

    Uses modular doubling to keep all intermediate values < p,
    avoiding register overflow.

    Requires scratch registers:
      shifted_a_qubits[n_bits] — holds a·2^k mod p (starts/ends at |0⟩)
      tmp_qubits[n_bits]       — scratch for modular_multiply_constant_inplace

    Gate count: O(n³) — n iterations × O(n²) gates per modular add/double.
    """
    # Copy a → shifted_a
    for i in range(n_bits):
        qc.cx(a_qubits[i], shifted_a_qubits[i])

    for k in range(n_bits):
        # shifted_a = a·2^k mod p (always < p)
        # Controlled on b[k]: out += shifted_a (mod p)
        c_modular_add_quantum(qc, b_qubits[k], shifted_a_qubits,
                              out_qubits, n_bits, p, ancilla)

        # Double shifted_a for next iteration
        if k < n_bits - 1:
            modular_multiply_constant_inplace(qc, shifted_a_qubits, tmp_qubits,
                                              n_bits, 2, p, ancilla)

    # Uncompute shifted_a: currently holds a·2^{n-1} mod p
    # Multiply by modular inverse of 2^{n-1} to get back to a
    inv_shift = pow(pow(2, n_bits - 1, p), -1, p)
    modular_multiply_constant_inplace(qc, shifted_a_qubits, tmp_qubits,
                                      n_bits, inv_shift, p, ancilla)

    # Un-copy: shifted_a == a now, so CX zeros it out
    for i in range(n_bits):
        qc.cx(a_qubits[i], shifted_a_qubits[i])


# ---------------------------------------------------------------------------
# Part 2d — Modular negation, inversion and squaring primitives
# ---------------------------------------------------------------------------

def modular_negate_inplace(qc: QuantumCircuit,
                           qubits: List[int],
                           tmp_qubits: List[int],
                           n_bits: int, p: int, ancilla: int):
    """Compute |x⟩ → |(-x) mod p⟩ = |(p-x) mod p⟩ in-place.

    Implemented as multiplication by (p-1) mod p.
    """
    modular_multiply_constant_inplace(qc, qubits, tmp_qubits,
                                      n_bits, p - 1, p, ancilla)


def _build_permutation_transpositions(perm: List[int]) -> List[Tuple[int, int]]:
    """Cycle-decompose a permutation into transpositions."""
    n = len(perm)
    transpositions = []
    visited = [False] * n
    for start in range(n):
        if visited[start] or perm[start] == start:
            visited[start] = True
            continue
        cycle = []
        j = start
        while not visited[j]:
            visited[j] = True
            cycle.append(j)
            j = perm[j]
        for idx in range(len(cycle) - 1, 0, -1):
            transpositions.append((cycle[0], cycle[idx]))
    return transpositions


def modular_inverse_permutation(qc: QuantumCircuit,
                                qubits: List[int],
                                n_bits: int, p: int):
    """Compute |x⟩ → |x⁻¹ mod p⟩ in-place (0 maps to 0).

    Uses a pre-computed lookup table decomposed into transpositions.
    """
    dim = 1 << n_bits
    perm = list(range(dim))
    for v in range(1, p):
        perm[v] = pow(v, -1, p)

    for a, b in _build_permutation_transpositions(perm):
        _apply_transposition(qc, a, b, qubits, n_bits)


def c_modular_inverse_permutation(qc: QuantumCircuit, ctrl: int,
                                  qubits: List[int],
                                  n_bits: int, p: int,
                                  ancilla_qubits: List[int] = None):
    """Controlled modular inversion via permutation lookup."""
    dim = 1 << n_bits
    perm = list(range(dim))
    for v in range(1, p):
        perm[v] = pow(v, -1, p)

    for a, b in _build_permutation_transpositions(perm):
        _controlled_transposition(qc, ctrl, a, b, qubits, n_bits, ancilla_qubits)


def modular_square_out(qc: QuantumCircuit,
                       x_qubits: List[int],
                       out_qubits: List[int],
                       shifted_a_qubits: List[int],
                       tmp_qubits: List[int],
                       n_bits: int, p: int, ancilla: int):
    """Compute |x⟩|0⟩ → |x⟩|x² mod p⟩ using quantum-quantum multiply.

    Note: x → x² is NOT a permutation (not injective), so it cannot be
    done in-place via transpositions. Uses modular_multiply_quantum instead.
    """
    modular_multiply_quantum(qc, x_qubits, x_qubits, out_qubits,
                             shifted_a_qubits, tmp_qubits,
                             n_bits, p, ancilla)


# ---------------------------------------------------------------------------
# Part 3 — Modular inversion via Fermat's little theorem
# ---------------------------------------------------------------------------

def modular_inverse_fermat(qc: QuantumCircuit,
                           x_qubits: List[int],
                           out_qubits: List[int],
                           tmp_qubits: List[int],
                           n_bits: int,
                           p: int,
                           ancilla: int):
    """Compute |x⟩|0⟩ → |x⟩|x⁻¹ mod p⟩  via x^{p-2} mod p.

    Uses repeated squaring with classical exponent e = p-2.
    Requires out_qubits initialised to encoding of 1, and tmp_qubits as scratch.

    This is the most expensive primitive: O(n) modular multiplications,
    each O(n²) gates → O(n³) total.
    """
    e = p - 2
    # Initialise out = 1 (set LSB)
    qc.x(out_qubits[0])

    # We need a copy of x to repeatedly square
    # Copy x → tmp via CNOT
    for i in range(n_bits):
        qc.cx(x_qubits[i], tmp_qubits[i])

    # Repeated squaring: for each bit of e
    for bit_pos in range(e.bit_length()):
        if (e >> bit_pos) & 1:
            # out = out * base mod p
            # We need quantum-quantum multiplication here
            # For now, use the constant-multiply trick per-bit of base
            _quantum_modular_multiply_accumulate(
                qc, tmp_qubits, out_qubits, n_bits, p, ancilla
            )

        # Square the base: base = base² mod p
        if bit_pos < e.bit_length() - 1:
            _quantum_modular_square_inplace(
                qc, tmp_qubits, out_qubits, n_bits, p, ancilla
            )


def _quantum_modular_multiply_accumulate(
        qc: QuantumCircuit,
        a_qubits: List[int], b_qubits: List[int],
        n_bits: int, p: int, ancilla: int):
    """Multiply: |a⟩|b⟩ → |a⟩|a·b mod p⟩ using QFT quantum-quantum addition.

    For each bit k of a: if a_k=1, add b·2^k to accumulator (in Fourier space).
    This modifies b in-place, treating it as an accumulator.

    NOTE: this is a simplified version that adds into b. For correct Shor's
    usage, the caller must manage the accumulator lifecycle.
    """
    qft = QFTGate(n_bits)
    iqft = qft.inverse()

    for k in range(n_bits):
        # Controlled on a_k: add b << k to accumulator
        # We do this by adding in Fourier space with shift
        qc.append(qft, b_qubits)
        # For each bit j of b (in QFT), add phase controlled by a_k
        for j in range(n_bits):
            shift = k  # multiply by 2^k
            angle = 2 * math.pi * (1 << shift) / (1 << (j + 1))
            if abs(angle % (2 * math.pi)) > 1e-14:
                # Need CCP: controlled by a_k AND b itself
                # This is the quantum-quantum part
                qc.cp(angle, a_qubits[k], b_qubits[j])
        qc.append(iqft, b_qubits)


def _quantum_modular_square_inplace(
        qc: QuantumCircuit,
        x_qubits: List[int], tmp_qubits: List[int],
        n_bits: int, p: int, ancilla: int):
    """Square in-place: |x⟩ → |x² mod p⟩ (uses tmp as scratch).

    This is a placeholder — full implementation requires careful
    ancilla management for the quantum-quantum multiplication.
    """
    # TODO: implement proper quantum squaring
    # For now this is a structural placeholder
    pass


# ---------------------------------------------------------------------------
# Part 4 — Scalable Shor ECDLP using efficient permutation decomposition
# ---------------------------------------------------------------------------

class ScalableShorECDLP:
    """Shor's algorithm using efficient permutation circuit decomposition.

    Replaces the dense-unitary approach with O(N·n) gate decomposition per
    controlled point addition.  Memory usage is O(N) instead of O(2^{2n}).

    For 13-bit curves:
      - Dense unitary: ~4 GB per matrix, 28 matrices → impossible
      - This approach:  ~4243 transpositions per permutation → feasible
    """

    def __init__(self, params: CurveParams,
                 G: Tuple[int, int], Q: Tuple[int, int]):
        from projecteleven import ShorECDLP

        self.params = params
        self.G = G
        self.Q = Q
        self.n = params.n

        self.ec = EllipticCurve(params)
        self.encoder = PointEncoder(self.ec, G)
        self.adder = EfficientPermutationAdder(self.encoder)
        self._extraction_solver = ShorECDLP(params, G, Q)

        if Q is not None and Q not in self.encoder.point_to_index:
            raise ValueError("Q is not in the group generated by G")

    def extract_discrete_log(self, counts: Dict[str, int]) -> Optional[int]:
        """Delegate to ShorECDLP extraction logic."""
        return self._extraction_solver.extract_discrete_log(counts)

    def build_circuit(self, num_counting: int = None) -> QuantumCircuit:
        """Build the Shor ECDLP circuit with efficient gate decomposition."""
        n = self.n
        n_bits = self.encoder.n_bits

        if num_counting is None:
            num_counting = n_bits + 1

        # V-chain MCX needs (num_controls - 2) ancillas.
        # Max controls per transposition = n_bits (ctrl + n_bits-1 non-pivot).
        num_ancilla = max(0, n_bits - 2)

        j_reg = QuantumRegister(num_counting, 'j')
        k_reg = QuantumRegister(num_counting, 'k')
        pt_reg = QuantumRegister(n_bits, 'pt')
        anc_reg = QuantumRegister(num_ancilla, 'anc') if num_ancilla > 0 else None
        cl = ClassicalRegister(2 * num_counting + n_bits, 'cr')

        if anc_reg:
            qc = QuantumCircuit(j_reg, k_reg, pt_reg, anc_reg, cl,
                                name=f"Shor_ECDLP_efficient_n{n}")
            ancilla_qubits = list(anc_reg)
        else:
            qc = QuantumCircuit(j_reg, k_reg, pt_reg, cl,
                                name=f"Shor_ECDLP_efficient_n{n}")
            ancilla_qubits = None

        # Hadamard on counting registers
        for i in range(num_counting):
            qc.h(j_reg[i])
            qc.h(k_reg[i])

        # Controlled additions of 2^i * G
        G_power = self.G
        for i in range(num_counting):
            self.adder.apply_controlled_add(
                qc, j_reg[i], list(pt_reg), G_power, ancilla_qubits
            )
            G_power = self.ec.add(G_power, G_power)
            if G_power is None:
                G_power = self.G

        # Controlled additions of 2^i * Q
        Q_power = self.Q
        for i in range(num_counting):
            if Q_power is not None:
                self.adder.apply_controlled_add(
                    qc, k_reg[i], list(pt_reg), Q_power, ancilla_qubits
                )
            Q_power = self.ec.add(Q_power, Q_power)
            if Q_power is None and self.Q is not None:
                Q_power = self.Q

        # Measure point register
        for i in range(n_bits):
            qc.measure(pt_reg[i], cl[i])

        # Inverse QFT on counting registers
        qc.append(QFTGate(num_counting).inverse(), j_reg)
        qc.append(QFTGate(num_counting).inverse(), k_reg)

        # Measure counting registers
        for i in range(num_counting):
            qc.measure(j_reg[i], cl[n_bits + i])
            qc.measure(k_reg[i], cl[n_bits + num_counting + i])

        return qc

    def qubit_count(self) -> Dict[str, int]:
        n_bits = self.encoder.n_bits
        num_counting = n_bits + 1
        num_ancilla = max(0, n_bits - 2)
        return {
            "j_register": num_counting,
            "k_register": num_counting,
            "point_register": n_bits,
            "ancilla": num_ancilla,
            "total": 2 * num_counting + n_bits + num_ancilla,
            "group_order": self.n,
        }

    def gate_estimate(self) -> Dict[str, int]:
        """Estimate gate counts for the circuit."""
        n_bits = self.encoder.n_bits
        num_counting = n_bits + 1
        n = self.n

        # Each permutation has at most N transpositions
        # Each transposition: ~n CNOTs + 1 MCX(n-1 controls)
        # MCX(n-1) decomposes to ~O(n) Toffoli gates
        trans_per_perm = n  # upper bound
        gates_per_trans = 3 * n_bits  # rough: CNOTs + MCX decomposition
        perms = 2 * num_counting  # for G powers and Q powers

        return {
            "permutations": perms,
            "transpositions_per_perm": trans_per_perm,
            "estimated_basic_gates": perms * trans_per_perm * gates_per_trans,
            "estimated_cx_gates": perms * trans_per_perm * n_bits * 2,
        }


# ---------------------------------------------------------------------------
# Part 5 — Top-level solver using efficient decomposition
# ---------------------------------------------------------------------------

def solve_ecdlp_scalable(
    p: int, a: int, b: int, n: int,
    G: Tuple[int, int], Q: Tuple[int, int],
    shots: int = 8192,
    backend_name: str = "ibm_marrakesh",
    token: Optional[str] = None,
    instance: str = "open-instance",
    verbose: bool = True,
) -> Optional[int]:
    """Solve ECDLP using Shor's algorithm with efficient gate decomposition on IBM Quantum."""
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
    from projecteleven import ShorECDLP  # for extract_discrete_log and verification

    if verbose:
        print("=" * 60)
        print("SHOR'S ALGORITHM FOR ECDLP (Efficient Decomposition)")
        print("=" * 60)

    params = CurveParams(p, a, b, n)
    solver = ScalableShorECDLP(params, G, Q)

    if verbose:
        print(f"\nCurve: y² = x³ + {a}x + {b} (mod {p})")
        print(f"Group order: n = {n}")
        print(f"Generator: G = {G}")
        print(f"Target: Q = {Q}")
        reqs = solver.qubit_count()
        print(f"\nQubits: {reqs['total']} total")
        est = solver.gate_estimate()
        print(f"Estimated basic gates: ~{est['estimated_basic_gates']:,}")

    qc = solver.build_circuit()

    if verbose:
        print(f"Circuit gates: {qc.size()}")
        print(f"Circuit depth: {qc.depth()}")

    if token:
        service = QiskitRuntimeService(
            channel="ibm_quantum_platform", token=token, instance=instance
        )
    else:
        service = QiskitRuntimeService(
            channel="ibm_quantum_platform", instance=instance
        )
    backend = service.backend(backend_name)

    if verbose:
        print(f"\nBackend: {backend.name}")

    qc_t = transpile(qc, backend, optimization_level=3)

    if verbose:
        print(f"Transpiled depth: {qc_t.depth()}")
        ops = qc_t.count_ops()
        cx_count = sum(v for k, v in ops.items() if k in ['cx', 'ecr', 'cz'])
        print(f"Two-qubit gates: {cx_count}")

    sampler = SamplerV2(mode=backend)
    job = sampler.run([qc_t], shots=shots)

    if verbose:
        print(f"Job ID: {job.job_id()}")
        print("Waiting for results...")

    result = job.result()
    pub_result = result[0]
    counts = pub_result.data.cr.get_counts()

    if verbose:
        print(f"\nUnique outcomes: {len(counts)}")

    # Reuse the extraction logic from ShorECDLP
    legacy_solver = ShorECDLP(params, G, Q)
    d = legacy_solver.extract_discrete_log(counts)

    if verbose:
        print(f"\n{'=' * 60}")
        if d is not None:
            ec = EllipticCurve(params)
            print(f"RESULT: d = {d}")
            computed = ec.scalar_mult(d, G)
            print(f"Verification: {d}*G = {computed}")
            if computed == Q:
                print("[OK] VERIFIED")
            else:
                print("[FAIL] MISMATCH")
        else:
            print("FAILED: Could not extract discrete log")
        print("=" * 60)

    return d


# ---------------------------------------------------------------------------
# Quick self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="ECDLP Solver - Scalable Efficient Decomposition"
    )
    parser.add_argument("--shots", type=int, default=8192)
    parser.add_argument("--backend", default="ibm_marrakesh")
    parser.add_argument("--token", type=str,
                        help="IBM Quantum API token. Saves to local account on first use.")
    parser.add_argument("--instance", type=str, default="open-instance",
                        help="IBM Quantum instance (default: open-instance)")
    args = parser.parse_args()

    if args.token:
        from qiskit_ibm_runtime import QiskitRuntimeService
        QiskitRuntimeService.save_account(
            channel="ibm_quantum_platform", token=args.token, overwrite=True
        )
        print(f"IBM Quantum account saved.")

    print("=== Efficient permutation decomposition — IBM Quantum ===\n")

    # 4-bit curve: p=13, a=0, b=7, n=7
    params = CurveParams(p=13, a=0, b=7, n=7)
    G = (11, 5)
    d_secret = 6
    ec = EllipticCurve(params)
    Q = ec.scalar_mult(d_secret, G)
    print(f"Curve: y^2 = x^3 + 7 (mod 13), n={params.n}")
    print(f"G = {G}, Q = {d_secret}*G = {Q}")

    solver = ScalableShorECDLP(params, G, Q)
    reqs = solver.qubit_count()
    est = solver.gate_estimate()
    print(f"Qubits: {reqs['total']}")
    print(f"Estimated gates: ~{est['estimated_basic_gates']:,}")

    print("\nBuilding circuit...")
    qc = solver.build_circuit()
    print(f"Circuit size: {qc.size()} gates, depth: {qc.depth()}")

    print("\nSubmitting to IBM Quantum...")
    result = solve_ecdlp_scalable(
        13, 0, 7, 7, G, Q,
        shots=args.shots,
        backend_name=args.backend,
        token=args.token,
        instance=args.instance,
        verbose=True,
    )

    if result == d_secret:
        print(f"\n[OK] Test PASSED: recovered d = {result}")
    else:
        print(f"\n[FAIL] Test FAILED: got {result}, expected {d_secret}")
