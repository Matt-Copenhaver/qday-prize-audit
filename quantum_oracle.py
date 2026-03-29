"""
Quantum Oracle for ECDLP — Coordinate-Based Encoding
=====================================================

Implements a full quantum oracle for Shor's ECDLP algorithm using
coordinate-based point encoding.  Instead of mapping EC points to
group indices (0..n-1), the quantum register holds actual (x, y)
field-element values in binary plus an identity flag.

Each controlled "add classical point S" operation is realised as a
permutation on the coordinate register, computed from the EC addition
formula and decomposed into transpositions via the proven infrastructure
in quantum_arithmetic.py.

Targets curves up to ~6-bit group order (p <= 43).  Enabled by
``--oracle coordinate`` on the CLI.
"""

import math
from typing import Tuple, Optional, List, Dict
from dataclasses import dataclass

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.circuit.library import QFTGate

from projecteleven import CurveParams, EllipticCurve, PointEncoder, ShorECDLP
from quantum_arithmetic import _controlled_transposition


# ---------------------------------------------------------------------------
# Coordinate-based point encoder
# ---------------------------------------------------------------------------

class CoordinatePointEncoder:
    """Encode EC points as (x, y, identity_flag) in binary.

    Bit layout of the encoded integer (LSB-first registers):
        bits [0 .. f_bits-1]          = x coordinate
        bits [f_bits .. 2*f_bits-1]   = y coordinate
        bit  [2*f_bits]               = identity flag (1 = point at infinity)

    Total width: 2*f_bits + 1  qubits.
    """

    def __init__(self, params: CurveParams):
        self.p = params.p
        self.a = params.a
        self.b = params.b
        self.f_bits = params.p.bit_length()
        self.pt_bits = 2 * self.f_bits + 1

    def encode(self, P: Optional[Tuple[int, int]]) -> int:
        """Encode a curve point to a basis-state index."""
        if P is None:
            return 1 << (2 * self.f_bits)          # id_flag = 1
        x, y = P
        return (y << self.f_bits) | x               # id_flag = 0

    def decode(self, val: int) -> Optional[Tuple[int, int]]:
        """Decode a basis-state index to a curve point (or None for O)."""
        id_flag = (val >> (2 * self.f_bits)) & 1
        if id_flag:
            return None
        mask = (1 << self.f_bits) - 1
        x = val & mask
        y = (val >> self.f_bits) & mask
        return (x, y)

    def is_valid(self, val: int) -> bool:
        """Check whether *val* encodes a valid curve point or identity."""
        P = self.decode(val)
        if P is None:
            # Valid only if it is the canonical identity encoding
            return val == self.encode(None)
        x, y = P
        if x >= self.p or y >= self.p:
            return False
        lhs = (y * y) % self.p
        rhs = (x * x * x + self.a * x + self.b) % self.p
        return lhs == rhs


# ---------------------------------------------------------------------------
# Coordinate-based controlled point adder
# ---------------------------------------------------------------------------

class CoordinatePointAdder:
    """Controlled addition of a classical point S to a quantum coordinate register.

    Builds the "add S" permutation over the 2*f_bits+1 qubit point register,
    then decomposes it into transpositions using ``_controlled_transposition``
    from quantum_arithmetic.py.  Each controlled addition is encapsulated as
    an opaque sub-circuit gate to keep Qiskit's DAG small.
    """

    def __init__(self, ec: EllipticCurve, params: CurveParams):
        self.ec = ec
        self.params = params
        self.encoder = CoordinatePointEncoder(params)
        self.f_bits = self.encoder.f_bits
        self.pt_bits = self.encoder.pt_bits
        self._perm_cache: Dict[Optional[Tuple[int, int]], List[Tuple[int, int]]] = {}
        self._circuit_cache: Dict[Tuple[Optional[Tuple[int, int]], int], QuantumCircuit] = {}

    # -- permutation construction -------------------------------------------

    def _get_transpositions(self, S: Optional[Tuple[int, int]]) -> List[Tuple[int, int]]:
        """Cycle-decompose the 'add S' permutation into transpositions."""
        if S in self._perm_cache:
            return self._perm_cache[S]

        pt_dim = 1 << self.pt_bits
        perm = list(range(pt_dim))

        for val in range(pt_dim):
            if not self.encoder.is_valid(val):
                continue  # invalid encoding -> identity (unmapped)
            P = self.encoder.decode(val)
            result = self.ec.add(P, S)
            perm[val] = self.encoder.encode(result)

        # Cycle decomposition -> transpositions
        transpositions: List[Tuple[int, int]] = []
        visited = [False] * pt_dim
        for start in range(pt_dim):
            if visited[start] or perm[start] == start:
                visited[start] = True
                continue
            cycle: List[int] = []
            j = start
            while not visited[j]:
                visited[j] = True
                cycle.append(j)
                j = perm[j]
            for idx in range(len(cycle) - 1, 0, -1):
                transpositions.append((cycle[0], cycle[idx]))

        self._perm_cache[S] = transpositions
        return transpositions

    # -- sub-circuit construction -------------------------------------------

    def _build_controlled_add_circuit(
        self, S: Optional[Tuple[int, int]], num_ancilla: int
    ) -> QuantumCircuit:
        """Build an isolated sub-circuit for controlled addition of S.

        Qubit layout inside the sub-circuit:
            [0 .. pt_bits-1]                        = point register (x|y|id)
            [pt_bits]                               = control qubit
            [pt_bits+1 .. pt_bits+num_ancilla]      = V-chain ancilla
        """
        cache_key = (S, num_ancilla)
        if cache_key in self._circuit_cache:
            return self._circuit_cache[cache_key]

        n_qubits = self.pt_bits + 1 + num_ancilla
        sub = QuantumCircuit(n_qubits, name=f"CAddCoord_{S}")

        pt_qubits = list(range(self.pt_bits))
        ctrl = self.pt_bits
        anc = list(range(self.pt_bits + 1, n_qubits)) if num_ancilla > 0 else None

        for a, b in self._get_transpositions(S):
            _controlled_transposition(sub, ctrl, a, b, pt_qubits, self.pt_bits, anc)

        self._circuit_cache[cache_key] = sub
        return sub

    # -- public interface ---------------------------------------------------

    def apply_controlled_add(
        self,
        qc: QuantumCircuit,
        ctrl_qubit: int,
        pt_qubits: List[int],
        S: Optional[Tuple[int, int]],
        ancilla_qubits: List[int] = None,
    ):
        """Append controlled coordinate-based point addition to *qc*."""
        if S is None:
            return

        num_ancilla = len(ancilla_qubits) if ancilla_qubits else 0
        sub = self._build_controlled_add_circuit(S, num_ancilla)

        qubit_map = pt_qubits + [ctrl_qubit]
        if ancilla_qubits:
            qubit_map.extend(ancilla_qubits)

        qc.append(sub.to_gate(), qubit_map)


# ---------------------------------------------------------------------------
# Shor ECDLP solver with coordinate-based oracle
# ---------------------------------------------------------------------------

class QuantumOracleShorECDLP:
    """Shor's algorithm for ECDLP using a coordinate-based quantum oracle.

    Drop-in replacement for ``ShorECDLP`` / ``ScalableShorECDLP``.
    The point register stores (x, y, id_flag) in binary instead of a
    group index.
    """

    def __init__(self, params: CurveParams,
                 G: Tuple[int, int], Q: Tuple[int, int]):
        self.params = params
        self.G = G
        self.Q = Q
        self.n = params.n

        self.ec = EllipticCurve(params)
        self.coord_encoder = CoordinatePointEncoder(params)
        self.adder = CoordinatePointAdder(self.ec, params)

        # Index-based encoder for translating measurements back to group indices
        self.index_encoder = PointEncoder(self.ec, G)

        # For extract_discrete_log reuse
        self._extraction_solver = ShorECDLP(params, G, Q)

        if Q is not None and Q not in self.index_encoder.point_to_index:
            raise ValueError("Q is not in the group generated by G")

    def build_circuit(self, num_counting: int = None) -> QuantumCircuit:
        """Build the Shor ECDLP circuit with coordinate-based oracle."""
        n = self.n
        n_bits = self.index_encoder.n_bits
        f_bits = self.coord_encoder.f_bits
        pt_bits = self.coord_encoder.pt_bits

        if num_counting is None:
            num_counting = n_bits + 1

        # V-chain MCX needs (num_controls - 2) ancillas.
        # Max controls per transposition = pt_bits (ctrl + pt_bits-1 non-pivot).
        num_ancilla = max(0, pt_bits - 2)

        j_reg = QuantumRegister(num_counting, 'j')
        k_reg = QuantumRegister(num_counting, 'k')
        x_reg = QuantumRegister(f_bits, 'x')
        y_reg = QuantumRegister(f_bits, 'y')
        id_reg = QuantumRegister(1, 'id')
        anc_reg = (QuantumRegister(num_ancilla, 'anc')
                   if num_ancilla > 0 else None)
        cl = ClassicalRegister(2 * num_counting + pt_bits, 'cr')

        regs = [j_reg, k_reg, x_reg, y_reg, id_reg]
        if anc_reg:
            regs.append(anc_reg)
        regs.append(cl)

        qc = QuantumCircuit(*regs, name=f"Shor_ECDLP_coord_n{n}")

        ancilla_qubits = list(anc_reg) if anc_reg else None

        # Point register qubit list: x | y | id  (matches encoder bit layout)
        pt_qubits = list(x_reg) + list(y_reg) + list(id_reg)

        # Initialise point register to identity O  (id_flag = 1)
        qc.x(id_reg[0])

        # Hadamard on counting registers
        for i in range(num_counting):
            qc.h(j_reg[i])
            qc.h(k_reg[i])

        # Controlled additions of 2^i * G
        G_power = self.G
        for i in range(num_counting):
            self.adder.apply_controlled_add(
                qc, j_reg[i], pt_qubits, G_power, ancilla_qubits
            )
            G_power = self.ec.add(G_power, G_power)
            if G_power is None:
                G_power = self.G

        # Controlled additions of 2^i * Q
        Q_power = self.Q
        for i in range(num_counting):
            if Q_power is not None:
                self.adder.apply_controlled_add(
                    qc, k_reg[i], pt_qubits, Q_power, ancilla_qubits
                )
            Q_power = self.ec.add(Q_power, Q_power)
            if Q_power is None and self.Q is not None:
                Q_power = self.Q

        # Measure point register (x | y | id)
        for i in range(pt_bits):
            qc.measure(pt_qubits[i], cl[i])

        # Inverse QFT on counting registers
        qc.append(QFTGate(num_counting).inverse(), j_reg)
        qc.append(QFTGate(num_counting).inverse(), k_reg)

        # Measure counting registers
        for i in range(num_counting):
            qc.measure(j_reg[i], cl[pt_bits + i])
            qc.measure(k_reg[i], cl[pt_bits + num_counting + i])

        return qc

    def extract_discrete_log(self, counts: Dict[str, int]) -> Optional[int]:
        """Extract d from measurement results.

        Translates coordinate-encoded point measurements to group indices,
        then delegates to the proven ShorECDLP extraction logic.
        """
        n_bits = self.index_encoder.n_bits
        pt_bits = self.coord_encoder.pt_bits
        num_counting = n_bits + 1

        expected_len = 2 * num_counting + pt_bits

        translated: Dict[str, int] = {}
        for bitstring, count in counts.items():
            if len(bitstring) != expected_len:
                continue

            # Qiskit MSB-left: k_bits | j_bits | pt_bits
            k_bits = bitstring[:num_counting]
            j_bits = bitstring[num_counting:2 * num_counting]
            coord_bits = bitstring[2 * num_counting:]

            # Decode coordinate measurement to a point
            coord_val = int(coord_bits, 2)
            P = self.coord_encoder.decode(coord_val)

            # Convert to group index
            r = self.index_encoder.point_to_index.get(P, None)
            if r is None:
                continue  # measurement landed on invalid encoding, skip

            # Rebuild bitstring with index-based point register for ShorECDLP
            r_bits = format(r, f'0{n_bits}b')
            new_bitstring = k_bits + j_bits + r_bits
            translated[new_bitstring] = translated.get(new_bitstring, 0) + count

        if not translated:
            return None

        return self._extraction_solver.extract_discrete_log(translated)

    def qubit_count(self) -> Dict[str, int]:
        n_bits = self.index_encoder.n_bits
        f_bits = self.coord_encoder.f_bits
        pt_bits = self.coord_encoder.pt_bits
        num_counting = n_bits + 1
        num_ancilla = max(0, pt_bits - 2)
        return {
            "j_register": num_counting,
            "k_register": num_counting,
            "x_register": f_bits,
            "y_register": f_bits,
            "id_flag": 1,
            "point_register_total": pt_bits,
            "ancilla": num_ancilla,
            "total": 2 * num_counting + pt_bits + num_ancilla,
            "group_order": self.n,
        }
