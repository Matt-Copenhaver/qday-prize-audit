"""
ECDLP Quantum Solver - Shor's Algorithm
========================================

Implements Shor's algorithm for the Elliptic Curve Discrete Logarithm Problem
without classical pre-computation of the discrete logarithm.

For Q-Day Prize Challenge: https://www.qdayprize.org/

Algorithm:
----------
1. Prepare: |j⟩|k⟩|O⟩ in superposition
2. Compute: |j⟩|k⟩|jG + kQ⟩ using quantum point addition
3. Measure point register → collapses to R
4. Apply inverse QFT to j, k registers
5. Measure j, k → extract d from j + kd ≡ r (mod n)
"""

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.circuit.library import QFTGate
import numpy as np
import math
from typing import Tuple, Optional, List, Dict
from dataclasses import dataclass


@dataclass(frozen=True)
class CurveParams:
    """Elliptic curve parameters: y² = x³ + ax + b (mod p)"""
    p: int
    a: int
    b: int
    n: int  # Group order
    
    @property
    def n_bits(self) -> int:
        return max(1, (self.n - 1).bit_length())


class EllipticCurve:
    """Elliptic curve arithmetic"""
    
    def __init__(self, params: CurveParams):
        self.p = params.p
        self.a = params.a
        self.b = params.b
        self.n = params.n
    
    def add(self, P: Optional[Tuple[int, int]], 
            Q: Optional[Tuple[int, int]]) -> Optional[Tuple[int, int]]:
        """Add two points on the curve"""
        if P is None:
            return Q
        if Q is None:
            return P
        
        x1, y1 = P
        x2, y2 = Q
        
        if x1 == x2:
            if (y1 + y2) % self.p == 0:
                return None
            num = (3 * x1 * x1 + self.a) % self.p
            denom = (2 * y1) % self.p
        else:
            num = (y2 - y1) % self.p
            denom = (x2 - x1) % self.p
        
        if denom == 0:
            return None
        
        lam = (num * pow(denom, -1, self.p)) % self.p
        x3 = (lam * lam - x1 - x2) % self.p
        y3 = (lam * (x1 - x3) - y1) % self.p
        return (x3, y3)
    
    def scalar_mult(self, k: int, P: Optional[Tuple[int, int]]) -> Optional[Tuple[int, int]]:
        """Compute k*P using double-and-add"""
        if k == 0 or P is None:
            return None
        k = k % self.n
        if k == 0:
            return None
        
        result = None
        addend = P
        while k:
            if k & 1:
                result = self.add(result, addend)
            addend = self.add(addend, addend)
            k >>= 1
        return result
    
    def enumerate_group(self, G: Tuple[int, int]) -> List[Optional[Tuple[int, int]]]:
        """Return [O, G, 2G, ..., (n-1)G]"""
        elements = [None]
        current = G
        for _ in range(1, self.n):
            elements.append(current)
            current = self.add(current, G)
        return elements


class PointEncoder:
    """Encode EC points as quantum basis state indices"""
    
    def __init__(self, ec: EllipticCurve, G: Tuple[int, int]):
        self.ec = ec
        self.n = ec.n
        self.n_bits = max(1, (self.n - 1).bit_length())
        self.elements = ec.enumerate_group(G)
        self.point_to_index = {pt: i for i, pt in enumerate(self.elements)}
    
    def encode(self, P: Optional[Tuple[int, int]]) -> int:
        return self.point_to_index.get(P, 0)
    
    def decode(self, index: int) -> Optional[Tuple[int, int]]:
        if 0 <= index < len(self.elements):
            return self.elements[index]
        return None


class QuantumPointAdder:
    """Quantum gates for EC point addition"""
    
    def __init__(self, encoder: PointEncoder):
        self.encoder = encoder
        self.n = encoder.n
        self.n_bits = encoder.n_bits
        self.ec = encoder.ec
        self._gate_cache: Dict[int, any] = {}
    
    def _build_controlled_unitary(self, S: Optional[Tuple[int, int]]) -> np.ndarray:
        """
        Build controlled unitary matrix directly.
        
        |0⟩|P⟩ → |0⟩|P⟩
        |1⟩|P⟩ → |1⟩|P + S⟩
        """
        s_idx = self.encoder.encode(S)
        
        if s_idx in self._gate_cache:
            return self._gate_cache[s_idx]
        
        n_bits = self.n_bits
        pt_dim = 2 ** n_bits
        total_dim = 2 * pt_dim  # control + point register
        
        # Build permutation for "add S"
        perm = list(range(pt_dim))
        for i in range(self.n):
            P = self.encoder.decode(i)
            P_plus_S = self.ec.add(P, S)
            perm[i] = self.encoder.encode(P_plus_S)
        
        # Build controlled unitary
        # |0⟩⊗|x⟩ → |0⟩⊗|x⟩ (identity on point register)
        # |1⟩⊗|x⟩ → |1⟩⊗|perm(x)⟩ (permutation on point register)
        U = np.zeros((total_dim, total_dim), dtype=complex)
        
        # Control = 0: identity
        for i in range(pt_dim):
            U[i, i] = 1.0
        
        # Control = 1: apply permutation
        for i in range(pt_dim):
            U[pt_dim + perm[i], pt_dim + i] = 1.0
        
        self._gate_cache[s_idx] = U
        return U
    
    def apply_controlled_add(self, qc: QuantumCircuit, 
                              ctrl_qubit: int,
                              pt_qubits: List[int],
                              S: Optional[Tuple[int, int]]):
        """Apply controlled point addition to circuit"""
        if S is None:
            return  # Adding identity does nothing
        
        U = self._build_controlled_unitary(S)
        
        # Apply as unitary on [point_register] + [control]
        # Control is MSB in the matrix (index = ctrl*pt_dim + pt_value),
        # so it must be the last qubit (Qiskit treats qubits[0] as LSB).
        qubits = pt_qubits + [ctrl_qubit]
        qc.unitary(U, qubits, label=f"CAdd{self.encoder.encode(S)}")


class ShorECDLP:
    """Shor's algorithm for ECDLP"""
    
    def __init__(self, params: CurveParams, G: Tuple[int, int], Q: Tuple[int, int]):
        self.params = params
        self.G = G
        self.Q = Q
        self.n = params.n
        
        self.ec = EllipticCurve(params)
        self.encoder = PointEncoder(self.ec, G)
        self.adder = QuantumPointAdder(self.encoder)
        
        # Validate Q is in the group
        if Q is not None and Q not in self.encoder.point_to_index:
            raise ValueError("Q is not in the group generated by G")
    
    def build_circuit(self, num_counting: int = None) -> QuantumCircuit:
        """
        Build Shor ECDLP circuit.
        
        The circuit uses G and Q as points - no knowledge of d.
        """
        n = self.n
        n_bits = self.encoder.n_bits
        
        if num_counting is None:
            num_counting = n_bits + 1
        
        j_reg = QuantumRegister(num_counting, 'j')
        k_reg = QuantumRegister(num_counting, 'k')
        pt_reg = QuantumRegister(n_bits, 'pt')
        
        cl = ClassicalRegister(2 * num_counting + n_bits, 'cr')

        qc = QuantumCircuit(j_reg, k_reg, pt_reg, cl,
                            name=f"Shor_ECDLP_n{n}")
        
        # Hadamard on counting registers
        for i in range(num_counting):
            qc.h(j_reg[i])
            qc.h(k_reg[i])
        
        # Controlled additions of 2^i * G
        G_power = self.G
        for i in range(num_counting):
            self.adder.apply_controlled_add(qc, j_reg[i], list(pt_reg), G_power)
            G_power = self.ec.add(G_power, G_power)
            if G_power is None:
                G_power = self.G
        
        # Controlled additions of 2^i * Q
        Q_power = self.Q
        for i in range(num_counting):
            if Q_power is not None:
                self.adder.apply_controlled_add(qc, k_reg[i], list(pt_reg), Q_power)
            Q_power = self.ec.add(Q_power, Q_power)
            if Q_power is None and self.Q is not None:
                Q_power = self.Q
        
        # Measure point register first (bits 0..n_bits-1)
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
    
    def extract_discrete_log(self, counts: Dict[str, int]) -> Optional[int]:
        """Extract d from measurement results"""
        n = self.n
        candidates: Dict[int, int] = {}
        
        parsed = []
        for bitstring, count in counts.items():
            result = self._parse_bitstring(bitstring)
            if result:
                parsed.append((result, count))
        
        # Direct extraction: d = (r - j) * k^(-1) mod n
        for (j, k, r), count in parsed:
            if k == 0 or math.gcd(k, n) != 1:
                continue
            try:
                k_inv = pow(k, -1, n)
                d_cand = ((r - j) * k_inv) % n
                if self._verify(d_cand):
                    candidates[d_cand] = candidates.get(d_cand, 0) + count
            except ValueError:
                pass
        
        # Pair-based extraction
        by_r: Dict[int, List[Tuple[int, int, int]]] = {}
        for (j, k, r), count in parsed:
            by_r.setdefault(r, []).append((j, k, count))
        
        for r, meas in by_r.items():
            for i, (j1, k1, c1) in enumerate(meas):
                for j2, k2, c2 in meas[i+1:]:
                    dk = (k1 - k2) % n
                    dj = (j2 - j1) % n
                    if dk == 0 or math.gcd(dk, n) != 1:
                        continue
                    try:
                        d_cand = (dj * pow(dk, -1, n)) % n
                        if self._verify(d_cand):
                            candidates[d_cand] = candidates.get(d_cand, 0) + c1 + c2
                    except ValueError:
                        pass
        
        if candidates:
            return max(candidates, key=candidates.get)
        return None
    
    def _parse_bitstring(self, bitstring: str) -> Optional[Tuple[int, int, int]]:
        """Parse measurement bitstring to (j, k, r).

        Single classical register layout (cr), Qiskit MSB-left convention:
          bitstring = k_bits | j_bits | pt_bits
        where each segment is MSB-first, so int(segment, 2) gives the value.
        """
        n_bits = self.encoder.n_bits
        num_counting = (len(bitstring) - n_bits) // 2

        if len(bitstring) != 2 * num_counting + n_bits:
            return None

        k_bits = bitstring[:num_counting]
        j_bits = bitstring[num_counting:2 * num_counting]
        pt_bits = bitstring[2 * num_counting:]

        j = int(j_bits, 2) % self.n
        k = int(k_bits, 2) % self.n
        r = int(pt_bits, 2) % self.n

        return (j, k, r)
    
    def _verify(self, d: int) -> bool:
        """Verify d is correct"""
        if d == 0:
            return self.Q is None
        if d < 0 or d >= self.n:
            return False
        return self.ec.scalar_mult(d, self.G) == self.Q
    
    def qubit_count(self) -> Dict[str, int]:
        n_bits = self.encoder.n_bits
        num_counting = n_bits + 1
        return {
            "j_register": num_counting,
            "k_register": num_counting,
            "point_register": n_bits,
            "total": 2 * num_counting + n_bits,
            "group_order": self.n,
        }


def solve_ecdlp(
    p: int, a: int, b: int, n: int,
    G: Tuple[int, int], Q: Tuple[int, int],
    shots: int = 8192,
    backend_name: str = "ibm_marrakesh",
    token: Optional[str] = None,
    instance: str = "open-instance",
    verbose: bool = True,
    oracle: Optional[str] = None,
    optimization_level: int = 3,
) -> Optional[int]:
    """Solve ECDLP using Shor's algorithm on IBM Quantum hardware."""
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2

    if verbose:
        print("=" * 60)
        print("SHOR'S ALGORITHM FOR ECDLP")
        print("Q-Day Prize: https://www.qdayprize.org/")
        print("=" * 60)

    params = CurveParams(p, a, b, n)
    n_bits = max(1, (n - 1).bit_length())

    # Strategy selection: explicit --oracle flag overrides auto-selection.
    if oracle == "ripple":
        from ripple_carry_shor import RippleCarryShorECDLP
        solver = RippleCarryShorECDLP(params, G, Q)
        strategy = "ripple-carry modular addition (CDKM)"
    elif oracle == "google":
        from google_semiclassical import SemiclassicalShorECDLP
        solver = SemiclassicalShorECDLP(params, G, Q)
        strategy = "Google semiclassical PE (qubit-recycled)"
    elif oracle == "arithmetic":
        from quantum_oracle import ArithmeticShorECDLP
        solver = ArithmeticShorECDLP(params, G, Q)
        strategy = "arithmetic (coordinate + QFT primitives)"
    elif oracle == "coordinate":
        from quantum_oracle import QuantumOracleShorECDLP
        solver = QuantumOracleShorECDLP(params, G, Q)
        strategy = "coordinate oracle"
    elif oracle == "permutation" or (oracle is None and n_bits > 6):
        from quantum_arithmetic import ScalableShorECDLP
        solver = ScalableShorECDLP(params, G, Q)
        strategy = "efficient permutation"
    else:
        solver = ShorECDLP(params, G, Q)
        strategy = "dense unitary"

    if verbose:
        print(f"\nCurve: y^2 = x^3 + {a}x + {b} (mod {p})")
        print(f"Group order: n = {n}")
        print(f"Generator: G = {G}")
        print(f"Target: Q = {Q}")
        print(f"Strategy: {strategy}")
        reqs = solver.qubit_count()
        print(f"\nQubits: {reqs['total']} total")

    qc = solver.build_circuit()

    if verbose:
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

    qc_t = transpile(qc, backend, optimization_level=optimization_level)

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

    print(f"\nUnique outcomes: {len(counts)}")
    d = solver.extract_discrete_log(counts)

    if verbose:
        print(f"\n{'=' * 60}")
        if d is not None:
            print(f"RESULT: d = {d}")
            computed = solver.ec.scalar_mult(d, G)
            print(f"Verification: {d}*G = {computed}")
            if computed == Q:
                print("[OK] VERIFIED")
            else:
                print("[FAIL] MISMATCH")
        else:
            print("FAILED: Could not extract discrete log")
        print("=" * 60)

    return d


# Test curves from qday challenge input
CURVES = {
    "curve_4": {
        "p": 13,
        "a": 0,
        "b": 7,
        "n": 7,
        "G": (11, 5),
        "Q": (11, 8)
    }
}


def verify_curve(curve: dict) -> bool:
    """Verify curve parameters"""
    params = CurveParams(curve["p"], curve["a"], curve["b"], curve["n"])
    ec = EllipticCurve(params)
    G = curve["G"]
    
    # Check G is on curve
    x, y = G
    left = (y * y) % params.p
    right = (x**3 + params.a * x + params.b) % params.p
    if left != right:
        print(f"G = {G} not on curve")
        return False
    
    # Check order
    test = ec.scalar_mult(params.n, G)
    if test is not None:
        print(f"{params.n}*G = {test}, expected O")
        return False
    
    return True


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description="ECDLP Solver - Shor's Algorithm")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--curve", choices=list(CURVES.keys()),
                       help="Use a built-in test curve")
    group.add_argument("--challenge", type=int, metavar="BIT_LENGTH",
                       help="Solve challenge curve from input_curves.json by bit length")
    parser.add_argument("--d", type=int,
                        help="Secret key for testing (ignored with --challenge)")
    parser.add_argument("--shots", type=int, default=8192)
    parser.add_argument("--backend", default="ibm_marrakesh")
    parser.add_argument("--token", type=str,
                        help="IBM Quantum API token. Saves to local account on first use.")
    parser.add_argument("--instance", type=str, default="open-instance",
                        help="IBM Quantum instance (default: open-instance)")
    parser.add_argument("--oracle", choices=["dense", "permutation", "coordinate", "arithmetic", "google", "ripple"],
                        default=None,
                        help="Oracle strategy (default: auto-select based on curve size)")
    parser.add_argument("--optimization-level", type=int, default=3, choices=[0, 1, 2, 3],
                        help="Qiskit transpilation optimization level (default: 3)")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()

    if args.token:
        from qiskit_ibm_runtime import QiskitRuntimeService
        QiskitRuntimeService.save_account(
            channel="ibm_quantum_platform", token=args.token, overwrite=True
        )
        print(f"IBM Quantum account saved.")

    if args.challenge:
        with open("input_curves.json") as f:
            all_curves = json.load(f)
        curve_data = next(
            (c for c in all_curves if c["bit_length"] == args.challenge), None
        )
        if curve_data is None:
            print(f"No challenge curve found for {args.challenge}-bit")
            exit(1)

        p = curve_data["prime"]
        a, b = 0, 7
        n = curve_data["subgroup_order"]
        G = tuple(curve_data["generator_point"])
        Q = tuple(curve_data["public_key"])
        d_expected = curve_data.get("private_key")

        print(f"\nChallenge curve: {args.challenge}-bit")
        print(f"p = {p}, n = {n}")
        print(f"G = {G}")
        print(f"Q = {Q}")
    else:
        curve = CURVES[args.curve].copy()
        p, a, b, n = curve["p"], curve["a"], curve["b"], curve["n"]
        G = curve["G"]

        params = CurveParams(p, a, b, n)
        ec = EllipticCurve(params)
        d_secret = (args.d or 6) % n or 1
        Q = ec.scalar_mult(d_secret, G)
        d_expected = d_secret

        print(f"\nCurve: {args.curve}")
        print(f"Secret: d = {d_secret}, Q = {Q}")

    if not verify_curve({"p": p, "a": a, "b": b, "n": n, "G": G}):
        print("Curve verification failed!")
        exit(1)

    if args.verify_only:
        exit(0)

    # Empty line for formatting
    print(f"")

    result = solve_ecdlp(
        p, a, b, n, G, Q,
        shots=args.shots,
        backend_name=args.backend,
        token=args.token,
        instance=args.instance,
        oracle=args.oracle,
        optimization_level=args.optimization_level,
    )

    if result is not None and d_expected is not None:
        if result == d_expected:
            print("\n[OK] SUCCESS: Recovered correct secret key")
        else:
            print(f"\n~ Found d = {result} (expected {d_expected})")
    elif result is None:
        print("\n[FAIL] FAILED")