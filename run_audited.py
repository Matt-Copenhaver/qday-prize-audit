#!/usr/bin/env python3
"""
Audited Runner for GiancarloLelli's Q-Day Prize Code
=====================================================

Wraps his exact solve_ecdlp() flow but captures artifacts that his code
only prints to stdout:

  - Raw measurement counts (JSON)
  - Transpiled circuit metrics (depth, 2Q gates, qubit count, gate breakdown)
  - Job ID, backend name, timestamps
  - Extraction result and verification
  - Source code hash for provenance

All results are saved to results/<run_id>/ with full provenance.

Usage:
  python run_audited.py --challenge 16 --oracle ripple --backend ibm_fez --shots 20000
  python run_audited.py --challenge 4 --backend ibm_torino --shots 8192

To replay extraction on previously saved counts (no hardware needed):
  python run_audited.py --challenge 16 --oracle ripple --from-counts results/<run_id>/raw_counts.json
"""

import sys
import os
import json
import math
import datetime
import hashlib
import argparse
from pathlib import Path
from typing import Dict, Optional, Tuple

from projecteleven import (
    CurveParams, EllipticCurve, ShorECDLP,
    verify_curve
)


def main():
    parser = argparse.ArgumentParser(description="Audited Q-Day Prize Run")
    parser.add_argument("--challenge", type=int, required=True,
                        help="Challenge bit length (e.g., 4, 6, 10, 15, 16, 17)")
    parser.add_argument("--oracle", choices=["dense", "permutation", "coordinate",
                                              "arithmetic", "google", "ripple"],
                        default=None, help="Oracle strategy (same as projecteleven.py)")
    parser.add_argument("--backend", default="ibm_fez", help="IBM Quantum backend")
    parser.add_argument("--shots", type=int, default=8192)
    parser.add_argument("--token", type=str, default=None)
    parser.add_argument("--instance", type=str, default="open-instance")
    parser.add_argument("--optimization-level", type=int, default=3)
    parser.add_argument("--run-id", type=str, default=None,
                        help="Custom run ID (default: auto-generated)")
    parser.add_argument("--from-counts", type=str, default=None,
                        help="Skip hardware; load raw counts from JSON file")
    args = parser.parse_args()

    # ── Load challenge curve ──
    with open("input_curves.json") as f:
        all_curves = json.load(f)
    curve_data = next(
        (c for c in all_curves if c["bit_length"] == args.challenge), None
    )
    if curve_data is None:
        print(f"No challenge curve found for {args.challenge}-bit")
        sys.exit(1)

    p = curve_data["prime"]
    a, b = 0, 7
    n = curve_data["subgroup_order"]
    G = tuple(curve_data["generator_point"])
    Q = tuple(curve_data["public_key"])
    true_d = curve_data.get("private_key")

    # ── Create results directory ──
    timestamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%d_%H%M%S")
    run_id = args.run_id or f"audit_{args.challenge}bit_{timestamp}"
    results_dir = Path("results") / run_id
    results_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print(f"AUDITED Q-DAY PRIZE RUN: {args.challenge}-bit challenge")
    print(f"Run ID: {run_id}")
    print("=" * 80)
    print(f"\nCurve: y^2 = x^3 + 0x + 7 (mod {p})")
    print(f"Group order: n = {n}")
    print(f"Generator: G = {G}")
    print(f"Target: Q = {Q}")
    if true_d is not None:
        print(f"Known d: {true_d}")

    # ── Build solver (his exact code path) ──
    params = CurveParams(p, a, b, n)

    if args.oracle == "ripple":
        from ripple_carry_shor import RippleCarryShorECDLP
        solver = RippleCarryShorECDLP(params, G, Q)
        strategy = "ripple-carry modular addition (CDKM)"
    elif args.oracle == "google":
        from google_semiclassical import SemiclassicalShorECDLP
        solver = SemiclassicalShorECDLP(params, G, Q)
        strategy = "Google semiclassical PE"
    elif args.oracle == "arithmetic":
        from quantum_oracle import ArithmeticShorECDLP
        solver = ArithmeticShorECDLP(params, G, Q)
        strategy = "arithmetic"
    elif args.oracle == "coordinate":
        from quantum_oracle import QuantumOracleShorECDLP
        solver = QuantumOracleShorECDLP(params, G, Q)
        strategy = "coordinate oracle"
    elif args.oracle == "permutation" or (args.oracle is None and params.n_bits > 6):
        from quantum_arithmetic import ScalableShorECDLP
        solver = ScalableShorECDLP(params, G, Q)
        strategy = "efficient permutation"
    else:
        solver = ShorECDLP(params, G, Q)
        strategy = "dense unitary"

    print(f"Strategy: {strategy}")
    reqs = solver.qubit_count()
    print(f"Qubits: {reqs['total']} total")

    # ── Get counts: either from hardware or from file ──
    if args.from_counts:
        print(f"\nLoading counts from: {args.from_counts}")
        with open(args.from_counts) as f:
            counts = json.load(f)
        job_id = f"from_file:{args.from_counts}"
        backend_name = "from_file"
        transpiled_depth = None
        two_qubit_gates = None
        total_qubits_transpiled = None
        circuit_ops = {}
    else:
        # Build and transpile circuit (his exact flow)
        qc = solver.build_circuit()
        pre_transpile_depth = qc.depth()
        print(f"Pre-transpile depth: {pre_transpile_depth}")

        from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
        from qiskit import transpile

        if args.token:
            service = QiskitRuntimeService(
                channel="ibm_quantum_platform", token=args.token, overwrite=True
            )
        else:
            service = QiskitRuntimeService(
                channel="ibm_quantum_platform", instance=args.instance
            )
        backend = service.backend(args.backend)
        backend_name = backend.name

        print(f"\nBackend: {backend_name}")

        qc_t = transpile(qc, backend, optimization_level=args.optimization_level)
        transpiled_depth = qc_t.depth()
        total_qubits_transpiled = qc_t.num_qubits
        circuit_ops = dict(qc_t.count_ops())
        two_qubit_gates = sum(v for k, v in circuit_ops.items()
                              if k in ['cx', 'ecr', 'cz'])

        print(f"Transpiled depth: {transpiled_depth}")
        print(f"Two-qubit gates: {two_qubit_gates}")
        print(f"Total qubits (transpiled): {total_qubits_transpiled}")

        # Run on hardware
        sampler = SamplerV2(mode=backend)
        job = sampler.run([qc_t], shots=args.shots)
        job_id = job.job_id()

        print(f"\nJob ID: {job_id}")
        print("Waiting for results...")

        result = job.result()
        pub_result = result[0]
        counts = pub_result.data.cr.get_counts()

    print(f"\nUnique outcomes: {len(counts)}")
    print(f"Total shots: {sum(counts.values())}")

    # ── Save raw counts ──
    counts_path = results_dir / "raw_counts.json"
    with open(counts_path, "w") as f:
        json.dump(counts, f, indent=2)
    print(f"Raw counts saved to: {counts_path}")

    # ── Run his extraction (unmodified) ──
    d_result = solver.extract_discrete_log(counts)

    print(f"\n{'=' * 60}")
    if d_result is not None:
        print(f"RESULT: d = {d_result}")
        computed = solver.ec.scalar_mult(d_result, G)
        print(f"Verification: {d_result}*G = {computed}")
        print(f"Match Q: {computed == Q}")
        if true_d is not None:
            print(f"Matches known d={true_d}: {d_result == true_d}")
    else:
        print("RESULT: FAILED — could not extract discrete log")
    print("=" * 60)

    # ── Compute fidelity estimate ──
    fidelity = None
    if two_qubit_gates is not None and two_qubit_gates > 0:
        fidelity = math.exp(-two_qubit_gates * 0.0005)

    # ── Save structured results ──
    results_data = {
        "run_id": run_id,
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "challenge": {
            "bit_length": args.challenge,
            "p": p, "a": a, "b": b, "n": n,
            "G": list(G), "Q": list(Q),
            "true_d": true_d,
        },
        "hardware": {
            "backend": backend_name,
            "job_id": job_id,
            "shots": args.shots,
            "optimization_level": args.optimization_level,
        },
        "circuit": {
            "oracle_strategy": strategy,
            "logical_qubits": reqs["total"],
            "transpiled_depth": transpiled_depth,
            "two_qubit_gates": two_qubit_gates,
            "transpiled_qubits": total_qubits_transpiled,
            "gate_ops": circuit_ops,
            "estimated_fidelity": fidelity,
        },
        "measurement": {
            "unique_outcomes": len(counts),
            "total_shots": sum(counts.values()),
            "counts_file": "raw_counts.json",
        },
        "extraction": {
            "method": "extract_discrete_log (original, with _verify)",
            "result_d": d_result,
            "correct": d_result == true_d if true_d is not None else None,
        },
        "provenance": {
            "source_code_hash_sha256": hashlib.sha256(
                Path("projecteleven.py").read_bytes()
            ).hexdigest(),
            "audit_script": "run_audited.py",
        },
    }

    results_path = results_dir / "results.json"
    with open(results_path, "w") as f:
        json.dump(results_data, f, indent=2)
    print(f"\nFull results saved to: {results_path}")


if __name__ == "__main__":
    main()
