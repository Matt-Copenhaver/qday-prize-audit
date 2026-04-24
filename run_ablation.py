#!/usr/bin/env python3
"""
Ablation Experiment: _verify() vs pure vote ranking
====================================================

Takes raw counts from a previous audited run and runs extraction two ways:

  A) WITH _verify()    — GiancarloLelli's original method (lines 282/301 of projecteleven.py)
  B) WITHOUT _verify() — identical logic, but no classical d*G==Q filter

If method A finds the correct d but method B does not rank it near the top,
then _verify() is doing all the work and there is no quantum signal.

This script uses NO quantum hardware — it operates entirely on saved counts.

Usage:
  python run_ablation.py --challenge 15 --oracle ripple --counts results/audit_15bit_*/raw_counts.json
"""

import sys
import json
import math
import datetime
import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional, Tuple, List

from projecteleven import CurveParams, EllipticCurve, ShorECDLP


def parse_bitstrings(solver, counts: Dict[str, int]) -> List[Tuple[Tuple[int, int, int], int]]:
    """Parse measurement bitstrings into (j, k, r) tuples.

    Handles both ShorECDLP (point register = n_bits) and
    RippleCarryShorECDLP (accumulator = m+1 bits) register layouts.
    Both use Qiskit MSB-left: k_bits | j_bits | point/acc_bits.
    """
    n = solver.n
    parsed = []

    # Determine point/accumulator register width
    if hasattr(solver, 'm'):
        # RippleCarryShorECDLP: accumulator is m+1 bits
        pt_width = solver.m + 1
    elif hasattr(solver, 'encoder'):
        pt_width = solver.encoder.n_bits
    else:
        pt_width = max(1, (n - 1).bit_length())

    for bitstring, count in counts.items():
        num_counting = (len(bitstring) - pt_width) // 2
        if len(bitstring) != 2 * num_counting + pt_width:
            continue

        k_bits = bitstring[:num_counting]
        j_bits = bitstring[num_counting:2 * num_counting]
        pt_bits = bitstring[2 * num_counting:]

        j = int(j_bits, 2) % n
        k = int(k_bits, 2) % n
        r = int(pt_bits, 2) % n

        parsed.append(((j, k, r), count))

    return parsed


def extract_with_verify(solver, counts: Dict[str, int]) -> Tuple[Optional[int], Dict[int, int]]:
    """GiancarloLelli's exact extraction — _verify() on every candidate.

    Uses only the direct path for RippleCarryShorECDLP (matching his code),
    and both direct + pair paths for ShorECDLP.
    """
    n = solver.n
    candidates: Dict[int, int] = {}
    parsed = parse_bitstrings(solver, counts)

    # Direct extraction
    for (j, k, r), count in parsed:
        if k == 0 or math.gcd(k, n) != 1:
            continue
        try:
            k_inv = pow(k, -1, n)
            d_cand = ((r - j) * k_inv) % n
            if solver._verify(d_cand):
                candidates[d_cand] = candidates.get(d_cand, 0) + count
        except (ValueError, ZeroDivisionError):
            pass

    # Pair-based extraction (ShorECDLP only — RippleCarry doesn't use it)
    if not hasattr(solver, 'm'):
        by_r: Dict[int, List] = {}
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
                        if solver._verify(d_cand):
                            candidates[d_cand] = candidates.get(d_cand, 0) + c1 + c2
                    except (ValueError, ZeroDivisionError):
                        pass

    if candidates:
        return max(candidates, key=candidates.get), candidates
    return None, {}


def extract_without_verify(solver, counts: Dict[str, int]) -> Tuple[Optional[int], Dict[int, int]]:
    """Same extraction logic, NO _verify(). All candidates ranked by raw votes."""
    n = solver.n
    candidates: Dict[int, int] = {}
    parsed = parse_bitstrings(solver, counts)

    # Direct extraction
    for (j, k, r), count in parsed:
        if k == 0 or math.gcd(k, n) != 1:
            continue
        try:
            k_inv = pow(k, -1, n)
            d_cand = ((r - j) * k_inv) % n
            if d_cand > 0:
                candidates[d_cand] = candidates.get(d_cand, 0) + count
        except (ValueError, ZeroDivisionError):
            pass

    # Pair-based extraction (ShorECDLP only)
    if not hasattr(solver, 'm'):
        by_r: Dict[int, List] = {}
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
                        if d_cand > 0:
                            candidates[d_cand] = candidates.get(d_cand, 0) + c1 + c2
                    except (ValueError, ZeroDivisionError):
                        pass

    if candidates:
        return max(candidates, key=candidates.get), candidates
    return None, {}


def main():
    parser = argparse.ArgumentParser(description="Ablation: _verify() vs pure vote ranking")
    parser.add_argument("--challenge", type=int, required=True)
    parser.add_argument("--oracle", choices=["dense", "permutation", "coordinate",
                                              "arithmetic", "google", "ripple"],
                        default=None)
    parser.add_argument("--counts", type=str, required=True,
                        help="Path to raw_counts.json from a previous audited run")
    parser.add_argument("--run-id", type=str, default=None)
    args = parser.parse_args()

    # Load challenge curve
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

    # Build solver (needed for _parse_bitstring and _verify)
    params = CurveParams(p, a, b, n)

    if args.oracle == "ripple":
        from ripple_carry_shor import RippleCarryShorECDLP
        solver = RippleCarryShorECDLP(params, G, Q)
    elif args.oracle == "google":
        from google_semiclassical import SemiclassicalShorECDLP
        solver = SemiclassicalShorECDLP(params, G, Q)
    elif args.oracle == "arithmetic":
        from quantum_oracle import ArithmeticShorECDLP
        solver = ArithmeticShorECDLP(params, G, Q)
    elif args.oracle == "coordinate":
        from quantum_oracle import QuantumOracleShorECDLP
        solver = QuantumOracleShorECDLP(params, G, Q)
    elif args.oracle == "permutation" or (args.oracle is None and params.n_bits > 6):
        from quantum_arithmetic import ScalableShorECDLP
        solver = ScalableShorECDLP(params, G, Q)
    else:
        solver = ShorECDLP(params, G, Q)

    # Load counts
    print(f"Loading counts from: {args.counts}")
    with open(args.counts) as f:
        counts = json.load(f)
    print(f"Unique outcomes: {len(counts)}")
    print(f"Total shots: {sum(counts.values())}")

    # ── Extraction A: WITH _verify() ──
    print("\n" + "=" * 70)
    print("EXTRACTION A: WITH _verify() (GiancarloLelli's original method)")
    print("=" * 70)

    result_v, cands_v = extract_with_verify(solver, counts)

    if result_v is not None:
        print(f"Result: d = {result_v}")
        if true_d is not None:
            print(f"Correct: {result_v == true_d}")
    else:
        print("Result: FAILED (no verified candidates)")
    print(f"Verified candidates found: {len(cands_v)}")
    if cands_v:
        sorted_v = sorted(cands_v.items(), key=lambda x: -x[1])
        for d_cand, votes in sorted_v[:10]:
            marker = " <-- true d" if d_cand == true_d else ""
            print(f"  d={d_cand}: {votes} votes{marker}")

    # ── Extraction B: WITHOUT _verify() ──
    print("\n" + "=" * 70)
    print("EXTRACTION B: WITHOUT _verify() (pure vote ranking)")
    print("=" * 70)

    result_nv, cands_nv = extract_without_verify(solver, counts)

    print(f"Total distinct candidates: {len(cands_nv)}")

    true_d_rank = None
    true_d_votes = 0
    top1_d = None
    top1_votes = 0

    if cands_nv:
        sorted_nv = sorted(cands_nv.items(), key=lambda x: -x[1])
        top1_d, top1_votes = sorted_nv[0]

        for i, (d_cand, votes) in enumerate(sorted_nv):
            if d_cand == true_d:
                true_d_rank = i + 1
                true_d_votes = votes
                break

        total_votes = sum(cands_nv.values())
        expected = total_votes / n
        z = (true_d_votes - expected) / math.sqrt(expected) if expected > 0 else 0

        print(f"\nTop-1 candidate: d = {top1_d} ({top1_votes} votes)")
        if true_d is not None:
            print(f"Top-1 is correct d: {top1_d == true_d}")

        if true_d is not None:
            print(f"\nTrue d = {true_d}:")
            print(f"  Rank: {true_d_rank} of {len(cands_nv)}")
            print(f"  Votes: {true_d_votes}")
            print(f"  Expected votes under noise: {expected:.1f}")
            print(f"  Z-score: {z:+.2f}")

        print(f"\nTop 20 candidates:")
        for i, (d_cand, votes) in enumerate(sorted_nv[:20]):
            marker = " <-- true d" if d_cand == true_d else ""
            print(f"  #{i+1}: d={d_cand}, votes={votes}{marker}")

    # ── Conclusion ──
    print("\n" + "=" * 70)
    print("CONCLUSION")
    print("=" * 70)

    if result_v is not None and result_v == true_d:
        print(f"With _verify():    d = {true_d} FOUND (correct)")
    else:
        print(f"With _verify():    FAILED")

    if true_d_rank is not None:
        print(f"Without _verify(): d = {true_d} at rank {true_d_rank} of {len(cands_nv)}")
        total_votes = sum(cands_nv.values())
        expected = total_votes / n
        z = (true_d_votes - expected) / math.sqrt(expected) if expected > 0 else 0
        print(f"                   z-score = {z:+.2f}")
        if z < 3:
            print(f"\n>>> NO QUANTUM SIGNAL: true d does not rank above noise (z < 3)")
            print(f">>> _verify() is doing all the work — this is classical brute-force")
    else:
        print(f"Without _verify(): d = {true_d} NOT FOUND in candidates")

    # ── Save results ──
    timestamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%d_%H%M%S")
    run_id = args.run_id or f"ablation_{args.challenge}bit_{timestamp}"
    results_dir = Path("results") / run_id
    results_dir.mkdir(parents=True, exist_ok=True)

    total_votes = sum(cands_nv.values()) if cands_nv else 0
    expected = total_votes / n if n > 0 else 0
    z = (true_d_votes - expected) / math.sqrt(expected) if expected > 0 else 0

    results_data = {
        "run_id": run_id,
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "source_counts": args.counts,
        "challenge": {
            "bit_length": args.challenge,
            "n": n,
            "true_d": true_d,
        },
        "with_verify": {
            "result_d": result_v,
            "correct": result_v == true_d if true_d is not None else None,
            "num_candidates": len(cands_v),
        },
        "without_verify": {
            "top1_d": top1_d,
            "top1_votes": top1_votes,
            "top1_correct": top1_d == true_d if true_d is not None else None,
            "true_d_rank": true_d_rank,
            "true_d_votes": true_d_votes,
            "num_candidates": len(cands_nv),
            "expected_votes_noise": expected,
            "z_score": z,
            "top_20": [{"d": d_cand, "votes": votes}
                       for d_cand, votes in sorted(cands_nv.items(), key=lambda x: -x[1])[:20]],
        },
    }

    results_path = results_dir / "results.json"
    with open(results_path, "w") as f:
        json.dump(results_data, f, indent=2)
    print(f"\nResults saved to: {results_path}")


if __name__ == "__main__":
    main()
