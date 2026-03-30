# Shor's Algorithm for ECDLP — Q-Day Prize Submission

Quantum solver for the Elliptic Curve Discrete Logarithm Problem (ECDLP), built for the [Q-Day Prize Challenge](https://www.qdayprize.org/) by [Project Eleven](https://www.projecteleven.com/). The goal: recover ECC private keys on real quantum hardware using Shor's algorithm.

- **Author**: Giancarlo Lelli
- **Contact**: [gcarlo.lelli@gmail.com](mailto:gcarlo.lelli@gmail.com)
- **LinkedIn**: https://www.linkedin.com/in/giancarlolelli
- **Background**: Technology leader with 10+ years in enterprise software, full-stack architecture, and cloud-native development. Background in computer science with hands-on experience across .NET, Python, Rust, and Cloud ecosystems. Currently working as Cloud GTM Specialist focused on solution architecture and sales engineering.

## Approach

All challenge curves use **y^2 = x^3 + 7** over F_p (a = 0, b = 7), matching the secp256k1 family. The solver implements the two-register variant of Shor's algorithm for ECDLP:

1. Prepare counting registers |j>, |k> in uniform superposition (Hadamard)
2. Compute |j>|k>|jG + kQ> via 2t controlled point additions (t = num_counting qubits)
3. Measure the point register, collapsing it to some group element R
4. Apply inverse QFT to the counting registers
5. Measure j, k and extract d from the relation j + kd = r (mod n)

The private key d is recovered by collecting multiple (j, k) samples that satisfy the same linear relation modulo the group order n. The solver supports four oracle strategies for the controlled point additions, selected automatically based on curve size or manually via `--oracle`.

## Oracle Strategies

### Strategy 1: Dense Unitary (default for n_bits <= 6)

Used for curves with group order up to ~6 bits. Implemented in `projecteleven.py`.

Each controlled point addition "add S" is represented as a **2^(n+1) x 2^(n+1) permutation matrix** applied via `qc.unitary()`. The matrix encodes the full group action: the upper-left block is identity (control=0), the lower-right block permutes basis states according to the map P -> P+S (control=1).

- **Encoding**: Group index (0..n-1)
- **Memory**: O(2^{2n}) per matrix
- **Qubits**: 2t + n (two counting registers + point register)
- **Limitation**: Qiskit's unitary decomposition is O(4^n), making this infeasible beyond ~6-bit

### Strategy 2: Efficient Permutation Decomposition (default for n_bits > 6)

Used for larger curves. Implemented in `quantum_arithmetic.py`.

Instead of building dense matrices, each "add S" permutation is **cycle-decomposed into transpositions**. Each transposition (swap of two basis states |a> <-> |b>) is implemented with:

1. **CNOT reduction** -- CNOTs from a pivot bit to all other differing bits, reducing the multi-bit difference to a single-bit difference
2. **Multi-controlled X** -- An MCX gate on the pivot bit, conditioned on all other bits matching the target pattern
3. **Undo CNOTs** -- Reverse step 1 to restore the non-pivot bits

The MCX uses **V-chain decomposition** with (n-2) dedicated ancilla qubits, giving O(n) Toffoli gates per MCX instead of O(n^2) without ancillas. Each controlled addition is built as an **isolated sub-circuit** and appended as a single opaque gate, avoiding quadratic DAG growth in Qiskit.

- **Encoding**: Group index (0..n-1)
- **Memory**: O(N) per addition (N = group order)
- **Qubits**: 2t + n + (n-2) ancillas
- **Gates per addition**: O(N * n)

### Strategy 3: Coordinate-Based Quantum Oracle (`--oracle coordinate`)

Available for curves up to ~6-bit. Implemented in `quantum_oracle.py`.

Instead of encoding points as group indices, the quantum register holds actual **(x, y) field-element coordinates** in binary plus an identity flag. The point register layout is:

- `x_reg`: f_bits qubits (f_bits = ceil(log2(p)))
- `y_reg`: f_bits qubits
- `id_flag`: 1 qubit (1 = point at infinity)

Each controlled "add S" is computed from the EC addition formula over all valid coordinate encodings, producing a permutation on the coordinate register. This permutation is cycle-decomposed into transpositions using the same CNOT-reduction + MCX infrastructure as Strategy 2.

- **Encoding**: (x, y, id_flag) coordinates
- **Qubits**: 2t + 2*f_bits + 1 + max(0, 2*f_bits - 1) ancillas
- **Gates per addition**: O(N * f_bits)

### Strategy 4: Arithmetic Oracle (`--oracle arithmetic`)

Framework for polynomial-scaling point addition. Implemented in `quantum_oracle.py` and `quantum_arithmetic.py`.

Uses coordinate encoding (same as Strategy 3) with QFT-based modular arithmetic primitives as building blocks toward fully arithmetic point addition. The codebase includes tested implementations of:

- **Beauregard modular adder** -- QFT-based (target + constant) mod p with proper ancilla uncompute
- **Quantum-quantum modular multiply** -- |a>|b>|0> -> |a>|b>|a*b mod p> via shift-and-add with explicit modular doubling, O(n^3) gates
- **Modular inverse permutation** -- |x> -> |x^{-1} mod p> via lookup table transpositions
- **Controlled quantum-quantum modular add** -- Controlled |a> -> |a + b mod p> with Beauregard reduction

The arithmetic primitives achieve O(n^3) scaling per point addition vs O(N*n) for the permutation approach. However, the QFT-based operations carry a ~150x larger constant factor, making the arithmetic approach more efficient only for curves above ~20-bit group order. For current challenge sizes (up to 12-bit), the permutation-based adder remains faster and is used by default.

### Comparison

| Metric | Dense Unitary | Efficient Permutation | Coordinate Oracle | Arithmetic Oracle |
|--------|--------------|----------------------|-------------------|-------------------|
| Point encoding | Group index | Group index | (x, y, id_flag) | (x, y, id_flag) |
| Scaling per addition | O(4^n) decomp. | O(N * n) | O(N * f_bits) | O(n^3) asymptotic |
| Qubits (4-bit) | 11 | 13 | 24 | 24 |
| Qubits (6-bit) | 17 | 21 | 36 | 36 |
| 2Q gates (4-bit) | 774 | ~1,200 | 6,449 | 6,449 |
| 2Q gates (6-bit) | 23,471 | ~38,000 | 95,254 | 95,254 |
| Practical range | <= 6-bit | <= ~16-bit | <= 6-bit | >= 20-bit (future) |

### QFT Arithmetic Primitives

The codebase includes QFT-based modular arithmetic building blocks (Beauregard/Draper adders, quantum-quantum modular multiplication, modular inverse/negation) as a foundation toward fully arithmetic coordinate-encoding at 256-bit. These primitives have been verified correct via Statevector simulation for primes up to p=13.

## Results

Successfully recovered private keys on IBM Quantum hardware for challenge curves up to 10-bit:

| Challenge | p | n | Strategy | Qubits | 2Q Gates | Transpiled Depth | Shots | Backend | Recovered d | Job ID |
|-----------|-----|------|----------|--------|----------|------------------|-------|-------------|-------------|--------|
| 4-bit | 13 | 7 | Dense unitary | 11 | 774 | 2,425 | 8,192 | ibm_torino | 6 | d73u28kvllmc73anvi90 |
| 4-bit | 13 | 7 | Coordinate oracle | 24 | 6,449 | 13,125 | 8,192 | ibm_kingston | 6 | d74ht798qmgc73fm32c0 |
| 4-bit | 13 | 7 | Arithmetic oracle | 24 | 6,477 | 13,452 | 8,192 | ibm_torino | 6 | d75648lbjrds73ec0eng |
| 6-bit | 43 | 31 | Dense unitary | 17 | 23,471 | 72,475 | 8,192 | ibm_torino | 18 | d73u2l5koquc73e24u8g |
| 6-bit | 43 | 31 | Coordinate oracle | 36 | 95,254 | 169,766 | 8,192 | ibm_kingston | 18 | d74hu918qmgc73fm33g0 |
| 8-bit | 163 | 139 | Efficient permutation | 32 | 294,628 | 599,517 | 8,192 | ibm_kingston | 103 | d73ui15koquc73e25e4g |
| 9-bit | 349 | 313 | Efficient permutation | 36 | 887,544 | 1,764,266 | 8,192 | ibm_torino | 135 | d73ua2h8qmgc73flei9g |
| 10-bit | 547 | 547 | Efficient permutation | 40 | 2,049,138 | 3,948,250 | 1,024 | ibm_torino | 165 | d752vfu8faus73evhovg |

All runs were executed on the IBM Quantum open-instance plan, which grants 10 minutes of free quantum computation per month. The 10-bit challenge required reducing shots to 1,024 to fit within the QPU time budget. Full execution logs are in the `executions/` folder.

> **Important**: I believe that going beyond 10-bit with our implementation is feasible, at least up until 12-bit. But due to limited resources, I wasn't able to verify this claim.

## Quick Start

```bash
git clone https://github.com/GiancarloLelli/quantum.git
cd quantum

python -m venv .
Scripts\Activate.ps1    # For Windows only

pip install -r requirements.txt
```

### How to run

You need an [IBM Quantum](https://quantum.ibm.com/) account. Pass your API token on the first run and it will be saved locally:

```bash
# Solve the 4-bit challenge curve:
python projecteleven.py --challenge 4 --token YOUR_IBM_TOKEN --backend ibm_marrakesh

# Subsequent runs (token already saved):
python projecteleven.py --challenge 4 --backend ibm_marrakesh

# Use the coordinate-based quantum oracle:
python projecteleven.py --challenge 4 --oracle coordinate --backend ibm_marrakesh

# Use the arithmetic oracle (coordinate encoding + QFT primitives):
python projecteleven.py --challenge 4 --oracle arithmetic --backend ibm_marrakesh

# Use a specific IBM Quantum instance:
python projecteleven.py --challenge 4 --instance ibm-q/open/main --backend ibm_marrakesh

# Verify curve parameters without quantum execution:
python projecteleven.py --curve curve_4 --verify-only
```

### CLI Options

| Flag | Description | Default |
|------|-------------|---------|
| `--challenge N` | Solve the N-bit challenge curve from `input_curves.json` | — |
| `--curve NAME` | Use a built-in test curve (`curve_4`) | — |
| `--token TOKEN` | IBM Quantum API token (saved locally on first use) | — |
| `--backend NAME` | IBM Quantum backend | `ibm_marrakesh` |
| `--instance ID` | IBM Quantum instance | `open-instance` |
| `--shots N` | Number of measurement shots | `8192` |
| `--oracle TYPE` | Oracle strategy: `dense`, `permutation`, `coordinate`, or `arithmetic` | auto |
| `--optimization-level N` | Qiskit transpilation optimization level (0-3) | `3` |
| `--d N` | Known secret key for testing (with `--curve`) | — |
| `--verify-only` | Validate curve parameters and exit | — |

## Project Structure

```
projecteleven.py        # Shor solver — dense unitary approach + CLI entry point
quantum_arithmetic.py   # Efficient permutation decomposition + QFT arithmetic primitives
quantum_oracle.py       # Coordinate-based oracle + arithmetic oracle framework
input_curves.json       # Challenge curves (4-bit to 30-bit)
problem/curves.py       # Curve generation utility
requirements.txt        # qiskit, qiskit-ibm-runtime
```

## References

- P. Shor, ["Algorithms for Quantum Computation: Discrete Logarithms and Factoring"](https://arxiv.org/abs/quant-ph/9508027) (1994)
- S. Beauregard, ["Circuit for Shor's algorithm using 2n+3 qubits"](https://arxiv.org/abs/quant-ph/0205095) (2003)
- M. Roetteler, M. Naehrig, K. Svore, K. Lauter, ["Quantum resource estimates for computing elliptic curve discrete logarithms"](https://arxiv.org/abs/1706.06752) (2017)

## License

This project is a submission to the Q-Day Prize Challenge release under [MIT LICENSE](LICENSE)
