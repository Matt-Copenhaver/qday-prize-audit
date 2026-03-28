# Shor's Algorithm for ECDLP — Q-Day Prize Submission

Quantum solver for the Elliptic Curve Discrete Logarithm Problem (ECDLP), built for the [Q-Day Prize Challenge](https://www.qdayprize.org/) by [Project Eleven](https://www.projecteleven.com/). The goal: recover ECC private keys on real quantum hardware using Shor's algorithm.

- **Author**: Giancarlo Lelli
- **Contact**: [gcarlo.lelli@gmail.com](mailto:gcarlo.lelli@gmail.com)
- **LinkedIn**: https://www.linkedin.com/in/giancarlolelli

## Approach description

All challenge curves use **y^2 = x^3 + 7** over F_p (a = 0, b = 7), matching the secp256k1 family. The solver selects the appropriate strategy automatically based on curve size.

Both strategies implement the same Shor circuit structure:
1. Prepare counting registers |j>, |k> in uniform superposition (Hadamard)
2. Compute |j>|k>|jG + kQ> via 2t controlled point additions (t = num_counting qubits)
3. Measure the point register, collapsing it to some group element R
4. Apply inverse QFT to the counting registers
5. Measure j, k and extract d from the relation j + kd = r (mod n)

The private key d is recovered by collecting multiple (j, k) samples that satisfy the same linear relation modulo the group order n.

### Strategy 1: Dense Unitary (n_bits <= 6)

Used for curves with group order up to ~6 bits. Implemented in `projecteleven.py`.

Each controlled point addition "add S" is represented as a **2^(n+1) x 2^(n+1) permutation matrix** applied via `qc.unitary()`. The matrix encodes the full group action: the upper-left block is identity (control=0), the lower-right block permutes basis states according to the map P -> P+S (control=1).

- **Memory**: O(2^{2n}) per matrix -- a 6-bit curve produces 128x128 matrices (~130 KB each), but a 13-bit curve would need 16384x16384 matrices (~2 GB each)
- **Qubits**: 2t + n (two counting registers + point register)
- **Limitation**: Qiskit's unitary decomposition is O(4^n), making this infeasible beyond ~6-bit

### Strategy 2: Efficient Permutation Decomposition (n_bits > 6)

Used for larger curves. Implemented in `quantum_arithmetic.py`.

Instead of building dense matrices, each "add S" permutation is **cycle-decomposed into transpositions**. Each transposition (swap of two basis states |a> <-> |b>) is implemented with:

1. **CNOT reduction** -- CNOTs from a pivot bit to all other differing bits between a and b, reducing the multi-bit difference to a single-bit difference
2. **Multi-controlled X** -- An MCX gate on the pivot bit, conditioned on all other bits matching the target pattern. Uses X gates to flip control bits that need to be |1> for the condition
3. **Undo CNOTs** -- Reverse step 1 to restore the non-pivot bits

The MCX uses **V-chain decomposition** with (n-2) dedicated ancilla qubits, giving O(n) Toffoli gates per MCX instead of O(n^2) without ancillas.

Each controlled addition is built as an **isolated sub-circuit** and appended to the main circuit as a single opaque gate. This avoids quadratic DAG growth in Qiskit -- the main circuit sees ~28 top-level gates regardless of how many transpositions each contains internally.

| Metric | Dense Unitary | Efficient Permutation |
|--------|--------------|----------------------|
| Memory per addition | O(2^{2n}) | O(N) |
| Gates per addition | O(4^n) decomposition | O(N * n) |
| Qubits | 2t + n | 2t + n + (n-2) ancillas |
| 4-bit (n=7) | 0.01s build | 0.01s build |
| 9-bit (n=313) | infeasible | 36 qubits, ~143K gates |
| 13-bit (n=4243) | infeasible | 52 qubits, ~4.1M gates |

Where N = group order, n = bit length, t = n+1 counting qubits.

### QFT Arithmetic Primitives

The codebase also includes QFT-based modular arithmetic building blocks (Beauregard/Draper adders, modular multiplication, Fermat inversion) as a foundation toward a fully arithmetic coordinate-encoding approach that would scale polynomially to 256-bit keys.

All execution targets **IBM Quantum hardware** via `qiskit-ibm-runtime`. There is no local simulator path.

## References

- P. Shor, ["Algorithms for Quantum Computation: Discrete Logarithms and Factoring"](https://arxiv.org/abs/quant-ph/9508027) (1994)
- S. Beauregard, ["Circuit for Shor's algorithm using 2n+3 qubits"](https://arxiv.org/abs/quant-ph/0205095) (2003)
- M. Roetteler, M. Naehrig, K. Svore, K. Lauter, ["Quantum resource estimates for computing elliptic curve discrete logarithms"](https://arxiv.org/abs/1706.06752) (2017)

## Results

Successfully recovered private keys on IBM Quantum hardware for challenge curves up to 9-bit:

| Challenge | p | n | Strategy | Qubits | 2Q Gates | Transpiled Depth | Backend | Recovered d | Job ID |
|-----------|-----|------|----------|--------|----------|------------------|-------------|-------------|--------|
| 4-bit | 13 | 7 | Dense unitary | 11 | 774 | 2,425 | ibm_torino | 6 | d73u28kvllmc73anvi90 |
| 6-bit | 43 | 31 | Dense unitary | 17 | 23,471 | 72,475 | ibm_torino | 18 | d73u2l5koquc73e24u8g |
| 8-bit | 163 | 139 | Efficient permutation | 32 | 294,628 | 599,517 | ibm_kingston | 103 | d73ui15koquc73e25e4g |
| 9-bit | 349 | 313 | Efficient permutation | 36 | 887,544 | 1,764,266 | ibm_torino | 135 | d73ua2h8qmgc73flei9g |

All runs used 8,192 shots and were executed on the IBM Quantum open-instance plan, which grants 10 minutes of free quantum computation per month. Full execution logs are in the `executions/` folder.

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
| `--d N` | Known secret key for testing (with `--curve`) | — |
| `--verify-only` | Validate curve parameters and exit | — |

## Project Structure

```
projecteleven.py        # Shor solver — dense unitary approach
quantum_arithmetic.py   # Shor solver — efficient permutation decomposition + QFT arithmetic
input_curves.json       # Challenge curves (4-bit to 30-bit)
problem/curves.py       # Curve generation utility
requirements.txt        # qiskit, qiskit-ibm-runtime
```

## License

This project is a submission to the Q-Day Prize Challenge release under [MIT LICENSE](LICENSE)
