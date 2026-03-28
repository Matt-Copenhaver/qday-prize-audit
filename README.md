# Shor's Algorithm for ECDLP — Q-Day Prize Submission

Quantum solver for the Elliptic Curve Discrete Logarithm Problem (ECDLP), built for the [Q-Day Prize Challenge](https://www.qdayprize.org/) by [Project Eleven](https://www.projecteleven.com/). The goal: recover ECC private keys on real quantum hardware using Shor's algorithm.

- **Author**: Giancarlo Lelli
- **Contact**: [gcarlo.lelli@gmail.com](mailto:gcarlo.lelli@gmail.com)
- **LinkedIn**: https://www.linkedin.com/in/giancarlolelli

## Approach

All challenge curves use **y² = x³ + 7** over F_p (a = 0, b = 7), matching the secp256k1 family.

The solver implements Shor's algorithm with two quantum register strategies:

- **Index-based encoding** (`projecteleven.py`) — Enumerates the elliptic curve group and encodes points as quantum basis state indices. Controlled point additions are applied as permutation unitaries. Practical for small curves (up to ~6-bit).

- **Efficient permutation decomposition** (`quantum_arithmetic.py`) — Replaces dense unitary matrices with cycle decomposition into controlled transpositions, reducing memory from O(2^{2n}) to O(N) and gate count to O(N·n) per controlled addition. Builds 13-bit circuits (41 qubits, ~4.1M gates) in ~15 seconds.

Both strategies follow the same Shor circuit structure:
1. Prepare counting registers |j⟩, |k⟩ in uniform superposition
2. Compute |j⟩|k⟩|jG + kQ⟩ via controlled point additions
3. Measure the point register
4. Apply inverse QFT to the counting registers
5. Extract d from the relationship j + kd ≡ r (mod n)

The module also includes QFT-based modular arithmetic primitives (Beauregard/Draper adders, modular multiplication, Fermat inversion) as building blocks toward a fully arithmetic coordinate-encoding approach that would scale polynomially to 256-bit keys.

All execution targets **IBM Quantum hardware** via `qiskit-ibm-runtime`. There is no local simulator path.

## References

- P. Shor, "Algorithms for Quantum Computation: Discrete Logarithms and Factoring" (1994)
- S. Beauregard, "Circuit for Shor's algorithm using 2n+3 qubits" (2003)
- M. Roetteler, M. Naehrig, K. Svore, K. Lauter, "Quantum resource estimates for computing elliptic curve discrete logarithms" (2017)

## Quick Start

```bash
git clone https://github.com/user/quantum.git
cd quantum

python -m venv .
source Scripts/activate        # Windows: Scripts\Activate.ps1

pip install -r requirements.txt
```

### Run

You need an [IBM Quantum](https://quantum.ibm.com/) account. Pass your API token on the first run and it will be saved locally:

```bash
# Solve the 4-bit challenge curve:
python projecteleven.py --challenge 4 --token YOUR_IBM_TOKEN --backend ibm_marrakesh

# Subsequent runs (token already saved):
python projecteleven.py --challenge 4 --backend ibm_marrakesh

# Use a specific IBM Quantum instance:
python projecteleven.py --challenge 4 --instance ibm-q/open/main --backend ibm_marrakesh

# Scalable solver (efficient decomposition):
python quantum_arithmetic.py --backend ibm_marrakesh --token YOUR_IBM_TOKEN

# Verify curve parameters without quantum execution:
python projecteleven.py --curve curve_4 --verify-only
```

### CLI Options

| Flag | Description | Default |
|------|-------------|---------|
| `--challenge N` | Solve the N-bit challenge curve from `input_curves.json` | — |
| `--curve NAME` | Use a built-in test curve (`curve_4`, `curve_13`) | — |
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
