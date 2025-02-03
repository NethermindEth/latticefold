# LatticeFold
A proof-of-concept implementation of the LatticeFold folding scheme, engineered by [Nethermind](https://nethermind.io), based on the work [LatticeFold: A Lattice-based Folding Scheme and its Applications to Succinct Proof Systems](https://eprint.iacr.org/2024/257) by Dan Boneh and Binyi Chen.

## Project Overview

LatticeFold implements a lattice-based folding scheme for succinct proof systems. It combines advanced cryptographic techniques to reduce the complexity of proofs while maintaining security and efficiency. This implementation demonstrates the potential of lattice-based cryptography for succinct proof systems.
### Applications
- **Recursive SNARKs**: Enables efficient recursive succinct non-interactive arguments of knowledge (SNARKs), essential for scaling blockchains and verifying large computations.  
- **Proof-Carrying Data (PCD)**: Facilitates the construction of PCD schemes, verifying sequences of computations to ensure correctness at each step.  
- **Post-Quantum Security**: Built on the Module Short Integer Solution (MSIS) problem, providing security against quantum attacks, making it ideal for post-quantum applications.  
- **Zero-Knowledge Proofs**: Constructs proofs that allow demonstrating knowledge of a statement without revealing the statement itself, crucial for privacy-preserving applications.  
- **Verifiable Computations**: Verifies computational correctness without re-executing the computations, ensuring computational integrity.  
- **Customizable Constraint Systems (CCS)**: Supports both low-degree (e.g., R1CS) and high-degree constraints, offering flexibility in computational proofs.  
- **Fully Homomorphic Encryption (FHE) and Lattice Signatures**: Operates on the same module structure as FHE and lattice-based signatures, benefiting from optimizations and hardware advancements in these areas.  


The project is supported by the Ethereum Foundation ZK Grant and is intended for research and experimentation, not production use.

### Key Features
- **Ajtai Commitment Scheme**: Provides cryptographic commitments.
- **R1CS/CCS Structures**: Facilitates succinct circuit representations.
- **Fiat-Shamir Transcript**: Adds non-interactive proof capabilities.

---

**DISCLAIMER:** This is a proof-of-concept prototype, and in particular has not received careful code review. This implementation is provided "as is" and NOT ready for production use. Use at your own risk.

## Building
### Prerequisites
- [Install Rust](https://www.rust-lang.org/tools/install) to manage toolchains.
- Use `rustup` to install the nightly toolchain pinned to `nightly-2024-11-05`:
  
   ```bash
   rustup install nightly-2024-11-05
  ```
- The [rust-toolchain](https://github.com/NethermindEth/latticefold/blob/main/rust-toolchain) file pins the version of the Rust toolchain, which the LatticeFold library builds with, to the specific version `nightly-2024-11-05`.
### Instructions
- Clone the repository:
```bash
git clone https://github.com/NethermindEth/latticefold.git
cd latticefold
```
- After that, use `cargo`, the standard Rust build tool, to build the library:
```bash
cargo build --release
```
## Troubleshooting
- Verify the correct toolchain is set using `rustup show`.
- Run `cargo update` to fetch the latest compatible dependencies.
- For platform-specific issues, consult the Rust documentation.
## Usage
Import the library to your `Cargo.toml`:
```toml
[dependencies]
latticefold = { git = "https://github.com/NethermindEth/latticefold.git", package = "latticefold" }
```

Available packages:
- `latticefold`: main crate, contains the non-interactive folding scheme implementation, together with the Ajtai commitment scheme, R1CS/CCS structures, Fiat-Shamir transcript machinery, etc.
- `cyclotomic-rings`: contains the trait definition of a ring suitable to be used in the LatticeFold protocol, a few ready-to-use rings and short challenge set machinery.

## Examples

- Check [latticefold/examples/README.md](latticefold/examples/README.md) for examples.
- Run an example using Cargo:
```bash
cargo run --example example_name
```


## Frontends

Currently, the only way to define a circuit to be folded is by specifying it as a [rank-1 constraint system (R1CS)](https://github.com/NethermindEth/latticefold/blob/main/latticefold/src/arith/r1cs.rs) - a representation of arithmetic circuits using degree-2 constraints or a [customizable constraint system (CCS)](https://github.com/NethermindEth/latticefold/blob/main/latticefold/src/arith.rs) - flexible circuit definitions that can handle constraints of any degree.

## License
The crates in this repository are licensed under either of the following licenses, at your discretion.

* Apache License Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
* MIT license ([LICENSE-MIT](LICENSE-MIT))

Unless you explicitly state otherwise, any contribution submitted for inclusion in this library by you shall be dual licensed as above (as defined in the Apache v2 License), without any additional terms or conditions.

## Acknowledgments

- This project is built on top of [our fork](https://github.com/NethermindEth/stark-rings) of [lattirust library](https://github.com/cknabs/lattirust) originally developed by [Christian Knabenhans](https://github.com/cknabs) and [Giacomo Fenzi](https://github.com/WizardOfMenlo). 
- We adapted [the sumcheck protocol from Jolt](https://github.com/a16z/jolt/blob/fa45507aaddb1815bafd54332e4b14173a7f8699/jolt-core/src/subprotocols/sumcheck.rs#L35) to the ring setting. 
- A lot of definitions are directly transferred from [sonobe](https://github.com/privacy-scaling-explorations/sonobe) library. 
- The implementation is supported by Ethereum Foundation [ZK Grant](https://blog.ethereum.org/2024/06/25/zk-grants-round-announce).

