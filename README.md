# LatticeFold

A proof-of-concept implementation of the LatticeFold folding scheme engineered by [Nethermind](nethermind.io) based on the work 
[LatticeFold: A Lattice-based Folding Scheme and its Applications to Succinct Proof Systems](https://eprint.iacr.org/2024/257) by Dan Boneh and Binyi Chen.

**DISCLAIMER:** This is a proof-of-concept prototype, and in particular has not received careful code review. This implementation is provided "as is" and NOT ready for production use. Use at your own risk.

## Building

The [rust-toolchain](https://github.com/NethermindEth/latticefold/blob/main/rust-toolchain) file pins the version of the Rust toolchain, which the LatticeFold library builds with, to the specific version `nightly-2024-06-05`. This is mainly because we are substantially using an incomplete rustc feature `generic_const_exprs`.

One can install the `nightly-2024-06-05` toolchain by invoking:
```bash
rustup install nightly-2024-06-05
```

After that, use `cargo`, the standard Rust build tool, to build the libraries:

```bash
git clone https://github.com/NethermindEth/latticefold.git
cd latticefold
cargo build --release
```

## Usage
Import the library:
```toml
[dependencies]
latticefold = { git = "https://github.com/NethermindEth/latticefold.git", package = "latticefold"}
```

Available packages:
- `latticefold`: main crate, contains the non-interactive folding scheme implementation, together with the Ajtai commitment scheme, R1CS/CCS structures, Fiat-Shamir transcript machinery, etc.
- `cyclotomic-rings`: contains the trait definition of a ring suitable to be used in the LatticeFold protocol, a few ready-to-use rings and short challenge set machinery.

## Frontends

Currently, the only way to define a circuit to be folded is by specifying it as a [rank-1 constraint system (R1CS)](https://github.com/NethermindEth/latticefold/blob/main/latticefold/src/arith/r1cs.rs) or a [customizable constraint system (CCS)](https://github.com/NethermindEth/latticefold/blob/main/latticefold/src/arith.rs).

## License

LatticeFold is [MIT Licensed](https://github.com/NethermindEth/latticefold/blob/main/LICENSE).

## Acknowledgments

- This project is built on top of [our fork](https://github.com/NethermindEth/lattirust) of [lattirust library](https://github.com/cknabs/lattirust) originally developed by [Christian Knabenhans](https://github.com/cknabs) and [Giacomo Fenzi](https://github.com/WizardOfMenlo). 
- We almost verbatim adapt [the sumcheck protocol from Nexus zkVM](https://github.com/nexus-xyz/nexus-zkvm/blob/f37401c477b680ce5334b2ca523ded8a7273d8c8/nova/src/folding/hypernova/ml_sumcheck/mod.rs) to the ring setting. 
- A lot of definitions are directly transferred from [sonobe](https://github.com/privacy-scaling-explorations/sonobe) library. 
- The implementation is supported by Ethereum Foundation [ZK Grant](https://blog.ethereum.org/2024/06/25/zk-grants-round-announce).