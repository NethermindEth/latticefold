# LatticeFold

A proof-of-concept implementation of the LatticeFold folding scheme engineered by [Nethermind](nethermind.io) based on the work 
[LatticeFold: A Lattice-based Folding Scheme and its Applications to Succinct Proof Systems](https://eprint.iacr.org/2024/257) by Dan Boneh and Binyi Chen.

**DISCLAIMER:** This is a proof-of-concept prototype, and in particular has not received careful code review. This implementation is provided "as is" and NOT ready for production use. Use at your own risk.

## Building

The [rust-toolchain](https://github.com/NethermindEth/latticefold/blob/main/rust-toolchain) file pins the version of the Rust toolchain, which the LatticeFold library builds with, to the specific version `nightly-2024-11-05`.

One can install the `nightly-2024-11-05` toolchain by invoking:
```bash
rustup install nightly-2024-11-05
```

After that, use `cargo`, the standard Rust build tool, to build the library:

```bash
git clone https://github.com/NethermindEth/latticefold.git
cd latticefold
cargo build --release
```

## Usage
Import the library:
```toml
[dependencies]
latticefold = { git = "https://github.com/NethermindEth/latticefold.git", package = "latticefold" }
```

Available packages:
- `latticefold`: main crate, contains the non-interactive folding scheme implementation, together with the Ajtai commitment scheme, R1CS/CCS structures, Fiat-Shamir transcript machinery, etc.
- `cyclotomic-rings`: contains the trait definition of a ring suitable to be used in the LatticeFold protocol, a few ready-to-use rings and short challenge set machinery.

## Examples

Check [latticefold-examples/README.md](examples/README.md) for examples.

## Frontends

Currently, the only way to define a circuit to be folded is by specifying it as a [rank-1 constraint system (R1CS)](https://github.com/NethermindEth/latticefold/blob/main/latticefold/src/arith/r1cs.rs) or a [customizable constraint system (CCS)](https://github.com/NethermindEth/latticefold/blob/main/latticefold/src/arith.rs).

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
