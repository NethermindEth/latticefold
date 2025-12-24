# LatticeFold+

A developing implementation of [LatticeFold+: Faster, Simpler, Shorter Lattice-Based Folding for Succinct Proof Systems](https://eprint.iacr.org/2025/247), a more performant<sup>1<sup> version of LatticeFold.

**DISCLAIMER:** This is a proof-of-concept prototype, and in particular has not received careful code review. This implementation is provided "as is" and NOT ready for production use. Use at your own risk.

<sup>1<sup> This is currently a work-in-progress and current performance is not final.

## Benches
Run benches with,

```sh
cargo bench
```

## TODO
- Add / use 128-bit ring / or support 128-bit security through extensions. |S| = 128 bit, |C| = 64 bit;
- Parameter tuning (B, kappa, etc);
- Some sumchecks over Zq (avoid Rq promotions);
- Various paper-provided performance improvements;
- Explore NTT (e.g. `SuitableRing`);
- Support CCS;
- Enhance main (Prover/Verifier) and more internal interfaces (sub-protocols): target real usability, reduce cloning;
- Add error types / error handling. Remove `unwrap`s;
- Docs;
- Add examples;
- Add more benches for each subprotocol. Add benches targeting LF comparison (with equivalent parameters).
