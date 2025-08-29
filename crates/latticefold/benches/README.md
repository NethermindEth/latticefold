# Benchmarking with Environment Variables

This benchmark suite allows users to selectively enable or disable certain benchmarks using environment variables. By leveraging environment variables, users can tailor the benchmark runs to specific configurations.

---

## Environment Variables

The following environment variables can be set to **enable or filter benchmarks**:

### Boolean Flags (Enable/Disable Filtering)
- `GOLDILOCKS`
- `STARK`
- `BABYBEAR`
- `FROG`

- `PROVER`
- `VERIFIER`
- `AJTAI`

- `LINEARIZATION`
- `DECOMPOSITION`
- `FOLDING`
- `E2E`

> **Default Behavior**:  
If none of the flags for a group are set, all flags are enabled by default.

### Optional Variables (Numeric Values)
The following variables control specific parameters for benchmarks. If not set, they will be ignored:
- `DURATION` (default: `50.0` seconds as floating point number) – Duration for benchmarks
- `WARMUP` (default: `1.0` seconds as floating point number) – Warmup time

---

## How It Works
- If none of `LINEARIZATION`, `DECOMPOSITION`, `FOLDING`, and `E2E` is set, all benchmarks would be run.
- If some of them is set, only set benchmarks would be run.
- Similarly, if none of `PROVER`, `VERIFIER` and `AJTAI` is set, all benchmarks would be run, otherwise,
  only set benchmarks would be run.
- `GOLDILOCKS`, `STARK`, `BABYBEAR` and `FROG` are used similarly.

---

## Examples

### **Bash**

1. **Run only `LINEARIZATION` benchmarks:**
   ```bash
   LINEARIZATION=1 cargo bench
   ```

2. **Run only `PROVER` `DECOMPOSITION` benchmarks:**
   ```bash
   PROVER=1 DECOMPOSITION=1 cargo bench
   ```

3. **Run benchmarks matching specific parameters for GOLDILOCKS ring:**
   ```bash
   KAPPA=10 GOLDILOCKS=1 L=8 GOLDILOCKS=1 cargo bench
   ```

---

### **PowerShell**
1. **Run only `LINEARIZATION` benchmarks:**
   ```powershell
   $env:LINEARIZATION=1; cargo bench
   ```

2. **Run only `PROVER` `DECOMPOSITION` benchmarks:**
   ```powershell
   $env:PROVER=1; $env:DECOMPOSITION=1; cargo bench
   ```

3. **Run benchmarks matching specific parameters for GOLDILOCKS ring:**
   ```powershell
   $env:KAPPA=10; $env:GOLDILOCKS=1; $env:L=8; $env:GOLDILOCKS=1; cargo bench
   ```

---

## Numeric Parameters

Numeric parameters can be used to filter benchmarks, if they are set they only benchmarks with matching parameter would be run

- X_LEN 
- KAPPA 
- W 
- WIT_LEN 
- B 
- L 
- B_SMALL 
- K

