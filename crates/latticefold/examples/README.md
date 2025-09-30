# Examples README

This file explains how to use the examples in this repository. Examples demonstrate functionality and can be customized using environment variables. Instructions are provided for Linux/MacOS (bash/zsh) and Windows (PowerShell).

## Implemented examples

- e2e

## Customization with Environment Variables

The examples in this repository support customization via environment variables. Most examples use the following parameters, except for `starkprime`, which has its own set of parameters detailed below, to tailor their behavior:

- **`RING`**: Establishes the used ring. Ring configurations are: `Goldilocks`, `BabyBear`, `Frog`, and `StarkPrime`.
    - Default: `Goldilocks`
- **`PARAM_B`**: Sets the value of `B` in `DecompositionParams`.
    - Default: `32768` (`1 << 15`)
- **`PARAM_L`**: Sets the value of `L` in `DecompositionParams`.
    - Default: `5`
- **`PARAM_B_SMALL`**: Sets the value of `B_SMALL` in `DecompositionParams`. Must be even.
    - Default: `2`
- **`PARAM_K`**: Sets the value of `K` in `DecompositionParams`.
    - Default: `15`
- **`PARAM_KAPPA`**: Sets the value of `KAPPA`, controlling challenge set parameters.
    - Default: `4`
- **`PARAM_WIT_LEN`**: Sets the witness length.
    - Default: `4`

### starkprime default values

- **`PARAM_B`**: Sets the value of `B` in `DecompositionParams`.
    - Default: `1073741824u128`
- **`PARAM_L`**: Sets the value of `L` in `DecompositionParams`.
    - Default: `9`
- **`PARAM_B_SMALL`**: Sets the value of `B_SMALL` in `DecompositionParams`.
    - Default: `2`
- **`PARAM_K`**: Sets the value of `K` in `DecompositionParams`.
    - Default: `30`
- **`PARAM_KAPPA`**: Sets the value of `KAPPA`, controlling challenge set parameters.
    - Default: `4`
- **`PARAM_WIT_LEN`**: Sets the witness length.
    - Default: `4`

These parameters influence the behavior and output of the examples.

## Setting Environment Variables

### Linux/MacOS (bash/zsh)

1. Open a terminal.

2. Export the desired environment variables before running the examples. For example:

   ```bash
   export PARAM_B=65536
   export PARAM_L=6
   export PARAM_B_SMALL=2
   export PARAM_K=16
   export PARAM_KAPPA=5
   export PARAM_WIT_LEN=5

   cargo run --example <example_name>
   ```

3. Replace `<example_name>` with the name of the example you want to run.

### Windows (PowerShell)

1. Open PowerShell.

2. Set the desired environment variables before running the examples. For example:

   ```powershell
   $env:PARAM_B=65536
   $env:PARAM_L=6
   $env:PARAM_B_SMALL=2
   $env:PARAM_K=16
   $env:PARAM_KAPPA=5
   $env:PARAM_WIT_LEN=5

   cargo run --example <example_name>
   ```

3. Replace `<example_name>` with the name of the example you want to run.

## Example Output

When you modify environment variables, the generated parameters are automatically updated in the example's output. This allows for testing different configurations and validating results under various conditions.

## Default Values

If no environment variables are specified, the examples will run with the following defaults:

- `PARAM_B`: `32768`
- `PARAM_L`: `5`
- `PARAM_B_SMALL`: `2`
- `PARAM_K`: `15`
- `PARAM_KAPPA`: `4`
- `PARAM_WIT_LEN`: `4`

## Notes

- Ensure you rebuild the examples after modifying environment variables to see the changes.

  ```bash
  cargo clean && cargo run --example <example_name>
  ```

- For detailed instructions on each example, refer to the example's source code or inline comments.

