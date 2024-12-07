#!/bin/bash

# Run the single_operations benchmark
cargo bench --bench single_operations || echo "single_operations benchmark failed"

# Run the linearization_ops benchmark
cargo bench --bench linearization_ops || echo "linearization_ops benchmark failed" 

# Run the decomposition_ops benchmark
cargo bench --bench decomposition_ops || echo "decomposition_ops benchmark failed" 

# Run the folding_ops benchmark
cargo bench --bench folding_ops || echo "folding_ops benchmark failed"
