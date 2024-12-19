from sage.all import *
from collections import defaultdict
from textwrap import dedent

# Define bound_2 function
def bound_2(d, kappa, p):
    return (2**(2 * sqrt(log(1.0045, 2) * d * kappa * log(p, 2)))).n()
# Define bound_inf function
def bound_inf(d, kappa, p, n):
    L = 1
    bound_value = floor(bound_2(d, kappa, p) / sqrt(d * (n * L)).n())
    # Ensure bound_value is a power of two
    if bound_value & (bound_value - 1) != 0:  # Check if not a power of two
        bound_value = 2**floor(log(bound_value, 2))  # Use log with base 2
    # Iterate until bound_value^L > p/2 or L exceeds 50
    while bound_value**L <= p:
        if L > 8:
            return "unpractical", "unpractical"
        L += 1
        bound_value = floor(bound_2(d, kappa, p) / sqrt(d * (n * L)).n())
        if bound_value & (bound_value - 1) != 0:  # Check if not a power of two
            # Find the largest power of two less than or equal to bound_value
            power_of_two = 2**floor(log(bound_value, 2))  # Use log with base 2
            bound_value = power_of_two  # Reduce to previous power of two
            if find_smallest_L_log(bound_value, p) != L:
                continue
    return bound_value, L
# Function to find the smallest L such that B^L > p using logarithms
def find_smallest_L_log(B, p):
    if B <= 0:
        return "unpractical"
    return ceil(log(p) / log(B))
# Function to find all (b, k) pairs such that b^k = B
def find_b_k_pairs(B):
    # Check if B is "unpractical"
    if B == "unpractical":
        return [("unpractical", "unpractical", "unpractical")]
    # Check if B is a power of two
    if B <= 0 or (B & (B - 1)) != 0:
        print("B is not a power of two")
        return [("unpractical", "unpractical", "unpractical")]
    k = int(log(B, 2))  # Calculate k such that 2^k = B using log with base 2
    return [2, k, B]
# Primes with their corresponding d values
params = {
    "BabyBear": {"p": 15 * 2^27 + 1, "d": 72},
    "Goldilocks": {"p": 2^64 - 2^32 + 1, "d": 24},
    "StarkPrime": {"p": 2^251 + (17 * 2^192) + 1, "d": 16},
    "Frog": {"p": 159120925213255836417, "d": 16},
}

# Range of num_cols values
num_cols_values = [2^9, 2^10, 2^11, 2^12, 2^13, 2^14, 2^15]

# Function to generate Rust macros and write to a file
def generate_macros_file():
    with open("latticefold/benches/config.toml", "w") as f:
        f.write("[benchmarks]\n")

        # Add more macros based on calculated parameters
        for prime_name, param in params.items():
            p = param["p"]
            d = param["d"]
            # f.write(f"\n//--- {prime_name} cyclotomic ring (modulus p = {p}, degree = {d}) ---\n")
                # Find the maximum kappa for which bound_2 < p / 2
            kappa = 1
            while bound_2(d, kappa, p) < p / 2:
                kappa += 1
            max_kappa = kappa - 1  # The last kappa where bound_2 was less than p / 2
            # f.write(f"//\tMaximum kappa for which bound_{{l_2}} < p/2: {max_kappa}")
            prime_lowercase = prime_name.lower()
            scalar_bench = f"{prime_lowercase}"
            non_scalar_bench = f"{prime_lowercase}_non_scalar"
            degree_three_non_scalar_bench = f"{prime_lowercase}_degree_three_non_scalar"
            bench_types = [scalar_bench, non_scalar_bench, degree_three_non_scalar_bench]
            for bench_type in bench_types:
                f.write(f"{bench_type} = [\n")
                kappa_values = range(1, max_kappa + 1)
                # Initialize a list to store all entries across all kappa values
                all_entries = []
                for kappa in kappa_values:
                    for n in num_cols_values:
                        # Calculate bound_inf for the current kappa and n
                        current_bound_inf, L = bound_inf(d, kappa, p, n)
                        # If the current bound is "unpractical", skip to the next kappa
                        if current_bound_inf == "unpractical":
                            continue
                        # Find all previous powers of two such that B^L > p/2
                        previous_powers_of_two = []
                        B = current_bound_inf
                        while B > 1:
                            if B**L > p:
                                previous_powers_of_two.append(B)
                            B //= 2  # Move to the previous power of two
                        # Display the results for each valid power of two
                        if previous_powers_of_two:
                            B_pow2 = min(previous_powers_of_two)
                            (b, k, B_pow2_in_pair) = find_b_k_pairs(B_pow2)
                            L = find_smallest_L_log(b**k, p)
                            all_entries.append((kappa, n, B_pow2_in_pair, L, b, k))
                # Sort all entries across all kappa values first by n, then by kappa
                all_entries.sort(key=lambda x: (x[1], x[0]))
                # Group entries by n
                entries_by_n = defaultdict(list)
                for entry in all_entries:
                    entries_by_n[entry[1]].append(entry)
                # Print results for each n
                for n, entries in entries_by_n.items():
                    min_by_kappa_entries = {entry[2:]: entry for entry in entries}
                    for entry in min_by_kappa_entries.values():
                        # f.write(f"\n\t\trun_single_{bench_type}_benchmark!(&mut $group, 1, {entry[0]}, {n}, {entry[2]}, {entry[3]}, {entry[4]}, {entry[5]});")
                        f.write(f"    {{ x_len = 1, c = {entry[0]}, w = {n}, b = \"{entry[2]}\", l = {entry[3]}, b_small = {entry[4]}, k = {entry[5]} }},\n")

                f.write("]\n\n")
                # f.write("   };\n")
                # f.write("}\n")

        f.write(dedent("""
            [ajtai]
            babybear = [
                { c = 1, w = 32768},
                { c = 2, w = 32768},
                { c = 3, w = 32768},
                { c = 4, w = 32768},
                { c = 5, w = 32768},
                { c = 6, w = 32768},
                { c = 1, w = 65536},
                { c = 2, w = 65536},
                { c = 3, w = 65536},
                { c = 4, w = 65536},
                { c = 5, w = 65536},
                { c = 6, w = 65536},
                { c = 1, w = 131072},
                { c = 2, w = 131072},
                { c = 3, w = 131072},
                { c = 4, w = 131072},
                { c = 5, w = 131072},
                { c = 6, w = 131072},
                { c = 1, w = 262144},
                { c = 2, w = 262144},
                { c = 3, w = 262144},
                { c = 4, w = 262144},
                { c = 5, w = 262144},
                { c = 6, w = 262144},
                { c = 1, w = 524288},
                { c = 2, w = 524288},
                { c = 3, w = 524288},
                { c = 4, w = 524288},
                { c = 5, w = 524288},
                { c = 6, w = 524288},
                { c = 1, w = 1048576},
                { c = 2, w = 1048576},
                { c = 3, w = 1048576},
                { c = 4, w = 1048576},
                { c = 5, w = 1048576},
                { c = 6, w = 1048576},
                { c = 12, w = 32768},
                { c = 13, w = 65536},
                { c = 13, w = 131072},
                { c = 14, w = 262144},
                { c = 14, w = 524288},
                { c = 15, w = 1048576},
            ]

            goldilocks = [
                { c = 1, w = 32768},
                { c = 2, w = 32768},
                { c = 3, w = 32768},
                { c = 4, w = 32768},
                { c = 5, w = 32768},
                { c = 6, w = 32768},
                { c = 1, w = 65536},
                { c = 2, w = 65536},
                { c = 3, w = 65536},
                { c = 4, w = 65536},
                { c = 5, w = 65536},
                { c = 6, w = 65536},
                { c = 1, w = 131072},
                { c = 2, w = 131072},
                { c = 3, w = 131072},
                { c = 4, w = 131072},
                { c = 5, w = 131072},
                { c = 6, w = 131072},
                { c = 1, w = 262144},
                { c = 2, w = 262144},
                { c = 3, w = 262144},
                { c = 4, w = 262144},
                { c = 5, w = 262144},
                { c = 6, w = 262144},
                { c = 1, w = 524288},
                { c = 2, w = 524288},
                { c = 3, w = 524288},
                { c = 4, w = 524288},
                { c = 5, w = 524288},
                { c = 6, w = 524288},
                { c = 1, w = 1048576},
                { c = 2, w = 1048576},
                { c = 3, w = 1048576},
                { c = 4, w = 1048576},
                { c = 5, w = 1048576},
                { c = 6, w = 1048576},
                { c = 17, w = 32768},
                { c = 17, w = 65536},
                { c = 18, w = 131072},
                { c = 19, w = 262144},
                { c = 19, w = 524288},
                { c = 20, w = 1048576},
            ]

            starkprime = [
                { c = 1, w = 32768},
                { c = 2, w = 32768},
                { c = 3, w = 32768},
                { c = 4, w = 32768},
                { c = 5, w = 32768},
                { c = 6, w = 32768},
                { c = 1, w = 65536},
                { c = 2, w = 65536},
                { c = 3, w = 65536},
                { c = 4, w = 65536},
                { c = 5, w = 65536},
                { c = 6, w = 65536},
                { c = 1, w = 131072},
                { c = 2, w = 131072},
                { c = 3, w = 131072},
                { c = 4, w = 131072},
                { c = 5, w = 131072},
                { c = 6, w = 131072},
                { c = 1, w = 262144},
                { c = 2, w = 262144},
                { c = 3, w = 262144},
                { c = 4, w = 262144},
                { c = 5, w = 262144},
                { c = 6, w = 262144},
                { c = 1, w = 524288},
                { c = 2, w = 524288},
                { c = 3, w = 524288},
                { c = 4, w = 524288},
                { c = 5, w = 524288},
                { c = 6, w = 524288},
                { c = 1, w = 1048576},
                { c = 2, w = 1048576},
                { c = 3, w = 1048576},
                { c = 4, w = 1048576},
                { c = 5, w = 1048576},
                { c = 6, w = 1048576},
                { c = 7, w = 131072},
                { c = 7, w = 262144},
                { c = 7, w = 524288},
                { c = 7, w = 1048576},
            ]

            frog = [
                { c = 1, w = 32768},
                { c = 2, w = 32768},
                { c = 3, w = 32768},
                { c = 4, w = 32768},
                { c = 5, w = 32768},
                { c = 6, w = 32768},
                { c = 1, w = 65536},
                { c = 2, w = 65536},
                { c = 3, w = 65536},
                { c = 4, w = 65536},
                { c = 5, w = 65536},
                { c = 6, w = 65536},
                { c = 1, w = 131072},
                { c = 2, w = 131072},
                { c = 3, w = 131072},
                { c = 4, w = 131072},
                { c = 5, w = 131072},
                { c = 6, w = 131072},
                { c = 1, w = 262144},
                { c = 2, w = 262144},
                { c = 3, w = 262144},
                { c = 4, w = 262144},
                { c = 5, w = 262144},
                { c = 6, w = 262144},
                { c = 1, w = 524288},
                { c = 2, w = 524288},
                { c = 3, w = 524288},
                { c = 4, w = 524288},
                { c = 5, w = 524288},
                { c = 6, w = 524288},
                { c = 1, w = 1048576},
                { c = 2, w = 1048576},
                { c = 3, w = 1048576},
                { c = 4, w = 1048576},
                { c = 5, w = 1048576},
                { c = 6, w = 1048576},
                { c = 17, w = 32768},
                { c = 18, w = 65536},
                { c = 19, w = 131072},
                { c = 19, w = 262144},
                { c = 20, w = 524288},
                { c = 21, w = 1048576},
            ]
            """))

# Call the function to generate the macros file
generate_macros_file()