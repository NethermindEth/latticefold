from sage.all import *
from collections import defaultdict
# Define bound_2 function
def bound_2(d, kappa, p):
    return (2**(2 * sqrt(log(1.0045, 2) * d * kappa * log(p, 2)))).n()
# Define bound_inf function
def bound_inf(d, kappa, p, n):
    L = 1
    bound_value = floor(bound_2(d, kappa, p) / sqrt(d * (n * L)).n())
    # Ensure bound_value is a power of two
    bound_value = 2**(bound_value.bit_length() - 1)
    # Iterate until bound_value^L > p/2 or L exceeds 50
    while bound_value^L <= p:
        if L > 8:
            return "unpractical", "unpractical"
        L += 1
        bound_value = floor(bound_2(d, kappa, p) / sqrt(d * (n * L)).n())
        if bound_value & (bound_value - 1) != 0:  # Check if not a power of two
            bound_value = 2**(bound_value.bit_length() - 1)  # Reduce to previous power of two
            if find_smallest_L_log(bound_value, p) != L:
                continue
    return bound_value, L
# Function to find the smallest L such that B^L > p/2 using logarithms
def find_smallest_L_log(B, p):
    if B <= 0:
        return "unpractical"
    return ceil(log(p) / log(B))
# Function to find all (b, k) pairs such that b^k = B
def find_b_k_pairs(B):
    # Check if B is "unpractical"
    if B == "unpractical":
        return [("unpractical", "unpractical", "unpractical")]
    if 2^(B.bit_length()-1) != B:
        print("B is not a power of two")
        return [("unpractical", "unpractical", "unpractical")]
    return [2, B.bit_length() - 1, B]
# Primes with their corresponding d values
params = {
#    "BabyBear": {"p": 15 * 2^27 + 1, "d": 72},
    "Goldilocks": {"p": 2^64 - 2^32 + 1, "d": 24},
#    "StarkPrime": {"p": 2^251 + (17 * 2^192) + 1, "d": 16},
#    "Frog": {"p": 159120925213255836417, "d": 16},
#    "Dilithium": {"p": 2^23 - 2^13 + 1, "d": 256}
}
# Range of num_cols values
num_cols_values = [2^9, 2^10, 2^11, 2^12, 2^13, 2^14]
# Iterate over each prime and calculate the maximum kappa and perform bound calculations
for prime_name, param in params.items():
    p = param["p"]
    d = param["d"]
    print("")
    print(f"--- {prime_name} prime modulus (d = {d}) ---")
    # Find the maximum kappa for which bound_2 < p / 2
    kappa = 1
    while bound_2(d, kappa, p) < p / 2:
        kappa += 1
    max_kappa = kappa - 1  # The last kappa where bound_2 was less than p / 2
    print(f"\tMaximum kappa for which bound_2 < p/2: {max_kappa}")
    # Define kappa_values from 1 to min(50, max_kappa)
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
                if B**L > p / 2:
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
        min_entries = {}
        for entry in entries:
            key = (entry[2], entry[3], entry[4], entry[5])  # B, L, b, k
            if key not in min_entries or entry[0] < min_entries[key][0]:
                min_entries[key] = entry
        for entry in min_entries.values():
            benchmark_function = f"run_single_{prime_name.lower()}_benchmark!"
            print(f"\t{benchmark_function}(&mut $group, 1, {entry[0]}, {n}, {entry[2]}, {entry[3]}, {entry[4]}, {entry[5]});")
    print("")
    for n, entries in entries_by_n.items():
        min_entries = {}
        for entry in entries:
            key = (entry[2], entry[3], entry[4], entry[5])  # B, L, b, k
            if key not in min_entries or entry[0] < min_entries[key][0]:
                min_entries[key] = entry
        for entry in min_entries.values():
            benchmark_function = f"run_single_{prime_name.lower()}_non_scalar_benchmark!"
            print(f"\t{benchmark_function}(&mut $group, 1, {entry[0]}, {n}, {entry[2]}, {entry[3]}, {entry[4]}, {entry[5]});")