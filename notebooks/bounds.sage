from sage.all import *

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
    while bound_value^L <= p / 2:
        if L > 50:
            return "unpractical", "unpractical"
        L += 1
        bound_value = floor(bound_2(d, kappa, p) / sqrt(d * (n * L)).n())
        if bound_value % 2 == 1:
            bound_value -= 1
            if find_smallest_L_log(bound_value, p) != L:
                continue
    
    return bound_value, L

# Function to find the smallest L such that B^L > p/2 using logarithms
def find_smallest_L_log(B, p):
    if B <= 0:
        return "unpractical"
    return ceil(log(p / 2) / log(B))

# Function to find all (b, k) pairs such that b^k = B
def find_b_k_pairs(B, original_L, target_pairs=10):
    pairs = []
    
    # Check if B is "unpractical"
    if B == "unpractical":
        return [("unpractical", "unpractical", "unpractical")]
    
    # Ensure B is a power of two
    B = 2**(B.bit_length() - 1)
    
    # Start from the target B and decrement until finding the required number of pairs
    while len(pairs) < target_pairs and B > 1:
        # Iterate over even values of b only
        b = 2
        while b <= B:
            if b > 64:
                break
            k = 1
            power = b
            while power <= B:
                L = find_smallest_L_log(b**k, p)
                if L == original_L:
                    pairs.append((b, k, B))  # Add B to the pair if L matches
                k += 1
                power *= b
            b *= 2 # Move to the next power of two
        
        # Sort pairs by how close b^k is to B, in descending order
        pairs.sort(key=lambda pair: abs(pair[0]**pair[1] - B))
        
        # If we have enough pairs, break the loop
        if len(pairs) >= target_pairs:
            break
        
        # Otherwise, decrement B by 2 to ensure it remains even
        B -= 2
    
    # Use a dictionary to track the maximum B for each k
    max_b_for_k = {}

    for b, k, B in pairs:
        if k not in max_b_for_k or B > max_b_for_k[k][2]:
            max_b_for_k[k] = (b, k, B)

    # Extract the pairs with the greatest B for each k
    unique_pairs = list(max_b_for_k.values())

    # Return the closest pairs found, limited to target_pairs
    return unique_pairs[:target_pairs]

# Primes with their corresponding d values
params = {
#    "BabyBear": {"p": 15 * 2^27 + 1, "d": 72},
    "Goldilocks": {"p": 2^64 - 2^32 + 1, "d": 24},
#    "Stark": {"p": 2^251 + (17 * 2^192) + 1, "d": 16},
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
    
    # Limit kappa to 16 for a specific case study
    if prime_name == "Goldilocks":  # Replace "Goldilocks" with your specific case
        max_kappa = min(max_kappa, 16)
    
    print(f"\tMaximum kappa for which bound_2 < p/2: {max_kappa}")
    
    # Define kappa_values from 1 to min(50, max_kappa)
    kappa_values = range(1, min(50, max_kappa) + 1)

    # Iterate over each kappa value
    for kappa in kappa_values:
        for n in num_cols_values:
            # Calculate bound_inf for the current kappa and n
            current_bound_inf, L = bound_inf(d, kappa, p, n)
            
            # If the current bound is "unpractical", skip to the next kappa
            if current_bound_inf == "unpractical":
                continue
            
            # Find the pairs (b, k) such that b^k = B
            pairs = find_b_k_pairs(current_bound_inf, L)
            
            # Display the results
            print(f"\tkappa = {kappa}, n = {n}: B = {current_bound_inf}")
            for b, k, B in pairs:
                # Recalculate L for each pair
                L = find_smallest_L_log(b**k, p)
                
                # Discard any L equal to 1 or less and where L*n > 2^14
                if L > 1 and L * n <= 2^14:
                    print(f"run_single_{prime_name.lower()}_benchmark!(&mut group, 1, {kappa}, {n}, {current_bound_inf}, 1, {b}, {k});")
