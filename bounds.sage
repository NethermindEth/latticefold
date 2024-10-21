from sage.all import *

# Define bound_2 function
def bound_2(d, kappa, p):
    return (2**(2 * sqrt(log(1.0045, 2) * d * kappa * log(p, 2)))).n()

# Define bound_inf function
def bound_inf(d, kappa, p, n):
    return (bound_2(d, kappa, p) / sqrt(d * n)).n()

# Define B function (bound_inf ceiling with adjustment for odd values)
def adjusted_B(bound_inf_value):
    B = ceil(bound_inf_value)
    if B % 2 == 1:  # Check if B is odd
        B -= 1  # Decrease by 1 if odd
    return B

# Function to find the smallest L such that B^L > p/2 using logarithms
def find_smallest_L_log(B, p):
    if B <= 0:
        return "unpractical"
    return ceil(log(p / 2) / log(B))

# Primes with their corresponding d values
params = {
    "BabyBear": {"p": 15 * 2^27 + 1, "d": 72},
    "Goldilocks": {"p": 2^64 - 2^32 + 1, "d": 24},
    "Stark": {"p": 2^251 + (17 * 2^192) + 1, "d": 16},
    "Frog": {"p": 159120925213255836417, "d": 16},
    "Dilithium": {"p": 2^23 - 2^13 + 1, "d": 256}
}

# Range of num_cols values
num_cols_values = [2^15, 2^16, 2^17, 2^18, 2^19, 2^20]

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
    print(f"Maximum kappa for which bound_2 < p/2: {max_kappa}")
    
    # Define kappa_values up to max_kappa
    kappa_values = range(1, max_kappa + 1)
    
    # Compute values of bound_2 and log_2(bound_2) for kappa from 1 to max_kappa
    bound_2_values = [(kappa, bound_2(d, kappa, p), log(bound_2(d, kappa, p), 2)) for kappa in kappa_values]
    
    # Display bound_2 results
    print("Values of bound_2 and log_2(bound_2):")
    for kappa, b2, log_b2 in bound_2_values:
        print(f"kappa = {kappa}, bound_2 = {b2}, log_2(bound_2) = {log_b2}")

    # Compute values of bound_inf, log_2(bound_inf), B, and L for kappa from 1 to max_kappa and num_cols values
    bound_inf_values = [
        (kappa, num_cols, bound_inf(d, kappa, p, num_cols), log(bound_inf(d, kappa, p, num_cols), 2), adjusted_B(bound_inf(d, kappa, p, num_cols)))
        for kappa in kappa_values for num_cols in num_cols_values
    ]
    
    # Display bound_inf results with B and L
    print("\nValues of bound_inf, log_2(bound_inf), B, and L:")
    for kappa, num_cols, b_inf, log_b_inf, B in bound_inf_values:
        if B == 0:
            print(f"kappa = {kappa}, num_cols = {num_cols}, B = {B}, L = unpractical")
        else:
            L = find_smallest_L_log(B, p)
            print(f"kappa = {kappa}, num_cols = {num_cols}, bound_inf = {b_inf}, log_2(bound_inf) = {log_b_inf}, B = {B}, L = {L}")

    # Find kappa for B = 2^16
    target_B = 2^16
    margin = 2^16
    print(f"\nKappa values for B ≈ 2^16 (within a margin of {margin}):")
    for num_cols in num_cols_values:
        closest_kappa = None
        closest_B = None
        smallest_difference = float('inf')
        
        for kappa in range(1, max_kappa + 1):
            B = adjusted_B(bound_inf(d, kappa, p, num_cols))
            difference = abs(B - target_B)
            
            if difference <= margin and difference < smallest_difference:
                closest_kappa = kappa
                closest_B = B
                smallest_difference = difference
        
        if closest_kappa is not None:
            L = find_smallest_L_log(closest_B, p)
            print(f"num_cols = {num_cols}: kappa = {closest_kappa}, B = {closest_B} = 2^{log(closest_B, 2).n()}, L = {L}")
        else:
            print(f"num_cols = {num_cols}: No suitable kappa found")
        
    # Calculate the largest B for max_kappa where L is an integer for all num_cols
    print("\nCalculating largest B for max_kappa where L is an integer for all num_cols:")
    for num_cols in num_cols_values:
        largest_B_with_integer_L = None
        B = adjusted_B(bound_inf(d, max_kappa, p, num_cols))
        L = find_smallest_L_log(B, p)
        
        largest_B_with_integer_L = B

        if largest_B_with_integer_L is not None:
            print(f"num_cols = {num_cols}: Largest B with integer L,  B = {largest_B_with_integer_L}, L = {L}")
        else:
            print(f"num_cols = {num_cols}: No B found with integer L")

    # Calculate the range of B for which L = 2, 3, 4, 5
    for L_target in range(2, 6):
        print(f"\nCalculating range of B for which L = {L_target} for all num_cols:")
        for num_cols in num_cols_values:
            min_B = None
            max_B = None
            min_kappa = None
            # Check max_kappa
            p_div_L = p / (2**L_target)
            
            # Start with B = 1 and increment to find the minimum B
            for kappa in range(1, max_kappa + 1):
                B = adjusted_B(bound_inf(d, kappa, p, num_cols))
                if B**L_target > p_div_L:
                    min_B = B
                    min_kappa = kappa
                    break
            
            # Use the largest B logic to find the maximum B for which L = L_target
            for kappa in range(max_kappa, 0, -1):
                B = adjusted_B(bound_inf(d, kappa, p, num_cols))
                if B**L_target > p_div_L:
                    max_B = B
                    max_kappa = kappa
                    break
            
            if min_B is not None and max_B is not None:
                print(f"num_cols = {num_cols}: Range of B for L = {L_target} is [{min_B} (kappa = {min_kappa}), {max_B} (kappa = {max_kappa})]")
            else:
                print(f"num_cols = {num_cols}: No range of B found for L = {L_target}")

    
