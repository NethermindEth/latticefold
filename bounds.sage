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

# Primes with their corresponding d values
params = {
    "BabyBear": {"p": 15 * 2^27 + 1, "d": 72},
    "Goldilocks": {"p": 2^64 - 2^32 + 1, "d": 24},
    "Stark": {"p": 2^251 + (17 * 2^192) + 1, "d": 16},
    "Frog": {"p": 159120925213255836417 * 2^24 + 1, "d": 16},
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

    # Compute values of bound_inf, log_2(bound_inf), and B for kappa from 1 to max_kappa and num_cols values
    bound_inf_values = [
        (kappa, num_cols, bound_inf(d, kappa, p, num_cols), log(bound_inf(d, kappa, p, num_cols), 2), adjusted_B(bound_inf(d, kappa, p, num_cols)))
        for kappa in kappa_values for num_cols in num_cols_values
    ]
    
    # Display bound_inf results with B
    print("\nValues of bound_inf, log_2(bound_inf), and B:")
    for kappa, num_cols, b_inf, log_b_inf, B in bound_inf_values:
        print(f"kappa = {kappa}, num_cols = {num_cols}, bound_inf = {b_inf}, log_2(bound_inf) = {log_b_inf}, B = {B}")

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
            print(f"num_cols = {num_cols}: kappa = {closest_kappa}, B = {closest_B} = 2^{log(closest_B, 2).n()}")
        else:
            print(f"num_cols = {num_cols}: No suitable kappa found")
