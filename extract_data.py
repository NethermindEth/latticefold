import re
import csv

        #r'(Linearization|Decomposition|Folding) Godlilocks/(Linearization|Decomposition|Folding) '
        #r'(Prover|Verifier)/Param\. Kappa=(\d+),\s*Cols=(\d+),\s*B=(\d+),\s*L=(\d+),\s*B_small=(\d+),\s*K=(\d+)\s+'
        #r'(\d+\.\d+)\s+(\d+\.\d+±\d+\.\d+)(ms|s)'
def parse_line(line):
    # Updated pattern to match the provided lines
    pattern = (
             r'(Linearization|Decomposition|Folding) Godlilocks/(Linearization|Decomposition|Folding) '
             r'(Prover|Verifier)/Param\. Kappa=(\d+),\s*Cols=(\d+),\s*B=(\d+),\s*L=(\d+),\s*B_small=(\d+),\s*K=(\d+)\s+'
             r'(\d+\.\d+)\s+(\d+\.\d+±\d+\.\d+)(ms|s)'
    )

    match = re.search(pattern, line)

    if match:
        # Extract components from the match
        component = match.group(2)  # Second occurrence
        prover_or_verifier = match.group(3)  # Prover or Verifier
        kappa = int(match.group(4))
        cols = int(match.group(5))
        B = int(match.group(6))
        L = int(match.group(7))
        B_small = int(match.group(8))
        K = int(match.group(9))
        base = float(match.group(11).split('±')[0])
        time_measurement = match.group(11)  # Time with uncertainty (e.g., 6.3±0.19)
        time_unit = match.group(12)  # 's' or 'ms'

        # Convert to ms if unit is seconds
        if time_unit == 's':
            base *= 1000
        if time_unit == 'µs':
            base /= 1000

        result = [kappa, cols, B, L, B_small, K, "?", "?", "?", "?", "?", "?"]
        if component == "Linearization" and prover_or_verifier == "Prover":
            result[6] = base
        elif component == "Linearization" and prover_or_verifier == "Verifier":
            result[7] = base
        elif component == "Decomposition" and prover_or_verifier == "Prover":
            result[8] = base
        elif component == "Decomposition" and prover_or_verifier == "Verifier":
            result[9] = base
        elif component == "Folding" and prover_or_verifier == "Prover":
            result[10] = base
        elif component == "Folding" and prover_or_verifier == "Verifier":
            result[11] = base

        return result

    return None


def process_file(input_file, output_file):
    results = {}
    
    with open(input_file, 'r') as f:
        for line in f:
            result = parse_line(line)
            if result:
                key = tuple(result[:6])
                if key not in results:
                    results[key] = result
                else:
                    existing = results[key]
                    for i in range(6, 12):  # Updated range to cover all timing columns
                        if result[i] != "?":
                            existing[i] = result[i]
    
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Kappa', '|z|', 'B', 'L', 'B_small', 'K', 
                        'Linearization Prover', 'Linearization Verifier',
                        'Decomposition Prover', 'Decomposition Verifier',
                        'Folding Prover', 'Folding Verifier'])  # Added Folding columns
        writer.writerows(results.values())

# Usage
process_file('output.txt', 'results.csv')
