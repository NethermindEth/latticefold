import re
import csv

def parse_line(line):
    # Updated pattern to capture the time unit (s or ms)
    pattern = r'(?:Linearization|Decomposition|Folding) \w+ with decompose witness/(?:(?:Linearization|Decomposition|Folding) (?:Prover|Verifier))/Param\. Kappa=(\d+),\s*Cols=(\d+),\s*B=(\d+),\s*L=(\d+),\s*B_small=(\d+),\s*K=(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)±\d+\.\d+(m?s)'
    match = re.search(pattern, line)
    
    if match:
        component = re.search(r'with decompose witness/(.*?)/', line).group(1)
        kappa = int(match.group(1))
        cols = int(match.group(2))
        B = int(match.group(3))
        L = int(match.group(4))
        B_small = int(match.group(5))
        K = int(match.group(6))
        base = float(match.group(8))
        time_unit = match.group(9)  # 's' or 'ms'
        
        # Convert to milliseconds if the unit is seconds // Beware of other units
        if time_unit == 's':
            base *= 1000
        
        result = [kappa, cols, B, L, B_small, K, "?", "?", "?", "?", "?", "?"]
        
        # Update the appropriate column based on which component we found
        if component == "Linearization Prover":
            result[6] = base
        elif component == "Linearization Verifier":
            result[7] = base
        elif component == "Decomposition Prover":
            result[8] = base
        elif component == "Decomposition Verifier":
            result[9] = base
        elif component == "Folding Prover":
            result[10] = base
        elif component == "Folding Verifier":
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
