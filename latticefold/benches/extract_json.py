import json
import csv
import os

def parse_value_str(value_str):
    """Parse the value_str to extract Kappa, W_CCS, W, L, B, B_SMALL, and K."""
    params = {}
    
    for param in value_str.split(', '):
        key, value = param.split('=')
        
        key = key.strip().upper()
        
        if key.startswith('PARAM. '):
            key = key.replace('PARAM. ', '')
        if key == 'COLS':
            key = 'W_CCS'
        
        params[key] = int(value)
    
    if 'W' not in params and 'W_CCS' in params and 'L' in params:
        params['W'] = params['W_CCS'] * params['L']
    
    return params

def create_csv_files_from_json(json_file):
    with open(json_file, 'r') as file:
        data = json.load(file)

    # Access the 'benchmarks' key in the JSON
    benchmarks = data.get('benchmarks', {})

    # Create a directory to store the CSV files
    output_dir = 'csv_output'
    os.makedirs(output_dir, exist_ok=True)

    # Dictionary to hold data for each group_id
    group_data = {}
    fieldnames_set = set()  # Set to track added fieldnames

    # List to hold skipped entries
    skipped_entries = []

    for key, value in benchmarks.items():  # Iterate over the benchmarks
        group_id = value['criterion_benchmark_v1']['group_id']
        function_id = value['criterion_benchmark_v1']['function_id']
        value_str = value['criterion_benchmark_v1'].get('value_str', '')

        # Skip entries without the required parameters
        if not value_str:
            print(f"Skipping entry {key} due to missing value_str")
            skipped_entries.append({
                'function_id': function_id,
                'value_str': value_str,
                'point_estimate': value['criterion_estimates_v1']['mean']['point_estimate'] / 1_000_000
            })
            continue

        # Get mean value and convert it to ms
        point_estimate = value['criterion_estimates_v1']['mean']['point_estimate'] / 1_000_000

        # Parse the value_str to get the parameters
        params = parse_value_str(value_str)

        # Check if all required parameters are present
        required_params = ['KAPPA', 'W_CCS', 'W', 'L', 'B', 'B_SMALL', 'K']
        if not all(param in params for param in required_params):
            print(f"Skipping entry {key} due to missing required parameters: {required_params}")
            skipped_entries.append({
                'function_id': function_id,
                'value_str': value_str,
                'point_estimate': point_estimate
            })
            continue

        # Add function_id as a column if not already added
        if function_id not in fieldnames_set:
            fieldnames_set.add(function_id)

        # Initialize group_id in group_data if not already present
        if group_id not in group_data:
            group_data[group_id] = []

        # Check if a row with the same parameters already exists
        existing_entry = next((entry for entry in group_data[group_id] if all(entry[param] == params[param] for param in required_params)), None)

        if existing_entry:
            # Update the existing entry with the new point_estimate
            existing_entry[function_id] = point_estimate
        else:
            # Add a new entry with the point_estimate
            params[function_id] = point_estimate
            group_data[group_id].append(params)

    # Write each group_id's data to a separate CSV file
    for group_id, entries in group_data.items():
        # Determine which columns have all zero values
        zero_columns = set()
        for field in required_params + sorted(fieldnames_set):
            if all(entry.get(field, 0) == 0 for entry in entries):
                zero_columns.add(field)

        # Filter out zero columns from fieldnames
        fieldnames = [field for field in required_params + sorted(fieldnames_set) if field not in zero_columns]

        # Append "(ms)" to each fieldname except for the specified ones
        fieldnames_with_units = [
            f"{field} (ms)" if field not in ['KAPPA', 'W_CCS', 'W', 'L', 'B', 'B_SMALL', 'K'] else field
            for field in fieldnames
        ]

        csv_file = os.path.join(output_dir, f"{group_id.replace(' ', '_')}.csv")
        with open(csv_file, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames_with_units)

            writer.writeheader()
            for entry in entries:
                # Use get method to avoid KeyError
                filtered_entry = {
                    (f"{field} (ms)" if field not in ['KAPPA', 'W_CCS', 'W', 'L', 'B', 'B_SMALL', 'K'] else field): entry.get(field, 'N/A')
                    for field in fieldnames
                }
                writer.writerow(filtered_entry)

    # Write skipped entries to a separate CSV file
    skipped_csv_file = os.path.join(output_dir, "skipped_entries.csv")
    with open(skipped_csv_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['function_id', 'value_str', 'point_estimate'])
        writer.writeheader()
        for entry in skipped_entries:
            writer.writerow(entry)

# Usage
create_csv_files_from_json('results.json')
