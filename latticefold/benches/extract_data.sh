#!/bin/bash

# Check if critcmp is installed
if ! command -v critcmp &> /dev/null
then
    echo "critcmp could not be found, please install it first."
    exit
fi

# Run critcmp and export results to results.json
critcmp --export base > results.json

# Run the extract_data.py script
python3 extract_json.py