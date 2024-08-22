# DIRECTORY="/tahoma/emsl61033/Spectralis/spec_input/testing"
# OUTPUT_DIR="/tahoma/emsl61033/spectralis_output"
# FILES="$DIRECTORY/*.mgf"

# if ls $FILES 1> /dev/null 2>&1; then
#     for filepath in $FILES; do
#         echo "Processing $filepath file..."
#         # Extract the filename from the path
#         spec_filename=$(basename "$filepath")

#         echo "Output file: $OUTPUT_DIR/$spec_filename.csv"

#         # Run the casanovo command
#         python -m spectralis.spectralis_master --mode="ea" --input_path="$filepath" --output="$OUTPUT_DIR/$spec_filename.csv" --config="/home/leej741/spectralis/spectralis_config.yaml"
#     done
# else
#     echo "No .mgf files found in $DIRECTORY"
# fi
#!/bin/bash

# Set the paths
CONFIG_PATH="/home/leej741/spectralis/spectralis_config.yaml"
MGF_DIR="/tahoma/emsl61033/Spectralis/spec_input/testing"
OUT_DIR="/tahoma/emsl61033/spectralis_output"

# Activate the virtual environment if needed
# source /path/to/your/venv/bin/activate

# Iterate over each .mgf file in the directory
for mgf_file in "$MGF_DIR"/*.mgf; do
    # Extract the file name without the directory and extension
    mgf_file_name=$(basename "$mgf_file" .mgf)
    
    # Run the Spectralis rescoring command using Python
    python3 -c "
import re
from spectralis.spectralis_master import Spectralis

spectralis = Spectralis(config_path='$CONFIG_PATH')
spectralis.rescoring_from_mgf(mgf_path='$mgf_file', out_path='$OUT_DIR/$mgf_file_name.csv')
"
done

# Deactivate the virtual environment if needed
# deactivate
