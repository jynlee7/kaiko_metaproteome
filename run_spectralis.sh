DIRECTORY="/tahoma/emsl61033/mgf_data/casanovo_validation_set_individual_files"
MODEL_PATH="/tahoma/emsl61033/casanovo_checkpoints/epoch=1-step=160000.ckpt"
OUTPUT_DIR="/tahoma/emsl61033/spectralis_output"
FILES="$DIRECTORY/*.mgf"

if ls $FILES 1> /dev/null 2>&1; then
    for filepath in $FILES; do
        echo "Processing $filepath file..."
        # Extract the filename from the path
        spec_filename=$(basename "$filepath")

        echo "Output file: $OUTPUT_DIR/$spec_filename.csv"

        # Run the casanovo command
        python -m spectralis.spectralis_master --mode="ea" --input_path="$filepath" --output="$OUTPUT_DIR/$spec_filename.csv" --config="/home/leej741/spectralis/spectralis_config.yaml"
    done
else
    echo "No .mgf files found in $DIRECTORY"
fi
