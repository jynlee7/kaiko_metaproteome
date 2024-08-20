# conda activate casanovo

DIRECTORY="/tahoma/emsl61033/mgf_data/casanovo_validation_set_individual_files"
MODEL_PATH="/tahoma/emsl61033/casanovo_checkpoints/epoch=1-step=160000.ckpt"
OUTPUT_DIR="/tahoma/emsl61033/casanovo_output_files/160000_trained"
FILES="$DIRECTORY/*.mgf"

if ls $FILES 1> /dev/null 2>&1; then
    for filepath in $FILES; do
        echo "Processing $filepath file..."
        # Extract the filename from the path
        mgf_filename=$(basename "$filepath")

        echo "Output file: $OUTPUT_DIR/$mgf_filename.mztab"

        # Run the casanovo command
        python -m casanovo.casanovo sequence "$filepath" --model="$MODEL_PATH" --output="$OUTPUT_DIR/$mgf_filename.mztab" --config="/tahoma/emsl61033/casanovo_checkpoints/casanovo_defaults_7_21_24/config.yaml"
    done
else
    echo "No .mgf files found in $DIRECTORY"
fi
