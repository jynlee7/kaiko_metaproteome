from pyteomics import mgf
import pandas as pd
import re
# import os

mgf_file_path = "/Users/leej741/Desktop/validation_set/testing/01_ECTABPP_NP_1_9Jan17_Wally_17-12-02.mzML_0.001_qvaluecutoff.mgf"
mgf_file_name = re.sub(r'^.*/|\.mgf$', '', mgf_file_path)
tsv_file_path = "/Users/leej741/Desktop/git/comparison.tsv"
output_mgf_file_path = "/Users/leej741/Desktop/validation_set/testing/testing_output.mgf"

def check_mgf_files_in_tsv(mgf_file_name, tsv_file_path):
    tsv_data = pd.read_csv(tsv_file_path, sep='\t')
    return tsv_data[tsv_data['mgf'].str.contains(mgf_file_name, na=False)]

def modify_sequence(sequence):
    if isinstance(sequence, float) and pd.isna(sequence):
        print(f"Skipping modification for NaN sequence: {sequence}")
        return ''
    
    pattern = r'\+\d+(\.\d+)?'  
    modified_sequence = re.sub(pattern, '(ox)', sequence)
    return modified_sequence

def create_new_mgf_with_predicted_sequences(mgf_file_name, tsv_file_path, output_mgf_file_path):
    mgf_matches = check_mgf_files_in_tsv(mgf_file_name, tsv_file_path)
    if mgf_matches.empty:
        print(f"No matching sequences found for {mgf_file_name} in {tsv_file_path}")
        return
    
    scan_to_sequence = mgf_matches.set_index('scans')['casanovo_seq'].to_dict()
    
    spectra = []
    with mgf.read(mgf_file_path) as reader:
        for spectrum in reader:
            original_seq = scan_to_sequence.get(int(spectrum['params']['title']), '')
            if not original_seq:
                print(f"No sequence found for scan {spectrum['params']['title']}")
            else:
                modified_seq = modify_sequence(original_seq)
                if len(modified_seq) == 0:
                    print(f"Skipping empty or invalid sequence for scan {spectrum['params']['title']}")
                    continue
                spectrum['params']['seq'] = modified_seq
            spectra.append(spectrum)

    with open(output_mgf_file_path, 'w') as output_file:
        for spectrum in spectra:
            output_file.write("BEGIN IONS\n")
            output_file.write(f"TITLE={spectrum['params'].get('TITLE', spectrum['params'].get('scans', ''))}\n")
            output_file.write(f"PEPMASS={spectrum['params']['pepmass'][0]}\n")
            output_file.write(f"CHARGE={spectrum['params']['charge']}\n")
            output_file.write(f"SCANS={spectrum['params']['scans']}\n")
            output_file.write(f"RTINSECONDS={spectrum['params'].get('rtinseconds', '')}\n")
            if 'seq' in spectrum['params']:
                output_file.write(f"SEQ={spectrum['params']['seq']}\n")
            for mz, intensity in zip(spectrum['m/z array'], spectrum['intensity array']):
                output_file.write(f"{mz} {intensity}\n")
            output_file.write("END IONS\n")
    
    print(f"New MGF file created with predicted sequences at {output_mgf_file_path}")

create_new_mgf_with_predicted_sequences(mgf_file_name, tsv_file_path, output_mgf_file_path)