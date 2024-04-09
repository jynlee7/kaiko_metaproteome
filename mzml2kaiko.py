### python mzml2kaiko.py --mzml_dir mint/ --out_dir mint/

from pyteomics import mzml,auxiliary
import time
import gzip
import glob
import os
import sys

import argparse

########################################################
parser = argparse.ArgumentParser()
parser.add_argument(
    '--mzml_dir', type=str,
    help='mzML directory')
parser.add_argument(
    '--out_dir', type=str,
    help='output directory')
parser.add_argument(
    '--gz', action='store_true',
    help='mzML.gz?')

FLAGS = parser.parse_args()
########################################################

def inspect_mzML_file(fpath, gzipped=True):
    spectra = []
    if gzipped:
        f = gzip.open(fpath, 'rb')
    else:
        f = fpath
    for obj in mzml.read(f):
        spectra.append(obj)
    if gzipped: f.close()
    return spectra

# start_time = time.time()

# spectra = inspect_mzML_file('mint/MinT_frac_Kans_2D_01_01_2D_27Feb17_Tiger_16-09-25.mzML', False)
# end_time = time.time()
# print('Num:{0}, time:{1}'.format(len(spectra), end_time-start_time))

def generate_mgf_without_annotation(mzml_spectra, file_index=0, ntops=500, out_file='out.mgf'):
    num_spectra = 0
    with open(out_file, 'w') as f:
        for spectrum in mzml_spectra:
            if spectrum['ms level'] != 2:
                continue
            scan = int(spectrum['id'].split('scan=')[1])
            try:
            
                mz_arr = spectrum['m/z array']
                int_arr = spectrum['intensity array']
                rtsec = 60.0*(spectrum['scanList']['scan'][0]['scan start time'])
                selectedIon = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]

                assert len(mz_arr) == len(int_arr), "[ERR] Wrong data format: len(mz_arr) != len(int_arr)"

                print("BEGIN IONS", file=f)
                print("TITLE={0}.{1}".format(file_index, scan), file=f)
                print("PEPMASS={0}".format(selectedIon['selected ion m/z']), file=f)
                # sometimes they don't have a charge info. if so, we use the annotation file
                if 'charge state' in selectedIon:
                    print("CHARGE={0:d}+".format(int(selectedIon['charge state'])), file=f)
                else:
                    print("CHARGE={0:d}+".format(999), file=f)
                print("SCANS={0}:{1}".format(file_index, scan), file=f)
                print("RTINSECONDS={0}".format(rtsec), file=f)
                print("SEQ=UNKNOWN", file=f)
                for i in range(len(mz_arr)):
                    print("{0} {1}".format(mz_arr[i], int_arr[i]), file=f)
                print("END IONS", file=f)
                num_spectra += 1
            except:
                print('[ERR]', scan, spectrum)
                continue
                    
        return num_spectra

def generate_mgf_files(data_dir, dest_dir='./', gzipped=True):
    # collect mzML.gz files
    if gzipped:
        mzML_files = glob.glob(data_dir + "/*.mzML.gz")
    else:
        mzML_files = glob.glob(data_dir + "/*.mzML")
    
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)

    mzML_log_handler = open(dest_dir + '/mgf_list.log', 'w+')
    print("id\tmgf_file\tnum_scans\ttotal_scans", file=mzML_log_handler)
    
    start_time = time.time()
    num_mzML_files = len(mzML_files)
    total_scans = 0
    for i, mzML_file in enumerate(mzML_files):
        try:
            if gzipped:
                common_name = os.path.basename(mzML_file).rsplit('.mzML.gz')[0]
            else:
                common_name = os.path.basename(mzML_file).rsplit('.mzML')[0]
        
            if os.path.exists(dest_dir + '/' + common_name + '.mgf'):
                print('[{0:3d}/{1:3d}] {2}, Already exists' \
                      .format(i+1,
                              num_mzML_files,
                              common_name))
                continue
                
        
            seq_file = glob.glob(data_dir + '/' + common_name + "*.txt")
            msg = ""
            scan_ids = []
            num_spectra = 0
            mzml_spectra = inspect_mzML_file(mzML_file, gzipped)
            
            if len(seq_file) == 1:    
                annotated = get_annotated_pepseq(seq_file[0])
                scan_ids = list(annotated.Scan)
                pepseqs = list(annotated.pepseq)
                charges = list(annotated.Charge)
                num_scans = len(scan_ids)
                num_spectra = generate_mgf(mzml_spectra,
                                           scan_ids,
                                           pepseqs,
                                           charges,
                                           file_index = i,
                                           out_file=dest_dir + '/' + common_name + '.mgf')
            else:
                num_spectra = generate_mgf_without_annotation(mzml_spectra,
                                                              file_index=i,
                                                              out_file=dest_dir + '/' + common_name + '.mgf')
                num_scans = num_spectra
            total_scans += num_spectra
            msg = "SUCCESS"
            print('[{0:3d}/{1:3d}] {2}, {3:d}/{4:d}/{5:d}, {6:.2f}sec' \
                      .format(i+1,
                              num_mzML_files,
                              common_name,
                              num_spectra,
                              num_scans,
                              total_scans,
                              time.time()-start_time))
            print("{0}\t{1}\t{2}\t{3}".format(i, common_name, num_spectra, total_scans), file=mzML_log_handler)
            sys.stdout.flush()
        except Exception as e:
            print('[ERR] {}'.format(mzML_file))
            print('[ERR]', e)
        

if __name__ == "__main__":
    generate_mgf_files(FLAGS.mzml_dir, FLAGS.out_dir, FLAGS.gz)
