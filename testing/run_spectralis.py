import glob
import re
from spectralis.spectralis_master import Spectralis

spectralis = Spectralis(config_path="../spectralis/spectralis_config.yaml")

MGF_FILES = glob.glob("/Users/leej741/Desktop/validation_set/testing/*.mgf")


for mgf_file in MGF_FILES:
    mgf_file_name = re.sub(r'^.*/|\.mgf$', '', mgf_file)
    spectralis.rescoring_from_mgf(mgf_path=mgf_file, out_path=f"/Users/leej741/Desktop/validation_set/spectralis_out/{mgf_file_name}.csv")
exit()