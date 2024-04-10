FROM camiloposso15/tensorflow2.12.0-py310

WORKDIR /Kaiko_metaproteome
ADD ./Kaiko_denovo /Kaiko_metaproteome/Kaiko_denovo
ADD ./Kaiko_2.py /Kaiko_metaproteome
ADD ./Kaiko_3.py /Kaiko_metaproteome
ADD ./Kaiko_4.py /Kaiko_metaproteome
ADD ./Kaiko_parse_uniref.py /Kaiko_metaproteome
ADD ./Kaiko_pipeline_main.py /Kaiko_metaproteome
ADD ./kaiko_unit_test.py /Kaiko_metaproteome
ADD ./unit_test_util.py /Kaiko_metaproteome
ADD ./kaiko_defaults.yaml /Kaiko_metaproteome

RUN pip install --no-cache-dir --upgrade pip && \ 
    pip install --no-cache-dir indexed_gzip==1.8.5 llvmlite==0.39.1 biopython==1.81 numba==0.56.4 pyteomics sigopt==3.2.0 memory-profiler pyyaml pathlib s3path pyodbc openpyxl xlsxwriter
