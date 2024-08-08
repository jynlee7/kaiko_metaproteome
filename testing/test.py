from pyteomics import mztab, mgf
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import difflib
import re
import pyarrow.parquet as pq
import pyarrow as pa
import time
import glob
import os

mgf_files = glob.glob("/Users/leej741/Desktop/validation_set/mgf_individual_files/*.mgf")
casanovo_files = glob.glob("/Users/leej741/Desktop/validation_set/casanovo_output_files/*.mztab")
output_path = "/Users/leej741/Desktop/git/comparison.tsv"


def get_mgf_path_by_filename(mgf_paths):
    mgf_dict = {}
    for mgf_path in mgf_paths:
        mgf_dict[os.path.splitext(os.path.basename(mgf_path))[0]] = mgf_path
    return mgf_dict

def get_casanovo_path_by_mgf_filename(casanovo_paths):
    casanovo_dict = {}
    for casanovo_path in casanovo_paths:
        tables = mztab.MzTab(casanovo_path)
        mgf_name_path = tables.metadata['ms_run[1]-location'][8:]
        mgf_name = os.path.splitext(os.path.basename(mgf_name_path))[0]
        casanovo_dict[mgf_name] = casanovo_path
    return casanovo_dict


def get_common_paths(mgf_dict, casanovo_dict):
    common_paths = []
    for mgf_name in mgf_dict.keys():
        if mgf_name in mgf_dict and mgf_name in casanovo_dict:
            common_paths.append((
                mgf_dict[mgf_name],
                casanovo_dict[mgf_name]
            ))
    return common_paths

def calc_length(seq):
    if '+' in seq or any(char.isdigit() for char in seq):
        count = sum(1 for char in seq if char.isupper())
        return count
    else:
        return len(seq)

def forgive_IL_mismatch(true_seq, casanovo_seq):
    """
    Forgive mismatches between Isoleucine (I) and Leucine (L).
    """
    if len(true_seq) != len(casanovo_seq):
        return False

    for t, p in zip(true_seq, casanovo_seq):
        if t == p:
            continue
        if (t == 'I' and p == 'L') or (t == 'L' and p == 'I'):
            continue
        return False
    return True

def get_score(true_seq, casanovo_seq, true_pepmass, casanovo_pepmass):
    """
    Calculate the score forgiving I/L mismatches and comparing pepmasses.
    """
    if forgive_IL_mismatch(true_seq, casanovo_seq):
        return 1.0 if abs(true_pepmass - casanovo_pepmass) < 0.001 else 0.0
    return difflib.SequenceMatcher(None, true_seq, casanovo_seq).ratio()

def get_lenient_match(true_seq, casanovo_seq, true_pepmass, casanovo_pepmass):
    if get_score(true_seq, casanovo_seq, true_pepmass, casanovo_pepmass) == 1.0:
        return True
    else:
        return true_seq == casanovo_seq

mgf_paths = get_mgf_path_by_filename(mgf_files)
casanovo_paths = get_casanovo_path_by_mgf_filename(casanovo_files)
common_paths = get_common_paths(mgf_paths, casanovo_paths)

def aggregate_mgf_casanovo(mgf_path, casanovo_path, output_path):
    """
    read a mgf file and create a dataframe for it
    """
    print(f"Reading mgf file...: {mgf_path}")
    
    start_time = time.perf_counter()
    mgf_data = mgf.read(mgf_path, use_index=True)

    # extract the precursor info from mgf_data as a dataframe
    precursors = [d['params'] for d in mgf_data]
    df = pd.DataFrame.from_records(precursors)
    
    # select the first element in a pepmass tuple
    df['pepmass'] = df['pepmass'].apply(lambda x: x[0])
    # get a integer for a charge value from the charge list
    df['charge'] = df['charge'].apply(lambda x: int(x[0]))

    print(f"Finish reading the mgf file: {time.perf_counter() - start_time:.3f} seconds")

    """
    read a mztab file and create a dataframe for it
    """
    tables = mztab.MzTab(casanovo_path)
    psms = tables.spectrum_match_table
    psms['scan_index'] = psms["spectra_ref"].str.extract("index=(\\d+)")
    # convert the scan index as an object to integer
    psms['scan_index'] = psms['scan_index'].astype(int)
    # create the mztab dataframe with scan index and sequence
    casanovo_df = psms[['scan_index', 'sequence']]
    # rename the column: sequnce to casanovo 
    casanovo_df = casanovo_df.rename(columns={'sequence': 'casanovo_seq'})
    casanovo_df = casanovo_df.set_index('scan_index')

    # merge the mgf df and casanovo df with index number
    result = pd.concat([df, casanovo_df], axis=1)
    # get the specific columns that are necessary, remove title
    del result['title']
    result = result.rename(
        columns={
            'seq': 'true_seq'
        }
    )
    result['true_seq'] = result['true_seq'].astype(str)
    result['casanovo_seq'] = result['casanovo_seq'].astype(str)

    result['score'] = result.apply(lambda x: get_score(x.true_seq, x.casanovo_seq, x.pepmass, x.pepmass), axis=1)
    result['lenient_match'] = result.apply(lambda x: get_lenient_match(x.true_seq, x.casanovo_seq, x.pepmass, x.pepmass), axis=1)
    result['exact_match'] = result.apply(lambda x: x.true_seq == x.casanovo_seq, axis=1)
    result["peptide_length"] = result.true_seq.apply(calc_length)

    parquet_file_path = mgf_path + '.parquet'
    print(f"Reading parquet file...: {mgf_path}")
    start_time = time.perf_counter()
    if not os.path.exists(parquet_file_path):
        result.to_parquet(parquet_file_path, engine='pyarrow')
    table = pq.read_table(parquet_file_path)
    df = table.to_pandas()
    print(f"Finish reading the parquet file: {time.perf_counter() - start_time:.3f} seconds")

    return result

def get_filename_without_extension(file_path, extension):
    pattern = rf'[^/]+(?=\.{extension})'
    match = re.search(pattern, file_path)
    return match.group(0) if match else None

def compare_multiple(common_paths):
    dfs = []
    for (mgf_path, casanovo_path) in common_paths:
        temp_df = aggregate_mgf_casanovo(mgf_path, casanovo_path, output_path)
        temp_df['mgf'] = get_filename_without_extension(mgf_path, 'mgf')
        # temp_df['casanovo'] = get_filename_without_extension(casanovo_path, 'mztab')
        dfs.append(temp_df)
    compare_df = pd.concat(dfs, ignore_index=True)
    return compare_df.to_csv(output_path, index=None, sep='\t')


def graph_and_create_table(tsv_path):
    # Read the data
    result = pd.read_csv(tsv_path, sep='\t')
    
    # Group by peptide length and calculate the mean of lenient matches
    grouped = result.groupby('peptide_length')['lenient_match'].mean() * 100

    # Find specific values for certain peptide lengths (e.g., 5, 6, 7)
    lengths_of_interest = [i for i in range(6, 51)]
    values = {length: grouped.get(length, None) for length in lengths_of_interest}

    # Create a table (DataFrame) for specific values vs. length of the peptide
    table = pd.DataFrame(list(values.items()), columns=['Peptide Length', 'Lenient Match Percentage'])
    
    # Print the table
    print("Table of Specific Values vs. Length of the Peptide:")
    print(table.to_string(index=False))

    # Plot the data using seaborn or matplotlib
    grouped_df = grouped.reset_index()

    plt.figure(figsize=(10, 6))
    sns.barplot(data=grouped_df, x='peptide_length', y='lenient_match', hue='peptide_length', palette='viridis', dodge=False, legend=False)

    # Set the title and labels
    plt.title('Lenient Match Percentage vs. Peptide Length')
    plt.xlabel('Peptide Length')
    plt.ylabel('Lenient Match Percentage (%)')

    # Display the plot
    plt.show()

compare_multiple(common_paths)
# graph_and_create_table(output_path)

exit()

def mgf_to_dataframe(mgf_path):
    """
    read a mgf file and create a dataframe for it
    """
    print(f"Reading mgf file...: {mgf_path}")
    
    start_time = time.perf_counter()
    mgf_data = mgf.read(mgf_path, use_index=True)

    # extract the precursor info from mgf_data as a dataframe
    precursors = [d['params'] for d in mgf_data]
    df = pd.DataFrame.from_records(precursors)
    
    # select the first element in a pepmass tuple
    df['pepmass'] = df['pepmass'].apply(lambda x: x[0])
    # get a integer for a charge value from the charge list
    df['charge'] = df['charge'].apply(lambda x: int(x[0]))

    print(f"Finish reading the mgf file: {time.perf_counter() - start_time:.3f} seconds")
    return df


# mgf_path = "/Users/leej741/Desktop/camilo/HA_Biocide_T3_06_Run2_1Sep17_Oak_Jup-17-09-01.mzML_0.001_qvaluecutoff.mgf"
# df = mgf_to_dataframe(mgf_path)

# # store a pyarrow.parquet file
# print(f"Saving mgf data to a parquet file...")
# table = pa.Table.from_pandas(df)
# pq.write_table(table, mgf_path + '.parquet')
# print("Done")


# read a parquet and create a pandas df
mgf_path = "/Users/leej741/Desktop/camilo/HA_Biocide_T3_06_Run2_1Sep17_Oak_Jup-17-09-01.mzML_0.001_qvaluecutoff.mgf"
print(f"Reading parquet file...: {mgf_path}")
start_time = time.perf_counter()
table = pq.read_table(mgf_path + '.parquet')
df = table.to_pandas()
print(f"Finish reading the parquet file: {time.perf_counter() - start_time:.3f} seconds")


casanovo_path = "/Users/leej741/Desktop/git/casanovo_20240708085126.mztab"
tables = mztab.MzTab(casanovo_path)
psms = tables.spectrum_match_table
psms['scan_index'] = psms["spectra_ref"].str.extract("index=(\\d+)")
# convert the scan index as an object to integer
psms['scan_index'] = psms['scan_index'].astype(int)
# create the mztab dataframe with scan index and sequence
casanovo_df = psms[['scan_index', 'sequence']]
# rename the column: sequnce to casanovo 
casanovo_df = casanovo_df.rename(columns={'sequence': 'casanovo_seq'})
casanovo_df = casanovo_df.set_index('scan_index')

# merge the mgf df and casanovo df with index number
result = pd.concat([df, casanovo_df], axis=1)

# merge the merged df with the casanovo df with scan index
# comparison_df = result.merge(casanovo_df, left_on="index", right_on="scan_index")
# get the specific columns that are necessary, remove title
del result['title']
result = result.rename(
    columns={
        'seq': 'true_seq'
    }
)
# rename the columns to specific names
# comparison_df = comparison_df.rename(
#     columns={
#         'scans': 'scan_num',
#         'output_seq': 'kaiko_seq',
#         'casanovo_sequence': 'casanovo_seq'
#     }
# )
# def calc_length(seq):
#     if '+' in seq or any(char.isdigit() for char in seq):
#         count = sum(1 for char in seq if char.isupper())
#         return count
#     else:
#         return len(seq)

# def get_score(true_seq, casanovo_seq):
#     comp_score = difflib.SequenceMatcher(None, true_seq, casanovo_seq).ratio()
#     return comp_score


# result['score'] = result.apply(lambda x: get_score(x.true_seq, x.casanovo_seq), axis=1)
# result['exact_match'] = result.apply(lambda x: x.true_seq == x.casanovo_seq, axis=1)
# result["peptide_length"] = result.true_seq.apply(calc_length, axis=1)

# res_graph = result.plot(kind='area',
#             x='seq',
#             y='score',
#             color='red')
# res_graph.title("original vs trained")
# plt.pyplot.show()
pattern = r'[^/]+(?=\.mgf)'
match = re.search(pattern, mgf_path)
mgf_filename = match.group(0) if match else None

mgf_filenames = [mgf_filename for i in range(len(result))]
result["dataset"] = mgf_filenames

# Group by peptide length and calculate the exact match percentage
grouped = result.groupby('peptide_length')['exact_match'].mean() * 100

# Reset the index to turn it back into a DataFrame
grouped_df = grouped.reset_index()

# Plot the data using seaborn or matplotlib
plt.figure(figsize=(10, 6))
sns.barplot(data=grouped_df, x='peptide_length', y='exact_match', hue='peptide_length', palette='viridis', dodge=False, legend=False)

# Set the title and labels
plt.title('Exact Match Percentage vs. Peptide Length')
plt.xlabel('Peptide Length')
plt.ylabel('Exact Match Percentage (%)')

# Display the plot
plt.show()

print(result)
# result.to_csv(f"/Users/leej741/Desktop/{mgf_filename}.csv", sep="\t")