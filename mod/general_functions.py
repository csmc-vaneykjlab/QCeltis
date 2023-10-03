"""
General Functions used for extracting ID-Based and ID-Free Metrics
"""

import pandas as pd
import numpy as np
import os
import sys
import logging
import argparse
from collections import Counter

#-------------------------------------------------------------------------- VARIABLES --------------------------------------------------------------------------

#50 unique colors
color_list = [
    "#1f77b4", "#2ca02c", "#d62728", "#ff7f0e", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
    "#393b79", "#e7ba52", "#ad494a", "#8c6d31", "#7b4173", "#843c39", "#d6616b", "#7b4173", "#ce6dbd", "#c6c6c6",
    "#cedb9c", "#9c9ede", "#5254a3", "#e6ab02", "#bd9e39", "#6b6ecf", "#6b6ecf", "#b5cf6b", "#7fc97f", "#e7cb94",
    "#e7969c", "#fddbc7", "#c7eae5", "#9edae5", "#fdae6b", "#fdd0a2", "#636363", "#989898", "#bdbdbd", "#f7f7f7"
]

#-------------------------------------------------------------------------- FUNCTIONS ---------------------------------------------------------------------------

def check_path(filepath):
    if not os.path.exists(str(filepath)):
        raise FileNotFoundError(f"{filepath} doesn't exist, please specify another output path")
    return None

def check_file(txtfile, level):

    #adding a check if the input file is tab delimited or not.
    if not is_tab_delimited(txtfile):
        print(f"The {txtfile} is not tab delimited, Please check the input")
        logging.error(f"The {txtfile} is not tab delimited, Please check the input")
        sys.exit(1)

    df = pd.read_csv(txtfile, sep="\t")

    if level == "Protein":
        if "Protein" not in df.columns.tolist():
            print(f"'Protein' column doesn't exist in the given file {txtfile}")
            logging.error(f"'Protein' column doesn't exist in the given file {txtfile}")
            sys.exit(1)

        cols = df.columns.tolist()
        cols.remove("Protein")

        if len(cols) < 2:
            print(f"The number of samples provided in {txtfile} is less than 2")
            logging.error(f"The number of samples provided in {txtfile} is less than 2")
            sys.exit(1)

        logging.info(f"{len(cols)} samples will be analyzed from {txtfile}")

    if level == "Peptide":
        if not all(ele in df.columns.tolist() for ele in ['Protein','Peptide']):
            print(f"'Protein' or 'Peptide' column doesn't exist in the given file {txtfile}")
            logging.error(f"'Protein' or 'Peptide' column doesn't exist in the given file {txtfile}")
            sys.exit(1)

        cols = df.columns.tolist()
        cols.remove("Protein")
        cols.remove("Peptide")

        if len(cols) < 2:
            print(f"The number of samples provided in {txtfile} is less than 2")
            logging.error(f"The number of samples provided in {txtfile} is less than 2")
            sys.exit(1)

        logging.info(f"{len(cols)} samples will be analyzed from {txtfile}")

    if level == "Precursor":
        if not all(ele in df.columns.tolist() for ele in ['Protein','Peptide','Precursor']):
            print(f"'Protein','Peptide' or 'Precursor' column doesn't exist in the given file {txtfile}")
            logging.error(f"'Protein','Peptide' or 'Precursor' column doesn't exist in the given file {txtfile}")
            sys.exit(1)

        cols = df.columns.tolist()
        cols.remove("Protein")
        cols.remove("Peptide")
        cols.remove("Precursor")

        if len(cols) < 2:
            print(f"The number of samples provided in {txtfile} is less than 2")
            logging.error(f"The number of samples provided in {txtfile} is less than 2")
            sys.exit(1)

        logging.info(f"{len(cols)} samples will be analyzed from {txtfile}")

    if level == "Peptide List":
        if "Peptide" not in df.columns.tolist():
            print(f"'Peptide' column doesn't exist in the given file: {txtfile}")
            logging.error(f"'Peptide' column doesn't exist in the given file: {txtfile}")
            sys.exit(1)

    return None

def check_grouping_file(txtfile, grouping_file):

    logging.info(f"Checking provided grouping file: {grouping_file}")

    #adding a check if the given input file is tab delimited or not.
    if not is_tab_delimited(grouping_file):
        print(f"The {grouping_file} is not tab delimited, Please check the input")
        logging.error(f"The {grouping_file} is not tab delimited, Please check the input")
        sys.exit(1)

    df = pd.read_csv(txtfile, sep="\t")

    #add error check for group df
    group_df = pd.read_csv(grouping_file, sep="\t")

    if not all(ele in group_df.columns.tolist() for ele in ['Filename','Group']):
        print(f"'Filename' or 'Group' column not present in {grouping_file}")
        logging.error(f"'Filename' or 'Group' column not present in {grouping_file}")
        sys.exit(1)

    #checking if samples between given protein/peptide/precursor file exist in grouping file
    expected_columns = ['Protein', 'Peptide', 'Precursor']
    samples = [x for x in group_df['Filename'].tolist() if x not in expected_columns]

    samples_not_present_in_grouping = []

    for filename in group_df['Filename'].tolist():
        if filename not in samples:
            samples_not_present_in_grouping.append(filename)

    if not len(samples_not_present_in_grouping) == 0:
        samples_string = ", ".join(samples_not_present_in_grouping)
        print(f"{samples_string} are not present in the given grouping file: {grouping_file}")
        logging.error(f"{samples_string} are not present in the given grouping file: {grouping_file}")
        sys.exit(1)

    groups_from_df = list(set(group_df['Group'].tolist()))

    if len(groups_from_df) < 2:
        print("At least 2 groups should be provided")
        logging.error("At least 2 groups should be provided")
        sys.exit(1)

    for group in groups_from_df:
        group_files = group_df[group_df['Group'] == group]['Filename'].tolist()
        if len(group_files) < 2:
            logging.error("For each group, atleast 2 filenames should be present")
            sys.exit(1)

    groups = [str(group) for group in groups_from_df]

    return ", ".join(groups)

def check_samples(mzml_dir, protein_level, peptide_level, precursor_level, grouping_file):

    filename_lists = []

    #getting filenames from mzml_dir
    if mzml_dir:
        mzml_filenames = []
        list_dir = os.listdir(mzml_dir)
        for mzml_file in list_dir:
            if not mzml_file.endswith(".mzML"):
                continue
            mzml_filenames.append(mzml_file)
        filename_lists.append(mzml_filenames)


    #getting filenames from protein file
    if protein_level:
        pt_level = pd.read_csv(protein_level, sep="\t")
        pt_filenames = pt_level.columns.tolist()
        pt_filenames.remove('Protein')
        filename_lists.append(pt_filenames)

    #getting filenames from peptide file
    if peptide_level:
        pep_level = pd.read_csv(peptide_level, sep="\t")
        pep_filenames = pep_level.columns.tolist()
        pep_filenames.remove('Protein')
        pep_filenames.remove('Peptide')
        filename_lists.append(pep_filenames)

    #getting filenames from precursor file
    if precursor_level:
        pre_level = pd.read_csv(precursor_level, sep="\t")
        pre_filenames = pre_level.columns.tolist()
        pre_filenames.remove('Protein')
        pre_filenames.remove('Peptide')
        pre_filenames.remove('Precursor')
        filename_lists.append(pre_filenames)

    #getting filenames from grouping file
    if grouping_file:
        groupdf = pd.read_csv(grouping_file, sep="\t")
        group_filenames = groupdf['Filename'].tolist()
        filename_lists.append(group_filenames)


    sorted_filenames_across_inputs = [sorted(lt) for lt in filename_lists]
    if all(sorted_filenames_across_inputs[0] == ele for ele in sorted_filenames_across_inputs): #comparing lists
        logging.info("Checked filenames across provided input files")
    else:
        logging.error("Filenames across provided input files are not consistent. Please provide the same filenames across all input files")
        print("Filenames across provided input files are not consistent. Please provide the same filenames across all input files")
        sys.exit(1)

    return None

def get_grouping_dict(grouping_file):

    df = pd.read_csv(grouping_file, sep="\t")

    groups = list(set(df['Group'].tolist()))

    grouping_dict = {}

    for group in groups:
        subset_df = df[df['Group'] == group]
        files = list(set(subset_df['Filename'].tolist()))
        grouping_dict[str(group)] = files

    return grouping_dict

cv = lambda x: np.std(x, ddof=1) / np.mean(x) * 100

def cv_status(row, thres):
    if row > thres:
        return "FAIL"
    else:
        return "PASS"

def check_threshold(row, threshold):
    if threshold > float(row):
        return "FAIL"
    else:
        return "PASS"

def groupname(filename, groups):

    for group in groups:
        if filename in groups[group]:
            return group

def quant_status(number, threshold):
    if number >= threshold:
        return "PASS"
    else:
        return "FAIL"

def transpose_DF(input_DF):
    input_DF_T = input_DF.T
    input_DF_T.index.name = 'newhead'
    input_DF_T.reset_index(inplace=True)
    return (input_DF_T)

def perc_qc(row, perc_thres):
    if float(row) < perc_thres:
        return "FAIL"
    else:
        return "PASS"

def label_outlier(value, outliers):
        if value in outliers:
            return 1
        else:
            return 0

def get_outlier_and_cv_status(series):
    outlier = series[0]
    threshold = series[1]

    if outlier == 1:
        return "FAIL"
    elif threshold == "FAIL":
        return "FAIL"
    else:
        return "PASS"

def only_outlier_status(outlier_value):
    if outlier_value == 1:
        return "FAIL"
    else:
        return "PASS"

def get_series_status(series):

    if series[0] == "FAIL":
        return "FAIL"
    elif series[1] == "FAIL":
        return "FAIL"
    else:
        return "PASS"

def get_overall_qc_status(series, cols_len):

    status_list = []
    for i in range(cols_len):
        status_list.append(series[i])

    if list(set(status_list)) == ['PASS']:
        return pd.Series(['PASS',f'0 out of {cols_len} metrics'])
    else:
        failed = status_list.count('FAIL')
        return pd.Series(['FAIL',f'{failed} out of {cols_len} metrics'])

#adding a function to help check the range of the input values
def int_range(min_value, max_value):
    def _type_checker(value):
        ivalue = int(value)
        if ivalue < min_value or ivalue > max_value:
            raise argparse.ArgumentTypeError("Value must be between {} and {}".format(min_value, max_value))
        return ivalue
    return _type_checker

#adding a function to help check if the file is tab delimited or not.
def is_tab_delimited(filename):
    with open(filename, 'r') as file:
        first_line = file.readline()
        #print(first_line)
        return '\t' in first_line


#adding a function to help check for duplicate values in the input file
def check_duplicates(*args):
    counts = Counter(args)
    duplicates = [item for item, count in counts.items() if count > 1]
    if duplicates == [False]:
        return False
    return duplicates
