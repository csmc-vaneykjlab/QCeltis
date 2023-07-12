"""
General Functions used for extracting ID-Based and ID-Free Metrics
"""

import pandas as pd
import numpy as np
import os
import sys
import logging
import argparse

#-------------------------------------------------------------------------- FUNCTIONS ---------------------------------------------------------------------------

def check_path(filepath):
    if not os.path.exists(str(filepath)):
        raise FileNotFoundError(f"{filepath} doesn't exist, please specify another output path")
    return None

def check_file(txtfile, level):

    #add error handling here - what if the file is not txt or whatever
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

    groups = list(set(group_df['Group'].tolist()))

    if len(groups) < 2:
        print("At least 2 groups should be provided")
        logging.error("At least 2 groups should be provided")
        sys.exit(1)

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
        grouping_dict[group] = files

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

def get_idfree_sample_qc_status(series):
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

#adding a function to help check the range of the input values
def int_range(min_value, max_value):
    def _type_checker(value):
        ivalue = int(value)
        if ivalue < min_value or ivalue > max_value:
            raise argparse.ArgumentTypeError("Value must be between {} and {}".format(min_value, max_value))
        return ivalue
    return _type_checker