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

# This function checks if a given file path exists in the file system.
# It takes one argument, 'filepath', which is expected to be a string representing the path to a file.
# If the file path does not exist, it raises a FileNotFoundError with a message including the provided path.
# This is useful for validating file paths before attempting to access files.
def check_path(filepath):
    """
    Checks if the specified file path exists.

    Args:
    filepath (str): The path to the file that needs to be verified.

    Raises:
    FileNotFoundError: If the specified file path does not exist.
    
    Returns:
    None: Returns nothing if the file path exists.
    """
    if not os.path.exists(str(filepath)):
        raise FileNotFoundError(f"{filepath} doesn't exist, please specify another output path")
    return None

def check_file(txtfile, level):
    """
    Verifies the format and content of a provided text file based on the specified analysis level.

    Args:
    txtfile (str): Path to the text file to be checked.
    level (str): The analysis level, which determines the required columns in the file.

    Raises:
    SystemExit: If the file is not in the expected format or lacks required columns, an error is logged, and the program exits.

    Note:
    This function checks for tab-delimitation and the presence of specific columns based on the 'level' argument.
    It also ensures there are enough samples for analysis.
    """
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
    """
    Verifies the format and content of the grouping file.

    Args:
    txtfile (str): Path to the main text file used for reference.
    grouping_file (str): Path to the grouping file to be checked.

    Raises:
    SystemExit: Exits the program if the grouping file is not in the expected format or lacks required columns.

    Note:
    This function checks for tab-delimitation, the presence of 'Filename' and 'Group' columns, and ensures that all samples in the grouping file are listed in the main file.
    """

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
    """
    Checks consistency of filenames across multiple input sources.

    Args:
    mzml_dir (str): Directory containing mzML files.
    protein_level (str): Path to the protein-level data file.
    peptide_level (str): Path to the peptide-level data file.
    precursor_level (str): Path to the precursor-level data file.
    grouping_file (str): Path to the grouping file.

    Raises:
    SystemExit: If filenames across input files are inconsistent, an error message is logged, and the program exits.

    Note:
    This function aggregates filenames from different sources and checks for consistency across all inputs.
    """

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
    """
    Creates a dictionary mapping groups to their corresponding filenames.

    Args:
    grouping_file (str): Path to the grouping file.

    Returns:
    dict: A dictionary where keys are group names and values are lists of filenames belonging to each group.
    """

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
    """
    Determines the pass/fail status based on the coefficient of variation.

    Args:
    row (float): The calculated CV value.
    thres (float): Threshold for the CV.

    Returns:
    str: "PASS" if CV is below the threshold, otherwise "FAIL".
    """

    if row > thres:
        return "FAIL"
    else:
        return "PASS"

def check_threshold(row, threshold):
    """
    Evaluates if a given value meets a specified threshold.

    Args:
    row (float): The value to be checked.
    threshold (float): The threshold to compare against.

    Returns:
    str: "FAIL" if the value is below the threshold, otherwise "PASS".
    """
    if threshold > float(row):
        return "FAIL"
    else:
        return "PASS"

def groupname(filename, groups):
    """
    Finds the group name for a given filename.

    Args:
    filename (str): The filename to be checked.
    groups (dict): Dictionary mapping groups to their filenames.

    Returns:
    str: The name of the group the filename belongs to.
    """
    for group in groups:
        if filename in groups[group]:
            return group

def quant_status(number, threshold):
    """
    Determines if a numeric value meets a specified quantitative threshold.

    Args:
    number (float): The numeric value to be checked.
    threshold (float): The threshold to compare against.

    Returns:
    str: "PASS" if the number is equal to or greater than the threshold, otherwise "FAIL".
    """
    if number >= threshold:
        return "PASS"
    else:
        return "FAIL"

def transpose_DF(input_DF):
    """
    Transposes a DataFrame and resets the index.

    Args:
    input_DF (DataFrame): The DataFrame to be transposed.

    Returns:
    DataFrame: The transposed DataFrame with the index reset.
    """
    input_DF_T = input_DF.T
    input_DF_T.index.name = 'newhead'
    input_DF_T.reset_index(inplace=True)
    return (input_DF_T)

def perc_qc(row, perc_thres):
    """
    Evaluates if a percentage value meets a specific threshold.

    Args:
    row (float): The percentage value to be evaluated.
    perc_thres (float): The threshold percentage.

    Returns:
    str: "FAIL" if the value is below the threshold, otherwise "PASS".
    """
    if float(row) < perc_thres:
        return "FAIL"
    else:
        return "PASS"

def label_outlier(value, outliers):
    """
    Labels a value as an outlier.

    Args:
    value: The value to be checked.
    outliers (list): A list of outlier values.

    Returns:
    int: 1 if the value is an outlier, otherwise 0.
    """
    if value in outliers:
        return 1
    else:
        return 0

def get_outlier_and_cv_status(series):
    """
    Determines the status based on outlier and CV threshold values.

    Args:
    series (iterable): Contains the outlier flag and CV status.

    Returns:
    str: "FAIL" if any condition fails, otherwise "PASS".
    """
    outlier = series[0]
    threshold = series[1]

    if outlier == 1:
        return "FAIL"
    elif threshold == "FAIL":
        return "FAIL"
    else:
        return "PASS"

def only_outlier_status(outlier_value):
    """
    Determines if a value is an outlier.

    Args:
    outlier_value (int): The outlier flag value.

    Returns:
    str: "FAIL" if the value is an outlier, otherwise "PASS".
    """
    if outlier_value == 1:
        return "FAIL"
    else:
        return "PASS"

def get_series_status(series):
    """
    Determines the overall status from a series of pass/fail conditions.

    Args:
    series (iterable): Contains multiple status values.

    Returns:
    str: "FAIL" if any status is "FAIL", otherwise "PASS".
    """
    if series[0] == "FAIL":
        return "FAIL"
    elif series[1] == "FAIL":
        return "FAIL"
    else:
        return "PASS"

def get_overall_qc_status(series, cols_len):
    """
    Computes the overall QC status based on multiple metrics.

    Args:
    series (iterable): Series of individual status values.
    cols_len (int): Total number of metrics.

    Returns:
    pd.Series: A series containing the overall status and the failure count.
    """
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
    """
    Creates a type checker function for argparse that validates if a value is within a specified range.

    Args:
    min_value (int): The minimum acceptable value.
    max_value (int): The maximum acceptable value.

    Returns:
    Function: A function that checks if a value is within the specified range.
    """
    def _type_checker(value):
        ivalue = int(value)
        if ivalue < min_value or ivalue > max_value:
            raise argparse.ArgumentTypeError("Value must be between {} and {}".format(min_value, max_value))
        return ivalue
    return _type_checker

#adding a function to help check if the file is tab delimited or not.
def is_tab_delimited(filename):
    """
    Checks if the given file is tab-delimited.

    Args:
    filename (str): The file to be checked.

    Returns:
    bool: True if the file is tab-delimited, False otherwise.
    """
    with open(filename, 'r') as file:
        first_line = file.readline()
        #print(first_line)
        return '\t' in first_line


#adding a function to help check for duplicate values in the input file
def check_duplicates(*args):
    """
    Checks for duplicate values in the provided arguments.

    Args:
    *args: Variable length argument list to check for duplicates.

    Returns:
    list or bool: List of duplicates if present, otherwise False.
    """
    counts = Counter(args)
    duplicates = [item for item, count in counts.items() if count > 1]
    if duplicates == [False]:
        return False
    return duplicates
