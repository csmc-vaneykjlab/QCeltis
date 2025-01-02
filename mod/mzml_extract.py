import pymzml
import statistics
import pandas as pd
import numpy as np
import os
import threading
import logging
import xlsxwriter
import plotly.express as px
import plotly
import plotly.graph_objects as go
import plotly.offline as offline
from jinja2 import Environment, FileSystemLoader
from scipy.stats import shapiro
import time

from mod.general_functions import cv, cv_status, check_threshold, groupname, label_outlier, get_outlier_and_cv_status, only_outlier_status, get_series_status, color_list

#-------------------------------------------------------------------------- FUNCTIONS ---------------------------------------------------------------------------

def get_mzml_list(mzml_dir):

    """
    Retrieves a list of .mzML files from the specified directory.

    Args:
    mzml_dir (str): Directory containing mzML files.

    Returns:
    list: List of paths to mzML files in the specified directory.
    """
    mzml_list = os.listdir(mzml_dir)

    final_mzml_list = []

    logging.info(f"Getting list of mzML files under provided directory path : {mzml_dir}")
    for mzml_file in mzml_list:
        if not mzml_file.endswith(".mzML"):
            logging.info(f"{mzml_file} is not a mzML file and will not be used for data extraction")
            continue
        full_filename = f"{mzml_dir}/{mzml_file}"
        final_mzml_list.append(full_filename)

    logging.info(f"{len(final_mzml_list)} mzML files have been found under {mzml_dir}")

    return final_mzml_list

def mzml_extract(mzml_path, mzml_data):

    """
    Extracts data from an mzML file and appends it to a provided list.

    Args:
    mzml_path (str): Path to the mzML file.
    mzml_data (list): List to which extracted data is appended.

    Returns:
    list: Updated list with data extracted from the mzML file.
    """

    logging.info(f"Extracting file: {mzml_path}")

    data_dict = {}
    basepeak_intensity_list = []
    ms1_spectra = 0
    ms2_spectra = 0
    ms1_tic = 0
    ms2_tic = 0

    msrun = pymzml.run.Reader(mzml_path)

    for spectrum in msrun:
        #getting basepeak intensity
        if 'base peak intensity' in spectrum: #MS:1000505 --> accession for basepeak intensity
            basepeak_intensity_value = spectrum['base peak intensity']
            basepeak_intensity_list.append(basepeak_intensity_value)
        #getting spectra count + tic
        if 'ms level' and 'total ion current' in spectrum:
            #getting spectra and tic from proteowizard converted mzML files
            if isinstance(spectrum['ms level'], list):
                if int(list(set(spectrum['ms level']))[0]) == 1:
                    ms1_spectra += 1
                    ms1_tic += spectrum['total ion current']
                if int(list(set(spectrum['ms level']))[0]) == 2:
                    ms2_spectra += 1
                    ms2_tic += spectrum['total ion current']
            else: #for thermoconverter converted mzML files
                if spectrum['ms level'] == 1:
                    ms1_spectra += 1
                    ms1_tic += spectrum['total ion current']
                if spectrum['ms level'] == 2:
                    ms2_spectra += 1
                    ms2_tic += spectrum['total ion current']

    data_dict['Filename'] = os.path.split(mzml_path)[1]

    #throwing an error if no spectra can be extracted from the mzML file
    if ms1_spectra == 0 and ms2_spectra == 0:
        logging.error(f"Not able to read spectra in the mzML file. Please remove this file: {mzml_path} from all inputs.")
        sys.exit(1)

    else:
        data_dict['MS1 Spectra'] = ms1_spectra
        data_dict['MS2 Spectra'] = ms2_spectra
        data_dict['MS2/MS1 Spectra'] = (ms2_spectra/ms1_spectra)

    if ms1_tic != 0:
        data_dict['Log MS1 TIC'] = np.log2(ms1_tic)

    if ms2_tic != 0:
        data_dict['Log MS2 TIC'] = np.log2(ms2_tic)

    if len(basepeak_intensity_list) != 0:
        data_dict['Log Max Basepeak Intensity'] = np.log2(max(basepeak_intensity_list))

    mzml_data.append(data_dict)

    return mzml_data

def get_mzml_info_dataframe(mzml_list):

    """
    Processes a list of mzML files in parallel and compiles extracted data into a DataFrame.

    Args:
    mzml_list (list): List of mzML file paths.

    Returns:
    DataFrame: DataFrame containing compiled data from all mzML files.
    """

    mzml_data = []
    threads = []

    for filename in mzml_list:

        job = threading.Thread(target=mzml_extract, args=(filename, mzml_data))
        threads.append(job)
        job.start()

    #Finish all threads
    for job in threads:
        job.join()

    mzml_dataframe = pd.DataFrame(mzml_data)
    mzml_dataframe = mzml_dataframe.sort_values("Filename")

    return mzml_dataframe

def apply_idfree_thresholds(mzml_df, mzml_threshold_dict):

    """
    Applies specified thresholds to mzML data and updates the DataFrame with QC status.

    Args:
    mzml_df (DataFrame): DataFrame containing mzML data.
    mzml_threshold_dict (dict): Dictionary containing threshold values for various metrics.

    Returns:
    DataFrame: Updated mzML DataFrame with added columns for QC status based on thresholds.
    """

    if mzml_threshold_dict['MS1 TIC Threshold'] and 'Log MS1 TIC' in mzml_df.columns.tolist():
        mzml_df[f"MS1TIC QC Threshold = {mzml_threshold_dict['MS1 TIC Threshold']}"] = mzml_df['Log MS1 TIC'].apply(check_threshold, args=[mzml_threshold_dict['MS1 TIC Threshold'],])

    if mzml_threshold_dict['MS2 TIC Threshold'] and 'Log MS2 TIC' in mzml_df.columns.tolist():
        mzml_df[f"MS2TIC QC Threshold = {mzml_threshold_dict['MS2 TIC Threshold']}"] = mzml_df['Log MS2 TIC'].apply(check_threshold, args=[mzml_threshold_dict['MS2 TIC Threshold'],])

    if mzml_threshold_dict['MS1 Spectra Threshold']:
        mzml_df[f"MS1Spectra QC Threshold = {mzml_threshold_dict['MS1 Spectra Threshold']}"] = mzml_df['MS1 Spectra'].apply(check_threshold, args=[mzml_threshold_dict['MS1 Spectra Threshold'],])

    if mzml_threshold_dict['MS2 Spectra Threshold']:
        mzml_df[f"MS2Spectra QC Threshold = {mzml_threshold_dict['MS2 Spectra Threshold']}"] = mzml_df['MS2 Spectra'].apply(check_threshold, args=[mzml_threshold_dict['MS2 Spectra Threshold'],])

    if mzml_threshold_dict['Max Basepeak Intensity Threshold'] and 'Log Max Basepeak Intensity' in mzml_df.columns.tolist():
        mzml_df[f"Max Basepeak Intensity QC Threshold = {mzml_threshold_dict['Max Basepeak Intensity Threshold']}"] = mzml_df['Log Max Basepeak Intensity'].apply(check_threshold, args=[mzml_threshold_dict['Max Basepeak Intensity Threshold'],])

    return mzml_df

def iqr_outliers(df, colname, iqr_sensitivity):

    """
    Identifies outliers in a DataFrame column based on the Interquartile Range (IQR).

    Args:
    df (DataFrame): DataFrame containing the data.
    colname (str): Name of the column to check for outliers.
    iqr_sensitivity (float): Sensitivity value for determining IQR decision range for outlier detection. Default is 1.5

    Returns:
    tuple: Updated DataFrame with outliers marked and the number of outliers found.
    """

    #calculating quantiles
    q1=df[colname].quantile(0.25)
    q3=df[colname].quantile(0.75)

    #iqr
    IQR=q3-q1

    outliers = df[((df[colname]<(q1-(iqr_sensitivity*IQR))) | (df[colname]>(q3+(iqr_sensitivity*IQR))))][colname].tolist()

    df[f"{colname} Outliers"] = df[colname].apply(label_outlier, args=[outliers, ])

    col_iqr_range = (q1-(iqr_sensitivity*IQR), q3+(iqr_sensitivity*IQR))

    return (df, len(outliers), col_iqr_range)

def outlier_detection(mzml_df, iqr_sensitivity):

    """
    Detects outliers in mzML data using either z-score or IQR methods based on data distribution.

    Args:
    mzml_df (DataFrame): DataFrame containing mzML data.
    iqr_sensitivity (float): Sensitivity value for determining IQR decision range for outlier detection. Default is 1.5

    Returns:
    DataFrame: mzML DataFrame updated with outlier detection results.
    """

    idfree_metrics = ['Log MS1 TIC', 'Log MS2 TIC', 'MS2/MS1 Spectra', 'Log Max Basepeak Intensity']
    iqr_ranges = {}

    for colname in mzml_df.columns.tolist():
        if colname in idfree_metrics:

            logging.info(f"Checking outliers for {colname}")

            mzml_df, num_outliers, col_iqr_range = iqr_outliers(mzml_df, colname, iqr_sensitivity)

            if num_outliers == 0:
                logging.info(f"No outliers found for {colname}")
            else:
                logging.info(f"{num_outliers} outliers were found for {colname}")

            iqr_ranges[colname] = col_iqr_range

    return (mzml_df, iqr_ranges)

def calculate_tic_cv(mzml_df, groups, tic_cv_threshold):

    """
    Calculates the Coefficient of Variation (CV%) for MS1 and MS2 TIC values across different groups.

    Args:
    mzml_df (DataFrame): DataFrame containing mzML data.
    groups (dict): Dictionary mapping groups to filenames.
    tic_cv_threshold (float): Threshold for the TIC CV%.

    Returns:
    DataFrame: DataFrame containing CV% for MS1 and MS2 TIC values.
    """

    tic_cv = mzml_df[['Filename','Log MS1 TIC','Log MS2 TIC']]

    group_cv = {}

    for group in groups:
        group_subset = tic_cv[tic_cv['Filename'].isin(groups[group])]
        ms1tic_cv = round(cv(group_subset['Log MS1 TIC'].tolist()),2)
        ms2tic_cv = round(cv(group_subset['Log MS2 TIC'].tolist()),2)
        group_cv[group] = {"Log MS1 TIC CV%": ms1tic_cv, "Log MS2 TIC CV%": ms2tic_cv}

    tic_cv = pd.DataFrame.from_dict(group_cv, orient="index")
    tic_cv.index = tic_cv.index.set_names(['Group'])
    tic_cv.reset_index(drop=False, inplace=True)
    tic_cv = tic_cv.sort_values('Group')

    tic_cv[f'MS1 TIC CV% Threshold = {int(tic_cv_threshold)}'] = tic_cv['Log MS1 TIC CV%'].apply(cv_status, args=[tic_cv_threshold,])
    tic_cv[f'MS2 TIC CV% Threshold = {int(tic_cv_threshold)}'] = tic_cv['Log MS2 TIC CV%'].apply(cv_status, args=[tic_cv_threshold,])

    return tic_cv

def get_sample_qc(mzml_df, mzml_threshold_dict, groupwise_comparison, groups):

    """
    Applies QC status checks to mzML data based on thresholds and outlier detection.

    Args:
    mzml_df (DataFrame): DataFrame containing mzML data.
    mzml_threshold_dict (dict): Dictionary containing QC thresholds.
    groupwise_comparison (bool): Indicates if group-wise comparison is being used.
    groups (dict): Dictionary mapping groups for analysis.

    Returns:
    DataFrame: Updated mzML DataFrame with sample QC status.
    """

    if 'Log MS1 TIC' in mzml_df.columns.tolist():
        if mzml_threshold_dict['MS1 TIC Threshold']:
            mzml_df['MS1 TIC Sample QC Status'] = mzml_df[['Log MS1 TIC Outliers',f"MS1TIC QC Threshold = {mzml_threshold_dict['MS1 TIC Threshold']}"]].apply(get_outlier_and_cv_status, axis=1)
        else:
            mzml_df['MS1 TIC Sample QC Status'] =  mzml_df['Log MS1 TIC Outliers'].apply(only_outlier_status)

    if 'Log MS2 TIC' in mzml_df.columns.tolist():
        if mzml_threshold_dict['MS2 TIC Threshold']:
            mzml_df['MS2 TIC Sample QC Status'] = mzml_df[['Log MS2 TIC Outliers',f"MS2TIC QC Threshold = {mzml_threshold_dict['MS2 TIC Threshold']}"]].apply(get_outlier_and_cv_status, axis=1)
        else:
            mzml_df['MS2 TIC Sample QC Status'] = mzml_df['Log MS2 TIC Outliers'].apply(only_outlier_status)

    if mzml_threshold_dict['MS1 Spectra Threshold']:
        mzml_df['MS1 Spectra QC Status'] = mzml_df[['MS2/MS1 Spectra Outliers',f"MS1Spectra QC Threshold = {mzml_threshold_dict['MS1 Spectra Threshold']}"]].apply(get_outlier_and_cv_status, axis=1)
    else:
        mzml_df['MS1 Spectra QC Status'] = mzml_df['MS2/MS1 Spectra Outliers'].apply(only_outlier_status)

    if mzml_threshold_dict['MS2 Spectra Threshold']:
        mzml_df['MS2 Spectra QC Status'] = mzml_df[['MS2/MS1 Spectra Outliers',f"MS2Spectra QC Threshold = {mzml_threshold_dict['MS2 Spectra Threshold']}"]].apply(get_outlier_and_cv_status, axis=1)
    else:
        mzml_df['MS2 Spectra QC Status'] = mzml_df['MS2/MS1 Spectra Outliers'].apply(only_outlier_status)

    if 'Max Basepeak Intensity' in mzml_df.columns.tolist():
        if mzml_threshold_dict['Max Basepeak Intensity Threshold']:
            mzml_df['Max Basepeak Intensity QC Status'] = mzml_df[['Log Max Basepeak Intensity Outliers', f"Max Basepeak Intensity QC Threshold = {mzml_threshold_dict['Max Basepeak Intensity Threshold']}"]].apply(get_outlier_and_cv_status, axis=1)
        else:
            mzml_df['Max Basepeak Intensity QC Status'] = mzml_df['Log Max Basepeak Intensity Outliers'].apply(only_outlier_status)

    sample_qc_cols = ['MS1 TIC Sample QC Status', 'MS2 TIC Sample QC Status', 'MS1 Spectra QC Status', 'MS2 Spectra QC Status', 'Max Basepeak Intensity QC Status']
    matched_sample_qc_cols = []
    for colname in mzml_df.columns.tolist():
        if colname in sample_qc_cols:
            matched_sample_qc_cols.append(colname)

    mzml_df = mzml_df[['Filename'] + matched_sample_qc_cols]
    mzml_df = mzml_df.sort_values('Filename')

    return mzml_df

def get_idfree_grouped_df(mzml_sample_df, tic_cv, tic_cv_threshold, groups):
    """
    Compiles a grouped DataFrame for ID-free data, merging TIC CV% and sample QC status.

    Args:
    mzml_sample_df (DataFrame): DataFrame containing sample QC status.
    tic_cv (DataFrame): DataFrame containing TIC CV% data.
    tic_cv_threshold (float): Threshold for TIC CV%.
    groups (dict): Dictionary mapping groups to filenames.

    Returns:
    DataFrame: Compiled DataFrame with grouped QC status for ID-free data.
    """

    tic_group_df = tic_cv[['Group',f'MS1 TIC CV% Threshold = {int(tic_cv_threshold)}', f'MS2 TIC CV% Threshold = {int(tic_cv_threshold)}']]
    idfree_status_params = ['MS1 TIC Sample QC Status', 'MS2 TIC Sample QC Status', 'MS1 Spectra QC Status', 'MS2 Spectra QC Status', 'Max Basepeak Intensity QC Status']
    mzml_sample_df['Group'] = mzml_sample_df['Filename'].apply(groupname, args=[groups, ])

    group_status_dict = {}

    for group in list(set(mzml_sample_df['Group'].tolist())):
        group_subset = mzml_sample_df[mzml_sample_df['Group'] == group]
        col_dict = {}

        for colname in ['MS1 TIC Sample QC Status', 'MS2 TIC Sample QC Status', 'MS1 Spectra QC Status', 'MS2 Spectra QC Status', 'Max Basepeak Intensity QC Status']:
            if colname in group_subset.columns.tolist():
                group_colname = colname.replace("Sample", "Group")
                if list(set(group_subset[colname].tolist())) == ["PASS"]:
                    col_dict[group_colname] = "PASS"
                else:
                    col_dict[group_colname] = "FAIL"

        group_status_dict[group] = col_dict

    grouped_df = pd.DataFrame.from_dict(group_status_dict, orient="columns")
    grouped_df = grouped_df.T
    grouped_df.reset_index(drop=False, inplace=True)
    grouped_df.rename(columns={'index':'Group'}, inplace=True)

    grouped_df = pd.merge(grouped_df, tic_group_df, on='Group')

    if 'MS1 TIC Group QC Status' in grouped_df.columns.tolist():
        grouped_df['MS1 TIC Group QC Status'] = grouped_df[['MS1 TIC Group QC Status', f'MS1 TIC CV% Threshold = {int(tic_cv_threshold)}']].apply(get_series_status, axis=1)
        grouped_df = grouped_df.drop(f'MS1 TIC CV% Threshold = {int(tic_cv_threshold)}', axis=1)

    if 'MS2 TIC Group QC Status' in grouped_df.columns.tolist():
        grouped_df['MS2 TIC Group QC Status'] = grouped_df[['MS2 TIC Group QC Status', f'MS2 TIC CV% Threshold = {int(tic_cv_threshold)}']].apply(get_series_status, axis=1)
        grouped_df = grouped_df.drop(f'MS2 TIC CV% Threshold = {int(tic_cv_threshold)}', axis=1)

    grouped_df = grouped_df.sort_values('Group')

    return grouped_df

#------------------------------------------------------------------------ PLOT FUNCTIONS ----------------------------------------------------------------------------

def tic_plots(mzml_df, tic_cv, ms1_tic_threshold, ms2_tic_threshold, tic_cv_threshold, groupwise_comparison, color_list, iqr_ranges):

    df = mzml_df[['Filename','Log MS1 TIC','Log MS2 TIC']]

    df = df.melt(id_vars=["Filename"],
        var_name="Label",
        value_name="TIC")

    tic_line = px.line(df, x='Filename', y="TIC", title="Total Ion Current", color="Label", line_shape="spline", markers=True)
    tic_line.update_xaxes(tickfont_size=6)
    tic_line.update_layout(title={'font': {'size': 9}})
    tic_line.update_layout(
            margin=dict(l=20, r=20, t=20, b=20)
    )

    if ms1_tic_threshold:
        tic_line.add_hline(y=ms1_tic_threshold, line_dash="dot", annotation_text=f"MS1 TIC Threshold = {ms1_tic_threshold}")

    if ms2_tic_threshold:
        tic_line.add_hline(y=ms2_tic_threshold, line_dash="dot", annotation_text=f"MS2 TIC Threshold = {ms2_tic_threshold}")

    tic_plot = plotly.io.to_html(tic_line, include_plotlyjs=False, full_html=False, default_width='900px', default_height='450px')
    #tic_plot = offline.plot(tic_line, output_type='div', include_plotlyjs=False)

    tic_report_params = {'total_ion_current': True,
                    'tic_plot': tic_plot,
                    'tic_ms_plot_description': 'The total ion current (TIC) is the summed intensity across the entire range of masses being detected in each sample. MS1 and MS2 Total Ion Current Values extracted from spectra within the given mzML files'}

    if not list(set(mzml_df['Log MS1 TIC Outliers'].tolist())) == [0]:
        ms1_outliers = mzml_df[mzml_df['Log MS1 TIC Outliers'] == 1]['Filename'].tolist()
        ms1_outliers_filenames = ", ".join(ms1_outliers)
        tic_report_params['tic_ms1_outlier_description'] = f"{mzml_df['Log MS1 TIC Outliers'].tolist().count(1)} outliers were found. The following files have been detected as outliers: {ms1_outliers_filenames}"

        tic_ms1_outlier = px.scatter(mzml_df, x='Filename', y='Log MS1 TIC', title="Log MS1 TIC Outliers", color='Log MS1 TIC Outliers')
        tic_ms1_outlier.update_layout(title={'font': {'size': 9}})
        tic_ms1_outlier.add_hline(y=iqr_ranges['Log MS1 TIC'][1], line_width=1.5, line_dash='dash', line_color="red")
        tic_ms1_outlier.add_hline(y=iqr_ranges['Log MS1 TIC'][0], line_width=1.5, line_dash='dash', line_color="red")
        tic_ms1_outlier.update_xaxes(tickfont_size=6)
        tic_ms1_outlier.update_layout(
                margin=dict(l=20, r=20, t=20, b=20)
        )
        tic_ms1_outlier.update_traces(marker=dict(line=dict(color='black', width=1)))

        tic_ms1_outlier_plot = plotly.io.to_html(tic_ms1_outlier, include_plotlyjs=False, full_html=False, default_width='900px', default_height='450px')

        tic_report_params['tic_ms1_outlier_plot'] = tic_ms1_outlier_plot

    if not list(set(mzml_df['Log MS2 TIC Outliers'].tolist())) == [0]:
        ms2_outliers = mzml_df[mzml_df['Log MS2 TIC Outliers'] == 1]['Filename'].tolist()
        ms2_outliers_filenames = ", ".join(ms2_outliers)
        tic_report_params['tic_ms2_outlier_description'] = f"{mzml_df['Log MS2 TIC Outliers'].tolist().count(1)} outliers were found. The following files have been detected as outliers: {ms2_outliers_filenames}"

        tic_ms2_outlier = px.scatter(mzml_df, x='Filename', y='Log MS2 TIC', title="Log MS2 TIC Outliers", color='Log MS2 TIC Outliers')
        tic_ms2_outlier.update_layout(title={'font': {'size': 9}})
        tic_ms2_outlier.add_hline(y=iqr_ranges['Log MS2 TIC'][1], line_width=1.5, line_dash='dash', line_color="red")
        tic_ms2_outlier.add_hline(y=iqr_ranges['Log MS2 TIC'][0], line_width=1.5, line_dash='dash', line_color="red")
        tic_ms2_outlier.update_xaxes(tickfont_size=6)
        tic_ms2_outlier.update_layout(
                margin=dict(l=20, r=20, t=20, b=20)
        )
        tic_ms2_outlier.update_traces(marker=dict(line=dict(color='black', width=1)))

        tic_ms2_outlier_plot = plotly.io.to_html(tic_ms2_outlier, include_plotlyjs=False, full_html=False, default_width='900px', default_height='450px')

        tic_report_params['tic_ms2_outlier_plot'] = tic_ms2_outlier_plot

    if groupwise_comparison:
        tic_report_params['tic_ms_cv_description'] = "When a grouping file is provided, CV% for TIC values across samples in each group is calculated. This provides an insignt into how consistent the samples are within each group."

        ms1tic_bar = px.bar(tic_cv, x='Group', y="Log MS1 TIC CV%", title="MS1 Total Ion Current - CV%", color="Group", color_discrete_sequence=color_list)
        ms1tic_bar.update_xaxes(tickfont_size=6)
        ms1tic_bar.update_layout(title={'font': {'size': 9}})
        ms1tic_bar.add_hline(y=tic_cv_threshold, line_dash="dot", annotation_text=f"TIC CV Threshold = {tic_cv_threshold}")
        ms1tic_bar.update_layout(
            margin=dict(l=20, r=20, t=20, b=20),
        )

        ms1_tic = plotly.io.to_html(ms1tic_bar, include_plotlyjs=False, full_html=False, default_width='900px', default_height='450px')
        tic_report_params['tic_ms1_cv_plot'] = ms1_tic

        if list(set(tic_cv[f'MS1 TIC CV% Threshold = {int(tic_cv_threshold)}'].tolist())) == ['PASS']:
            tic_report_params['tic_ms1_cv_description'] = 'CV% for MS1 TIC was calculated using TIC values from each sample within a provided group. All groups have passed the CV Threshold'
        else:
            failed_ms1_groups = ", ".join(tic_cv[tic_cv[f'MS1 TIC CV% Threshold = {int(tic_cv_threshold)}'] == 'FAIL']['Group'].tolist())
            tic_report_params['tic_ms1_cv_description'] = f'CV% for MS1 TIC was calculated using TIC values from each sample within a provided group. The following groups have not met the CV Threshold: {failed_ms1_groups}. This represents an inconsistent TIC pattern, please check the samples within the failed groups.'

        ms2tic_bar = px.bar(tic_cv, x='Group', y="Log MS2 TIC CV%", title="MS2 Total Ion Current - CV%", color="Group", color_discrete_sequence=color_list)
        ms2tic_bar.update_xaxes(tickfont_size=6)
        ms2tic_bar.update_layout(title={'font': {'size': 9}})
        ms2tic_bar.add_hline(y=tic_cv_threshold, line_dash="dot", annotation_text=f"TIC CV Threshold = {tic_cv_threshold}")
        ms2tic_bar.update_layout(
            margin=dict(l=20, r=20, t=20, b=20),
        )

        ms2_tic = plotly.io.to_html(ms2tic_bar, include_plotlyjs=False, full_html=False, default_width='900px', default_height='450px')
        tic_report_params['tic_ms2_cv_plot'] = ms2_tic

        if list(set(tic_cv[f'MS2 TIC CV% Threshold = {int(tic_cv_threshold)}'].tolist())) == ['PASS']:
            tic_report_params['tic_ms2_cv_description'] = 'CV% for MS2 TIC was calculated using TIC values from each sample within a provided group. All groups have passed the CV Threshold'
        else:
            failed_ms2_groups = ", ".join(tic_cv[tic_cv[f'MS2 TIC CV% Threshold = {int(tic_cv_threshold)}'] == 'FAIL']['Group'].tolist())
            tic_report_params['tic_ms2_cv_description'] = f'CV% for MS1 TIC was calculated using TIC values from each sample within a provided group. The following groups have not met the CV Threshold: {failed_ms2_groups}. This represents an inconsistent TIC pattern, please check the samples within the failed groups.'

    return tic_report_params

def spectral_plot(mzml_df, iqr_ranges):

    df = mzml_df[['Filename','MS2/MS1 Spectra', 'MS2/MS1 Spectra Outliers']]

    count_line = px.line(df, x='Filename', y="MS2/MS1 Spectra", title="MS2/MS1 Spectral Ratio", line_shape="spline", markers=True)
    count_line.update_xaxes(tickfont_size=6)
    count_line.update_layout(title={'font': {'size': 9}})
    count_line.update_layout(
        margin=dict(l=20, r=20, t=20, b=20),
    )

    spectral_count = plotly.io.to_html(count_line, include_plotlyjs=False, full_html=False, default_width='900px', default_height='450px')
    spectra_report_params = {'ms2_ms1_spectral_ratio': True,
                            'ms2_ms1_spectral_ratio_plot': spectral_count,
                            'ms2_ms1_spectral_ratio_description': 'MS2/MS1 Spectral Count Ratio extracted from given mzML files'}

    if not list(set(mzml_df['MS2/MS1 Spectra Outliers'].tolist())) == [0]:
        spectra_outliers = mzml_df[mzml_df['MS2/MS1 Spectra Outliers'] == 1]['Filename'].tolist()
        spectra_outliers_filenames = ", ".join(spectra_outliers)
        spectra_report_params['ms2_ms1_spectral_ratio_outlier_description'] = f"{mzml_df['MS2/MS1 Spectra Outliers'].tolist().count(1)} outliers were found. The following files have been detected as outliers: {spectra_outliers_filenames}"

        ms2_ms1_spectral_ratio_outlier = px.scatter(mzml_df, x='Filename', y='MS2/MS1 Spectra', title="MS2/MS1 Spectra Outliers", color='MS2/MS1 Spectra Outliers')
        ms2_ms1_spectral_ratio_outlier.update_xaxes(tickfont_size=6)
        ms2_ms1_spectral_ratio_outlier.add_hline(y=iqr_ranges['MS2/MS1 Spectra'][1], line_width=1.5, line_dash='dash', line_color="red")
        ms2_ms1_spectral_ratio_outlier.add_hline(y=iqr_ranges['MS2/MS1 Spectra'][0], line_width=1.5, line_dash='dash', line_color="red")
        ms2_ms1_spectral_ratio_outlier.update_layout(title={'font': {'size': 9}})
        ms2_ms1_spectral_ratio_outlier.update_layout(
                margin=dict(l=20, r=20, t=20, b=20)
        )
        ms2_ms1_spectral_ratio_outlier.update_traces(marker=dict(line=dict(color='black', width=1)))

        ms2_ms1_spectral_ratio_plot = plotly.io.to_html(ms2_ms1_spectral_ratio_outlier, include_plotlyjs=False, full_html=False, default_width='900px', default_height='450px')

        spectra_report_params['ms2_ms1_spectral_ratio_outlier_plot'] = ms2_ms1_spectral_ratio_plot

    return spectra_report_params

def basepeak_graph(mzml_df, max_basepeak_intensity_threshold, groups, groupwise_comparison, color_list, iqr_ranges):

    if groupwise_comparison:
        mzml_df['Group'] = mzml_df['Filename'].apply(groupname, args=[groups,])
        mzml_df = mzml_df.sort_values('Group')
        bp_bar = px.bar(mzml_df, x='Filename', y="Log Max Basepeak Intensity", title="Log Max Base Peak Intensity", color="Group", color_discrete_sequence=color_list)
    else:
        bp_bar = px.bar(mzml_df, x='Filename', y="Log Max Basepeak Intensity", title="Log Max Base Peak Intensity")

    bp_bar.update_layout(title={'font': {'size': 9}})
    bp_bar.update_xaxes(tickfont_size=6)
    bp_bar.update_layout(
        margin=dict(l=20, r=20, t=20, b=20),
    )

    if max_basepeak_intensity_threshold:
        bp_bar.add_hline(y=max_basepeak_intensity_threshold, line_dash="dot", annotation_text=f"Max Base Peak Intensity Threshold = {max_basepeak_intensity_threshold}")

    bp_plot = plotly.io.to_html(bp_bar, include_plotlyjs=False, full_html=False, default_width='900px', default_height='450px')

    basepeak_report_params = {'max_basepeak_intensity' : True,
                              'max_basepeak_intensity_plot': bp_plot}

    if not list(set(mzml_df['Log Max Basepeak Intensity Outliers'].tolist())) == [0]:
        bp_outliers = mzml_df[mzml_df['Log Max Basepeak Intensity Outliers'] == 1]['Filename'].tolist()
        bp_outliers_filenames = ", ".join(bp_outliers)
        basepeak_report_params['max_basepeak_intensity_outlier_description'] = f"{mzml_df['Log Max Basepeak Intensity Outliers'].tolist().count(1)} outliers were found. The following files have been detected as outliers: {bp_outliers_filenames}"

        max_basepeak_intensity_outlier = px.scatter(mzml_df, x='Filename', y='Log Max Basepeak Intensity', title="Log Max Base Peak Intensity Outliers", color='Log Max Basepeak Intensity Outliers')
        max_basepeak_intensity_outlier.update_layout(title={'font': {'size': 9}})
        max_basepeak_intensity_outlier.add_hline(y=iqr_ranges['Log Max Basepeak Intensity'][1], line_width=1.5, line_dash='dash', line_color="red")
        max_basepeak_intensity_outlier.add_hline(y=iqr_ranges['Log Max Basepeak Intensity'][0], line_width=1.5, line_dash='dash', line_color="red")
        max_basepeak_intensity_outlier.update_xaxes(tickfont_size=6)
        max_basepeak_intensity_outlier.update_layout(
                margin=dict(l=20, r=20, t=20, b=20)
        )
        max_basepeak_intensity_outlier.update_traces(marker=dict(line=dict(color='black', width=1)))

        max_basepeak_intensity_outlier_plot = plotly.io.to_html(max_basepeak_intensity_outlier, include_plotlyjs=False, full_html=False, default_width='900px', default_height='450px')

        basepeak_report_params['max_basepeak_intensity_outlier_plot'] = max_basepeak_intensity_outlier_plot

    return basepeak_report_params

def create_graphs(mzml_df, tic_cv, groupwise_comparison, groups, mzml_threshold_dict, iqr_ranges):

    if 'Log MS1 TIC' and 'Log MS2 TIC' in mzml_df.columns.tolist():
        tic_report_params = tic_plots(mzml_df, tic_cv, mzml_threshold_dict['MS1 TIC Threshold'], mzml_threshold_dict['MS2 TIC Threshold'], mzml_threshold_dict['TIC CV Threshold'], groupwise_comparison, color_list, iqr_ranges)
    else:
        logging.info("No TIC information was extracted from provided mzML files, no plots for TIC will be generated")
        tic_report_params = {}

    if 'Log Max Basepeak Intensity' in mzml_df.columns.tolist():
        basepeak_report_params = basepeak_graph(mzml_df, mzml_threshold_dict['Max Basepeak Intensity Threshold'], groups, groupwise_comparison, color_list, iqr_ranges)
    else:
        logging.info("No Basepeak Intensity information was extracted from provided mzML files, no plots for Max Basepeak Intensity will be generated")
        basepeak_report_params = {}

    spectra_report_params = spectral_plot(mzml_df, iqr_ranges)

    idfree_report_parameters = dict(tuple(tic_report_params.items()) + tuple(spectra_report_params.items()) + tuple(basepeak_report_params.items()))

    return idfree_report_parameters

#---------------------------------------------------------------------- MAIN FUNCTION CALL -------------------------------------------------------------------------

def calculate_idfree_metrics(out_dir, reportname, mzml_dir, groupwise_comparison, groups, mzml_threshold_dict):

    #getting list of mzML files
    mzml_list = get_mzml_list(mzml_dir)

    if len(mzml_list) > 30:
        mzml_list_chunks = [mzml_list[x:x+30] for x in range(0, len(mzml_list), 30)]
    else:
        mzml_list_chunks = [mzml_list]

    mzml_df = pd.DataFrame()

    for lt in mzml_list_chunks:
        #extracting data from mzml files
        extracted_df = get_mzml_info_dataframe(lt)
        mzml_df = pd.concat([mzml_df, extracted_df], ignore_index=True)
        time.sleep(100)

    #applying thresholds + outlier detection
    mzml_df = apply_idfree_thresholds(mzml_df, mzml_threshold_dict)
    mzml_df, iqr_ranges = outlier_detection(mzml_df, mzml_threshold_dict['IQR Sensitivity'])

    if groupwise_comparison:
        mzml_df['Group'] = mzml_df['Filename'].apply(groupname, args=[groups, ])
        mzml_df = mzml_df.sort_values("Group")
        tic_cv = calculate_tic_cv(mzml_df, groups, mzml_threshold_dict['TIC CV Threshold'])
    else:
        tic_cv = ""

    logging.info(f"Saving ID-Free QC Report to {out_dir}/{reportname}_ID-Free_QC_Report.xlsx")

    #saving dataframes to excel document
    writer = pd.ExcelWriter(f"{out_dir}/{reportname}_ID-Free_QC_Report.xlsx", engine='xlsxwriter')
    mzml_df.to_excel(writer, index=False, sheet_name="ID-Free Metrics Summary")
    if groupwise_comparison:
        tic_cv.to_excel(writer, index=False, sheet_name='Group TIC CV')
    writer.close()

    idfree_report_parameters = create_graphs(mzml_df, tic_cv, groupwise_comparison, groups, mzml_threshold_dict, iqr_ranges)

    mzml_sample_df = get_sample_qc(mzml_df, mzml_threshold_dict, groupwise_comparison, groups)
    if groupwise_comparison:
        idfree_grouped_df = get_idfree_grouped_df(mzml_sample_df, tic_cv, mzml_threshold_dict['TIC CV Threshold'], groups)
    else:
        idfree_grouped_df = ""

    return (mzml_sample_df, idfree_grouped_df, idfree_report_parameters)
