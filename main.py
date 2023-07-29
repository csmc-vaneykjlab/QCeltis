import argparse
from argparse import RawTextHelpFormatter
import pandas as pd
import os
import sys
import logging
import time
from jinja2 import Environment, FileSystemLoader
from mod.mzml_extract import calculate_idfree_metrics
from mod.idbased_metrics import calculate_idbased_metrics
from mod.general_functions import check_path, check_file, check_grouping_file, get_grouping_dict, check_samples, int_range, check_duplicates

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

start_time = time.time()

def main():

    logging.basicConfig(filename="qc_script.log", level=logging.INFO)

    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)

    #output parameters
    parser.add_argument('-o', '--outdirectory', type=str, required=True, help='[Required] Output Directory Path')
    parser.add_argument('-r', '--reportname', type=str, required=True, help='[Required] Report Name for HTML and Excel Reports')

    #id-free metrics - mzml extraction
    #ADD - checks for thresholds (like range check, type check etc)
    parser.add_argument('-m', '--mzml_directory', type=str, default=False, help='[Optional] Path to directory where mzML files are present')
    parser.add_argument('-t1', '--ms1_tic_threshold', type=float, default=False, help='[Optional] MS1 TIC Threshold')
    parser.add_argument('-t2', '--ms2_tic_threshold', type=float, default=False, help='[Optional] MS2 TIC Threshold')
    parser.add_argument('-s1', '--ms1_spectra_threshold', type=float, default=False, help='[Optional] MS2 Spectra Threshold')
    parser.add_argument('-s2', '--ms2_spectra_threshold', type=float, default=False, help='[Optional] MS2 Spectra Threshold')
    parser.add_argument('-bp', '--max_basepeak_intensity', type=float, default=False, help='[Optional] Maximum Basepeak Intensity Threshold')

    #id-based inputs
    parser.add_argument('-pt', '--protein_level', type=str, default=False, help='[Optional] Path to Protein Intensity File')
    parser.add_argument('-pep', '--peptide_level', type=str, default=False, help='[Optional] Path to Peptide Intensity File')
    parser.add_argument('-pre', '--precursor_level', type=str, default=False, help='[Optional] Path to Precursor Intensity File')
    parser.add_argument('-peprt', '--peptide_rt',  type=str, default=False, help='[Optional] Path to Peptide Retention Time File')
    parser.add_argument('-prert', '--precursor_rt',type=str,  default=False, help='[Optional] Path to Precursor Retention Time File')
    parser.add_argument('-g', '--grouping_file', type=str, default=False, help='[Optional] Path to Grouping File')
    parser.add_argument('-peplt', '--peptide_list', type=str, default=False, help='[Optional] Path to file containing list of peptides to monitor intensity and RT distribution across samples')

    #id-based metric thresholds
    parser.add_argument('-x', '--protein_threshold', type=int, default=False, help='[Optional] Protein Threshold for each sample')
    parser.add_argument('-y', '--peptide_threshold', type=int, default=False, help='[Optional] Peptide Threshold for each sample')
    parser.add_argument('-z', '--precursor_threshold', type=int, default=False, help='[Optional] Precursor Threshold for each sample')
    parser.add_argument('-e', '--enzyme', type=str, default=False, help="[Optional] User input enzyme. Available input of enzymes will be:\n"+
        'Asp-N\t\tcleaves at N terminus of D. \tNo exceptions\n'+
        'Lys-C\t\tcleaves at C terminus of K. \tException: KP\n'+
        'Lys-N\t\tcleaves at N terminus of K. \tNo exception\n'+
        'Thermolysin\tcleaves at N terminus of A,F,I,L,M,V. \tExceptions: DA,DF,DI,DL,DM,DV,EA,EF,EI,EL,EM,EV\n'+
        'ProteinaseK\tcleaves at C terminus A,F,Y,W,L,I,V. \tNo exception\n'+
        'MFH\t\t(Microwave-assisted formic acid hydrolysis), cleaves at C terminus of D. \tNo exception\n'+
        'Gluc\t\tcleaves at C terminus of D,E. \tExceptions: DP,EP,EE,DE\n'+
        'GlucBicarb\tcleaves at C terminus of E. \tExceptions: EP,EE\n'+
        'ChymotrypsinFLMWY\tcleaves at C terminus of F,L,M,W,Y. \tExceptions: FP,LP,MP,WP,YP,PY\n'+
        'Trypsin\t\tcleaves at C terminus of K,R. \tExceptions: KP,RP\n'+
        'TrypsinP\tcleaves at C terminus of K,R. \tNo exception\n'+
        'CNBr\t\tcleaves at C terminus of M. \tNo exceptions\n'+
        'Arg-C\t\tcleaves at C terminus of R. \tException: RP\n')
    parser.add_argument('-c', '--miscleavage_threshold', type=int_range(0, 100), default=False, help='[Optional] 0% Missed Cleavage Threshold for each sample')
    parser.add_argument('-t', '--tic_cv_threshold', type=int_range(0, 100), default=False, help='[Optional] TIC CV Threshold for groupwise QC status - Percentage between 0 and 100 - will only be used if grouping file is provided')
    parser.add_argument('-s', '--cv_percent_threshold', type=int_range(0, 100), default=False, help='[Optional] Intensity CV Threshold - Percentage between 0 and 100')
    parser.add_argument('-d', '--data_percent_threshold', type=int_range(0, 100), default=False, help='[Optional] Data Threshold for Intensity CV - will only be used if platewise comparison is selected')
    parser.add_argument('-irt', '--irtlabel', default=False, type=str, help='[Optional] If iRT peptides are present in your peptide intensity file, please provide how the iRT proteins are labelled in your dataset')
    parser.add_argument('-v', '--coverage_threshold', type=int_range(0, 100), default=False, help='[Optional] Intensity or Retention Time Coverage % Threshold in each sample')

    args = parser.parse_args()

    # -------------------------------------------- ASSIGNING ARGUMENTS -------------------------------------------

    out_dir = str(args.outdirectory)
    reportname = str(args.reportname)

    mzml_dir = args.mzml_directory
    ms1_tic_threshold = float(args.ms1_tic_threshold)
    ms2_tic_threshold = float(args.ms2_tic_threshold)
    ms1_spectra_threshold = int(args.ms1_spectra_threshold)
    ms2_spectra_threshold = int(args.ms2_spectra_threshold)
    max_basepeak_intensity_threshold = float(args.max_basepeak_intensity)

    protein_level = args.protein_level
    peptide_level = args.peptide_level
    precursor_level = args.precursor_level
    peptide_rt = args.peptide_rt
    precursor_rt = args.precursor_rt
    grouping_file = args.grouping_file
    peptide_list = args.peptide_list

    protein_threshold = int(args.protein_threshold)
    peptide_threshold = int(args.peptide_threshold)
    precursor_threshold = int(args.precursor_threshold)
    digestion_enzyme = args.enzyme
    miscleavage_threshold = float(args.miscleavage_threshold)
    tic_cv_threshold = float(args.tic_cv_threshold)
    cv_percent_threshold = float(args.cv_percent_threshold)
    data_percent_threshold = float(args.data_percent_threshold)
    irtlabel = args.irtlabel
    coverage_threshold = float(args.coverage_threshold)

    # ------------------------------------------- CHECKING INPUTS AND THRESHOLDS -----------------------------------
    logging.info("----------------------------- CHECKING PROVIDED INPUTS AND THRESHOLDS -----------------------------\n")

    logging.info("----------------------------- Checking provided inputs for duplicates -----------------------------\n")


    duplicates = check_duplicates(protein_level, peptide_level, precursor_level, grouping_file)
    if duplicates:
        print(f"ERROR: Duplicate filenames have been given {','.join(duplicates)} ")
        logging.error(f"ERROR: Duplicate filenames have been given {','.join(duplicates)}")
        sys.exit(1)


    check_path(out_dir)
    logging.info(f"All outputs will be saved to {out_dir}")

    if mzml_dir:
        logging.info("--------------------------------------- Checking mzML Directory ------------------------------------------\n")

        #checking if directory exists
        if not os.path.exists(mzml_dir):
            print(f"Given MzML Directory path doesn't exist = {mzml_dir}")
            logging.error(f"Given MzML Directory path doesn't exist = {mzml_dir}")
            sys.exit(1)
        else:
            logging.info(f"mzML files found under {mzml_dir} will be used to extract ID-Free QC Metrics")

        if grouping_file:
            if not tic_cv_threshold:
                logging.error("TIC CV threshold is not provided, no MS1 and MS2 TIC comparison between groups can be performed. Please provide TIC CV threshold")
                sys.exit(1)

    else:
        logging.info("mzML directory not provided, no id-free metrics will be calculated")


    if protein_level:
        logging.info("---------------------------------------- Checking Protein Intensity File ----------------------------------------\n")

        check_file(protein_level, "Protein")

        if grouping_file:
            groups = check_grouping_file(protein_level, grouping_file)
            logging.info(f"{grouping_file} will be used for comparing samples across the following groups: {groups}")

        if protein_threshold:
            logging.info(f"Protein Threshold: {protein_threshold} will be applied")
        else:
            logging.info("No Protein Threshold will be applied")

        if cv_percent_threshold:
            logging.info(f"Intensity CV Threshold: {cv_percent_threshold} will be applied for protein intensities")
        else:
            logging.info("No Intensity CVs will be calculated, no threshold provided")

        if data_percent_threshold:
            if not cv_percent_threshold:
                logging.info(f"Intensity CV and Data Percent Threshold both need to be provided")
            else:
                logging.info(f"Intensity CV Threshold: {cv_percent_threshold} will be applied along with data percent threshold: {data_percent_threshold}")
        else:
            logging.info("Protein Intensities over data percentage will not be calculated, no data percent threshold is provided")

    else:
        logging.info("Protein Level file not provided, protein level metrics will not be calculated \n")


    if peptide_level:

        logging.info("--------------------------------------------- Checking Peptide Intensity File --------------------------------------------- ")
        check_file(peptide_level, "Peptide")

        if grouping_file:
            groups = check_grouping_file(peptide_level, grouping_file)
            logging.info(f"{grouping_file} will be used for comparing samples across the following groups: {groups}")

            if not tic_cv_threshold:
                logging.error("TIC CV threshold is not provided, no peptide TIC comparison between groups can be performed. Please provide TIC CV threshold")
                sys.exit(1)

        if protein_threshold:
            logging.info(f"Peptide Threshold: {peptide_threshold} will be applied")
        else:
            logging.info("No Peptide Threshold will be applied")

        if cv_percent_threshold:
            logging.info(f"Intensity CV Threshold: {cv_percent_threshold} will be applied for peptide intensities")
        else:
            logging.info("No Intensity CVs will be calculated, no threshold provided")

        if data_percent_threshold:
            if not cv_percent_threshold:
                logging.info(f"Intensity CV and Data Percent Threshold both need to be provided")
            else:
                logging.info(f"Intensity CV Threshold: {cv_percent_threshold} will be applied along with data percent threshold: {data_percent_threshold}")
        else:
            logging.info("Peptide Intensities over data percentage will not be calculated, no data percent threshold is provided")

        if digestion_enzyme:
            if miscleavage_threshold:
                logging.info(f"Enzyme: {digestion_enzyme} and Missed cleavage Percentage Threshold: {miscleavage_threshold} will be used to calculated missed cleavage percentage for each sample")
            if not miscleavage_threshold:
                logging.info(f"Missed Cleavage Percentage Threshold will not be applied, no threshold provided")
        else:
            logging.info(f"Enzyme not provided, missed cleavage percentage will not be calculated")

    else:
        logging.info("Peptide Level file not provided, peptide level metrics will not be calculated")


    if precursor_level:

        logging.info("--------------------------------------------- Checking Precursor Intensity File --------------------------------------------- ")
        check_file(precursor_level, "precursor")

        if grouping_file:
            groups = check_grouping_file(precursor_level, grouping_file)
            logging.info(f"{grouping_file} will be used for comparing samples across the following groups: {groups}")

            if not tic_cv_threshold:
                logging.error("TIC CV threshold is not provided, no precursor TIC comparison between groups can be performed. Please provide TIC CV threshold")
                sys.exit(1)

        if precursor_threshold:
            logging.info(f"Precursor Threshold: {precursor_threshold} will be applied")
        else:
            logging.info("No Precursor Threshold will be applied")

        if cv_percent_threshold:
            logging.info(f"Intensity CV Threshold: {cv_percent_threshold} will be applied for precursor intensities")
        else:
            logging.info("No Intensity CVs will be calculated, no threshold provided")

        if data_percent_threshold:
            if not cv_percent_threshold:
                logging.info(f"Intensity CV and Data Percent Threshold both need to be provided")
            else:
                logging.info(f"Intensity CV Threshold: {cv_percent_threshold} will be applied along with data percent threshold: {data_percent_threshold}")
        else:
            logging.info("Precursor Intensities over data percentage will not be calculated, no data percent threshold is provided")

        if digestion_enzyme:
            if miscleavage_threshold:
                logging.info(f"Enzyme: {digestion_enzyme} and Missed cleavage Percentage Threshold: {miscleavage_threshold} will be used to calculated missed cleavage percentage for each sample")
            if not miscleavage_threshold:
                logging.info(f"Missed Cleavage Percentage Threshold will not be applied, no threshold provided")
        else:
            logging.info(f"Enzyme not provided, missed cleavage percentage will not be calculated")

        if grouping_file:
            if not tic_cv_threshold:
                logging.error("TIC CV threshold is not provided, no precursor TIC comparison between groups can be performed. Please provide TIC CV threshold")
                sys.exit(1)

    else:
        logging.info("Precursor Level file not provided, precursor level metrics will not be calculated")


    if peptide_level or precursor_level:

        if irtlabel:
            logging.info(f"iRT peptides labelled with {irtlabel} will be used for analysis")
        else:
            logging.info("No iRT QC analysis will be performed")

        if peptide_list:
            check_file(peptide_list, "Peptide List")

    if precursor_rt or peptide_rt:
        logging.info("---------------------------------------- Checking Retention Time Inputs --------------------------------------------------------")

    if precursor_rt:
        check_file(precursor_rt, "precursor")
        #coverage_threshold
    else:
        logging.info("Precursor RT file not provided, no precursor RT metrics will be calculated")

    if peptide_rt:
        check_file(peptide_rt, "peptide")
        #coverage_threshold
    else:
        logging.info("Peptide RT file not provided, no peptide RT metrics will be calculated")

    #check samples across all provided inputs
    logging.info("--------------------------------------------- Checking Samples in Inputs --------------------------------------------- ")
    check_samples(mzml_dir, protein_level, peptide_level, precursor_level, grouping_file)

    if grouping_file:
        logging.info("--------------------------------------------- Checking for Groups -----------------------------------------------------")
        groupwise_comparison = True
        groups = get_grouping_dict(grouping_file)
    else:
        groupwise_comparison = False
        groups = ""

    if mzml_dir:
        logging.info("-------------------------------------- Calculating ID Free Metrics ------------------------------------------------ ")

        mzml_threshold_dict = {}
        mzml_threshold_dict['MS1 TIC Threshold'] = ms1_tic_threshold
        mzml_threshold_dict['MS2 TIC Threshold'] = ms2_tic_threshold
        mzml_threshold_dict['MS1 Spectra Threshold'] = ms1_spectra_threshold
        mzml_threshold_dict['MS2 Spectra Threshold'] = ms2_spectra_threshold
        mzml_threshold_dict['Max Basepeak Intensity Threshold'] = max_basepeak_intensity_threshold
        mzml_threshold_dict['TIC CV Threshold'] = tic_cv_threshold

        mzml_sample_df, mzml_group_df, idfree_report_parameters = calculate_idfree_metrics(out_dir, reportname, mzml_dir, groupwise_comparison, groups, mzml_threshold_dict)

    if protein_level or peptide_level or precursor_level:
        logging.info("-------------------------------------- Calculating ID Based Metrics ------------------------------------------------ ")

        input_dict = {}
        input_dict['Protein Level'] = protein_level
        input_dict['Peptide Level'] = peptide_level
        input_dict['Precursor Level'] = precursor_level
        input_dict['Peptide RT'] = peptide_rt
        input_dict['Precursor RT'] = precursor_rt
        input_dict['Peptide List'] = peptide_list

        threshold_dict = {}
        threshold_dict['Protein Threshold'] = protein_threshold
        threshold_dict['Peptide Threshold'] = peptide_threshold
        threshold_dict['Precursor Threshold'] = precursor_threshold
        threshold_dict['Enzyme'] = digestion_enzyme
        threshold_dict['Miscleavage Threshold'] = miscleavage_threshold
        threshold_dict['TIC CV Threshold'] = tic_cv_threshold
        threshold_dict['CV Percent Threshold'] = cv_percent_threshold
        threshold_dict['Data Percent Threshold'] = data_percent_threshold
        threshold_dict['iRT Label'] = irtlabel
        threshold_dict['Coverage Threshold'] = coverage_threshold

        #calculating_idbased_metrics()
        idbased_sample_df, idbased_group_df, idbased_report_parameters = calculate_idbased_metrics(out_dir, reportname, input_dict, threshold_dict, groups, groupwise_comparison)

    if mzml_dir:
        if protein_level or peptide_level or precursor_level:
            sample_df = pd.merge(mzml_sample_df, idbased_sample_df, on="Filename")
            if groupwise_comparison:
                grouped_df = pd.merge(mzml_group_df, idbased_group_df, on="Group")
        else:
            sample_df = mzml_sample_df
            if groupwise_comparison:
                grouped_df = mzml_group_df
    elif protein_level or peptide_level or precursor_level:
        sample_df = idbased_sample_df
        if groupwise_comparison:
            grouped_df = idbased_group_df


    logging.info(f"Saving Overall QC Report to {out_dir}/{reportname}_QC_Status_Report.xlsx")
    #saving dataframes to excel document
    writer = pd.ExcelWriter(f"{out_dir}/{reportname}_QC_Status_Report.xlsx", engine='xlsxwriter')
    sample_df.to_excel(writer, index=False, sheet_name="Samplewise QC Metrics")
    if groupwise_comparison:
        grouped_df.to_excel(writer, index=False, sheet_name='Groupwise QC Metrics')
    writer.save()

    #creating report
    all_report_params = {}

    if mzml_dir:
        if protein_level or peptide_level or precursor_level:
            all_report_params = dict(tuple(idfree_report_parameters.items()) + tuple(idbased_report_parameters.items()))
        else:
            all_report_params = idfree_report_parameters
    elif protein_level or peptide_level or precursor_level:
        all_report_params = idbased_report_parameters

    all_report_params['groupwise_comparison'] = groupwise_comparison
    all_report_params['mzml_dir'] = mzml_dir

    logging.info(f"Saving HTML QC Report to {out_dir}/{reportname}.html")
    if all_report_params:
        env = Environment(loader=FileSystemLoader(str("/common/vegesnam/qcpackage/github/QCPackage/templates"))) #remove hardcoding - added this for now
        template = env.get_template(str("report_template.html"))

        output_from_parsed_template = template.render(all_report_params)

        with open(f'{out_dir}/{reportname}.html', 'w',encoding="utf-8") as f:
            f.write(output_from_parsed_template)


if __name__ == '__main__':
    main()

print("--- %s seconds ---" % (time.time() - start_time))
logging.info("--- %s seconds ---" % (time.time() - start_time))
