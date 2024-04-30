# QCeltis
A python package developed for performing quality control analysis on large-scale DIA proteomics datasets. It was designed to enable detection of technical variability due to several factors such as sample collection, transportation, storage, preparation, and/or instrument performance, thus helping improve the accuracy of any biological interpretations of large-scale mass spectrometry-based proteomics data. It allows users to monitor quality control (QC) samples within and across different batches and helps create metrics and plots to not only easily depict outliers, but also to tease out potential causes of these outliers.

## Table of Contents
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Usage](#Usage)
4. [Description of Parameters](#description-of-parameters)
5. [Input File Description](#input-file-description)
6. [Outputs and Explanations](#outputs-and-explanations)
7. [Cite](#Cite)
8. [Contributions](#contributions)
9. [Release Notes](#release-notes)

## Introduction

### Why is Quality Control Required
Quality control is essential in large-scale DIA proteomics projects to minimize variability and ensure the accuracy and reproducibilty of measurements and biological results. Technical variability, stemming from factors like sample handling and MS instrument performance, can be addressed through regular and periodic measurement of QC samples used in the experimental dataset. Teasing out the cause of outliers in the category of technical variation can be done by evaluating some data metrics which vary depending on MS instrument, LC column, instrument cleaning, length of proteolytic digestion, and sample cleanup, and can improve downstream analysis. QCeltis aims to aid in performing comprehensive quality control analysis, further enhancing the reliability of proteomics research results. 

### QCeltis Overview

QCeltis a command-line python package that provides users with quality control assessement along with easy-to-intepret visualizations to help idenitfy issues at sample preparation, digestion and acquisition steps of the MS experiment process. To do so, QCeltis provides a way to perform quality control analysis through 2 primary steps: 

1) ID-Free: parses input consisting of mzML files and extracts raw data metrics such as: MS1 and MS2 TIC, MS2/MS1 Spectral Ratio, Max Base Peak Intensities. During this step, the script applies user-defined thresholds, if given and outlier analysis to identify "bad" quality control samples that are not consistent with the rest of the dataset. 
2) ID-Based: parses input from data processed through a search engine, such as Dia-NN or OpenSWATH. The input file contain protein, peptide or precursor-level intensities, or all 3 can be provided. Metrics such as quantification information, intensity CVs, common peptide/precursor TIC, miscleaved peptides count, etc are extracted and user-defined thresholds are applied to identify "bad" samples.

After extracting ID-Free and ID-Based metrics and identifying outliers, the script provides samplewise "PASS/FAIL" labels which are reported in the output Excel Reports. Additionally, a sample grouping file can be provided by the user which is used to perform groupwise (either plates or batches) quality control analysis across the provided quality control samples. Additional plots and reports are generated when the grouping file is provided and a groupwise "PASS/FAIL" label is reported in the outputs. 

- Refer to the [Input File Description](#input-file-description) section for further details on required input files. 
- Refer to the [Outputs and Explanations](#outputs-and-explanations) section for further details on metrics evaluated and graphical representation of the metrics assessed through QCeltis 

![image](https://github.com/vegesnam/QCeltis/assets/87665957/3d084e98-8f3c-472a-959a-af331ef8947f)

#### Metrics: 

| Metric                                | Category | Explanation                                                                                                           |
|---------------------------------------|----------|--------------------------------------------------------------------------------------------------------------------------|
| MS1 TIC                               | ID-Free  | Total Ion Current extracted from MS1 Spectra. MS1 TIC values are expected to be consistent across replicate quality control samples, inconsistencies could point to issues with data acquisition and improper LC-MS instrument performance.                            |                        
| MS2 TIC                               | ID-Free  | Total Ion Current extracted from MS2 Spectra. MS2 TIC values are expected to be consistent across replicate quality control samples, inconsistencies could point to issues with data acquisition and improper LC-MS instrument performance.                            |
| MS2/MS1 Spectral Ratio                | ID-Free  | The ratio of MS2-to-MS2 Spectra provides insights into the depth of proteome coverage and fragmentation efficiency. Consistent or expected numbers of spectra across samples indicates stable instrument performance, any significant deviations can indicate issues with the MS instrument and deviation from the experimental protocol.                         |
| Max Base peak Intensity               | ID-Free  | Represents the highest base peak intensity in each replicate sample. Extreme maximum base peak intensities within samples can indicate issues with the data acquisition and LC-MS instrument performance.                     |
| Quantification                        | ID-Based | The number of identified proteins, peptides and precursors is a fundamental metric for assessing the data quality from the experiment. Inconsistent numbers across replicate samples points to issues with sample preparation, digestion protocols, experimental reproducibility, and data acquisition                                                        |
| Intensity CVs                         | ID-Based | Coefficient of Variation (CV%) across the dataset or within provided groups reveals the degree of variability. Higher CVs indicate greater variation in the intensity measurements across replicate samples stemming from sample preparation, data acquisition or instrument performance. |
| PCA                                   | ID-Based | PCA plots can reveal the percentage of variability of provided replicate sample intensities and potential batch effects across the provided batches or plates. |
| Common TIC                            | ID-Based | TIC values from common peptides or precursors across the replicate samples aids in assessing the reproducibility. Variations in TIC for common peptides or precursors can highlight potential sources of variation, such as sample preparation or contamination issues. |
| 0 Missed Cleavage Percentage          | ID-Based | Used to assess the efficiency and specificity of the enzymatic digestion during the sample preparation step.                                           |
| iRT Peptide Intensity Distribution    | ID-Based | iRT (Indexed Retention Time) peptides can be used as a reference for retention time normalization and calibration in MS proteomics. Monitoring the intensities and presence of these iRT peptides can reveal any sample preparation or data acquisition issues.           |
| Selected Peptide Intensity Distribution | ID-Based | Monitoring the intensities and presence of user-defined peptides can reveal any sample preparation or data acquisition issues.         |

## Installation
QCeltis can be installed directly from [PyPI](https://pypi.org/project/qceltis/) using `pip` as follows. 

~~~bash
pip install qceltis
~~~
Alternatively, if you would like the development version from github. You can perform the below steps to use the source. 

First, clone the repository from GitHub:

~~~bash
git clone https://github.com/csmc-vaneykjlab/QCeltis.git
~~~

Then, navigate to the downloaded directory and install the package requirements using `pip`:

~~~bash
cd QCeltis
pip install -r requirements.txt
~~~

Check your installation:
~~~bash
python3 main.py --help
~~~



## Usage

If downloaded using github: 

```python
python  main.py [-h] --outdirectory OUTDIRECTORY --reportname REPORTNAME [--mzml_directory MZML_DIRECTORY]
               [--ms1_tic_threshold MS1_TIC_THRESHOLD] [--ms2_tic_threshold MS2_TIC_THRESHOLD]
               [--ms1_spectra_threshold MS1_SPECTRA_THRESHOLD] [--ms2_spectra_threshold MS2_SPECTRA_THRESHOLD]
               [--max_basepeak_intensity MAX_BASEPEAK_INTENSITY] [--protein_level PROTEIN_LEVEL]
               [--peptide_level PEPTIDE_LEVEL] [--precursor_level PRECURSOR_LEVEL]
               [--grouping_file GROUPING_FILE] [--peptide_list PEPTIDE_LIST]
               [--protein_threshold PROTEIN_THRESHOLD] [--peptide_threshold PEPTIDE_THRESHOLD]
               [--precursor_threshold PRECURSOR_THRESHOLD] [--enzyme ENZYME] [--miscleavage_threshold MISCLEAVAGE_THRESHOLD]
               [--tic_cv_threshold TIC_CV_THRESHOLD] [--cv_percent_threshold CV_PERCENT_THRESHOLD]
               [--data_percent_threshold DATA_PERCENT_THRESHOLD] [--irt IRTLABEL]
               [--coverage_threshold COVERAGE_THRESHOLD]
```

If downloaded as a library through pip install:

```python
qceltis [-h] --outdirectory OUTDIRECTORY --reportname REPORTNAME [--mzml_directory MZML_DIRECTORY]
               [--ms1_tic_threshold MS1_TIC_THRESHOLD] [--ms2_tic_threshold MS2_TIC_THRESHOLD]
               [--ms1_spectra_threshold MS1_SPECTRA_THRESHOLD] [--ms2_spectra_threshold MS2_SPECTRA_THRESHOLD]
               [--max_basepeak_intensity MAX_BASEPEAK_INTENSITY] [--protein_level PROTEIN_LEVEL]
               [--peptide_level PEPTIDE_LEVEL] [--precursor_level PRECURSOR_LEVEL]
               [--grouping_file GROUPING_FILE] [--peptide_list PEPTIDE_LIST]
               [--protein_threshold PROTEIN_THRESHOLD] [--peptide_threshold PEPTIDE_THRESHOLD]
               [--precursor_threshold PRECURSOR_THRESHOLD] [--enzyme ENZYME] [--miscleavage_threshold MISCLEAVAGE_THRESHOLD]
               [--tic_cv_threshold TIC_CV_THRESHOLD] [--cv_percent_threshold CV_PERCENT_THRESHOLD]
               [--data_percent_threshold DATA_PERCENT_THRESHOLD] [--irt IRTLABEL]
               [--coverage_threshold COVERAGE_THRESHOLD]
```

## Parameters

| Parameter                | Short Form | Description                                      | Default Value |
|--------------------------|------------|--------------------------------------------------|---------------|
| --outdirectory           | -o         | Output directory path                            | None          |
| --reportname             | -r         | Report name for HTML and Excel reports           | None          |
| --mzml_directory         | -m         | Path to the directory where mzML files are present | None          |
| --ms1_tic_threshold      | -t1        | MS1 TIC threshold                                | False         |
| --ms2_tic_threshold      | -t2        | MS2 TIC threshold                                | False         |
| --ms1_spectra_threshold  | -s1        | MS1 spectra threshold                            | False         |
| --ms2_spectra_threshold  | -s2        | MS2 spectra threshold                            | False         |
| --max_basepeak_intensity | -bp        | Maximum base peak intensity threshold             | False         |
| --protein_level          | -pt        | Path to protein intensity file                   | None          |
| --peptide_level          | -pep       | Path to peptide intensity file                   | None          |
| --precursor_level        | -pre       | Path to precursor intensity file                 | None          |
| --grouping_file          | -g         | Path to grouping file                            | None          |
| --peptide_list           | -peplt     | Path to file containing list of peptides to monitor intensity and RT distribution across samples | None          |
| --protein_threshold      | -x         | Protein threshold for each sample                | False         |
| --peptide_threshold      | -y         | Peptide threshold for each sample                | False         |
| --precursor_threshold    | -z         | Precursor threshold for each sample              | False         |
| --enzyme                 | -e         | User input enzyme                                | None          |
| --miscleavage_threshold  | -c         | Missed cleavage threshold for each sample        | False         |
| --tic_cv_threshold       | -t         | TIC CV threshold for groupwise QC status         | 40         |
| --cv_percent_threshold   | -s         | Intensity CV threshold                           | 40         |
| --data_percent_threshold | -d         | Data threshold for intensity CV                  | 70         |
| --irtlabel               | -irt       | Label for iRT peptides present in your peptide intensity file | None          |
| --coverage_threshold     | -v         | Intensity or retention time coverage % threshold in each sample | False         |

## Input Files Descriptions

### Input for Raw Data quality (ID-Free) Assessment: 

#### mzml_directory

A directory containing mzML files needs to be given as input. Please refer https://github.com/HUPO-PSI/mzML for specifications of the mzML file format. The mzML files, should contain data as spectra, and can be converted from .raw (for data from a Thermo-Fisher instrument) or .wiff (for raw data from a Sciex instrument). The script will only pick the files with .mzML extention from the input directory.

### Input for Search Results (ID-Based) Assessment: 

#### protein_level input file ([Example](https://github.com/csmc-vaneykjlab/QCeltis/blob/main/example-dataset/protein_level.txt))

| Column                | Description |
|--------------------------|------------|
| Protein           | A column containing Protein Ids/Names      |
| < list of samples >             | a list of samples followed by their intensity values     |

#### peptide_level input file

| Column                | Description |
|--------------------------|------------|
| Protein           | A column containing Protein Ids/Names      |
| Peptide           | A column containing Peptide Ids      |
| < list of samples >             | a list of samples followed by their intensity values     |

#### precursor_level input file ([Example](https://github.com/csmc-vaneykjlab/QCeltis/blob/main/example-dataset/precursor_level.txt))

| Column                | Description |
|--------------------------|------------|
| Protein           | A column containing Protein Ids/Names      |
| Peptide           | A column containing Peptide Ids      |
| Precursor           | A column containing Precursor Ids      |
| < list of samples >             | a list of samples followed by their intensity values     |

### For Batch-wise QC Comparison: 

#### Sample Grouping File ([Example](https://github.com/csmc-vaneykjlab/QCeltis/blob/main/example-dataset/grouping_file.txt))

| Column                | Description |
|--------------------------|------------|
| Filename           | A column containing file names      |
| Group            | The batch the sample belongs to     |

## Outputs and Explanations

### ID-Free and ID-Based Metrics:

| Metric                                | Category | Expected Output                                                                                                           |
|---------------------------------------|----------|--------------------------------------------------------------------------------------------------------------------------|
| MS1 TIC                               | ID-Free  | TIC Values, MS1 + MS2 TIC Graph, Outlier Detection, Threshold FAIL/PASS, TIC CV% across groups                            |
| MS2 TIC                               | ID-Free  | TIC Values, MS1 + MS2 TIC Graph, Outlier Detection, Threshold FAIL/PASS, TIC CV% across groups                            |
| MS2/MS1 Spectra                       | ID-Free  | TIC Values, MS2/MS1 Graph, Outlier Detection on MS2/MS1, Threshold FAIL/PASS                                             |
| Max Base peak Intensity                | ID-Free  | Max Base peak Values, Max Base peak Intensity Graph + Outlier Detection                                                        |
| Quantification                        | ID-Based | Protein/Peptide/Precursor Quantification + Threshold PASS/FAIL + graph                                                            |
| Intensity CVs                         | ID-Based | Protein/Peptide/Precursor Overall CV + Cumulative CV% + Groupwise CV % + Cumulative/Groupwise Graphs + Threshold PASS/FAIL (if groupwise) |
| Common TIC                            | ID-Based | Peptide/Precursor Common TIC values + TIC CV% across groups + graphs for samples + groups + Threshold PASS/FAIL (if groupwise) |
| 0 Missed Cleavage Percentage          | ID-Based | Missed Cleavage Summary + Threshold PASS/FAIL                                                                           |
| iRT Peptide Intensity Distribution    | ID-Based | Peptide/Precursor Level Intensities Distribution + Graph + Coverage Summary across samples + Threshold PASS/FAIL         |
| Selected Peptide Intensity Distribution | ID-Based | Peptide/Precursor Level Intensities Distribution + Graph + Coverage Summary across samples + Threshold PASS/FAIL         |

### Excel reports

| Report Type        | Sheets                                                  |
|--------------------|---------------------------------------------------------|
| ID-Free            | ID-Free Metrics Summary, Group TIC CV                   |
| Protein Level      | Protein Quant Summary, Groupwise Protein Quant, Protein Level CV, Protein CV Group Summary |
| Peptide Level      | Peptide Quant Summary, Groupwise Peptide Quant, Peptide Level CV, Peptide CV Group Summary, Common Peptide TIC, Common Peptide TIC Group CV, Miscleavage Threshold, iRT Peptide Intensity/Selected Peptide Intensity | 
| Precursor Level    | Precursor Quant Summary, Groupwise Precursor Quant, Precursor Level CV, Precursor CV Group Summary, Common Precursor TIC, Common Precursor TIC Group CV, Miscleavage Threshold, iRT Precursor Intensity/Selected Precursor Intensity | 
| Status Report      | Samplewise QC Metrics, Groupwise QC Metrics             |

For examples of Excel Reports: [ID-Free + ID-Based Example Dataset Reports](https://github.com/vegesnam/QCeltis/tree/main/example-dataset/results/IDFree_IDBased_GroupComparison_QCResult)
(Manasa: change this once the repo is changed)

### HTML Report

The HTML report is generated at the end of the QCeltis analysis. The report is divided into 2 tabs: ID-Free Metrics and ID-Based Metrics. If mzML directory is provided as input, ID-Free Tab will be populated with plots using the Id-Free metrics. The Id-Based Tab will contain metrics and results from your Search Engine Results dataset - Protein, Peptide or Precursor Level intensity files.

Note: If grouping file is provided, additional groupwise plots are produced within the report. The colors within each graph will represent the groups provided. 

#### ID-Free Metrics Tab

QCeltis extracts the following ID-Free metrics from the provided mzML files and detects samples that are not consistent across the dataset or are outliers using outlier analysis methods and user-defined thresholds. 

#### Total Ion Current

<ins>Total Ion Current Line Graph:</ins>

![TIC Line Graph](https://github.com/vegesnam/QCeltis/assets/87665957/862a19ee-12e0-4fd4-bd55-e3b9305e882d)

MS1 and MS2 Total Ion Current Values are extracted from spectra present in the mzML files and are plotted using a line graph. Total Ion Current values from MS1 and MS2 Spectra are expected to be consistent across the replicate quality control samples. If 'MS1 TIC Threshold' or 'MS2 TIC Threshold' is provided by the user, a dotted line is used to display the threshold in the line graph. Any samples above that threshold are considered as failed samples. Apart from threshold-defined identification of failed samples, outlier analysis is also performed. If extreme values are found across the samples, they are labelled as outliers. 

<ins>TIC Outlier Plot:</ins> 

![MS1 TIC Outlier Plot](https://github.com/vegesnam/QCeltis/assets/87665957/b73575cd-c6da-48cc-9444-8397b7a681c4)

When outliers are identified, a scatter plot of the TIC values is produced, with the outlier highlighted in yellow. Outliers/Failed Samples found using this metric indicate issues with with data acquisition and LC-MS instrument performance such as improper autosampler sample pickup. 
 
<ins>TIC CV% Across Groups Plot:</ins> 

![MS1 TIC CV% Plot](https://github.com/vegesnam/QCeltis/assets/87665957/6447f8ed-c446-4772-8693-333eef06d746)

When a grouping file is provided, TIC CV% is calculated across samples within each provided group and the 'TIC CV Threshold' is applied. If any group isn't within the threshold, the group is labelled as a 'FAIL' and this indicates an inconsistent TIC pattern within the samples of the group. 

#### Spectral Ratio

<ins>Spectral Ratio Line Graph:</ins> 

![Spectral Ratio](https://github.com/vegesnam/QCeltis/assets/87665957/911c43df-38d9-4d01-af67-537848dc038c)

The number of MS1 and MS2 spectra are counted from each provided mzML file and the MS2-to-MS1 spectral ratio is calculated. Consistent and expected numbers of spectra across samples indicates stable instrument performance. Outlier analysis is performed to check for any extreme values and identify outliers. If outliers are identified, a scatter plot (similar to the above TIC Outlier Plot) is generated, with the outliers highlighted in yellow. Outliers/Failed Samples found using this metric indicate issues with with LC-MS instrument performance or experiment protocol changes. 

#### Max Base Peak Intensity

<ins>Max Base Peak Intensity Bar Graph:</ins>

![Max Base Peak Intensity](https://github.com/vegesnam/QCeltis/assets/87665957/4e8881ee-95a0-4318-96fe-13781fae0e59)

When the provided mzML files have base peak intensity information, the maximum base peak intensity is extracted from the files and plotted as a bar graph. The Max Base Peak Intensity represents the highest recorded base peak intensity in each mzML file. If 'Max Base Peak Intensity Threshold' is provided by the user, a dotted line is used to display the threshold in the line graph. Any samples above that threshold are considered as failed samples. Apart from threshold-defined identification of failed samples, outlier analysis is also performed. If extreme values are found across the samples, they are labelled as outliers through outlier analysis. If outliers are identified, a scatter plot (similar to the above TIC Outlier Plot) is generated, with the outliers highlighted in yellow. Outliers/Failed Samples found using this metric indicate issues with with data acquisition or LC-MS instrument performance. Outliers identified within this metric usually correlate with outliers identified within the MS1 or MS2 TIC metrics.  

#### ID-Based Metrics Tab 

QCeltis extracted the following ID-Based metrics from the provided search engine results, applies user-defined thresholds and performs groupwise comparison if the grouping file is provided. Depending on the type of search engine input, different plots are generated. 

#### Quantification Plots

Based on the type of search engine result, the number of proteins, peptides or precursors is extracted. Here, examples of number of proteins and precursors identified is provided. 

<ins>Quantification Bar Graphs:<ins>

<ins>Protein Bar Graph:</ins> 

![Number of Proteins](https://github.com/vegesnam/QCeltis/assets/87665957/009c4de5-cb94-4daf-8bc9-9c33375a96ab)

Number of Proteins identified in each sample is plotted. Here, a 'Protein Threshold = 200' was provided and it is represented using a dotted line. All the samples are considered as "PASS" since the number of proteins is above the provided threshold. If there are any samples that have proteins numbers below the threshold, they are considered as "FAIL" and this indicates issues with sample preparation, digestion protocols, experimental reproducibility

<ins>Precursor Bar Graph:</ins> 

![Number of Precursors](https://github.com/vegesnam/QCeltis/assets/87665957/de40c702-105b-4f33-a7c6-eba1d1fbd2c5)

Number of Precursors identified in each sample is plotted. Here, a 'Precursor Threshold = 3000' was provided and it is represented using a dotted line. All the samples are considered as "PASS" since the number of precursors is above the provided threshold. If there are any samples that have precursors numbers below the threshold, they are considered as "FAIL" and this indicates issues with sample preparation, digestion protocols, experimental reproducibility. 

When a peptide-level input is provided, a bar plot similar to the above is produced. 

#### Intensity CV% Plots

<ins>Intensity CV% Cumulative Frequency Line Graph:</ins>

![CV Cumulative Frequency](https://github.com/vegesnam/QCeltis/assets/87665957/8c341781-5627-467e-aa64-9bf7106df54b)

Using the intensity values provided, CV% (Coefficient of Variation %) is calculated for each protein, peptide or precursor and the CV% cumulative frequency is plotted using a line graph. Here, the protein and precursor CV% cumulative frequency is plotted. Higher CVs indicate greater variation in the intensity measurements across replicate samples stemming from sample preparation, data acquisition or instrument performance. Lower CVs indicate higher reproducibility of protein, peptide or precursor intensities across replicate samples. 

<ins>Intensity CV% Bar Graph:</ins> 

![Protein CV%](https://github.com/vegesnam/QCeltis/assets/87665957/393fe0cc-8e2b-47f1-83c9-c1d058f098f4)

![Precursor CV%](https://github.com/vegesnam/QCeltis/assets/87665957/a09a7943-d8b2-462b-91c3-7cad97a157cf)

When a grouping file is provided, CV% is calculated across samples within each provided group and a bar graph is generated. Each bar represents the percentage of proteins, peptide or precursors under the 'CV Percent Threshold' (By default, cv percent threshold is set to 40%, this can be changed by the user). This provides information about the consistency and reproducibility of intensity values within each provided group (batch/plate). Groupwise CV% helps identify batch/plate specific trends and patterns. A 'Data Percent Threshold' (By default, data percent threshold is set to 70%, this can be changed by the user) is applied and any groups not passing this threshold are labelled as "FAIL". 

Here, protein and precursor CV% across groups are shown. A similar plot will be generated for peptides as well (when peptide-level input is provided). 

#### Common TIC

<ins>Common Peptide/Precursor TIC Bar Graph:</ins> 

![Common Precursor TIC](https://github.com/vegesnam/QCeltis/assets/87665957/4f74e742-f41a-410b-b725-b8f3d130779d)

<ins>Common Peptide/Precursor TIC CV%:</ins> 

![Common Precursor TIC CV%](https://github.com/vegesnam/QCeltis/assets/87665957/e10960a0-50b9-4d02-9030-6d86da117177)

Common Peptide/Precursor TIC is the summed intensity of all the common peptides or precursors found in all samples. If both peptide and precursor files are provided, only common peptide will be calculated.

If TIC is consistent with other samples here but the TIC is labelled as an outlier in ID-free, then it could be a contamination issue (if TIC is high in id-free)
CVs are calculated with intensity values from samples within each group.

#### Miscleavages Plot

<ins>Number of No Miscleavages Bar Graph:</ins> 

![image](https://github.com/csmc-vaneykjlab/QCeltis/assets/87665957/617e6bcc-53fe-4e65-b5c1-3fe13a65bca4)

The number of miscleavages present in each peptide is calculated and the number of peptides with no miscleavages in each sample is plotted in a bar graph. If the 'Miscleavage Threshold' is provided, a dotted line is plotted across the bars. Any samples not passing the threshold are labelled as "FAIL". Having a high number of peptides with 0 miscleavages indicates a good digestion protocol. Any samples not meeting the threshold point to the inefficiency of the enzymetic digestion during the sample preparation step. 

##### Other Plots:

##### PCA:
When a grouping file is given, a PCA plot is generated using the protein, peptide or precursor intensity values. Replicate quality control samples are expected to cluster together, if separation is observed among the provided groups, this could indicate possible batch effects across the dataset.  

##### iRT / Selected Peptide Intensity Distribution + Coverage Summary:

If your dataset contains iRT peptides, the '--irtlabel' parameter can be used to plot the intensity distribution of the iRT peptides/precursors (from the iRT Biognosys) across the dataset using a line graph. Apart from the intensity distribution, a coverage summary bar graph is generated, where each bar represents the percentage of samples each iRT peptide/precursor is present in. If the 'Coverage Threshold' is provided by the user, a dotted line is also plotted in the bar graph. Any samples not meeting the threshold are labelled as "FAIL". Low coverage of iRT peptides across the datasets could indicate issues with sample preparation. 

Using the '--peptide_list' parameter, a user-defined list of peptides can be monitored instead of the iRT peptides. A similar intensity distribution line graph and coverage summary bar graphs will be plotted and the 'Coverage Threshold' is applied if it is provided by the user. We recommend monitoring [CiRT peptides](https://www.sciencedirect.com/science/article/pii/S1535947620326335) for eukaryotic datasets. 

## Cite

## Support

If you encounter any bugs or issues, please help us improve QCeltis by creating a new issue at: https://github.com/csmc-vaneykjlab/QCeltis/issues. For any other queries, email us at GroupHeartBioinformaticsSupport@cshs.org.

## Release Notes

Version 0.0.1
 - Initial version released.

## License
QCeltis © 2024 by Jennifer Van Eyk is licensed under [MIT License](https://github.com/csmc-vaneykjlab/QCeltis/blob/main/LICENSE.txt).
