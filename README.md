# QCPackage
Initial Code Source for the QC Package

## Table of Contents
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Usage](#Usage)
4. [Description of Parameters](#description-of-parameters)
5. [Input File Description](#input-file-description)
6. [Outputs and explanations](#outputs-and-explanations)
7. [Cite](#Cite)
8. [Contributions](#contributions)
9. [Release Notes](#release-notes)

## Introduction

### Why is Quality Control Required
The sources of variability in a proteomics experiment that are addressed by QC can be categorized into two groups: biological and technical. Technical variability is derived from sample collection, transportation, storage, preparation, and/or instrument performance. Teasing out the cause of outliers in the category of technical variability can be done by evaluating some data parameters which can vary depending on the mass spectrometer, LC column, instrument cleaning, length of proteolytic digestion, and sample cleanup, and can improve downstream analysis. This is the aim of the QCPackage.


## Installation
QCP can be installed from the command line using `git` and `pip`.

First, clone the repository from GitHub:

git clone https://github.com/vegesnam/QCPackage.git

Then, navigate to the downloaded directory and install the package using `pip`:

~~~bash
cd qcp
pip install -r requirements.txt
~~~

Alternatively, the package can be installed directly from PyPI:

`pip install -i < pypi instance will go here >`

## Usage

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
| --tic_cv_threshold       | -t         | TIC CV threshold for groupwise QC status         | False         |
| --cv_percent_threshold   | -s         | Intensity CV threshold                           | False         |
| --data_percent_threshold | -d         | Data threshold for intensity CV                  | False         |
| --irtlabel               | -irt       | Label for iRT peptides present in your peptide intensity file | None          |
| --coverage_threshold     | -v         | Intensity or retention time coverage % threshold in each sample | False         |


## Input Files Descriptions

### protein_level input file

| Column                | Description |
|--------------------------|------------|
| Protein           | A column containing Protein names      |
| < list of samples >             | a list of samples followed by their intensity values     |

see example: and link

### peptide_level input file

| Column                | Description |
|--------------------------|------------|
| Peptides           | A column containing Peptides names      |
| < list of samples >             | a list of samples followed by their intensity values     |

see example: and link

### precursor_level input file

| Column                | Description |
|--------------------------|------------|
| Precursors           | A column containing Precursors names      |
| < list of samples >             | a list of samples followed by their intensity values     |

see example: and link

### Grouping input file

| Column                | Description |
|--------------------------|------------|
| Filename           | A column containing file names      |
| Group            | The group the sample belongs to     |

see example: and link

### mzml_directory

A directory containing mzML files needs to be given as input. Please refer https://github.com/HUPO-PSI/mzML for specifications of the mzML file format.
The script will only pick the files with .mzML extention from the input directory.

## Outputs and explanantions

### ID-Free and ID-Based Metrics:

| Metric                                | Category | Expected Output                                                                                                           |
|---------------------------------------|----------|--------------------------------------------------------------------------------------------------------------------------|
| MS1 TIC                               | ID-Free  | TIC Values, MS1 + MS2 TIC Graph, Outlier Detection, Threshold FAIL/PASS, TIC CV% across groups                            |
| MS2 TIC                               | ID-Free  | TIC Values, MS1 + MS2 TIC Graph, Outlier Detection, Threshold FAIL/PASS, TIC CV% across groups                            |
| MS2/MS1 Spectra                       | ID-Free  | TIC Values, MS2/MS1 Graph, Outlier Detection on MS2/MS1, Threshold FAIL/PASS                                             |
| Max Basepeak Intensity                | ID-Free  | Max Basepeak Values, Basepeak Intensity Graph + Outlier Detection                                                        |
| Quant                                 | ID-Based | Protein/Peptide/Precursor Quant + Threshold PASS/FAIL + graph                                                            |
| Intensity CVs                         | ID-Based | Protein/Peptide/Precursor Overall CV + Cumulative CV% + Groupwise CV % + Cumulative/Groupwise Graphs + Threshold PASS/FAIL (if groupwise) |
| Common TIC                            | ID-Based | Peptide/Precursor Common TIC values + TIC CV% across groups + graphs for samples + groups + Threshold PASS/FAIL (if groupwise) |
| 0 Missed Cleavage Percentage          | ID-Based | Missed Cleavage Summary + Threshold PASS/FAIL                                                                           |
| iRT Peptide Intensity Distribution    | ID-Based | Peptide/Precursor Level Intensities Distribution + Graph + Coverage Summary across samples + Threshold PASS/FAIL         |
| Selected Peptide Intensity Distribution | ID-Based | Peptide/Precursor Level Intensities Distribution + Graph + Coverage Summary across samples + Threshold PASS/FAIL         |
| RT Coverage Summary Distribution      | ID-Based | if RT files are given, iRT or selected peptide RT distribution + coverage summary + Threshold PASS/FAIL                  |

### Excel reports

| Report Type        | Sheets                                                  |
|--------------------|---------------------------------------------------------|
| ID-Free            | ID-Free Metrics Summary, Group TIC CV                   |
| Protein Level      | Protein Quant Summary, Groupwise Protein Quant, Protein Level CV, Protein CV Group Summary |
| Precursor Level    | Precursor Quant Summary, Groupwise Precursor Quant, Precursor Level CV, Precursor CV Group Summary, Common Precursor TIC, Common Precursor TIC Group CV, Miscleavage Threshold, iRT Precursor Intensity |
| Status Report      | Samplewise QC Metrics, Groupwise QC Metrics             |

### HTML Reports

#### ID-Free plots

##### Total Ion Current

The total ion current (TIC) is the summed intensity across the entire range of masses being detected in each sample. MS1 and MS2 Total Ion Current Values extracted from spectra within the given mzML files. Total Ion Current values from MS1 and MS2 Spectra are expected to be consistent across the QC samples. If extreme values are found, they are labelled as outliers. Outliers found in TIC can indicate an issue If groups are provided, TIC CV% is calculated across samples within each group and the tic cv threshold is applied if provided.  If any group doesn’t pass the threshold, this indicates an inconsistent TIC pattern within the samples of the group.

![TIC_CV](https://github.com/vegesnam/QCPackage/assets/32958585/960d1ac5-491b-43fc-a0a5-43d91e51f766)

![TIC_Outlier](https://github.com/vegesnam/QCPackage/assets/32958585/f0102ca2-5adb-4d56-b443-ec7b4bd09f49)

![total_ion_current](https://github.com/vegesnam/QCPackage/assets/32958585/2daa3e5c-e486-4547-bacf-500802184db9)

If extreme values are found, they are labelled as outliers. Outliers found in TIC can indicate an issue If groups are provided, TIC CV% is calculated across samples within each group and the tic cv threshold is applied if provided.  If any group doesn’t pass the threshold, this indicates an inconsistent TIC pattern within the samples of the group.

##### Spectral Ratio

![Spectral_ratio](https://github.com/vegesnam/QCPackage/assets/32958585/79bef057-7082-47d0-9823-8e24879ef853)

The spectral ratio indicates the count of MS2 Spectra relative to MS1 Spectra in each file.

The mass spectrometer consistently tallies both MS1 and MS2 Counts; thus, if there are discrepancies in the numbers, it's not a sample problem but rather an issue with the MS acquisition. If the ratio isn't stable, evaluate the instrument's performance on samples identified as anomalies.

##### Max Basepeak Intensity

![max_basepeak_intensity](https://github.com/vegesnam/QCPackage/assets/32958585/3921d0e6-d8be-4f9a-8f54-3d6f41fee59b)

The base peak intensity is the recorded intensity of the most intense peak from each spectrum in the mzML file. The Max Base Peak Intensity represents the highest recorded base peak intensity in each mzML file.

![Uploading max_basepeak_intensity_outlier_plot.png…]()

Max Base Peak Intensity is expected to be consistent across QC samples. Any outliers detected are highlighted in yellow. Outliers detected could point to issues with sample pickup or samples being dried out. Usually correlates with ID-free TIC.
Note: The outlier here will also be labelled as an outlier in the TIC graph.


#### ID-Based plots

##### Quantification Plots

Displays the number of proteins and peptides or precursors. When the threshold is provided, any samples not meeting the threshold is flagged as a ‘FAIL’, this could indicate an issue with the sample prep or digestion protocol.

<p align="center">
  <!-- First Row -->
  <img src="https://github.com/vegesnam/QCPackage/assets/32958585/5dbd3332-6dba-46e2-93d9-20ea261ca8b2" width="50%" />
  <img src="https://github.com/vegesnam/QCPackage/assets/32958585/35f1128e-1d24-4951-8b3e-ebc748d03095" width="50%" />

  <!-- Second Row -->
  <img src="[URL_to_Image3](https://github.com/vegesnam/QCPackage/assets/32958585/1780abe6-b0a4-4f20-a775-0ebaa031b5ea)" width="50%" />
  <img src="[URL_to_Image4](https://github.com/vegesnam/QCPackage/assets/32958585/2fcd40f4-3b57-4fcd-8909-7405c13f0948)" width="50%" />
</p>

In this case, even though the threshold is met, the last 3 groups have significantly less proteins and precursors. Further investigation showed a high percentage of missingness in the samples and indicated an issue with sample prep.
So these plots can be a good way to visualize any varying patterns within samples

##### CV Plots

![CV_IDbased_Protein](https://github.com/vegesnam/QCPackage/assets/32958585/72cc5220-ca01-49c0-997d-8c1e1bae96c1)

![CV_IDbased_Precursor](https://github.com/vegesnam/QCPackage/assets/32958585/57060297-0e84-48e3-b5ae-6ebc7f26c6d2)

CVs are calculated with intensity values from samples within each group [ Please note: groupwise input must be given ]. Data and CV Thershold: simple explanation.

If inconsistent ( data threshold failure ), then the samples were poorly digested (unless TIC suggests otherwise). If TIC is high but the CV isn't, then digestion issue or if TIC is also low, then the sample didn't get picked up properly.
[ write better explanations for this ]

##### Common TIC

![IDbased_Common_TIC](https://github.com/vegesnam/QCPackage/assets/32958585/48f71adf-ae7f-471b-8607-2218da6bd0a8)

Common Peptide/Precursor TIC is the summed intensity of all the common peptides or precursors found in all samples. If both peptide and precursor files are provided, only common peptide will be calculated.

If TIC is consistent with other samples here but the TIC is labelled as an outlier in ID-free, then it could be a contamination issue (if TIC is high in id-free)

![IDbased_common_TIC_CV](https://github.com/vegesnam/QCPackage/assets/32958585/e28bb8dd-660c-46f5-8109-1fb63917a6e8)

CVs are calculated with intensity values from samples within each group.

##### Miscleavage Plot

![miscleavage_TR](https://github.com/vegesnam/QCPackage/assets/32958585/217d48be-94cf-479c-9c8d-1370f4d6b450)

Having a high number of peptides with 0 miscleavages indicates a good digestion protocol. If the samples don’t have enough 0 missed cleaved peptides, (doesn’t meet threshold), then this indicates a digestion issue.

![miscleavage_DR](https://github.com/vegesnam/QCPackage/assets/32958585/98746b87-073c-4f4d-93eb-da77731e7131)

Here, there is a difference in numbers and this is due to the missingness in the last 3 plates (also indicated in the protein and peptide quant) – further indicates that it is a digestion/sampleprep issue.

##### Other Plots:

###### PCA:
if groups are given, QC Samples from groups are expected to cluster together, if separation is observed, this could indicate an issue with mass spec (if using TR samples) or sampleprep/digestion protocols (if using DR samples) or could be batch effects

###### Irt / Selected Peptide Intensity Distribution + Coverage Summary:

Intensity Distribution for peptides is plotted. If iRT is selected, the script will look for peptides from the iRT Biognosys kit. If precursor file is given, it will select iRT precursors with charge 2. If both peptide and precursor input files are given, then only peptides will be plotted.If no iRTs were used, users can provide a list of peptides to plot. Intensity coverage of these peptides is also plotted. CiRT peptide list can be used for eukaryotic datasets: link to Identification of a Set of Conserved Eukaryotic Internal Retention Time Standards for Data-independent Acquisition Mass Spectrometry - PubMed (nih.gov) 



## Cite

## Support

## release notes
