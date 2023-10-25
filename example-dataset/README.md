# Using the QCPackage using an example dataset

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

## Parameters and Commands

Parameters (link to parameters in the main readme)
Input File Descriptions (link to input file descriptions)

To begin with analysis, you can download the mzML files, protein and precursor level files from (example_dataset link). You can provide these as input parameters and execute the following command:

ID-Free + ID-Based + Groupwise Comparison Command:

```python
python  main.py --outdirectory ./output --reportname ExampleDataset-QCReport --mzml_directory ./mzML_files
               --protein_level ./protein_level.txt --precursor_level ./precursor_level.txt --grouping_file ./grouping_file.txt
               --protein_threshold 200 --precursor_threshold 3000 --enzyme 'trypsin'
               --miscleavage_threshold 80 --tic_cv_threshold 40 --cv_percent_threshold 40
               --data_percent_threshold 70 --irtlabel 'iRT'--coverage_threshold 80
```

ID-Free + Groupwise Comparison Command:

```python
python  main.py --outdirectory ./output --reportname ExampleDataset-QCReport --mzml_directory ./mzML_files
               --grouping_file ./grouping_file.txt --tic_cv_threshold 40 --cv_percent_threshold 40
```

ID-Based + Groupwise Comparison Command:

```python
python  main.py --outdirectory ./output --reportname ExampleDataset-QCReport
               --protein_level ./protein_level.txt --precursor_level ./precursor_level.txt --grouping_file ./grouping_file.txt
               --protein_threshold 200 --precursor_threshold 3000 --enzyme 'trypsin'
               --miscleavage_threshold 80 --tic_cv_threshold 40 --cv_percent_threshold 40
               --data_percent_threshold 70 --irtlabel 'iRT'--coverage_threshold 80
```

## Dataset

Dataset description + Groups

## Outputs

After the command is executed, the corresponding excel reports and html report are generated in the given output directory.
Results included here: (link to results)

For the ID-Free + ID-Based + Groupwise Comparison Command, the following outputs will be generated:

1. ExampleDataset-QCReport.html (both ID-Free and ID-Based tabs will be populated)
2. ExampleDataset-QCReport_ID-Free_QC_Report.xlsx
3. ExampleDataset-QCReport_ProteinLevel_QC_Report.xlsx
4. ExampleDataset-QCReport_PrecursorLevel_QC_Report.xlsx
5. ExampleDataset-QCReport_QC_Status_Report.xlsx (Combined overall QC PASS/FAIL status is provided from ID-Free + ID-Based metrics)

For the ID-Free + Groupwise Comparison Command, the following outputs will be generated:

1. ExampleDataset-QCReport.html (only ID-Free tab will be populated)
2. ExampleDataset-QCReport_ID-Free_QC_Report.xlsx
5. ExampleDataset-QCReport_QC_Status_Report.xlsx (Overall QC PASS/FAIL status is provided based on ID-Free metrics only)

For the ID-Based + Groupwise Comparison Command, the following outputs will be generated:

1. ExampleDataset-QCReport.html (only ID-Based tab will be populated)
3. ExampleDataset-QCReport_ProteinLevel_QC_Report.xlsx
4. ExampleDataset-QCReport_PrecursorLevel_QC_Report.xlsx
5. ExampleDataset-QCReport_QC_Status_Report.xlsx (Overall QC PASS/FAIL status is provided based on ID-Based metrics only)

### Detailed Explanation for ID-Free + ID-Based + Groupwise Comparison QC Analysis

### HTML Report

The report is divided into 2 tabs:

1. ID-Free - When a mzML directory is provided as input, ID-Free Tab is populated with plots from ID-Free parameters extracted from your mzML files.
2. ID-Based - When a protein, peptide or precursor level intensity file is provided, ID-based Tab will be populated with ID-based parameters extracted from your quantified input.

In the case of the example dataset, both ID-Free and ID-Based tabs will be populated since we have provided the mzML directory, protein level, precursor level and grouping file parameters.

#### ID-Free Tab

##### Total Ion Current

MS1 and MS2 Total Ion Current Values is extracted from spectra within the given <number> mzML files.

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

### Excel Reports
