# Using QCeltis to perform Quality Control using an Example Dataset

## Table of Contents
1. [Installation](#installation)
2. [Parameters and Commands](#parameters-and-commands)
4. [Dataset](#dataset)
6. [Outputs](#outputs-and-explanations)
7. [Results](#results-from-qc-analysis)

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

Click for [Parameters](https://github.com/vegesnam/QCPackage/blob/main/README.md#parameters)
Click for [Input File Descriptions](https://github.com/csmc-vaneykjlab/QCeltis/tree/main?tab=readme-ov-file#input-file-description)

To begin with analysis, you can download the mzML files, protein and precursor level files from [here](https://github.com/csmc-vaneykjlab/QCeltis/tree/main/example-dataset). You can provide these as input parameters and execute the following command via command-line: 

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

(Manasa: Dataset description + Batches : DBS DR) 
(How the data was acquired, number of samples + batches, where is it available)

## Inputs 

(Manasa: mzML dir - number of files, protein, precursor and grouping file)

## Outputs

After the command is executed, the corresponding excel reports and html report are generated in the given output directory.
Results included here: (Manasa: link to results)

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

## Results from QC Analysis

### HTML Report

The report is divided into 2 tabs:

1. ID-Free - When a mzML directory is provided as input, ID-Free Tab is populated with plots from ID-Free parameters extracted from your mzML files.
2. ID-Based - When a protein, peptide or precursor level intensity file is provided, ID-based Tab will be populated with ID-based parameters extracted from your quantified input.

In the case of the example dataset, both ID-Free and ID-Based tabs will be populated since we have provided the mzML directory, protein level, precursor level and grouping file parameters.

#### ID-Free Tab

##### Total Ion Current

MS1 and MS2 Total Ion Current Values are extracted from spectra within the given 47 mzML files.

![TIC_Plot](https://github.com/vegesnam/QCPackage/assets/87665957/bafb01fb-00cd-4a36-8e9e-ab7d33ad4ae5)

If extreme values are found, they are labelled as outliers.

![MS1 TIC Outliers](https://github.com/vegesnam/QCPackage/assets/87665957/c3904df1-3c67-4589-b49f-068ba174ccda)

![MS2 TIC Outliers](https://github.com/vegesnam/QCPackage/assets/87665957/ceb42e10-d10c-48cc-a8f9-91baf64dc8ed)

Here, 3 outliers were found. The following files have been detected as outliers: 
* 20211018_Seroconversion_DBS_Plate1_DR4_centroid.mzML
* 20211018_Seroconversion_DBS_Plate6_DR4_centroid.mzML
* 20211018_Seroconversion_DBS_Plate8_DR2_centroid.mzML

(Manasa: Reason for the 3 outliers)  

Groupwise comparison using CV% is performed with MS1 and MS2 TIC values:

![MS1 TIC CV%](https://github.com/vegesnam/QCPackage/assets/87665957/6e6ffe91-a1f1-4e99-816c-4494bb342e6c)

![MS2 TIC CV%](https://github.com/vegesnam/QCPackage/assets/87665957/e6cd281f-79af-4e4a-b69d-878d69772c80)

(Manasa: All groups are under given tic cv threshold, so no issues) 

##### Spectral Ratio
(Manasa: fix this)
The spectral ratio indicates the count of MS2 Spectra relative to MS1 Spectra in each file. The mass spectrometer consistently tallies both MS1 and MS2 Counts; thus, if there are discrepancies in the numbers, it's not a sample problem but rather an issue with the MS acquisition. If the ratio isn't stable, evaluate the instrument's performance on samples identified as anomalies.

![Spectral Ratio](https://github.com/vegesnam/QCPackage/assets/87665957/5184c1ef-6d22-4f8d-82bd-120b44ead9eb)

Pretty consistent since the values all lay around 50 

##### Max Basepeak Intensity

The base peak intensity is the recorded intensity of the most intense peak from each spectrum in the mzML file. The Max Base Peak Intensity represents the highest recorded base peak intensity in each mzML file.

![Max BP Plot](https://github.com/vegesnam/QCPackage/assets/87665957/c9018489-3f43-4a6e-a113-8dfa45d90293)

Max Base Peak Intensity is expected to be consistent across QC samples. Any outliers detected are highlighted in yellow. Outliers detected could point to issues with sample pickup or samples being dried out. Usually correlates with ID-free TIC.

![Max BP Outlier](https://github.com/vegesnam/QCPackage/assets/87665957/986eb8ec-1df0-42f8-8bba-17f3fdca9733)

Here, 1 outlier was found. The sample "20211018_Seroconversion_DBS_Plate1_DR3_centroid.mzML" has an extremely high max base peak intensity compared to the other samples. 

Reason

#### ID-Based plots

##### Quantification Plots

Displays the number of proteins and peptides or precursors. When the threshold is provided, any samples not meeting the threshold is flagged as a ‘FAIL’, this could indicate an issue with the sample prep or digestion protocol.

In this case, even though the threshold is met, the last 3 groups have significantly less proteins and precursors. Further investigation showed a high percentage of missingness in the samples and indicated an issue with sample prep.
So these plots can be a good way to visualize any varying patterns within samples

Protein Quant 

![Protein Quant](https://github.com/vegesnam/QCPackage/assets/87665957/d1605eeb-6bd8-4942-b7a6-792caf262952)

Precursor Quant 

![image](https://github.com/vegesnam/QCPackage/assets/87665957/7d8128c9-f67f-4b6d-ac09-6b2148eb93d2)


##### CV Plots

CVs are calculated with intensity values from samples within each group [ Please note: groupwise input must be given ]. Data and CV Thershold: simple explanation.

If inconsistent ( data threshold failure ), then the samples were poorly digested (unless TIC suggests otherwise). If TIC is high but the CV isn't, then digestion issue or if TIC is also low, then the sample didn't get picked up properly.
[ write better explanations for this ]

![Cumulative CV](https://github.com/vegesnam/QCPackage/assets/87665957/e792951d-97c9-45e4-98ca-5086fe68aa9f)

![Proteins under CV%](https://github.com/vegesnam/QCPackage/assets/87665957/03a29fa8-d061-4377-9fad-72d15efafce3)

![Precursors under CV%](https://github.com/vegesnam/QCPackage/assets/87665957/0bb211ae-fc08-4858-b9ee-e3af3e6dcae7)

Explain reasons 

##### PCA 

Protein PCA: 

![Protein PCA](https://github.com/vegesnam/QCPackage/assets/87665957/875b8281-6152-4d85-b853-0f81ff8e795c)

Precursor PCA: 

![image](https://github.com/vegesnam/QCPackage/assets/87665957/998f01de-af69-49f1-8bf6-544805dd8eec)

Clustering/Grouping and what was done. 

##### Common TIC

Common Peptide/Precursor TIC is the summed intensity of all the common peptides or precursors found in all samples. If both peptide and precursor files are provided, only common peptide will be calculated.

![CommonTIC](https://github.com/vegesnam/QCPackage/assets/87665957/aff8cec3-bbde-41c2-9ecd-0ffe6875d24f)

If TIC is consistent with other samples here but the TIC is labelled as an outlier in ID-free, then it could be a contamination issue (if TIC is high in id-free)

![CommonTICCV](https://github.com/vegesnam/QCPackage/assets/87665957/2b02ad88-0922-40b8-9632-1a4e3fc1997a)

CVs are calculated with intensity values from samples within each group.

##### Miscleavage Plot

![image](https://github.com/vegesnam/QCPackage/assets/87665957/ad7f6047-3754-4590-afba-1d27d5835bb8)

Having a high number of peptides with 0 miscleavages indicates a good digestion protocol. If the samples don’t have enough 0 missed cleaved peptides, (doesn’t meet threshold), then this indicates a digestion issue.

Here, there is a difference in numbers and this is due to the missingness in the last 3 plates (also indicated in the protein and peptide quant) – further indicates that it is a digestion/sampleprep issue.

### Excel Reports

(Manasa: what each excel file contains - sheet names and what the sheets contain) 
