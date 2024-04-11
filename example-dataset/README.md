# QCeltis analysis using an example dataset

## Table of Contents
1. [Installation](#installation)
2. [Dataset](#dataset)
3. [Inputs](#inputs)
4. [Parameters and Commands](#parameters-and-commands)
5. [Outputs](#outputs-and-explanations)
6. [Results](#results-from-qc-analysis)

## Installation
QCeltis can be installed from the command line using `git` and `pip`.

First, clone the repository from GitHub:

git clone [https://github.com/vegesnam/QCPackage.git](https://github.com/csmc-vaneykjlab/QCeltis.git)

Then, navigate to the downloaded directory and install the package using `pip`:

~~~bash
cd qcp
pip install -r requirements.txt
~~~

Alternatively, the package can be installed directly from PyPI:

`pip install -i < pypi instance will go here >`

## Dataset

The cohort consisted of dried blood samples spotted on Whatman cards, samples were prepped across 8 digestion plates and acquired across 15 MS batches. Along with cohort dried blood samples, an aliquot of pooled, mixed-gender dried blood samples were used as the quality control (QC) samples while performing the mass-spectrometry-based assays. Around 6 QC samples, known as digestion replicates (DRs) were embedded between cohort samples across each digestion plate. Technical Replicate (TR) quality control samples were digested once and ran before every MS batch. In this example dataset, we are going to perform analysis on the digestion replicate samples from the 8 digestion plates to monitor data quality and reproducibility. 

## Inputs 

To begin with analysis, you can download the mzML files, protein and precursor level files from [here](https://github.com/csmc-vaneykjlab/QCeltis/tree/main/example-dataset). 

- mzML files can be downloaded from <PXID>
- Protein Intensity File containing Protein IDs, 47 sample columns with protein intensity values
- Precursor Intensity File containing Protein IDs, Peptide IDs, Precursor IDs, 47 sample columns with precursor intensity values
- Grouping File containing sample filename and the respective plate information

## Parameters and Commands

A complete list of input file descriptions can be found [here](https://github.com/csmc-vaneykjlab/QCeltis/tree/main?tab=readme-ov-file#input-file-description), along with the [parameters](https://github.com/vegesnam/QCPackage/blob/main/README.md#parameters) available. 

You can provide the following input parameters and execute the following commands via command-line for this example dataset: 

ID-Free + ID-Based + Groupwise Comparison Command:

```python
python  main.py --outdirectory ./output --reportname ExampleDataset-QCReport --mzml_directory ./mzML_files
               --protein_level ./protein_level.txt --precursor_level ./precursor_level.txt --grouping_file ./grouping_file.txt
               --protein_threshold 200 --precursor_threshold 3000 --enzyme 'trypsin'
               --miscleavage_threshold 60 --tic_cv_threshold 40 --cv_percent_threshold 40
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
               --miscleavage_threshold 60 --tic_cv_threshold 40 --cv_percent_threshold 40
               --data_percent_threshold 70 --irtlabel 'iRT'--coverage_threshold 80
```

## Outputs

After the command is executed, the corresponding excel reports and html report are generated in the given output directory.
Results included [here](https://github.com/csmc-vaneykjlab/QCeltis/tree/main/example-dataset/results): 

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

### ID-Free Tab

#### Total Ion Current

MS1 and MS2 Total Ion Current values are extracted from spectra within the given 47 mzML files.

<ins>TIC Plot</ins>:

![TIC Plot](https://github.com/csmc-vaneykjlab/QCeltis/assets/87665957/a9d88777-eb00-40c0-b659-b2e0810b3834)

If extreme values are found, they are labelled as outliers.

![MS1 TIC Outliers](https://github.com/vegesnam/QCPackage/assets/87665957/c3904df1-3c67-4589-b49f-068ba174ccda)

![MS2 TIC Outliers](https://github.com/vegesnam/QCPackage/assets/87665957/ceb42e10-d10c-48cc-a8f9-91baf64dc8ed)

Here, 3 outliers were found. The following files have been detected as outliers: 
* 20211018_Seroconversion_DBS_Plate1_DR4_centroid.mzML
* 20211018_Seroconversion_DBS_Plate6_DR4_centroid.mzML
* 20211018_Seroconversion_DBS_Plate8_DR2_centroid.mzML

These 3 outlier samples were rerun. 

<ins>TIC CV% Plot</ins>:

Groupwise comparison using CV% is performed with MS1 and MS2 TIC values:

![MS1 TIC CV%](https://github.com/vegesnam/QCPackage/assets/87665957/6e6ffe91-a1f1-4e99-816c-4494bb342e6c)

![MS2 TIC CV%](https://github.com/vegesnam/QCPackage/assets/87665957/e6cd281f-79af-4e4a-b69d-878d69772c80)

The intra-plate CV% is lower than the set threshold, which indicates consistent TIC values within the DRs of each plate. All the plates have passed.  

#### Spectral Ratio

MS1 and MS2 spectral counts are extracted from the input mzML files. The MS2-to-MS1 spectral ratio indicates the count of MS2 Spectra relative to MS1 Spectra in each file. The mass spectrometer consistently tallies both MS1 and MS2 Counts; thus, if there are discrepancies in the numbers, it's not a sample preparation or digestion problem but rather an issue with the MS acquisition. If the ratio isn't stable, evaluate the instrument's performance on samples identified as anomalies.

![Spectral Ratio](https://github.com/vegesnam/QCPackage/assets/87665957/5184c1ef-6d22-4f8d-82bd-120b44ead9eb)

In this case, the spectral ratio is consistent across the samples, with values that all lay around 50.  

#### Max Basepeak Intensity

The base peak intensity is the recorded intensity of the most intense peak from each spectrum in the mzML file. The Max Base Peak Intensity represents the highest recorded base peak intensity in each mzML file.

![Max BP Plot](https://github.com/vegesnam/QCPackage/assets/87665957/c9018489-3f43-4a6e-a113-8dfa45d90293)

Max Base Peak Intensity is expected to be consistent across replicate QC samples. Any outliers detected are highlighted in yellow. Outliers detected could point to issues with sample pickup or samples being dried out. 

![Max BP Outlier](https://github.com/vegesnam/QCPackage/assets/87665957/986eb8ec-1df0-42f8-8bba-17f3fdca9733)

Here, 1 outlier was found. The sample "20211018_Seroconversion_DBS_Plate1_DR3_centroid.mzML" has an extremely high max base peak intensity compared to the other samples. 

### ID-Based plots

#### Quantification Plots

Number of proteins, peptides and precursors found in each sample are displayed using a bar plot. When the threshold is provided, any samples not meeting the threshold is flagged as a ‘FAIL’, this could indicate an issue with the sample preparation or digestion protocol. These plots are a good way to visualize any varying patterns within samples and detect sample preparation or digestion issues. 

<ins>Protein Quant</ins>:

![Protein Quant](https://github.com/vegesnam/QCPackage/assets/87665957/d1605eeb-6bd8-4942-b7a6-792caf262952)

<ins>Precursor Quant</ins>:

![image](https://github.com/vegesnam/QCPackage/assets/87665957/7d8128c9-f67f-4b6d-ac09-6b2148eb93d2)


#### CV Plots

Using the protein, peptide or precursor intensity values provided, coefficient of variation (CV) percent values are calculated across the entire dataset and across samples within each group, in this case, within each plate.  

<ins>Cumulative CV% Plot</ins>:

Cumulative frequency percentage of calculated CV% of intensity values across all samples is plotted. This plot reveals the degree of variability across the dataset, higher CVs indicate greater variation that could be stemming from sample preparation, data acquisition or instrument performance. Lower CVs indicate higher reproducibility of protein, peptide or precursor intensities across replicate samples.

![Cumulative CV](https://github.com/vegesnam/QCPackage/assets/87665957/e792951d-97c9-45e4-98ca-5086fe68aa9f)

<ins>Percentage of Proteins under CV%</ins>: 

Percent of proteins under the provided CV threshold is plotted. A data threshold, which represents the minimum percent of proteins under CV threshold, is applied and any groups not meeting the threshold are marked as "FAIL"

![Proteins under CV%](https://github.com/vegesnam/QCPackage/assets/87665957/03a29fa8-d061-4377-9fad-72d15efafce3)

<ins>Percentage of Precursors under CV%</ins>:

Percent of precursors under the provided CV threshold is plotted. A data threshold, which represents the minimum percent of precursors under CV threshold, is applied and any groups not meeting the threshold are marked as "FAIL"

![Precursors under CV%](https://github.com/vegesnam/QCPackage/assets/87665957/0bb211ae-fc08-4858-b9ee-e3af3e6dcae7)

Here, in this example, all the groups have passed the intensity CV and data thresholds. 

#### PCA 

Protein, peptide or precursor intensities are used to perform principal component analysis across provided groups. Any clustering observed across replicate samples can be an indication of batch effects. 

<ins>Protein PCA</ins>: 

![Protein PCA](https://github.com/vegesnam/QCPackage/assets/87665957/875b8281-6152-4d85-b853-0f81ff8e795c)

<ins>Precursor PCA</ins>: 

![image](https://github.com/vegesnam/QCPackage/assets/87665957/998f01de-af69-49f1-8bf6-544805dd8eec)

#### Common TIC

Common Peptide/Precursor TIC is the summed intensity of all the common peptides or precursors found in all samples. If both peptide and precursor files are provided, only common peptide will be calculated.

<ins>Common Precursor TIC</ins>: 

![CommonTIC](https://github.com/vegesnam/QCPackage/assets/87665957/aff8cec3-bbde-41c2-9ecd-0ffe6875d24f)

<ins>Common Precursor TIC CV%</ins>:

CVs are calculated with intensity values from common precursors across the samples within each group. All groups within the TIC CV threshold are considered to have passed.

![CommonTICCV](https://github.com/vegesnam/QCPackage/assets/87665957/2b02ad88-0922-40b8-9632-1a4e3fc1997a)


#### Miscleavage Plot

Number of peptides with 0 miscleavages are plotted from each sample. Having a high number of peptides with 0 miscleavages indicates a good digestion protocol. If the samples don’t have enough 0 missed cleaved peptides, (doesn’t meet threshold), then this indicates a digestion issue.

![No Miscleavages](https://github.com/csmc-vaneykjlab/QCeltis/assets/87665957/412eadc6-e056-486d-8bc5-4d1381d0dc20)


### Excel Reports

1. ExampleDataset-QCReport_ID-Free_QC_Report.xlsx: contains id-free metric summary and groupwise TIC CV% calculations. 
2. ExampleDataset-QCReport_ProteinLevel_QC_Report.xlsx: contains protein-level id-based metrics summary, such as quant summary and groupwise CV% calculations.
3. ExampleDataset-QCReport_PrecursorLevel_QC_Report.xlsx: contains precursor-level id-based metrics summary, such as quant summary, CV% calculations, common precursor TIC and miscleavage summary.
4. ExampleDataset-QCReport_QC_Status_Report.xlsx: contains combined overall QC PASS/FAIL status from id-free and id-based metrics for each sample and group 
