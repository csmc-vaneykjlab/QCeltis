# QCPackage
Initial Code Source for the QC Package

## Table of Contents
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Usage](#Usage)
4. [Description of Parameters](#description-of-parameters)
5. [Input File Description](#input-file-description)
6. [Contributions](#contributions)
7. [Example Output](#example-output)

## Introduction
Quality Control for Proteomics (QCP) is a Python package designed to...

## Installation
QCP can be installed from the command line using `git` and `pip`. 

First, clone the repository from GitHub: git clone https://github.com/vegesnam/QCPackage.git

Then, navigate to the downloaded directory and install the package using `pip`:

~~~bash
cd qcp
pip install -r requirements.txt
~~~

Alternatively, the package can be installed directly from PyPI:

`pip install -i < pypi instance will go here >`

## Usage 

```python
python  main.py [-h] -o OUTDIRECTORY -r REPORTNAME [-m MZML_DIRECTORY]
               [-t1 MS1_TIC_THRESHOLD] [-t2 MS2_TIC_THRESHOLD]
               [-s1 MS1_SPECTRA_THRESHOLD] [-s2 MS2_SPECTRA_THRESHOLD]
               [-bp MAX_BASEPEAK_INTENSITY] [-pt PROTEIN_LEVEL]
               [-pep PEPTIDE_LEVEL] [-pre PRECURSOR_LEVEL] [-peprt PEPTIDE_RT]
               [-prert PRECURSOR_RT] [-g GROUPING_FILE] [-peplt PEPTIDE_LIST]
               [-x PROTEIN_THRESHOLD] [-y PEPTIDE_THRESHOLD]
               [-z PRECURSOR_THRESHOLD] [-e ENZYME] [-c MISCLEAVAGE_THRESHOLD]
               [-t TIC_CV_THRESHOLD] [-s CV_PERCENT_THRESHOLD]
               [-d DATA_PERCENT_THRESHOLD] [-irt IRTLABEL]
               [-v COVERAGE_THRESHOLD]
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
| --peptide_rt             | -peprt     | Path to peptide retention time file              | None          |
| --precursor_rt           | -prert     | Path to precursor retention time file            | None          |
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

## outputs and explanantions

[ Place holder ]

## Cite 

## Support

## release notes 



