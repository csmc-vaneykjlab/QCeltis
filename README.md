# QCPackage
Initial Code Source for the QC Package

## Table of Contents
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Description of Parameters](#description-of-parameters)
4. [Example Run Command](#example-run-command)
5. [Contributions](#contributions)
6. [Example Output](#example-output)

## Introduction
Quality Control for Proteomics (QCP) is a Python package designed to...

## Installation
QCP can be installed from the command line using `git` and `pip`. 

First, clone the repository from GitHub: git clone https://github.com/vegesnam/QCPackage.git

Then, navigate to the downloaded directory and install the package using `pip`:

cd qcp
pip install -r requirements.txt
pip install .

Alternatively, the package can be installed directly from PyPI:
pip install -i < pypi instance will go here > 

## Usage 

## Parameters 

--outdirectory or -o: Output directory path

--reportname or -r: Report name for HTML and Excel reports

--mzml_directory or -m: Path to the directory where mzML files are present

--ms1_tic_threshold or -t1: MS1 TIC threshold

--ms2_tic_threshold or -t2: MS2 TIC threshold

--ms1_spectra_threshold or -s1: MS1 spectra threshold

--ms2_spectra_threshold or -s2: MS2 spectra threshold

--max_basepeak_intensity or -bp: Maximum base peak intensity threshold

--protein_level or -pt: Path to protein intensity file

--peptide_level or -pep: Path to peptide intensity file

--precursor_level or -pre: Path to precursor intensity file

--peptide_rt or -peprt: Path to peptide retention time file

--precursor_rt or -prert: Path to precursor retention time file

--grouping_file or -g: Path to grouping file

--peptide_list or -peplt: Path to file containing list of peptides to monitor intensity and RT distribution across samples

--protein_threshold or -x: Protein threshold for each sample

--peptide_threshold or -y: Peptide threshold for each sample

--precursor_threshold or -z: Precursor threshold for each sample

--enzyme or -e: User input enzyme

--miscleavage_threshold or -c: Missed cleavage threshold for each sample

--tic_cv_threshold or -t: TIC CV threshold for groupwise QC status

--cv_percent_threshold or -s: Intensity CV threshold

--data_percent_threshold or -d: Data threshold for intensity CV

--irtlabel or -irt: Label for iRT peptides present in your peptide intensity file

--coverage_threshold or -v: Intensity or retention time coverage % threshold in each sample

## Example run command

python main.py \ [ to be changed with pip execution command ]
--outdirectory /path/to/output/directory \
--reportname ReportName \
[Optional] --mzml_directory /path/to/mzml/directory \
[Optional] --ms1_tic_threshold 1000 \
[Optional] --ms2_tic_threshold 2000 \
[Optional] --ms1_spectra_threshold 3000 \
[Optional] --ms2_spectra_threshold 4000 \
[Optional] --max_basepeak_intensity 5000 \
[Optional] --protein_level /path/to/protein/intensity/file \
[Optional] --peptide_level /path/to/peptide/intensity/file \
[Optional] --precursor_level /path/to/precursor/intensity/file \
[Optional] --peptide_rt /path/to/peptide/retention/time/file \
[Optional] --precursor_rt /path/to/precursor/retention/time/file \
[Optional] --grouping_file /path/to/grouping/file \
[Optional] --peptide_list /path/to/peptide/list/file \
[Optional] --protein_threshold 100 \
[Optional] --peptide_threshold 200 \
[Optional] --precursor_threshold 300 \
[Optional] --enzyme enzyme_name \
[Optional] --miscleavage_threshold 10 \
[Optional] --tic_cv_threshold 15 \
[Optional] --cv_percent_threshold 20 \
[Optional] --data_percent_threshold 25 \
[Optional] --irtlabel iRT_Label \
[Optional] --coverage_threshold 30

## outputs and explanantions

[ Place holder ]

## Cite 

## Support

## release notes 



