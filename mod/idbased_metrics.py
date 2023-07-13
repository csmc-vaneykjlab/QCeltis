"""
Functions to calculate QC metrics from protein, peptide and precursor level data
"""
import pandas as pd
import sys
import os
import numpy as np
import pymzml
import statistics
import pandas as pd
import os
import re
import threading
import logging
import xlsxwriter
import collections
from collections import defaultdict
import plotly
import plotly.express as px
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

from mod.general_functions import cv, groupname, quant_status, check_threshold, transpose_DF, perc_qc, cv_status

enzyme_info =  {'asp-n':{'terminus' : 'N' , 'cleave' : ['D'], 'exceptions' : []},
                'lys-c' : {'terminus' : 'C' , 'cleave' : ['K'], 'exceptions' : ['KP']},
                'lys-n':{'terminus' : 'N' , 'cleave' : ['K'], 'exceptions' : []},
                'thermolysin':{'terminus' : 'N' , 'cleave' : ['A','F','I','L','M','V'],
                                'exceptions' : ['DA','DF','DI','DL','DM','DV','EA','EF','EI','EL','EM','EV']},
                'proteinasek':{'terminus' : 'C' , 'cleave' : ['A','F','Y','W','L','I','V'], 'exceptions' : []},
                'mfh':{'terminus' : 'C' , 'cleave' : ['D'], 'exceptions' : []},
                'gluc':{'terminus' : 'C','cleave':['D','E'],'exceptions':['DP','EP',"EE","DE"]},
                'glucbicarb':{'terminus' : 'C','cleave':['E'],'exceptions':['EP','EE']},
                'chymotrypsinflmwy':{'terminus' : 'C' , 'cleave' : ['F','L','M','W','Y'], 'exceptions' : ['FP','LP','MP','WP','YP','PY']},
                'trypsin' : {'terminus' : 'C' , 'cleave' : ['K','R'], 'exceptions' : ['KP','RP']},
                'trypsinp':{'terminus' : 'C' , 'cleave' : ['K','R'], 'exceptions' : []},
                'cnbr':{'terminus' : 'C' , 'cleave' : ['M'], 'exceptions' : []},
                'arg-c':{'terminus' : 'C' , 'cleave' : ['R'], 'exceptions' : ['RP']}}

combined_pat = r'|'.join(('\[.*?\]', '\(.*?\)','\{.*?\}'))

#------------------------------------------------------------------------- FUNCTIONS ---------------------------------------------------------------------------------

def calc_miscleavage(each_peptide, enzyme):
    ''' Calculate miscleavages for a peptide sequence, based on enzyme selected by user. Considerations to exclude exceptions to the cleavage sites are made. If no end cleavage site found, add 1000 to the original MC_score as punishment'''
    cleavage_score = 0
    FOUND = False
    cleave_site = set(enzyme_info[enzyme]['cleave'])
    exceptions = set(enzyme_info[enzyme]['exceptions'])

    peptide = re.sub(combined_pat,'',each_peptide)

    # adjust and calculate miscleavage score
    if enzyme_info[enzyme.lower()]['terminus']=="C":
        for i in range (len(peptide)):
            if peptide[i] in cleave_site:
                FOUND = True
                if peptide[i:(i+2)] in exceptions:
                    pass
                else:
                    cleavage_score += 1

        if peptide[-1] in cleave_site:
            cleavage_score -= 1
        else:
            cleavage_score += 1000

    elif enzyme_info[enzyme.lower()]['terminus']=="N":
        for i in range (len(peptide)):
            if peptide[i] in cleave_site:
                FOUND = True

                if peptide[(i-1):(i+1)] in exceptions:
                    pass
                else:
                    cleavage_score += 1

        if peptide[0] in cleave_site:
            cleavage_score -= 1
        else:
            cleavage_score += 1000
    else:
        pass

    if FOUND == False:
        return (1000)

    return(cleavage_score)

def get_quant(df, filenames, threshold, level, groupwise_comparison, groups):

    #getting quant df
    quant = pd.DataFrame(df[filenames].count())
    quant = quant.reset_index()
    quant.columns = ['Filename', f'{level} Number']

    if threshold:
        threshold_col = f'{level} Threshold = {threshold}'
        quant[threshold_col] = quant[f'{level} Number'].apply(check_threshold, args=[threshold,])

        if groupwise_comparison:
            #get group based PASS or FAIL
            group_status_dict = {}
            group_quant = quant[['Filename', threshold_col]]
            group_quant['Group'] = group_quant['Filename'].apply(groupname, args=[groups,])
            group_quant = group_quant[['Group',threshold_col]]

            for group in list(set(group_quant['Group'].tolist())):
                subset_quant = group_quant[group_quant['Group'] == group]
                group_status = list(set(subset_quant[threshold_col].tolist()))

                if group_status == ['PASS']:
                    group_status_dict[group] = "PASS"
                else:
                    group_status_dict[group] = "FAIL"


            group_status_df = pd.DataFrame(list(group_status_dict.items()), columns=['Group',f'{level} Threshold QC Status'])

        else:
            group_status_df = ""
    else:
        group_status_df = ""

    return (quant, group_status_df)

def intensity_cvs(df, intensity_cv_threshold, data_percent_threshold, filenames, level, groupwise_comparison, groups):

    overall_cv_df = df
    overall_cv_df['Overall Average'] = overall_cv_df[filenames].mean(axis=1)
    overall_cv_df['Overall Standard Deviation'] = overall_cv_df[filenames].std(axis=1)
    overall_cv_df['Overall CV %'] = (overall_cv_df['Overall Standard Deviation']/overall_cv_df['Overall Average'])*100

    cv_bins = [10,20,30,40,50,60,70,80,90,100]

    cumfeq_dict = {}
    cv_level = overall_cv_df
    cv_level.dropna(subset=['Overall CV %'], inplace=True)

    for num in cv_bins:
        feature_count = len(cv_level[cv_level['Overall CV %'] <= num]['Overall CV %'].tolist())
        cumfeq_dict[num] = {f"{level} Cumulative Frequency %" :(feature_count/len(cv_level['Overall CV %'].tolist()))*100}

    cv_sum = pd.DataFrame.from_dict(cumfeq_dict, orient="index")
    cv_sum.reset_index(drop=False, inplace=True)
    cv_sum.rename({"index":"CV%"}, axis=1, inplace=True)

    if groupwise_comparison:
        if intensity_cv_threshold and data_percent_threshold:

            grouped_cv_dict = {}

            for group in groups:
                df[f'{group}-Average'] = df[groups[group]].mean(axis=1)
                df[f'{group}-Standard Deviation'] = df[groups[group]].std(axis=1)
                df[f'{group}-CV%'] = (df[f'{group}-Standard Deviation']/df[f'{group}-Average'])*100

                group_subset_cols = [level] + groups[group] + [f'{group}-CV%']
                group_subset = df[group_subset_cols].dropna(subset=groups[group], how='all')
                total = len(list(group_subset[level].tolist()))
                total_under_cv = len(list(set(group_subset[group_subset[f'{group}-CV%'] <= intensity_cv_threshold][level].tolist())))
                total_perc_under_cv = round((total_under_cv/total)*100, 2)

                if total_perc_under_cv >= data_percent_threshold:
                    group_cv_status = "PASS"
                else:
                    group_cv_status = "FAIL"

                grouped_cv_dict[group] = {f"{level} Number": total,
                                         f"{level} CV% < {intensity_cv_threshold}":total_under_cv,
                                         f"{data_percent_threshold}% {level}s <= {intensity_cv_threshold}% CV": group_cv_status,}

                grouped_df = pd.DataFrame.from_dict(grouped_cv_dict, orient="index")
                grouped_df.reset_index(drop=False, inplace=True)
                grouped_df.rename(columns = {"index":"Group"}, inplace=True)

    else:
        grouped_df = ""

    return (overall_cv_df, cv_sum, grouped_df)

def miscleavage(pep_level, enzyme, miscleavage_threshold, filenames, groups, groupwise_comparison):

    peptide_level = pep_level[['Peptide'] + filenames]
    peptide_level.fillna(0, inplace=True)
    peptide_level.drop_duplicates(subset=['Peptide'], inplace=True)
    peptide_level.index = peptide_level['Peptide']
    peptide_level.drop(['Peptide'], axis=1, inplace=True)

    digestion_info = peptide_level.to_dict()
    filesNames = sorted(digestion_info.keys())

    digestionDF = pd.DataFrame()
    digestionDF["Digestion"] = ["No end cleavage site found","0 missed cleavage","1 missed cleavage","2 missed cleavage",
                           "3 missed cleavage","4 missed cleavage",">=5 missed cleavage","Total Peptides"]

    cleavageScore_dict = {}
    for file in filesNames:
        mcScoreCount_list = []
        for pepSeqs in digestion_info[file]:
            if not int(digestion_info[file][pepSeqs]) == 0:
                mcScore = calc_miscleavage(pepSeqs,enzyme)
                cleavageScore_dict[pepSeqs] = mcScore
                mcScoreCount_list.append(mcScore)

        mcScore_Counter = collections.Counter(mcScoreCount_list)
        ms_list = [0]*(len(digestionDF["Digestion"])-1)

        for mcScore in mcScore_Counter:
            if mcScore >= 1000:
                ms_list[0]+= mcScore_Counter[mcScore]
            elif mcScore >= 0 and mcScore <=4:
                ms_list[mcScore+1] = mcScore_Counter[mcScore]
            elif mcScore >=5 and mcScore < 1000:
                ms_list[-1] += mcScore_Counter[mcScore]
            else:
                print ("Error in digestion")

        fileSum = sum(ms_list)

        ms_list.append(fileSum)
        digestionDF[file] = ms_list

    dig_df = transpose_DF(digestionDF)
    headers = dig_df.iloc[0].values
    dig_df.columns = headers
    dig_df.drop(index=0, axis=0, inplace=True)
    dig_df['0 missed cleavage percentage'] = (dig_df['0 missed cleavage']/dig_df['Total Peptides'])*100

    if miscleavage_threshold:
        dig_df['0 missed cleavage QC Status'] = dig_df['0 missed cleavage percentage'].apply(perc_qc, args=[int(miscleavage_threshold),])

        if groupwise_comparison and miscleavage_threshold:
            group_dig_df = dig_df[['Digestion','0 missed cleavage QC Status']]
            group_dig_df['Group'] = group_dig_df['Digestion'].apply(groupname, args=[groups,])
            group_dig_df = group_dig_df[['Group','0 missed cleavage QC Status']]

            group_level = {}
            for group in list(set(group_dig_df['Group'].tolist())):
                group_subset = group_dig_df[group_dig_df['Group'] == group]
                status = list(set(group_subset['0 missed cleavage QC Status'].tolist()))
                if status == ['PASS']:
                    group_level[group] = "PASS"
                else:
                    group_level[group] = "FAIL"

            group_df = pd.DataFrame(list(group_level.items()), columns=['Group','0 Miscleaved Peptides QC Status'])
            group_df['0 Miscleaved Peptides QC Status'] = group_df['0 Miscleaved Peptides QC Status'].astype(str)

        else:
            group_df = ""

    else:
        group_df = ""

    dig_df.rename(columns={'Digestion':'Filename'}, inplace=True)

    return (dig_df, group_df)

def common_tic(df_level, level, tic_cv_threshold, filenames, groups, groupwise_comparison):

    if level == "Peptide":
        df = df_level[['Peptide'] + filenames]
    elif level == "Precursor":
        df = df_level[['Precursor'] + filenames]

    df.replace({0.0:np.nan,0:np.nan}, inplace=True)
    df.dropna(inplace=True)

    #calculating TIC
    df.loc['TIC'] = df.sum(numeric_only=True, axis=0)
    df_tic = df.T.reset_index(drop=False)[['index','TIC']]
    df_tic.rename({'index':'Filename'}, axis=1, inplace=True)
    df_tic = df_tic[df_tic['Filename'] != level]

    if groupwise_comparison and tic_cv_threshold:
        #getting group level information
        df_tic['Group'] = df_tic['Filename'].apply(groupname, args=[groups, ])

        #calculating cv across plates
        group_level = {}
        for group in list(set(df_tic['Group'].tolist())):
            group_subset = df_tic[df_tic['Group'] == group]
            tic_list = group_subset['TIC'].tolist()
            group_level[group] = round(cv(tic_list),2)

        group_df = pd.DataFrame(list(group_level.items()), columns=['Group','CV %'])
        group_df[f'Common {level} TIC QC Status'] = group_df['CV %'].apply(cv_status, args=[float(tic_cv_threshold),])

        df_tic.drop(['Group'], axis=1, inplace=True)
        group_df[f'Common {level} TIC QC Status'] = group_df[f'Common {level} TIC QC Status'].astype(str)

    else:
        group_df = ""

    return (df_tic, group_df)

def selected_peps(df_level, level, coverage_threshold, filenames, irtlabel, peptide_list, peptide_list_df):

    #iRT kit
    irt_peptides = ["LGGNEQVTR","GAGSSEPVTGLDAK","VEATFGVDESNAK","YILAGVENSK","TPVISGGPYEYR","TPVITGAPYEYR","DGLDAASYYAPVR","ADVTPADFSEWSK","GTFIIDPGGVIR","GTFIIDPAAVIR","LFLQFGAQGSPFLK"]
    irt_precursors = ["LGGNEQVTR2","GAGSSEPVTGLDAK2","VEATFGVDESNAK2","YILAGVENSK2","TPVISGGPYEYR2","TPVITGAPYEYR2","DGLDAASYYAPVR2","ADVTPADFSEWSK2","GTFIIDPGGVIR2","GTFIIDPAAVIR2","LFLQFGAQGSPFLK2"]

    if not irtlabel:
        irt_plots = False
        irt_level = ""
        logging.info("No iRT Label provided, iRT QC Analysis will not be performed")

    if irtlabel:
        irt_plots = True

        irt_present = False
        for protein in list(set(df_level['Protein'].tolist())):
            if irtlabel in protein:
                irt_present = True
                break

        if not irt_present:
            irt_plots = False
            irt_level = ""
            logging.info("No labelled iRT peptides found, no iRT QC analysis will be performed")

        if irt_present:
            if level == "Peptide":
                irt_level = df_level[df_level['Peptide'].isin(irt_peptides)]
            if level == "Precursor":
                irt_level = df_level[df_level['Precursor'].isin(irt_precursors)]

            irt_level['Coverage %'] = round((irt_level[filenames].count(axis=1)/len(filenames))*100,2)

            if coverage_threshold:
                irt_level[f'Coverage Threshold = {coverage_threshold}'] = irt_level['Coverage %'].apply(check_threshold,  args=[coverage_threshold,])

    #selected peptide list
    if peptide_list:
        peptide_list = list(set(peptide_list_df['Peptide'].tolist()))

        peptides_not_found = []
        for peptide in peptide_list:
            if peptide not in list(set(df_level['Peptide'].tolist())):
                peptides_not_found.append(peptide)

        if len(peptides_not_found) > 0:
            #logging.info(f"The following peptides provided in the peptide list were not found in the input files: {", ".join(peptides_not_found)}")
            print(f"The following peptides provided in the peptide list were not found in the input files: {', '.join(peptides_not_found)}")

        selected_pep_df = df_level[df_level['Peptide'].isin(peptide_list)]
        selected_pep_df['Coverage %'] = round((selected_pep_df[filenames].count(axis=1)/len(filenames))*100,2)

        if coverage_threshold:
            selected_pep_df[f'Coverage Threshold = {coverage_threshold}'] = selected_pep_df['Coverage %'].apply(check_threshold,  args=[coverage_threshold,])

    else:
        selected_pep_df = ""

    return (irt_level, irt_plots, selected_pep_df)

def get_sample_df(protein_level, peptide_level, precursor_level, pt_sample_df, pep_sample_df, pre_sample_df, threshold_dict):

    if protein_level and peptide_level and precursor_level and threshold_dict['Protein Threshold'] and threshold_dict['Peptide Threshold'] and threshold_dict['Precursor Threshold']:
        overall_sample_df = pd.merge(pt_sample_df, pep_sample_df, on="Filename")
        overall_sample_df = pd.merge(overall_sample_df, pre_sample_df, on="Filename")

    elif protein_level and peptide_level and threshold_dict['Protein Threshold'] and threshold_dict['Peptide Threshold']:
        overall_sample_df = pd.merge(pt_sample_df, pep_sample_df, on="Filename")

    elif protein_level and precursor_level and threshold_dict['Protein Threshold'] and threshold_dict['Precursor Threshold']:
        overall_sample_df = pd.merge(pt_sample_df, pre_sample_df, on="Filename")

    elif peptide_level and precursor_level and threshold_dict['Peptide Threshold'] and threshold_dict['Precursor Threshold']:
        overall_sample_df = pd.merge(pep_sample_df, pre_sample_df, on="Filename")

    elif protein_level and threshold_dict['Protein Threshold']:
        return pt_sample_df

    elif peptide_level and threshold_dict['Peptide Threshold']:
        return pep_sample_df

    elif precursor_level and threshold_dict['Precursor Threshold']:
        return pre_sample_df

    return overall_sample_df

def get_overall_df(protein_level, peptide_level, precursor_level, pt_group_df, pep_group_df, pre_group_df):

    if protein_level and peptide_level and precursor_level:
        overall_group_df = pd.merge(pt_group_df, pep_group_df, on="Group")
        overall_group_df = pd.merge(overall_group_df, pre_group_df, on="Group")

    elif protein_level and peptide_level:
        overall_group_df = pd.merge(pt_group_df, pep_group_df, on="Group")

    elif protein_level and precursor_level:
        overall_group_df = pd.merge(pt_group_df, pre_group_df, on="Group")

    elif peptide_level and precursor_level:
        overall_group_df = pd.merge(pep_group_df, pre_group_df, on="Group")

    elif protein_level:
        return pt_group_df

    elif peptide_level:
        return pep_group_df

    elif precursor_level:
        return pre_group_df

    return overall_group_df

#---------------------------------------------------------------------- GRAPH FUNCTIONS -----------------------------------------------------------------------------

def get_quant_plot(quant_df, threshold, level, groupwise_comparison, groups):

    #getting protein quant graph
    if groupwise_comparison:
        quant_df['Group'] = quant_df['Filename'].apply(groupname, args=[groups,])
        quant_plot = px.bar(quant_df, x='Filename', y=f'{level} Number', title=f"Number of {level}s Identified", color="Group")
    else:
        quant_plot = px.bar(quant_df, x='Filename', y=f'{level} Number', title=f"Number of {level}s Identified")

    if threshold:
        quant_plot.add_hline(y=threshold, line_dash="dot", annotation_text=f"{level} Threshold = {threshold}")

    quant_plot.update_xaxes(tickfont_size=8)
    quant_plot.update_layout(margin=dict(l=20, r=20, t=20, b=20))
    quant_plot.update_layout(title={'font': {'size': 9}})

    quant_div = plotly.io.to_html(quant_plot, include_plotlyjs=True, full_html=False)

    if level == "Protein":
        quant_report_params = {'protein_file': True,
                            'protein_quant_plot': quant_div,
                            'protein_quant_description': "Number of Proteins Identified from each sample" }
    if level == "Peptide":
        quant_report_params = {'peptide_file': True,
                            'peptide_quant_plot': quant_div,
                            'peptide_quant_description': "Number of Peptides Identified from each sample" }

    if level == "Precursor":
        quant_report_params = {'precursor_file': True,
                            'precursor_quant_plot': quant_div,
                            'precursor_quant_description': "Number of Precursors Identified from each sample" }

    return quant_report_params

def intensity_cv_graphs(cv_sum, grouped_cv, level, groupwise_comparison, cv_percent_threshold, data_percent_threshold):

    feature_col_name = f"{level} CV% < {cv_percent_threshold}"
    feature_column = f"% {level} under CV% < {cv_percent_threshold}"

    grouped_cv[feature_column] = (grouped_cv[feature_col_name]/grouped_cv[f'{level} Number'])*100

    grouped_cv_graph = px.bar(grouped_cv, x='Group', y=feature_column, title=f"Number of {level}s Identified under {cv_percent_threshold}% CV", color='Group')
    grouped_cv_graph.add_hline(y=data_percent_threshold, line_dash="dot", annotation_text=f"Data Percent Threshold = {data_percent_threshold}")
    grouped_cv_graph.update_xaxes(tickfont_size=8)
    grouped_cv_graph.update_layout(
            margin=dict(l=20, r=20, t=20, b=20),
        )
    grouped_cv_graph.update_layout(title={'font': {'size': 9}})
    grouped_cv_plot = plotly.io.to_html(grouped_cv_graph, include_plotlyjs=False, full_html=False)

    if level == "Protein":
        intensity_cv_report_params = {'protein_file': True,
                            'percentage_proteins_undercv_plot': grouped_cv_plot,
                            'percentage_proteins_undercv_description': f"Number of Proteins under {cv_percent_threshold}" }

    if level == "Peptide":
        intensity_cv_report_params = {'peptide_file': True,
                            'percentage_peptides_undercv_plot': grouped_cv_plot,
                            'percentage_peptides_undercv_description': f"Number of Peptides under {cv_percent_threshold}" }

    if level == "Precursor":
        intensity_cv_report_params = {'precursor_file': True,
                            'percentage_precursors_undercv_plot': grouped_cv_plot,
                            'percentage_precursors_undercv_description': f"Number of Peptides under {cv_percent_threshold}" }

    return intensity_cv_report_params

def pca_plot(df_level, level,filenames, groups): #call it only if groupwise is given

    #keeping only required columns and filling NAs with 0
    df = df_level[[level]+filenames]
    df.fillna(0, inplace=True)

    samples = df.drop(level, axis=1).columns.tolist()
    df.set_index(level, inplace=True)
    df = df.T
    df['Samples'] = df.index
    df.reset_index(inplace=True, drop=True)

    df['Group'] = df['Samples'].apply(groupname, args=[groups, ])
    df.sort_values(by=['Group'], inplace=True)
    samples = df.drop(['Samples', 'Group'], axis=1).columns.tolist()
    features = df.loc[:, samples]
    features = StandardScaler().fit_transform(features)

    pca3 = PCA(n_components=3)
    pC3 = pca3.fit_transform(features)
    total_var = pca3.explained_variance_ratio_.sum() * 100

    pca_fig3 = px.scatter_3d(
        pC3, x=0, y=1, z=2, color=df['Group'],
        title=f'Total Explained Variance: {total_var:.2f}%',
        labels={'0': 'PC 1', '1': 'PC 2', '2': 'PC 3'}
    )

    pca3_graph = plotly.io.to_html(pca_fig3, include_plotlyjs=False, full_html=False)

    if level == "Protein":
        pca_report_params = {'protein_pca_plot':pca3_graph,
                            'protein_pca_description': "PCA for Protein Intensities Across Groups"}

    if level == "Peptide":
        pca_report_params = {'peptide_pca_plot':pca3_graph,
                            'peptide_pca_description': "PCA for Peptide Intensities Across Groups"}

    if level == "Precursor":
        pca_report_params = {'precursor_pca_plot':pca3_graph,
                            'precursor_pca_description': "PCA for Precursor Intensities Across Groups"}

    return pca_report_params

def common_tic_plot(df_tic, group_tic, level, tic_cv_threshold, groupwise_comparison, groups):

    if groupwise_comparison:
        df_tic['Group'] = df_tic['Filename'].apply(groupname, args=[groups, ])
        tic_bar = px.bar(df_tic, x='Filename', y='TIC', title=f"Common {level} TIC", color=df_tic['Group'].tolist())
    else:
        tic_bar = px.bar(df_tic, x='Filename', y='TIC', title=f"Common {level} TIC")

    tic_bar.update_xaxes(tickfont_size=8)
    tic_bar.update_layout(
        margin=dict(l=20, r=20, t=20, b=20),
    )
    tic_bar.update_layout(title={'font': {'size': 9}})
    tic_bar_graph = plotly.io.to_html(tic_bar, include_plotlyjs=True, full_html=False)

    if level == "Peptide":
        common_tic_report_params = {'common_peptide_tic_plot': tic_bar_graph,
                                'common_peptide_tic_description': "Total Ion Current calculated from peptide intensities found in all given samples"}

    if level == "Precursor":
        common_tic_report_params = {'common_precursor_tic_plot': tic_bar_graph,
                                'common_precursor_tic_description': "Total Ion Current calculated from precursor intensities found in all given samples"}

    if groupwise_comparison:

        cv_bar = px.bar(group_tic, x='Group', y='CV %', title=f'Common {level} TIC Group CV%', color=group_tic['Group'].tolist())

        if tic_cv_threshold:
            cv_bar.add_hline(y=tic_cv_threshold, line_dash="dot", annotation_text=f"TIC CV Threshold = {tic_cv_threshold}")

        cv_bar.update_xaxes(tickfont_size=8)
        cv_bar.update_layout(
            margin=dict(l=20, r=20, t=20, b=20),
        )
        cv_bar.update_layout(title={'font': {'size': 9}})
        cv_bar_graph = plotly.io.to_html(cv_bar, include_plotlyjs=True, full_html=False)

        if level == "Peptide":
            common_tic_report_params['common_peptide_tic_group_cv'] = cv_bar_graph
            common_tic_report_params['common_peptide_tic_group_cv_description'] = "Common Peptide TIC CV% calculated across groups"

        if level == "Precursor":
            common_tic_report_params['common_precursor_tic_group_cv'] = cv_bar_graph
            common_tic_report_params['common_precursor_tic_group_cv_description'] = "Common Precursor TIC CV% calculated across groups"

    return common_tic_report_params

def miscleavage_plot(dig_df, miscleavage_threshold, groupwise_comparison, groups):

    if groupwise_comparison:
        dig_df['Group'] = dig_df['Filename'].apply(groupname, args=[groups, ])
        dig = px.bar(dig_df, x='Filename', y='0 missed cleavage percentage', title="Percentage of No Missed Cleavages", color=dig_df['Group'].tolist())
    else:
        dig = px.bar(dig_df, x='Filename', y='0 missed cleavage percentage', title="Percentage of No Missed Cleavages")

    if miscleavage_threshold:
        dig.add_hline(y=miscleavage_threshold, line_dash="dot", annotation_text=f"No Missed Cleavage Percentage Threshold = {miscleavage_threshold}")

    dig.update_xaxes(tickfont_size=8)
    dig.update_layout(margin=dict(l=20, r=20, t=20, b=20))
    dig.update_layout(title={'font': {'size': 9}})

    dig_graph = plotly.io.to_html(dig, include_plotlyjs=True, full_html=False)

    miscleavage_report_params = {'percent_miscleavage_plot': dig_graph,
                                'percent_miscleavage_description': "0 miscleaved peptides found in each sample."}

    return miscleavage_report_params

def selected_peptide_plots(df_level, filenames, level, level2, coverage_threshold):

    #Intensity Distribution Plot
    int_level = df_level[[level] + filenames]
    df_int = pd.melt(int_level, id_vars=[level], var_name="Filename", value_name="Intensity")

    if level2 == "iRT":
        int_plot = px.line(df_int, x='Filename', y="Intensity", title=f"Intensity Distribution Across iRT {level}s", color=level)
    else:
        int_plot = px.line(df_int, x='Filename', y="Intensity", title=f"Intensity Distribution Across {level}s", color=level)

    int_plot.update_xaxes(tickfont_size=8)
    int_plot.update_layout(margin=dict(l=20, r=20, t=20, b=20))
    int_plot.update_layout(title={'font': {'size': 9}})
    int_dist = plotly.io.to_html(int_plot, include_plotlyjs=True, full_html=False)

    #Peptide/Precursor Coverage Across Samples
    int_cov = df_level[[level, 'Coverage %']]
    if level2 == "iRT":
        cov = px.bar(int_cov, x=level, y="Coverage %", title=f"Coverage Percentage of iRT {level}s", color=int_cov[level].tolist())
    else:
        cov = px.bar(int_cov, x=level, y="Coverage %", title=f"Coverage Percentage of {level}s", color=int_cov[level].tolist())
    cov.add_hline(y=coverage_threshold)
    cov.update_xaxes(tickfont_size=8)
    cov.update_layout(
        margin=dict(l=20, r=20, t=20, b=20),
    )
    cov_plot = plotly.io.to_html(cov, include_plotlyjs=True, full_html=False)

    if level2 == "iRT":
        selected_peptide_report_params = {'irt_intensity_plot': int_dist,
                                          'irt_intensity_coverage_plot': cov_plot,
                                          'irt_intensity_description': 'iRT Intensity Distribution'}
    else:
        selected_peptide_report_params = {'selected_peptide_intensity_plot': int_dist,
                                          'selected_peptide_intensity_coverage_plot': cov_plot,
                                          'selected_peptide_intensity_description': 'User-Given Peptides Intensity Distribution'}

    return selected_peptide_report_params

def cumulative_freq_graph(protein_level, peptide_level, precursor_level, pt_cv_sum, pep_cv_sum, pre_cv_sum):

    if protein_level and peptide_level and precursor_level:
        df = pd.merge(pt_cv_sum, pep_cv_sum, on="CV%")
        df = pd.merge(df, pre_cv_sum, on="CV%")

        cv_sum = df.melt(id_vars=["CV%"],
                    var_name="Label",
                    value_name="Cumulative Frequency %")

        cv_line = px.line(cv_sum, x='CV%', y="Cumulative Frequency %", title="Number of Proteins, Peptides and Precursors under CV% (Across all Samples)", color="Label", line_shape="spline")
        cv_line.update_xaxes(tickfont_size=8)
        cv_line.update_layout(
                margin=dict(l=20, r=20, t=20, b=20)
        )

    elif protein_level and peptide_level:
        df = pd.merge(pt_cv_sum, pep_cv_sum, on="CV%")

        cv_sum = df.melt(id_vars=["CV%"],
                    var_name="Label",
                    value_name="Cumulative Frequency %")

        cv_line = px.line(cv_sum, x='CV%', y="Cumulative Frequency %", title="Number of Proteins and Peptides under CV% (Across all Samples)", color="Label", line_shape="spline")
        cv_line.update_xaxes(tickfont_size=8)
        cv_line.update_layout(
                margin=dict(l=20, r=20, t=20, b=20)
        )

    elif protein_level and precursor_level:
        df = pd.merge(pt_cv_sum, pre_cv_sum, on="CV%")

        cv_sum = df.melt(id_vars=["CV%"],
                    var_name="Label",
                    value_name="Cumulative Frequency %")

        cv_line = px.line(cv_sum, x='CV%', y="Cumulative Frequency %", title="Number of Proteins and Precursors under CV% (Across all Samples)", color="Label", line_shape="spline")
        cv_line.update_xaxes(tickfont_size=8)
        cv_line.update_layout(
                margin=dict(l=20, r=20, t=20, b=20)
        )

    elif peptide_level and precursor_level:
        df = pd.merge(pep_cv_sum, pre_cv_sum, on="CV%")

        cv_sum = df.melt(id_vars=["CV%"],
                    var_name="Label",
                    value_name="Cumulative Frequency %")

        cv_line = px.line(cv_sum, x='CV%', y="Cumulative Frequency %", title="Number of Peptides and Precursors under CV% (Across all Samples)", color="Label", line_shape="spline")
        cv_line.update_xaxes(tickfont_size=8)
        cv_line.update_layout(
                margin=dict(l=20, r=20, t=20, b=20)
        )

    elif protein_level:
        cv_line = px.line(pt_cv_sum, x='CV%', y="Protein Cumulative Frequency %", title="Number of Proteins under CV% (Across all Samples)", line_shape="spline")
        cv_line.update_xaxes(tickfont_size=8)
        cv_line.update_layout(
                margin=dict(l=20, r=20, t=20, b=20)
        )

    elif peptide_level:
        cv_line = px.line(cv_sum, x='CV%', y="Peptide Cumulative Frequency %", title="Number of Peptides under CV% (Across all Samples)", line_shape="spline")
        cv_line.update_xaxes(tickfont_size=8)
        cv_line.update_layout(
                margin=dict(l=20, r=20, t=20, b=20)
        )

    elif precursor_level:
        cv_line = px.line(cv_sum, x='CV%', y="Precursor Cumulative Frequency %", title="Number of Precursors under CV% (Across all Samples)", line_shape="spline")
        cv_line.update_xaxes(tickfont_size=8)
        cv_line.update_layout(
                margin=dict(l=20, r=20, t=20, b=20)
        )

    cumfreq_cv_line = plotly.io.to_html(cv_line, include_plotlyjs=True, full_html=False)

    cumfreq_report_params = {'cumulative_frequency_plot': cumfreq_cv_line,
                             'cumulative_frequency_description': "Cumulative Frequency % of CV%"}

    return cumfreq_report_params

#-------------------------------------------------------------------- MAIN FUNCTION --------------------------------------------------------------------------------

def calculate_idbased_metrics(out_dir, reportname, input_dict, threshold_dict, groups, groupwise_comparison):

    protein_report_params = {}
    peptide_report_params = {}
    precursor_report_params = {}

    if input_dict['Protein Level']:

        pt_level = pd.read_csv(input_dict['Protein Level'], sep="\t")
        pt_level = pt_level[pt_level.columns.tolist()].replace({'0':np.nan, 0:np.nan})

        filenames = pt_level.columns.tolist()
        filenames.remove("Protein")

        #protein quant and intensity CVs
        pt_quant, pt_grouped_quant = get_quant(pt_level, filenames, threshold_dict['Protein Threshold'], "Protein", groupwise_comparison, groups)
        pt_level_cv, pt_cv_sum, pt_grouped_cv = intensity_cvs(pt_level, threshold_dict['CV Percent Threshold'], threshold_dict['Data Percent Threshold'], filenames, "Protein", groupwise_comparison, groups)

        #getting protein level report parameters - plots + descriptions
        pt_quant_report_params = get_quant_plot(pt_quant, threshold_dict['Protein Threshold'], "Protein", groupwise_comparison, groups)
        #pt_intensity_cv_report_params = intensity_cv_graphs(pt_cv_sum, pt_grouped_cv, "Protein", groupwise_comparison, threshold_dict['CV Percent Threshold'], threshold_dict['Data Percent Threshold'])
        if groupwise_comparison:
            pt_intensity_cv_report_params = intensity_cv_graphs(pt_cv_sum, pt_grouped_cv, "Protein", groupwise_comparison, threshold_dict['CV Percent Threshold'], threshold_dict['Data Percent Threshold'])
            protein_pca_report_params = pca_plot(pt_level,"Protein",filenames, groups)
        else:
            pt_intensity_cv_report_params = {}
            protein_pca_report_params = {}

        protein_report_params = dict(tuple(pt_quant_report_params.items()) + tuple(pt_intensity_cv_report_params.items()) + tuple(protein_pca_report_params.items()))

        #saving dataframes to excel document
        protein_report_writer = pd.ExcelWriter(f"{out_dir}/{reportname}_ProteinLevel_QC_Report.xlsx", engine='xlsxwriter')
        pt_quant.to_excel(protein_report_writer, index=False, sheet_name='Protein Quant Summary')
        if groupwise_comparison and threshold_dict['Protein Threshold']:
            pt_grouped_quant.to_excel(protein_report_writer, index=False, sheet_name='Groupwise Protein Quant')
        pt_level_cv.to_excel(protein_report_writer, index=False, sheet_name='Protein Level CV')
        if groupwise_comparison:
            pt_grouped_cv.to_excel(protein_report_writer, index=False, sheet_name='Protein CV Group Summary')
        protein_report_writer.save()

        #getting protein overall sample dataframe
        if threshold_dict["Protein Threshold"]:
            pt_sample_df = pt_quant[['Filename', f'Protein Threshold = {threshold_dict["Protein Threshold"]}']]
        else:
            pt_sample_df = ""

        #getting protein group overall dataframe
        if groupwise_comparison:
            pt_group_df = pd.merge(pt_grouped_quant, pt_grouped_cv, on="Group")
            pt_group_df = pt_group_df[["Group",
                               "Protein Threshold QC Status",
                               f"{threshold_dict['Data Percent Threshold']}% Proteins <= {threshold_dict['CV Percent Threshold']}% CV"]]


    if input_dict['Peptide Level']:

        pep_level = pd.read_csv(input_dict['Peptide Level'], sep="\t")
        pep_level = pep_level[pep_level.columns.tolist()].replace({'0':np.nan, 0:np.nan})

        filenames = pep_level.columns.tolist()
        filenames.remove("Protein")
        filenames.remove("Peptide")

        #peptide quant, intensity_cvs, tic, missed cleavage percentage, peptide intensity distribution
        pep_quant, pep_grouped_quant = get_quant(pep_level, filenames, threshold_dict['Peptide Threshold'], "Peptide", groupwise_comparison, groups)
        pep_level_cv, pep_cv_sum, pep_grouped_cv = intensity_cvs(pep_level, threshold_dict['CV Percent Threshold'], threshold_dict['Data Percent Threshold'], filenames, "Peptide", groupwise_comparison, groups)
        pep_tic, pep_grouped_tic = common_tic(pep_level, "Peptide", threshold_dict['TIC CV Threshold'], filenames, groups, groupwise_comparison)

        if groupwise_comparison:
            pep_group_df = pd.merge(pep_grouped_quant, pep_grouped_cv, on="Group")
            pep_group_df = pd.merge(pep_group_df, pep_grouped_tic, on="Group")

        if threshold_dict['Enzyme']:
            dig_df, dig_grouped = miscleavage(pep_level, threshold_dict['Enzyme'], threshold_dict['Miscleavage Threshold'], filenames, groups, groupwise_comparison)
            if groupwise_comparison:
                pep_group_df = pd.merge(pep_group_df, dig_grouped, on="Group")
        else:
            dig_df, dig_grouped = ""

        if threshold_dict['iRT Label'] or input_dict['Peptide List']:
            if input_dict['Peptide List']:
                peptide_list_df = pd.read_csv(input_dict['Peptide List'], sep="\t")
            else:
                peptide_list_df = ""
            irt_level, irt_plots, selected_pep_df = selected_peps(pep_level, "Peptide", threshold_dict['Coverage Threshold'], filenames, threshold_dict['iRT Label'], input_dict['Peptide List'], peptide_list_df)

        #peptide sample dataframe
        if threshold_dict['Peptide Threshold'] and threshold_dict['Enzyme'] and threshold_dict['Miscleavage Threshold']:
            pep_sample_df = pd.merge(pep_quant, dig_df, on="Filename")
            pep_sample_df = pep_sample_df[['Filename', f'Peptide Threshold = {threshold_dict["Peptide Threshold"]}', '0 missed cleavage QC Status']]
        else:
            pep_sample_df = ""

        #getting peptide level report parameters - plots + descriptions
        pep_quant_report_params = get_quant_plot(pep_quant, threshold_dict['Peptide Threshold'], "Peptide", groupwise_comparison, groups)
        #pep_intensity_cv_report_params = intensity_cv_graphs(pep_cv_sum, pep_grouped_cv, "Peptide", groupwise_comparison, threshold_dict['CV Percent Threshold'], threshold_dict['Data Percent Threshold'])
        if groupwise_comparison:
            pep_intensity_cv_report_params = intensity_cv_graphs(pep_cv_sum, pep_grouped_cv, "Peptide", groupwise_comparison, threshold_dict['CV Percent Threshold'], threshold_dict['Data Percent Threshold'])
            peptide_pca_report_params = pca_plot(pep_level, "Peptide", filenames, groups)
        else:
            pep_intensity_cv_report_params = {}
            peptide_pca_report_params = {}

        pep_common_tic_report_params = common_tic_plot(pep_tic, pep_grouped_tic, "Peptide", threshold_dict['TIC CV Threshold'], groupwise_comparison, groups)
        if threshold_dict['Enzyme']:
            miscleavage_report_params = miscleavage_plot(dig_df, threshold_dict['Miscleavage Threshold'], groupwise_comparison, groups)
        else:
            miscleavage_report_params = {}

        if threshold_dict['iRT Label']:
            irt_report_params = selected_peptide_plots(irt_level, filenames, "Peptide", "iRT", threshold_dict['Coverage Threshold'])
        else:
            irt_report_params = {}

        if input_dict['Peptide List']:
            selected_peptide_report_params = selected_peptide_plots(selected_peptide_df, filenames, "Peptide", "Peptide List", threshold_dict['Coverage Threshold'])
        else:
            selected_peptide_report_params = {}

        peptide_report_params = dict(tuple(pep_quant_report_params.items()) +
                                tuple(pep_intensity_cv_report_params.items()) +
                                tuple(peptide_pca_report_params.items()) +
                                tuple(pep_common_tic_report_params.items()) +
                                tuple(miscleavage_report_params.items()) +
                                tuple(irt_report_params.items()) +
                                tuple(selected_peptide_report_params.items()))

        #saving dataframes to excel document
        peptide_report_writer = pd.ExcelWriter(f"{out_dir}/{reportname}_PeptideLevel_QC_Report.xlsx", engine='xlsxwriter')
        pep_quant.to_excel(peptide_report_writer, index=False, sheet_name='Peptide Quant Summary')
        if groupwise_comparison and threshold_dict['Peptide Threshold']:
            pep_grouped_quant.to_excel(peptide_report_writer, index=False, sheet_name='Groupwise Peptide Quant')
        pep_level_cv.to_excel(peptide_report_writer, index=False, sheet_name='Peptide Level CV')
        if groupwise_comparison:
            pep_grouped_cv.to_excel(peptide_report_writer, index=False, sheet_name='Peptide CV Group Summary')
        pep_tic.to_excel(peptide_report_writer, index=False, sheet_name='Common Peptide TIC')
        if groupwise_comparison:
            pep_grouped_tic.to_excel(peptide_report_writer, index=False, sheet_name='Common Peptide TIC Group CV')
        if threshold_dict['Enzyme']:
            dig_df.to_excel(peptide_report_writer, index=False, sheet_name='Miscleavage Threshold')
        if threshold_dict['iRT Label']:
            irt_level.to_excel(peptide_report_writer, index=False, sheet_name='iRT Peptide Intensity')
        if input_dict['Peptide List']:
            selected_pep_df.to_excel(peptide_threshold_dict, index=False, sheet_name='Selected Peptide Intensity')
        peptide_report_writer.save()

        #getting peptide group overall dataframe
        if groupwise_comparison:
            if threshold_dict['Enzyme']:
                pep_group_df = pep_group_df[["Group", "Peptide Threshold QC Status",
                            f"{threshold_dict['Data Percent Threshold']}% Peptides <= {threshold_dict['CV Percent Threshold']}% CV",
                            "Common Peptide TIC QC Status", "0 Miscleaved Peptides QC Status"]]
            else:
                pep_group_df = pep_group_df[["Group", "Peptide Threshold QC Status",
                            f"{threshold_dict['Data Percent Threshold']}% Peptides <= {threshold_dict['CV Percent Threshold']}% CV",
                            "Common Peptide TIC QC Status"]]

    if input_dict['Precursor Level']:

        pre_level = pd.read_csv(input_dict['Precursor Level'], sep="\t")
        pre_level = pre_level[pre_level.columns.tolist()].replace({'0':np.nan, 0:np.nan})

        filenames = pre_level.columns.tolist()
        filenames.remove("Protein")
        filenames.remove("Peptide")
        filenames.remove("Precursor")

        #precursor quant, intensity_cvs, tic, missed cleavage percentage, precursor intensity distribution
        pre_quant, pre_grouped_quant = get_quant(pre_level, filenames, threshold_dict['Precursor Threshold'], "Precursor", groupwise_comparison, groups)
        pre_level_cv, pre_cv_sum, pre_grouped_cv = intensity_cvs(pre_level, threshold_dict['CV Percent Threshold'], threshold_dict['Data Percent Threshold'], filenames, "Precursor", groupwise_comparison, groups)
        pre_tic, pre_grouped_tic = common_tic(pre_level, "Precursor", threshold_dict['TIC CV Threshold'], filenames, groups, groupwise_comparison)

        if groupwise_comparison:
            pre_group_df = pd.merge(pre_grouped_quant, pre_grouped_cv, on="Group")

        if not input_dict['Peptide Level']:
            if groupwise_comparison:
                pre_group_df = pd.merge(pre_group_df, pre_grouped_tic, on="Group")

            if threshold_dict['Enzyme']:
                dig_df, dig_grouped = miscleavage(pre_level, threshold_dict['Enzyme'], threshold_dict['Miscleavage Threshold'], filenames, groups, groupwise_comparison)
                if groupwise_comparison:
                    pre_group_df = pd.merge(pre_group_df, dig_grouped, on="Group")
            else:
                dig_df, dig_grouped = ""

            if threshold_dict['iRT Label'] or input_dict['Peptide List']:

                if input_dict['Peptide List']:
                    peptide_list_df = pd.read_csv(input_dict['Peptide List'], sep="\t")
                else:
                    peptide_list_df = ""

                irt_level, irt_plots, selected_pep_df = selected_peps(pre_level, "Precursor", threshold_dict['Coverage Threshold'], filenames, threshold_dict['iRT Label'], input_dict['Peptide List'], peptide_list_df)

        #peptide sample dataframe
        if threshold_dict['Precursor Threshold']:
            pre_sample_df = pre_quant[['Filename', f'Precursor Threshold = {threshold_dict["Precursor Threshold"]}']]
            if not input_dict['Peptide Level']:
                if threshold_dict['Enzyme'] and threshold_dict['Miscleavage Threshold']:
                    pre_sample_df = pd.merge(pre_sample_df, dig_df, on="Filename")
                    pre_sample_df = pre_sample_df[['Filename', f'Precursor Threshold = {threshold_dict["Precursor Threshold"]}', '0 missed cleavage QC Status']]
        else:
            pre_sample_df = ""

        #getting precursor level report parameters - plots + descriptions
        precursor_quant_report_params = get_quant_plot(pre_quant, threshold_dict['Precursor Threshold'], "Precursor", groupwise_comparison, groups)
        #precursor_intensity_cv_report_params = intensity_cv_graphs(pre_cv_sum, pre_grouped_cv, "Precursor", groupwise_comparison, threshold_dict['CV Percent Threshold'], threshold_dict['Data Percent Threshold'])
        if groupwise_comparison:
            precursor_intensity_cv_report_params = intensity_cv_graphs(pre_cv_sum, pre_grouped_cv, "Precursor", groupwise_comparison, threshold_dict['CV Percent Threshold'], threshold_dict['Data Percent Threshold'])
            precursor_pca_report_params = pca_plot(pre_level, "Precursor", filenames, groups)
        else:
            precursor_intensity_cv_report_params = {}
            precursor_pca_report_params = {}

        if not input_dict['Peptide Level']:

            pre_common_tic_report_params = common_tic_plot(pre_tic, pre_grouped_tic, "Precursor", threshold_dict['TIC CV Threshold'], groupwise_comparison, groups)

            if threshold_dict['Enzyme']:
                miscleavage_report_params = miscleavage_plot(dig_df, threshold_dict['Miscleavage Threshold'], groupwise_comparison, groups)
            else:
                miscleavage_report_params = {}

            if threshold_dict['iRT Label']:
                irt_report_params = selected_peptide_plots(irt_level, filenames, "Precursor", "iRT", threshold_dict['Coverage Threshold'])
            else:
                irt_report_params = {}

            if input_dict['Peptide List']:
                selected_peptide_report_params = selected_peptide_plots(selected_peptide_df, filenames, "Precursor", "Peptide List", threshold_dict['Coverage Threshold'])
            else:
                selected_peptide_report_params = {}

        precursor_report_params = dict(tuple(precursor_quant_report_params.items()) +
                                tuple(precursor_intensity_cv_report_params.items()) +
                                tuple(precursor_pca_report_params.items()) +
                                tuple(pre_common_tic_report_params.items()) +
                                tuple(miscleavage_report_params.items()) +
                                tuple(irt_report_params.items()) +
                                tuple(selected_peptide_report_params.items()))

        #saving dataframes to excel document
        precursor_report_writer = pd.ExcelWriter(f"{out_dir}/{reportname}_PrecursorLevel_QC_Report.xlsx", engine='xlsxwriter')
        pre_quant.to_excel(precursor_report_writer, index=False, sheet_name='Precursor Quant Summary')
        if groupwise_comparison and threshold_dict['Precursor Threshold']:
            pre_grouped_quant.to_excel(precursor_report_writer, index=False, sheet_name='Groupwise Precursor Quant')
        pre_level_cv.to_excel(precursor_report_writer, index=False, sheet_name='Precursor Level CV')
        if groupwise_comparison:
            pre_grouped_cv.to_excel(precursor_report_writer, index=False, sheet_name='Precursor CV Group Summary')
        pre_tic.to_excel(precursor_report_writer, index=False, sheet_name='Common Precursor TIC')
        if groupwise_comparison:
            pre_grouped_tic.to_excel(precursor_report_writer, index=False, sheet_name='Common Precursor TIC Group CV')

        if not input_dict['Peptide Level']:
            if threshold_dict['Enzyme']:
                dig_df.to_excel(precursor_report_writer, index=False, sheet_name='Miscleavage Threshold')
            if threshold_dict['iRT Label']:
                irt_level.to_excel(precursor_report_writer, index=False, sheet_name='iRT Precursor Intensity')
            if input_dict['Peptide List']:
                selected_pep_df.to_excel(precursor_report_writer, index=False, sheet_name='Selected Precursor Intensity')

        precursor_report_writer.save()

        if groupwise_comparison:
            if not input_dict['Peptide Level']:
                if threshold_dict['Enzyme']:
                    pre_group_df = pre_group_df[["Group", "Precursor Threshold QC Status",f"{threshold_dict['Data Percent Threshold']}% Precursors <= {threshold_dict['CV Percent Threshold']}% CV",
                    "Common Precursor TIC QC Status", "0 Miscleaved Peptides QC Status"]]
                else:
                    pre_group_df = pre_group_df[["Group", "Precursor Threshold QC Status",f"{threshold_dict['Data Percent Threshold']}% Precursors <= {threshold_dict['CV Percent Threshold']}% CV",
                    "Common Precursor TIC QC Status"]]

            else:
                pre_group_df = pre_group_df[["Group", "Precursor Threshold QC Status",f"{threshold_dict['Data Percent Threshold']}% Precursors <= {threshold_dict['CV Percent Threshold']}% CV"]]

    #adding cumulative frequency graph
    if not input_dict['Protein Level']:
        pt_group_df = ""
        pt_sample_df = ""
        pt_cv_sum = ""

    if not input_dict['Peptide Level']:
        pep_group_df = ""
        pep_sample_df = ""
        pep_cv_sum = ""

    if not input_dict['Precursor Level']:
        pre_group_df = ""
        pre_sample_df = ""
        pre_cv_sum = ""


    #overall sample dataframe
    overall_sample_df = get_sample_df(input_dict['Protein Level'], input_dict['Peptide Level'], input_dict['Precursor Level'], pt_sample_df, pep_sample_df, pre_sample_df, threshold_dict)

    #overall grouped dataframe
    if groupwise_comparison:
        overall_group_df = get_overall_df(input_dict['Protein Level'], input_dict['Peptide Level'], input_dict['Precursor Level'], pt_group_df, pep_group_df, pre_group_df)
    else:
        overall_group_df = ""

    #Report Parameters
    cumfreq_report_params = cumulative_freq_graph(input_dict['Protein Level'], input_dict['Peptide Level'], input_dict['Precursor Level'], pt_cv_sum, pep_cv_sum, pre_cv_sum)
    report_parameters = dict(tuple(protein_report_params.items())
                        + tuple(peptide_report_params.items())
                        + tuple(precursor_report_params.items()) +
                        tuple(cumfreq_report_params.items()))

    return (overall_sample_df, overall_group_df, report_parameters)
