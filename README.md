# Scripts and data from Griffiths (2018) in review

# Physiological Analysis

Script: Respiration_analysis.R

Input files: Resp_all_data.csv, Resp_day9_data.csv, Resp_day29_data.csv

Description: Script contains code for normalizing respiration data and performing an ANOVA and a Tukey’s post-hoc test


Script: lipid_protein_ananlysis.R

Input files: lipid_protein_data.csv

Description: Script contains code for normalizing lipid and protein data and performing an ANOVA



# Differentially Expressed Gene Analysis

Script: DeSeq2.R

Input files: Orthoblast_RSEM_merged_matrix, column_2transcriptomes.txt

Output files: output_DEG_GOLtp1.txt, output_DEG_GOLtp2.txt, output_DEG_PACtp1.txt, output_DEG_PACtp2.txt

Description: Script contains analysis using DeSeq2. The merged matrix file contains the total counts of contains for each “Orthogroup”. The column file gives an explanation for the headers in the matrix file (the treatment conditions and population names for each read file. Output files contain the pvalues for log fold changes in Orthogroups in response to pH for each population timepoint (day 9 and 29).



# Functional Enrichment Analysis

Script: Please see scripts and explanations located at: https://github.com/z0on/GO_MWU\

Input files: Orthoblast_interproresults_nonredun.csv, pvalue_GOL_tp1, pvalue_GOL_tp2, pvalue_PAC_tp1, pvalue_PAC_tp2

Output files: BP_pvalue_GOL_tp1, BP_pvalue_GOL_tp2, BP_pvalue_PAC_tp1, BP_pvalue_PAC_tp2, MF_pvalue_GOL_tp1, MF_pvalue_GOL_tp2, MF_pvalue_PAC_tp1, MF_pvalue_PAC_tp2

Description: Script contains analysis for functional enrichment of logfold contig changes. The Orthoblast_interproresults_nonredun.csv input file contains the GO terms associated with each Orthogroup. The pvalue input files contain the signed logfold pvalues for expression changes derived from the DEG analysis using DeSeq2. Output files correspind to each input file name and whether the Biological Processes (BP) or Molecular Functions (MF) were analyzed.
