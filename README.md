# estimating_admixture

----------------------------------------------------------------------------------

Title: DNA, Self-Reported Ancestry, and Social Scientific Inquiry

Journal: Demography

Authors: Zhang and Trejo

Date: Jan 2026

----------------------------------------------------------------------------------

This replication package contains data and code to reproduce the analysis, tables, 
and figures included in the main text and supplementary material. 

Instructions: Prior to running the scripts, users will need to change the path 
information in line 111 of 02-merge_all.R and at the top of 03-create_ipw.do.

----------------------------------------------------------------------------------

Description of the script files: 


1. 00-create_GSPs: prepare genotype data and run ADMIXTURE to create GSPs
2. 00-functions.R: functions to load into 04-tables_figures.R
3. 01-create_geoanc.R: prepare self-reported geographic ancestry
4. 02-merge_all.R: prepare other survey data and merge with GSPs
5. 03-create_ipw.do: perform IPW weighting
6. 04-tables_figures.R: main analysis and create tables and figures
7. 05-create_maps.R: create maps
8. 06-PCA.R: perform PCA on genotype data
9. 07-compare_gsp&pc2.R: compare GSPs and PCs

----------------------------------------------------------------------------------

Description of the data files: 


1. geoanc_codebook.xlsx: codebook for the self-reported geographic ancestry variable
2. igsr_populations.tsv: coordinates information for the reference panel, downloaded from
https://www.internationalgenome.org/data-portal/population
3. multilevel_regions.xlsx: map the subregion categories to the corresponding region categories
4. reference_composition.txt: composition of the reference panel
