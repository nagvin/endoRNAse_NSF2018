# endoRNAse_NSF2018
Scripts used in the RNA degradome analysis to identify 5'P sites that are dependent on Arabidopsis DNE1 endoRNase

Publication:
Nagarajan et al.,
RNA degradome analysis reveals DNE1 endoribonuclease is required for the turnover of diverse mRNA substrates in Arabidopsis
Plant Cell 2023 May 29;35(6):1936-1955. doi: 10.1093/plcell/koad085.


Analysis details:
---------------------
RNA degradome libraries (using the PARE and GMUCT methods) were analyzed as follows
Pre-processing steps: 
1. QC checks with FastQC

2. Trim adapter, crop to 50 nt (GMUCT) and 20 nt (PARE) and remove low quality reads with Trimmomatic 

Overview of the ComPARE computational pipeline followed by Hurtig et al., 2021 (PNAS) – deviations are shown with an asterisk(*) 

Data analysis steps:
1. Align reads to TAIR10 genome (v49) using Tophat v2.0 

2. Calculate 5’P and transcript read counts per million (CPM) using BamCoverage (Deeptools)

3. Normalize 5’P (first nt position) reads to depth

4. Assign gene and transcript attributes to 5’P sites from GFF files (Bedtools intersect)*

5. Calculate log2 fold changes between 5’P site CPM values from Col-0, xrn4 and dne1xrn4 (or dbl)

6. Extract 20 nt sequences flanking 5’P site (Bedtools getfasta)*

##############################
List of scripts:
##############################
1. PARE_score_filter_v11.pl: Reads bedgraph file from the comPARE analysis pipeline to extract 5'P sites >= 1 CPM from all replicates. Removes sites within chloroplast and mitochondria. Output is a text file containing, 5'P sites, coordinates, orientation, 5'P site abundance (CPM) for Col-0, xrn4 and dne1xrn4  and log2 fold-changes at 5'P sites between xrn4 vs Col-0 and xrn4 vs dne1xrn4.   

2. PARE_annotate_filter.pl: Reads output file of bedtools intersect that includes genomic and cDNA coordinates of 5'P sites. Output includes coordinates and sequence of the region flanking the 5'P sites from Bedtools.  

3. comPARE_datasets_v3.pl: 5'P site coordinates and corresponding sequences are annotated with transcript and genes names, Cap PARE (mRNA decapping) sites (Nagarajan et al., 2019),  transcript abundances (RNA-seq) results, and datasets from multiple studies (e.g. NMD mutants). The output dataset is presented in Supplemental Data Set S1 (Seedling) and S3 (Leaf).  

4. metagene_v1.pl. Metagene analysis to identify the positional information of DNE1-dependent 5'P sites across the transcript based on TAIR10 transcript annotations. Output was used in Figure 3B in the manuscript 


