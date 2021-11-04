# Workflow analysis for leiomyosarcoma origin and oncogenesis

This workflow is composed of 11 R scripts organized in 4 main parts. The scripts are ordered and need to be run sequentially. The MASTER_script.R script allows to run the whole pipeline at once. The majority of the scripts depends on the previous script, no intermediate results are kept. Functions used in each part are written in a separate file _functions.R_.

The aim of the workflow is to reproduce the results presented in the manuscript.

## Data
### From public databases

- functional annotations from MSigDB in a _gmt_ format
- probes annotations for Affymetrix and Agilent platforms
- cytoband genomic position and chromosome size from hg19 genome version from UCSC
- miRNA-mRNA interaction from miRTarbase and miRecords formatted and saved in RData format
- TCGA PANCAN batch corrected miRNA expression
- TCGA PANCAN batch corrected miRNA clinical data formated and saved as RData
- TCGA SARC miRNA expression saved as RData
- TCGA RNA-seq expression
- TCGA clinical data formated saved in a RData object.
- TCGA gene level copy number (CN) combined saved in a RData object.
- TCGA + GTEX RNA-seq combined data saved as RData. This object contains also a vector with tissue of origin of the samples ordered according to the samples in the data.frame and a vector containing the colours associated with each tissue type.
- List of micro-array data from both Affymetrix and agilent platforms saved in a RData object
- Gene level copy number estimation from SNP6 and BAC array technologies saved in RData objects

### From this study 

- Annotation file for all patients analysed on micro-arrays
- ICGC clinical data
- ICGC RNA-seq gene level log2(FPKM+1)
- ICGC DMD isoform expression
- ICGC TP53 and RB1 genetic status and PTEN, DMD and ATRX location
- ICGC breakpoint prediction
- ICGC somatic mutation prediction formated to be pass to _MutationPatterns_ package

### Pre-computed intermediate results

- ICGC CN estimation data: precomputed from cn.MOPS output. Used script is available in pre_computed.R file
- Correlation in step 1.1: Consistency of gene expression profiles is measured by computing the correlation of each gene with itself between the 2 platforms and with all the other genes in both platforms. This analysis is time consuming (around 1h30) and a precomputed object is used by default. To re-compute it, change "TRUE" to "FALSE" in the _run_ file.
- Functional and regulatory feature enrichment: GSEA and iCistarget results. All folders contain html reports allowing to visualize the results. _report.html_ in icistarget\_\* folders and _pos\__ and _neg\_snapshot.html_ in GSEA/Affy\_\* folders. 
- Functional annotations by GSAn for top hLMS expressed genes (ribosomal genes excluded, list available in text file) and miRNA targets. A json file is available and can be loaded on GSAn web server to visualize the results.
- Cell lines log2(FPKM+1) transcriptome files \*\__FPKM\_log2.txt_


## Scripts
### Step 1: Transcriptome analysis
#### 1.1. Detection of LMS groups: harmonization and combination of Affymetrix and Agilent micro-arrays data 

This script aims at identify genes with consistent expression profiles between Agilent and Affymetrix platform through 87 patients analyzed on both platforms and harmonized the data to be able to combine them. All datasets are then merged and harmonized. Modules of co-expressed genes are defined and used in a hierarchical clustering. LMSs grouped in a single cluster are compared to others. Definition of hLMS and oLMS in Affymetrix cohort. Differentially expressed genes are defined and used to compute centroids.

* Produced plots: 
    - Figure 1A / S1B: _Figure1A_all_samples_clustering_
    - Figure 1B: _Figure1B_DE_Affymetrix_
    - Figure S1B: _FigureS1B_ Venn diagram 
    - Figure S1C: _FS1_example_transformation_87_ line plots
    - Figure S1C: _FS1_raw_heatmap_87_samples_ normalized 1st round
    - Figure S1C: _FS1_Qnorm_heatmap_87_samples_ normalized 2nd round
    - Figure S1C: _FS1_harmonized_heatmap_87_samples_ harmonized

* Text files produced:
    - functional enrichment of module of co-expression: hypergeometric test: _enrichment\_cluster\_multiplatform\_\*_
    - table with 2 columns: gene symbols and t-scores in decreasing order. This file can be used to reproduce GSEA analysis _Affy\_t\_scores\_OCEANCODE.rnk_ from command line: 

    java -cp gsea-3.0.jar -Xmx2048m xtools.gsea.GseaPreranked -gmx /data/ext/KEGG.symbols.gmt.txt -norm meandiv -nperm 1000 -rnk /results/Affy_t_scores_OCEANCODE.rnk -scoring_scheme weighted -rpt_label Affy_t_scores_OCEANCODE -create_svgs false -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out "/results/GSEA_LMS_DE" -gui false; done

    - Two files, with significantly up- and down- regulated genes _hLMS\_genes\_\*.txt_ which can be used to reproduce iCistarget results from web site: https://gbiomed.kuleuven.be/apps/lcb/i-cisTarget/

#### 1.2. Visualize results obtained with GSEA and iCistarget

This script allows to parse _html_ reports from GSEA and retrieve the normalized enrichment scores and FDR associated to enriched terms. The median t-scores (hLMS/oLMS) for all genes belonging to a given term is computed and the number of significantly differentially expressed genes are counted. Related enriched terms were manually grouped in _GSEA\_\*\_table\_group_ files in the _precomp_ folder. In iCistarget analysis, the _report.html_ are parsed. Manual grouping is also made at this step.

* produced plot:
    - Figure 1D: _Figure1D_GSEA_functional_annots_
    - Figure 1E: _Figure1E_cistarget_results_

#### 1.3. Classify patients from ICGC and TCGA cohorts

Correlation to centroids obtained in the previous step is computed for each TCGA and ICGC patients. A Gaussian mixture model is computed on distances to the hLMS centroid from all patients (Affymetric cohort included) and patients are classified. Clinical associations are then tested and metastasis free survival is evaluated.

* Produced plots: 
    - Figure 1F: _Figure1F_centroid_classification_mclust_
    - Figure 1C: _Figure1C_Affy_MFS_surv_
    - Figure 5A: _Figure5A_OCs_centroids.pdf_

* Text files produced:
    - Table 1: _Table1_Affymetrix_clinical.tab_
    - Table 2: _Table2_ICGC_TCGA_clinical.tab_
    - Table S2A: _TS2A_stat_table_all_cohorts.tab_

#### 1.4. Evaluate to which normal tissue each LMS group is the closest

Using TCGA and GTEX combined data, median gene expression is computed for the hLMS group and the top 100 are selected to cluster patients using tSNE method. Median gene expression of those genes is computed for all tissue groups and visualized with an heatmap.

* Produced plots: 
    - Figure S2A: _FS2A_tSNE_top100_hLMS_GTEX_
    - Figure S2B: _FS2B_heatmap_top100_median_tissue_


#### Functions
- _stat.probes_ : compute Pearson's correlations from 2 matrices containing the same patients and genes as input
_select.genes.by.cor_ : select genes having a correlation above >0.8 with itself between platforms or better then with any other gene.
- _harmonize_ : median normalization for the 87 patients
_plot.couple.gene_ : plot gene expression profiles for 2 genes for a given expression matrix
- _harmonize.2_ : generalized version of median normalization between platforms
- _import.clinical_ : Parse and format clinical annotations for micro-arrays analyzed patients
- _readGSEA_ : parse HTML GSEA output
- _import.icistarget_ : parse HTML iCistarget output
- _clinical.fisher_ : compute Fisher's exact test for categorical data, one _versus_ all other
- _t.score_ : compute t-scores and t-test p-values from gene expression data between 2 groups
- _computeCentroids_: compute centroids based on _cit.dfAggregate_
- _classify_LMS_: compute Spearman's correlation to centroids


### Step 2: Copy number variation analysis
#### 2.1. Analyze copy number alteration recurrences in hLMS and oLMS and focus on MYOCD locus

During this step, frequencies of homozygous (deletion), heterozygous (loss) losses, gains of one copy or amplifications (gains of more than one copy) are computed. These frequencies were computed in separated and merged cohorts, frequency profiles were compared between cohorts to evaluate the consistency using a linear regression model. Enrichment of gene alteration events on the merged cohort is analyzed using a Fisher's exact test comparing hLMS and oLMS. Enrichment of significant events per cytogenic band is estimated using Fisher's exact test. MYOCD expression level is categorized by cohort according to the whole transcriptome signal distribution.


* Produced plot:
    - Figure 3A: _Figure3A_heatmap_CN_puce_available_genes_ 
    - Figure 3B: _Figure3B_puces_avail_all_hLMS_
    - Figure 3B: _Figure3B_puces_avail_all_oLMS_
    - Figure 4A: _Figure4A_zoom_MYOCD_
    - Figure 4B: _Figure4B_MYOCD_expression_per_cohort_
    - Figure S4A: _FS4A_heatmap_rsquared_CN_
    - Figure S4B: _FS4B_MYOCD_CN_amplified_expression_
    - Figure S4C: _FS4C_cumul_MYOCD_expr_

* Text files produced:
    - Table S3A: _TS3_CN_per_cytoband.tab_
    - Table S3B: _TS3_CN_per_genes.tab_

#### Functions

- _copy.test_: Assign numeric values to predicted copy numbers (A: 4, amplification; G: 3, gain; N: 2, normal; P and L: 1, heterozygous loss; D: 0, homozygous loss). If a gene locus is assigned with more than one status, the lowest copy number is kept.
- _comp.rsquared_: compute linear regression per event type and return R-squared values
- _percentage.events_: compute percentage of events in hLMS and oLMS. The total number of patients per gene do not count missing values.
- _percentage.all_: compute global percentage of event. The total number of patients per gene do not count missing values.
- _import.copyNumber.by.genes_: compute counting of events per group and Fisher's exact test if parameter _fisher_ set to true.
- _enrich.CN_: function using _import.copyNumber.by.genes_ and allowing to choose the reference group with parameter _ref_.

### Step 3: miRNA transcriptome profiling
#### 3.1. Analysis of SARC TCGA normalized mature miRNA data

Differential expression is computed with miRComb package, using limma method and Holm's multi-testing correction. miRNA with less than one normalized reads in at least one group are discarded.

* Produced plot:
    - Figure 2A: _Figure2A_PCA_TCGA_miRNA_
    - Figure 2C: _Figure2C_heatmap_diff_miRNA_cr_TCGA_ 

#### 3.2. Analysis of ICGC miRNA data

Raw counts are submitted to edgeR package for filtering, normalization and differential expression analysis.

* Produced plot:
    - Figure 2A: _Figure2A_PCA_ICGC_miRNA_
    - Figure 2C: _Figure2C_heatmap_diff_miRNA_centered_reduced_ICGC_ 

#### 3.3. Comparison between TCGA and ICGC miRNA differential expression

Log Fol Changes (logFC) estimated during the previous steps are compared using a linear regression model. Results of each analysis are merged.
logFC and median expression of miRNA from the DIO3-DLK1 miRNA cluster on chromosome 14 which are not significantly differentially expressed (DE) are evaluated.

* Produced plot:
    - Figure 2A: _Fgure2A_XY_logFC_ICGC_TCGA_miRNA_
    - Figure S3AB: _FigureS3AB_boxplots_DIO3_DLK1_no_sig_

* Text files produced:
    - Table S2B: _miRNA_ICGC_TCGA_DE.tab_

#### 3.4. Pan cancer miRNA profiling

Differentially expressed miRNAs in both cohort and DIO3-DLK1 miRNAs are used to cluster patients from TCGA PANCAN data. Distribution of cancer histotypes in clusters is evaluated.

* Produced plot:
    - Figure 2B: _Figure2B_heatmap_diff_miRNA_all_TCGA_
    - Figure 2B: _Figure2B_PANCAN_cluster_diff_histotype_
    - Figure S3C: _FS3C_heatmap_DIO3_miRNA_all_TCGA_
    - Figure S3C: _FS3C_PANCAN_cluster_DIO3_histotype_
    - Figure S3C: _points_cluster_DIO3_DLK1_all_TCGA_

#### 3.5. Visualization of miRNA target functional mapping (GSAn)

miRComb algorithm is applied to TCGA and ICGC cohorts independently together with validated miRNA-mRNA interactions from public databases. Results were merged. Validated targets from databases which were found to significantly anti-correlate in both cohort were analyzed with GSAn method. This script helps in visualize this result.

* Produced plot:
    - Figure2D: _Figure2D_GSAN_miRNA_targets_

* Text files produced:
    - Table S2C: _miRComb_table.tab_
    - _test.mRNA.miRNA.txt_: selected targets that has been analyzed with GSAn

#### Functions
- _run_mirComb_: Apply classic miRComb pipeline: apply differential expression analysis on mRNAs and miRNAs (limma), select on significance (p-value < 0.01 and absolute logFC > 1), compute correlation between mRNA and miRNA expression and compute a p-value on Pearson's correlation score and correct it with FDR approach, add public database information. 


### Step 4: Single nucleotide variant analysis

The VCFs from ICGC containing somatic variants were analysed with MutationPattarns R package following their suggested pipeline: compute 96 trinucleotides matrix, compare with COSMIC signatures (cosine similarity) and estimate COSMIC signature contribution for each patients. Only those with a cosine similarity above 0.75 are reported.

Enrichment for specific allele status (TP53 and RB1) and protein cellular localization (DMD, PTEN, ATRX) was estimated using a Fisher's exact test comparing one category to all other possible categories. 

#### 4.1 Single nucleotide variant analysis

* Produced plot:
    - Figure 4C: _Figure4C_cat_alteration_categories_
    - Figure 4D: _Figure4D_DMD_isoforms_
    - Figure 4E: _Figure4E_loc_alteration_categories_
    - Figure S4D: _FigureS4D_TBM_BP_
    - Figure S4E: _FigureS4E_cosine_cosmic_signature_
    - Figure S4E: _FigureS4E_contribution_cosmic_signature_barplot_

#### Functions
- _plot.prop.mutation_: plot percentages of allele status or protein localization occurrence per LMS group.


