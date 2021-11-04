#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
run.preProcessed <- as.logical(args[1])

#run.preProcessed <- TRUE

data.dir <- "../data"
res.dir <- "../results"

message("Loading packages")
suppressMessages(suppressWarnings(source("0.utils.R")))
start_time <- Sys.time()
message("Start 1.1.harmonize_extract_classify_Affymetrix_LMS.R")
source("1.1.harmonize_extract_classify_Affymetrix_LMS.R")
end_time <- Sys.time()
print(end_time - start_time)

message("Start 1.2.parse_GSEA_iCistarget_results.R")
source("1.2.parse_GSEA_iCistarget_results.R")

end_time <- Sys.time()
print(end_time - start_time)


message("Start 1.3.classify_patients_clinical_enrichment.R")
source("1.3.classify_patients_clinical_enrichment.R")

end_time <- Sys.time()
print(end_time - start_time)

message("Start 1.4.GTEX_analysis.R")
source("1.4.GTEX_analysis.R")

end_time <- Sys.time()
print(end_time - start_time)


message("Start 2.1.copynumber_analysis.R")
source("2.1.copynumber_analysis.R")

end_time <- Sys.time()
print(end_time - start_time)


message("Start 3.1.TCGA_miRNA_DE_analysis.R")
source("3.1.TCGA_miRNA_DE_analysis.R")

end_time <- Sys.time()
print(end_time - start_time)


message("Start 3.2.ICGC_miRNA_DE_analysis.R")
source("3.2.ICGC_miRNA_DE_analysis.R")

end_time <- Sys.time()
print(end_time - start_time)

message("Start 3.3.Compare_miRNA_cohorts.R")
source("3.3.Compare_miRNA_cohorts.R")

end_time <- Sys.time()
print(end_time - start_time)

message("Start 3.4.PANCAN_miRNA_analysis.R")
source("3.4.PANCAN_miRNA_analysis.R")

end_time <- Sys.time()
print(end_time - start_time)

message("Start 3.5.miRNA_mRNA_interaction_analysis.R")
source("3.5.miRNA_mRNA_interaction_analysis.R")

end_time <- Sys.time()
print(end_time - start_time)

end_time <- Sys.time()
print(end_time - start_time)
message("Start 4.1.mutations_analysis.R")
source("4.1.mutations_analysis.R")
end_time <- Sys.time()
print(end_time - start_time)
