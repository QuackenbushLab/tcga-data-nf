library(TCGAbiolinks)
library(optparse)
library(SummarizedExperiment)

option_list = list(
  make_option(c("-p", "--project_id"), type="character", default=NULL, 
              help="project_id", metavar="character"),
  make_option(c("-d", "--downloads_dir"), type="character", default=NULL, 
              help="downloads_dir", metavar="character"),
  make_option(c("--analysis_workflow_type"), type="character", default="ASCAT3", 
              help="", metavar="character"),  
  make_option(c("--output_table"), type="character", default=NULL, 
              help="", metavar="character"),
  make_option(c("--output_rds"), type="character", default=NULL, 
              help="", metavar="character"),
  make_option(c("--sample_list"), type="character", default=' ', 
              help="test file with sample names", metavar="character")    
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
sample_list = opt$sample_list
print(opt)

project_id = opt$project_id # args[1]
analysis_workflow_type = opt$analysis_workflow_type # always methylation_beta_value at this point
print(paste('Downloading',project_id,',with workflow:',analysis_workflow_type))
output_table = opt$output_table
output_rds = opt$output_rds

#barcode = listSamples
if (nchar(sample_list) > 3 && !grepl("NA$", sample_list)){
  barcodes = read.table(opt$sample_list, header = FALSE, sep = "", dec = ".")
  
  queried <- GDCquery(
    project = project_id, 
    data.category = "Copy Number Variation", 
    data.type = "Gene Level Copy Number",
    workflow.type = analysis_workflow_type,
    barcode = barcodes$V1)
  } else {
  queried <- GDCquery(
    project = project_id, 
    data.category = "Copy Number Variation", 
    data.type = "Gene Level Copy Number",
    workflow.type = analysis_workflow_type
  ) }


GDCdownload(queried,files.per.chunk = 50)
cn_rds <- GDCprepare(queried)



saveRDS(cn_rds, output_rds)
write.csv( assays(cn_rds)$copy_number, file = output_table)