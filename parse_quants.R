library(tximport)
library(biomaRt)


# modifies the pheno data to include the path to the relevant quant files
get_pheno_data = function(){
  
  # get list of folders containing quants
  mcf7_quant_folders = Sys.glob(file.path(MCF7_QUANTS, 'SRR*'))
  t47d_quant_folders = Sys.glob(file.path(T47D_QUANTS, 'SRR*'))
  
  # get list of accessions from quant_folders (i.e. the last part of file path)
  mcf7_accessions = unlist(lapply(mcf7_quant_folders, basename))
  t47d_accessions = unlist(lapply(t47d_quant_folders, basename))
  
  # get list of quant files 
  mcf7_quant_files = Sys.glob(file.path(MCF7_QUANTS, '*/quant.sf' ))
  t47d_quant_files = Sys.glob(file.path(T47D_QUANTS, '*/quant.sf' ))
  
  mcf7_df = data.frame(mcf7_accession=mcf7_accessions, quant_file=mcf7_quant_files)
  t47d_df = data.frame(t47d_accession=t47d_accessions, quant_file=t47d_quant_files)
  
  pd_mcf7 = merge(PD_MCF7, mcf7_df, by.x = 'Run', by.y = 'mcf7_accession')
  pd_t47d = merge(PD_T47D, t47d_df, by.x = 'Run', by.y = 't47d_accession')
  
  return (list(mcf7=pd_mcf7, t47d=pd_t47d))
}



#get annotation
get_annotation = function(ensemble_ids){
  ensemble_ids_no_dot = gsub("\\.[0-9]$", "", ensemble_ids)
  

  
  if ( file.exists(ANNOTATION_DATA_FILE)){
    annotLookup = readRDS(ANNOTATION_DATA_FILE)
  } else {
    mart <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl', 
                    host = "www.ensembl.org",
                    ensemblRedirect = FALSE)
    annotLookup <- getBM(
      mart=mart,
      attributes=c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
      filter="ensembl_transcript_id",
      values=ensemble_ids_no_dot,
      # values="ENST00000464927",
      uniqueRows=TRUE)
  }
    return(annotLookup)
}

# read data using tximport. Input is the phenoData dataframe (i.e. SraRunTable) 
#  which includes additional column for path to file containing salmon quant file
read_data = function(pd){
  
  #isolate files again and name them with accession number

  t47d_files = as.vector(pd[['t47d']][,'quant_file'])
  names(t47d_files) = as.vector(pd[['t47d']][,'Run'])
  
  mcf7 = tximport(files = mcf7_files, type='salmon', txOut = T, tx2gene = annotLookup)
  t47d = tximport(files = t47d_files, type='salmon', txOut = T, tx2gene = annotLookup)
  
  data = list(mcf7=mcf7, t47d=t47d)
  return(data)
}



#######################################################################################################
# manage directories and file paths
source('global_variables.R')

# read pheno data into R
PD_MCF7 = read.csv(MCF7_RUN_TABLE, sep = '\t')
PD_T47D = read.csv(T47D_RUN_TABLE, sep = '\t')



# get pd with correct path to quant files
pd = get_pheno_data()

# write.csv(ensemble_ids, file = ensembl_id_csv)
ensembl_ids = read.csv(ENSEMBL_ID_CSV, row.names = 1)

# get annotation data
annotLookup = get_annotation(ensembl_ids)

# parse the data into r
data = read_data(pd)

saveRDS(data, file = RDS_FILE )
saveRDS(annotLookup, file = ANNOTATION_DATA_FILE)

