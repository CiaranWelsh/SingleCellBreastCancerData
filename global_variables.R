# this file contains variables shared among all scripts 

validate_directories = function(){
  # check that directories exist
  for (i in c(WORKING_DIRECTORY, DATA_DIR, QUANT_DIR, MCF7_QUANTS, T47D_QUANTS)){
    assertthat::is.dir(i)
  }
  
  # check that required files exist
  for (i in c(MCF7_RUN_TABLE, T47D_RUN_TABLE)){
    assertthat::is.readable(i)
  }
}

WORKING_DIRECTORY = '/media/ncw135/DATA/SingleCellBreastCancer'
DATA_DIR = file.path(WORKING_DIRECTORY, 'data')
QUANT_DIR = file.path(DATA_DIR, 'quants')
MCF7_QUANTS = file.path(QUANT_DIR, 'MCF7s')
T47D_QUANTS = file.path(QUANT_DIR, 'T47Ds')

MCF7_RUN_TABLE = file.path(DATA_DIR, 'SraRunTableMCF7.txt')
T47D_RUN_TABLE = file.path(DATA_DIR, 'SraRunTableT47D.txt')

# directories for storing some data
SAVED_OBJECTS = file.path(WORKING_DIRECTORY, 'SavedObjects')
RDS_FILE = file.path(SAVED_OBJECTS, 'POST_SALMON_DATA.rds')
MCF7_ANNOTATION_DATA_FILE = file.path(SAVED_OBJECTS, 'MCF7_ANNOTATION_DATA.rds')
T47D_ANNOTATION_DATA_FILE = file.path(SAVED_OBJECTS, 'T47D_ANNOTATION_DATA.rds')

# get ensemble ids for input into annotation function
ENSEMBL_ID_CSV = file.path(SAVED_OBJECTS, 'ensembl_ids.csv')

# list of bulk data IDS
MCF7_BULK_DATA_IDS = c('SRR6301063', 'SRR6301064', 
                       'SRR6301065', 'SRR6301066', 
                       'SRR6301067', 'SRR6301068', 
                       'SRR6301069', 'SRR6301070')
T47D_BULK_DATA_IDS = c('SRR6301155', 'SRR6301201', 
                       'SRR6301212', 'SRR6301126')

# will throw error if any of these files/folders do not exist
validate_directories()