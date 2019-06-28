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

WORKING_DIRECTORY = '/media/ncw135/DATA/PublicRepositories/geodatasets'
DATA_DIR = file.path(WORKING_DIRECTORY, 'data/GSE107864')
QUANT_DIR = file.path(DATA_DIR, 'quants')
MCF7_QUANTS = file.path(QUANT_DIR, 'MCF7s')
T47D_QUANTS = file.path(QUANT_DIR, 'T47Ds')

MCF7_RUN_TABLE = file.path(DATA_DIR, 'SraRunTableMCF7.txt')
T47D_RUN_TABLE = file.path(DATA_DIR, 'SraRunTableT47D.txt')

# directories for storing some data
SAVED_OBJECTS = file.path(WORKING_DIRECTORY, 'SavedObjects')
RDS_FILE = file.path(SAVED_OBJECTS, 'POST_SALMON_DATA.rds')
ANNOTATION_DATA_FILE = file.path(SAVED_OBJECTS, 'ANNOTATION_DATA.rds')

# get ensemble ids for input into annotation function
ENSEMBL_ID_CSV = file.path(SAVED_OBJECTS, 'ensembl_ids.csv')

# will throw error if any of these files/folders do not exist
validate_directories()