library(GEOquery)

WD = getwd()
DATA_DIR = file.path(WD, 'data')
DATA_DIR

accessions.file = file.path(DATA_DIR, 'accessions')
# These two experiments were conducted in https://www.cell.com/cell-reports/pdf/S2211-1247(18)31714-5.pdf
# time series  
acc1 = 'GSE107864'
# not time series
acc2 = 'GSE119455'

accessions = read.csv(accessions.file)

# get a single GDS by accession number
get1GDS = function(accession_number){
  if (length(accession_number) != 1){
    stop('dont be a dick')
  }
    if (class(accession_number) == 'factor'){
      accession_number = as.character(accession_number)
    }
    # get path for storing the data and create
    path = file.path(DATA_DIR, accession_number)
    dir.create(path, showWarnings = T,  recursive = T)
    print(path)
    print(dir.exists(path))
    # if already exists, use this
    existing_path = file.path(path, paste0(accession_number, '.soft.gz'))
    # get data
    if (file.exists(existing_path)){
      print('exist')
      gds = getGEO(accession_number, destdir = path, filename = existing_path,
                   AnnotGPL = T, getGPL = T)
    } else{
      print('not exist')
      gds = getGEO(accession_number, destdir = path,
                   AnnotGPL = T, getGPL = T)
    }
    
    
    return(gds[[1]])
}

getListOfGDS = function(accessions){
  l = list()
  for (i in 1:dim(accessions)[1]){
    message(paste('getting accesion number', i))
    gds = get1GDS(accession_number = accessions[i, 1])
    l[[accessions[i, 1]]] = gds
  }
  return (l)
}

# gds = get1GDS(accessions[1, 1])
# 
# gds.list = getListOfGDS(accessions)
# 
# names(gds.list)
# Meta(gds)



# gds = getGEO(acc1, destdir = file.path(DATA_DIR, acc1), GSEMatrix = T)
# gds = getGEO(acc1)
# gds = get1GDS(acc1)
# gds

# getGEOSuppFiles(acc1, baseDir = file.path(DATA_DIR, acc1))


# getGEO(filename = file.path(DATA_DIR, paste0(acc1, '/', acc1, '.soft.gz')))
# getGEO(filename = file.path(DATA_DIR, paste0(acc1, '/', acc1, '.soft.gz')))
mcf7_acc = 'GSE107858' 
t47d_acc = 'GSE107863'
# getGEO(mcf7_acc)
mcf7 = get1GDS(mcf7_acc)
mcf7 

pd = phenoData(mcf7)

varLabels(pd)


class(mcf7)

names(mcf7)

eset = exprs(mcf7)
eset

?exprs
class(eset)
dim(eset)

'SRP125153'


