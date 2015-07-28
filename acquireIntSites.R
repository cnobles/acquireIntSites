## check for presence of R packages 
rPackages <- c("stats", "methods", "RMySQL") 
stopifnot(all(sapply(rPackages, require, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE))) 

options(stringsAsFactors=F) 
#' increase output width to console width 
#wideScreen <- function(howWide=as.numeric(strsplit(system('stty size', intern=T), ' ')[[1]])[2]) { 
#    options(width=as.integer(howWide)) 
# } 
# wideScreen() 

## get argument 
args <- commandArgs(trailingOnly=TRUE) 
gtspid <- args[1] 
query_by <- "GTSP"
sample_query <- gtspid

if(grepl(".csv", gtspid)){
  queryInfo <- read.csv(file = gtspid)
  if("patient" %in% colnames(queryInfo) & "GTSP" %in% colnames(queryInfo)){
    stop("Both patient and GTSP columns found in ", gtspid,".\n\tPlease only use one.")
  }else if("patient" %in% colnames(queryInfo)){
    sample_query <- unique(queryInfo$patient)
    
    junk <- sapply(dbListConnections(MySQL()), dbDisconnect) 
    dbConn <- dbConnect(MySQL(), group="specimen_management")  
    stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)
    
    query_selection <- "SELECT Patient,SpecimenAccNum FROM specimen_management.gtsp "
    string <- paste(sample_query, collapse="' OR Patient = '")
    query_request <- paste0("WHERE Patient = '", string, "'")
    query <- paste0(query_selection, query_request)
    patient_GTSP <- dbGetQuery(dbConn, query)
    
    dbDisconnect(dbConn)
    rm(string, query, query_request, query_selection, dbConn, junk)
    
    gtspid <- unique(patient_GTSP$SpecimenAccNum)
    query_by <- "patient"
  }else if("GTSP" %in% colnames(queryInfo)){
    sample_query <- unique(queryInfo$GTSP)
    gtspid <- unique(queryInfo$GTSP)
    query_by <- "GTSP"
  }else{
    stop("Neither patient or GTSP column found in csv file.")
  }
}else if(!grepl("GTSP", gtspid)){
  sample_query <- gtspid
  
  junk <- sapply(dbListConnections(MySQL()), dbDisconnect) 
  dbConn <- dbConnect(MySQL(), group="specimen_management")  
  stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)
  
  query_selection <- "SELECT Patient,SpecimenAccNum FROM specimen_management.gtsp "
  string <- paste(sample_query, collapse="' OR Patient = '")
  query_request <- paste0("WHERE Patient = '", string, "'")
  query <- paste0(query_selection, query_request)
  patient_GTSP <- dbGetQuery(dbConn, query)
  
  dbDisconnect(dbConn)
  rm(string, query, query_request, query_selection, dbConn, junk)
  
  gtspid <- unique(patient_GTSP$SpecimenAccNum)
  query_by <- "patient"
  
}

if( is.na(gtspid) ) { 
   message("Usage:\n\tRscript acquireSites.R GTSP????\n\tor\n\tRscript acquireSites.R *.csv")
   q(status=1) 
} 


## check if file exist and permission .my.cnf  
stopifnot(file.exists("~/.my.cnf")) 
stopifnot(file.info("~/.my.cnf")$mode == as.octmode("600")) 

## initialize connection to database 
## ~/.my.cnf must be present 
junk <- sapply(dbListConnections(MySQL()), dbDisconnect) 
dbConn <- dbConnect(MySQL(), group="intSitesDev237")  
stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)

allSampleName <- suppressWarnings( dbGetQuery(dbConn, "SELECT * FROM samples") ) 
dbDisconnect(dbConn)
rm(dbConn, junk)

## get all replicates 
rep_string <- do.call(c, lapply(gtspid, function(x){
  grep(x, allSampleName$sampleName, value = TRUE)
  }))
replicates <- allSampleName[match(rep_string, allSampleName$sampleName),] 
stopifnot( nrow(replicates)>0 ) 

## get sampleInfo from specimen_management table
junk <- sapply(dbListConnections(MySQL()), dbDisconnect) 
dbConn <- dbConnect(MySQL(), group="specimen_management")  
stopifnot(dbGetQuery(dbConn, "SELECT 1")==1) 

string <- paste(gtspid, collapse="' OR SpecimenAccNum = '")
query_request <- paste0("WHERE SpecimenAccNum = '", string, "'")
sql <- paste0("SELECT * FROM specimen_management.gtsp ", query_request) 
sampleInfo <- suppressWarnings( dbGetQuery(dbConn, sql) ) 
dbDisconnect(dbConn) #Disconnect from specimen_management after query
rm(junk, dbConn, string, query_request) #Cleaning up so there is no interference or mistakes with the second query

if(query_by == "patient"){
  isThere <- unique(sampleInfo$Patient)
}else if(query_by == "GTSP"){
  isThere <- unique(sampleInfo$SpecimenAccNum)
}
stopifnot( (NA %in% match(sample_query, isThere)) == FALSE) 
colnames(sampleInfo) <- tolower(colnames(sampleInfo)) 

sampleIDin <- sprintf("(%s)", paste(unique(replicates$sampleID), collapse=",")) 

##get unique sites 
junk <- sapply(dbListConnections(MySQL()), dbDisconnect) 
dbConn <- dbConnect(MySQL(), group="intSitesDev237")  
stopifnot(dbGetQuery(dbConn, "SELECT 1")==1) 

sql <- paste("SELECT DISTINCT * 
             FROM samples JOIN sites 
             ON samples.sampleID = sites.sampleID 
             JOIN pcrbreakpoints 
             ON pcrbreakpoints.siteID = sites.siteID  
             WHERE samples.sampleID in ", sampleIDin ) 
sites.uniq <- suppressWarnings( dbGetQuery(dbConn, sql) )  
sites.uniq <- sites.uniq[, !duplicated(colnames(sites.uniq))] 

##get multihit sites 
sql <- paste("SELECT * 
             FROM samples JOIN multihitpositions 
             ON samples.sampleID = multihitpositions.sampleID 
             JOIN multihitlengths 
             ON multihitpositions.multihitID = multihitlengths.multihitID ", 
             "WHERE samples.sampleID in ",  sampleIDin ) 
sites.multi <- suppressWarnings( dbGetQuery(dbConn, sql) )  
sites.multi <- sites.multi[, !duplicated(colnames(sites.multi))] 
sites.multi$breakpoint <- ifelse(sites.multi$strand=="+", 
sites.multi$position+sites.multi$length-1, 
sites.multi$position-sites.multi$length+1)

dbDisconnect(dbConn)
rm(dbConn, junk)

#Save data of interest
save(sampleInfo, file = "sampleInfo.RData")
save(sites.uniq, file = "sites.uniq.RData")
save(sites.multi, file = "sites.multi.RData")