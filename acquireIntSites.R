## check for presence of R packages 
rPackages <- c("stats", "methods", "RMySQL") 
stopifnot(all(sapply(rPackages, require, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE))) 

options(stringsAsFactors=F) 
#' increase output width to console width 
wideScreen <- function(howWide=as.numeric(strsplit(system('stty size', intern=T), ' ')[[1]])[2]) { 
    options(width=as.integer(howWide)) 
 } 
 wideScreen() 

## get argument 
args <- commandArgs(trailingOnly=TRUE) 
gtspid <- args[1] 
##gtspid <- "GTSP0568" 

if( is.na(gtspid) ) { 
   message("Usage:\n\tRscript acquireSites.R GTSP####")  
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

## get all replicates 
replicates <- subset(allSampleName, grepl(paste0("^",gtspid), sampleName) ) 
stopifnot( nrow(replicates)>0 ) 

## get sampleInfo from specimen_management table
sql <- sprintf("SELECT * FROM specimen_management.gtsp WHERE SpecimenAccNum = '%s'", gtspid) 
sampleInfo <- suppressWarnings( dbGetQuery(dbConn, sql) ) 
stopifnot( nrow(sampleInfo)==1) 
colnames(sampleInfo) <- tolower(colnames(sampleInfo)) 

sampleIDin <- sprintf("(%s)", paste(unique(replicates$sampleID), collapse=",")) 

##get unique sites 
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

##output bed file 
needed <- c("chr", "position", "breakpoint",
            "sampleName", "refGenome", "strand") 

bed <- rbind(subset(sites.uniq, select=needed), 
             subset(sites.multi, select=needed)) 

bed <- plyr::arrange(bed, chr, position, strand) 

##bed$name <- gtspid 
bed$name <- paste(sampleInfo$patient, sampleInfo$timepoint, sampleInfo$celltype, sep=":") 
bed$score <- 400 

bed <- plyr::arrange(bed, chr, position, strand) 

bed <- data.frame(chrom=as.character(bed$chr), 
                  chromStart=as.integer(pmin(bed$position, bed$breakpoint)), 
                  chromEnd=as.integer(pmax(bed$position, bed$breakpoint)), 
                  name=as.character(bed$name), 
                  score=500, 
                  strand=as.character(bed$strand)) 

fileName <- paste(gtspid, sampleInfo$patient, sampleInfo$timepoint, sampleInfo$celltype, "bed", sep=".") 
fileName <- paste(gtspid, "bed", sep=".") 

write.table(bed, file=fileName, 
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE) 
