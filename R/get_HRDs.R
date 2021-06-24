# New HRD scar score

GetGzFromUrl <- function(url, ...) {
  # http://stackoverflow.com/questions/9548630/read-gzipped-csv-directly-from-a-url-in-r
  con <- gzcon(url(url))
  txt <- readLines(con)
  dat <- read.delim(textConnection(txt), ...)
  return(dat)
}

GetChrominfo <- function() {
  # Get chromInfo table from UCSC
  #chrom <- GetGzFromUrl("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz", header = FALSE)
  chrom <- as.data.frame(read.delim("/home/vijay/my_projects/ovaHRDscar/chromInfo.txt"))
  gaps <- as.data.frame(read.delim("/home/vijay/my_projects/ovaHRDscar/gap.txt"))
  chrom <- subset(chrom, grepl("^chr[0-9XY]{1,2}$", chrom[,1]))
  # Get gap table from UCSC
  #gaps <- GetGzFromUrl("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz", header = FALSE)
  centro <- subset(gaps, gaps[,8] == "centromere")
  # Merge the relevant info from the two tables
  chrominfo <- merge(chrom[,1:2], centro[,2:4], by.x = 1, by.y = 1) # merge size and centromere location
  chrominfo$centromere <- rowMeans(chrominfo[,3:4]) # convert centromere start and end into one location (the mean)
  chrominfo <- chrominfo[,c(1,2,5,3,4)] # keep chromosome, size and centromere location
  colnames(chrominfo) <- c("chr", "size", "centromere", "centstart", "centend")
  chrominfo[,1] <- as.character(chrominfo[,1])
  chrominfo$chr <- sub("chr", "", chrominfo$chr)
  chrominfo$chr <- sub("X", "23", chrominfo$chr)
  chrominfo$chr <- sub("Y", "24", chrominfo$chr)
  chrominfo[,1] <- as.numeric(chrominfo[,1])
  chrominfo <- chrominfo[order(chrominfo$chr), ]  # order by chromosome number
  rownames(chrominfo) <- as.character(chrominfo[,1])
  chrominfo <- as.matrix(chrominfo)
  return(chrominfo)
}

 # hg19

get.ovaHRDscars <- function(seg, chrominfo = grch37, LOH_windos=c(10,50), LST_segSize=12e6, LST_mindistance=1e6){
  seg <- preparing.input(seg)
  
  #if (chrominfo == "grch38"){
  #  chrom = chrominfo_grch38
  #}
  #if (chrominfo == "grch37"){
  #  chrom = chrominfo_grch37
  #}
  grch37 <- GetChrominfo()
  chrom = grch37

  #Calculating HRD-LOH
  HRD_LOHs <- features.LOH(seg, MbSizes=LOH_windos)
  HRD_LOHs <- HRD_LOHs[,1]

  #Calculating LSTs
  LSTs <- LSTs(seg, chrominfo= chrom, segsizes=LST_segSize, mindistance=LST_mindistance)

  #Calculating nTAIs
  res_ai<- calc.ai_new(seg,chrom)
  nTAIs <- res_ai[,1]

  HRDScore <- HRD_LOHs + LSTs + nTAIs

  #Concatenating results

  HRDresults <- cbind(HRD_LOHs, LSTs, nTAIs, HRDScore)
  colnames(HRDresults) <- cbind("nLOH","LSTs","nTAIs", "HRDScore")
  assign("HRDresults",as.data.frame(HRDresults),envir = .GlobalEnv)
  return(HRDresults)
}
