#' this script prepares fastq files downloaded from SRADb R package or files downloaded from European Nucleic Acids ftp fastq downloads, which have sra fastq headers, and changes the sra headers to illumina standard headers as a preparatory step for uploading fastq files into the  cloud system,into an open account and to a project under an authorized user.  We assume that the user can import files from SRA or ENA, and use fastq-dump to convert SRA files into fastq files.  arkas does support fastq header conversions into BaseSpace standard nomenclature. this method will check the header of the fastq, if it is an sra header, then arkas will transform the header to illumina standard; further this script will also rename the fastq file itself to illumina standard, and prepare the fastq file for an upload to the basespace cloud using the script fastqFileUploadToBaseSpace.R , inputs must be gzipped.
#' @param headerStyle ,  character string, where SRA is the header style obtained from importing from SRAdb; the Normal header format is a standard two field record delimited by a space, with each field delimited by several colons.
#' @param fastqPath,  character string path to the fastq directory
#' @param fastqFile a vector of sample files downloaded from SRA to convert to illumina std
#' @param sraOutputDir character vector where to spit out the converted fastq
#' @param fastqReadNumber  integer 1, or 2, this is the read number that will be written in the fastq file, for SRAdb defaults to 1, but we give the option to have 1 or 2 here.
#' @param hasBarcode   boolean , some SRA files have a barcode and some sra do not have a barcode, this flag specifies which case.  if the raw SRA fastq has a barcode (run ID) it looks like @SRR1564893.1 HWI-ST972:1180:D225DACXX:7:1101:1247:2104 length=50, if the run ID is missing the raw SRA header looks like @SRR892995.1 HWI-ST601:8:1101:1219:2076 length=100.  if the run ID is missing, we make a fake runIDNumber along with a runID tag as the 2nd and 3rd column.
#' @importFrom nat.utils is.gzip
#' @export
#' @return a integer defining success or failure
sraFastqHeaderToIlluminaStandard<-function(headerFormat=c("SRA","Normal"), fastqPath,fastqFile, sraOutputDir, fastqReadNumber=1,hasBarcode=FALSE  ) {
#this function is not exported, and sraFastqHeaderToIlluminaStandard handles
hFormatted <- match.arg(headerFormat, c("SRA", "Normal"))


# requirements: the fastq file must be gzipped
stopifnot(file.exists(paste0(fastqPath,"/",fastqFile))=="TRUE")

stopifnot( is.gzip(paste0(fastqPath,"/",fastqFile))=="TRUE")


message("files appear gzipped .. O.K!, checking format")



if(hFormatted =="SRA") {
#CASE 1:  SRA headers  @SRR3173882.sra.1 HWI-ST1209-LAB:323:HA9TPADXX:1:1101:1408:2086 length=50
fastqHeaderConvert<-system.file("bin","fastqHeaderConvert.sh",package="arkas")
fastqConvert2<-system.file("bin","fastqConvert2.sh",package="arkas")
convertFastqToIlluminaStandard<-system.file("bin","convertFastqFileToIlluminaStd.sh",package="arkas")
sraSingleUpload<-system.file("bin","sraSingleFastqBaseSpaceUpload.sh", package="arkas")

setwd(fastqPath)
    if(hasBarcode==TRUE){
message("running SRA fastq header conversion e.g.  @SRR3173882.sra.1 HWI-ST1209-LAB:323:HA9TPADXX:1:1101:1408:2086 length=50...")
fastqConvert3<-system.file("bin","fastqConvert3.sh",package="arkas")
command<-paste0(fastqHeaderConvert," ",fastqFile," ", fastqReadNumber," | ",fastqConvert2," | ", fastqConvert3," | gzip -c > ",sraOutputDir,"/",fastqFile)
system(command)
  } #has barcode like @SRR1564893.1 HWI-ST972:1180:D225DACXX:7:1101:1247:2104 length=50

   if(hasBarcode==FALSE){
message("running SRA fastq header conversion e.g.  @SRR3173882.sra.1 HWI-ST1209-LAB:323:1:1101:1408:2086 length=50...")

fastqConvert3<-system.file("bin","fastqConvert_NoRunId.sh",package="arkas")
command<-paste0(fastqHeaderConvert," ",fastqFile," ", fastqReadNumber," | ",fastqConvert2," | ", fastqConvert3," | gzip -c > ",sraOutputDir,"/",fastqFile)
system(command)
}

#now to change the name to illumina std and upload to basespace
command2<-paste0(convertFastqToIlluminaStandard," ",paste0(sraOutputDir,"/"),fastqFile)
system(command2) 
message(paste0("converted ",fastqFile," ready to upload to illumina basespace")) 
}

if(hFormatted =="Normal") {
#CASE 2:  Normal Fastq Headers @HWI-ST1209:323:HA9TPADXX:1:1101:1408:2086 1:N:0:1
message("running upload for normal headers e.g.  @HWI-ST1209:323:HA9TPADXX:1:1101:1408:2086 1:N:0:1 ..")

}


} #{{{ main
