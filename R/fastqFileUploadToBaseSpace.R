#' once the sraFastqHeaderToIlluminaStandard has been ran on SRA imported samples, the basespace fileupload can be ran. this assumes that a user has a basespace account and has a project directory listed correctly.  We assume that the basespace CLI is installed please see (https://help.basespace.illumina.com/articles/descriptive/basespace-cli/).  this R interface to BaseSpace cloud system is mostly useful for numerous samples targeting the uses of single cell sequencing with samples ranging from upwards 800; thus making fastq uploading difficult; hence the automation of it within Arkas. Fastq Headers must be in illumina standards. FIX ME: we would also like to add the execution of arkas via the basespace CLI to run cloud applications via R.
#' @param illuminaDirPath  a illumina directory path to illumina standard fastq files with illumina header and naming conventions.
#' @param illuminafastqFile  fastq files with illumina headers and naming conventions, a vector of file names in illumina standard, this is optional parameter for large directories where hte names are numerious, in this case use the file signature for multi uploads
#' @param basespaceProject   character string of the basespace project name, this must exist on basespace 
#' @param fastqFileSignature character that is unique to the fastq file directory where upon grep'ing the desired files will get matched.  the default is the illumina standard suffix _001.fastq.gz which should pick out the illumina files in the case where the fastq directory has multiple raw files.
#' @param illuminaDirs a character of illumina sample directories which contain illumina fastqs, if SE, only 1 fastq in this dir, if PE, there should be 2 or more.
#' @param paired boolean   if the illuminaDirs contains a paired fastq pair
#' @param dryRun boolena   if true then for bs CLI will simulate an upload, if false bs CLI will upload a fastq to the real account.
#' @return nothing, a successful indication that the files were uploaded to basespace
#' @examples 
#' \dontrun{
#' fastqFileUploadToBaseSpace(illuminaDirPath,
#'                           illuminafastqFile=illuminaFastqFiles,
#'                           basespaceProject=basespaceProject,
#'                           illuminaDirs=illuminaDirs)
#' } 
#' @return integer for success or failure
#' 
#' @export
fastqFileUploadToBaseSpace<-function(illuminaDirPath=NULL, illuminafastqFile=NULL, basespaceProject=NULL,fastqFileSignature="_001.fastq.gz",illuminaDirs=NULL,paired=FALSE, dryRun=FALSE) { # {{{

  # FIXME: if a directory of fastqs exist with multiple R1, and multiple R2, there may be errors.  This was designed for 1 R1 and 1 R2 (if PE)
  # dependencies: https://help.basespace.illumina.com/articles/descriptive/basespace-cli/ must be installed
    dirs<-paste0(illuminaDirPath,"/",illuminaDirs,"/")
    illuminafastqFile<-dir(dirs)[grepl(fastqFileSignature,dir(dirs))]
    x0<-paste0(dirs,"/",dir(dirs)[grepl("R1",dir(dirs))]) #a vector , 2 entries if PE
    stopifnot(file.exists(x0)==TRUE)
    if(paired==TRUE){
      for(j in 2:length(illuminafastqFile)){
        x1<-paste0(dirs,"/",illuminafastqFile[j])
        stopifnot(file.exists(x1)==TRUE)
        x0<-paste(x0,x1)
      } # looping over illuminafastqFile vector, should be 2 here if PE==TRUE
    }
    if(dryRun==FALSE){
      command<-paste0("bs upload sample ",x0," -p ",basespaceProject)
    } else { 
      command<-paste0("bs upload sample --dry-run ",x0," -p ",basespaceProject)
    } #simulated upload for unit testing
    
    system(command) 
  
}# }}} main

