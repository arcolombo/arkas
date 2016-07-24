#' Often times it is needed to cross compare edgeR results, limma/voom results across various filtering criteria and comparing normalized results and unnormalized results, which is not a good idea.
#' @param crossLevel character option of tx_id, or gene_id which will compare at the transcript or gene level collapse. 
#' @param cutoffMax integer, this will be the maximum read.cutoff that will compare each read.cutoff up to the max, i.e. from 1<=cutoffMax thresholding.
#' @param dataType character either normalized or unnormalized data to compare. if normalized is selected, then ruv is ran to only compare across normalized results.
#' @param outputDir a character path to save all the pdfs printed, includes limma volcano plots, heatmaps.
#' @param design a matrix specifying contrasts to compare.
#' @param p.val  numeric pvalue used for linear fitting
#' @param adjustBy character for FDR filters
#' @import edgeR
#' @import limma
#' @import ComplexHeatmap
#' @import circlize
#' @import pvclust
#' @import dendsort
#' @import EDASeq
#' @export
#' @return returns a data frame and several images to plot.
crossModels<-function(kexp,crossLevel=c("tx_id","gene_id"),cutoffMax=3, dataType=c("normalized","unnormalized"), outputDir=".", design=NULL, p.val=0.05, adjustBy="BH",species=c("Homo.sapiens","Mus.musculus")){
  #FIX ME: have this output a single PDF of all the data at the end for each flag, you do not want to print out many PDFs. the pdf should have the EdgeR, limma, GWA/TWA, and RWA, with RUV/without RUV for a fixed cutoffMax.  then you can call this method for various cutoffMaxes.    perhaps do an overlap analysis with a beeSwarm of top Genes, top repeats all of this in 1 PDF. 
  dataType<-match.arg(dataType,c("normalized","unnormalized"))
  crossLevel<-match.arg(crossLevel,c("tx_id","gn_id"))
  adjustBy<-match.arg(adjustBy,c("BH","none","BY","holm"))
  species<-match.arg(species,c("Homo.sapiens","Mus.musculus"))
  
 if(design==NULL){
 stopifnot(design!=NULL)
 }

  annoKexp<-annotateFeatures(kexp,"transcript")

 if(crossLevel=="gene_id"){
 #run collapse bundles by gene_id for an incrementation to cutoffMatx
 g.bundles<-collapseBundles(kexp,"gene_id",read.cutoff=cutoffMax)
 #print out the RLE for gene, and repeat only regions, and print a PCA of the whole data set, no RUV here.
 if(dataType=="unnormalized"){
 pdf(paste0(outputDir,"/GeneBundles_readCutoff_",cutoffMax,".pdf"))
 plotRLE(g.bundles,outline=FALSE)
 title("Gene Bundles Raw CPM")
 plotPCA(g.bundles,cex=1.2)
 title("PCA Gene Bundles Raw CPM")
 dev.off()
 }
 else { #ruv selected
 ruv.gwa<-geneWiseAnalysis(kexp,design=design,how="cpm",p.cutoff=2,read.cutoff=cutoffMax,species=species,fitOnly=TRUE,adjustBy=adjustBy)
  tSel<-ruv.gwa$top
  sili<-tSel[order(tSel$adj.P.Val,decreasing=TRUE),]
  in.sili<-rownames(sili[1:100]) #take the top 100 highest pval
  #run RUV here and recast design matrix
  gwaRuv<-RUVg(round(g.bundles),in.sili,k=1)
  ruvDesign<-cbind(design,gwaRuv$W)
  rownames(ruvDesign)<-colnames(kexp) 
  #plotting pca and RLE
  normalizedCounts<-gwaRuv$normalizedCounts
  pdf(paste0(outputDir,"./NormalizedGeneBundles_readcutoff_",cutoffMax,".pdf"))
  plotRLE(normalizedCounts,outline=FALSE,ylim=c(-2,2))
  title(main="Normalized Gene Bundles CPM")
  plotPCA(normalizedCounts,cex=1.2)
  title("Normalized Gene Bundles PCA")
  dev.off()
  #Run GWA with RUV limma
  #print pdf volcano plots of repeats
  #perform RWA with RUV limma
  #
  #run edgeR GWA and RWA with RUV
  }
 } #gene id 

 else {
 # run collapse bundles by tx_id for an incrementation to cutoffMax.
 t.bundles<-collapseBundles(kexp,"tx_id",read.cutoff=cutoffMax)
 #print out RLE of gene and repeat regions, print out a PCA plot, no RUV here. only CPM data
 if(dataType=="unnormalized"){
 pdf(paste0(outputDir,"/TranscriptBundles_readCutoff_",cutoffMax,".pdf"))
 plotRLE(t.bundles,outline=FALSE)
 title("Transcript Bundles Raw CPM")
 plotPCA(t.bundles,cex=1.2)
 title("PCA Transcript Bundles Raw CPM")
 dev.off()
 }
 else {
 ruv.twa<-transcriptWiseAnalysis(kexp,design=design,p.cutoff=2,read.cutoff=cutoffMax,species=species,adjustBy=adjustBy)
  tSel<-ruv.twa$top
  sili<-rownames(tSel[order(tSel$adj.P.Val,decreasing=TRUE),])
  in.sili<-sili[1:100] #take the top 100 highest pval
  #run RUV here recreate design matrix
  twaRuv<-RUVg(round(t.bundles),in.sili,k=1) 
  ruvDesign<-cbind(design,twaRuv$W)
  rownames(ruvDesign)<-colnames(kexp)
  #plotting pca and RLE
  normalizedCounts<-twaRuv$normalizedCounts
  pdf(paste0(outputDir,"./NormalizedTranscriptBundles_readcutoff_",cutoffMax,".pdf"))
  plotRLE(normalizedCounts,outline=FALSE,ylim=c(-2,2))
  title(main="Normalized Transcript Bundles CPM")
  plotPCA(normalizedCounts,cex=1.2)
  title("Normalized Transcript Bundles PCA")
  dev.off()
  #run TWA with RUV 
  } 

 }
 
 #output Heatmaps, volcano plots for each incrementations


} ###{{{ 
