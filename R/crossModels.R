#' Often times it is needed to cross compare edgeR results, limma/voom results across various filtering criteria and comparing normalized results and unnormalized results. This function can run a cross model generation set across various filtering criteria. by default the design is a group-means factorization where the coeficient is the second column.
#' @param crossLevel character option of tx_id, or gene_id which will compare at the transcript or gene level collapse. 
#' @param cutoffMax integer, this will be the maximum read.cutoff that will compare each read.cutoff up to the max, i.e. from 1<=cutoffMax thresholding.
#' @param dataType character either normalized or unnormalized data to compare. if normalized is selected, then ruv is ran to only compare across normalized results.
#' @param outputDir a character path to save all the pdfs printed, includes limma volcano plots, heatmaps.
#' @param numberSelected integer, this is the number of the highest ranked adj.P.Val genes to print into a heatmap, the max amount is the number of genes returned from an analysis.
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
#' @return returns several images plotted.
crossModels<-function(kexp,crossLevel=c("tx_id","gene_id"),cutoffMax=3, dataType=c("normalized","unnormalized"), outputDir=".", design=NULL, p.val=0.05, adjustBy="BH",species=c("Homo.sapiens","Mus.musculus"),numberSelected=200){
  #FIX ME: have this output a single PDF of all the data at the end for each flag, you do not want to print out many PDFs. the pdf should have the EdgeR, limma, GWA/TWA, and RWA, with RUV/without RUV for a fixed cutoffMax.  then you can call this method for various cutoffMaxes.    perhaps do an overlap analysis with a beeSwarm of top Genes, top repeats all of this in 1 PDF. 
  dataType<-match.arg(dataType,c("normalized","unnormalized"))
  crossLevel<-match.arg(crossLevel,c("tx_id","gn_id"))
  adjustBy<-match.arg(adjustBy,c("BH","none","BY","holm"))
  species<-match.arg(species,c("Homo.sapiens","Mus.musculus"))
  
 if(design==NULL){
 stopifnot(design!=NULL)
 }

  kexp<-annotateFeatures(kexp,"transcript")

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
  gwa<-geneWiseAnalysis(kexp,design=design,how="cpm",p.cutoff=p.val,read.cutoff=cutoffMax,species=species,fitOnly=FALSE,adjustBy=adjustBy)
  message("Currently only support normalized analysis")
  } else { #ruv selected
 ruv.gwa<-geneWiseAnalysis(kexp,design=design,how="cpm",p.cutoff=2,read.cutoff=cutoffMax,species=species,fitOnly=TRUE,adjustBy=adjustBy)
  tSel<-ruv.gwa$top
  sili<-tSel[order(tSel$adj.P.Val,decreasing=TRUE),]
  in.sili<-rownames(sili)[1:numberSelected] #take the top numberSelected highest pval
  #run RUV here and recast design matrix
  gwaRuv<-RUVg(round(g.bundles),in.sili,k=1)
  ruvDesign<-cbind(design,gwaRuv$W)
  rownames(ruvDesign)<-colnames(kexp) 
  #plotting pca and RLE
  normalizedCounts<-gwaRuv$normalizedCounts
  #Run GWA with RUV limma
  gwa.RUV<-geneWiseAnalysis(kexp,
                            design=ruvDesign,
                            how="cpm",
                            p.cutoff=p.val,
                            read.cutoff=cutoffMax,
                            species=species,
                            fitOnly=FALSE,
                            adjustBy=adjustBy)
  #gwa RUV heatmap
  write.csv(gwa.RUV,file=paste0(outputDir,"/gwaRUV_pval_",p.val,"_readCutoff_",cutoffMax,".csv"),
  row.names=FALSE)
  write.table(ruvDesign,
              file=paste0(outputDir,"./ruvDesign.txt"),
              quote=FALSE,sep="\t")
  #prepare gwa.RUV into heatmap, by gene name
  gwa.heat<-gwa.RUV$limmaWithMeta
  gwa.heat<-gwa.heat[order(gwa.heat$adj.P.Val,decreasing=FALSE),]
 
  ruv.targets<-rownames(gwa.heat)[1:numberSelected]  
 
  tpm<-collapseTpm(kexp,"gene_id")
  df.kexp<-as.data.frame(tpm,stringsAsFactors=FALSE)
  ruv.idxTpm<-rownames(df.kexp) %in% ruv.targets
  ruv.mt<-as.matrix(df.kexp[ruv.idxTpm,])
    

  for( i in 1:nrow(ruv.mt)){
 rg.idx<-which(rownames(ruv.mt)[i]==rownames(gwa.heat))
   replace<-as.character(gwa.heat$Gene.symbol[rg.idx])
   
    if(!is.na(replace)==TRUE){
    rownames(ruv.mt)[i]<-replace
    }#if gene symbol exists
   else {
   rownames(ruv.mt)[i]<-rownames(gwa.heat)[rg.idx]
   }
 }#for

  #FIX ME: sometimes the cluster doesn't work for shitty data
  x.pv<-pvclust(log(1+ruv.mt),nboot=100)
  dend<-dendsort(hclust(dist(log(1+ruv.mt))),isReverse=TRUE)
  ruv.gwa.heatmap<-Heatmap(log(1+ruv.mt),
                   name="log(1+tpm)",
                   cluster_rows=dend,
                   cluster_columns=x.pv$hclust,
                   column_title=paste0("Limma Normalized P.Val ",p.val," cutoff ",cutoffMax),
                   row_names_gp=gpar(fontsize=6))




  #perform RWA with RUV limma
  rpt.tnx<-collapseTranscripts(kexp,read.cutoff=cutoffMax)
  rpt.idx<-!grepl("^ENS",rownames(rpt.tnx))
  rpts<-rpt.tnx[rpt.idx,]
  rpt.id<- !grepl("^ERCC",rownames(rpts))
  rpts<-rpts[rpt.id,]
  Rpts<-rpts[,colnames(rpts) %in% rownames(design)]
  
  
  ruvRWA<-repeatWiseAnalysis(kexp,
                      design=ruvDesign,
                       how="cpm",
                       p.cutoff=p.val,
                       fold.cutoff=1,
                       read.cutoff=cutoffMax,
                       species=species,
                       adjustBy=adjustBy)
  write.table(ruvRWA$top,
              file=paste0(outputDir,"/NormalizedrepeatWiseAnalysis_pval_",p.val,"_cutoff_",cutoffMax,".txt"),
              quote=FALSE,
              row.names=FALSE,
              col.names=TRUE,
              sep="\t")
   #create a heatmap of normalized repeats
  rpt.tpm<-collapseTpm(kexp,"tx_id")
  df.rpt<-as.data.frame(rpt.tpm,stringsAsFactors=FALSE)
   
  if(nrow(ruvRWA$top)<numberSelected){
    rpt.targets<-rownames(ruvRWA$top)[1:nrow(ruvRWA$top)]
   } else {
    rpt.targets<-rownames(ruvRWA$top)[1:numberSelected]
   }
  df.rpt<-df.rpt[rownames(df.rpt)%in%rpt.targets,]
  rpt.mt<-df.rpt
  rpt.mt<-as.matrix(rpt.mt)

  if(nrow(rpt.mt)>=20){
rpt.pv<-pvclust(log(1+rpt.mt),nboot=50)
rpt.dend<-dendsort(hclust(dist(log(1+rpt.mt))),isReverse=TRUE)


  rh.rpt<-Heatmap(log(1+rpt.mt),
                  name="log(1+tpm)",
                  cluster_rows=rpt.dend,
                  cluster_columns=rpt.pv$hclust,
                  column_title=paste0("RUV Repeats P.Val ",p.val," ",adjustBy," cut ",cutoffMax),
                 row_names_gp=gpar(fontsize=6)) +
  Heatmap(features(kexp)[rownames(rpt.mt)]$tx_biotype,
          name="tx_biotype",
          width=unit(5,"mm")) +
  Heatmap(features(kexp)[rownames(rpt.mt)]$gene_biotype,
          name="rpt_family",
          width=unit(5,"mm"))
  } else {
  rh.rpt<-Heatmap(log(1+rpt.mt),
          name="log(1+tpm)",
          column_title=paste0("RUV Repeats P.Val ",p.val," cut at ",cutoffMax," ",adjustBy),
          row_names_gp=gpar(fontsize=6)) +
  Heatmap(features(kexp)[rownames(rpt.mt)]$tx_biotype,
          name="tx_biotype",
          width=unit(5,"mm")) +
  Heatmap(features(kexp)[rownames(rpt.mt)]$gene_biotype,
          name="rpt_family",
          width=unit(5,"mm"))
  }


 
  #run edgeR GWA and RWA with RUV at gene level
  geneCounts<-collapseBundles(kexp,"gene_id",read.cutoff=cutoffMax)
  d2 <- DGEList(counts=geneCounts, group=phenoDat$type)
  d2 <- calcNormFactors(d2)
  d2 <- estimateGLMTrendedDisp(d2, ruvDesign)
  d2 <- estimateGLMTagwiseDisp(d2, ruvDesign)
  d2 <- estimateGLMRobustDisp(d2, ruvDesign) # superb, but slow!
  #for repeat dispersions
  txCounts<-collapseTranscripts(kexp,read.cutoff=cutoffMax)
  idx<-!grepl("^ENS",rownames(txCounts))
  rpts<-txCounts[idx,]
  id<- !grepl("^ERCC",rownames(rpts))
  rpts<-rpts[id,]
  d <- DGEList(counts=rpts, group=phenoDat$type)
  d <- calcNormFactors(d)
  d <- estimateGLMTrendedDisp(d, ruvDesign)
  d <- estimateGLMTagwiseDisp(d, ruvDesign)
  d <- estimateGLMRobustDisp(d, ruvDesign) # superb, but slow!
  
  gn.fit<-glmFit(d2,ruvDesign)
  gn.lrt<-glmLRT(gn.fit,coef=2)
  gn.Tags<-topTags(gn.lrt,p=p.val,n=nrow(kexp),adjust.method=adjustBy)[[1]]
  gn.feats<-mcols(features(kexp)[features(kexp)$gene_id %in% rownames(gn.Tags)])
  gn.key<-gn.feats[,c(4,5)]
  gn.key<-gn.key[!duplicated(gn.key$gene_id),]

  write.table(gn.Tags,file=paste0(outputDir,"/Normalized.edgeR.pval_",p.val,"_cut_",cutoffMax,"_",adjustBy),quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")
  #edgeR heatmap

  ruv.targets<-rownames(gn.Tags)[1:numberSelected]
  eg.tpm<-collapseTpm(kexp,"gene_id")
  eg.df.kexp<-as.data.frame(eg.tpm,stringsAsFactors=FALSE)
  eg.ruv.idxTpm<-rownames(eg.df.kexp) %in% ruv.targets
  eg.ruv.mt<-as.matrix(eg.df.kexp[eg.ruv.idxTpm,])

  for(i in 1:nrow(eg.ruv.mt)){
   rg.idx<-which(rownames(eg.ruv.mt)[i]==gn.key$gene_id)
   replace<-unique(as.character(gn.key$gene_name[rg.idx]))
    if(!is.na(replace)==TRUE){
      if(any(rownames(eg.ruv.mt)==replace)){
      replace<-paste0(replace,".",i)
      rownames(eg.ruv.mt)[i]<-replace
     } else {
     rownames(eg.ruv.mt)[i]<-replace
     }
    }#if gene symbol exists
  }#for
  eg.x.pv<-pvclust(log(1+eg.ruv.mt),nboot=100)
  eg.dend<-dendsort(hclust(dist(log(1+eg.ruv.mt))),isReverse=TRUE)
  eg.ruv.gwa.heatmap<-Heatmap(log(1+eg.ruv.mt),
                   name="log(1+tpm)",
                   cluster_rows=eg.dend,
                   cluster_columns=eg.x.pv$hclust,
                   column_title=paste0("EdgeR Normalized P.Val ",p.val," cutoff ",cutoffMax," top ",numberSelected),
                   row_names_gp=gpar(fontsize=6))

  ##edgeR rwa analysis csv and heatmap

  

 
  #print only 1 pdf with everything
  pdf(paste0(outputDir,"./NormalizedGeneBundles_readcutoff_",cutoffMax,".pdf"))
  plotRLE(normalizedCounts,outline=FALSE,ylim=c(-2,2))
  title(main="Normalized Gene Bundles CPM")
  plotPCA(normalizedCounts,cex=1.2)
  title("Normalized Gene Bundles PCA")
  volcanoplot(gwa.RUV$fit,coef=2,highlight=20,names=(gwa.RUV$limmaWithMeta$Gene.symbol))
   title("Top Normalized Limma Gene CPM")
  #plot heatmap
   ruv.gwa.heatmap
   rh.rpt
  #include edgeR RUV
   plotBCV(d2)
   title("All Gene RUV Dispersions")
   plotBCV(d)
   title("All Repeat RUV Dispersions")
   eg.ruv.gwa.heatmap
  # include some beeswarm
  dev.off()

  } #normalized
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
 else { #for transcript level
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
 
} ###{{{ 
