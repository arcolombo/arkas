#' Often times it is needed to cross compare edgeR results, limma/voom results across various filtering criteria and comparing normalized results and unnormalized results, comparing gene-level, and / or transcript level at various filtering criteria reproducibly. This function can run a cross model generation set across various filtering criteria. by default the design is a treatment factorization where the coeficient is the second column with an intercept term.
#' @param crossLevel character option of tx_id, or gene_id which will compare at the transcript or gene level collapse. 
#' @param cutoffMax integer, this will be the maximum read.cutoff that will compare each read.cutoff up to the max, i.e. from 1<=cutoffMax thresholding.
#' @param dataType character either normalized or unnormalized data to compare. if normalized is selected, then ruv is ran to only compare across normalized results.
#' @param outputDir a character path to save all the pdfs printed, includes limma volcano plots, heatmaps.
#' @param design  a matrix with a treatment level contrasts, does not yet support group-means factorization
#' @param setP numeric for linear fitting
#' @param adjustBy character either BH,none,BY,holm
#' @param species character Homo.sapiens or Mus.musculus
#' @param numberSelected integer, this is the number of the highest ranked adj.P.Val genes to print into a heatmap, the max amount is the number of genes returned from an analysis.
#' @param saveReport boolean, if true then a txt and csv files are printed out to file, if false, then no report is printed out
#' @import edgeR
#' @import limma
#' @import ComplexHeatmap 
#' @import circlize
#' @importFrom pvclust pvclust
#' @importFrom dendsort dendsort
#' @importFrom EDASeq plotRLE
#' @importFrom EDASeq plotPCA
#' @importFrom RUVSeq RUVg
#' @importFrom grid gpar
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics title
#' @importFrom utils data read.delim write.csv write.table
#' @export
#' @return returns several images plotted.
crossModels<-function(kexp,
                      crossLevel=c("tx_id","gene_id"),
                      cutoffMax=3, 
                      dataType=c("normalized","unnormalized"),
                      outputDir=".",
                      design=NULL, 
                      setP=0.05, 
                      adjustBy="BH",
                      species=c("Homo.sapiens","Mus.musculus"),
                      numberSelected=200,
                      saveReport=FALSE){
   
 readkey<-function()
{
    cat ("Press [enter] to continue")
    line <- readline()
}

bySd <- function(x, k=500) { 
   sds<-vector() 
  for(i in 1:nrow(x)){
   sds[i]<-sd(x[i,],na.rm=TRUE)
   }
  x[rev(order(sds))[seq_len(k)],]
}



#FIX ME: add a flag for doing a multiple group contrast. the user must input a good working contrast.matrix, then the output if flagged would go through each coefficient.

  dataType<-match.arg(dataType,c("normalized","unnormalized"))
  crossLevel<-match.arg(crossLevel,c("tx_id","gene_id"))
  adjustBy<-match.arg(adjustBy,c("BH","none","BY","holm"))
  species<-match.arg(species,c("Homo.sapiens","Mus.musculus"))
  
 if(is.null(design)==TRUE){
 stop("design must be input")
 }

  kexp<-annotateFeatures(kexp,"transcript")

 if(crossLevel=="gene_id"){
 #run collapse bundles by gene_id for an incrementation to cutoffMatx
 g.bundles<-collapseBundles(kexp,"gene_id",read.cutoff=cutoffMax)
 #print out the RLE for gene, and repeat only regions, and print a PCA of the whole data set, no RUV here.
 if(dataType=="unnormalized"){
 print(plotRLE(g.bundles,outline=FALSE))
 title("Gene Bundles Raw CPM")
 readkey()
 print(plotPCA(g.bundles,cex=1.2))
 title("PCA Gene Bundles Raw CPM")
 readkey()
 
  gwa<-geneWiseAnalysis(kexp,design=design,how="cpm",p.cutoff=setP,read.cutoff=cutoffMax,species=species,fitOnly=FALSE,adjustBy=adjustBy)
   #gwa NO RUV heatmap
  if(saveReport==TRUE){
  write.csv(gwa$limmaWithMeta[order(gwa$limmaWithMeta$adj.P.Val),],file=paste0(outputDir,"/gwaLimmaTMM_pval_",setP,"_readCutoff_",cutoffMax,".csv"),
  row.names=FALSE)
  }#print to file
  
  #prepare gwa into heatmap, by gene name
  gwa.heat<-gwa$limmaWithMeta
  gwa.heat<-gwa.heat[order(gwa.heat$adj.P.Val,decreasing=FALSE),]
  targets<-rownames(gwa.heat)[1:numberSelected]

  tpm<-collapseTpm(kexp,"gene_id")
  df.kexp<-as.data.frame(tpm,stringsAsFactors=FALSE)
  idxTpm<-rownames(df.kexp) %in% targets
  mt<-as.matrix(df.kexp[idxTpm,])

 for( i in 1:nrow(mt)){
  rg.idx<-which(rownames(mt)[i]==rownames(gwa.heat))
   replace<-as.character(gwa.heat$Gene.symbol[rg.idx])
    if(replace==""){
    replace<-rownames(gwa.heat)[rg.idx]
    } else if(!is.na(replace)==TRUE){
    rownames(mt)[i]<-replace
    } else {
   rownames(mt)[i]<-rownames(gwa.heat)[rg.idx]
   }
 }#for

  x.pv<-pvclust(log(1+mt),nboot=100)
  dend<-dendsort(hclust(dist(log(1+mt))),isReverse=TRUE)
  gwa.heatmap<-Heatmap(log(1+mt),
                   name="log(1+tpm)",
                   cluster_rows=dend,
                   cluster_columns=x.pv$hclust,
                   column_title=paste0("Limma Gene TMM Scaled P.Val ",setP," cutoff ",cutoffMax),
                   row_names_gp=gpar(fontsize=6),
                   column_names_gp=gpar(fontsize=8))

   draw(gwa.heatmap)
   readkey()


  #perform RWA with RUV limma
  rpt.tnx<-collapseTranscripts(kexp,read.cutoff=cutoffMax)
  rpt.idx<-!grepl("^ENS",rownames(rpt.tnx))
  rpts<-rpt.tnx[rpt.idx,]
  rpt.id<- !grepl("^ERCC",rownames(rpts))
  rpts<-rpts[rpt.id,]
  Rpts<-rpts[,colnames(rpts) %in% rownames(design)]

  RWA<-repeatWiseAnalysis(kexp,
                      design=design,
                       how="cpm",
                       p.cutoff=setP,
                       fold.cutoff=1,
                       read.cutoff=cutoffMax,
                       species=species,
                       adjustBy=adjustBy)
  if(saveReport==TRUE){
  write.csv(RWA$top,
              file=paste0(outputDir,"/limmaTMMrepeatWiseAnalysis_pval_",setP,"_cutoff_",cutoffMax,".csv"),
              row.names=FALSE)
   } #print to file

   rpt.tpm<-collapseTpm(kexp,"tx_id")
  df.rpt<-as.data.frame(rpt.tpm,stringsAsFactors=FALSE)
   
  if(nrow(RWA$top)<numberSelected){
    rpt.targets<-rownames(RWA$top)[1:nrow(RWA$top)]
   } else {
    rpt.targets<-rownames(RWA$top)[1:numberSelected]
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
                  column_title=paste0("Top Repeats TMM P.Val ",setP," ",adjustBy," cut ",cutoffMax),
                 row_names_gp=gpar(fontsize=6),
                 column_names_gp=gpar(fontsize=8)) 
 rh.rpt2<-Heatmap(features(kexp)[rownames(rpt.mt)]$tx_biotype,
          name="tx_biotype",
          width=unit(5,"mm"))
 rh.rpt3<- Heatmap(features(kexp)[rownames(rpt.mt)]$gene_biotype,
          name="rpt_family",
          width=unit(5,"mm"))
  } else {
  rh.rpt<-Heatmap(log(1+rpt.mt),
          name="log(1+tpm)",
          column_title=paste0("Top Repeats TMM P.Val ",setP," cut at ",cutoffMax," ",adjustBy),

          row_names_gp=gpar(fontsize=6))
  rh.rpt2<-Heatmap(features(kexp)[rownames(rpt.mt)]$tx_biotype,
          name="tx_biotype",
          width=unit(5,"mm"))
  rh.rpt3<-Heatmap(features(kexp)[rownames(rpt.mt)]$gene_biotype,
          name="rpt_family",
          width=unit(5,"mm"))
  }

  draw(rh.rpt+rh.rpt2+rh.rpt3)
  readkey()


  #edgeR too at the moment
   geneCounts<-collapseBundles(kexp,"gene_id",read.cutoff=cutoffMax)
  #phenoDat object assumes the coeff of A.vs.B is the second column, may encounter some order issues with rownames and the group A or B
  phenoDat<-data.frame(term=rownames(design),type=c(rep("A",nrow(design[design[,2]>0,])),rep("B",nrow(design[design[,2]<1,] ))))
  d2 <- DGEList(counts=geneCounts, group=phenoDat$type)
  d2 <- calcNormFactors(d2)
  d2 <- estimateGLMCommonDisp(d2,design)
  d2 <- estimateGLMTrendedDisp(d2, design)
  d2 <- estimateGLMTagwiseDisp(d2, design)
  d2 <- estimateGLMRobustDisp(d2, design) # superb, but slow!


  #for repeat dispersions
  txCounts<-collapseTranscripts(kexp,read.cutoff=cutoffMax)
  idx<-!grepl("^ENS",rownames(txCounts))
  rpts<-txCounts[idx,]
  id<- !grepl("^ERCC",rownames(rpts))
  rpts<-rpts[id,]
  d <- DGEList(counts=rpts, group=phenoDat$type)
  d <- calcNormFactors(d)
  d<- estimateGLMCommonDisp(d,design)
  d <- estimateGLMTrendedDisp(d, design)
  d <- estimateGLMTagwiseDisp(d, design)
  d <- estimateGLMRobustDisp(d, design) # superb, but slow!

  gn.fit<-glmFit(d2,design)
  gn.lrt<-glmLRT(gn.fit,coef=2) #assumes the coef is 2nd column, no multigroups supported yet
  gn.Tags<-topTags(gn.lrt,p=setP,n=nrow(kexp))[[1]]
  gn.feats<-mcols(features(kexp)[features(kexp)$gene_id %in% rownames(gn.Tags)])
  gn.key<-gn.feats[,c(4,5)]
  gn.key<-gn.key[!duplicated(gn.key$gene_id),]

  rpt.fit<-glmFit(d,design)
  rpt.lrt<-glmLRT(rpt.fit,coef=2)
  rpt.Tags<-topTags(rpt.lrt,p=setP,n=nrow(kexp))[[1]]
  if(saveReport==TRUE){
  write.table(gn.Tags,file=paste0(outputDir,"/TMM.gwa.edgeR.pval_",setP,"_cut_",cutoffMax,"_",adjustBy,".txt"),quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")
   write.table(rpt.Tags,file=paste0(outputDir,"/TMM.rwa.edgeR.pval_",setP,"_cut_",cutoffMax,"_",adjustBy,".txt"),quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")
  }
 
   edgeR.targets<-rownames(gn.Tags)[1:numberSelected]
  eg.tpm<-collapseTpm(kexp,"gene_id")
  eg.df.kexp<-as.data.frame(eg.tpm,stringsAsFactors=FALSE)
  eg.idxTpm<-rownames(eg.df.kexp) %in% edgeR.targets
  eg.mt<-as.matrix(eg.df.kexp[eg.idxTpm,])

  for(i in 1:nrow(eg.mt)){
   rg.idx<-which(rownames(eg.mt)[i]==gn.key$gene_id)
   replace<-unique(as.character(gn.key$gene_name[rg.idx]))
    if(!is.na(replace)==TRUE){
      if(any(rownames(eg.mt)==replace)){
      replace<-paste0(replace,".",i)
      rownames(eg.mt)[i]<-replace
     } else {
     rownames(eg.mt)[i]<-replace
     }
    }#if gene symbol exists
  }#for
  eg.x.pv<-pvclust(log(1+eg.mt),nboot=100)
  eg.dend<-dendsort(hclust(dist(log(1+eg.mt))),isReverse=TRUE)
  eg.gwa.heatmap<-Heatmap(log(1+eg.mt),
                   name="log(1+tpm)",
                   cluster_rows=eg.dend,
                   cluster_columns=eg.x.pv$hclust,
                   column_title=paste0("EdgeR Gene TMM P.Val ",setP," cutoff ",cutoffMax," top ",numberSelected), 
                   row_names_gp=gpar(fontsize=6),
                   column_names_gp=gpar(fontsize=8))
                   draw(eg.gwa.heatmap)
                   readkey()




 ##edgeR rwa analysis csv and heatmap
  rpt.targets<-rownames(rpt.Tags)
  rpt.tpm<-collapseTpm(kexp,"tx_id")
    
  rpt.idx<-!grepl("^ENS",rownames(rpt.tpm))
  eg.rpts.tpm<-rpt.tpm[rpt.idx,]
  tpm.id<- !grepl("^ERCC",rownames(eg.rpts.tpm))
  all.rpts.tpm<-eg.rpts.tpm[tpm.id,]

  eg.df.rpt<-as.data.frame(rpt.tpm,stringsAsFactors=FALSE)
  eg.df.rpt<-eg.df.rpt[rownames(eg.df.rpt) %in%rpt.targets,]
  eg.rpt.mt<-as.matrix(eg.df.rpt)

  eg.rwa.heatmap<-Heatmap(log(1+eg.rpt.mt),name="log(1+tpm)",
                column_title=paste0("EdgeR Top TMM Repeats TPM  P.Val ",setP," cut ",cutoffMax),
                row_names_gp=gpar(fontsize=8),
                column_names_gp=gpar(fontsize=8))
  eg.rwa.heatmap2<-Heatmap(features(kexp)[rownames(eg.rpt.mt)]$tx_biotype,name="tx_biotype",width=unit(5,"mm"))
 eg.rwa.heatmap3<- Heatmap(features(kexp)[rownames(eg.rpt.mt)]$gene_biotype,name="rpt_family",width=unit(5,"mm"))

  draw(eg.rwa.heatmap+eg.rwa.heatmap2+eg.rwa.heatmap3)
  readkey()
  #a single PDF with all images
  if(saveReport==TRUE){
  pdf(paste0(outputDir,"/GeneBundles_readCutoff_",cutoffMax,".pdf"))
  plotRLE(g.bundles,outline=FALSE)
  title("Gene Bundles Raw CPM")
  plotPCA(g.bundles,cex=1.2)
  title("PCA Gene Bundles Raw CPM")
  draw(gwa.heatmap)
  draw(rh.rpt+rh.rpt2+rh.rpt3)
  draw(eg.gwa.heatmap)
  draw(eg.rwa.heatmap+eg.rwa.heatmap2+eg.rwa.heatmap3)
  dev.off()
  }


  } else { #ruv selected
 ruv.gwa<-geneWiseAnalysis(kexp,design=design,how="cpm",p.cutoff=2,read.cutoff=cutoffMax,species=species,fitOnly=TRUE,adjustBy=adjustBy)
  tSel<-ruv.gwa$top
  sili<-tSel[order(tSel$adj.P.Val,decreasing=TRUE),]
  in.sili<-rownames(sili)[1:numberSelected] #take the top numberSelected highest pval
  #run RUV here and recast design matrix
print(plotRLE(g.bundles,outline=FALSE))
 title("Normalized Gene Bundles CPM")
  readkey()
 print(plotPCA(g.bundles,cex=1.2))
 title("Normalized PCA Gene Bundles CPM")
 readkey()
  gwaRuv<-RUVg(round(g.bundles),in.sili,k=1)
  ruvDesign<-cbind(design,gwaRuv$W)
  rownames(ruvDesign)<-colnames(kexp) 
  #plotting pca and RLE
  normalizedCounts<-gwaRuv$normalizedCounts
  #Run GWA with RUV limma
  gwa.RUV<-geneWiseAnalysis(kexp,
                            design=ruvDesign,
                            how="cpm",
                            p.cutoff=setP,
                            read.cutoff=cutoffMax,
                            species=species,
                            fitOnly=FALSE,
                            adjustBy=adjustBy)
  #gwa RUV heatmap
  if(saveReport==TRUE){
  write.table(in.sili,file=paste0(outputDir,"/neg.insilico_",setP,"_cut_",cutoffMax,".txt"),quote=FALSE,sep="\t")
  write.csv(gwa.RUV$limmaWithMeta[order(gwa.RUV$limmaWithMeta$adj.P.Val),],file=paste0(outputDir,"/gwaRUV_pval_",setP,"_readCutoff_",cutoffMax,".csv"),
  row.names=FALSE)
  write.table(ruvDesign,
              file=paste0(outputDir,"/ruvDesign_pval_",setP,"_cut_",cutoffMax,".txt"),
              quote=FALSE,sep="\t")
  } #print to file
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
    if(replace==""){
    replace<-rownames(gwa.heat)[rg.idx]
    } else if(!is.na(replace)==TRUE){
    rownames(ruv.mt)[i]<-replace
    } else {
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
                   column_title=paste0("Limma Gene Normalized P.Val ",setP," cutoff ",cutoffMax),
                   row_names_gp=gpar(fontsize=6),
                   column_names_gp=gpar(fontsize=8))
     
   draw(ruv.gwa.heatmap)
   readkey()


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
                       p.cutoff=setP,
                       fold.cutoff=1,
                       read.cutoff=cutoffMax,
                       species=species,
                       adjustBy=adjustBy)
  
  if(saveReport==TRUE){
  write.csv(ruvRWA$top,
              file=paste0(outputDir,"/NormalizedrepeatWiseAnalysis_pval_",setP,"_cutoff_",cutoffMax,".csv"),
              row.names=FALSE)
  }
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
                  column_title=paste0("RUV Repeats P.Val ",setP," ",adjustBy," cut ",cutoffMax),
                 row_names_gp=gpar(fontsize=6),
                 column_names_gp=gpar(fontsize=8))
  rh.rpt2<-Heatmap(features(kexp)[rownames(rpt.mt)]$tx_biotype,
          name="tx_biotype",
          width=unit(5,"mm"))
  rh.rpt3<-Heatmap(features(kexp)[rownames(rpt.mt)]$gene_biotype,
          name="rpt_family",
          width=unit(5,"mm"))
  } else {
  rh.rpt<-Heatmap(log(1+rpt.mt),
          name="log(1+tpm)",
          column_title=paste0("RUV Repeats P.Val ",setP," cut at ",cutoffMax," ",adjustBy),
          row_names_gp=gpar(fontsize=6))
  rh.rpt2<-Heatmap(features(kexp)[rownames(rpt.mt)]$tx_biotype,
          name="tx_biotype",
          width=unit(5,"mm"))
  rh.rpt3<-Heatmap(features(kexp)[rownames(rpt.mt)]$gene_biotype,
          name="rpt_family",
          width=unit(5,"mm"))
  }

  draw(rh.rpt+rh.rpt2+rh.rpt3)
  readkey()
 
  #run edgeR GWA and RWA with RUV at gene level
  geneCounts<-collapseBundles(kexp,"gene_id",read.cutoff=cutoffMax)
  #phenoDat object assumes the coeff of A.vs.B is the second column, may encounter some order issues with rownames and the group A or B
  phenoDat<-data.frame(term=rownames(design),type=c(rep("A",nrow(design[design[,2]>0,])),rep("B",nrow(design[design[,2]<1,] ))))
  d2 <- DGEList(counts=geneCounts, group=phenoDat$type)
  d2 <- calcNormFactors(d2)
  d2 <- estimateGLMCommonDisp(d2,ruvDesign)
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
  d<- estimateGLMCommonDisp(d,ruvDesign)
  d <- estimateGLMTrendedDisp(d, ruvDesign)
  d <- estimateGLMTagwiseDisp(d, ruvDesign)
  d <- estimateGLMRobustDisp(d, ruvDesign) # superb, but slow!
 
  gn.fit<-glmFit(d2,ruvDesign)
  gn.lrt<-glmLRT(gn.fit,coef=2) #assumes the coef is 2nd column, no multigroups supported yet
  gn.Tags<-topTags(gn.lrt,p=setP,n=nrow(kexp),adjust.method=adjustBy)[[1]]
  gn.feats<-mcols(features(kexp)[features(kexp)$gene_id %in% rownames(gn.Tags)])
  gn.key<-gn.feats[,c(4,5)]
  gn.key<-gn.key[!duplicated(gn.key$gene_id),]

  rpt.fit<-glmFit(d,ruvDesign)
  rpt.lrt<-glmLRT(rpt.fit,coef=2)
  rpt.Tags<-topTags(rpt.lrt,p=setP,n=nrow(kexp))[[1]]
  
  if(saveReport==TRUE){
  write.table(gn.Tags,file=paste0(outputDir,"/Normalized.gwa.edgeR.pval_",setP,"_cut_",cutoffMax,"_",adjustBy,".txt"),quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")
   write.table(rpt.Tags,file=paste0(outputDir,"/Normalized.rwa.edgeR.pval_",setP,"_cut_",cutoffMax,"_",adjustBy,".txt"),quote=FALSE,row.names=TRUE,col.names=TRUE,sep="\t")
  }
 
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
                   column_title=paste0("EdgeR Gene Normalized P.Val ",setP," cutoff ",cutoffMax," top ",numberSelected),
                   row_names_gp=gpar(fontsize=6),
                   column_names_gp=gpar(fontsize=8))

  draw(eg.ruv.gwa.heatmap)
  readkey()

  ##edgeR rwa analysis csv and heatmap
  rpt.targets<-rownames(rpt.Tags)
  rpt.tpm<-collapseTpm(kexp,"tx_id")
    
  rpt.idx<-!grepl("^ENS",rownames(rpt.tpm))
  eg.rpts.tpm<-rpt.tpm[rpt.idx,]
  tpm.id<- !grepl("^ERCC",rownames(eg.rpts.tpm))
  all.rpts.tpm<-eg.rpts.tpm[tpm.id,]

  eg.df.rpt<-as.data.frame(rpt.tpm,stringsAsFactors=FALSE)
  eg.df.rpt<-eg.df.rpt[rownames(eg.df.rpt) %in%rpt.targets,] 
  eg.rpt.mt<-as.matrix(eg.df.rpt)

  eg.rpt<-Heatmap(log(1+eg.rpt.mt),name="log(1+tpm)",
                column_title=paste0("EdgeR Top Repeats TPM  P.Val ",setP," cut ",cutoffMax),
                row_names_gp=gpar(fontsize=8),
                column_names_gp=gpar(fontsize=8))
  eg.rpt2<-Heatmap(features(kexp)[rownames(eg.rpt.mt)]$tx_biotype,name="tx_biotype",width=unit(5,"mm")) 
  eg.rpt3<-Heatmap(features(kexp)[rownames(eg.rpt.mt)]$gene_biotype,name="rpt_family",width=unit(5,"mm"))

  draw(eg.rpt+eg.rpt2+eg.rpt3)
  readkey()
  ###all repeats TPM
  all.rpt.mt<-as.matrix(all.rpts.tpm)
  all.rpt.x.pv<-pvclust(log(1+all.rpt.mt),nboot=100)
  all.rpt.dend<-dendsort(hclust(dist(log(1+all.rpt.mt))),isReverse=TRUE)


  all.eg.rpt<-Heatmap(log(1+all.rpt.mt),name="log(1+tpm)",
                cluster_rows=dend,
                column_title=paste0("EdgeR All Repeats TPM  P.Val ",setP," cut ",cutoffMax),
                row_names_gp=gpar(fontsize=8),
                column_names_gp=gpar(fontsize=8))
  all.eg.rpt2<-Heatmap(features(kexp)[rownames(all.rpt.mt)]$tx_biotype,name="tx_biotype",width=unit(5,"mm"))
 all.eg.rpt3<-Heatmap(features(kexp)[rownames(all.rpt.mt)]$gene_biotype,name="rpt_family",width=unit(5,"mm"))
  draw(all.eg.rpt+all.eg.rpt2+all.eg.rpt3)
  readkey()

  if(saveReport==TRUE){
  #print only 1 pdf with everything
  pdf(paste0(outputDir,"/NormalizedBundles_readcutoff_",cutoffMax,".pdf"))
  plotRLE(normalizedCounts,outline=FALSE,ylim=c(-2,2))
  title(main="Normalized Gene Bundles CPM")
  plotPCA(normalizedCounts,cex=1.2)
  title("Normalized Gene Bundles PCA")
  volcanoplot(gwa.RUV$fit,coef=2,highlight=20,names=(gwa.RUV$limmaWithMeta$Gene.symbol))
   title("Top Normalized Limma Gene CPM")
  #plot heatmap
   draw(ruv.gwa.heatmap)
   draw(rh.rpt+rh.rpt2+rh.rpt3)
  #include edgeR RUV
   plotBCV(d2)
   title("All Gene RUV EdgeR Dispersions")
   plotBCV(d)
   title("All Repeat RUV EdgeR Dispersions")
   draw(eg.ruv.gwa.heatmap)
   draw(eg.rpt+eg.rpt2+eg.rpt3)
   draw(all.eg.rpt+all.eg.rpt2+all.eg.rpt3)
  # include some beeswarm
  dev.off()
  }
  } #normalized
 } #gene id 
 else {
 # run collapse bundles by tx_id for an incrementation to cutoffMax.
 t.bundles<-collapseBundles(kexp,"tx_id",read.cutoff=cutoffMax)
 #print out RLE of gene and repeat regions, print out a PCA plot, no RUV here. only CPM data
 if(dataType=="unnormalized"){
 if(saveReport==TRUE){
 pdf(paste0(outputDir,"/TranscriptBundles_readCutoff_",cutoffMax,".pdf"))
 plotRLE(t.bundles,outline=FALSE)
 title("Transcript Bundles Raw CPM")
 plotPCA(t.bundles,cex=1.2)
 title("PCA Transcript Bundles Raw CPM")
 dev.off()
 } else {
 print(plotRLE(t.bundles,outline=FALSE))
 title("Transcript Bundles Raw CPM")
 readkey()
 print(plotPCA(t.bundles,cex=1.2))
 title("PCA Transcript Bundles Raw CPM")
 readkey()
    } 
  } else{ #for transcript level
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
  if(saveReport==TRUE){
  pdf(paste0(outputDir,"/NormalizedTranscriptBundles_readcutoff_",cutoffMax,".pdf"))
  plotRLE(normalizedCounts,outline=FALSE,ylim=c(-2,2))
  title(main="Normalized Transcript Bundles CPM")
  plotPCA(normalizedCounts,cex=1.2)
  title("Normalized Transcript Bundles PCA")
  dev.off()
  } else {
  print(plotRLE(normalizedCounts,outline=FALSE,ylim=c(-2,2)))
  title(main="Normalized Transcript Bundles CPM")
  readkey()
  print(plotPCA(normalizedCounts,cex=1.2))
  title("Normalized Transcript Bundles PCA")
  readkey()
   }
  #run TWA with RUV 
   
 }#dataType normalized
} 

} ###{{{ main
