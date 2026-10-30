
# Run: Rscript R/compartment_MAST_GO_adapted.R INPUT.rds OUTPUT_DIR

suppressPackageStartupMessages({library(Seurat); library(clusterProfiler); library(org.Hs.eg.db)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop('Supply input RDS and output directory')
obj <- readRDS(args[[1]])
out <- args[[2]]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
md <- obj[[]]
required <- c('condition','sample','celltype_res05')
if(!all(required %in% colnames(md))) stop('Missing metadata: ',paste(setdiff(required,colnames(md)),collapse=', '))
if(!all(c('CTRL','FUS','TDP') %in% unique(md$condition))) stop('Expected CTRL, FUS, TDP')
if(anyNA(md[,required])) stop('Missing values in required metadata')
if(!'SCT' %in% Assays(obj)) stop('Missing SCT assay')
if(!inherits(obj[['SCT']], 'SCTAssay')) stop('SCT must be a fitted SCTAssay with its models')
# Detect this before any costly model fitting; do not silently change the model.
for (comp in unique(md$celltype_res05)) for (mut in c('FUS','TDP')) {
 d <- droplevels(data.frame(condition=factor(md$condition),sample=factor(md$sample)) [md$celltype_res05==comp & md$condition %in% c('CTRL',mut),,drop=FALSE])
 if(nrow(d)==0 || nlevels(d$condition)<2 || nlevels(d$sample)<2) stop('Insufficient samples in ',comp,' / ',mut)
 design <- model.matrix(~sample+condition,data=d)
 if(qr(design)$rank < ncol(design)) stop('Condition is confounded with sample in ',comp,' / ',mut, '. Sample cannot be used as a fixed latent variable to estimate the condition effect; use a suitable sample-level or random-effect model.')
}
set.seed(42)
md <- md[md$condition %in% c('CTRL','FUS','TDP'),,drop=FALSE]
selected <- unlist(lapply(split(rownames(md),md$condition,drop=TRUE),function(cells) sample(cells,min(table(md$condition)))),use.names=FALSE)
obj <- subset(obj,cells=selected)
DefaultAssay(obj) <- 'SCT'
obj <- PrepSCTFindMarkers(obj, assay='SCT', verbose=TRUE)
Idents(obj) <- 'condition'
for (comp in unique(obj$celltype_res05)) for (mut in c('FUS','TDP')) {
 part <- subset(obj,subset=celltype_res05==comp)
 result <- FindMarkers(part,assay='SCT',slot='data',ident.1=mut,ident.2='CTRL',test.use='MAST',logfc.threshold=0,min.pct=0,latent.vars='sample',recorrect_umi=FALSE)
 result$gene <- rownames(result)
 result$BH_adjusted_P <- p.adjust(result$p_val,method='BH')
 stem <- paste(gsub('[^A-Za-z0-9_-]','_',comp),mut,sep='_')
 write.csv(result,file.path(out,paste0(stem,'_DEG.csv')),row.names=FALSE)
 fc <- if('avg_log2FC' %in% names(result)) 'avg_log2FC' else stop('Missing avg_log2FC')
 for (direction in c('up','down')) {
  genes <- unique(result$gene[result$BH_adjusted_P<=0.05 & if(direction=='up') result[[fc]]>0 else result[[fc]]<0])
  if(length(genes)<5L) {message(stem,' ',direction,': fewer than 5 genes; GO skipped');next}
  ego <- enrichGO(gene=genes,universe=unique(result$gene),OrgDb=org.Hs.eg.db,keyType='SYMBOL',ont='BP',pAdjustMethod='BH',pvalueCutoff=0.05,qvalueCutoff=0.10,readable=TRUE)
  write.csv(as.data.frame(ego),file.path(out,paste0(stem,'_GO_BP_',direction,'.csv')),row.names=FALSE)
 }
}
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
