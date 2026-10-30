
# Rscript R/mn_SCT_MAST_GO_GSEA.R INPUT.rds OUTPUT_DIR [SAMPLE_COLUMN]
# The input must already have celltype_res05 (meta data column for cell compartments) and original SCT counts/data.
suppressPackageStartupMessages({
  library(Seurat); library(MAST); library(Matrix); library(lme4)
  library(clusterProfiler); library(org.Hs.eg.db); library(msigdbr); library(fgsea)
})

select_mn <- function(obj) {
  md <- obj[[]]; genes <- c('SLC5A7','ISL1','ISL2','MNX1','STMN2','VSX2','SIM1')
  if (!'celltype_res05' %in% names(md)) stop('Missing celltype_res05')
  if (!'SCT' %in% SeuratObject::Assays(obj)) stop('Missing SCT assay')
  missing <- setdiff(genes,rownames(obj[['SCT']]))
  if(length(missing)) stop('Missing marker genes: ',paste(missing,collapse=', '))
  x <- SeuratObject::GetAssayData(obj,assay='SCT',layer='counts')[genes,,drop=FALSE]
  keep <- !is.na(md$celltype_res05) & md$celltype_res05=='Neuron' & x['SLC5A7',]>0 &
    (x['ISL1',]>0 | x['ISL2',]>0 | x['MNX1',]>0 | x['STMN2',]>5) &
    x['VSX2',]<1 & x['SIM1',]<1
  colnames(obj)[keep]
}

rank_genes <- function(d) {
  d <- d[!is.na(d$gene) & nzchar(d$gene) & is.finite(d$z_score),,drop=FALSE]
  d <- d[order(-abs(d$z_score),d$gene),,drop=FALSE]
  d <- d[!duplicated(d$gene),,drop=FALSE]
  d <- d[order(-d$z_score,d$gene),,drop=FALSE]
  setNames(d$z_score,d$gene)
}

filter_pathways <- function(pathways,genes) {
  pathways <- lapply(pathways,unique)
  pathways <- pathways[lengths(pathways)>=10L & lengths(pathways)<=500L]
  n <- vapply(pathways,function(x) length(intersect(x,genes)),integer(1))
  pathways[n>=10L & n<=500L]
}

write_table <- function(d,path) {
  d <- as.data.frame(d)
  for(nm in names(d)) if(is.list(d[[nm]])) d[[nm]] <- vapply(d[[nm]],paste,collapse='/',FUN.VALUE=character(1))
  write.csv(d,path,row.names=FALSE,na='')
}

fit_mast <- function(x,meta,mutant) {
  wanted <- meta$condition %in% c('CTRL',mutant)
  meta <- meta[wanted,,drop=FALSE]; x <- x[,rownames(meta),drop=FALSE]
  meta$condition <- factor(meta$condition,levels=c('CTRL',mutant))
  meta$sample_id <- factor(meta$sample_id)
  if(anyNA(meta$condition) || anyNA(meta$sample_id)) stop('Missing condition/sample')
  samples <- tapply(as.character(meta$condition),meta$sample_id,function(v) length(unique(v)))
  if(any(samples!=1L)) stop('Each sample must belong to one condition')
  n_per_group <- tapply(meta$sample_id,meta$condition,function(v) length(unique(v)))
  if(any(n_per_group<2L)) stop('At least two samples per condition are required for a random intercept')
  detection_ctrl <- Matrix::rowMeans(x[,meta$condition=='CTRL',drop=FALSE]>0)
  detection_mut <- Matrix::rowMeans(x[,meta$condition==mutant,drop=FALSE]>0)
  keep <- detection_ctrl>=.10 | detection_mut>=.10
  if(sum(keep)<2L) stop('Fewer than two genes pass the 10% detection filter')
  x <- x[keep,,drop=FALSE]
  # CDR is standardized number of detected genes after the 10% filter, per comparison.
  cdr <- Matrix::colSums(x>0)
  if(sd(cdr)==0) stop('Cellular detection rate has zero variance')
  cd <- data.frame(wellKey=colnames(x),condition=meta$condition,sample_id=meta$sample_id,
                   cngeneson=as.numeric(scale(cdr)),row.names=colnames(x))
  fd <- data.frame(primerid=rownames(x),row.names=rownames(x))
  sca <- MAST::FromMatrix(exprsArray=as.matrix(x),cData=cd,fData=fd,check_sanity=FALSE)
  fit <- MAST::zlm(~condition+cngeneson+(1|sample_id),sca,method='glmer',ebayes=FALSE,
    strictConvergence=FALSE,parallel=FALSE,
    fitArgsD=list(control=lme4::glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=200000))))
  contrast <- paste0('condition',mutant)
  s <- as.data.frame(summary(fit,doLRT=contrast)$datatable)
  h <- s[s$contrast==contrast & s$component=='H',c('primerid','Pr(>Chisq)'),drop=FALSE]
  z <- s[s$contrast==contrast & s$component=='logFC',c('primerid','coef','z'),drop=FALSE]
  if(anyDuplicated(h$primerid) || anyDuplicated(z$primerid)) stop('Duplicate MAST gene result rows')
  res <- data.frame(gene=rownames(x),pct_CTRL=unname(detection_ctrl[keep]),pct_mutant=unname(detection_mut[keep]))
  res$p_val <- h[['Pr(>Chisq)']][match(res$gene,h$primerid)]
  res$logFC <- z$coef[match(res$gene,z$primerid)]
  res$z_score <- z$z[match(res$gene,z$primerid)]
  res$p_adj <- NA_real_; ok <- is.finite(res$p_val)
  res$p_adj[ok] <- p.adjust(res$p_val[ok],method='BH',n=nrow(res))
  res$direction <- ifelse(is.finite(res$logFC) & res$logFC>0,'up',ifelse(is.finite(res$logFC) & res$logFC<0,'down','undetermined'))
  list(table=res,fit=fit,cdata=cd)
}

run_go <- function(deg,out,mutant) {
 for(direction in c('up','down')) {
  genes <- unique(deg$gene[!is.na(deg$p_adj) & deg$p_adj<=.05 & deg$direction==direction])
  dest <- file.path(out,paste0(mutant,'_GO_BP_',direction,'.csv'))
  if(length(genes)<5L) {writeLines('GO skipped: fewer than five input genes.',sub('[.]csv$','.txt',dest));next}
  universe <- unique(deg$gene[is.finite(deg$p_val)])
  go <- clusterProfiler::enrichGO(gene=genes,universe=universe,OrgDb=org.Hs.eg.db,keyType='SYMBOL',ont='BP',pAdjustMethod='BH',pvalueCutoff=.05,qvalueCutoff=.10,readable=TRUE)
  write_table(if(is.null(go)) data.frame() else as.data.frame(go),dest)
 }
}

run_gsea <- function(deg,pathways,out,mutant,collection) {
 ranks <- rank_genes(deg)
 if(length(ranks)<10L) stop('Too few finite z-scores for GSEA')
 write_table(data.frame(gene=names(ranks),z_score=unname(ranks)),file.path(out,paste0(mutant,'_GSEA_ranks.csv')))
 sets <- filter_pathways(pathways,names(ranks))
 if(!length(sets)) {writeLines('No eligible pathways after size filtering.',file.path(out,paste0(mutant,'_GSEA_',collection,'_skipped.txt')));return(invisible(NULL))}
 set.seed(42)
 g <- as.data.frame(fgsea::fgseaMultilevel(pathways=sets,stats=ranks,minSize=10,maxSize=500,eps=0,BPPARAM=BiocParallel::SerialParam()))
 g <- g[order(g$padj,g$pval,na.last=TRUE),,drop=FALSE]
 write_table(g,file.path(out,paste0(mutant,'_GSEA_',collection,'_all.csv')))
 write_table(g[is.finite(g$padj) & g$padj<=.10,,drop=FALSE],file.path(out,paste0(mutant,'_GSEA_',collection,'_FDR10.csv')))
}

main <- function(input,out,sample_column='sample') {
  if(as.character(packageVersion('fgsea'))!='1.36.2' || as.character(packageVersion('msigdbr'))!='26.1.1')
    warning('Installed fgsea/msigdbr versions differ from manuscript; see sessionInfo.txt')
  obj <- readRDS(input); md <- obj[[]]
  if(!all(c('condition',sample_column) %in% names(md))) stop('Missing condition or sample metadata')
  if(anyNA(md[,c('condition',sample_column),drop=FALSE])) stop('Missing condition or sample metadata')
  mn <- select_mn(obj) # fixed before any SCT preparation
  md <- md[mn,,drop=FALSE]; md$sample_id <- as.character(md[[sample_column]])
  if(!all(c('CTRL','FUS','TDP') %in% md$condition)) stop('MN selection lacks one or more conditions')
  x <- SeuratObject::GetAssayData(obj,assay='SCT',layer='data')[,mn,drop=FALSE]
  if(anyDuplicated(rownames(x)) || anyNA(rownames(x))) stop('SCT gene names must be unique and nonmissing')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  write_table(as.data.frame(table(condition=md$condition,sample=md$sample_id)),file.path(out,'MN_cell_counts.csv'))
  writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
  bp <- msigdbr::msigdbr(species='Homo sapiens',collection='C5',subcollection='GO:BP')
  hall <- msigdbr::msigdbr(species='Homo sapiens',collection='H')
  write_table(unique(bp[,intersect(c('db_version','db_species','gs_collection','gs_subcollection'),names(bp)),drop=FALSE]),file.path(out,'MSigDB_GO_metadata.csv'))
  write_table(unique(hall[,intersect(c('db_version','db_species','gs_collection','gs_subcollection'),names(hall)),drop=FALSE]),file.path(out,'MSigDB_Hallmark_metadata.csv'))
  saveRDS(list(GO_BP=bp,Hallmark=hall),file.path(out,'MSigDB_gene_sets.rds'))
  sets <- list(GO_BP=split(bp$gene_symbol,bp$gs_name),Hallmark=split(hall$gene_symbol,hall$gs_name))
  for(mutant in c('FUS','TDP')) {
    message('Fitting ',mutant,' vs CTRL')
    fit <- fit_mast(x,md,mutant)
    write_table(fit$table,file.path(out,paste0(mutant,'_MN_MAST_all_tested.csv')))
    write_table(fit$table[is.finite(fit$table$p_adj) & fit$table$p_adj<=.05,,drop=FALSE],file.path(out,paste0(mutant,'_MN_MAST_BH05.csv')))
    saveRDS(fit$fit,file.path(out,paste0(mutant,'_MAST_fit.rds')))
    write_table(fit$cdata,file.path(out,paste0(mutant,'_cell_detection_rate.csv')))
    run_go(fit$table,out,mutant)
    for(collection in names(sets)) run_gsea(fit$table,sets[[collection]],out,mutant,collection)
  }
  writeLines('Completed both comparisons',file.path(out,'COMPLETE.txt'))
}
if(sys.nframe()==0L) {
 args <- commandArgs(trailingOnly=TRUE)
 if(length(args)<2L || length(args)>3L) stop('Usage: Rscript R/mn_SCT_MAST_GO_GSEA.R INPUT.rds OUTPUT_DIR [SAMPLE_COLUMN]')
 main(args[1],args[2],if(length(args)==3L) args[3] else 'sample')
}
