library(plyr)
library(dplyr)

## quick function to simplify gene annotation provided by 450k gene annotation
## which often repeats the same gene symbol between semi-colons
simple.gene <- function(gene.vec)
  laply(gene.vec, function(x) paste(unique(strsplit(x, ';')[[1]]), collapse = ';'))


## from betas or mvals, using rownames basically
## unless specifically changed, removes all fields in rm.fld
getCleanAnno <- function(mat, rm.fld = c('AddressA', 'AddressB', 'ProbeSeqA', 'ProbeSeqB', 'NextBase', 'Color', 'Forward_Sequence', 'SourceSeq', 'Random_Loci', 'Methyl27_Loci', 'UCSC_RefGene_Accession', 'Phantom', 'HMM_Island', 'Islands_Name', 'Regulatory_Feature_Name')){
  library(minfi)

  anno <- RatioSet(mat, annotation = c(array = 'IlluminaHumanMethylation450k',
                                       annotation = 'ilmn12.hg19'))
  anno <- getAnnotation(anno)
  stopifnot(all(rm.fld %in% names(anno)))
  anno <- anno[, !names(anno) %in% rm.fld]
  anno$Gene <- simple.gene(anno$UCSC_RefGene_Name)
  anno
}

pca.theme <- function(i = 1){
  require(ggplot2)
  theme(text = element_text(colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = NA),
        panel.border = element_rect(fill = NA, colour = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(size = rel(i), colour = 'black'),
        axis.title = element_text(size = rel(i)),
        plot.title = element_text(size = rel(i)))
}

my.grepl <- function(vec, ..., fun = c('any', 'all')){
  lis <- lapply(vec, function(x) grepl(x, ...))
  foo <- do.call(rbind, lis)
  fun <- get(match.arg(fun))
  apply(foo, 2, fun)
}

## limma function:
limma.tt <- function(mval.mat, betas.mat, pheno, var.interest, covars = NULL, fDat){
  library(limma)
  stopifnot(all(c(covars, var.interest) %in% names(pheno)),
            all(colnames(betas.mat) == colnames(mval.mat)),
            all(c('GeneSymbol', 'cpgSite') %in% names(fDat)))
  form <- paste('~', var.interest)
  if(!is.null(covars))
    form <- paste(form, '+', paste(covars, collapse = ' + '))
  des <- model.matrix(as.formula(form), data = pheno)
  cat('Limma on m-values... \n')
  fit <- lmFit(mval.mat, des)
  eb.fit <- eBayes(fit)#, trend = TRUE, robust = TRUE
  tt <- topTable(eb.fit, coef = 2, number = nrow(mval.mat), adjust.method = 'BH')
  tt$Gene <- fDat$GeneSymbol[match(rownames(tt), fDat$cpgSite)]
  cat('Limma on beta values... \n')
  fit.beta <- lmFit(betas.mat, des)
  eb.fit <- eBayes(fit.beta)
  tt.beta <- topTable(eb.fit, coef = 2, number = nrow(betas.mat), adjust.method = 'BH')
  tt$betaFC <- tt.beta[rownames(tt), 'logFC']
  tt
}

## returns a list of toptables, one for each column of contr.mat
limma.tt.contrasts <- function(mval.mat, betas.mat, des.mat, contr.mat, fDat){
  options(warn = -1)
  library(limma)
  stopifnot(all(colnames(betas.mat) == colnames(mval.mat)),
            all(c('GeneSymbol', 'cpgSite') %in% names(fDat)))

  cat('\n', 'Limma on m-values... \n ')
  fit <- lmFit(mval.mat, des.mat)
  fit <- contrasts.fit(fit, contr.mat)
  eb.fit <- eBayes(fit)
  tt.list <- llply(1:ncol(contr.mat), function(x)
    topTable(eb.fit, coef = x, number = nrow(mval.mat)), .progress = 'text')

  cat('\n', 'Adding gene names... \n')
  tt.list <- llply(tt.list, function(tt) {
    tt$Gene <- fDat$GeneSymbol[match(rownames(tt), fDat$cpgSite)]
    tt}, .progress = 'text')

  cat('\n Limma on beta values... \n')
  fit <- lmFit(betas.mat, des.mat)
  fit <- contrasts.fit(fit, contr.mat)
  eb.fit <- eBayes(fit)
  tt.list.beta <- llply(1:ncol(contr.mat), function(x)
    topTable(eb.fit, coef = x, number = nrow(betas.mat)), .progress = 'text')

  tt.list <- mapply(function(x, y){
    x$betaFC <- y[rownames(x), 'logFC']
    x}, tt.list, tt.list.beta, SIMPLIFY = F)
  names(tt.list) <- colnames(contr.mat)
  options(warn = 0)
  tt.list
}

## function which takes the differences of mvals and betas
## for each PAIR of samples for each SUBJECT as denoted by subj.fld in pheno,
## taken corresponding to the var.interest which must have two levels for each subject
## limma applied with now difference of values wrt var.interest, and controlled for
## covars, which mean control for difference now instead of methylation value
## the
limma.tt.paired <- function(mval.mat, betas.mat, pheno, subj.fld, var.interest,
                            covars.diff = NULL, covars.fixed = NULL, fDat){
  library(limma)
  library(plyr)
  library(dplyr)
  stopifnot(all(c(covars.diff, covars.fixed, var.interest, subj.fld) %in% names(pheno)),
            all(colnames(betas.mat) == colnames(mval.mat)),
            all(c('GeneSymbol', 'cpgSite') %in% names(fDat)),
            ncol(mval.mat) == nrow(pheno))
  id.fld <- names(pheno)[which(laply(pheno, function(x) all(x == colnames(betas.mat))))[1]]
  phen <- pheno[c(id.fld, subj.fld, var.interest, covars.diff, covars.fixed)]
  names(phen)[1:3] <- c('id', 'subj', 'var.int')
  if(!is.null(covars.diff)){
    dups <- phen %>% group_by(subj) %>% dplyr::select(-one_of(covars.diff)) %>%
      summarise(n.samp = n_distinct(var.int))
  } else dups <- phen %>% group_by(subj) %>% summarise(n.samp = n_distinct(var.int))
  stopifnot(all(dups$n.samp == 2))
  levs <- sort(unique(phen$var.int))
  mval.diff <- t(daply(phen, .(subj), function(x){
    samps <- x$id[match(levs, x$var.int)]
    mval.mat[,samps[2]] - mval.mat[,samps[1]]}))
  rownames(mval.diff) = rownames(mval.mat)
  beta.diff <- t(daply(phen, .(subj), function(x){
    samps <- x$id[match(levs, x$var.int)]
    betas.mat[,samps[2]] - betas.mat[,samps[1]]}))
  rownames(beta.diff) = rownames(betas.mat)

  ## create difference for each covars.diff:
  for(v in covars.diff){
    tmp <- phen[c('subj', 'var.int', v)]; names(tmp)[3] <- 'v'
    tmp <- tmp %>% group_by(subj) %>% arrange(var.int) %>% summarise(v = diff(v))
    phen <- phen %>% dplyr::select(-one_of(v)) %>% merge(tmp, all.x = T)
    names(phen)[ncol(phen)] <- paste0(v, '_Diff')
  }
  covars <- setdiff(names(phen), c('subj', 'id', 'var.int'))
  if(length(covars) == 0) covars <- NULL

  phen <- phen %>% tbl_df %>% dplyr::select(-id, -var.int) %>% group_by(subj) %>%
    distinct() %>% arrange(subj)
  form <- '~ 1'
  if(!is.null(covars))
    form <- paste(form, '+', paste(covars, collapse = ' + '))
  des <- model.matrix(as.formula(form), data = phen)
  stopifnot(all(phen$subj == colnames(beta.diff)),
            all(phen$subj == colnames(mval.diff)))
  cat('\n', 'Limma on m-values...')
  fit <- lmFit(mval.diff, des)
  eb.fit <- eBayes(fit)
  tt <- topTable(eb.fit, coef = 1, number = nrow(mval.diff), adjust.method = 'BH')
  cat('\n Re-annotating gene names...')
  tt$Gene <- fDat$GeneSymbol[match(rownames(tt), fDat$cpgSite)]
  cat('\n Limma on beta values... \n')
  fit.beta <- lmFit(beta.diff, des)
  eb.fit <- eBayes(fit.beta)
  tt.beta <- topTable(eb.fit, coef = 1, number = nrow(beta.diff), adjust.method = 'BH')
  tt$betaFC <- tt.beta[rownames(tt), 'logFC']
  tt
}

## linear modeling functions for gene expression## linear modeling functions
limma.tt.exprs <- function(mat, pheno, var.interest, form, fDat = NULL){
  library(limma)
  #stopifnot(c('Probe','GeneSymbol') %in% colnames(fDat))
  form <- as.formula(form)

  des <- model.matrix(form, data=pheno)
  fit <- lmFit(mat, des)
  eb.fit <- eBayes(fit)

  tt <- topTable(eb.fit, coef = var.interest, number = nrow(mat), adjust.method = 'BH')
  #cat(paste0('\nreturning coef of ', var.interest,' \n'))

  if(!is.null(fDat)){
    tt$Gene <- fDat$GeneSymbol[match(rownames(tt), fDat$Probe)]
    tt$Probe <- rownames(tt)
    tt %>% dplyr::select(Gene, Probe, logFC, AveExpr, t, P.Value, adj.P.Val)
  } else{
    tt
  }
}

lm.tt.exprs <- function(mat, pheno, var.interest, form, fDat=NULL){

  form <- as.formula(form)
  res = data.frame(probe_id = rownames(mat), estimate=0, se=0, t=0, P=0)

  for(i in 1:nrow(mat)){
    pheno[,exprs := mat[i,]]
    mod = summary(lm(as.formula(form), data = pheno))
    res[i,2:5]=mod$coefficients[var.interest,]
    pheno[, exprs := NULL]
    if(i%%500==0) cat(paste0('LM of row ',i,' finished! \n\n'))
  }
  if(!is.null(fDat)){
    res %>% arrange(P) %>%
      mutate(FDR = p.adjust(P, method = 'BH'),
             gene_symbol = fDat$gene_symbol[match(probe_id, fDat$probe_id)]) %>%
      dplyr::select(gene_symbol, probe_id, estimate, se, t, P, FDR)
  }else{
    res %>% arrange(P) %>%
      mutate(FDR = p.adjust(P, method = 'BH')) %>%
      dplyr::select(probe_id, estimate, se, t, P, FDR)
  }
}


limma.tt.contrasts.exprs <- function(exprs.mat, des.mat, contr.mat, fDat){
  options(warn = -1)
  library(limma)
  stopifnot(all(c('GeneSymbol', 'Probe') %in% names(fDat)))

  cat('\n', 'Limma on expression... \n ')
  fit <- lmFit(exprs.mat, des.mat)
  fit <- contrasts.fit(fit, contr.mat)
  eb.fit <- eBayes(fit)
  tt.list <- llply(1:ncol(contr.mat), function(x)
    topTable(eb.fit, coef = x, number = nrow(exprs.mat)), .progress = 'text')

  cat('\n', 'Adding gene names... \n')
  tt.list <- llply(tt.list, function(tt) {
    tt$Gene <- fDat$GeneSymbol[match(rownames(tt), fDat$Probe)]
    tt}, .progress = 'text')

  names(tt.list) <- colnames(contr.mat)
  options(warn = 0)
  tt.list
}



## helper functions for DMRcate:

## this should plot the top 5 DMRs in the dmrcate.output object with meaningful
## file suffix and plot title
plotTop_n <- function(res.obj, filePrefix, betas, colors, n = 5){
  tab <- res.obj$results
  tab <- subset(tab, no.probes > 1)
  res.obj$results <- subset(res.obj$results, no.probes > 1)
  for(i in 1:n){
    print(i)
    gene <- ifelse(tab$gene_assoc[i] == '', 'noGeneAssoc',
                   gsub(',', '_', tab$gene_assoc[i]))
    title <- sprintf('DMR consisting of %d probes with mean p-value %.4e',
                     tab$no.probes[i], tab$meanpval[i])
    filename <- sprintf('%s_%s_%s.png', filePrefix,
                        gsub(':.*', '', tab$hg19coord[i]), gene)
    png(filename, width = 9, height = 8, units = 'in', res = 100)
    DMR.plot(res.obj, i, betas = betas, phen.col = colors,
             toscale = T, plotmedians = T, pch = 16, main = title)
    dev.off()
  }
}

plotGenes <- function(res.obj, genes, filePrefix, betas, colors){
  tab <- res.obj$results
  for(g in genes){
    if(!g %in% tab$gene_assoc){
      cat(sprintf('\n %s not found in results \n', g))
      if(grepl(g, tab$gene_assoc)){
        cat(sprintf('perhaps you meant %s?'),
            paste(grep(g, tab$gene_assoc, value = T), collapse = ' or '))
      }
    } else{
      print(g)
      i <- which(tab$gene_assoc == g)
      title <- sprintf('DMR consisting of %d probes with mean p-value %.4e',
                       tab$no.probes[i], tab$meanpval[i])
      filename <- sprintf('%s_%s_%s.png', filePrefix,
                          gsub(':.*', '', tab$hg19coord[i]), g)
      png(filename, width = 9, height = 8, units = 'in', res = 100)
      DMR.plot(res.obj, i, betas = betas, phen.col = colors,
               toscale = T, plotmedians = T, pch = 16, main = title)
      dev.off()
    }
  }
}

## for comparisons using a continuous predictor, this function plots scatterplots
## across a DMR region for each probe

DMR.plot.cont <- function(dmrcoutput, dmr, filePrefix, betas,
                          phen, phenoDat, id.fld = 'STUDY_ID',
                          annotation = c(array = 'IlluminaHumanMethylation450k',
                                         annotation = 'ilmn12.hg19')){
  library(ggplot2)
  stopifnot(is(dmrcoutput, "dmrcate.output"))
  stopifnot(is.matrix(betas))
  stopifnot((length(dmr) == 1) && (dmr %in% 1:nrow(dmrcoutput$results)))
  coords <- dmrcoutput$results$hg19coord[dmr]
  chr <- sub(":.*", "", coords)
  bookends <- sub(".*:", "", coords)
  startcpg <- as.integer(sub("-.*", "", bookends))
  stopcpg <- as.integer(sub(".*-", "", bookends))
  RSobject <- RatioSet(betas, annotation = annotation)
  RSanno <- getAnnotation(RSobject)
  cpgs <- rownames(RSanno)[RSanno$chr %in% chr & RSanno$pos >=
                             startcpg & RSanno$pos <= stopcpg]
  cpgs <- cpgs[order(RSanno[cpgs, "pos"])]
  betas <- betas[as.character(cpgs),]
  m <- match(as.character(cpgs), rownames(RSanno))
  clusta <- data.frame(gene = simple.gene(RSanno$UCSC_RefGene_Name[m]),
                       group = RSanno$UCSC_RefGene_Group[m], pos = RSanno$pos[m])
  clusta <- cbind(clusta, cpgSite = rownames(betas), betas)
  plotDat <- melt(clusta, id.vars = c('gene', 'group', 'pos', 'cpgSite'))
  names(plotDat)[5:6] <- c('ID', 'beta')
  plotDat$pheno <- phenoDat[[phen]][match(plotDat$ID, phenoDat[[id.fld]])]
  plotDat$pos2 <- with(plotDat, paste(gene, group, pos, sep = '_'))
  g <- ggplot(plotDat, aes(x = pheno, y = beta)) + geom_point() +
    facet_wrap(~pos2, nrow = 1) + ylim(0,1) +
    geom_hline(yintercept = c(0,0.5, 1), color = 'lightgrey', alpha = 0.5) +
    stat_smooth(method = 'lm') + theme_bw() +
    xlab(phen)
  ggsave(paste0(filePrefix, '.png'), g, height = 5, width = 2*nrow(clusta),
         limitsize = F)

}

## plot top 5 for coninuous variable:
plotTop5.cont <- function(res.obj, filePrefix, betas, phen, phenoDat,
                          id.fld = 'STUDY_ID'){
  tab <- res.obj$results
  tab <- tab[order(tab$meanpval),]
  for(i in 1:5){
    print(i)
    gene <- ifelse(tab$gene_assoc[i] == '', 'noGeneAssoc',
                   gsub(',', '_', tab$gene_assoc[i]))
    title <- sprintf('DMR consisting of %d probes with mean p-value %.4e',
                     tab$no.probes[i], tab$meanpval[i])
    filepref <- sprintf('%s_%s_%s', filePrefix,
                        gsub(':.*', '', res.obj$results$hg19coord[i]), gene)
    DMR.plot.cont(res.obj, i, filepref, betas, phen, phenoDat, id.fld)
  }
}


## writes all the gene_assocs found in the dmrcate.output object OR a top-table
## to a simple text list
write.genes <- function(res.obj, file, thresh = 0.05){
  if(is(res.obj, 'data.frame')){
    stopifnot(all(c('Gene', 'adj.P.Val') %in% names(res.obj)))
    genes <- subset(res.obj, adj.P.Val < thresh)$Gene
    genes <- unique(unlist(strsplit(genes, ';')))
  } else if(is(res.obj, 'dmrcate.output')){
    cat(sprintf('\n %d DMR regions in the results file \n', nrow(res.obj$results)))
    genes <- res.obj$results$gene_assoc
    genes <- unique(unlist(strsplit(genes, ',')))
    cat(sprintf('%d unique DMR-associated genes identified \n', length(genes)))
  } else stop('res.obj is the wrong kind of format')
  write.table(unique(genes), file = file, col.names = F, row.names = F,
              quote = F, sep = '\n')
}

write.results <- function(res.obj, file){
  tab <- res.obj$results
  write.csv(tab, file, row.names = F)
}

getBadCG <- function(dist = 2, mafcut = 0.05, crosshyb = T){
  library(DMRcatedata)
  library(plyr)
  env <- new.env(parent = emptyenv())
  data(dmrcatedata, envir = env)
  cross.hyb <- as.character(env$crosshyb)
  snp.probes <- env$illuminaSNPs
  len0 <- vapply(snp.probes$Distance, length, integer(1), USE.NAMES = F)
  len1 <- vapply(snp.probes$MinorAlleleFrequency, length, integer(1), USE.NAMES = F)
  snp.probes <- snp.probes[len0 == len1,]
  dist0 <- as.integer(unlist(snp.probes$Distance, use.names = FALSE))
  distrange <- range(dist0)
  stopifnot(dist >= min(distrange) && dist <= max(distrange))
  test0 <- (dist0 >= -1) & (dist0 <= dist)
  test1 <- unlist(snp.probes$MinorAlleleFrequency, use.names = FALSE) > mafcut
  test <- (test0 & test1)
  ntrue <- cumsum(test)[cumsum(len0)]
  badidxs <- ntrue - c(0, head(ntrue, -1)) != 0
  out <- data.frame(CG = row.names(snp.probes)[badidxs],
                    Detail = laply(snp.probes$Distance[badidxs],
                                   function(x) paste(x, collapse = ',')))
  rbind(out, data.frame(CG = cross.hyb,
                        Detail = 'CrossHybridized'))
}

subset.dmrcate <- function(res.obj, rmBlankAssoc = T, min.no.probes = 5){
  foo <- res.obj
  if(rmBlankAssoc)
    foo$results <- subset(res.obj$results, gene_assoc != '')
  foo$results <- subset(foo$results, no.probes >= min.no.probes)
  foo
}

summarize.tt <- function(tt){
  stopifnot(all(c('adj.P.Val', 'betaFC') %in% names(tt)))
  sig <- tt$adj.P.Val < 0.05
  sig.eff <- sum(abs(tt$betaFC[sig]) > 0.2 )
  c(Signif = sum(sig), Sig.Eff = sig.eff)
}

summarize.dmrcate <- function(res.obj){
  stopifnot(is(res.obj, 'dmrcate.output'))
  tab <- subset(res.obj$results, no.probes >= 3)
  sig.eff <- sum(abs(tab$maxbetafc) > 0.2)
  ann <- sum(tab$gene_assoc != '')
  c(min3 = nrow(tab), min3_beta02 = sig.eff, min3_anno = ann)
}

## produces boxplots for CG sites stratified by plotVar (must be categorical)
## pheno rows must correspond to mval columns
sanity.check <- function(cgSites, plotVar, colorVar, shapeVar = NULL, beta, pheno){
  library(reshape2)
  cat('\n Make sure the rows in pheno correspond to the columns in beta \n ')
  stopifnot(all(c(plotVar, colorVar, shapeVar) %in% names(pheno)),
            all(cgSites %in% rownames(beta)),
            !is.numeric(pheno[[plotVar]]))
  if(!is.null(shapeVar)) stopifnot(!is.numeric(pheno[[shapeVar]]))
  plotDat <- cbind(as.data.frame(t(beta[cgSites,])), var = pheno[[plotVar]],
                   color = pheno[[colorVar]])
  if(!is.null(shapeVar)) plotDat$shape = pheno[[shapeVar]]
  idvars <- c('var', 'color')
  if(!is.null(shapeVar)) idvars <- c(idvars, 'shape')
  plotDat <- melt(plotDat, id.vars = idvars, variable.name ='CGsite', value.name = 'Beta')
  g <- ggplot(plotDat, aes(x = var, y = Beta)) +
    geom_boxplot(outlier.size = 0) +
    xlab(plotVar) + pca.theme() + facet_wrap(~CGsite, scales = 'free_y')
  if(is.null(shapeVar)){
    g <- g + geom_point(aes(color = color),
                        position = position_jitter(width = 0.2, height = 0))
  } else{
    g <- g + geom_point(aes(shape = shape, color = color),
                        position = position_jitter(width = 0.2, height = 0)) +
      scale_shape_discrete(name = shapeVar)
  }
  if(is.numeric(plotDat$color))
    g <- g + scale_color_gradient(low = 'blue', high = 'red', name = colorVar) else
      g <- g + scale_color_discrete(name = colorVar)
  g + ggtitle(sprintf('Beta-values for various CpG sites grouped by %s', plotVar))
}

### same as above, but uses paired samples for each group denoted by plotVar
### i.e. in pheno, each subject in subjVar must have exactly two levels of plotVar
sanity.check.paired <- function(cgSites, plotVar, colorVar, beta, pheno, subjVar){
  library(reshape2)
  cat('\n Make sure the rows in pheno correspond to the columns in beta \n ')
  stopifnot(all(c(plotVar, colorVar, subjVar) %in% names(pheno)),
            all(cgSites %in% rownames(beta)),
            !is.numeric(pheno[[plotVar]]),
            !is.numeric(pheno[[subjVar]]))
  plotDat <- cbind(as.data.frame(t(beta[cgSites,])), var = pheno[[plotVar]],
                   color = pheno[[colorVar]], subj = pheno[[subjVar]])
  stopifnot(all(table(plotDat$subj) == 2),
            plotDat %>% group_by(subj) %>% summarise(n.samp = n_distinct(var)) %>%
              with(all(n.samp == 2)))
  plotDat <- melt(plotDat, id.vars = c('var', 'color', 'subj'),
                  variable.name ='CGsite', value.name = 'Beta')
  g <- ggplot(plotDat, aes(x = var, y = Beta, color = color)) +
    geom_boxplot(outlier.size = 0) +
    xlab(plotVar) + pca.theme() + facet_wrap(~CGsite, scales = 'free_y') +
    geom_line(aes(group = subj))
  if(is.numeric(plotDat$color))
    g <- g + scale_color_gradient(low = 'blue', high = 'red', name = colorVar) else
      g <- g + scale_color_discrete(name = colorVar)
  g + ggtitle(sprintf('Beta-values for various CpG sites paired by %s', plotVar))
}

## generates "Table 1" for two-group (case/control) comparisons
table1.cc <- function(pheno, var.int, cat.fld, cont.fld){
  library(plyr)
  stopifnot(all(c(var.int, cat.fld, cont.fld) %in% names(pheno)),
            length(unique(pheno[[var.int]])) == 2)
  chisq.fn <- function(fld){
    tab <- table(pheno[[fld]], pheno[[var.int]])
    res <- chisq.test(tab, correct = F)
    cbind(tab, c(rep('', nrow(tab) - 1), format(res$p.value, digits = 3)))
  }
  ttest.fn <- function(fld){
    p <- pheno[c(var.int, fld)]
    names(p) <- c('var.int', 'fld')
    tab <- unlist(daply(p, .(var.int), summarise,
                        paste(format(mean(fld, na.rm = T), digits = 2),
                              format(sd(fld, na.rm = T)/sqrt(sum(!is.na(fld))),
                                     digits = 2), sep = '+/-')))
    res <- t.test(fld ~ var.int, data = p)
    matrix(c(tab, format(res$p.value, digits = 3)), nrow = 1)
  }
  cat.foo <- do.call('rbind', llply(cat.fld, chisq.fn))
  cont.foo <- do.call('rbind', llply(cont.fld, ttest.fn))
  rownames(cont.foo) <- cont.fld
  foo <- rbind(cat.foo, cont.foo)
  colnames(foo)[3] <- 'pvalue'
  foo
}

## these use sanity.check.paired for each group of CpG in the top list of a
## DMRcate output object, res.obj
plotTop_n_paired <- function(res.obj, filePrefix, plotVar, colorVar, beta, pheno, subjVar, fdat= fDat, n =5, height = 10, width = 15){
  tab <- res.obj$results
  stopifnot('cpgSite' %in% names(fdat))
  for(i in 1:n){
    print(i)
    gene <- ifelse(tab$gene_assoc[i] == '', 'noGeneAssoc',
                   gsub(',', '_', tab$gene_assoc[i]))
    title <- sprintf('DMR consisting of %d probes with mean p-value %.4e',
                     tab$no.probes[i], tab$meanpval[i])
    chr <- gsub(':.*', '', tab$hg19coord[i])
    pos1 <- as.numeric(gsub('.*:', '', gsub('-.*', '', tab$hg19coord[i])))
    pos2 <- as.numeric(gsub('.*-', '', tab$hg19coord[i]))
    filename <- sprintf('%s_%s_%s.png', filePrefix,
                        chr, gene)
    cpgSites <- fDat %>% filter(chr == chr, pos >= pos1, pos <= pos2) %>%
      arrange(pos) %>% with(cpgSite)
    g <- sanity.check.paired(cpgSites, 'VIDUS_Time', 'CD4T_Diff', beta, pheno, 'VIDUS_ID') +
      ggtitle(title)
    ggsave(file = filename, g, height = 10, width = 15)
  }
}

## volcano plot function
### Label points with label.fld (if not NULL)
### use label.thresh and label.thresh.fld to determine what gets labeled.
volcano_plot <- function(toptable, contrast.fld = 'betaFC', pval.fld = 'P.Value',
                         pthresh = 0.01, contrast.thresh = 0.25,
                         xlabel = 'beta difference',
                         ylabel = '-log10 p-value',
                         title = '',
                         line.label='FDR=0.1',
                         label.fld = 'Gene',
                         label.thresh = 0.05, label.thresh.fld = 'adj.P.Val',
                         colors = c(neut = '#8b8b8b', nonsig='grey', low = 'blue', high = 'red')){
  require(ggplot2);require(ggrepel);require(ggrastr)

  df <- data.frame(contrast = toptable[,contrast.fld],
                   pval = toptable[,pval.fld])
  df$colorClass <- 'nonsig'
  df$colorClass[df$pval <= pthresh] <- 'neut'
  df$colorClass[df$pval <= pthresh & df$contrast <= -contrast.thresh] <- 'Hypomethylated'
  df$colorClass[df$pval <= pthresh & df$contrast >= contrast.thresh] <- 'Hypermethylated'
  names(colors)[3:4] <- c('Hypomethylated', 'Hypermethylated')

  ### labeling:
  df$label <- ''
  df$hjust <- as.numeric(df$contrast < 0)
  if(!is.null(label.fld)){
    if(label.thresh.fld == 'betaFC'){
      ix <- which(abs(toptable[[label.thresh.fld]]) >= label.thresh)
      df$label[ix] <- toptable[[label.fld]][ix]
    } else{
      ix <- which(toptable[[label.thresh.fld]] <= label.thresh & abs(toptable[[contrast.fld]]) >= contrast.thresh)
      df$label[ix] <- toptable[[label.fld]][ix]

    }
  }

  ### add label to p-value threshold line:
  pval.line.label <- line.label
  pval.x <- min(toptable[[contrast.fld]])

  ggplot(df, aes(x = contrast, y = -log10(pval), color = colorClass)) +
    geom_point_rast(size = 2) +
    xlab(xlabel) + ylab(ylabel) +
    scale_color_manual(name = 'Effect Direction', values = colors,
                       limits = c('Hypomethylated', 'Hypermethylated')) +
    geom_vline(xintercept = contrast.thresh * c(-1,1), linetype = 'dotdash') +
    geom_hline(yintercept = -log10(pthresh), linetype = 'dotted') +
    ggtitle(title) +
    geom_text_repel(aes(label=label), size = 4) +
    theme(axis.text = element_text(color = '#000000', size = 12),
          axis.title = element_text(color = 'black', size = 15),
          plot.title = element_text(size = 15, face = 'bold'),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(color = 'black', fill = NA),
          legend.position = c(0,0), legend.justification = c(0,0),
          legend.background = element_rect(colour = "black"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    annotate('text', x = pval.x, y = -log10(pthresh), label = pval.line.label,
             hjust = 0, vjust = 0, size = rel(4), fontface = 3)
  #guides(fill = guides_legend())
}

## volcano plot function - slightly modified for gene expression
### Label points with label.fld (if not NULL)
### use label.thresh and label.thresh.fld to determine what gets labeled.
volcano_plot_exprs <- function(toptable, contrast.fld = 'betaFC', pval.fld = 'P.Value',
                               pthresh = 0.01, contrast.thresh = 0.25,
                               xlabel = 'log2 fold change',
                               ylabel = '-log10 p-value',
                               title = '',
                               label.fld = 'Gene',
                               line.label='',
                               label.thresh = 0.05, label.thresh.fld = 'adj.P.Val',
                               colors = c(neut = '#8b8b8b', nonsig='grey', low = 'blue', high = 'red')){
  require(ggplot2);require(ggrepel);require(ggrastr)

  df <- data.frame(contrast = toptable[,contrast.fld],
                   pval = toptable[,pval.fld])
  df$colorClass <- 'nonsig'
  df$colorClass[df$pval < pthresh] <- 'neut'
  df$colorClass[df$pval < pthresh & df$contrast < -contrast.thresh] <- 'Down-Regulated'
  df$colorClass[df$pval < pthresh & df$contrast > contrast.thresh] <- 'Up-Regulated'
  names(colors)[3:4] <- c('Down-Regulated', 'Up-Regulated')

  ### labeling:
  df$label <- ''
  df$hjust <- as.numeric(df$contrast < 0)
  if(!is.null(label.fld)){
    ix <- which(toptable[[label.thresh.fld]] < label.thresh)
    df$label[ix] <- toptable[[label.fld]][ix]
  }

  ### add label to p-value threshold line:
  pval.line.label <- line.label
  pval.x <- min(toptable[[contrast.fld]])

  ggplot(df, aes(x = contrast, y = -log10(pval), color = colorClass)) +
    geom_point_rast(size = 2) +
    xlab(xlabel) + ylab(ylabel) +
    scale_color_manual(name = 'Effect Direction', values = colors,
                       limits = c('Down-Regulated', 'Up-Regulated')) +
    geom_vline(xintercept = contrast.thresh * c(-1,1), linetype = 'dotdash') +
    geom_hline(yintercept = -log10(pthresh), linetype = 'dotted') +
    ggtitle(title) +
    geom_text_repel(aes(label=label), size = 3.5) +
    theme(axis.text = element_text(color = '#000000', size = 12),
          axis.title = element_text(color = 'black', size = 15, face = 'bold'),
          plot.title = element_text(size = 15, face = 'bold'),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(color = 'black', fill = NA),
          legend.position = c(0,0), legend.justification = c(0,0),
          legend.background = element_rect(colour = "black"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    annotate('text', x = pval.x, y = -log10(pthresh), label = pval.line.label,
             hjust = 0, vjust = 0, size = rel(4), fontface = 3)
  #guides(fill = guides_legend())
}


## This is my own version of the DMR.plot, using ggplot
## plot of betas, colored by a class, subset of CpG sites according to a subset
## of features, fSub, which should contain the fields gene.fld, cpg.fld,
## pos.fld, and chr.fld
## attempt to plot medians
## make sure columns of betas correspond to rows of pheno
## To add: gene names and genomic locations (e.g. 5'UTR, Body, etc.)
gene.plot <- function(fSub, betas, pheno, groupVar = NULL,
                      gene.fld = 'GeneSymbol', cpg.fld = 'cpgSite',
                      pos.fld = 'pos', chr.fld = 'chr',
                      cols = c('red', 'blue')){
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  stopifnot(ncol(betas) == nrow(pheno),
            groupVar %in% names(pheno),
            all(c(gene.fld, cpg.fld, pos.fld, chr.fld) %in% names(fSub)))

  ### get CpG sites relating to gene, using features
  feat <- fSub[c(gene.fld, cpg.fld, pos.fld, chr.fld)]
  names(feat) <- c('Gene', 'CpG', 'POS', 'CHR')
  feat <- filter(feat, CpG %in% rownames(betas))
  stopifnot(n_distinct(feat$CHR) == 1)
  if(!is.null(groupVar)){
    plotDat <- cbind(as.data.frame(t(betas[feat$CpG, ])),
                     varb = pheno[[groupVar]])

    plotDat <- melt(plotDat, id.var = 'varb', variable.name = 'cgSite',
                    value.name = 'Beta')
    plotDat <- merge(plotDat, feat, all.x = T, by.x = 'cgSite', by.y = 'CpG')
    plotDat.med <- plotDat %>% group_by(cgSite, varb) %>%
      summarise(medianBeta = median(Beta), POS = unique(POS))

    d <- min(dist(unique(plotDat$POS)))
    ggplot(plotDat, aes(x = POS, y = Beta, color = varb, fill = varb)) +
      geom_point(alpha = 0.5,
                 position = position_jitterdodge(jitter.width = 0,
                                                 dodge.width = d*0.8)) +
      geom_line(aes(x = POS, y = medianBeta, color = varb), plotDat.med,
                position = position_jitterdodge(jitter.width = 0,
                                                dodge.width = d*0.8)) +
      scale_color_manual(values = cols, name = groupVar) +
      xlab(paste(unique(feat$CHR), 'position')) + ylab('Methylation beta') +
      pca.theme() +
      scale_fill_manual(values = cols, guide = 'none')
  } else{
    plotDat <- t(betas[feat$CpG, ])
    plotDat <- plotDat %>%
      melt(varnames = c('subj', 'cgSite'), value.name = 'Beta') %>%
      merge(feat, all.x = T, by.x = 'cgSite', by.y = 'CpG') %>%
      mutate(POSeven = factor(match(POS, sort(unique(POS))) %% 2))
    plotDat.med <- plotDat %>% group_by(cgSite) %>%
      summarise(medianBeta = median(Beta), POS = unique(POS))

    ggplot(plotDat, aes(x = POS, y = Beta, color = POSeven)) +
      geom_point(alpha = 0.5) +
      geom_line(aes(x = POS, y = medianBeta), plotDat.med, color = 'black') +
      scale_color_manual(values = cols) +
      xlab(paste(unique(feat$CHR), 'position')) + ylab('Methylation beta') +
      pca.theme() +
      theme(legend.position = 'none') +
      geom_hline(yintercept = c(-1,0,1), color = 'grey20') +
      ylim(range(plotDat$Beta))

  }
  ###############
  ####
  ###!!!!!!!!!! TO BE COMPLETED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}





### This is a wrapper for gene.plot (above) which takes in either a region e.g.
### 'chr5:122433740-122435550' like in DMRcate results OR a gene symbol
### ... other arguments are passed to gene.plot (betas, pheno, groupVar)
### if both are provided, region overrides gene
my.DMR.plot <- function(region = NULL, gene = NULL, f = fDat, ...){
  stopifnot(!is.null(region) | !is.null(gene))
  if(!is.null(region)){
    chrom <- gsub(':.*', '', region)
    sta.pos <- as.numeric(gsub('.*:', '', gsub('-.*', '', region)))
    sto.pos <- as.numeric(gsub('.*-', '', region))
    f <- f %>% filter(chr == chrom, pos <= sto.pos, pos >= sta.pos)
  } else{
    f <- filter(f, grepl(paste0('(^|;)', gene, '(;|$)'), GeneSymbol))

  }
  gene.plot(f, ...)

}


## 2016/07/05: modify gene.plot by adding UCSC_RefGene_group annotation regions
gene.plot <- function(fSub, betas, pheno, groupVar = NULL,
                      gene.fld = 'GeneSymbol', cpg.fld = 'cpgSite',
                      pos.fld = 'pos', chr.fld = 'chr',
                      gene.2.fld = 'UCSC_RefGene_Name', annot.2.fld = 'UCSC_RefGene_Group',
                      group.cols = c('red', 'blue'),
                      annot.cols = terrain.colors(12)[seq(1, 11, 2)]){
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  stopifnot(ncol(betas) == nrow(pheno),
            groupVar %in% names(pheno),
            all(c(gene.fld, cpg.fld, pos.fld, chr.fld) %in% names(fSub)))

  ### get CpG sites relating to gene, using features
  feat <- fSub[c(gene.fld, cpg.fld, pos.fld, chr.fld, gene.2.fld, annot.2.fld)]
  names(feat) <- c('Gene', 'CpG', 'POS', 'CHR', 'Gene2', 'Annot2')
  feat <- filter(feat, CpG %in% rownames(betas))
  stopifnot(n_distinct(feat$CHR) == 1)



  if(!is.null(groupVar)){
    plotDat <- cbind(as.data.frame(t(betas[feat$CpG, ])),
                     varb = pheno[[groupVar]])

    plotDat <- gather(plotDat, 'cgSite','Beta', -varb)
    #     plotDat <- melt(plotDat, id.var = 'varb', variable.name = 'cgSite',
    #                     value.name = 'Beta')
    plotDat <- merge(plotDat, feat, all.x = T, by.x = 'cgSite', by.y = 'CpG')
    plotDat.med <- plotDat %>% group_by(cgSite, varb) %>%
      summarise(medianBeta = median(Beta), POS = unique(POS))


    d <- min(dist(unique(plotDat$POS)))
    G <- ggplot(plotDat, aes(x = POS, y = Beta)) +
      geom_point(alpha = 0.5, mapping = aes(color = varb),
                 position = position_jitterdodge(jitter.width = 0,
                                                 dodge.width = d*0.8)) +
      geom_line(aes(x = POS, y = medianBeta, color = varb), plotDat.med,
                position = position_jitterdodge(jitter.width = 0,
                                                dodge.width = d*0.8)) +
      scale_color_manual(values = group.cols, name = groupVar) +
      xlab(paste(unique(feat$CHR), 'position')) + ylab('Methylation beta') +
      pca.theme()

  } else{
    plotDat <- t(betas[feat$CpG, ])
    plotDat <- plotDat %>%
      melt(varnames = c('subj', 'cgSite'), value.name = 'Beta') %>%
      merge(feat, all.x = T, by.x = 'cgSite', by.y = 'CpG') %>%
      mutate(POSeven = factor(match(POS, sort(unique(POS))) %% 2))
    plotDat.med <- plotDat %>% group_by(cgSite) %>%
      summarise(medianBeta = median(Beta), POS = unique(POS))

    G <- ggplot(plotDat, aes(x = POS, y = Beta)) +
      geom_point(aes(color = POSeven), alpha = 0.5) +
      geom_line(aes(x = POS, y = medianBeta), plotDat.med, color = 'black') +
      scale_color_manual(values = group.cols, guide = 'none') +
      xlab(paste(unique(feat$CHR), 'position')) + ylab('Methylation beta') +
      pca.theme() +
      geom_hline(yintercept = 0, color = 'grey20', linetype = 'dashed', size = 0.5) #+
    #ylim(range(plotDat$Beta))

  }
  ## GENE ANNOTATION:
  geneDat <- filter(feat, Gene != '')
  if(nrow(geneDat) > 0){
    geneDat <- geneDat %>%
      group_by(CpG) %>%
      do(data.frame(POS = .$POS, Gene2 = strsplit(.$Gene2, ';')[[1]],
                    Annot2 = strsplit(.$Annot2, ';')[[1]])) %>%
      ungroup %>% distinct() %>%
      mutate(Annot2 = factor(Annot2, levels = c('TSS1500', 'TSS200', '1stExon',
                                                "5'UTR", "Body", "3'UTR")))
    gene.pos <- group_by(geneDat, Gene2) %>%
      summarise(minPos = min(POS), maxPos = max(POS)) %>%
      arrange(minPos)

    overlappy <- arrange(gene.pos, minPos)
    negs <- -1
    n <- 0

    while(!is.na(negs) && any(negs < 0, na.rm = T)){
      n <- n + 1
      negs <- with(overlappy, lead(minPos, n = n) - maxPos)
    }

    minmax <- range(plotDat$Beta)
    incr <- diff(minmax)*0.1
    ymax <- minmax[2] + incr

    ## add x-pos, xend, xlabel pos, ypos ylabelpos
    gene.pos <- mutate(gene.pos, xlab = (maxPos + minPos)/2,
                       yposBase = rep(ymax + incr*(1:n), nrow(gene.pos))[1:nrow(gene.pos)],
                       ylab = yposBase - incr/4)
    ## merge this to geneDat for positions of annotations
    geneDat <- left_join(geneDat, gene.pos, by = 'Gene2') %>%
      mutate(ypos = yposBase + incr*(as.numeric(Annot2)/9 - 1/18))

    G <- G + #geom_segment(aes(x = minPos, xend = maxPos, y = yposBase, yend = yposBase),
      #                   gene.pos, color = 'black') +
      #geom_text(aes(x = xlab, y = ylab, label = Gene2), gene.pos, color = 'black',
      #          size = 2) +
      #geom_point(aes(x = POS, y = ypos, fill = Annot2), geneDat, shape = 22, color = 'transparent', size = 3) +
      scale_fill_manual(values = annot.cols, name = 'UCSC Annotation') +
      guides(fill = guide_legend(override.aes = list(size = 5)))
  }
  G

}




## plot for a DMRcate result object:
## volcano plot but have variable size per point, based on number of probes in gene
## volcano plot function
### Label points with label.fld (if not NULL)
### use label.thresh and label.thresh.fld to determine what gets labeled.
### automatically replaces p-values of 0 with 1e-300
volcano_plot.DMRcate1 <- function(res.obj, pval.fld = c('mean', 'min'),
                                  pthresh = 0.01, contrast.thresh = 0.25, min.p = 1e-100,
                                  xlabel = 'max beta difference',
                                  ylabel = '-log10 p-value',
                                  title = '',
                                  label.thresh = 0.05,
                                  colors = c(neut = 'darkgrey', low = 'red', high = 'blue')){
  require(ggplot2)

  stopifnot(is(res.obj, 'dmrcate.output'))
  p.fld <- match.arg(pval.fld)
  df <- res.obj$results[c('gene_assoc', paste0(p.fld, 'pval'), 'maxbetafc', 'no.probes')]
  names(df) <- c('Gene', 'pval', 'contrast', 'no.probes')
  df$colorClass <- 'neut'
  df$colorClass[df$pval < pthresh & df$contrast < -contrast.thresh] <- 'hypomethylated'
  df$colorClass[df$pval < pthresh & df$contrast > contrast.thresh] <- 'hypermethylated'
  names(colors)[2:3] <- c('hypomethylated', 'hypermethylated')

  ### labeling:
  df$labl <- ''
  df$hjust <- as.numeric(df$contrast < 0)

  ix <- which(df$pval < label.thresh)
  df$labl[ix] <- gsub(';.*', '', df$Gene[ix])

  ## add p-value to label for points below min.p
  ix <- which(df$pval < min.p)
  df$labl[ix] <- sprintf('%s(%.0e)', df$labl[ix], df$pval[ix])
  # set low p-values to min.p + jitter
  ix <- ix[order(df$pval[ix])]
  df$pval[ix] <- 10^(log10(min.p) + sort(rnorm(length(ix), 0, 3)))

  ### add label to p-value threshold line:
  pval.line.label <- paste('p =', pthresh)
  pval.x <- min(df$contrast)

  ## for x-limit:
  xmax <- max(abs(df$contrast))*1.2

  ggplot(df, aes(x = contrast, y = -log10(pval), color = colorClass, size = no.probes)) +
    geom_point(alpha = 0.33) +
    xlab(xlabel) + ylab(ylabel) +
    ylim(c(0, max(-log10(df$pval)) + 1)) + xlim(c(-xmax, xmax)) +
    scale_color_manual(name = 'effect direction', values = colors,
                       limits = c('hypomethylated', 'hypermethylated')) +
    scale_size_area(name = '#probes', max_size = 15) +
    geom_vline(xintercept = contrast.thresh * c(-1,1), linetype = 'dotdash') +
    geom_hline(yintercept = -log10(pthresh), linetype = 'dotted') +
    ggtitle(title) +
    geom_text(aes(x = contrast, y = -log10(pval), label = labl, hjust = hjust),
              data = df, vjust = 0, color = 'black', size = 3) +
    theme(axis.text = element_text(color = '#000000', size = 12),
          axis.title = element_text(color = 'black', size = 15, face = 'bold'),
          plot.title = element_text(size = 15, face = 'bold'),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(color = 'black', fill = NA),
          legend.position = 'right',
          legend.background = element_rect(colour = "black"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    annotate('text', x = pval.x, y = -log10(pthresh), label = pval.line.label,
             hjust = 0, vjust = 0, size = rel(4), fontface = 3)
  #guides(fill = guides_legend())
}

# slightly modified from Nick's code -  label region with >= 3 cpgs
volcano_plot.DMRcate2 <- function(res.obj,
                                  pthresh = 0.001, contrast.thresh = 0.5,
                                  xlabel = 'max beta difference',
                                  ylabel = '- log10 min FDR',
                                  title = '',
                                  colors = c(neut = 'darkgrey', low = 'blue', high = 'red')){
  require(ggplot2)
  df <- res.obj[,c('Gene.Symbols', 'minfdr', 'Stouffer', 'maxbetafc', 'no.cpgs')]
  names(df) <- c('Gene', 'minfdr', 'Stouffer', 'contrast', 'no.cpgs')
  df$colorClass <- 'neut'
  df$colorClass[df$minfdr < pthresh & df$contrast < -contrast.thresh] <- 'Hypomethylated'
  df$colorClass[df$minfdr < pthresh & df$contrast > contrast.thresh] <- 'Hypermethylated'
  names(colors)[2:3] <- c('Hypomethylated', 'Hypermethylated')

  ### labeling:
  df$labl <- ''
  df$hjust <- as.numeric(df$contrast < 0)

  ix <- which(df$no.cpgs > 2 & abs(df$contrast) > contrast.thresh)
  df$labl[ix] <- gsub(';.*', '', df$Gene[ix])

  ### add label to p-value threshold line:
  pval.line.label <- paste('min FDR =', pthresh)
  pval.x <- min(df$contrast)

  ## for x-limit:
  xmax <- max(abs(df$contrast))*1.2

  ggplot(df, aes(x = contrast, y = -log10(minfdr), color = colorClass, size = no.cpgs)) +
    geom_point(alpha = 0.33) +
    xlab(xlabel) + ylab(ylabel) +
    ylim(c(0, max(-log10(df$minfdr)) + 1)) + xlim(c(-xmax, xmax)) +
    scale_color_manual(name = 'Effect Direction', values = colors,
                       limits = c('Hypomethylated', 'Hypermethylated')) +
    scale_size_area(name = 'Number of Probes', max_size = 15, breaks = c(1,4,7,10,13)) +
    geom_vline(xintercept = contrast.thresh * c(-1,1), linetype = 'dotdash') +
    #geom_hline(yintercept = -log10(pthresh), linetype = 'dotted') +
    ggtitle(title) +
    geom_text_repel(aes(label=labl), size = 7, point.padding = unit(0.7, 'lines'),min.segment.length = unit(0.1, 'lines'))+
    theme_bw()+
    theme(axis.text = element_text(color = '#000000', size = 12.5),
          axis.title = element_text(color = 'black', size = 15, face = 'bold'))#+
  #annotate('text', x = pval.x, y = -log10(pthresh), label = pval.line.label,
  #         hjust = 0, vjust = 0, size = 5, fontface = 3)
}

##### Manhattan-like plot for DMRcate results... will plot p-value vs. location
# dmrcate.manhattan = function(res.obj, title=NULL, cols = c('red', 'skyblue'),
#                              pval.fld = c('mean', 'min'),
#                              chrs = 1:22, ymax = NULL,
#                              annotate=F, geneList=NULL) {
#   library(ggplot2)
#   library(dplyr)
#
#   stopifnot(is(res.obj, 'dmrcate.output'))
#   p.fld <- match.arg(pval.fld)
#   p.fld <- paste0(p.fld, 'pval')
#
#   res <- res.obj$results
#
#   if (annotate & is.null(geneList))
#     stop("You requested annotation but provided no SNPlist!")
#   if(annotate & any(! geneList %in% res$gene_assoc))
#     stop('No SNPs in SNPlist are found in the dataframe')
#
#   df <- res %>%
#     mutate(chr = as.numeric(gsub(':.*', '', gsub('chr', '', hg19coord))),
#            pos1 = as.numeric(gsub('.*:', '', gsub('-.*', '', hg19coord))),
#            pos2 = as.numeric(gsub('.*-', '', hg19coord))) %>%
#     mutate(pos = (pos1 + pos2)/2) %>%
#     select(gene = gene_assoc, pval = one_of(p.fld), chr, pos1, pos2, pos)
#
#   ### replace p=0 with something small:
#   df$pval[df$pval == 0] <- 1e-300
#
#   #limit to only chrs
#   df <- df %>% tbl_df %>% filter(chr %in% chrs) %>%
#     mutate(logp = -log10(pval)) %>% arrange(chr, pos)
#   df$chr <- factor(df$chr)
#
#   # get position references for each CHR:
#   posDat <- df %>% group_by(chr) %>%
#     summarise(minBP = min(pos1), maxBP = max(pos2)) %>%
#     mutate(BPrange = maxBP - minBP)
#
#   tmp <- data.frame(0, 0, 0, 0); names(tmp) <- names(posDat)
#   posDat <- rbind(tmp, posDat)
#   posDat$cumBP <- sapply(1:nrow(posDat), function(x) sum(posDat$BPrange[1:x]))
#   posDat <- data.frame(chr = posDat$chr[-1], addBP = posDat$cumBP[-nrow(posDat)])
#
#   df <- merge(df, posDat, all.x = T)
#
#   ## set plotting position for each chromosome based on BP and previous CHR BP:
#   df <- df %>% mutate(start.pos = pos1 + addBP, stop.pos = pos2 + addBP, BP = pos + addBP)
#
#   ## get tick location for each chr:
#   ticks <- df %>% group_by(chr) %>% summarise(x = median(BP))
#
#   mycols <- rep(cols, length(chrs))[1:length(unique(chrs))]
#
#   if(length(unique(chrs)) == 1){
#     plot <- ggplot(df, aes(x = BP, y = logp, color = chr)) +
#       xlab(paste('Chromosome', unique(chrs)))
#   } else {
#     plot <- ggplot(df, aes(x = BP, y = logp, color = chr)) +
#       scale_x_continuous(name = 'Chromosome', breaks = ticks$x, labels = ticks$chr)
#   }
#   plot <- plot + geom_point() +
#     geom_segment(aes(x = start.pos, xend = stop.pos, y = logp, yend = logp)) +
#     scale_color_manual(values = mycols) +
#     scale_y_continuous(breaks = seq(1, ceiling(max(df$logp)), by = 2)) +
#     ylab(expression(-log[10](italic(p)))) +
#     ggtitle(title) +
#     theme_bw() +
#     theme(axis.text.y = element_text(color = 'black'),
#           axis.text.x = element_text(angle = -20, color = 'black'),
#           legend.position = 'none',
#           panel.grid.minor = element_blank())
#   if(annotate){
#     df$labels <- ''
#     df$labels[match(geneList, df$gene)] <- geneList
#     plot <- plot + geom_text(aes(x = BP, y = logp, label = labels), df,
#                              vjust = 0, hjust = 0.5, color = 'black')
#   }
# #   if (!is.null(suggestiveline))
# #     plot <- plot + geom_hline(yintercept = suggestiveline, colour="brown", alpha = 0.33)
# #   if (!is.null(genomewideline))
# #     plot <- plot + geom_hline(yintercept = genomewideline, colour = "brown")
#   if(!is.null(ymax) && is.numeric(ymax)){
#     if(ymax <= max(d$logp) & ymax >= 0)
#       plot <- plot + ylim(0, ymax)
#   }
#   plot
# }

### plot density of no.probes for DMRcate results
### add information about #annotated regions...
noProbes.summPlot <- function(res.obj){
  library(ggplot2)
  library(gridExtra)
  stopifnot(is(res.obj, 'dmrcate.output'))
  res <- res.obj$results
  res$annotated <- c('annotated', 'unannotated')[(res$gene_assoc == '') + 1]
  tab <- table(res$annotated)
  tab <- tableGrob(matrix(tab, nrow = 1, dimnames = list(NULL, names(tab))))
  r <- density(log10(res$no.probes))
  ggplot(res, aes(no.probes, fill = annotated)) +
    scale_x_log10() + geom_density(alpha = 0.5) + theme_bw() +
    theme(legend.position = c(1, 1), legend.justification = c(1,1)) +
    annotation_custom(tab, median(log10(res$no.probes)), log10(max(res$no.probes)),
                      max(r$y)/2, max(r$y)/2)
}

summarize.tt.cutoffs <- function(tt, pval.thresh, p.fld = c('P.Value', 'adj.P.Val'),
                                 beta.thresh = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)){
  library(dplyr)
  library(reshape2)
  p.fld <- match.arg(p.fld)
  stopifnot(is.data.frame(tt),
            all(c('betaFC', p.fld, 'Gene') %in% names(tt)))
  tt <- dplyr::select(tt, pval = one_of(p.fld), betaFC, Gene)
  df <- data.frame(P = rep(sort(pval.thresh), each = length(beta.thresh)),
                   B = rep(sort(beta.thresh, decreasing = T), length(pval.thresh)))
  df2 <- df %>% group_by(P, B) %>%
    summarise(N.cpg = nrow(filter(tt, pval < P, abs(betaFC) > B)),
              N.gene = length(unique(unlist(strsplit(filter(tt, pval < P, abs(betaFC) > B, Gene != '')$Gene, ';'))))) %>%
    ungroup
  out <- list(cpgSites = acast(df2, P ~ B, value.var = 'N.cpg'),
              genes = acast(df2, P ~ B, value.var = 'N.gene'))
  lapply(out, function(x){
    names(dimnames(x)) <- c('pvalue cutoff', 'beta cutoff')
    x
  })
}


## function which plots a stacked bar plot of CpG site annotated features
## limma.res: limma top-table with cg sites in either rownames or in 'cgSite' field
## annoDat: annotation data.frame which needs at least 'UCSC_RefGene_Group' and
##         'Relation_to_Island' fields, and 'Name' for the CpG sites
##feature: choose which feature to plot
## if we're choosing gene groups, we first restrict our attention to CpG sites
##   annotated to genes
limma450k_stackBar <- function(limma.res, annoDat,
                               feature = c('CpG.Density', 'GeneGroup'),
                               p.fld = 'adj.P.Val',
                               p.thresh = 0.1){
  library(ggplot2)
  library(dplyr)
  stopifnot(all(c('Relation_to_Island', 'Name', 'UCSC_RefGene_Group') %in%
                  names(annoDat)))
  stopifnot(all(c(p.fld, 'betaFC') %in% names(limma.res)),
            p.thresh < 1, p.thresh > 0)
  cgSites <- rownames(limma.res)
  if('cgSite' %in% names(limma.res))
    cgSites <- limma.res$cgSite
  if('cpgSite' %in% names(limma.res))
    cgSites <- limma.res$cpgSite
  if('CpG' %in% names(limma.res))
    cgSites <- limma.res$CpG

  limma.res$cgSite <- cgSites
  feature <- match.arg(feature)

  limma.res$P <- limma.res[[p.fld]]

  limma.res.anno = merge(limma.res, annoDat, by.x = 'cgSite', by.y = 'Name', all.x = TRUE)
  limma.res.anno$UCSC_RefGene_Group[limma.res.anno$UCSC_RefGene_Group == ''] = 'Not Assigned'
  #  if(feature == 'GeneGroup'){
  #    tmp <- filter(annoDat, UCSC_RefGene_Group != '')
  #    limma.res <- filter(limma.res, cgSite %in% tmp$Name)
  #  }

  ## significant hits:
  cg.sig <- filter(limma.res, P < p.thresh) %>% with(cgSite)
  cg.sig.1 <- filter(limma.res, P < p.thresh, abs(betaFC) > 0.1) %>% with(cgSite)
  cg.sig.hypo <- filter(limma.res, P < p.thresh, betaFC < -0.1) %>% with(cgSite)
  cg.sig.hyper <- filter(limma.res, P < p.thresh, betaFC > 0.1) %>% with(cgSite)

  cg.list <- list('Total' = limma.res$cgSite,
                  'differential methylation' = cg.sig, '> 10% differential methylation' = cg.sig.1,
                  '> 10% hypomethylation' = cg.sig.hypo,
                  '> 10% hypermethylation' = cg.sig.hyper)

  annoDat <- filter(annoDat, Name %in% limma.res$cgSite)
  summDat <- plyr::ldply(cg.list, length) %>%
    dplyr::rename(CpG.class = .id, n = V1) %>%
    mutate(n = sprintf('n=%d', n))

  if(feature == 'CpG.Density'){
    ## transform Relation_to_Island to factor:
    trans <- c(Island = 'Island', N_Shore = 'Shore', S_Shore = 'Shore',
               N_Shelf = 'Shelf', S_Shelf = 'Shelf', OpenSea = 'OpenSea')
    annoDat <- mutate(annoDat,
                      CpG.density = factor(trans[Relation_to_Island],
                                           levels = c('Island', 'Shore',
                                                      'Shelf', 'OpenSea')))

    plotDat <- plyr::ldply(cg.list, function(x)
      filter(annoDat, Name %in% x) %>% dplyr::select(CpG.density)) %>%
      dplyr::rename(CpG.class = .id)

    plotDat$CpG.class = factor(plotDat$CpG.class, levels = c('Total',
                                                             'differential methylation', '> 10% differential methylation',
                                                             '> 10% hypomethylation',
                                                             '> 10% hypermethylation'))
    # plot object
    g <- ggplot(plotDat, aes(x = CpG.class)) +
      geom_bar(aes(fill = CpG.density),  color = 'black', position = 'fill') +
      scale_fill_brewer(name = 'CpG Density', type = 'seq',palette = 4)+
      scale_y_continuous(labels = scales::percent)+
      theme_bw() + ylab('% probes') + xlab('CpG class') +
      #geom_text(aes(y = 1.05, label = n), summDat, vjust = 0) +
      #theme(axis.text.x = element_text(angle = 20, hjust = 1))
      scale_x_discrete(labels=c('Total' = 'Total', 'differential methylation' = 'differential\nmethylation',
                                '> 10% differential methylation' = '> 10% differential\nmethylation',
                                '> 10% hypomethylation'='> 10% hypomethylation',
                                '> 10% hypermethylation'='> 10% hypermethylation'))

    # table object:
    tab <- group_by(plotDat, CpG.class) %>%
      do(broom::tidy(addmargins(table(.$CpG.density)))) %>%
      tidyr::spread(Var1, Freq)

  }
  if(feature == 'GeneGroup'){
    helper <- function(vec){
      unlist(lapply(strsplit(vec, ';'), unique))
    }
    plotDat <- lapply(cg.list, function(x)
      filter(limma.res.anno, cgSite %in% x) %>% with(helper(UCSC_RefGene_Group)))
    plotDat <- plyr::ldply(plotDat, function(x) data.frame(Gene.Feature = x)) %>%
      dplyr::rename(CpG.class = .id)
    plotDat$Gene.Feature <- gsub('1stExon', '1st Exon', plotDat$Gene.Feature)
    plotDat <- plotDat %>% mutate(Gene.Feature = factor(Gene.Feature,
                                                        levels = c('Not Assigned',"3'UTR", "Body",
                                                                   "5'UTR", '1st Exon', 'TSS200', 'TSS1500')))
    plotDat$CpG.class = factor(plotDat$CpG.class, levels = c('Total',
                                                             'differential methylation', '> 10% differential methylation',
                                                             '> 10% hypomethylation',
                                                             '> 10% hypermethylation'))
    # plot object:
    g <- ggplot(plotDat, aes(x = CpG.class)) +
      geom_bar(aes(fill = Gene.Feature), color = 'black', position = 'fill') +
      scale_fill_brewer(name = 'Gene Feature',type = 'seq',palette = 3)+
      scale_y_continuous(labels = scales::percent)+
      theme_bw() + ylab('% probes') +  xlab('CpG class') +
      #geom_text(aes(y = 1.05, label = n), summDat, vjust = 0) +
      #theme(axis.text.x = element_text())+
      scale_x_discrete(labels=c('Total' = 'Total', 'differential methylation' = 'differential\nmethylation',
                                '> 10% differential methylation' = '> 10% differential\nmethylation',
                                '> 10% hypomethylation'='> 10% hypomethylation',
                                '> 10% hypermethylation'='> 10% hypermethylation'))
    #guides(fill = guide_legend(reverse = T))


    # table object:
    tab <- group_by(plotDat, CpG.class) %>%
      do(broom::tidy(addmargins(table(.$Gene.Feature)))) %>%
      tidyr::spread(Var1, Freq)

  }
  return(list(plot = g, table = tab))
}

## divide the function into two parts upon Rachel's requests
limma450k_stackBar_part1=function(limma.res, annoDat,
                                  feature = c('CpG.Density', 'GeneGroup'),
                                  p.fld = 'adj.P.Val',
                                  p.thresh = 0.05){
  library(ggplot2)
  library(dplyr)
  stopifnot(all(c('Relation_to_Island', 'Name', 'UCSC_RefGene_Group') %in%
                  names(annoDat)))
  stopifnot(all(c(p.fld, 'betaFC') %in% names(limma.res)),
            p.thresh < 1, p.thresh > 0)
  cgSites <- rownames(limma.res)
  if('cgSite' %in% names(limma.res))
    cgSites <- limma.res$cgSite
  if('cpgSite' %in% names(limma.res))
    cgSites <- limma.res$cpgSite
  limma.res$cgSite <- cgSites
  feature <- match.arg(feature)

  limma.res$P <- limma.res[[p.fld]]

  limma.res.anno = merge(limma.res, annoDat, by.x = 'cgSite', by.y = 'Name', all.x = TRUE)
  limma.res.anno$UCSC_RefGene_Group[limma.res.anno$UCSC_RefGene_Group == ''] = 'Not Assigned'

  ## significant hits:
  cg.sig <- filter(limma.res, P < p.thresh) %>% with(cgSite)
  cg.sig.5 <- filter(limma.res, P < p.thresh, abs(betaFC) > 0.5) %>% with(cgSite)

  cg.list <- list('Total' = limma.res$cgSite,
                  'differential methylation' = cg.sig,
                  '> 50% differential methylation' = cg.sig.5)

  annoDat <- filter(annoDat, Name %in% limma.res$cgSite)
  summDat <- plyr::ldply(cg.list, length) %>%
    dplyr::rename(CpG.class = .id, n = V1) %>%
    mutate(n = sprintf('n=%d', n))

  if(feature == 'CpG.Density'){
    ## transform Relation_to_Island to factor:
    trans <- c(Island = 'Island', N_Shore = 'N_Shore', S_Shore = 'S_Shore',
               N_Shelf = 'N_Shelf', S_Shelf = 'S_Shelf', OpenSea = 'OpenSea')
    annoDat <- mutate(annoDat,
                      CpG.density = factor(trans[Relation_to_Island],
                                           levels = c('Island','N_Shore', 'S_Shore',
                                                      'N_Shelf','S_Shelf', 'OpenSea')))

    plotDat <- plyr::ldply(cg.list, function(x)
      filter(annoDat, Name %in% x) %>% dplyr::select(CpG.density)) %>%
      dplyr::rename(CpG.class = .id)

    plotDat$CpG.class = factor(plotDat$CpG.class, levels = c('Total',
                                                             'differential methylation',
                                                             '> 50% differential methylation'))
    # plot object
    g <- ggplot(plotDat, aes(x = CpG.class)) +
      geom_bar(aes(fill = CpG.density),  color = 'black', position = 'fill') +
      scale_fill_brewer(name = 'CpG Density', type = 'seq',palette = 4)+
      scale_y_continuous(labels = scales::percent)+
      theme_bw() + ylab('% probes') + xlab('CpG class') +
      #geom_text(aes(y = 1.05, label = n), summDat, vjust = 0) +
      #theme(axis.text.x = element_text(angle = 20, hjust = 1))
      scale_x_discrete(labels=c('Total' = 'Total', 'differential methylation' = 'differential\nmethylation',
                                '> 50% differential methylation' = '> 50% differential\nmethylation'))

    # table object:
    tab <- group_by(plotDat, CpG.class) %>%
      do(broom::tidy(addmargins(table(.$CpG.density)))) %>%
      tidyr::spread(Var1, n)

  }
  if(feature == 'GeneGroup'){
    helper <- function(vec){
      unlist(lapply(strsplit(vec, ';'), unique))
    }
    plotDat <- lapply(cg.list, function(x)
      filter(limma.res.anno, cgSite %in% x) %>% with(helper(UCSC_RefGene_Group)))
    plotDat <- plyr::ldply(plotDat, function(x) data.frame(Gene.Feature = x)) %>%
      dplyr::rename(CpG.class = .id)
    plotDat$Gene.Feature <- gsub('1stExon', '1st Exon', plotDat$Gene.Feature)
    plotDat <- plotDat %>% mutate(Gene.Feature = factor(Gene.Feature,
                                                        levels = c('Not Assigned',"3'UTR", "Body",
                                                                   "5'UTR", '1st Exon', 'TSS200', 'TSS1500')))
    plotDat$CpG.class = factor(plotDat$CpG.class, levels = c('Total',
                                                             'differential methylation',
                                                             '> 50% differential methylation'))
    # plot object:
    g <- ggplot(plotDat, aes(x = CpG.class)) +
      geom_bar(aes(fill = Gene.Feature), color = 'black', position = 'fill') +
      scale_fill_brewer(name = 'Gene Feature',type = 'seq',palette = 3)+
      scale_y_continuous(labels = scales::percent)+
      theme_bw() + ylab('% probes') +  xlab('CpG class') +
      #geom_text(aes(y = 1.05, label = n), summDat, vjust = 0) +
      #theme(axis.text.x = element_text())+
      scale_x_discrete(labels=c('Total' = 'Total', 'differential methylation' = 'differential\nmethylation',
                                '> 50% differential methylation' = '> 50% differential\nmethylation'))
    #guides(fill = guide_legend(reverse = T))


    # table object:
    tab <- group_by(plotDat, CpG.class) %>%
      do(broom::tidy(addmargins(table(.$Gene.Feature)))) %>%
      tidyr::spread(Var1, n)

  }
  return(list(plot = g, table = tab))
}



limma450k_stackBar_part2=function(limma.res, annoDat,
                                  feature = c('CpG.Density', 'GeneGroup'),
                                  p.fld = 'adj.P.Val',
                                  p.thresh = 0.05){
  library(ggplot2)
  library(dplyr)
  stopifnot(all(c('Relation_to_Island', 'Name', 'UCSC_RefGene_Group') %in%
                  names(annoDat)))
  stopifnot(all(c(p.fld, 'betaFC') %in% names(limma.res)),
            p.thresh < 1, p.thresh > 0)
  cgSites <- rownames(limma.res)
  if('cgSite' %in% names(limma.res))
    cgSites <- limma.res$cgSite
  if('cpgSite' %in% names(limma.res))
    cgSites <- limma.res$cpgSite
  limma.res$cgSite <- cgSites
  feature <- match.arg(feature)

  limma.res$P <- limma.res[[p.fld]]

  limma.res.anno = merge(limma.res, annoDat, by.x = 'cgSite', by.y = 'Name', all.x = TRUE)
  limma.res.anno$UCSC_RefGene_Group[limma.res.anno$UCSC_RefGene_Group == ''] = 'Not Assigned'
  #  if(feature == 'GeneGroup'){
  #    tmp <- filter(annoDat, UCSC_RefGene_Group != '')
  #    limma.res <- filter(limma.res, cgSite %in% tmp$Name)
  #  }

  ## significant hits:
  cg.sig <- filter(limma.res, P < p.thresh) %>% with(cgSite)
  cg.sig.hypo <- filter(limma.res, P < p.thresh, betaFC < 0) %>% with(cgSite)
  cg.sig.hyper <- filter(limma.res, P < p.thresh, betaFC > 0) %>% with(cgSite)

  cg.list <- list('Total' = limma.res$cgSite,
                  'differential methylation' = cg.sig,
                  'decrease methylation' = cg.sig.hypo,
                  'increase methylation' = cg.sig.hyper)

  annoDat <- filter(annoDat, Name %in% limma.res$cgSite)
  summDat <- plyr::ldply(cg.list, length) %>%
    dplyr::rename(CpG.class = .id, n = V1) %>%
    mutate(n = sprintf('n=%d', n))

  if(feature == 'CpG.Density'){
    ## transform Relation_to_Island to factor:
    trans <- c(Island = 'Island', N_Shore = 'N_Shore', S_Shore = 'S_Shore',
               N_Shelf = 'N_Shelf', S_Shelf = 'S_Shelf', OpenSea = 'OpenSea')
    annoDat <- mutate(annoDat,
                      CpG.density = factor(trans[Relation_to_Island],
                                           levels = c('Island','N_Shore', 'S_Shore',
                                                      'N_Shelf','S_Shelf', 'OpenSea')))

    plotDat <- plyr::ldply(cg.list, function(x)
      filter(annoDat, Name %in% x) %>% dplyr::select(CpG.density)) %>%
      dplyr::rename(CpG.class = .id)

    plotDat$CpG.class = factor(plotDat$CpG.class, levels = c('Total',
                                                             'differential methylation',
                                                             'decrease methylation',
                                                             'increase methylation'))
    # plot object
    g <- ggplot(plotDat, aes(x = CpG.class)) +
      geom_bar(aes(fill = CpG.density),  color = 'black', position = 'fill') +
      scale_fill_brewer(name = 'CpG Density', type = 'seq',palette = 4)+
      scale_y_continuous(labels = scales::percent)+
      theme_bw() + ylab('% probes') + xlab('CpG class') +
      #geom_text(aes(y = 1.05, label = n), summDat, vjust = 0) +
      #theme(axis.text.x = element_text(angle = 20, hjust = 1))
      scale_x_discrete(labels=c('Total' = 'Total',
                                'differential methylation' = 'differential\nmethylation',
                                'decrease methylation'='decrease\nmethylation',
                                'increase methylation'='increase\nmethylation'))

    # table object:
    tab <- group_by(plotDat, CpG.class) %>%
      do(broom::tidy(addmargins(table(.$CpG.density)))) %>%
      tidyr::spread(Var1, n)

  }
  if(feature == 'GeneGroup'){
    helper <- function(vec){
      unlist(lapply(strsplit(vec, ';'), unique))
    }
    plotDat <- lapply(cg.list, function(x)
      filter(limma.res.anno, cgSite %in% x) %>% with(helper(UCSC_RefGene_Group)))
    plotDat <- plyr::ldply(plotDat, function(x) data.frame(Gene.Feature = x)) %>%
      dplyr::rename(CpG.class = .id)
    plotDat$Gene.Feature <- gsub('1stExon', '1st Exon', plotDat$Gene.Feature)
    plotDat <- plotDat %>% mutate(Gene.Feature = factor(Gene.Feature,
                                                        levels = c('Not Assigned',"3'UTR", "Body",
                                                                   "5'UTR", '1st Exon', 'TSS200', 'TSS1500')))
    plotDat$CpG.class = factor(plotDat$CpG.class, levels = c('Total',
                                                             'differential methylation',
                                                             'decrease methylation',
                                                             'increase methylation'))
    # plot object:
    g <- ggplot(plotDat, aes(x = CpG.class)) +
      geom_bar(aes(fill = Gene.Feature), color = 'black', position = 'fill') +
      scale_fill_brewer(name = 'Gene Feature',type = 'seq',palette = 3)+
      scale_y_continuous(labels = scales::percent)+
      theme_bw() + ylab('% probes') +  xlab('CpG class') +
      #geom_text(aes(y = 1.05, label = n), summDat, vjust = 0) +
      #theme(axis.text.x = element_text())+
      scale_x_discrete(labels=c('Total' = 'Total',
                                'differential methylation' = 'differential\nmethylation',
                                'decrease methylation'='decrease\nmethylation',
                                'increase methylation'='increase\nmethylation'))
    #guides(fill = guide_legend(reverse = T))


    # table object:
    tab <- group_by(plotDat, CpG.class) %>%
      do(broom::tidy(addmargins(table(.$Gene.Feature)))) %>%
      tidyr::spread(Var1, n)

  }
  return(list(plot = g, table = tab))
}



# stackBar plot for REMP results
# (require to have column 'Probe','P','FDR','BetaFC' as well as the REMP annotations)
REMP_stackBar <- function(res,
                          p.fld = 'FDR',
                          p.thresh = 0.1){
  library(ggplot2)
  library(dplyr)
  res$P <- res[[p.fld]]

  ## significant hits:
  cg.sig <- filter(res, P < p.thresh) %>% with(Probe)
  cg.sig.hypo <- filter(res, P < p.thresh, BetaFC < 0) %>% with(Probe)
  cg.sig.hyper <- filter(res, P < p.thresh, BetaFC > 0) %>% with(Probe)

  cg.list <- list('Total' = res$Probe,
                  'differential methylation' = cg.sig,
                  'hypomethylation' = cg.sig.hypo,
                  'hypermethylation' = cg.sig.hyper)

  summDat <- plyr::ldply(cg.list, length) %>%
    dplyr::rename(CpG.class = .id, n = V1) %>%
    mutate(n = sprintf('n=%d', n))

  res.new = melt(res, id.vars=colnames(res)[1:14]) %>% filter(value != '')
  plotDat <- lapply(cg.list, function(x)
    filter(res.new, Probe %in% x) %>% with(variable))
  plotDat <- plyr::ldply(plotDat, function(x) data.frame(Gene.Feature = x)) %>%
    dplyr::rename(CpG.class = .id)
  plotDat$Gene.Feature = gsub('.symbol','',plotDat$Gene.Feature)
  plotDat <- plotDat %>% mutate(Gene.Feature = factor(Gene.Feature,
                                                      levels = c('InNM',"InNR", "InTSS",
                                                                 "In5UTR", 'In3UTR', 'InCDS', 'InExon')))
  plotDat$CpG.class = factor(plotDat$CpG.class, levels = c('Total',
                                                           'differential methylation',
                                                           'hypomethylation',
                                                           'hypermethylation'))
  # plot object:
  g <- ggplot(plotDat, aes(x = CpG.class)) +
    geom_bar(aes(fill = Gene.Feature), color = 'black', position = 'fill') +
    scale_fill_brewer(name = 'Gene Feature',type = 'seq',palette = 3)+
    scale_y_continuous(labels = scales::percent)+
    theme_bw() + ylab('% probes') +  xlab('CpG class') +
    geom_text(aes(y = 1.05, label = n), summDat, vjust = 0) +
    theme(axis.text.x = element_text())+
    scale_x_discrete(labels=c('Total' = 'Total', 'differential methylation' = 'differential\nmethylation',
                              'hypomethylation'='hypomethylation',
                              'hypermethylation'='hypermethylation'))
  #guides(fill = guide_legend(reverse = T))

  # table object:
  tab <- group_by(plotDat, CpG.class) %>%
    do(as.data.frame(addmargins(table(.$Gene.Feature)))) %>%
    tidyr::spread(Var1, Freq)
  return(list(plot = g, table = tab))
}



## Function for romer-style GSEA for 450K data (mvals). Key assumptions:
## 1) Restricts 450K data to CG sites with annotated genes
## 2) If CpG site annotated to multiple genes, takes the first one
## 3) Since genes have multiple CG sites annotated to them, they are represented that many times
## 4) Returns list of data.frame results
romer450k <- function(mval.mat, f, design, ctrst = ncol(design), geneSets = 1:7){
  library(limma)
  library(biomaRt)
  library(plyr)
  library(dplyr)

  stopifnot(all(c('cpgSite', 'GeneSymbol') %in% names(f)),
            all(rownames(mval.mat) %in% f$cpgSite),
            nrow(design) == ncol(mval.mat),
            all(geneSets %in% 1:7))

  cat('\n subsetting data to cpg sites with gene annotations... \n')
  fsub <- filter(f, GeneSymbol != '') %>%
    mutate(firstGene = laply(strsplit(GeneSymbol, ';'), '[', 1))
  mv <- mval.mat[intersect(fsub$cpgSite, rownames(mval.mat)), ]
  geneSymbols <- fsub$firstGene[match(rownames(mv), fsub$cpgSite)]
  cat(sprintf('%d of the original %d CpG sites map to %d unique genes \n',
              length(geneSymbols), nrow(mval.mat), n_distinct(geneSymbols)))

  cat('retrieving entrez gene IDs from ensembl... \n')
  ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  gn.1perDB <- getBM(c('entrezgene', 'hgnc_symbol'), filters = 'hgnc_symbol',
                     values = geneSymbols, mart = ensembl)


  entrezGeneIDs <- gn.1perDB$entrezgene[match(geneSymbols, gn.1perDB$hgnc_symbol)]
  cat(sprintf('%d (%d unique) of the %d gene symbols did not map \n',
              sum(is.na(entrezGeneIDs)), n_distinct(geneSymbols[is.na(entrezGeneIDs)]),
              length(entrezGeneIDs)))

  mv <- mv[!is.na(entrezGeneIDs),]
  ezIDs <- entrezGeneIDs[!is.na(entrezGeneIDs)]

  resList <- vector('list', length = length(geneSets))
  names(resList) <- paste0('C', geneSets)

  for(i in geneSets){
    cat(sprintf('Running romer for gene set C%d ... \n', i))
    gs <- get(load(sprintf(
      '/data/compute01/parelab/MSigDB_geneSets_human/human_c%d_v4.rdata', i)))
    ix <- ids2indices(gs, ezIDs)
    cat(sprintf('%d of the %d C%d gene sets exist in the dataset \n',
                length(ix), length(gs), i))
    r <- romer(mv, ix, design, contrast = ctrst, nrot = 499)
    r <- r %>% as.data.frame %>%
      mutate(GeneSet = rownames(r),
             Up.BH = p.adjust(Up, method = 'BH'),
             Down.BH = p.adjust(Down, method = 'BH'),
             Mixed.BH = p.adjust(Mixed, method = 'BH'))
    resList[[paste0('C', i)]] <- r
  }

  resList
}

## Function for roast-style GSEA for 450K data (mvals). Key assumptions:
## 1) Restricts 450K data to CG sites with annotated genes
## 2) If CpG site annotated to multiple genes, takes the first one
## 3) Since genes have multiple CG sites annotated to them, they are represented that many times
## 4) Returns list of data.frame results
roast450k <- function(mval.mat, f, design, ctrst = ncol(design), geneSets = 1:7){
  library(limma)
  library(biomaRt)
  library(plyr)
  library(dplyr)

  stopifnot(all(c('cpgSite', 'GeneSymbol') %in% names(f)),
            all(rownames(mval.mat) %in% f$cpgSite),
            nrow(design) == ncol(mval.mat),
            all(geneSets %in% 1:7))

  cat('\n subsetting data to cpg sites with gene annotations... \n')
  fsub <- filter(f, GeneSymbol != '') %>%
    mutate(firstGene = laply(strsplit(GeneSymbol, ';'), '[', 1))
  mv <- mval.mat[intersect(fsub$cpgSite, rownames(mval.mat)), ]
  geneSymbols <- fsub$firstGene[match(rownames(mv), fsub$cpgSite)]
  cat(sprintf('%d of the original %d CpG sites map to %d unique genes \n',
              length(geneSymbols), nrow(mval.mat), n_distinct(geneSymbols)))

  cat('retrieving entrez gene IDs from ensembl... \n')
  ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  gn.1perDB <- getBM(c('entrezgene', 'hgnc_symbol'), filters = 'hgnc_symbol',
                     values = geneSymbols, mart = ensembl)


  entrezGeneIDs <- gn.1perDB$entrezgene[match(geneSymbols, gn.1perDB$hgnc_symbol)]
  cat(sprintf('%d (%d unique) of the %d gene symbols did not map \n',
              sum(is.na(entrezGeneIDs)), n_distinct(geneSymbols[is.na(entrezGeneIDs)]),
              length(entrezGeneIDs)))

  mv <- mv[!is.na(entrezGeneIDs),]
  ezIDs <- entrezGeneIDs[!is.na(entrezGeneIDs)]

  resList <- vector('list', length = length(geneSets))
  names(resList) <- paste0('C', geneSets)

  for(i in geneSets){
    cat(sprintf('Running mroast for gene set C%d ... \n', i))
    gs <- get(load(sprintf(
      '/data/compute01/parelab/MSigDB_geneSets_human/human_c%d_v4.rdata', i)))
    ix <- ids2indices(gs, ezIDs)
    cat(sprintf('%d of the %d C%d gene sets exist in the dataset \n',
                length(ix), length(gs), i))
    r <- mroast(mv, ix, design, contrast = ctrst, nrot = 499)
    r <- mutate(r, GeneSet = rownames(r)) %>% dplyr::select(GeneSet, NGenes:FDR.Mixed)

    resList[[paste0('C', i)]] <- r
  }

  resList
}

# modification of DMRcate cpg.annotate() function to take my own test statistics as input
# the test statistics should look like a typical limma topTable
# rownames of statistics should be cpgnames, essential columns are: t, P.Value and adj.P.Val. betaFC can be missing but will become 0 if not provided.
# the CpGs in the object need to be same as in the output

# This function is useful if your differential methylation analysis model is not a typical limma model (eg. LME model etc.)
# Specify p if you want to use nominal p as differential cutoff, otherwise, specify fdr.
# eg. annot <- cpg.annotate.mod(betaMatrix, 'Beta', toptable, 'EPIC', fdr = 0.05)

cpg.annotate.mod = function (object, what = c("Beta", "M"), output, arraytype = c("EPIC", "450K"), fdr = NULL, p = NULL) {
  library(DMRcate)
  what <- match.arg(what)
  arraytype <- match.arg(arraytype)
  stopifnot(class(object) %in% c("matrix", "GenomicRatioSet"))
  if (arraytype == "450K") {
    grset <- makeGenomicRatioSetFromMatrix(object,
                                           array = "IlluminaHumanMethylation450k",
                                           annotation = "ilmn12.hg19",
                                           mergeManifest = TRUE, what = what)
  }
  if (arraytype == "EPIC") {
    grset <- makeGenomicRatioSetFromMatrix(object,
                                           array = "IlluminaHumanMethylationEPIC",
                                           annotation = "ilm10b2.hg19",
                                           mergeManifest = TRUE, what = what)
  }
  m <- match(rownames(grset), rownames(output))
  output <- output[m, ]
  anno <- getAnnotation(grset)
  stopifnot(all(rownames(anno) == rownames(output)))
  stat <- output$t
  if(!is.null(fdr)){
    annotated <- data.frame(ID = rownames(output), stat = stat,
                            CHR = anno$chr, pos = anno$pos, betafc = if(!is.null(output$betaFC)) output$betaFC else 0,
                            indfdr = output$adj.P.Val, is.sig = output$adj.P.Val < fdr)
  }
  if(!is.null(p)){
    annotated <- data.frame(ID = rownames(output), stat = stat,
                            CHR = anno$chr, pos = anno$pos, betafc = if(!is.null(output$betaFC)) output$betaFC else 0,
                            indfdr = output$adj.P.Val, is.sig = output$P.Value < p)
  }


  annotated <- annotated[order(annotated$CHR, annotated$pos), ]
  class(annotated) <- "annot"
  return(annotated)
}

# modification of DMRcate extractRanges() function to return extra information such as gene symbols and cpg IDs
# input should be a 'dmrcoutput' object
# eg. DMRranges <- extractRanges.mod(dmrcoutput, 'hg19')

extractRanges.mod = function (dmrcoutput, genome = c("hg19", "hg38", "mm10"))
{
  library(DMRcate)
  library(dplyr)
  env <- new.env(parent = emptyenv())
  data(dmrcatedata, envir = env)
  genome <- match.arg(genome)
  stopifnot(is(dmrcoutput, "dmrcate.output"))
  coords <- extractCoords(dmrcoutput$results$coord)
  coords <- cbind(coords, dmrcoutput$results[, c("no.cpgs",
                                                 "minfdr", "Stouffer", "maxbetafc", "meanbetafc")])
  coords$chromStart <- as.integer(as.character(coords$chromStart))
  coords$chromEnd <- as.integer(as.character(coords$chromEnd))
  ranges <- makeGRangesFromDataFrame(coords, keep.extra.columns = TRUE)
  switch(genome, hg19 = {
    tx = env$tx.hg19
  }, hg38 = {
    tx = env$tx.hg38
  }, mm10 = {
    tx = env$tx.mm10
  })
  promsidx <- as.data.frame(findOverlaps(ranges, promoters(tx,
                                                           2000, 2000)))
  proms <- tapply(promsidx$subjectHits, promsidx$queryHits,
                  function(x) tx[x])
  op.A <- sapply(proms, function(l) paste(l$tx_name, collapse = ", "))
  name.A <- names(proms)
  m.A <- as.numeric(name.A)
  M <- length(ranges)
  overlapping.promoters <- rep(NA_character_, M)
  overlapping.promoters[m.A] <- op.A
  ranges$overlapping.promoters <- overlapping.promoters

  op.B <- sapply(proms, function(l) paste(unique(l$gene_name), collapse = ", "))
  gene.symbols <- rep(NA_character_, M)
  gene.symbols[m.A] <- op.B
  ranges$Gene.Symbols <- gene.symbols

  dmrcoutput$input$CHR = as.character(dmrcoutput$input$CHR)
  dmrcoutput$input$ID = as.character(dmrcoutput$input$ID)

  update_cpgid <- function(region){
    chr <- region$seqnames
    d <- dmrcoutput$input %>% filter(CHR == chr, pos <= region$end, pos >= region$start)
    paste(d$ID, collapse = ' / ')
  }

  ranges = as.data.frame(ranges)
  ranges$seqnames = as.character(ranges$seqnames)
  ranges.update <- rowwise(ranges) %>%
    do(data.frame(., CpG.ID = update_cpgid(.))) %>%
    ungroup
  as.data.frame(ranges.update)
}

# modification of DMRcate extractRanges() function to return extra information such as gene symbols
# start from a result.range table
extractRanges.mod2 = function(dmrcoutput, genome = c("hg19", "hg38", "mm10"))
{
  env <- new.env(parent = emptyenv())
  data(dmrcatedata, envir = env)
  genome <- match.arg(genome)
  coords <- dmrcoutput[,c('seqnames', 'start', 'end')]
  colnames(coords) = c('chrom', 'chromStart', 'chromEnd')
  coords <- cbind(coords, dmrcoutput[, c("no.cpgs",
                                         "minfdr", "Stouffer", "maxbetafc", "meanbetafc")])
  coords$chromStart <- as.integer(as.character(coords$chromStart))
  coords$chromEnd <- as.integer(as.character(coords$chromEnd))
  ranges <- makeGRangesFromDataFrame(coords, keep.extra.columns = TRUE)
  switch(genome, hg19 = {
    tx = env$tx.hg19
  }, hg38 = {
    tx = env$tx.hg38
  }, mm10 = {
    tx = env$tx.mm10
  })
  promsidx <- as.data.frame(findOverlaps(ranges, promoters(tx,
                                                           2000, 2000)))
  proms <- tapply(promsidx$subjectHits, promsidx$queryHits,
                  function(x) tx[x])
  op.A <- sapply(proms, function(l) paste(l$tx_name, collapse = ", "))
  name.A <- names(proms)
  m.A <- as.numeric(name.A)
  M <- length(ranges)
  overlapping.promoters <- rep(NA_character_, M)
  overlapping.promoters[m.A] <- op.A
  ranges$overlapping.promoters <- overlapping.promoters

  op.B <- sapply(proms, function(l) paste(unique(l$gene_name), collapse = ", "))
  gene.symbols <- rep(NA_character_, M)
  gene.symbols[m.A] <- op.B
  ranges$gene.symbols <- gene.symbols

  dmrcoutput$input$CHR = as.character(dmrcoutput$input$CHR)
  dmrcoutput$input$ID = as.character(dmrcoutput$input$ID)

  update_cpgid <- function(region){
    chr <- region$seqnames
    d <- dmrcoutput$input %>% filter(CHR == chr, pos <= region$end, pos >= region$start)
    paste(d$ID, collapse = ' / ')
  }

  ranges = as.data.frame(ranges)
  ranges$seqnames = as.character(ranges$seqnames)
  ranges.update <- rowwise(ranges) %>%
    do(data.frame(., CpG.ID = update_cpgid(.))) %>%
    ungroup
  as.data.frame(ranges.update)
}

##########
# Gviz region plot for DMR (1 data track for 2 groups)

# prepare annotation for cpg islands as follows before using the function (hg19 data downloaded from UCSC)
#gunzip('data/cpgIslandExt.txt.gz')
#cpgIsland = fread('data/cpgIslandExt.txt')
#cpgIsland = cpgIsland %>% dplyr::select(V2,V3,V4,V7)
#colnames(cpgIsland) = c('chrom', 'chromStart', 'chromEnd', 'id')
#cpgIsland$id = paste0('Number of CpG = ',cpgIsland$id)
#cpgIsland = makeGRangesFromDataFrame(cpgIsland, keep.extra.columns = T)

plotRegion_2n = function(chr, from, to, grouplab, beta_g1, fDat, left=0, right=0, asize=1, grsize=1){
  library(Gviz)
  library(GenomicRanges)
  stopifnot('cpgSite' %in% colnames(fDat))
  # make annotation track
  atrack <- AnnotationTrack(cpgIsland, name = "CpG Island",
                            chromosome = chr, start = from, end = to,
                            #showFeatureId=T,
                            just.group = 'above',
                            size = asize,
                            background.panel = "#deebf7",
                            background.title = "#756bb1", fill='#bf812d')


  # make gene model track
  library(biomaRt)
  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  biomTrack = BiomartGeneRegionTrack(genome="hg19", biomart=mart,
                                     chromosome=chr, start=from, end=to, just.group = 'above',
                                     showId=T, geneSymbols=T, name = 'Gene Annotation',
                                     showFeatureId = T,size=grsize,
                                     col.line=NULL, transcriptAnnotation = "symbol",
                                     # filters=list(biotype="protein_coding"),
                                     collapseTranscripts = F,
                                     background.panel = "#FFFEDB",
                                     background.title = "#756bb1")
  biomTrack@range$symbol = paste0(biomTrack@range$symbol,' (',biomTrack@range$transcript,')')
  # make genome track and genome axis track
  gtrack <- GenomeAxisTrack() # genome track
  itrack <- IdeogramTrack(genome = 'hg19', chromosome = chr) # axis track

  # make my data track
  beta_g1fDat = fDat[match(rownames(beta_g1), fDat$cpgSite),]
  coordsg1 = as.data.frame(beta_g1fDat$chr)
  coordsg1 = coordsg1 %>% mutate(start = beta_g1fDat$pos, end = beta_g1fDat$pos)
  colnames(coordsg1)[1] = 'chr'
  betag1_coords = cbind(coordsg1,beta_g1)

  betag1_coords = makeGRangesFromDataFrame(na.omit(betag1_coords), keep.extra.columns = TRUE)
  dTrack1 <- DataTrack(betag1_coords, name = "Beta Value",chromosome = chr, genome = 'hg19',
                       background.title = "#756bb1",background.panel = "white",
                       groups = grouplab,type = c("a", "p","g"), legend=TRUE,
                       col = c('#2166ac','#b2182b')) # data track

  # plot all tracks
  plotTracks(list(itrack,gtrack,atrack,biomTrack,dTrack1),from = from, to = to,showId=TRUE, extend.left = left, extend.right = right)
}


# Gviz region plot for DMR (1 data track for paired data, 1 data track showing beta values and 1 data track showing beta differences)
plotRegion_pair = function(chr, from, to, group1, beta_g1, beta_g2, fDat, left=0, right=0, asize=1, grsize=1){
  library(Gviz)
  library(GenomicRanges)

  # make annotation track
  atrack <- AnnotationTrack(cpgIsland, name = "CpG Island",
                            chromosome = chr, start = from, end = to,
                            #showFeatureId=T,
                            just.group = 'above',
                            size = asize,
                            background.panel = "#deebf7",
                            background.title = "#756bb1", fill='#bf812d')


  # make gene model track
  #hg19db <-makeTxDbFromUCSC(genome = "hg19",
  #                          tablename = "knownGene") # can make txdb from databases,knownGene = UCSC,refGene=refseq
  #txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  #grtrack <- GeneRegionTrack(
  #  txdb,
  #  chromosome = chr, start = from, end = to,
  #  showId = TRUE,
  #  name = "Gene Annotation",
  #  showFeatureId=TRUE,showId = TRUE,geneSymbol=TRUE,
  #  transcriptAnnotation = "symbol",
  #  background.panel = "#FFFEDB",
  #  background.title = "#756bb1"
  #) # gene annotation track
  # alternatively use a biomart gene model track
  library(biomaRt)
  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  biomTrack = BiomartGeneRegionTrack(genome="hg19", biomart=mart,
                                     chromosome=chr, start=from, end=to, just.group = 'above',
                                     showId=T, geneSymbols=T, name = 'Gene Annotation',
                                     showFeatureId = T,size=grsize,
                                     col.line=NULL, transcriptAnnotation = "symbol",
                                     # filters=list(biotype="protein_coding"),
                                     collapseTranscripts = F,
                                     background.panel = "#FFFEDB",
                                     background.title = "#756bb1")
  biomTrack@range$symbol = paste0(biomTrack@range$symbol,' (',biomTrack@range$transcript,')')
  # make genome track and genome axis track
  gtrack <- GenomeAxisTrack() # genome track
  itrack <- IdeogramTrack(genome = 'hg19', chromosome = chr) # axis track

  # make my data track - BrF
  beta_g1fDat = fDat[match(rownames(beta_g1), rownames(fDat)),]
  coordsg1 = as.data.frame(beta_g1fDat$chr)
  coordsg1 = coordsg1 %>% mutate(start = beta_g1fDat$pos, end = beta_g1fDat$pos)
  colnames(coordsg1)[1] = 'chr'
  betag1_coords = cbind(coordsg1,beta_g1)

  betag1_coords = makeGRangesFromDataFrame(betag1_coords, keep.extra.columns = TRUE)
  dTrack1 <- DataTrack(betag1_coords, name = "Beta Value",chromosome = chr, genome = 'hg19',
                       background.title = "#756bb1",background.panel = "white",
                       groups = group1,type = c("a", "p","g"), legend=TRUE,
                       col = c('#b2182b','#2166ac')) # data track

  # make my data track - LgF
  beta_g2fDat = fDat[match(rownames(beta_g2), rownames(fDat)),]
  coordsg2 = as.data.frame(beta_g2fDat$chr)
  coordsg2 = coordsg2 %>% mutate(start = beta_g2fDat$pos, end = beta_g2fDat$pos)
  colnames(coordsg2)[1] = 'chr'
  betag2_coords = cbind(coordsg2,beta_g2)

  betag2_coords = makeGRangesFromDataFrame(betag2_coords, keep.extra.columns = TRUE)
  dTrack2 <- DataTrack(betag2_coords, name = "Difference in\nBeta Value\n(Parenchymal - Airway)",chromosome = chr, genome = 'hg19',
                       background.title = "#756bb1",background.panel = "white",
                       type = c("a", "p","g"), legend=TRUE, col = '#4d9221') # data track
  # plot all tracks
  plotTracks(list(itrack,gtrack,atrack,biomTrack,dTrack1, dTrack2),from = from, to = to,showId=TRUE, extend.left = left, extend.right = right)
}


##################
# Gviz region plot for DMR (2 data tracks for paired 2 groups)
plotRegion_double = function(chr, from, to, group1, group2, beta_g1, beta_g2, fDat, left=0, right=0, asize=1, grsize=1){
  library(Gviz)
  library(GenomicRanges)

  # make annotation track
  atrack <- AnnotationTrack(cpgIsland, name = "CpG Island",
                            chromosome = chr, start = from, end = to,
                            #showFeatureId=T,
                            just.group = 'above',
                            size = asize,
                            background.panel = "#deebf7",
                            background.title = "#756bb1", fill='#bf812d')

  # use a biomart gene model track
  library(biomaRt)
  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  biomTrack = BiomartGeneRegionTrack(genome="hg19", biomart=mart,
                                     chromosome=chr, start=from, end=to, just.group = 'above',
                                     showId=T, geneSymbols=T, name = 'Gene Annotation',
                                     showFeatureId = T,size=grsize,
                                     col.line=NULL, transcriptAnnotation = "symbol",
                                     # filters=list(biotype="protein_coding"),
                                     collapseTranscripts = F,
                                     background.panel = "#FFFEDB",
                                     background.title = "#756bb1")
  biomTrack@range$symbol = paste0(biomTrack@range$symbol,' (',biomTrack@range$transcript,')')

  # make genome track and genome axis track
  gtrack <- GenomeAxisTrack() # genome track
  itrack <- IdeogramTrack(genome = 'hg19', chromosome = chr) # axis track

  # make my data track - BrF
  beta_g1fDat = fDat[match(rownames(beta_g1), rownames(fDat)),]
  coordsg1 = as.data.frame(beta_g1fDat$chr)
  coordsg1 = coordsg1 %>% mutate(start = beta_g1fDat$pos, end = beta_g1fDat$pos)
  colnames(coordsg1)[1] = 'chr'
  betag1_coords = cbind(coordsg1,beta_g1)

  betag1_coords = makeGRangesFromDataFrame(betag1_coords, keep.extra.columns = TRUE)
  dTrack1 <- DataTrack(betag1_coords, name = "Beta Value\n(Airway)",chromosome = chr, genome = 'hg19',
                       background.title = "#756bb1",background.panel = "white",
                       groups = group1,type = c("a", "p","g"), legend=TRUE,
                       col = c('#b2182b','#2166ac')) # data track

  # make my data track - LgF
  beta_g2fDat = fDat[match(rownames(beta_g2), rownames(fDat)),]
  coordsg2 = as.data.frame(beta_g2fDat$chr)
  coordsg2 = coordsg2 %>% mutate(start = beta_g2fDat$pos, end = beta_g2fDat$pos)
  colnames(coordsg2)[1] = 'chr'
  betag2_coords = cbind(coordsg2,beta_g2)

  betag2_coords = makeGRangesFromDataFrame(betag2_coords, keep.extra.columns = TRUE)
  dTrack2 <- DataTrack(betag2_coords, name = "Beta Value\n(Parenchymal)",chromosome = chr, genome = 'hg19',
                       background.title = "#756bb1",background.panel = "white",
                       groups = group2,type = c("a", "p","g"), legend=TRUE,
                       col = c('#b2182b','#2166ac')) # data track

  # plot all tracks
  plotTracks(list(itrack,gtrack,atrack, biomTrack,dTrack1, dTrack2),from = from, to = to,showId=TRUE, extend.left = left, extend.right = right)
}

#cpgDensity = anno %>% mutate(chrom=chr, chromStart=pos, chromEnd=pos) %>% select(chrom, chromStart, chromEnd, Relation_to_Island) %>% makeGRangesFromDataFrame(., keep.extra.columns = T)
#anno$GeneFeature = sapply(strsplit(anno$UCSC_RefGene_Group,split=';'), function(x) unique(x)[1])
#geneFeature = anno %>% mutate(chrom=chr, chromStart=pos, chromEnd=pos) %>% select(chrom, chromStart, chromEnd, GeneFeature) %>% makeGRangesFromDataFrame(., keep.extra.columns = T)

plotRegion = function(chr, from, to, group1, beta_g1, beta_g2, fDat, left=0, right=0, asize=1){
  library(Gviz)
  library(GenomicRanges)
  # make annotation track
  chrom=paste0('chr',chr)
  atrack <- AnnotationTrack(cpgDensity[cpgDensity@seqnames==chrom,],
                            name = "CpG",
                            chromosome = chr, start = from, end = to,
                            group = cpgDensity[cpgDensity@seqnames==chrom,]$Relation_to_Island,
                            just.group = 'above',
                            size = asize,  background.title = "#756bb1",
                            genome = "hg19",min.width = 2, min.distance = 5,collapse = TRUE,fontsize=15)
  feature(atrack)= cpgDensity[cpgDensity@seqnames==chrom,]$Relation_to_Island

  atrack2 <- AnnotationTrack(geneFeature[geneFeature@seqnames==chrom,],
                             name = "Gene",
                             chromosome = chr, start = from, end = to,
                             group = geneFeature[geneFeature@seqnames==chrom,]$GeneFeature,
                             just.group = 'above',
                             size = asize,  fill='#bf812d',background.title = "#756bb1",
                             genome = "hg19",min.width = 2, min.distance = 5,collapse = TRUE,fontsize=15)
  feature(atrack2)= geneFeature[geneFeature@seqnames==chrom,]$GeneFeature

  # make genome track and genome axis track
  gtrack <- GenomeAxisTrack() # genome track
  #itrack <- IdeogramTrack(genome = 'hg19', chromosome = chr) # axis track

  # make my data track - BrF
  beta_g1fDat = fDat[match(rownames(beta_g1), rownames(fDat)),]
  coordsg1 = as.data.frame(beta_g1fDat$chr)
  coordsg1 = coordsg1 %>% mutate(start = beta_g1fDat$pos, end = beta_g1fDat$pos)
  colnames(coordsg1)[1] = 'chr'
  betag1_coords = cbind(coordsg1,beta_g1)

  betag1_coords = makeGRangesFromDataFrame(betag1_coords, keep.extra.columns = TRUE)
  dTrack1 <- DataTrack(betag1_coords, name = "Beta\n ",chromosome = chr, genome = 'hg19',
                       background.title = "#756bb1",background.panel = "white",
                       groups = group1,type = c("a", "p","g"), legend=TRUE,
                       col = c('#b2182b','#2166ac'),fontsize=17) # data track
  # make my data track - LgF
  beta_g2fDat = fDat[match(rownames(beta_g2), rownames(fDat)),]
  coordsg2 = as.data.frame(beta_g2fDat$chr)
  coordsg2 = coordsg2 %>% mutate(start = beta_g2fDat$pos, end = beta_g2fDat$pos)
  colnames(coordsg2)[1] = 'chr'
  betag2_coords = cbind(coordsg2,beta_g2)

  betag2_coords = makeGRangesFromDataFrame(betag2_coords, keep.extra.columns = TRUE)
  dTrack2 <- DataTrack(betag2_coords, name = "Beta Difference\n ",chromosome = chr, genome = 'hg19', background.title = "#756bb1",background.panel = "white",type = c("a", "p","g"), legend=TRUE, col = '#4d9221',fontsize=17) # data track
  # plot all tracks
  plotTracks(list(gtrack,atrack,atrack2,dTrack1,dTrack2),
             from = from, to = to,showId=T,
             extend.left = left, extend.right = right,stacking = "dense",
             groupAnnotation = "feature",
             OpenSea='#313695',N_Shelf='#74add1',N_Shore='#e0f3f8',
             Island='#8c510a',S_Shelf='#80cdc1', S_Shore='#01665e',
             Body='#e31a1c',TSS1500='#ff7f00',TSS200='#fb9a99',
             `1stExon`='#6a3d9a',`5'UTR`='#33a02c',`3'UTR`='#1f78b4')
}

# plot DMR region for a 2 by 2 factorial design
# anno$GeneFeature = sapply(strsplit(anno$UCSC_RefGene_Group,split=';'), function(x) unique(x)[1])
#anno$GeneFeature = factor(anno$GeneFeature,levels=c("3'UTR", "5'UTR", "1stExon", "TSS200", "TSS1500","Body"))
#anno$Relation_to_Island = factor(anno$Relation_to_Island, levels=c("N_Shelf", "N_Shore", "S_Shelf", "S_Shore","OpenSea","Island"))
# anno is the annotation data.frame which has columns: Name, chr, pos, Relation_to_Island, GeneFeature (as above)
# group1 is a single group label
# group2 is a single group label
# pheno_in_g1 is a vector of phenotypes in group1
# pheno_in_g2 is a vector of phenotypes in group2
# beta1 is a beta matrix of group1 subjects (cpg by subjects)
# beta2 is a beta matrix of group2 subjects (cpg by subjects)
plotDMR_2by2 = function(chr, from, to, anno, beta1, beta2, pheno_in_g1, pheno_in_g2, group1, group2, xlabel, ylabel){
  myColors <- c('#8dd3c7','#ffffb3','#fb8072','#bebada','#80b1d3','#fdb462',
                '#dfc27d', '#bf812d', '#c7eae5', '#35978f','#4575b4', '#f03b20')
  names(myColors) <- c(levels(anno$GeneFeature), levels(anno$Relation_to_Island))
  mycolScale <- scale_fill_manual(name = "Features",values = myColors)

  f = anno[anno$chr==paste0('chr',chr) & anno$pos >= from & anno$pos <= to,]
  plotDat1 = cbind(as.data.frame(t(beta1[f$Name,])), Pheno = pheno_in_g1, Group=group1)
  plotDat2 = cbind(as.data.frame(t(beta2[f$Name,])), Pheno = pheno_in_g2, Group=group2)
  plotDat = melt(rbind(plotDat1,plotDat2), id.vars=c('Group','Pheno')) %>% mutate(Position = as.numeric(as.factor(variable)))
  plotDat = merge(plotDat, f, by.x='variable',by.y='Name')
  plotDat = plotDat %>% arrange(pos) %>% mutate(Position = as.numeric(as.factor(pos)))

  p=ggplot(plotDat,aes(x=Position, y=value))+
    geom_point(aes(color=Pheno))+
    stat_smooth(aes(color=Pheno),se=F)+
    geom_rect(mapping=aes(xmin=Position-0.5, xmax=Position+0.5, ymin=-0.04, ymax=0, fill=Relation_to_Island),color = NA)+
    geom_rect(mapping=aes(xmin=Position-0.5, xmax=Position+0.5, ymin=-0.08, ymax=-0.04, fill=GeneFeature),color=NA)+
    facet_grid(Group ~.)+
    theme_bw()+
    scale_color_brewer('Groups',type = 'qual',palette = 6, direction = -1)+
    xlab(xlabel)+
    ylab(ylabel)+
    mycolScale+
    scale_x_continuous(breaks= function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
  return(p)
}

# plot DMR region for a 1-factor design
plotDMR = function(chr, from, to, anno, beta1, pheno_in_g1, group1, xlabel, ylabel){
  myColors <- c('#8dd3c7','#ffffb3','#fb8072','#bebada','#80b1d3','#fdb462',
                '#dfc27d', '#bf812d', '#c7eae5', '#35978f','#4575b4', '#f03b20')
  names(myColors) <- c(levels(anno$GeneFeature), levels(anno$Relation_to_Island))
  mycolScale <- scale_fill_manual(name = "Features",values = myColors)

  f = anno[anno$chr==paste0('chr',chr) & anno$pos >= from & anno$pos <= to,]
  plotDat = cbind(as.data.frame(t(beta1[f$Name,])), Pheno = pheno_in_g1)
  plotDat = melt(plotDat, id.vars=c('Pheno')) %>% mutate(Position = as.numeric(as.factor(variable)))
  plotDat = merge(plotDat, f, by.x='variable',by.y='Name')
  plotDat = plotDat %>% arrange(pos) %>% mutate(Position = as.numeric(as.factor(pos)))

  p=ggplot(plotDat,aes(x=Position, y=value))+
    geom_rect(mapping=aes(xmin=Position-0.5, xmax=Position+0.5, ymin=-0.04, ymax=0, fill=Relation_to_Island),color = NA)+
    geom_rect(mapping=aes(xmin=Position-0.5, xmax=Position+0.5, ymin=-0.08, ymax=-0.04, fill=GeneFeature),color=NA)+
    geom_point(aes(color=Pheno))+
    stat_smooth(aes(color=Pheno),se=F)+
    theme_bw()+
    scale_color_brewer('Groups',type = 'qual',palette = 6, direction = -1)+
    xlab(xlabel)+
    ylab(ylabel)+
    mycolScale+
    scale_x_continuous(breaks= function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
  return(p)
}

# enrichment analysis function using EnrichR (no background correction)
# commonly used databases: GO_Molecular_Function_2017b, GO_Biological_Process_2017b, GO_Cellular_Component_2017b,
#                          KEGG_2016, Reactome_2016, WikiPathways_2016
# example: enrichr.analysis(DEGenes)
enrichr.analysis = function(genes, databases=c('KEGG_2016', 'GO_Biological_Process_2017b', 'GO_Cellular_Component_2017b', 'GO_Molecular_Function_2017b')){
  library(enrichR)
  enrichment = enrichr(genes, databases)
  p = list()

  for (i in 1:length(enrichment)){
    p.cutoff.enrich = enrichment[[3]] %>% filter(Adjusted.P.value <= 0.1) %>%
      with(max(P.value)) %>% signif(3)
    p[[i]] = ggplot(enrichment[[i]] %>% arrange(P.value) %>% .[1:10,], aes(x=reorder(Term, -log10(P.value)), y=-log10(P.value))) +
      geom_bar(stat = "identity")+
      geom_hline(yintercept = p.cutoff.enrich, linetype = 'dotted')+
      theme(axis.text.x = element_text(angle = 0, hjust = 1))+
      xlab(gsub('_',' ',names(enrichment[i])))+
      theme(legend.position = "none")+
      coord_flip()+
      theme_bw()
  }
  names(p) = names(enrichment)
  return(list(enrichment = enrichment, plots = p))
}

# enrichment analysis function using WebGestalt (with background correction)
# commonly used databases: geneontology_Biological_Process, geneontology_Cellular_Component,
# 							geneontology_Molecular_Function, geneontology_Biological_Process_noRedundant,
#							geneontology_Cellular_Component_noRedundant, geneontology_Molecular_Function_noRedundant,
#							pathway_KEGG, pathway_Reactome, pathway_Wikipathway, disease_OMIM
# same as 2017 version
# example: webGestalt.analysis(DEGenes, referenceGenes)
webGestalt.analysis = function(genes, ref, databases = c('pathway_KEGG', 'geneontology_Biological_Process_noRedundant', 'geneontology_Cellular_Component_noRedundant', 'geneontology_Molecular_Function_noRedundant', 'pathway_KEGG')){
  library(WebGestaltR)
  library(ggplot2)
  library(dplyr)

  enrichment = list()
  for(n in 1:length(databases)){
    enrichment[[n]] <- WebGestaltR(enrichMethod="ORA",organism="hsapiens",
                                   enrichDatabase=databases[n],
                                   interestGene=genes, interestGeneType="genesymbol",
                                   referenceGene=ref, referenceGeneType="genesymbol",
                                   sigMethod='top', minNum = 5, maxNum = 2000, referenceSet='genome',
                                   isOutput=F, outputDirectory=getwd())
  }

  names(enrichment) = databases %>% gsub('_',' ',.) %>% gsub('geneontology','GO',.)
  #enrichment$PValue[enrichment$PValue==0]=2.22e-16 # p-value too small to display, change to the smallest p-value
  #enrichment$FDR[enrichment$FDR==0]=2.22e-16
  p = list()
  labl=list(expression('FDR '<' 0.1'), expression('FDR '>=' 0.1'))
  for (i in 1:length(enrichment)){
    p[[i]] = ggplot(enrichment[[i]] %>% dplyr::mutate(color = ifelse(FDR<0.1, 'FDR < 0.1', 'FDR \u2265 0.1')), aes(x=reorder(description, -log10(pValue)), y=-log10(pValue), fill = color, size = overlap)) +
      geom_point(shape = 21)+
      scale_size_area(name = 'Number of Overlapping Genes', max_size = 7)+
      scale_fill_brewer(type = 'qual', palette = 2, name = 'FDR', limits = c('FDR < 0.1', 'FDR \u2265 0.1'), labels=labl)+
      xlab(names(enrichment[i]))+
      ylab(expression(paste('-log'[10],' P')))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 0, hjust = 1))+
      coord_flip()
  }
  names(p) = names(enrichment)
  return(list(enrichment = enrichment, plots = p))
}

#noquote(unlist(format(.Machine)))
#           double.eps        double.neg.eps           double.xmin           double.xmax           double.base
#         2.220446e-16          1.110223e-16         2.225074e-308         1.797693e+308                     2
#        double.digits       double.rounding          double.guard     double.ulp.digits double.neg.ulp.digits
#                   53                     5                     0                   -52                   -53
#      double.exponent        double.min.exp        double.max.exp           integer.max           sizeof.long
#                   11                 -1022                  1024            2147483647                     4
#      sizeof.longlong     sizeof.longdouble        sizeof.pointer
#                    8                    16                     8

webGestalt.analysis2 = function(genes, ref, databases = c('pathway_KEGG', 'geneontology_Biological_Process_noRedundant', 'geneontology_Cellular_Component_noRedundant', 'geneontology_Molecular_Function_noRedundant', 'pathway_KEGG')){
  library(WebGestaltR)
  library(ggplot2)
  library(dplyr)

  enrichment = list()
  for(n in 1:length(databases)){
    enrichment[[n]] <- WebGestaltR(enrichMethod="ORA",organism="hsapiens",
                                   enrichDatabase=databases[n],
                                   interestGene=genes, interestGeneType="genesymbol",
                                   referenceGene=ref, referenceGeneType="genesymbol",
                                   sigMethod='fdr', minNum = 5, maxNum = 2000, referenceSet='genome',
                                   isOutput=F, outputDirectory = getwd())
  }

  names(enrichment) = databases %>% gsub('_',' ',.) %>% gsub('geneontology','GO',.)
  #enrichment$PValue[enrichment$PValue==0]=2.22e-16 # p-value too small to display, change to the smallest p-value
  #enrichment$FDR[enrichment$FDR==0]=2.22e-16
  p = list()
  for (i in 1:length(enrichment)){
    p[[i]] = ggplot(enrichment[[i]] %>% dplyr::mutate(color = ifelse(FDR<0.05, 'FDR < 0.05', 'FDR >= 0.05')), aes(x=reorder(description, -log10(PValue)), y=-log10(PValue), fill = color, size = O)) +
      geom_point(shape = 21)+
      scale_size_area(name = '# of Overlapping Genes', max_size = 7)+
      scale_fill_brewer(type = 'qual', palette = 2, direction = -1, name = 'FDR', limits = c('FDR < 0.05', 'FDR >= 0.05'))+
      xlab(names(enrichment[i]))+
      ylab('-log10 p-value')+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 0, hjust = 1))+
      coord_flip()
  }
  names(p) = names(enrichment)
  return(list(enrichment = enrichment, plots = p))
}


################## function to do model selection of an LME model
###  and to plot a PCs vs variables heatmap
###  the color of each grid indicates the p-value of model lme(PC ~ variable)
### adapted from Amrit Singh's compVar() function to do LME
### assuming each subject have multiple samples,
### the samples should have colname 'SampleID' and the subjects should have colname 'SubjectID'
compVar_lme = function (demo, eset, variables, ncomp = 10)
{
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(nlme)
  #demo = demo %>% filter(complete.cases(.))
  pcaX <- prcomp(eset[intersect(demo$SampleID,rownames(eset)),], scale. = TRUE, center = TRUE)
  demo <- demo[match(rownames(pcaX$x),demo$SampleID),]
  pval <- do.call(rbind, lapply(variables, function(i) {
    apply(pcaX$x, 2, function(j) {
      pred <- demo[, c(i,'SubjectID')]
      predictor <- demo[, i]
      dat <- cbind('x'= j, pred)
      form <- as.formula(paste0("as.numeric(x) ~ ", i))
      if (class(predictor) == "factor") {
        if (nlevels(predictor) == 2) {
          coef(summary(lme(form, random = ~ 1 | SubjectID, data = dat)))[2,"p-value"]
        }
        else {
          anova(lme(form, random = ~ 1 | SubjectID, data = dat))[2,"p-value"]
        }
      }
      else {
        coef(summary(lme(form, random = ~ 1 | SubjectID, data = dat)))[2,"p-value"]
      }
    })
  }))
  rownames(pval) <- variables
  colnames(pval) <- paste(colnames(pval), paste0(round(100 * (pcaX$sdev^2/sum(pcaX$sdev^2)), 1), "%"), sep = "-")
  pval <- pval[, 1:ncomp]
  pvalheatmap <- pval
  pvalheatmap[pvalheatmap < 0.01] <- 0.01
  pvalheatmap[pvalheatmap > 0.1] <- 1
  pvalheatmap[pvalheatmap > 0.01 & pvalheatmap < 0.05] <- 0.05
  pvalheatmap[pvalheatmap > 0.05 & pvalheatmap < 0.1] <- 0.1
  pvalheatmap[pvalheatmap == "0.01"] <- "p < 0.01"
  pvalheatmap[pvalheatmap == "0.05"] <- "0.01 < p < 0.05"
  pvalheatmap[pvalheatmap == "0.1"] <- "0.05 < p < 0.10"
  pvalheatmap[pvalheatmap == "1"] <- "p > 0.10"
  p <- pvalheatmap %>% as.data.frame %>% mutate(Variable = rownames(.)) %>%
    gather(Threshold, Value, -Variable) %>%
    mutate(Threshold = factor(Threshold, levels = unique(Threshold))) %>%
    mutate(Variable = factor(as.character(Variable), levels = variables)) %>%
    mutate(Value = factor(Value,levels = c("p < 0.01", "0.01 < p < 0.05", "0.05 < p < 0.10", "p > 0.10"))) %>%
    ggplot(aes(Threshold, Variable)) +
    geom_tile(aes(fill = Value), colour = "white") +
    scale_fill_manual(values = rev(brewer.pal(n = 8, name = "Blues")[c(2, 4, 6, 8)])) +
    #customTheme(sizeStripFont = 10, xAngle = 40, hjust = 1, vjust = 1, xSize = 10, ySize = 10,
    #xAxisSize = 10, yAxisSize = 10) +
    xlab("") + ylab("")
  return(list(pval = pval, pvalheatmap = pvalheatmap, p = p))
}

################## function to do model selection of an LME model
###  and to plot a PCs vs variables heatmap
###  the color of each grid indicates the p-value of model lme(PC ~ variable)
### adapted from Amrit Singh's compVar() function to do LME
### need to have 'SampleID' as the colname for the samples and 'SubjectID' as the colname for the subjects

compVar_lme = function (demo, eset, variables, ncomp = 10)
{
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(nlme)
  demo = demo %>% filter(complete.cases(.))
  pcaX <- prcomp(eset[intersect(demo$SampleID,rownames(eset)),], scale. = TRUE, center = TRUE)
  demo <- demo[match(rownames(pcaX$x),demo$SampleID),]
  pval <- do.call(rbind, lapply(variables, function(i) {
    apply(pcaX$x, 2, function(j) {
      pred <- demo[, c(i,'SubjectID')]
      predictor <- demo[, i]
      dat <- cbind('x'= j, pred)
      form <- as.formula(paste0("as.numeric(x) ~ ", i))
      if (class(predictor) == "factor") {
        if (nlevels(predictor) == 2) {
          coef(summary(lme(form, random = ~ 1 | SubjectID, data = dat)))[2,"p-value"]
        }
        else {
          anova(lme(form, random = ~ 1 | SubjectID, data = dat))[2,"p-value"]
        }
      }
      else {
        coef(summary(lme(form, random = ~ 1 | SubjectID, data = dat)))[2,"p-value"]
      }
    })
  }))
  rownames(pval) <- variables
  colnames(pval) <- paste(colnames(pval), paste0(round(100 * (pcaX$sdev^2/sum(pcaX$sdev^2)), 1), "%"), sep = "-")
  pval <- pval[, 1:ncomp]
  pvalheatmap <- pval
  pvalheatmap[pvalheatmap < 0.01] <- 0.01
  pvalheatmap[pvalheatmap > 0.1] <- 1
  pvalheatmap[pvalheatmap > 0.01 & pvalheatmap < 0.05] <- 0.05
  pvalheatmap[pvalheatmap > 0.05 & pvalheatmap < 0.1] <- 0.1
  pvalheatmap[pvalheatmap == "0.01"] <- "p < 0.01"
  pvalheatmap[pvalheatmap == "0.05"] <- "0.01 < p < 0.05"
  pvalheatmap[pvalheatmap == "0.1"] <- "0.05 < p < 0.10"
  pvalheatmap[pvalheatmap == "1"] <- "p > 0.10"
  p <- pvalheatmap %>% as.data.frame %>% mutate(Variable = rownames(.)) %>%
    gather(Threshold, Value, -Variable) %>%
    mutate(Threshold = factor(Threshold, levels = unique(Threshold))) %>%
    mutate(Variable = factor(as.character(Variable), levels = variables)) %>%
    mutate(Value = factor(Value,levels = c("p < 0.01", "0.01 < p < 0.05", "0.05 < p < 0.10", "p > 0.10"))) %>%
    ggplot(aes(Threshold, Variable)) +
    geom_tile(aes(fill = Value), colour = "white") +
    scale_fill_manual(values = rev(brewer.pal(n = 8, name = "Blues")[c(2, 4, 6, 8)])) +
    #customTheme(sizeStripFont = 10, xAngle = 40, hjust = 1, vjust = 1, xSize = 10, ySize = 10,
    #xAxisSize = 10, yAxisSize = 10) +
    xlab("") + ylab("")
  return(list(pval = pval, pvalheatmap = pvalheatmap, p = p))
}


################## function to do model selection of a linear model
compVar = function (demo, eset, variables, ncomp = 10)
{
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  demo = demo %>% filter(complete.cases(.))
  pcaX <- prcomp(eset[intersect(demo$SampleID,rownames(eset)),], scale. = TRUE, center = TRUE)
  demo <- demo[match(rownames(pcaX$x),demo$SampleID),]
  pval <- do.call(rbind, lapply(variables, function(i) {
    apply(pcaX$x, 2, function(j) {
      predictor <- demo[, i]
      dat <- cbind('x'= j, predictor %>% as.data.frame())
      form <- as.formula("as.numeric(x) ~ predictor")
      if (class(predictor) == "factor") {
        if (nlevels(predictor) == 2) {
          coef(summary(lm(form, data = dat)))[2,"Pr(>|t|)"]
        }
        else {
          anova(lm(form, data = dat))[2,"Pr(>F)"]
        }
      }
      else {
        coef(summary(lm(form, data = dat)))[2,"Pr(>|t|)"]
      }
    })
  }))
  rownames(pval) <- variables
  colnames(pval) <- paste(colnames(pval), paste0(round(100 * (pcaX$sdev^2/sum(pcaX$sdev^2)), 1), "%"), sep = "-")
  pval <- pval[, 1:ncomp]
  pvalheatmap <- pval
  pvalheatmap[pvalheatmap < 0.01] <- 0.01
  pvalheatmap[pvalheatmap > 0.1] <- 1
  pvalheatmap[pvalheatmap > 0.01 & pvalheatmap < 0.05] <- 0.05
  pvalheatmap[pvalheatmap > 0.05 & pvalheatmap < 0.1] <- 0.1
  pvalheatmap[pvalheatmap == "0.01"] <- "p < 0.01"
  pvalheatmap[pvalheatmap == "0.05"] <- "0.01 < p < 0.05"
  pvalheatmap[pvalheatmap == "0.1"] <- "0.05 < p < 0.10"
  pvalheatmap[pvalheatmap == "1"] <- "p > 0.10"
  p <- pvalheatmap %>% as.data.frame %>% mutate(Variable = rownames(.)) %>%
    gather(Threshold, Value, -Variable) %>%
    mutate(Threshold = factor(Threshold, levels = unique(Threshold))) %>%
    mutate(Variable = factor(as.character(Variable), levels = variables)) %>%
    mutate(Value = factor(Value,levels = c("p < 0.01", "0.01 < p < 0.05", "0.05 < p < 0.10", "p > 0.10"))) %>%
    ggplot(aes(Threshold, Variable)) +
    geom_tile(aes(fill = Value), colour = "white") +
    scale_fill_manual(values = rev(brewer.pal(n = 8, name = "Blues")[c(2, 4, 6, 8)])) +
    #customTheme(sizeStripFont = 10, xAngle = 40, hjust = 1, vjust = 1, xSize = 10, ySize = 10,
    #xAxisSize = 10, yAxisSize = 10) +
    xlab("") + ylab("")
  return(list(pval = pval, pvalheatmap = pvalheatmap, p = p))
}




# modified ggbiplot to get rid of equal coords
ggbiplot = function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE,
                     obs.scale = 1 - scale, var.scale = scale, groups = NULL,
                     ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3,
                     alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69,
                     varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE,
                     ...){
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord),
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale,
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)",
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) +
    ylab(u.axis.labs[2])
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi,
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r *
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"),
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0,
                                           xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2,
                                                                                                  "picas")), color = muted("red"))
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups),
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2,
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    g <- g + geom_text(data = df.v, aes(label = varname,
                                        x = xvar, y = yvar, angle = angle, hjust = hjust),
                       color = "darkred", size = varname.size)
  }
  return(g)
}









## manhattan plot
## used when y-axis should be in -log10 scale
gg.manhattan = function(dataframe, title=NULL, cols = c('red', 'skyblue'),
                        chr.fld = 'CHR', pos.fld = 'BP', p.fld = 'P', snp.fld = 'rsid', gene.fld=NULL,
                        chrs = 1:22, ymax = NULL,
                        suggestiveline=0, genomewideline=-log10(5e-8),
                        annotate=F, annotate.fld='SNP', SNPlist=NULL,lab2lim=NULL,lab3lim=NULL) {
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)

  if (annotate & is.null(SNPlist))
    stop("You requested annotation but provided no SNPlist!")
  #if(annotate & any(! unlist(SNPlist) %in% dataframe[[snp.fld]]))
  #  stop('No SNPs in SNPlist are found in the dataframe')

  if(is.data.table(dataframe)) dataframe <- as.data.frame(dataframe)

  stopifnot(all(c(chr.fld, pos.fld, p.fld, snp.fld) %in% names(dataframe)),
            is.data.frame(dataframe),
            all(!is.na(dataframe[[chr.fld]])),
            all(unlist(lapply(dataframe[c(chr.fld, pos.fld, p.fld)], is.numeric))))
  if(is.null(gene.fld)){
    d <- dataframe[c(chr.fld, pos.fld, p.fld, snp.fld)]
    names(d) <- c('CHR', 'BP', 'P', 'SNP')
  }else{
    d <- dataframe[c(chr.fld, pos.fld, p.fld, snp.fld, gene.fld)]
    names(d) <- c('CHR', 'BP', 'P', 'SNP','Gene')
  }

  #limit to only chrs 1-22
  d <- d %>% filter(CHR %in% chrs, P > 0, P <= 1, !is.na(BP)) %>%
    mutate(logp = -log10(P)) %>% arrange(CHR, BP)

  # get position references for each CHR:
  d$CHR <- factor(d$CHR)
  posDat <- d %>% group_by(CHR) %>% dplyr::summarise(maxBP = max(BP))

  tmp <- data.frame(0, 0); names(tmp) <- names(posDat)
  posDat <- rbind(tmp, posDat)
  posDat$cumBP <- sapply(1:nrow(posDat), function(x) sum(posDat$maxBP[1:x]))
  posDat <- data.frame(CHR = posDat$CHR[-1], addBP = posDat$cumBP[-nrow(posDat)])

  d <- merge(d, posDat, all.x = T)

  d <- d %>% mutate(pos = BP + addBP)
  ticks <- d %>% group_by(CHR) %>% dplyr::summarise(x = median(pos))

  mycols <- rep(cols, length(chrs))[1:length(unique(chrs))]

  if(length(unique(chrs)) == 1){
    plot <- ggplot(d, aes(x = pos, y = logp, color = CHR)) +
      xlab(paste('Chromosome', unique(chrs)))
  } else {
    plot <- ggplot(d, aes(x = pos, y = logp, color = CHR)) +
      scale_x_continuous(name = 'Chromosome', breaks = ticks$x, labels = ticks$CHR)
  }
  upper <- ceiling(max(d$logp))
  if(!is.null(ymax) && is.numeric(ymax)){
    upper <- ceiling(ymax)
  }
  plot <- plot + geom_point() + scale_color_manual(values = mycols) +
    scale_y_continuous(breaks = seq(1, upper, by = floor(upper/6))) +
    ylab(expression(-log[10](italic(p)))) +
    ggtitle(title) +
    theme_bw() +
    theme(axis.text.y = element_text(size = rel(2), color = 'black'),
          axis.text.x = element_text(size = rel(1.5), color = 'black'),
          legend.position = 'none', plot.title = element_text(size = rel(5), hjust = 0.5),
          panel.grid.minor = element_blank())

  if(annotate){
    if(annotate.fld=='SNP'){
      if(length(SNPlist)==1){
        d$labl1 <- ''
        d$labl1[match(SNPlist[[1]], d$SNP)] <- SNPlist[[1]]
        plot <- plot +
          geom_text_repel(aes(label=d$labl1), size = 5, color='black',ylim=c(5,ymax))
      }else{
        d$labl1 <- ''
        d$labl1[match(SNPlist[[1]], d$SNP)] <- SNPlist[[1]]
        d$labl2 <- ''
        d$labl2[match(SNPlist[[2]], d$SNP)] <- SNPlist[[2]]
        plot <- plot +
          geom_text_repel(aes(label=d$labl1), size = 5, color='black',ylim=lab2lim,min.segment.length=0.1)+
          geom_text_repel(aes(label=d$labl2), size = 5, color='orange',ylim=lab3lim,segment.alpha=0.4,min.segment.length=0.1)
      }
    }else if(annotate.fld == 'Gene'){
      if(length(SNPlist)==1){
        d$labl1 <- ''
        d$labl1[match(SNPlist[[1]], d$SNP)] <- d$Gene[match(SNPlist[[1]], d$SNP)]
        plot <- plot +
          geom_text_repel(aes(label=d$labl1), size =5, color='black',ylim=c(5,ymax))
      }else{
        d$labl1 <- ''
        d$labl1[match(SNPlist[[1]], d$SNP)] <- d$Gene[match(SNPlist[[1]], d$SNP)]
        d$labl2 <- ''
        d$labl2[match(SNPlist[[2]], d$SNP)] <- d$Gene[match(SNPlist[[2]], d$SNP)]
        plot <- plot +
          geom_text_repel(aes(label=d$labl1), size =5, color='black',ylim=lab2lim,min.segment.length=0.1)+
          geom_text_repel(aes(label=d$labl2), size = 5, color='orange',ylim=lab3lim,segment.alpha=0.4,min.segment.length=0.1)
      }
    }
  }

  if (!is.null(suggestiveline))
    plot <- plot + geom_hline(yintercept = suggestiveline, colour="brown", alpha = 0.33)
  if (!is.null(genomewideline))
    plot <- plot + geom_hline(yintercept = genomewideline, colour = "brown")
  if(!is.null(ymax) && is.numeric(ymax)){
    if(ymax <= max(d$logp) & ymax >= 0)
      plot <- plot + ylim(0, ymax)
  }
  plot
}


## manhattan plot
## used when y-axis should be in the regular scale
gg.manhattan2 = function(dataframe, title=NULL, cols = c('red', 'skyblue'),
                         chr.fld = 'CHR', pos.fld = 'BP', p.fld = 'P', snp.fld = 'rsid', gene.fld=NULL,
                         chrs = 1:22, ymax = NULL,
                         suggestiveline=0, genomewideline=-log10(5e-8),
                         annotate=F, annotate.fld='SNP', SNPlist=NULL, lab2lim=NULL, lab3lim=NULL) {
  library(data.table)
  library(ggplot2)
  library(dplyr)

  if (annotate & is.null(SNPlist))
    stop("You requested annotation but provided no SNPlist!")
  if(annotate & any(! unlist(SNPlist) %in% dataframe[[snp.fld]]))
    stop('No SNPs in SNPlist are found in the dataframe')

  if(is.data.table(dataframe)) dataframe <- as.data.frame(dataframe)

  stopifnot(all(c(chr.fld, pos.fld, p.fld, snp.fld) %in% names(dataframe)),
            is.data.frame(dataframe),
            all(!is.na(dataframe[[chr.fld]])),
            all(unlist(lapply(dataframe[c(chr.fld, pos.fld, p.fld)], is.numeric))))
  if(is.null(gene.fld)){
    d <- dataframe[c(chr.fld, pos.fld, p.fld, snp.fld)]
    names(d) <- c('CHR', 'BP', 'P', 'SNP')
  }else{
    d <- dataframe[c(chr.fld, pos.fld, p.fld, snp.fld, gene.fld)]
    names(d) <- c('CHR', 'BP', 'P', 'SNP','Gene')
  }

  #limit to only chrs 1-22
  d <- d %>% filter(CHR %in% chrs, !is.na(BP)) %>% arrange(CHR, BP)

  # get position references for each CHR:
  d$CHR <- factor(d$CHR)
  posDat <- d %>% group_by(CHR) %>% dplyr::summarise(maxBP = max(BP))

  tmp <- data.frame(0, 0); names(tmp) <- names(posDat)
  posDat <- rbind(tmp, posDat)
  posDat$cumBP <- sapply(1:nrow(posDat), function(x) sum(posDat$maxBP[1:x]))
  posDat <- data.frame(CHR = posDat$CHR[-1], addBP = posDat$cumBP[-nrow(posDat)])

  d <- merge(d, posDat, all.x = T)

  d <- d %>% mutate(pos = BP + addBP)
  ticks <- d %>% group_by(CHR) %>% dplyr::summarise(x = median(pos))

  mycols <- rep(cols, length(chrs))[1:length(unique(chrs))]

  if(length(unique(chrs)) == 1){
    plot <- ggplot(d, aes(x = pos, y = P, color = CHR)) +
      xlab(paste('Chromosome', unique(chrs)))
  } else {
    plot <- ggplot(d, aes(x = pos, y = P, color = CHR)) +
      scale_x_continuous(name = 'Chromosome', breaks = ticks$x, labels = ticks$CHR)
  }
  upper <- ceiling(max(d$P))
  if(!is.null(ymax) && is.numeric(ymax)){
    upper <- ceiling(ymax)
  }
  plot <- plot + geom_point() + scale_color_manual(values = mycols) +
    ylab(expression(-log[10](italic(p)))) +
    ggtitle(title) +
    theme_bw() +
    theme(axis.text.y = element_text(size = rel(2), color = 'black'),
          axis.text.x = element_text(size = rel(1.5), color = 'black'),
          legend.position = 'none', plot.title = element_text(size = rel(5), hjust = 0.5),
          panel.grid.minor = element_blank())

  if(annotate){
    if(annotate.fld=='SNP'){
      if(length(SNPlist)==1){
        d$labl1 <- ''
        d$labl1[match(SNPlist[[1]], d$SNP)] <- SNPlist[[1]]
        plot <- plot +
          geom_text_repel(aes(label=d$labl1), size = 6, color='black')
      }else{
        d$labl1 <- ''
        d$labl1[match(SNPlist[[1]], d$SNP)] <- SNPlist[[1]]
        d$labl2 <- ''
        d$labl2[match(SNPlist[[2]], d$SNP)] <- SNPlist[[2]]
        d$labl3 <- ''
        d$labl3[match(SNPlist[[3]], d$SNP)] <- SNPlist[[3]]
        plot <- plot +
          geom_text_repel(aes(label=d$labl1), size = 6, color='black')+
          geom_text_repel(aes(label=d$labl2), size = 6, color='orange',segment.alpha=0.4,ylim=lab2lim)+
          geom_text_repel(aes(label=d$labl3), size = 6, color='orange',segment.alpha=0.4,ylim=lab3lim)
      }
    }else if(annotate.fld == 'Gene'){
      if(length(SNPlist)==1){
        d$labl1 <- ''
        d$labl1[match(SNPlist[[1]], d$SNP)] <- d$Gene[match(SNPlist[[1]], d$SNP)]
        plot <- plot +
          geom_text_repel(aes(label=d$labl1), size = 6, color='black',box.padding = 1,force=1)
      }else{
        d$labl1 <- ''
        d$labl1[match(SNPlist[[1]], d$SNP)] <- d$Gene[match(SNPlist[[1]], d$SNP)]
        d$labl2 <- ''
        d$labl2[match(SNPlist[[2]], d$SNP)] <- d$Gene[match(SNPlist[[2]], d$SNP)]
        d$labl3 <- ''
        d$labl3[match(SNPlist[[3]], d$SNP)] <- d$Gene[match(SNPlist[[3]], d$SNP)]
        plot <- plot +
          geom_text_repel(aes(label=d$labl1), size = 6, color='black')+
          geom_text_repel(aes(label=d$labl2), size = 6, color='orange',segment.alpha=0.4,ylim=lab2lim)+
          geom_text_repel(aes(label=d$labl3), size = 6, color='orange',segment.alpha=0.4,ylim=lab3lim)
      }
    }
  }

  if (!is.null(suggestiveline))
    plot <- plot + geom_hline(yintercept = suggestiveline, colour="brown", alpha = 0.33)
  if (!is.null(genomewideline))
    plot <- plot + geom_hline(yintercept = genomewideline, colour = "brown")
  if(!is.null(ymax) && is.numeric(ymax)){
    if(ymax <= max(d$P) & ymax >= 0)
      plot <- plot + ylim(0, ymax)
  }
  plot
}

## qq plot
gg.qq = function(pvector, title=NULL, size.labels=2) {
  library(ggplot2)
  pvector <- pvector[!is.na(pvector) & pvector<1 & pvector>0]
  o <- -log10(sort(pvector, decreasing = F))
  e <- -log10( ppoints(length(o) ))
  plot <- qplot(e, o, xlim=c(0,max(e)), ylim=c(0,max(o))) +
    geom_abline(intercept=0,slope=1, col="red") +
    ggtitle(title) +
    scale_x_continuous(name = expression(Expected~~-log[10](italic(p)))) +
    scale_y_continuous(name = expression(Observed~~-log[10](italic(p)))) +
    theme_bw() +
    theme(axis.text = element_text(size = rel(1.5), color = 'black'),
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(size.labels), color = 'black'))
  plot
}


gg.chr.scatter = function(dataframe, cols = c('#ca0020', '#abd9e9'),
                          chr.fld1 = 'chr', pos.fld1 = 'pos', feature.fld = 'pair',
                          chr.fld2 = 'chromosome_name', pos.fld2 = 'middle_position',
                          loc.fld = 'location',gene.fld='Exprs', n.thresh=500,
                          chrs = 1:22){
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)

  if(is.data.table(dataframe)) dataframe <- as.data.frame(dataframe)
  d <- dataframe[c(feature.fld, chr.fld1, pos.fld1, chr.fld2, pos.fld2, loc.fld, gene.fld)]
  names(d) <- c('Feature','CHR', 'BP', 'chr', 'bp','location','Gene')


  #limit to only chrs 1-22
  d <- d %>% filter(CHR %in% chrs, chr %in% chrs, !is.na(BP), !is.na(bp))
  summ = d %>% group_by(Gene) %>% dplyr::summarize(N=length(Gene))
  hits = summ %>% filter(N>=n.thresh) %>% select(Gene) %>% unlist()
  d$Gene[!d$Gene %in% hits]=''
  d$Gene[duplicated(d$Gene)]=''


  # get position references for each CHR:
  d$CHR <- factor(d$CHR)
  d$chr <- factor(d$chr)
  posDat1 <- d %>% group_by(CHR) %>% dplyr::summarise(maxBP = max(BP))
  posDat2 <- d %>% group_by(chr) %>% dplyr::summarise(maxbp = max(bp))
  posDat=merge(posDat1,posDat2,by.x='CHR',by.y='chr',all=T)
  posDat=apply(posDat,2,as.numeric) %>% as.data.frame %>% arrange(CHR)
  tmp <- data.frame(0, 0, 0); names(tmp) <- c('CHR','maxBP','maxbp')
  posDat <- rbind(tmp, posDat)
  posDat$cumBP <- sapply(1:nrow(posDat), function(x) sum(posDat$maxBP[1:x]))
  posDat$cumbp <- sapply(1:nrow(posDat), function(x) sum(posDat$maxbp[1:x]))
  posDat <- data.frame(CHR = posDat$CHR[-1], addBP = posDat$cumBP[-nrow(posDat)], addbp=posDat$cumbp[-nrow(posDat)])

  d <- merge(d, posDat %>% dplyr::select(CHR, addbp), by.x='chr',by.y='CHR', all.x = T)
  d <- merge(d, posDat %>% dplyr::select(CHR, addBP), by.x='CHR',by.y='CHR', all.x = T)

  d <- d %>% mutate(Pos = BP + addBP, pos = bp + addbp)
  Ticks <- d %>% group_by(CHR) %>% dplyr::summarise(x = median(Pos))
  ticks <- d %>% group_by(chr) %>% dplyr::summarise(x = median(pos))
  #annot = d %>% filter(Gene %in% hits) %>% dplyr::select(Gene, chr, bp, pos) %>% distinct() %>% mutate(xpos=0)

  plot <- ggplot(d, aes(x = Pos, y = pos, color=location, label=Gene)) +
    scale_x_continuous(name = 'Chromosome (Methylation)', breaks = Ticks$x, labels = Ticks$CHR) +
    scale_y_continuous(name = 'Chromosome (Gene Expression)', breaks = ticks$x, labels = ticks$chr)

  plot <- plot +
    geom_point(size=0.2)+
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = 'black'),
          axis.title = element_text(size = 12, color = 'black'),
          legend.position = 'none')+
    scale_color_manual(values=cols)+
    #geom_text(aes(x=0), color='black',size=5)
    geom_text_repel(aes(x=100),color='black',force=0.1, size=5)
  #annotate('text', x=0, y=annot$pos, label=annot$Gene, size=3)
  plot
}


# enrichment using cluster profiler
enrich.analysis = function(genes, ref){
  library(clusterProfiler)
  library(DOSE)
  ids = bitr(as.character(na.omit(unique(genes))), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  Ref = bitr(as.character(na.omit(unique(ref))), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")

  cat('\nGOMF...\n')
  GOMF = enrichGO(gene=ids$ENTREZID, OrgDb="org.Hs.eg.db", universe=Ref$ENTREZID, ont = "MF", keyType = 'ENTREZID', pvalueCutoff = 0.1,qvalueCutoff=1)
  GOMF = setReadable(GOMF, "org.Hs.eg.db", "ENTREZID")
  cat('\nGOMF finished!\n')
  cat('\nGOCC...\n')
  GOCC = enrichGO(gene=ids$ENTREZID, OrgDb="org.Hs.eg.db", universe=Ref$ENTREZID, ont = "CC", keyType = 'ENTREZID',pvalueCutoff = 0.1,qvalueCutoff=1)
  GOCC = setReadable(GOCC, "org.Hs.eg.db", "ENTREZID")
  cat('\nGOCC finished!\n')
  cat('\nGOBP...\n')
  GOBP = enrichGO(gene=ids$ENTREZID, OrgDb="org.Hs.eg.db", universe=Ref$ENTREZID, ont = "BP", keyType = 'ENTREZID',pvalueCutoff = 0.1,qvalueCutoff=1)
  GOBP = setReadable(GOBP, "org.Hs.eg.db", "ENTREZID")
  cat('\nGOBP finished!\n')
  cat('\nKEGG...\n')
  KEGG = enrichKEGG(gene=ids$ENTREZID, organism = "hsa", universe=Ref$ENTREZID, keyType = "kegg",pvalueCutoff = 0.1,qvalueCutoff=1)
  KEGG = setReadable(KEGG, "org.Hs.eg.db", "ENTREZID")
  cat('\nKEGG finished!\n')
  cat('\nDO...\n')
  DO = enrichDO(gene=ids$ENTREZID, ont = "DO", universe=Ref$ENTREZID,pvalueCutoff = 0.1,qvalueCutoff=1, readable=T)
  DO = setReadable(DO, "org.Hs.eg.db", "ENTREZID")
  cat('\nDO finished!\n')

  enrichments = list(GOBP, GOMF, GOCC, KEGG, DO)
  names(enrichments)=c('GOBP', 'GOMF', 'GOCC', 'KEGG', 'DO')

  return(enrichments)
  cat('\nDone!\n')
}

enrich.plot = function(enrichres, method=c('top','FDR'),FDRcutoff=0.05,topNum=20){
  cat('\nNow plotting...\n')
  if(method=='FDR'){
    GOBP_plot = ggplot(enrichres$GOBP@result %>% dplyr::filter(p.adjust<FDRcutoff),
                       aes(x=reorder(Description, -log10(p.adjust)), y=-log10(p.adjust), fill='red')) +
      geom_bar(stat="identity")+
      xlab('GO Biological Process')+
      ylab(expression(paste('- log'[10],' BH-FDR')))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = 'none')+
      coord_flip()
    GOMF_plot = ggplot(enrichres$GOMF@result %>% dplyr::filter(p.adjust<FDRcutoff),
                       aes(x=reorder(Description, -log10(p.adjust)), y=-log10(p.adjust), fill='red')) +
      geom_bar(stat="identity")+
      xlab('GO Molecular Function')+
      ylab(expression(paste('- log'[10],' BH-FDR')))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = 'none')+
      coord_flip()
    GOCC_plot = ggplot(enrichres$GOCC@result %>% dplyr::filter(p.adjust<FDRcutoff),
                       aes(x=reorder(Description, -log10(p.adjust)), y=-log10(p.adjust), fill='red')) +
      geom_bar(stat="identity")+
      xlab('GO Cellular Component')+
      ylab(expression(paste('- log'[10],' BH-FDR')))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = 'none')+
      coord_flip()
    KEGG_plot = ggplot(enrichres$KEGG@result %>% dplyr::filter(p.adjust<FDRcutoff),
                       aes(x=reorder(Description, -log10(p.adjust)), y=-log10(p.adjust), fill='red')) +
      geom_bar(stat="identity")+
      xlab('KEGG Pathway')+
      ylab(expression(paste('- log'[10],' BH-FDR')))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = 'none')+
      coord_flip()
    DO_plot = ggplot(enrichres$DO@result %>% dplyr::filter(p.adjust<FDRcutoff),
                     aes(x=reorder(Description, -log10(p.adjust)), y=-log10(p.adjust), fill='red')) +
      geom_bar(stat="identity")+
      xlab('Disease Ontology')+
      ylab(expression(paste('- log'[10],' BH-FDR')))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = 'none')+
      coord_flip()
  }else if(method=='top'){
    GOBP_plot = ggplot(enrichres$GOBP@result %>% arrange(p.adjust) %>% .[1:min(nrow(.),topNum),],
                       aes(x=reorder(Description, -log10(p.adjust)), y=-log10(p.adjust), fill='red')) +
      geom_bar(stat="identity")+
      xlab('GO Biological Process')+
      ylab(expression(paste('- log'[10],' BH-FDR')))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = 'none')+
      coord_flip()
    GOMF_plot = ggplot(enrichres$GOMF@result %>% arrange(p.adjust) %>% .[1:min(nrow(.),topNum),],
                       aes(x=reorder(Description, -log10(p.adjust)), y=-log10(p.adjust), fill='red')) +
      geom_bar(stat="identity")+
      xlab('GO Molecular Function')+
      ylab(expression(paste('- log'[10],' BH-FDR')))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = 'none')+
      coord_flip()
    GOCC_plot = ggplot(enrichres$GOCC@result %>% arrange(p.adjust) %>% .[1:min(nrow(.),topNum),],
                       aes(x=reorder(Description, -log10(p.adjust)), y=-log10(p.adjust), fill='red')) +
      geom_bar(stat="identity")+
      xlab('GO Cellular Component')+
      ylab(expression(paste('- log'[10],' BH-FDR')))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = 'none')+
      coord_flip()
    KEGG_plot = ggplot(enrichres$KEGG@result %>% arrange(p.adjust) %>% .[1:min(nrow(.),topNum),],
                       aes(x=reorder(Description, -log10(p.adjust)), y=-log10(p.adjust), fill='red')) +
      geom_bar(stat="identity")+
      xlab('KEGG Pathway')+
      ylab(expression(paste('- log'[10],' BH-FDR')))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = 'none')+
      coord_flip()
    DO_plot = ggplot(enrichres$DO@result %>% arrange(p.adjust) %>% .[1:min(nrow(.),topNum),],
                     aes(x=reorder(Description, -log10(p.adjust)), y=-log10(p.adjust), fill='red')) +
      geom_bar(stat="identity")+
      xlab('Disease Ontology')+
      ylab(expression(paste('- log'[10],' BH-FDR')))+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.position = 'none')+
      coord_flip()
  }

  plots = list(GOBP_plot, GOMF_plot, GOCC_plot, KEGG_plot, DO_plot)
  names(plots)=c('GOBP', 'GOMF', 'GOCC', 'KEGG', 'DO')
  return(plots)
  cat('\nDone!\n')
}


gsea.analysis = function(genes){
  library(clusterProfiler)
  library(DOSE)
  cnvt = bitr(as.character(names(genes)), fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  genes = genes[names(genes) %in% cnvt$SYMBOL]
  names(genes) = as.character(cnvt$ENTREZID[match(names(genes), cnvt$SYMBOL)])

  cat('\nGOMF...\n')
  GOMF = gseGO(gene=genes, OrgDb="org.Hs.eg.db", ont = "MF", keyType = 'ENTREZID', pvalueCutoff = 0.1)
  GOMF = setReadable(GOMF, "org.Hs.eg.db", "ENTREZID")
  cat('\nGOMF finished!\n')
  cat('\nGOCC...\n')
  GOCC = gseGO(gene=genes, OrgDb="org.Hs.eg.db", ont = "CC", keyType = 'ENTREZID',pvalueCutoff = 0.1)
  GOCC = setReadable(GOCC, "org.Hs.eg.db", "ENTREZID")
  cat('\nGOCC finished!\n')
  cat('\nGOBP...\n')
  GOBP = gseGO(gene=genes, OrgDb="org.Hs.eg.db", ont = "BP", keyType = 'ENTREZID',pvalueCutoff = 0.1)
  GOBP = setReadable(GOBP, "org.Hs.eg.db", "ENTREZID")
  cat('\nGOBP finished!\n')
  cat('\nKEGG...\n')
  KEGG = gseKEGG(gene=genes, organism = "hsa", keyType = "kegg",pvalueCutoff = 0.1)
  KEGG = setReadable(KEGG, "org.Hs.eg.db", "ENTREZID")
  cat('\nKEGG finished!\n')
  cat('\nDO...\n')
  DO = gseDO(gene=genes,pvalueCutoff = 0.1)
  DO = setReadable(DO, "org.Hs.eg.db", "ENTREZID")
  cat('\nDO finished!\n')

  enrichments = list(GOBP, GOMF, GOCC, KEGG, DO)
  names(enrichments)=c('GOBP', 'GOMF', 'GOCC', 'KEGG', 'DO')

  return(enrichments)
  cat('\nDone!\n')
}





## function for colocalisation with COLOC
## gwas1.region.dat: a region of gwas1
## gwas2: full gwas2
## need to have columns [region_id, region_name, region_pos, chr, bp, A1, A2, MAF, Effect, StdErr, Pvalue, N] in gwas1
## need to have columns [rsid, chr, bp, A1, A2, beta, se, p, maf, n] in gwas2
## eg. for xQTL studies, region_id = probes, region_name = gene symbol, region_pos = probe position
coloc.region = function(gwas1.region.dat, gwas2, ...){
  test1 = gwas1.region.dat

  test1 = test1 %>% dplyr::select(region_id, region_name, region_pos, chr, bp, A1, A2, MAF, Effect, StdErr, Pvalue, N)
  test2 = test1
  test2$A1 = test1$A2
  test2$A2 = test1$A1
  test2$Effect = -test1$Effect
  test = rbind(test1, test2)

  d_coloc <- inner_join(gwas2 %>% select(rsid, chr, bp, A1, A2, beta, se, p, maf, n),
                        test)
  # data for coloc
  d_gwas1 <- list(snp = d_coloc$rsid, N = d_coloc$N, MAF = d_coloc$MAF,
                  beta = d_coloc$Effect, varbeta = d_coloc$StdErr^2,
                  type = 'quant')
  d_gwas2 <- list(snp = d_coloc$rsid, N = d_coloc$n, MAF = d_coloc$maf,
                  beta = d_coloc$beta, varbeta = d_coloc$se^2,
                  type = 'quant') ## replace the sample size when available

  # coloc
  abf <- try(coloc.abf(d_gwas1, d_gwas2, ...),silent = T)
  if(class(abf) != "try-error") {
    abf_summ <- as.data.frame(t(abf$summary))
    # most possible shared snp
    abf_snp <- arrange(abf$results, desc(SNP.PP.H4)) %>% slice(1) %>% with(snp)
    abf_snp <- filter(d_coloc, rsid == abf_snp) %>%
      mutate(coloc_snp_gwas1 = sprintf('%.2g (P=%.2g)', Effect, Pvalue),
             coloc_snp_gwas2 = sprintf('%.2g (P=%.2g)', beta, p)) %>%
      select(region_name, region_pos, coloc_region = region_id, rsid, chr, bp, A1, A2, coloc_snp_gwas1, coloc_snp_gwas2)
  } else {
    abf_summ <- data.frame(nsnps = NA, PP.H0.abf = NA, PP.H1.abf = NA, PP.H2.abf=NA, PP.H3.abf = NA, PP.H4.abf=NA)
    abf_snp <- data.frame(region_name=NA, region_pos=NA, coloc_region = NA, rsid = NA, chr=NA, bp=NA, effect_allele=NA, non_effect_allele=NA, coloc_snp_gwas1 = NA, coloc_snp_gwas2 = NA)
  }

  cbind(select(abf_summ, nsnps), abf_snp, select(abf_summ, -nsnps))
}


# region plot for how two traits colocalize with each other
# parameter requirements follow the above coloc.region function
# top.snp is a character string, which is the top coloc SNP in the region (what you want to label)
plot.specific.probe = function(gwas1.region.dat, gwas2, top.snp, mycolor=c('#5e4fa2','#fdae61'),trait1='t1',trait2='t2',trait1LayerAtTop=F){
  test1 = gwas1.region.dat

  test1 = test1 %>% dplyr::select(region_id, region_name, region_pos, chr, bp, A1, A2, MAF, Effect, StdErr, Pvalue, N)
  test2 = test1
  test2$A1 = test1$A2
  test2$A2 = test1$A1
  test2$Effect = -test1$Effect
  test = rbind(test1, test2)

  d_coloc <- inner_join(gwas2 %>% select(rsid, chr, bp, A1, A2, beta, se, p, maf, n), test)

  plotdat = d_coloc %>% select(chr, bp, rsid, Pvalue, p) %>% setnames(c('Pvalue', 'p'), c('tr1','tr2')) %>%
    arrange(chr,bp) %>% mutate(pos=as.numeric(as.factor(bp)), labl=ifelse(rsid==top.snp,top.snp,''))
  tmp = plotdat %>% filter(labl==top.snp)
  axisConvFctr =  (-log10(tmp$tr1))/(-log10(tmp$tr2))

  p=ggplot(plotdat, aes(x=pos))

  if(trait1LayerAtTop){
    p=p+
      geom_point(aes(y=-log10(tr2)*axisConvFctr, color=trait2),alpha=0.4)+
      geom_point(aes(y=-log10(tr1), color=trait1),alpha=0.4)
  }else{
    p=p+
      geom_point(aes(y=-log10(tr1), color=trait1),alpha=0.4)+
      geom_point(aes(y=-log10(tr2)*axisConvFctr, color=trait2),alpha=0.4)
  }
  p=p+
    scale_y_continuous(sec.axis = sec_axis(~./axisConvFctr,
                                           name = eval(bquote(expression(.(trait2) ~ paste(' [-log'[10],'P]'))))))+
    geom_text_repel(aes(label=labl,y=-log10(tr1)), size =5, color='black',min.segment.length=0.1,point.padding=0.7)+
    theme_bw()+
    scale_color_manual(values = mycolor)+
    theme(legend.position='top',
          axis.text = element_text(size = 13, color = 'black'),
          axis.title = element_text(size = 15, color = 'black'),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17))+
    labs(y=eval(bquote(expression(.(trait1) ~ paste(' [-log'[10],'P]')))),
         x='Genomic Position',
         color='')

  return(p)

}





## function for colocalisation with HyPrColoc
## gwas1.region.dat: a region of gwas1
## gwas.dat: dataset containing all the results of all the gwases
## traits: vector containing the trait of each gwas (starting with the trait of gwas1)
## need to have columns [region_id, region_name, region_pos, chr, bp, A1, A2, Effect, StdErr] in gwas1
## need to have columns [rsid, chr, bp, A1, A2, beta_1, se_1, beta_2, se_2, ..., beta_*, se_*] in gwas.dat
## eg. for eQTL studies, region_id = probes, region_name = gene symbol, region_pos = probe position
hyprcoloc.region.multi = function(gwas1.region.dat, gwas.dat, traits, ...){
  test1 = gwas1.region.dat

  test1 = test1 %>% dplyr::select(region_id, chr, bp, A1, A2, Effect, StdErr)
  test2 = test1
  test2$A1 = test1$A2
  test2$A2 = test1$A1
  test2$Effect = -test1$Effect
  test = rbind(test1, test2)

  d_coloc <- merge(gwas.dat,
                   test, by.x=c('chr','bp','A1','A2'),
                   by.y=c('chr','bp','A1','A2')) %>%
    select(chr,bp,A1,A2,rsid,starts_with('beta_'),starts_with('se_'),Effect,StdErr) %>%
    filter(complete.cases(.))
  # data for coloc
  betas = d_coloc %>% select(Effect, starts_with('beta_')) %>% setnames(colnames(.),traits) %>% as.matrix()
  rownames(betas) = d_coloc$rsid
  ses = d_coloc %>% select(StdErr, starts_with('se_')) %>% setnames(colnames(.),traits) %>% as.matrix()
  rownames(ses) = d_coloc$rsid

  marker = rownames(betas)
  out <- try(hyprcoloc(betas, ses, trait.names=traits, snp.id=marker, ...),silent = T)

  if(class(out) != "try-error") {
    res = out$results
  }else{
    res = data.frame(iteration = NA, traits = NA, posterior_prob = NA, regional_prob=NA,
                     candidate_snp = NA, posterior_explained_by_snp=NA,dropped_trait=NA)
  }

  print(paste0('done for ',unique(gwas1.region.dat$region_id),'!!'))
  return(cbind(res, region_id=unique(gwas1.region.dat$region_id),
               region_name=unique(gwas1.region.dat$region_name),
               region_pos=unique(gwas1.region.dat$region_pos)))
}


# make an upset plot in intersection mode
#listInput <- list(`TWAS (IPF)` = gene_twas, `Coloc (IPF)` = gene_coloc,
#				`SMR (IPF)` = gene_smr, `HyPrColoc (IPF)` = gene_hyprcoloc,
#				`TWAS (GTEx Lung)` = gene_twas_GTEx, `Coloc (GTEx Lung)` = gene_coloc_GTEx,
#				`SMR (GTEx Lung)` = gene_smr_GTEx, `HyPrColoc (GTEx Lung)` = gene_hyprcoloc_GTEx)
#m = make_comb_mat(listInput, mode = "intersect")
#df = m %>% t() %>% as.data.frame()
#df$freq=comb_size(m)
#df = df %>% uncount(freq)
#png('figs/IPF_GTEx_res_gene_overlap_intersectMode.png', res=600, width=6000,height=3000)
#upset(df,nsets = 4, nintersects = 400,
#sets = c("HyPrColoc (GTEx Lung)","SMR (GTEx Lung)","Coloc (GTEx Lung)","TWAS (GTEx Lung)",
#"HyPrColoc (IPF)", "SMR (IPF)", "Coloc (IPF)","TWAS (IPF)"),
#empty.intersections = "on",
#mainbar.y.label = "Intersection Size", sets.x.label = "Number of significant genes",keep.order = TRUE)
#dev.off()
