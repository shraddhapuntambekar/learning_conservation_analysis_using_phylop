FolderPaths <- c('out_map_roi/all/ACC/', 
                 'out_map_roi/all/CONACC/',
                 'out_map_roi/ON/ACC/',
                 'out_map_roi/ON/CONACC/',
                 'out_map_roi/PN/ACC/',
                 'out_map_roi/PN/CONACC/')
OutputFigureName <- c('ONmod_Allbranch_ACC','ONmod_Allbranch_CONACC',
                      'ONmod_ONbranch_ACC','ONmod_ONbranch_CONACC',
                      'ONmod_PNbranch_ACC','ONmod_PNbranch_CONACC') 

k<-2

statisticalTest <- c('welch')  ## welch or wilcox

dash3UTR <- data.frame(read.delim(paste0(FolderPaths[k],'meanAndSd_3dashUTR_ON_id_sorted_perbase.bed'), sep = '\t', col.names = c('chr', 'start','end','mean','sd'), as.is = TRUE, header = FALSE))
dash5UTR <- data.frame(read.delim(paste0(FolderPaths[k],'meanAndSd_5dashUTR_ON_id_sorted_perbase.bed'), sep = '\t', col.names = c('chr', 'start','end','mean','sd'), as.is = TRUE, header = FALSE))
introns <- data.frame(read.delim(paste0(FolderPaths[k],'meanAndSd_introns_ON_id_sorted_perbase.bed'), sep = '\t', col.names = c('chr', 'start','end','mean','sd'), as.is = TRUE, header = FALSE))
exons <- data.frame(read.delim(paste0(FolderPaths[k],'meanAndSd_exons_ON_id_sorted_perbase.bed'), sep = '\t', col.names = c('chr', 'start','end','mean','sd'), as.is = TRUE, header = FALSE))
cds <- data.frame(read.delim(paste0(FolderPaths[k],'meanAndSd_codingExons_ON_id_sorted_perbase.bed'), sep = '\t', col.names = c('chr', 'start','end','mean','sd'), as.is = TRUE, header = FALSE))
on_novelIntergenic <- data.frame(read.delim(paste0(FolderPaths[k],'meanAndSd_all_assemblies_ON_intergenic_sorted_perbase.bed'), sep = '\t', col.names = c('chr', 'start','end','mean','sd'), as.is = TRUE, header = FALSE))
on_novelIntronic <-data.frame(read.delim(paste0(FolderPaths[k],'meanAndSd_all_assemblies_ON_intronic_sorted_perbase.bed'), sep = '\t', col.names = c('chr', 'start','end','mean','sd'), as.is = TRUE, header = FALSE))
intergenes <- data.frame(read.delim(paste0(FolderPaths[k],'meanAndSd_ON_intergenic_complementof_wholegene_and_chrmSizes_sorted_perbase.bed'), sep = '\t', col.names = c('chr', 'start','end','mean','sd'), as.is = TRUE, header = FALSE)) 
ancestralRepeats <- data.frame(read.delim(paste0(FolderPaths[k],'meanAndSd_allAR_OreochromisNiloticus_fromCichlidBrowser_AR_overlap_8waymaf_8April_id_perbase.bed'), sep = '\t', col.names = c('chr', 'start','end','mean','sd'), as.is = TRUE, header = FALSE)) 

if (statisticalTest == 'wilcox'){
## wilcox test
if (k == 2 || k == 4 || k == 6 ) {
outtest_ar_5utr<- wilcox.test(ancestralRepeats$mean,dash3UTR$mean,paired=FALSE,alternative = "two.sided")
outtest_ar_3utr<- wilcox.test(ancestralRepeats$mean,dash5UTR$mean,paired=FALSE,alternative = "two.sided")
outtest_ar_introns<- wilcox.test(ancestralRepeats$mean,introns$mean,paired=FALSE,alternative = "two.sided")
outtest_ar_exons<- wilcox.test(ancestralRepeats$mean,exons$mean,paired=FALSE,alternative = "two.sided")
outtest_ar_cds<- wilcox.test(ancestralRepeats$mean,cds$mean,paired=FALSE,alternative = "two.sided")
outtest_ar_intergenes<- wilcox.test(ancestralRepeats$mean,intergenes$mean,paired=FALSE,alternative = "two.sided")
outtest_ar_novelIntergenic<- wilcox.test(ancestralRepeats$mean,on_novelIntergenic$mean,paired=FALSE,alternative = "two.sided")
outtest_ar_novelIntronic<- wilcox.test(ancestralRepeats$mean,on_novelIntronic$mean,paired=FALSE,alternative = "two.sided")
outtests_pvals <- c(outtest_ar_introns$p.value,outtest_ar_intergenes$p.value,outtest_ar_exons$p.value,outtest_ar_cds$p.value,outtest_ar_5utr$p.value,outtest_ar_3utr$p.value,outtest_ar_novelIntergenic$p.value,outtest_ar_novelIntronic$p.value)
outtest_names <- c('introns','intergenes','exons','cds','5utr','3utr','novelIntergenes','novelIntronic')
outtest <-data.frame("names"=outtest_names,"pval"=outtests_pvals)
outtest11<-outtest
outtest11[outtest$pval < 0.05,]
outtests_pvals1<-outtests_pvals[outtests_pvals<0.05]

outtest_novelIntronic_introns <- wilcox.test(on_novelIntronic$mean,introns$mean,paired=FALSE,alternative = "two.sided")
outtest_novelIntergenic_intergenes <-  wilcox.test(on_novelIntergenic$mean,intergenes$mean,paired=FALSE,alternative = "two.sided")

outtest_novelIntronic_cds <- wilcox.test(on_novelIntronic$mean,cds$mean,paired=FALSE,alternative = "two.sided")
outtest_novelIntergenic_cds <-  wilcox.test(on_novelIntergenic$mean,cds$mean,paired=FALSE,alternative = "two.sided")
}
#------------------------
} else if (statisticalTest == 'welch'){
## Welch t-test
if (k == 2 || k == 4 || k == 6 ) {
  outtest_ar_5utr<- t.test(ancestralRepeats$mean,dash3UTR$mean,paired=FALSE,alternative = "two.sided",var.equal = FALSE)
  outtest_ar_3utr<- t.test(ancestralRepeats$mean,dash5UTR$mean,paired=FALSE,alternative = "two.sided",var.equal = FALSE)
  outtest_ar_introns<- t.test(ancestralRepeats$mean,introns$mean,paired=FALSE,alternative = "two.sided",var.equal = FALSE)
  outtest_ar_exons<- t.test(ancestralRepeats$mean,exons$mean,paired=FALSE,alternative = "two.sided",var.equal = FALSE)
  outtest_ar_cds<- t.test(ancestralRepeats$mean,cds$mean,paired=FALSE,alternative = "two.sided",var.equal = FALSE)
  outtest_ar_intergenes<- t.test(ancestralRepeats$mean,intergenes$mean,paired=FALSE,alternative = "two.sided",var.equal = FALSE)
  outtest_ar_novelIntergenic<- t.test(ancestralRepeats$mean,on_novelIntergenic$mean,paired=FALSE,alternative = "two.sided",var.equal = FALSE)
  outtest_ar_novelIntronic<- t.test(ancestralRepeats$mean,on_novelIntronic$mean,paired=FALSE,alternative = "two.sided",var.equal = FALSE)
  outtests_pvals <- c(outtest_ar_introns$p.value,outtest_ar_intergenes$p.value,outtest_ar_exons$p.value,outtest_ar_cds$p.value,outtest_ar_5utr$p.value,outtest_ar_3utr$p.value,outtest_ar_novelIntergenic$p.value,outtest_ar_novelIntronic$p.value)
  outtest_names <- c('introns','intergenes','exons','cds','5utr','3utr','novelIntergenes','novelIntronic')
  outtest <-data.frame("names"=outtest_names,"pval"=outtests_pvals)
  outtest11<-outtest
  outtest11[outtest$pval < 0.05,]
  outtests_pvals1<-outtests_pvals[outtests_pvals<0.05]
  
  outtest_novelIntronic_introns <- t.test(on_novelIntronic$mean,introns$mean,paired=FALSE,alternative = "two.sided",var.equal = FALSE)
  outtest_novelIntergenic_intergenes <-  t.test(on_novelIntergenic$mean,intergenes$mean,paired=FALSE,alternative = "two.sided",var.equal = FALSE)
  
  outtest_novelIntronic_cds <- t.test(on_novelIntronic$mean,cds$mean,paired=FALSE,alternative = "two.sided",var.equal = FALSE)
  outtest_novelIntergenic_cds <-  t.test(on_novelIntergenic$mean,cds$mean,paired=FALSE,alternative = "two.sided",var.equal = FALSE)
}
}

#----------------------------------
analysis_introns <- data.frame("numacc" = length(introns$mean[introns$mean < 0]) ,"numcon"=length(introns$mean[introns$mean > 0]),"minScore"=min(introns$mean),"maxScore"=max(introns$mean))
analysis_intergenes <- data.frame("numacc" = length(intergenes$mean[intergenes$mean < 0]) ,"numcon"=length(intergenes$mean[intergenes$mean > 0]),"minScore"=min(intergenes$mean),"maxScore"=max(intergenes$mean))
analysis_cds <- data.frame("numacc" = length(cds$mean[cds$mean < 0]) ,"numcon"=length(cds$mean[cds$mean > 0]),"minScore"=min(cds$mean),"maxScore"=max(cds$mean))
analysis_dash3UTR <- data.frame("numacc" = length(dash3UTR$mean[dash3UTR$mean < 0]) ,"numcon"=length(dash3UTR$mean[dash3UTR$mean > 0]),"minScore"=min(dash3UTR$mean),"maxScore"=max(dash3UTR$mean))
analysis_dash5UTR <- data.frame("numacc" = length(dash5UTR$mean[dash5UTR$mean < 0]) ,"numcon"=length(dash5UTR$mean[dash5UTR$mean > 0]),"minScore"=min(dash5UTR$mean),"maxScore"=max(dash5UTR$mean))
analysis_on_novelIntergenic <- data.frame("numacc" = length(on_novelIntergenic$mean[on_novelIntergenic$mean < 0]) ,"numcon"=length(on_novelIntergenic$mean[on_novelIntergenic$mean > 0]),"minScore"=min(on_novelIntergenic$mean),"maxScore"=max(on_novelIntergenic$mean))
analysis_on_novelIntronic <- data.frame("numacc" = length(on_novelIntronic$mean[on_novelIntronic$mean < 0]) ,"numcon"=length(on_novelIntronic$mean[on_novelIntronic$mean > 0]),"minScore"=min(on_novelIntronic$mean),"maxScore"=max(on_novelIntronic$mean))
analysis_AR <- data.frame("numacc" = length(ancestralRepeats$mean[ancestralRepeats$mean < 0]) ,"numcon"=length(ancestralRepeats$mean[ancestralRepeats$mean > 0]),"minScore"=min(ancestralRepeats$mean),"maxScore"=max(ancestralRepeats$mean))
analysis_exons <- data.frame("numacc" = length(exons$mean[exons$mean < 0]) ,"numcon"=length(exons$mean[exons$mean > 0]),"minScore"=min(exons$mean),"maxScore"=max(exons$mean))
analysis <- rbind(analysis_introns,analysis_intergenes,analysis_exons,analysis_cds,analysis_dash5UTR,analysis_dash3UTR,analysis_on_novelIntergenic,analysis_on_novelIntronic,analysis_AR)
analysis1 <- analysis[!is.na(analysis[,3]),]
analysis2 <- analysis1[!is.na(analysis1[,4]),]

min_score <- min(analysis2[,3])
max_score <- max(analysis2[,4])

warnings()
#jpeg(paste0('figures/ONmodel/',OutputFigureName[k],'.jpg'), width = 480, height = 480,units='px')
#pdf(paste0('figures/ONmodel/',OutputFigureName[k],'.pdf'))
plot(ecdf(introns$mean), col = "green",lwd=3,xlim=range(min_score,1),xlab ='', ylab = '',main='', do.points=FALSE, verticals=TRUE)
lines(ecdf(intergenes$mean), col = "brown",lwd=3,xlim=range(min_score,1),xlab ='', ylab = '',main='', do.points=FALSE, verticals=TRUE)
lines(ecdf(ancestralRepeats$mean), col = "black",lwd=3,xlim=range(min_score,1),xlab ='', ylab = '',main='')
lines(ecdf(dash5UTR$mean), col = "darkorange2",lwd=3,xlim=range(min_score,1),xlab ='', ylab = '',main='')
lines(ecdf(dash3UTR$mean), col = "darkmagenta",lwd=3,xlim=range(min_score,1),xlab ='', ylab = '',main='', do.points=FALSE, verticals=TRUE)
lines(ecdf(exons$mean), col = "pink",lwd=3,xlim=range(min_score,1),xlab ='', ylab = '',main='', do.points=FALSE, verticals=TRUE)
lines(ecdf(cds$mean), col = "red",lwd=3,xlim=range(min_score,1),xlab ='', ylab = '',main='', do.points=FALSE, verticals=TRUE)
lines(ecdf(on_novelIntergenic$mean), col = "gray70",lwd=3,xlim=range(min_score,1),xlab ='', ylab = '', do.points=FALSE, verticals=TRUE,main='')
lines(ecdf(on_novelIntronic$mean), col = "blue4",lwd=3,xlim=range(min_score,1),xlab ='', ylab = '', do.points=FALSE, verticals=TRUE,main='')


if (k == 2 || k == 4 || k == 6 ) {
  legend('topleft',c('introns','intergenes','AR','5\'UTR','3\'UTR','exons','CDS','ON-novel-intergenic','ON-novel-intronic'),fill = c('green','brown','black','darkorange2','darkmagenta','pink','red','gray70','blue4'), border=NA)
    mtext(text = expression('CONACC score'), side = 1, line = 2.5,main=paste0(OutputFigureName[k]))
} else {
  legend('bottomright',c('introns','intergenes','AR','5\'UTR','3\'UTR','exons','CDS','ON-novel-intergenic','ON-novel-intronic'),fill = c('green','brown','black','darkorange2','darkmagenta','pink','red','gray70','blue4'), border=NA)
  mtext(text = expression('ACC score'), side = 1, line = 2.5,main=paste0(OutputFigureName[k]))
}    
mtext(text = expression(CDF), side = 2, line = 2.5)

dev.off()


