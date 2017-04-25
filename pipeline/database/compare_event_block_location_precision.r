#!/usr/bin/R
pdf('compare_event_block_location_presision.pdf', width=5, height=5)

difft = read.table('compare_transfer_event_block_location_presision.tab')
colnames(difft) = c("anc_block_id", "event_id", "nb_repsn", "nb_depsn", "nb_rbpsn", "nb_dbpsn")
difft$diff_r = difft$nb_repsn-difft$nb_rbpsn
difft$diff_d = difft$nb_depsn-difft$nb_dbpsn
ht = hist(difft$diff_r[difft$diff_r>=0 & difft$diff_r<=10], breaks=0:11,
 xlab='size reduction of the set of\npossible locations of transfer receiver', ylab='numer of single-gene events',
 main='Precision gain in event location inferrence\nwhen considering block events')
 
dim(difft[difft$nb_depsn>difft$nb_dbpsn,])[1]/dim(difft)[1]
dim(difft[difft$nb_repsn>difft$nb_rbpsn,])[1]/dim(difft)[1]
 
diffd = read.table('compare_duplication_event_block_location_presision.tab')
colnames(diffd) = c("anc_block_id", "event_id", "nb_repsn", "nb_rbpsn")
diffd$diff_r = diffd$nb_repsn-diffd$nb_rbpsn
hd = hist(diffd$diff_r[diffd$diff_r>=0 & diffd$diff_r<=10], breaks=0:11,
 xlab='reduction in size of duplication possible location set', ylab='numer of single-gene events',
 main='Precision gain in event location inferrence when considering block events')
dim(diffd[diffd$nb_repsn>diffd$nb_rbpsn,])[1]/dim(diffd)[1]

coul = c('#FFA620', '#FF0000')
barplot(rbind(ht$counts, hd$counts), names.arg=hd$breaks[1:11], beside=T, col=coul, ylim=c(0,5000),
 xlab='reduction in size of duplication possible location set', ylab='numer of single-gene events',
 main='Precision gain in event location inferrence when considering block events')
legend('topright', c('transfers', 'duplications'), fill=coul)

dev.off()

pdf('block_size_distrib.pdf', width=10, height=5)

trablocsizeditrib = read.table('trans_block_size_ditrib.tab')
dupblocsizeditrib = read.table('dup_block_size_ditrib.tab')
fillzeroclass = function(dataframe, classrange=1:25){
	sapply(classrange, function(s){ if (s %in% dataframe[,1]){ dataframe[dataframe[,1]==s,2] }else{0} })
}
fillzeroclasslog = function(dataframe, classrange=1:25){
	sapply(classrange, function(s){ if (s %in% dataframe[,1]){ log10(dataframe[dataframe[,1]==s,2]) }else{0} })
}
barplot(rbind(fillzeroclasslog(trablocsizeditrib), fillzeroclasslog(dupblocsizeditrib)),
 names.arg=1:25, col=coul, beside=T,
 xlab='block size (max. nb. genes in extant genomes)', ylab='log10(number of block events)',
 main='Distribution of block event sizes')
legend('topright', c('transfers', 'duplications'), fill=coul)
dev.off()
