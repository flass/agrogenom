#!/usr/bin/R
coul = c('#FFA620', '#FF0000')
coulstdg = c('#0000AA', '#FFA620', '#8800CC', '#FF0000', '#00CC00')
cr=1:50
fillzeroclass = function(dataframe, classrange=cr){
	sapply(classrange, function(s){ if (s %in% dataframe[,1]){ dataframe[dataframe[,1]==s,2] }else{0} })
}
poslog = function(x, n=10){
log(x+1, n)
}


###############################

pdf('uncertainty_distribution.pdf', width=6, height=5)

un = read.table('distr_uncertainty.tab', sep='\t')
colnames(un) = c('eventtype', 'carac', 'uncertainty', 'count')
un$logcount = log10(un$count+1)
m = t(sapply(list(c('S','location'), c('T', 'rec'), c('T', 'don'), c('D','location'), c('G','location')), function(l){
e = l[1]
c = l[2]
sapply(1:13, function(u){ un$logcount[un$uncertainty==u & un$eventtype==e & un$carac==c][1]  })
}))
barplot(m, beside=T, names.arg=1:13, col=coulstdg,
 xlab='number of possible locations in the species tree', ylab='number of single-gene events [log10(n+1)]',
 main='Distribution of the degree of uncertainty on\nthe location of events in the species tree')
legend('topright', c('speciations', 'transfers (receiver)', 'transfers (donor)', 'duplications', 'originations'), fill=coulstdg)

dev.off()


###############################

pdf('compare_event_block_location_precision.pdf', width=7, height=5)

difft = read.table('compare_transfer_event_block_location_precision.tab')
colnames(difft) = c("anc_block_id", "event_id", "nb_repsn", "nb_depsn", "nb_rbpsn", "nb_dbpsn", "nb_evt")
# do not count the events that cannot be refined
difft = difft[difft$nb_repsn>1,]
difft$diff_r = difft$nb_repsn-difft$nb_rbpsn
difft$diff_d = difft$nb_depsn-difft$nb_dbpsn
ht = hist(difft$diff_r[difft$diff_r>=0 & difft$diff_r<=10], breaks=0:11,
 xlab='reduction in the number of possible locations in the species tree', ylab='numer of single-gene events',
 main='Gain in precision of the inferrence of the location\nof transfer and duplication events\nwhen considering the regional signal')
 
dim(difft[difft$nb_depsn>difft$nb_dbpsn,])[1]/dim(difft)[1]
dim(difft[difft$nb_repsn>difft$nb_rbpsn,])[1]/dim(difft)[1]
 
diffd = read.table('compare_duplication_event_block_location_precision.tab')
colnames(diffd) = c("anc_block_id", "event_id", "nb_repsn", "nb_rbpsn", "nb_evt")
# do not count the events that cannot be refined
diffd = diffd[diffd$nb_repsn>1,]
diffd$diff_r = diffd$nb_repsn-diffd$nb_rbpsn
hd = hist(diffd$diff_r[diffd$diff_r>=0 & diffd$diff_r<=10], breaks=0:11,
 xlab='reduction in the number of possible locations in the species tree', ylab='number of single-gene events',
 main='Precision gain in event location inferrence when considering block events')
dim(diffd[diffd$nb_repsn>diffd$nb_rbpsn,])[1]/dim(diffd)[1]

#~ barplot(poslog(rbind(ht$counts, hd$counts))[, 2:length(ht$counts)], names.arg=hd$breaks[1:11], beside=T, col=coul, ylim=c(0,5000),
barplot(poslog(rbind(ht$counts, hd$counts)), names.arg=hd$breaks[1:11], beside=T, col=coul, 
 xlab='reduction in the number of possible locations in the species tree', ylab='number of single-gene events [log10(n+1)]',
 main='Gain in precision on the location of events in the species tree\nwhen considering the regional signal')
legend('topright', c('transfers (receiver)', 'duplications'), fill=coul)

dev.off()

###############################

pdf('block_size_distrib.pdf', width=12, height=5)

trablocsizedistrib = read.table('trans_leafblock_size_distrib.tab')
dupblocsizedistrib = read.table('dup_leafblock_size_distrib.tab')
barplot(poslog(rbind(fillzeroclass(trablocsizedistrib), fillzeroclass(dupblocsizedistrib))),
 names.arg=cr, col=coul, beside=T, 
 xlab='Leaf block size (max. nb. genes in extant genomes)', ylab='number of leaf block events [log10(n+1)]',
 main='Distribution of block event sizes')
legend('topright', c('transfers', 'duplications'), fill=coul)

trablocsizedistrib = read.table('trans_ancblock_size_distrib.tab')
dupblocsizedistrib = read.table('dup_ancblock_size_distrib.tab')
oriblocsizedistrib = read.table('ori_ancblock_size_distrib.tab')
noninvetblockevt = read.table('event_non_investigated_for_blocks.tab', row.names=1)

barplot(poslog(rbind(fillzeroclass(trablocsizedistrib), fillzeroclass(dupblocsizedistrib))),
 names.arg=cr, col=coul, beside=T, 
 xlab='Ancestral block size (number of unique gene families)', ylab='number of ancestral block events [log10(n+1)]',
 main='Distribution of sizes of blocks of co-evolved genes
')
legend('topright', c('transfers', 'duplications'), fill=coul)

dev.off()

# single-gene events = 1-event blocks + single events non investigated for block aggregation
#~ > oriblocsizedistrib[1,2]+noninvetblockevt['G',1]
#~ [1] 3600
#~ > dupblocsizedistrib[1,2]+noninvetblockevt['D',1]
#~ [1] 5041
#~ > trablocsizedistrib[1,2]+noninvetblockevt['T',1]
#~ [1] 26606

# multiple-gene events 
#~ > sum(oriblocsizedistrib[2:dim(oriblocsizedistrib)[1],2])
#~ [1] 667
#~ > sum(dupblocsizedistrib[2:dim(dupblocsizedistrib)[1],2])
#~ [1] 778
#~ > sum(trablocsizedistrib[2:dim(trablocsizedistrib)[1],2])
#~ [1] 5649

# total events
#~ > sum(oriblocsizedistrib[,2])+noninvetblockevt['G',1]
#~ [1] 4267
#~ > sum(dupblocsizedistrib[,2])+noninvetblockevt['D',1]
#~ [1] 5819
#~ > sum(trablocsizedistrib[,2])+noninvetblockevt['T',1]
#~ [1] 32255



###################################################################
# compare blocks length and functional homogeneity across phylogeny
###################################################################

tlineages = read.table('/pandata/lassalle/agrogenom/reftree/47RHIZOB.leaf_lineages.tab', sep='\t', stringsAsFactors=F)
lineages = sapply(tlineages[,1], function(cg){ c(cg ,strsplit(tlineages[tlineages[,1]==cg,2], split=' ')[[1]]) })
atustrains = sapply(tlineages[,1], function(cg){ 'N15' %in% lineages[[cg]] })
atulineages = sapply(tlineages[atustrains,1], function(cg){ c(cg ,strsplit(strsplit(tlineages[tlineages[,1]==cg,2], split=' N15 ')[[1]][1], split=' ')[[1]]) })
atugenomovars = c('N35', 'N38', 'N39', 'N42', 'N32', 'N27', 'N29', 'N37')
currgenomes = names(atulineages)
atustrainspecies = sapply(currgenomes, function(cg){ max(atulineages[[cg]][atulineages[[cg]] %in% atugenomovars]) })

distlbl = read.table('distrib_size_leafblocks_by_rec_by_currgenome.tab')
colnames(distlbl) = c('receiver', 'current_genome', 'block_size', 'count')
distlbfh = read.table('distrib_funcsim_leafblocks_by_rec_by_currgenome.tab')
colnames(distlbfh) = c('leafblock_id', 'receiver', 'current_genome', 'block_size', 'funsimP', 'funsimF', 'funsim')
cols=c('blue', 'green', 'red', 'cyan', 'grey', 'orange', 'purple')
#~ currgenomes = names(lineages)
# range of interest of block lengths
N = 8

###############################
# blocks lengths
###############################

pdf('distrib_size_leafblocks_by_rec_by_currgenome.pdf', width=12, height=5)

sapply(currgenomes, function(cg){
 td = distlbl[distlbl$current_genome==cg & (distlbl$receiver %in% lineages[[cg]]),]
#~  N = max(td$block_size)
 td$receiver = factor(td$receiver, levels=lineages[[cg]])
 td$block_size = factor(td$block_size, levels=1:N)
 td = td[!is.na(td$block_size),]
 mat = sapply(lineages[[cg]], function(rec){
  tdr = td[td$receiver==rec,]
  sapply(1:N, function(n){
   ifelse(n %in% tdr$block_size, tdr$count[tdr$block_size==n], 0)
  })
 })
 dmat = sapply(lineages[[cg]], function(rec){ mat[,rec]/sum(mat[,rec]) })
 bp=barplot(dmat, names.arg=lineages[[cg]], beside=T, ylab='proportion of total blocks', main=paste('Distribution of sizes of transferred block received in', cg, 'lineage'))
 mtext(text=c(1,N), at=as.numeric(bp[c(1,N),]), side=1)
 ct = chisq.test(mat)
#~  text(x=1, y=0.9, adj=0, paste('Chi-squared = ', ct$statistic, ', df = ', ct$parameter, ', p-value =', ct$p.value))
 mtext(line=2, side=1, paste('Chi-squared = ', ct$statistic, ', df = ', ct$parameter, ', p-value =', ct$p.value))
 return(ct)
})
dev.off()

###############################

pdf('distrib_size_leafblocks_strain_vs_genomovar.pdf', width=12, height=5)

sapply(currgenomes[!is.na(atustrainspecies)], function(cg){
 anc = atustrainspecies[[cg]]
 strainanc = c(cg, anc)
 td = distlbl[distlbl$current_genome==cg & (distlbl$receiver %in% strainanc),]
#~  N = max(td$block_size)
 td$receiver = factor(td$receiver, levels=strainanc)
 td$block_size = factor(td$block_size, levels=1:N)
 td = td[!is.na(td$block_size),]
 mat = sapply(strainanc, function(rec){
  tdr = td[td$receiver==rec,]
  sapply(1:N, function(n){
   ifelse(n %in% tdr$block_size, tdr$count[tdr$block_size==n], 0)
  })
 })
 print(mat)
 dmat = sapply(strainanc, function(rec){ mat[,rec]/sum(mat[,rec]) })
 bp=barplot(dmat, names.arg=strainanc, beside=T, ylab='proportion of total blocks', main=paste('Distribution of sizes of transferred block\nreceived in', cg, 'and its ancestor', anc))
 mtext(text=c(1,N), at=as.numeric(bp[c(1,N),]), side=1)
 ct = chisq.test(mat)
#~  text(x=1, y=0.9, adj=0, paste('Chi-squared = ', ct$statistic, ', df = ', ct$parameter, ', p-value =', ct$p.value))
 mtext(line=2, side=1, paste('Chi-squared = ', ct$statistic, ', df = ', ct$parameter, ', p-value =', ct$p.value))
 return(ct)
})
dev.off()

###############################

pdf('distrib_size_leafblocks_recent_vs_ancient.pdf', width=12, height=5)

multimat = lapply(currgenomes, function(cg){
 td = distlbl[distlbl$current_genome==cg & (distlbl$receiver %in% lineages[[cg]]),]
#~ N = max(distlbl$block_size)
 mat = t(sapply(lineages[[cg]], function(rec){
  tdr = td[td$receiver==rec,]
  sapply(1:N, function(n){
   ifelse(n %in% tdr$block_size, tdr$count[tdr$block_size==n], 0)
  })
 }))
 facs = data.frame(current_genome=rep(cg, length(lineages[[cg]])), receiver=lineages[[cg]])
 cgmat = cbind(facs, mat)
 colnames(cgmat) = c('current_genome', 'receiver', 1:N)
 return(cgmat)
})
gigamat = do.call(rbind, multimat)

leaves = as.character(gigamat$receiver)==as.character(gigamat$current_genome)
internals = !leaves & !(as.character(gigamat$receiver) %in% c('N1', 'N2', 'N3'))
atus = as.character(gigamat$current_genome[gigamat$receiver=='N15'])
atuleaves = leaves & (as.character(gigamat$current_genome) %in% atus) 
atuinternals = !leaves & (as.character(gigamat$current_genome) %in% atus) & as.character(gigamat$receiver)!='N15'

matleaves = apply(gigamat[leaves,as.character(1:8)], 2, sum)
matinternals = apply(gigamat[internals,as.character(1:8)], 2, sum)
ct = chisq.test(cbind(matleaves, matinternals))
bp=barplot(cbind(matleaves/sum(matleaves), matinternals/sum(matinternals)), beside=T, names.arg=c('leaf_nodes', 'internal_nodes'), ylab='proportion of total blocks', main='Distribution of sizes of transferred block received all Rhizobiaceae (excluding ancestors N1, N2, N3)')
mtext(text=1:8, at=as.numeric(bp), side=1)
text(x=2, y=0.7, adj=0, paste('Chi-squared = ', ct$statistic, ', df = ', ct$parameter, ', p-value =', ct$p.value))
(matleaves/sum(matleaves))/(matinternals/sum(matinternals))

atumatleaves = apply(gigamat[atuleaves,as.character(1:8)], 2, sum)
atumatinternals = apply(gigamat[atuinternals,as.character(1:8)], 2, sum)
ct = chisq.test(cbind(atumatleaves, atumatinternals))
bp=barplot(cbind(atumatleaves/sum(atumatleaves), atumatinternals/sum(atumatinternals)), beside=T, names.arg=c('leaf_nodes', 'internal_nodes'), ylab='proportion of total blocks', main='Distribution of sizes of transferred block received all A. tumefaciens (excluding ancestor N15)')
mtext(text=1:8, at=as.numeric(bp), side=1)
text(x=2, y=0.7, adj=0, paste('Chi-squared = ', ct$statistic, ', df = ', ct$parameter, ', p-value =', ct$p.value))
(atumatleaves/sum(atumatleaves))/(atumatinternals/sum(atumatinternals))
dev.off()

###############################
# blocks functional homogeneity
###############################

pdf('distrib_funcsim_leafblocks_by_rec_by_currgenome.pdf', width=12, height=5)

sapply(currgenomes, function(cg){
#~  td = distlbfh[distlbfh$current_genome==cg & (distlbfh$receiver %in% lineages[[cg]]),]
 td = distlbfh[distlbfh$current_genome==cg & (distlbfh$receiver %in% atulineages[[cg]]),]
 td$receiver = factor(td$receiver, levels=atulineages[[cg]])
 td$block_size = factor(td$block_size, levels=2:N)
 td = td[!is.na(td$block_size),]
 if (dim(td)[1] > (N-1)*length(atulineages[[cg]])*3){
 funsim.aov = aov(lm(funsim ~ receiver * block_size, data=td))
#~  tuk = TukeyHSD(funsim.aov, "receiver", order=T)
 b = boxplot(funsim ~ receiver * block_size, data=td, plot=F)
 boxplot(funsim ~ receiver * block_size, data=td, names=sapply(b$names, function(s){ strsplit(s, split='\\.')[[1]][2] }),
  main=paste('Distribution of functional similarities whithin transferred blocks in', cg),
  xlab="block size", ylab="Schlicker's funsim", col=cols[1:length(atulineages[[cg]])], ylim=c(0,1))
  abline(v=((1:N-2)*length(atulineages[[cg]]))+0.5)
  legend('topright', fill=c(NA,cols[1:length(atulineages[[cg]])]), legend=c('receiver genome', atulineages[[cg]]), bg='white') 
#~  return(tuk[tuk$p<0.05])
 return(summary(funsim.aov))
 }else{return(NA)}
})
dev.off()

###############################

pdf('distrib_funcsim_leafblocks_strain_vs_genomovar.pdf', width=12, height=5)


sapply(currgenomes, function(cg){
 anc = atustrainspecies[[cg]]
 strainanc = c(cg, anc)
#~  td = distlbfh[distlbfh$current_genome==cg & (distlbfh$receiver %in% lineages[[cg]]),]
 td = distlbfh[distlbfh$current_genome==cg & (distlbfh$receiver %in% strainanc),]
 td$receiver = factor(td$receiver, levels=strainanc)
 td$block_size = factor(td$block_size, levels=2:N)
 td = td[!is.na(td$block_size),]
 if (dim(td)[1] > (N-1)*length(strainanc)*3){
  funsim.aov = aov(lm(funsim ~ receiver * block_size, data=td))
#~  tuk = TukeyHSD(funsim.aov, "receiver", order=T)
  b = boxplot(funsim ~ receiver * block_size, data=td, plot=F)
  boxplot(funsim ~ receiver * block_size, data=td, names=sapply(b$names, function(s){ strsplit(s, split='\\.')[[1]][2] }),
   main=paste('Distribution of functional similarities whithin transferred blocks in', cg),
   xlab="block size", ylab="Schlicker's funsim", col=cols[1:length(strainanc)], ylim=c(0,1))
  abline(v=((1:N-2)*length(strainanc))+0.5)
  legend('topright', fill=c(NA,cols[1:length(strainanc)]), legend=c('receiver genome', strainanc), bg='white') 
  mtext(line=4, side=1, adj=1, text=paste(capture.output(summary(funsim.aov))[1:5], collapse='\n'))
#~  return(tuk[tuk$p<0.05])
  return(summary(funsim.aov))
 }else{return(NA)}
})
dev.off()

###############################

pdf('distrib_funcsim_leafblocks_recent_vs_ancient.pdf', width=12, height=5)

N = 12

multimat = lapply(currgenomes, function(cg){
 td = distlbfh[distlbfh$current_genome==cg & (distlbfh$receiver %in% atulineages[[cg]]),]
 td$block_size = factor(td$block_size, levels=2:N, ordered=T)
 td$current_genome = factor(td$current_genome, levels=currgenomes)
 td = td[!is.na(td$block_size),]
})
gigamat = do.call(rbind, multimat)
leaves = as.character(gigamat$receiver)==as.character(gigamat$current_genome)
internals = !leaves & !(as.character(gigamat$receiver) %in% c('N1', 'N2', 'N3'))

gigamat$node = ifelse(leaves, 'leaf', 'internal')
gigamat$receiver = factor(gigamat$receiver, levels=levels(as.factor(do.call(c, atulineages))))
atugigamat = gigamat[!(is.na(gigamat$node) | is.na(gigamat$receiver)),]

#~ t.test(gigamat$funsim[atuleaves], gigamat$funsim[atuinternals])
funsim.aov = aov(lm(funsim ~ 0 + node + block_size, data=atugigamat))
b = boxplot(funsim ~ 0 + node + block_size, data=atugigamat, plot=F)
boxplot(funsim ~ 0 + node + block_size, data=atugigamat, names=sapply(b$names, function(s){ strsplit(s, split='\\.')[[1]][2] }),
 main='Distribution of functional similarities whithin transferred block received along A. tumefaciens history',
 xlab="block size", ylab="Schlicker's funSim", col=cols[1:2], ylim=c(0,1))
#~ N = max(atugigamat$block_size)
abline(v=((1:N-2)*2+0.5))
legend('topleft', fill=c(NA,cols[1:2]), legend=c('receiver genome', c('internal', 'leaves')), bg='white') 
#~ mtext(line=4, side=1, adj=1, text=paste(capture.output(summary(funsim.aov))[1:5], collapse='\n'))
text(x=(N-1)*2, y=0.8, adj=1, labels=paste(capture.output(summary(funsim.aov))[1:5], collapse='\n'))

#~ 
#~ funsim.aov = aov(lm(funsimP ~ 0 + node + block_size, data=atugigamat))
#~ b = boxplot(funsimP ~ 0 + node + block_size, data=atugigamat, plot=F)
#~ boxplot(funsimP ~ 0 + node + block_size, data=atugigamat, names=sapply(b$names, function(s){ strsplit(s, split='\\.')[[1]][2] }),
#~  main='Distribution of functional similarities whithin transferred block received along A. tumefaciens history',
#~  xlab="block size", ylab="functional similarity on GO aspect \"Biological Process\"", col=cols[1:2], ylim=c(0,1))
#N = max(atugigamat$block_size)
#~ abline(v=((1:N-2)*2+0.5))
#~ legend('topleft', fill=c(NA,cols[1:2]), legend=c('receiver genome', c('internal', 'leaves')), bg='white') 
#mtext(line=4, side=1, adj=1, text=paste(capture.output(summary(funsim.aov))[1:5], collapse='\n'))
#~ text(x=(N-1)*2, y=0.8, adj=1, labels=paste(capture.output(summary(funsim.aov))[1:5], collapse='\n'))
#~ 
#~ 
#~ funsim.aov = aov(lm(funsimF ~ 0 + node + block_size, data=atugigamat))
#~ b = boxplot(funsimF ~ 0 + node + block_size, data=atugigamat, plot=F)
#~ boxplot(funsimF ~ 0 + node + block_size, data=atugigamat, names=sapply(b$names, function(s){ strsplit(s, split='\\.')[[1]][2] }),
#~  main='Distribution of functional similarities whithin transferred block received along A. tumefaciens history',
#~  xlab="block size", ylab="functional similarity on GO aspect \"Molecular Function\"", col=cols[1:2], ylim=c(0,1))
#N = max(atugigamat$block_size)
#~ abline(v=((1:N-2)*2+0.5))
#~ legend('topleft', fill=c(NA,cols[1:2]), legend=c('receiver genome', c('internal', 'leaves')), bg='white') 
#mtext(line=4, side=1, adj=1, text=paste(capture.output(summary(funsim.aov))[1:5], collapse='\n'))
#~ text(x=(N-1)*2, y=0.8, adj=1, labels=paste(capture.output(summary(funsim.aov))[1:5], collapse='\n'))

###############################

multimat = lapply(currgenomes, function(cg){
 anc = atustrainspecies[[cg]]
 strainanc = c(cg, anc)
 td = distlbfh[distlbfh$current_genome==cg & (distlbfh$receiver %in% strainanc),]
 td$block_size = factor(td$block_size, levels=2:N)
 td$current_genome = factor(td$current_genome, levels=currgenomes)
 td = td[!is.na(td$block_size),]
})
gigamat = do.call(rbind, multimat)
leaves = as.character(gigamat$receiver)==as.character(gigamat$current_genome)

gigamat$node = ifelse(leaves, 'leaf', 'internal')
atugigamat = gigamat[!is.na(gigamat$node),]

#~ t.test(gigamat$funsim[atuleaves], gigamat$funsim[atuinternals])
funsim.aov = aov(lm(funsim ~ 0 + node + block_size, data=atugigamat))
b = boxplot(funsim ~ 0 + node + block_size, data=atugigamat, plot=F)
boxplot(funsim ~ 0 + node + block_size, data=atugigamat, names=sapply(b$names, function(s){ strsplit(s, split='\\.')[[1]][2] }),
 main='Distribution of functional similarities whithin transferred blocks received by A. tumefaciens strains vs. species ancestors',
 xlab="block size", ylab="Schlicker's funSim", col=cols[1:2], ylim=c(0,1))
#~ N = max(atugigamat$block_size)
abline(v=((1:N-2)*2+0.5))
legend('topleft', fill=c(NA,cols[1:2]), legend=c('receiver genome', c('genomic species or\nclosest higher-order ancestor', 'leaves')), bg='white') 
#~ mtext(line=4, side=1, adj=1, text=paste(capture.output(summary(funsim.aov))[1:5], collapse='\n'))
text(x=(N-1)*2, y=0.8, adj=1, labels=paste(capture.output(summary(funsim.aov))[1:5], collapse='\n'))

#~ funsim.aov = aov(lm(funsimP ~ 0 + node + block_size, data=atugigamat))
#~ b = boxplot(funsimP ~ 0 + node + block_size, data=atugigamat, plot=F)
#~ boxplot(funsimP ~ 0 + node + block_size, data=atugigamat, names=sapply(b$names, function(s){ strsplit(s, split='\\.')[[1]][2] }),
#~  main='Distribution of functional similarities whithin transferred blocks received by A. tumefaciens strains vs. species ancestors',
#~  xlab="block size", ylab="functional similarity on GO aspect \"Biological Process\"", col=cols[1:2], ylim=c(0,1))
#~ #N = max(atugigamat$block_size)
#~ abline(v=((1:N-2)*2+0.5))
#~ legend('topleft', fill=c(NA,cols[1:2]), legend=c('receiver genome', c('genomic species or\nclosest higher-order ancestor', 'leaves')), bg='white') 
#~ #mtext(line=4, side=1, adj=1, text=paste(capture.output(summary(funsim.aov))[1:5], collapse='\n'))
#~ text(x=(N-1)*2, y=0.8, adj=1, labels=paste(capture.output(summary(funsim.aov))[1:5], collapse='\n'))
#~ 
#~ funsim.aov = aov(lm(funsimF ~ 0 + node + block_size, data=atugigamat))
#~ b = boxplot(funsimF ~ 0 + node + block_size, data=atugigamat, plot=F)
#~ boxplot(funsimF ~ 0 + node + block_size, data=atugigamat, names=sapply(b$names, function(s){ strsplit(s, split='\\.')[[1]][2] }),
#~  main='Distribution of functional similarities whithin transferred blocks received by A. tumefaciens strains vs. species ancestors',
#~  xlab="block size", ylab="functional similarity on GO aspect \"Molecular Function\"", col=cols[1:2], ylim=c(0,1))
#~ #N = max(atugigamat$block_size)
#~ abline(v=((1:N-2)*2+0.5))
#~ legend('topleft', fill=c(NA,cols[1:2]), legend=c('receiver genome', c('genomic species or\nclosest higher-order ancestor', 'leaves')), bg='white') 
#~ #mtext(line=4, side=1, adj=1, text=paste(capture.output(summary(funsim.aov))[1:5], collapse='\n'))
#~ text(x=(N-1)*2, y=0.8, adj=1, labels=paste(capture.output(summary(funsim.aov))[1:5], collapse='\n'))

dev.off()

#~ ##############################################################
#~ # Clustering of clade-specific genes
#~ ##############################################################
#~ 
#~ cladspegenecoord = read.table('coordinates_clade-specifc_genes_relaxed2')
#~ 



