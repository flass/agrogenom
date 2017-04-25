#~ library(MASS)
#~ orderLoci = function(loc, repli, repliOrder){
#~     orderLoc = c()
#~     repliBound = c()
#~     start = 0
#~     for (r in repliOrder){
#~         if (start != 0){ repliBound = c(repliBound, start) }
#~         rloc = loc[repli==r] + start
#~         orderLoc = c(orderLoc, rloc)
#~         start = max(rloc)
#~     }
#~     return(list(orderLoc, repliBound))
#~ }


combineAspects = function(gosim, aspects, metric, groupmetric){
		# combine mesures for different aspects of GO as in Schlicker et al. (BMC Bioinfo. 2006)
		aspectgosim = sapply(aspects, function(a){ gosim[gosim$metric==metric & gosim$GOaspect==a, groupmetric] })
		if (length(dim(aspectgosim))==2){
			aspectgosim = as.data.frame(aspectgosim)^2
			if (any(!is.na(aspectgosim))){ return(apply(aspectgosim, 1, sum)/length(aspects))
			}else{ return(NA) }
		}else{ return(NA) }
}

plotFunctionalHomogeneity = function(blockgosim, metric, groupmetric, aspects, blocksize, code, randgosim=NULL, consgosim=NULL, strainancs=NULL){
	if (blocksize > 0){
		bcondb = blockgosim$windowsize==blocksize
		if (!is.null(randgosim)){
			# random gene groups
			rcondb = randgosim$windowsize==blocksize
			condr = (rcondb & randgosim$loc==-1)
			rgs = combineAspects(randgosim[condr,], aspects, metric, groupmetric)
			nr = length(which(!is.na(rgs)))
			# systematic gene windows
			condw = (rcondb & randgosim$loc>=0)
			wgs = combineAspects(randgosim[condw,], aspects, metric, groupmetric)
			nw = length(which(!is.na(wgs)))
		}else{
			nr = nw = 0
		}
		if (!is.null(consgosim)){ 
			# conserved tranfer block event groups
			ccondb = consgosim$windowsize==blocksize
			cgs = combineAspects(consgosim[ccondb,], aspects, metric, groupmetric)
			nc = length(which(!is.na(cgs)))
		}else{
			nc = 0
		}
	}else{
		# all blocks
		bcondb = rcondb = TRUE
		blocksize = 'all'
		nc = nr = 0
	}
	# tranfer block event groups
	bgs = combineAspects(blockgosim[bcondb,], aspects, metric, groupmetric)
	nb = length(which(!is.na(bgs)))
	lwt = NA
	means = c(NA, NA)
	if (nb >= 10){
		print(c("# co-transfered blocks", nb))
		if (is.null(consgosim)){
			if ((nw >= 10) & (nr >= 10)){
			## plot the expected functional homogeneity distributions
			# plot the similarity distribution of random (uniformly sampled) gene groups
			plot(density(rgs, bw=0.05, na.rm=T), xlim=c(0, 1), ylim=c(0, 3), col='black', lwd=2,
			 xlab=paste('functional similarity ( metrics:', metric, '/', groupmetric, ')'),
			 main=paste('funSim measure (combination of', paste(aspects, collapse=', '), 'GO aspects)\nsimilarity in', code, 'for groups of', blocksize, 'genes')) 
			# plot the similarity distribution of all windows of gene group
			lines(density(wgs, bw=0.05, na.rm=T), col='blue', lwd=2)					
			## plot the functional similarity distribution of all observed block events
			lines(density(bgs, bw=0.05, na.rm=T), col='red', lwd=2)
			legtext = c("randomly sampled genes", "random gene windows", paste("block event genes", "\n(n =", nb, ')'))
			legcol = c('black', 'blue', 'red')
			legend('topright', legend=legtext, col=legcol, lwd=2, bg='white')
			wt =  wilcox.test(bgs, wgs) #, alternative='greater')
			means = c(mean(bgs, na.rm=T), mean(wgs, na.rm=T))
			}else{
				print('too few data points')
				print(c('nw', nw, 'nr', nr))
				wt = NULL
			}
		}else{	
			if (nc >= 10){ 			
			## plot the functional similarity distribution of (non-conserved) observed block events
			plot(density(bgs, bw=0.05, na.rm=T), xlim=c(0, 1), ylim=c(0, 3), col='red', lwd=2,
			 xlab=paste('functional similarity ( metrics:', metric, '/', groupmetric, ')'),
			 main=paste('funSim measure (combination of', paste(aspects, collapse=', '), 'GO aspects)\nsimilarity in', code, 'for groups of', blocksize, 'genes')) 
			## plot the functional similarity distribution of conserved observed block events
			print(c("# conserved co-transfered blocks", nc))
			lines(density(cgs, bw=0.05, na.rm=T), col='green', lwd=2)
			legtext = c(paste("block event genes\nnot conserved in", strainancs[code], "\n(n =", nb, ')'), paste("block event genes\nconserved in", strainancs[code], "\n(n =", nc, ')'))
			legcol = c('red', 'green')
			legend('topright', legend=legtext, col=legcol, lwd=2, bg='white')
			wt =  wilcox.test(cgs, bgs) #, alternative='greater')
			means = c(mean(cgs, na.rm=T), mean(bgs, na.rm=T))
			}else{
				print('too few data points')
				print(c('nc', nc))
				wt = NULL
			}
		}
		lwt = ifelse(is.null(wt), NA, -log(wt$p.value))
	}else{
		print('too few data points')
		print(c('nb', nb))
	}
	return(c(lwt, means))
}

distRecDon = function(blockgosim, spedist){
	dRD = sapply(1:dim(blockgosim)[1], function(i){
		d = spedist[blockgosim[i, 'code.rec'], blockgosim[i, 'code.don']]
	})
	return(dRD)
}



#~ metrics = c('funSimAverage', 'funSimMax')
#~ groupmetrics = c("meanMaxSim", "meanPWSim")
metrics = c('funSimMax')
groupmetrics = c("meanMaxSim")
#~ aspects = c('biological_process', 'molecular_function') # ,'cellular_component'
coul = c('red', 'blue', 'black', 'green')
#~ blocksizes = c(0,2:10) # , 20, 30
#~ minannotperwindow = 2

# information relative to the reference phylogeny
#~ agrodatadir = '/home/lassalle/agrogenom/'
agrodatadir = '/pandata/lassalle/agrogenom/'
GOdir = paste(agrodatadir, 'duplications/maxlossrate02/blockevents_260813/GeneOntology/', sep='')
sggfdir = paste(GOdir, 'sggf_out/', sep='')
lcodes = as.character(read.table(paste(GOdir, "Agrobacterium.codes", sep=''))[,1])
reftreedir = paste(agrodatadir, 'reftree/', sep='')
#~ speinfo = data.frame(read.table(paste(reftreedir, '47RHIZOB.age_lineages.tab', sep=''), h=F))
speinfo = data.frame(read.table(paste(reftreedir, '47RHIZOB.age_nodes.tab', sep=''), h=F))
colnames(speinfo) = c('code', 'sp_node_id', 'age')
spedist = data.frame(read.table(paste(reftreedir, '47RHIZOB.distance_nodes.mat', sep=''), h=T))
rownames(spedist) = colnames(spedist)
nfsqlout = paste(sggfdir, 'funsim_dump.tab', sep='')
write('', nfsqlout, ncolumns=1, sep="", append=FALSE)
spenames = sort(speinfo$code)
maphhb = rep(0, dim(speinfo)[1])
names(maphhb) = spenames

# correspondancies between strains and upper nodes, including genomic species ancestors
tlineages = read.table('/pandata/lassalle/agrogenom/reftree/47RHIZOB.leaf_lineages.tab', sep='\t', stringsAsFactors=F)
lineages = sapply(tlineages[,1], function(cg){ c(cg ,strsplit(tlineages[tlineages[,1]==cg,2], split=' ')[[1]]) })
atustrains = sapply(tlineages[,1], function(cg){ 'N15' %in% lineages[[cg]] })
atulineages = sapply(tlineages[atustrains,1], function(cg){ c(cg ,strsplit(strsplit(tlineages[tlineages[,1]==cg,2], split=' N15 ')[[1]][1], split=' ')[[1]]) })
atugenomovars = c('N35', 'N38', 'N39', 'N42', 'N32', 'N27', 'N29', 'N37')
currgenomes = names(atulineages)
atustrainspecies = sapply(currgenomes, function(cg){ max(atulineages[[cg]][atulineages[[cg]] %in% atugenomovars]) })

# Gene Ontology aspects
aspects = c('F', 'C', 'P')
names(aspects) = c('molecular_function', 'cellular_compartment', 'biological_process')
aspects = aspects[c('molecular_function', 'biological_process')]

bigs = NULL
layout(matrix(1:2, 1,2))
MWs = lapply(lcodes, function(code){
	nfrandgosim = paste(sggfdir, code, '.genegroup_funsim', sep='')
	nfbgosim = paste(sggfdir, code, '.eventblocks_funsim', sep='')
	nfblockinfos = paste(sggfdir, code, '.blockinfos', sep='')
	if (!is.na(atustrainspecies[code])){
		nfcgosim = paste(sggfdir, code, '.', atustrainspecies[code], '-conserved.eventblocks_funsim', sep='')
		nfconsinfos = paste(sggfdir, code, '.blockinfos', sep='')
		print(nfcgosim)
	}else{
		nfcgosim = nfconsinfos = NA	
	}
	if (file.exists(nfrandgosim) & file.exists(nfbgosim) & file.exists(nfblockinfos)){
		randgosim = read.table(nfrandgosim, h=T, sep='\t')
		randgosim$GOaspect = sapply(randgosim$GOaspect, function(a){ aspects[a] })
		bgosim = read.table(nfbgosim, h=T, sep='\t')
		blockinfos = read.table(nfblockinfos, h=T, sep='\t')
		if (dim(randgosim)[1]>0 & dim(bgosim)[1]>0 & dim(blockinfos)[1]>0){
			# collect data relative to transfer coordinates in species tree
			pdf(paste(sggfdir, code,'.distrib_funcsim.pdf', sep=''), width=5, height=5, onefile=T)
			blockgosim = merge(x=bgosim, y=blockinfos, by.x='leafblockid', by.y='leaf_block_id', all.x=TRUE)
			blockgosim = merge(x=blockgosim, y=speinfo, by.x='rec_sp_node_id', by.y='sp_node_id', all.x=TRUE)
			blockgosim = merge(x=blockgosim, y=speinfo, by.x='don_sp_node_id', by.y='sp_node_id', suffixes=c(".rec", ".don"), all.x=TRUE)
			blockgosim$GOaspect = sapply(blockgosim$GOaspect, function(a){ aspects[a] })
			print(c('blockgosim', code))
			print(table(blockgosim$event_type))
			blockgosim$dist_rec_don = distRecDon(blockgosim, spedist)
			if (is.null(bigs)){ bigs = blockgosim 
			}else{ bigs = rbind(bigs, blockgosim) }
			# insure independency of samples
			notransferinrandom = !(randgosim$labels %in% blockgosim$labels)
			randgosim = randgosim[notransferinrandom,]
			
			if (!is.na(nfcgosim)){
				if (file.exists(nfcgosim) & file.exists(nfconsinfos)){
					cgosim = read.table(nfcgosim, h=T, sep='\t')
					consinfos = read.table(nfconsinfos, h=T, sep='\t')
					consgosim = merge(x=cgosim, y=consinfos, by.x='leafblockid', by.y='leaf_block_id', all.x=TRUE)
					consgosim$GOaspect = sapply(consgosim$GOaspect, function(a){ aspects[a] })
					print(c('consgosim', code))
					print(table(consgosim$event_type))
					# insure independency of samples
					noconservervedintransfer = !(blockgosim$leafblockid %in% consgosim$leafblockid)
					noconsgosim = blockgosim[noconservervedintransfer,]
				}else{
					consgosim = NULL
				}
			}else{
				consgosim = NULL
			}
			
			mannwitneys = lapply(groupmetrics, function(groupmetric){
				bgos = bgosim[,c('leafblockid', 'GOaspect', groupmetric, 'metric', 'windowsize')]
				bgos$groupmetric = rep(groupmetric, dim(bgos)[1])
				write.table(bgos, nfsqlout, sep="\t", append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
		#~ 			for (aspect in aspects){
		#~ 				bgs = blockgosim[blockgosim$GOaspect==aspect & blockgosim$metric==metric & blockgosim$age.rec>0,]			
		#~ 				image(kde2d(bgs[['age.rec']], bgs[[groupmetric]]), xlab="receptor clade age",
		#~ 				 ylab=paste('block functional similarity ( metrics:', metric, '/', groupmetric, ')'), main=paste(aspect, 'similarity in', code))
		#~ 				image(kde2d(bgs[['dist_rec_don']], bgs[[groupmetric]]), xlab="receptor-to-donor clade distance",
		#~ 				 ylab=paste('block functional similarity ( metrics:', metric, '/', groupmetric, ')'), main=paste(aspect, 'similarity in', code))
		#~ 			}
				mw = lapply(metrics, function(metric){
					mwbw = as.data.frame(t(sapply(as.numeric(levels(as.factor(bgos$windowsize))), function(blocksize){
						print(paste(metric, groupmetric, blocksize, 'BvsW', code))
						mannwitney = plotFunctionalHomogeneity(blockgosim[blockgosim$event_type!='D',], metric, groupmetric, aspects, blocksize, code, randgosim=randgosim)
					return(c('BvsW', code, blocksize, mannwitney))
					})))
					if (!is.na(nfcgosim)){
						mwcnc = as.data.frame(t(sapply(as.numeric(levels(as.factor(bgos$windowsize))), function(blocksize){
							print(paste(metric, groupmetric, blocksize, 'CvsNC', code))
							mannwitney = plotFunctionalHomogeneity(noconsgosim[noconsgosim$event_type!='D',], metric, groupmetric, aspects, blocksize, code, consgosim=consgosim, strainancs=atustrainspecies)
						return(c('CvsNC', code, blocksize, mannwitney))
						})))
					}else{
						mwcnc = data.frame(comp=c(), code=c(), blocksize=c(), U=numeric(), testmean=numeric(), controlmean=numeric(), stringsAsFactors=F)
					}
					return(rbind(mwbw, mwcnc))
					
				})
				names(mw) = metrics
				return(mw)
			})
			names(mannwitneys) = groupmetrics
			dev.off()
			hhb = blockgosim[(blockgosim$windowsize > 2 & (blockgosim$meanMaxSim >=0.75 | blockgosim$meanPWSim >=0.75)),]
#~ 			write.table(hhb, paste(sggfdir, code, ".highly_homogeneous_blocks", sep=''), sep="\t", row.names=F)
			maphhb = maphhb + table(blockgosim$code.rec)[spenames]
		}else{ mannwitneys = data.frame(comp=c(), code=c(), blocksize=c(), U=numeric(), testmean=numeric(), controlmean=numeric(), stringsAsFactors=F) }
	}else{ mannwitneys = data.frame(comp=c(), code=c(), blocksize=c(), U=numeric(), testmean=numeric(), controlmean=numeric(), stringsAsFactors=F) }
	return(mannwitneys)
})
names(MWs) = lcodes
save(MWs, file=paste(sggfdir, "MWW_funsim.RData", sep=""))

nfmaphhb = paste(sggfdir, "map_to_nodes.highly_homogeneous_blocks", sep='')
write.table(maphhb, nfmaphhb, append=FALSE, sep="\t")

pdf(paste(sggfdir, "MWWdiff_distr_funsim.pdf", sep=''), height=5, width=5)
for (groupmetric in groupmetrics){
 	for (metric in metrics){
		blockU = data.frame(comp=c(), code=c(), blocksize=c(), U=numeric(), testmean=numeric(), controlmean=numeric(), stringsAsFactors=F) #, coul=c())
		for (code in lcodes){
			mws = MWs[[code]][[groupmetric]][[metric]]
			blockU = rbind(blockU, mws) #cbind(mws, rep(code, dim(mws)[1])))
		}
		colnames(blockU) = c('comp', 'code', 'blocksize','U', 'testmean', 'controlmean')
		blockU = blockU[!is.na(blockU$U),]
		blockU$blocksize = as.numeric(as.character(blockU$blocksize))
		blockU$U = as.numeric(as.character(blockU$U))
		blockU$testmean = as.numeric(as.character(blockU$testmean))
		blockU$controlmean = as.numeric(as.character(blockU$controlmean))
		print(blockU)
		
		diffmeans =  blockU$testmean - blockU$controlmean
		blockU$blocksize = blockU$blocksize + ((diffmeans>0)*0.1)
		blockU$blocksize = blockU$blocksize - ((diffmeans<=0)*0.1)
		plot(blockU$blocksize[blockU$comp=='BvsW'], blockU$U[blockU$comp=='BvsW'], xlim=c(2,6), ylim=c(0,25), pch=1, col=ifelse(diffmeans[blockU$comp=='BvsW']>0, 'red', 'blue'),
		 ylab="-log10(Mann-Witney U test p-value)", xlab="size of co-transfered blocks",
		 main="difference of functional homogeneity\nbetween random windows and co-transfered gene blocks")
		legend("topright", pch=1, col=c('blue', 'red'),
		 legend=c(paste("FH(random windows) > FH(transferred blocks) ( n =", length(which(diffmeans[blockU$comp=='BvsW']<0)), ")"),
		  paste("FH(random windows) < FH(transferred blocks) ( n =", length(which(diffmeans[blockU$comp=='BvsW']>0)), ")")))
		
		plot(blockU$blocksize[blockU$comp=='CvsNC'], blockU$U[blockU$comp=='CvsNC'], xlim=c(2,6), ylim=c(0,8), pch=1, col=ifelse(diffmeans[blockU$comp=='CvsNC']>0, 'green', 'purple'),
		 ylab="-log10(Mann-Witney U test p-value)", xlab="size of co-transfered blocks",
		 main="difference of functional homogeneity\nbetween co-transfered gene blocks\nconserved or not in genomic species")
		legend("topright", pch=1, col=c('purple', 'green'),
		 legend=c(paste("FH(non-conserved blocks) > FH(conserved blocks) ( n =", length(which(diffmeans[blockU$comp=='CvsNC']<0)), ")"),
		  paste("FH(non-conserved blocks) < FH(conserved blocks) ( n =", length(which(diffmeans[blockU$comp=='CvsNC']>0)), ")")))
	}
}
dev.off()
 	
save(bigs, file=paste(sggfdir, "score_genegroup_funsim.RData", sep=""))

#~ pdf(paste(sggfdir, 'leafVSold.distrib_funcsim.pdf', sep=''), width=10, height=10, onefile=T)
#~ # global interspecies plot
#~ for (groupmetric in groupmetrics){
#~ 	for (metric in metrics){
#~ 		for (aspect in aspects){
#~ 			big = bigs[bigs$GOaspect==aspect & bigs$metric==metric,]
#~ 			print(summary(big))
#~ 			for (blocksize in levels(as.factor(big$windowsize))){
#~ 				if (length(which(big$windowsize==blocksize)) > 25){
#~ 					bgs = big[big$windowsize==blocksize,]
#~ 					plot(density(bgs[[groupmetric]][bgs$age.rec==0], bw=0.05, na.rm=TRUE), xlim=c(0, 1), ylim=c(0, 3), col='purple', lwd=2,
#~ 					 xlab=paste('functional similarity ( metrics:', metric, '/', groupmetric, ')'),
#~ 					 main=paste(aspect, 'similarity in every species\nfor groups of', blocksize, 'genes')) 
#~ 					lines(density(bgs[[groupmetric]][bgs$age.rec>0], bw=0.05, na.rm=TRUE), col='green', lwd=2)	
#~ 					legend('topright', legend=c("recent block transfered genes", paste("(n =", length(which(bgs$age.rec==0)), ')'), "ancient block transfered genes", paste("(n =", length(which(bgs$age.rec>0)), ')')), 
#~ 					 col=c('purple', 'white', 'green', 'white'), lwd=2, bg='white')	
#~ 				}
#~ 			}
#~ 		}
#~ 	}
#~ }
#~ dev.off()
