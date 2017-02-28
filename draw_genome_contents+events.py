#!/usr/bin/python
# -*- coding: utf8 -*-

import sys
import tree2
import scipy

import pysvg.structure
import pysvg.builders
import pysvg.text
import pysvg.shape

def eventangles(devent, genomesize):
	devtangle = {}
	for evt in devent:
		devtangle[evt] = (2.0 * scipy.pi * float(devent[evt])/genomesize)
	devtangle['gainreplace'] = devtangle['gain'] + devtangle['replacement']
	return devtangle

#~ nfin = sys.argv[1]
#~ nfin = '/home/lassalle/agrogenom/ancestral_content_230413/Count.AsymmetricWagner.gain_15/synthesis/translocations/genome_synthesis.replicon_location.pickle'
#~ nfin = '/home/lassalle/agrogenom/synthesis_270813/genome_synthesis.pickle'
nfin = '/home/lassalle/agrogenom/synthesis_270813/block_aware/genome_synthesis.pickle'
#~ nftabin = '/home/lassalle/agrogenom/synthesis_270813/genome_synthesis.tab'
#~ ftabin = open(nftabin, 'r')
#~ head = ftabin.readline().rstrip('\n').split('\t')
#~ dfieldsi = dict(zip(head, range(len(head))))
#~ dnodefields = {}
#~ for line in ftabin:
	#~ lsp = line.rstrip('\n').split('\t')
	#~ dnodefields[lsp[0]] = lsp
#~ ftabin.close()
magnif = 1.5
#~ magnif = sys.argv[2]
topanccode = 'N15'
#~ topanccode = sys.argv[3]

reftree = tree2.load_pickle(nfin)

doc = pysvg.structure.svg()
sb = pysvg.builders.ShapeBuilder()

topanc = reftree[topanccode]
dnodecoord = topanc.getDictNodeCoords(phylofact=25000, interleaves=100)
#~ print dnodecoord

devtcol = {'loss':'red', 'gain':'blue', 'duplication':'cyan', 'replacement':'lime', 'gainreplace':'lime'}
devtstyle = {}
devtstyle['genome'] = pysvg.builders.StyleBuilder({"stroke":'black', "fill":'white' , "stroke-width":"1"})
devtstyle['core'] = pysvg.builders.StyleBuilder({"stroke":'black', "fill":'#ccccdc' , "stroke-width":"1"})
devtstyle['specific'] = pysvg.builders.StyleBuilder({"stroke":'black', "fill":'#ffcc00' , "stroke-width":"1"})
for evy in devtcol:
	devtstyle[evy] = pysvg.builders.StyleBuilder({"stroke":'black', "fill":devtcol[evy] , "stroke-width":"1"})

for node in reftree:
	if not node.is_childorself(topanc): continue
	x, y = dnodecoord[node.label()]
	code = node.label()
	genomesize = float(node.homolog_count())
	devents = node.getEvents()
	genome = pysvg.structure.g()
	#~ name = pysvg.text.text(node.label(), 0, y)
	#~ genome.addElement(name)
	
	# circle of the size of the genome
	r = magnif * scipy.sqrt(float(genomesize)/(scipy.pi))
	genecontent = sb.createCircle(x, y, r)
	#~ if node.is_leaf():
		#~ genecontent.set_style(devtstyle['core'].getStyle())
		#~ genome.addElement(genecontent)
	#~ else:
		#~ genecontent.set_style(devtstyle['genome'].getStyle())
		#~ genome.addElement(genecontent)
		#~ # circle of the size of the core genome
		#~ rc = magnif * scipy.sqrt(float(dnodefields[code][dfieldsi['coresize']])/(scipy.pi))
		#~ coregenome = sb.createCircle(x, y, rc)
		#~ coregenome.set_style(devtstyle['core'].getStyle())
		#~ genome.addElement(coregenome)
	#~ # circle of the size of the specifc genome
	#~ rs = magnif * scipy.sqrt(float(dnodefields[code][dfieldsi['cladespecific']])/(scipy.pi))
	#~ spegenome = sb.createCircle(x, y, rs)
	#~ spegenome.set_style(devtstyle['specific'].getStyle())
	#~ genome.addElement(spegenome)
	genecontent.set_style(devtstyle['genome'].getStyle())
	genome.addElement(genecontent)
	
	devtangles = eventangles(devents, genomesize)
	pdataM = "M %f %f "%(x-r, y)
	# make a circle arcs which angle is proportional to the nb of lost genes
	xl, yl = x-15, y+15
	pdataMl = "M %f %f "%(xl-r, yl)
	angle = devtangles['loss']
	pdataA = "a %f,%f 0 0,0 %f,%f"%(r, r, r*(1-scipy.cos(angle)), r*(+scipy.sin(angle)))
	arc = pysvg.shape.path(pathData=pdataMl+pdataA)
	
	arc.appendLineToPath(xl, yl, relative=False)
	arc.appendCloseCurve()
	arc.set_style(devtstyle['loss'].getStyle())
	genome.addElement(arc)
	
	
	# make two circle arcs which angle is proportional to the nb of gained genes + replaced genes, gained genes, duplicated genes
	for eventtype in ['gainreplace', 'gain', 'duplication']:
		angle = devtangles[eventtype]
		largearcflag = (0 if angle<scipy.pi else 1) 
		pdataA = "a %f,%f 0 %d,1 %f,%f"%(r, r, largearcflag, r*(1-scipy.cos(angle)), r*(-scipy.sin(angle)))
		#~ pdataA += " z"		# close the path
		arc = pysvg.shape.path(pathData=pdataM+pdataA)
		arc.appendLineToPath(x, y, relative=False)
		arc.appendCloseCurve()
		arc.set_style(devtstyle[eventtype].getStyle())
		genome.addElement(arc)
			
	doc.addElement(genome)
	
	
doc.save('genome_contents+events.svg')
#~ print doc.getXML()
fout = open('reftree.svg', 'w')
fout.write(topanc.svgTree(dnodecoord=dnodecoord, comment='gainlosscount', textorbit=50, defaultstyle='stroke:black; fill:none; stroke-width:4; '))
fout.close()
