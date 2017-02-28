#!/usr/bin/python
# -*- coding: utf8 -*-

import sys
import tree2
import scipy

import pysvg.structure
import pysvg.builders
import pysvg.text
import pysvg.shape


def eventradius(devent, repliconsize, magnif=1.0):
	devtradius = {}
	for evt in devent:
		devtradius[evt] = magnif * scipy.sqrt(2.0 * float(devent[evt])/scipy.pi)
	return devtradius

nfin = sys.argv[1]
# example: nfin = '/home/lassalle/agrogenom/synthesis_270813/translocations/genome_synthesis.replicon_location.pickle'
# this file is the main output from 'chromosomal_translocations.py' script
# example: magnif = 1.2
magnif = sys.argv[2]
# example: topanccode = 'N15'
topanccode = sys.argv[3]

reftree = tree2.load_pickle(nfin)

doc = pysvg.structure.svg()
sb = pysvg.builders.ShapeBuilder()

topanc = reftree[topanccode]
#~ dnodecoord = topanc.getDictNodeCoords(phylofact=25000, interleaves=100)
dnodecoord = topanc.getDictNodeCoords(cladogram=True, interleaves=150, cladofact=3)
#~ print dnodecoord

dreplicoord = {'I':(0,0), 'II':(100,0), 'pAt':(-20, 50), 'pTi':(40, 50), 'p':(80, 50), '?': (120,50)}
dreplitxt = {'I':'primary', 'II':'secondary', 'pAt':'pAt', 'pTi':'pTi', 'p':'plasmid(s)', '?':'unknown'}

devtcol = {'loss':'red', 'gain':'blue', 'duplication':'cyan', 'replacement':'green', 'gainreplace':'green'}
devtstyle = {}
devtstyle['genome'] = pysvg.builders.StyleBuilder({"stroke":'black', "fill":'white' , "stroke-width":"1"})
for evy in devtcol:
	devtstyle[evy] = pysvg.builders.StyleBuilder({"stroke":devtcol[evy], "fill":devtcol[evy] , "stroke-width":"1"})

for node in reftree:
	if not node.is_childorself(topanc): continue
	xg, yg = dnodecoord[node.label()]
	#~ print node.label(), xg, yg
	code = node.label()
	drepevents = node.misc()
	genome = pysvg.structure.g()
	cell = sb.createRect(xg-40*magnif, yg-40*magnif, 200*magnif, 120*magnif, 10, 10)
	cell.set_style(devtstyle['genome'].getStyle())
	genome.addElement(cell)
	
	for repli in dreplicoord:
		xr, yr = dreplicoord[repli]
		x = xg + xr*magnif
		y = yg + yr*magnif
		#~ print repli, x, y
		replicon = pysvg.structure.g()
		repliconsize = drepevents.get(repli,0)
		#~ name = pysvg.text.text(node.label(), 0, y)
		#~ replicon.addElement(name)
		
		r = magnif * scipy.sqrt(float(repliconsize)/(scipy.pi))
		genecontent = sb.createCircle(x, y, r)
		genecontent.set_style(devtstyle['genome'].getStyle())
		replicon.addElement(genecontent)
		devtradius = eventradius(drepevents, repliconsize, magnif)

		devtup = {'gain': ('-', repli), 'loss':(repli, '-')}
		# make two circle arcs which angle is proportional to the nb of gained genes + replaced genes, gained genes, duplicated genes
		for eventtype in devtup:
			repev = devtup[eventtype]
			re = devtradius.get(repev)
			#~ print eventtype, re
			if not re: continue
			pdataM = "M %f %f "%(x, y+re)
			sweepflag = (1 if eventtype=='gain' else 0)
			pdataA = "a %f,%f 0 1,%d %f,%f"%(re, re, sweepflag, 0, 0-2*re)
			arc = pysvg.shape.path(pathData=pdataM+pdataA)
			arc.appendCloseCurve()
			arc.set_style(devtstyle[eventtype].getStyle())
			replicon.addElement(arc)
		genome.addElement(replicon)
	
	for repli in dreplicoord:	
		xr, yr = dreplicoord[repli]
		x = xg + xr*magnif
		y = yg + yr*magnif		
		for replirec in dreplicoord:
			if (repli, replirec) in drepevents:
				# translocation
				xc, yc = dreplicoord[replirec]
				xrr = xg+xc*magnif
				yrr = yg+yc*magnif
				i = (2 if repli>replirec else -2)
				#~ rrl = sb.createLine(x+i, y+i, xrr+i, yrr+i)
				pdataM = "M %d %d "%(x+i, y+i)
				pdataL = "L %d %d "%(xrr+i, yrr+i)
				rrl = pysvg.shape.path(pathData=pdataM+pdataL)
				lstyle = "stroke:#00f800; stroke-width:2; marker-end:url(#Arrow2Mend)"
				rrl.set_style(lstyle)
				replicon.addElement(rrl)
				rrcount = pysvg.text.text( "%d %s->%s"%(drepevents[(repli, replirec)], repli, replirec), ((x+xrr)/2)+4*i, ((y+yrr)/2)+4*i)
				genome.addElement(rrcount)
			
	doc.addElement(genome)
	
	
doc.save('replicon_contents+events.svg')
#~ print doc.getXML()
fout = open('reftree.svg', 'w')
fout.write(topanc.writeSvgTree(dnodecoord=dnodecoord, labels=True, textorbit=50, defaultstyle='stroke:black; fill:none; stroke-width:5; '))
fout.close()
