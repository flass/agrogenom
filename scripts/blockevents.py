#!/usr/bin/python
# -*- coding: utf8 -*-


import getpass
import numpy as np
import tree2
import sys, getopt, os,  shutil
import copy, re, random
import rec_to_db


### utilitary functions
#~ def which(L, value):
	#~ w = []
	#~ i = -1
	#~ try:
		#~ while 1:
			#~ i = L.index(value, i+1)
			#~ w.append(i)
	#~ except ValueError:
		#~ pass
	#~ return w
	#~ 
#~ def whichmax(L):
	#~ return which(L, max(L))
			
### object classes for protein, block, replicon, ancestral block
class ProtInBlock(object):
	def __init__(self, dfields, eventnode, dicevent=None):
		self.__id = dfields['hogenom_gene_id']
		if not dicevent:
			dfields.update(eventnode.getdicevent())	# join informations about the evolutionary event above it
		else:
			dfields.update(dicevent)
		self.__dfields = dict(dfields)	# informations about the protein
		self.__recset = set(dfields['eventlocation'][1])			# set of possible receptors
		if len(dfields['eventlocation']) > 2:
			self.__donset = set(dfields['eventlocation'][3])			# set of possible donors
		else:
			self.__donset = set()
		self.__gteventnode = eventnode		# GeneTree node object where is located the considered event
		self.__status = None			# protein status in block : 'seed' of, 'compatible_set' or 'compatible_nni' with, 'against' or 'not_against' the common transfer scenario of the block
		self.__notagainstsets = []		# list of receptor/donor set pairs for which the tree containing the protein have no signal against a transfer event
		if self['genomic_beg'] > self['genomic_end']:
			self.editCoords(self['genomic_end'], self['genomic_beg'])
		
	def __getitem__(self, field):
		return self.__dfields[field]
		
	def __str__(self):
		return "%s (%d): %s"%(self.__id, list(self.eventids())[0], '\t'.join(self.dictToList(['locus_tag', 'family_accession', 'genomic_beg', 'genomic_end'])))
			
	def getid(self):
		return self.__id
		
	def setid(self, id):
		self.__id = id	
				
	def fields(self):
		return self.__dfields
		
	def editCoords(self, beg, end):
		self.__dfields['genomic_beg'] = min(beg, end)
		self.__dfields['genomic_end'] = max(beg, end)
	
	def modifRecSet(self, recset):
		self.__recset = set(recset)
		
	def modifDonSet(self, donset):
		self.__donset = set(donset)
		
	def getRecSet(self):
		return self.__recset
		
	def getDonSet(self):
		return self.__donset
		
	def eventnode(self):
		return self.__gteventnode
		
	def family(self):
		return self.__dfields.get('family_accession')
		
	#~ def eventfamnodes(self):
		#~ return set([(self.family(), self.eventnode().nodeid())])
		
	def eventids(self):
		return set([self['event_id']])
		
	def eventtype(self):
		return self.__dfields.get('eventtype')
		
	def getevent(self):
		return (self.__dfields.get('eventtype'), self.__dfields.get('eventlocation'))
		
	def duplication(self):
		if self.__dfields.get('eventtype')=='duplication':
			return self.__dfields.get('eventlocation')
		
	def speciation(self):
		if self.__dfields.get('eventtype')=='speciation':
			return self.__dfields.get('eventlocation')
		
	def gain(self):
		if self.__dfields.get('eventtype')=='gain':
			return self.__dfields.get('eventlocation')
			
	def transfer(self):
		if self.__dfields.get('eventtype')=='transfer':
			return self.__dfields.get('eventlocation')
			
	def getStatus(self):
		return self.__status
		
	def getEventSide(self):
		"""recognition of the side of the event"""
		side = None
		if self.eventtype() == 'transfer':
			trchildid = self['eventlocation'][4]
			genetree = self.__gteventnode.go_root()
			trchild = genetree.idgetnode(trchildid)
			if trchild:
				if genetree[self.getid()]==trchild or genetree[self.getid()].is_child(trchild):
					side = 'rec'
				else:
					side = 'don'
		return side
		
	def setStatus(self, status):
		known = ['seed', 'compatible_set', 'against', 'not_against', 'gap']	#, 'compatible_nni'
		if status not in known:
			raise ValueError, 'Unknown status : %s ; must be part of %s'%(status, str(known))
		self.__status = status
			
	def addNotAgainstSet(self, sets):
		self.__notagainstsets.append(sets)
		
	def getNotAgainstSets(self):
		return self.__notagainstsets
		
	def mergeNotAgainstSets(self):
		rs = set()
		ds = set()
		for notagainstset in self.__notagainstsets:
			rs |= set(notagainstset[0])
			ds |= set(notagainstset[1])
		self.__recset = rs
		self.__donset = ds
			
	def isSimilarTo(self, prot, fields=None):
		"""verifies the identity of both protein objects by comparing the provided list of fields (every field if none)."""
		if not fields:
			fields=self.__dfields.keys()
		for field in fields:
			#~ # ignore field for comparison if it is NULL (None in Python)
			#~ if prot[field]:
			if self[field] != prot[field]:
				return False
		else:
			return True

	#~ def possibleReceptorSet(self, reftree, receptorPair=None):
		#~ """find the node set of possible recptors in the species tree
		#~ 
		 #~ = the path between the two given cordinates
		 #~ """
		#~ if receptorPair:
			#~ receptor = receptorPair[0]
			#~ cautiousreceptor = receptorPair[1]
		#~ else:
			#~ coords = self['eventlocation']	
			#~ receptor = coords[0]
			#~ if len(coords)>1: cautiousreceptor = coords[1]
			#~ else: cautiousreceptor = coords[0]
		#~ if cautiousreceptor==receptor:
			#~ recSet = set([receptor])
		#~ else:
			#~ recSet = set( reftree[receptor].path_to(reftree[cautiousreceptor], returnLabels=True) )
		#~ self.__recset = recSet
		#~ 
	#~ def possibleDonorSet(self, reftree, donorflex=0, donorPair=None, genetree=None, minbs=0.9):
		#~ """find the node set of possible donors in the species tree
		#~ 
		 #~ = the difference between sets of nodes located under the second (cautiousdonor) and the first (donor) cordinate.
		 #~ if a gene tree is provided, will add nodes corresponding to any node reachable from donor without crossong a supported branch (support >= minbs).
		 #~ """
		#~ if donorPair:
			#~ donor = donorPair[0]
			#~ cautiousdonor = donorPair[1]
		#~ else:
			#~ coords = self['eventlocation']	
			#~ if len(coords)> 2:
				#~ donor = coords[2]
				#~ cautiousreceptor = coords[3]	
			#~ else:
				#~ self.__donset = set()
				#~ return
#~ 
		#~ if cautiousdonor==donor:
			#~ donSet = set([donor])
		#~ else:
			#~ largeset = set( reftree[cautiousdonor].get_children_labels() )
			#~ smallset = set( reftree[donor].get_children_labels() )
			#~ donSet = (largeset - smallset) | set([donor])
		#~ if genetree:
			#~ recid = coords[-1] # childid in transfer description
			#~ reccla = genetree.idgetnode(recid)
			#~ doncla = reccla.go_brother()
			#~ lposdon = doncla.explore_nonsupported_paths(minbs)
			#~ for node in lposdon:
				#~ # map possible nodes of the gene tree to the species tree
				#~ nlspe = node.dictSpeciesToLeafLabels().keys()
				#~ nref = reftree.map_to_node(nlspe)
				#~ # add those nodes
				#~ donSet = donSet | set([nref.label()])
				
		#~ ## extension of donor set to nodes at 'donor.flexibility' distance
		#~ if donorflex:
			#~ for i in range(donorflex):
				#~ for node in donSet:
					#~ f = reftree[node].go_father()
					#~ if f:
						#~ donSet = donSet | set([f.label()])
					#~ lc = reftree[node].get_children()
					#~ if lc:
						#~ for c in lc:
							#~ donSet = donSet | set([c.label()])
		#~ self.__donset = donSet

	def dictToList(self, keylist=None):
		"""converts components of a sequence into strings and return them into a list"""
		l = []
		if not keylist:
			keylist = self.__dfields.keys()
		for key in keylist:
			l.append(str(self.__dfields[key]))
		return l


class BlockEvent(object):
	def __init__(self, nblock, seed):
		self._id = nblock
		self._eventtype = seed.eventtype()
		self._seedid = seed.getid()
		self._recset = seed.getRecSet()			# set of possible locations in species tree common to block proteins (for all kind of events)
		self._donset = seed.getDonSet()			# set of possible donors in species tree common to block proteins (for transfers only)
		self._macrec = None
		self._macdon = None
		#~ self._eventfamnodes = seed.eventfamnodes()	# set of events in the block located at a node in a gene family tree
		self._eventid = seed.eventids()
				
	def getid(self):
		return self._id
				
	def edit_id(self, newid):
		self._id = newid
		
	def seed(self):
		return self[self._seedid]
		
	def eventtype(self):
		return self._eventtype
	
	def eventids(self, compute=False):
		if compute:
			self._eventid = set()
			for prot in self:
				if isinstance(prot, ProtInBlock):
					if not (prot.getStatus() in ['gap', 'against', 'not_against']):
						self.add_eventids(prot.eventids())
				else:
					self.add_eventids(prot.eventids(compute=compute))
		return self._eventid
		
	def add_eventids(self, eventid):
		if isinstance(eventid, int):
			self._eventid |= set([eventid])
		elif isinstance(eventid, set):
			self._eventid |= eventid
		else:
			raise TypeError, type(eventid)
		
	def eventidgetmembers(self, eventid):
		lp = []
		for prot in self:
			if eventid in prot.eventids():
				lp.append(prot)
		return lp
		
	def sharedeventids(self, block):
		sharedids = self.eventids() & block.eventids()
		return sharedids
				
	def getRecSet(self):
		return self._recset
		
	def getDonSet(self):
		return self._donset
			
	def commonRecDonSet(self, refine=False):
		"""defines common receptor and donor sets of proteins/blocks contained in the the block(self).
		
		By default, uses only information from proteins with transfer signal ; if refine=True, consider also the information from 'not_against' gap proteins. 
		"""
		recset = copy.copy(self.seed().getRecSet())
		donset = copy.copy(self.seed().getDonSet())
		# iters over object yielded by TransferedBlock subclass generator() function ; can be ProtInBlock or LeafBlock objects
		for prot in self:
			if isinstance(prot, ProtInBlock): # and (not prot.eventtype()=='transfer' ):
				if prot.getStatus() in ['gap', 'against']:
					# never consider undetermined 'gap' proteins or 'against' proteins (to be deleted)
					continue
				if (prot.getStatus()=='not_against') and (not refine):
					# do not consider 'not_against' gap proteins but for refining the block event coordinates
					continue
			recset &= prot.getRecSet()
			donset &= prot.getDonSet()		
		self._recset = recset
		self._donset = donset
				
	def getCommonRecDonSet(self, compute=True, refine=False):
		if compute:
			self.commonRecDonSet(refine=refine)			
		return (self._recset, self._donset)
		
	def enlargeCommonRecDonSet(self, suprec=None, supdon=None):
		"""adds nodes to receptor or donor sets (set union)"""
		if suprec:
			self._recset |= suprec
		if supdon:
			self._donset |= supdon
		
	def getMAC(self, reftree, coordtype, compute=True, returnLabel=False):
		"""returns the most ancient common receptor/donor of all the elements in the block"""
		if coordtype=='rec':
			selfmac = self._macrec
			selfset = self._recset
		elif coordtype=='don':
			selfmac = self._macdon
			selfset = self._donset
		else:
			raise ValueError, "wrong coordinate type value: %s"%coordtype
		if (not compute) and (selfmac):
			if returnLabel:
				return selfmac
			else:
				return reftree[selfmac]
		else:
			if (not selfset):
				mac = None
				maclab = None
			else:	
				mac = reftree.coalesce(selfset)
				maclab = mac.label()
				if (not mac.label() in selfset):
					if isinstance(self, LeafBlock):
						for prot in self:
							print "\t%s %s [%s / %s]"%(prot, prot.getStatus().upper(), " ".join(list(prot.getRecSet())), " ".join(list(prot.getDonSet())) )
						raise IndexError, "all %sor nodes are not in the same lineage : [%s]"%(coordtype, " ".join(list(self._recset)))
					elif isinstance(self, AncestralBlock):
						for block in self:
							for prot in block:
								print "\t%s %s [%s / %s]"%(prot, prot.getStatus().upper(), " ".join(list(prot.getRecSet())), " ".join(list(prot.getDonSet())) )
						raise IndexError, "all %sor nodes are not in the same lineage : [%s]"%(coordtype, " ".join(list(self._recset)))

			if coordtype=='rec': self._macrec = maclab
			elif coordtype=='don': self._macdon = maclab
			if returnLabel:
				return maclab
			else:
				return mac
		
	def getMACReceptor(self, reftree, compute=True, returnLabel=False):
		"""returns the most ancient common receptor of all the elements in the block"""
		return self.getMAC(reftree, 'rec', compute=compute, returnLabel=returnLabel)
			
	def getMACDonor(self, reftree, compute=True, returnLabel=False):
		"""returns the most ancient common donor of all blocks"""
		return self.getMAC(reftree, 'don', compute=compute, returnLabel=returnLabel)
		
	def compatibleRecDonSet(self, inrecset, indonset, returnRecDon=False):
		"""if Receptor and Donor sets of gene block (self) and input sets show non-null intersection, returns those intersections or boolean stating the compatibility"""
		commonRec = self._recset & inrecset
		if commonRec:
			if self._donset:
				# case of a transfer event, a second coordinate is needed
				commonDon = self._donset & indonset
				if commonDon:
					if not returnRecDon: return 1
					else: return (commonRec, commonDon)
			else:
				if not returnRecDon: return 1
				else: return (commonRec, self._donset)	
		if not returnRecDon: return 0
		else: return None
			
	def compatibleEvent(self, ngbr, returnRecDon=False, excludeTypes=True): 
		"""if Receptor and Donor sets of gene block (self) and a neighbour protein or another block show non-null intersection, returns those intersections or boolean stating the compatibility"""
		if self._eventtype == ngbr.eventtype():
			return  self.compatibleRecDonSet(ngbr.getRecSet(), ngbr.getDonSet(), returnRecDon=returnRecDon)
		elif excludeTypes:
			return 0
		else:
			return -1
			
	def getXMLTreeCollection(self, dirxmltrees):
		lnfxmltrees = []
		for fam in self.getFamList():
			nfxmltree = "%s/%s.xml"%(dirxmltrees, fam)
			if not os.access(nfxmltree, os.F_OK):
				raise ValueError, "missing phyloXML tree file at %s"%nfxmltree
			lnfxmltrees.append(nfxmltree)
		return lnfxmltrees
		
class LeafBlock(BlockEvent):
	def __init__(self, nblock, seed, side=None):
		super(LeafBlock, self).__init__(nblock, seed)
		self.__lapid = []				# list of identifiers of protein in the block
		self.__dprot = {}				# dictionaries of ProtInBlock objects
		if not side: self.__eventside = seed.getEventSide()
		else: self.__eventside = side			# side of the event : None (default) for non-directional events like duplication and speciation, or 'rec' or 'don' for transfers
		self.addProt(seed)
		
	def __getitem__(self, apid):
		return self.__dprot[apid]
		
	def __str__(self):
		return "[leaf block %s] %s [%s -> %s] :\t%s"%(str(self._id), self._eventtype, ' '.join(list(self._donset)), ' '.join(list(self._recset)), ' '.join(self.__lapid))
		
	def __iter__(self):
		return self.generator()
		
	def generator(self):
		for protid in self.__lapid:
			yield self.__dprot[protid]
			
	def eventside(self):
		return self.__eventside
		
	def addProt(self, protinblock, silent=True):
		"""adds a protein to the block"""
		if not isinstance(protinblock, ProtInBlock):
			raise TypeError, 'expected a ProtInBlock object'
		if not protinblock.getid() in self.__lapid:
			self.__lapid.append(protinblock.getid())
			self.__dprot[protinblock.getid()] = protinblock
		else:
			# can happen for duplications in series at the same location in the species tree (yield several similar events that may be compatible)
			if not self.seed().eventtype() == 'duplication':
				raise IndexError, 'Protein %s already recorded in block %s'%(protinblock.getid(), self )
		if not silent: print '\t\tadd', protinblock.getid(), protinblock.getStatus()
		if not protinblock.getStatus() in ['gap', 'not_against', 'against']:
			self.add_eventids(protinblock.eventids())
			if not silent: print self.getCommonRecDonSet()
			self.commonRecDonSet()
			if not silent: print '->', self.getCommonRecDonSet()
			
	def delProt(self, protid, pop=False, toEmpty=False):
		if protid==self._seedid:
			for apid in self.__lapid:
				if apid!=self._seedid:
					if self[apid].getStatus()=='compatible_set':
						print 'Change seed of block %s from %s to %s'%(self.getid(), self._seedid, apid)
						self._seedid = apid
						self[apid].setStatus('seed')
						break
			else:
				for apid in self.__lapid:
					if apid!=self._seedid:
						print 'Change seed of block %s from %s to %s'%(self.getid(), self._seedid, apid)
						self._seedid = apid
						self[apid].setStatus('seed')
						break
				else:
					if toEmpty:
						self._seedid = None
						self._recset = set([])
						self._donset = set([])
					else:
						raise IndexError, 'Cannot delete seed protein as it is the only one in the block'
		self.__lapid.remove(protid)
		if self.seed(): self.commonRecDonSet()
		self.eventids(compute=True)
		if not pop:
			del self.__dprot[protid]
		else:
			prot = self.__dprot[protid]
			del self.__dprot[protid]
			return prot
		
	def getProtDic(self):
		return self.__dprot
		
	def getProtList(self):
		return self.__lapid
		
	def getProts(self):
		lp = []
		for p in self.__dprot:
			lp.append(self[p])
		return lp
		
	def getDictProtToEventids(self):
		dprotevents = {}
		for prot in self:
			if not (prot.getStatus() in ['gap', 'against', 'not_against']):
				dprotevents[prot.getid()] = list(prot.eventids())
		return dprotevents
		
	def getDictProtToRecDonSets(self):
		dprotrds = {}
		for prot in self:
			recset = copy.copy(prot.getRecSet())
			donset = copy.copy(prot.getDonSet())
			#~ recset = copy.copy(prot['eventlocation'][1])
			#~ if len(prot['eventlocation'])>2: donset = copy.copy(prot['eventlocation'][3])
			#~ else: donset = set()
			dprotrds[prot] = (recset, donset)
		return dprotrds	
		
	def sortProtsByRecDonSet(self):
		dselfprotrds = self.getDictProtToRecDonSets()
		lpatterns = []
		lsamepatprots = []
		for prot in self:
			pat = dselfprotrds[prot]
			if not pat in lpatterns:
				lpatterns.append(pat)
				lsamepatprots.append([prot])
			else:
				# this pattern was already met with another prot, alredy listed
				i = lpatterns.index(pat)
				lsamepatprots[i].append(prot)
		return lpatterns, lsamepatprots				
		
	def sortByPos(self, referenceFirst=1):
		def cmp(x,y):
			# sorts by chromosomal 'genomic_beg' coordinate
			xbeg = self[x]['genomic_beg']
			ybeg = self[y]['genomic_beg']
			if xbeg<ybeg: return -1
			elif xbeg>ybeg: return 1
			else: return 0
		self.__lapid.sort(cmp)
		
			
	def eventInternality(self, eventids):
		"""return a score reflecting the degree to which events from input list are located internally or externally in the block
		
		example: if b a leaf block with events A, B, C, D arranged (A, B, C, D), 
		eventInternality(b, [B, C]) = sum_over_events(min(index(event), end(b) - index(event))) = 2
		eventInternality(b, [A, D]) = 0
		"""
		self.sortByPos()
		levt = [prot['event_id'] for prot in self]
		end = len(levt) - 1
		I = 0
		for evt in eventids:
			if evt in levt:
				i = levt.index(evt)
				I += min(i, end - i)
		return I
		
	def getBlockSize(self):
		"""Return the number of distinct protein in the block."""
		n = 0
		for prot in self:
			n+=1
		return n
		
	def getCoords(self):
		lbeg = []
		lend = []
		for prot in self:
			lbeg.append(int(prot['genomic_beg']))
			lend.append(int(prot['genomic_end']))
		return (min(lbeg), max(lend))
		
	def getFamList(self, filterStatus=[]):
		"""return the list of gene families involved in the block"""
		l = []	
		for prot in self:
			if prot.getStatus() not in filterStatus:
				l.append(prot.family())
		return list(set(l))
		
	def searchNgbrEvent(self, tprot, n, walkstep, fields, maxgapsize, reftree, cursor, dirtreepickles, nonshuffledtprot=None, silent=True):
		"""walks on replicon to find similar events by scanning all events linked to a protein and then moving to the next one
		
		only the more recent event of each type are considered, as successive events erase gene order signal.
		walkstep can be any non-null integer.
		"""
		if int(walkstep)==0:
			raise ValueError, "'walkstep' must be a non-null integer"
		if nonshuffledtprot:
			# get protein attribute from the ordered replicon
			dref = rec_to_db.getGeneInfo(nonshuffledtprot[n][0], fields, cursor)
		# moves one protein forward on the replicon
		nunknown = 0
		n += walkstep
		while (nunknown <= maxgapsize[self.eventtype()]) and n < len(tprot) :
			# gets information about the next protein
			apid = tprot[n][0]
			dgene = rec_to_db.getGeneInfo(apid, fields, cursor)
			apfam = dgene['family_accession']
			nfgenetreepickle = "%s/%s.pickle"%(dirtreepickles, apfam)
			if not os.access(nfgenetreepickle, os.F_OK):
				# may be a small (n<4) family with no tree computed
				tdfamily = rec_to_db.getFamilyInfo(apfam, '*', cursor, returnDict=False)
				if (len(tdfamily) < 4) :
					n += walkstep
					nunknown += 1
					if not silent: print "\t\t\tskip prot from small family", n, apid, apfam
					continue
				else:
					raise IndexError, "family %s has no tree computed according to database record:\n%s"%(apfam, str(tdfamily))
			if not silent: print "\t\tcheck prot", n, apid, apfam
			genetree = tree2.load_pickle(nfgenetreepickle)
			if nonshuffledtprot:
				# get protein attribute from the ordered replicon
				dref = rec_to_db.getGeneInfo(nonshuffledtprot[n][0], fields, cursor)
			# protein locus starts with unknown status relative to the event
			compat = -1
			#~ ddevents = genetree.getEvents(lineage=apid, returnDict=True)
			lineage = genetree[apid].lineage(value='id')
			recgiid = rec_to_db.getGeneTreeInfo(apfam, ['rec_gi_id'], cursor, returnDict=False)[0]
			for startnodeid in lineage:
				eventnode = genetree.idgetnode(startnodeid)
				tdevents = rec_to_db.getEventInfo(recgiid, startnodeid, ['event_id', 'event_type', 'rec_gi_end_node'], cursor, returnDict=True)
				for devent in tdevents:
				# for eventnodelab in ddevents:
					# eventnode = genetree[eventnodelab]
					# devent = ddevents[eventnodelab]
					eventtype = devent['eventtype']
					eventid = devent['event_id']
					devent['eventlocation'] = rec_to_db.getEventLocation(eventid, eventtype, cursor)
					if not silent: print "\t\t\tdevent:", devent
					ngbr = ProtInBlock( dfields=dgene, eventnode=eventnode, dicevent=devent)
					if nonshuffledtprot:
						# changes protein coordinates of the shuffled replicon into those of the ordered replicon (for compatibility with position dependent functions)
						ngbr.editCoords(dref['genomic_beg'], dref['genomic_end'])
					if ngbr.eventtype()==eventtype:
						# comparison of infered transfer patterns between block and neighbor gene
						compat = self.compatibleEvent(ngbr, excludeTypes=True)
						if compat:
							if not silent: print '\t\t-> compatible event'
							# the more recent event of this type is compatible 
							ngbr.setStatus('compatible_set')
							self.addProt(ngbr, silent=silent)
							break # the for devent loop
				if ngbr.getStatus()=='compatible_set':
					break # the for startnodeid loop
			else: 
				# if completed evaluation of all events (for devent / startnodeid loops not broken)
				if not silent: print '\t\t-> no compatible event: gap protein'
				ngbr.setStatus('gap')
				self.addProt(ngbr, silent=silent)
				nunknown += 1
				#  goes to next protein on replicon
			n += walkstep
					
					#~ compat = self.compatibleEvent(ngbr, excludeTypes=False)
					#~ if compat==1:
						#~ # the more recent event of this type is compatible 
						#~ ngbr.setStatus('compatible_set')
						#~ self.addProt(ngbr)
						#~ break # the for event loop
					#~ elif compat==0:
						#~ # the more recent event of this type is exclusive (e.g. two transfers from/to different locations)
						#~ break # the for event loop (enough to consider neighbor history as exclusive)
					#~ elif compat==-1:
						#~ # non-exculsive event (e.g. a duplication when considering a transfer)
						#~ continue # the for event loop
			#~ if compat: # == -1 or 1
				#~ # if completed evaluation of all events (for event loop not broken)
				#~ if compat == -1:
					#~ # without encountering an exclusive event, add the protein as a gap protein with unknown status and count the gap
					#~ self.addProt(ngbr)
					#~ nunknown += 1
				#~ #  goes to next protein on replicon
				#~ n += walkstep
			#~ else: # compat == 0
				#~ # if not (for event loop broken), stops block elongation (protein with uncompatible event was encountered)
				#~ # will be treated by next search (no need to increment n += walkstep)
				#~ break # the while loop
		return n
				
	def trimNonTransferedEnds(self, maxgapsize=None, silent=True, toEmpty=False):
		"""clean the block extremities of 'gap proteins' added while extending the block, but not legitimated by occurence of a further compatible transfered protein"""
		if maxgapsize: mgs = maxgapsize[self.eventtype()]
		else: mgs = len(self.getProtList())
		(ntrimleft, ntrimright) = (0, 0)
		# trimming left end
		self.sortByPos(referenceFirst=1)
		i = 0
		#~ while not self[self.__lapid[0]].eventtype()==self.eventtype() and (self[self.__lapid[0]].getStatus()):
		while not self[self.__lapid[0]].getStatus() in ['seed','compatible_set']:
			i += 1
			ntrimleft += 1
			if i > mgs+1:
				raise IndexError, "gap extensions should not be longer than maxgapsize=%d, protein out of range : %s"%(mgs, self.__lapid[0])
			if not silent: print 'trim', self.__lapid[0]
			self.delProt(self.__lapid[0], toEmpty=toEmpty)
			if toEmpty and not self.__lapid: break
		# trimming right end
		self.sortByPos(referenceFirst=-1)
		i = 0
		#~ while not self[self.__lapid[-1]].eventtype()==self.eventtype() and (self[self.__lapid[-1]].getStatus() in ['seed','compatible_set']):
		while not self[self.__lapid[-1]].getStatus() in ['seed','compatible_set']:
			i += 1
			ntrimright += 1
			if i > mgs+1:
				raise IndexError, "gap extensions should not be longer than maxgapsize=%d, protein out of range : %s"%(mgs, self.__lapid[0])	
			if not silent: print 'trim', self.__lapid[-1]	
			self.delProt(self.__lapid[-1], toEmpty=toEmpty)
			if toEmpty and not self.__lapid: break
		return (ntrimleft, ntrimright)
			
	def splitBlock(self, splitprot, newid, maxgapsize):
		"""splits the block into a shorter version of itself (left part) and a new one (right part)"""
		self.sortByPos()
		if isinstance(splitprot, ProtInBlock):
			splitprot = splitprot.getid()
		else:
			splitprot = str(splitprot)
		i = self.__lapid.index(splitprot)
		todelete = []
		while self[self.__lapid[i]].getStatus() not in ['seed', 'compatible_set']:
			# stores the names of the split protein and the possible right trail of proteins without signal for the event for further deletion
			todelete.append(self.__lapid[i])
			i += 1
		# creates a new block with the first encountered protein with signal for the right type of event
		newblock = LeafBlock(newid, self[self.__lapid[i]])
		todelete.append(self.__lapid[i])
		# walks rightward on the block
		i += 1
		while i < len(self.__lapid):
			# adds proteins to new block and stores their name for further deletion of self block
			newblock.addProt(self[self.__lapid[i]])
			todelete.append(self.__lapid[i])
			i += 1
		for apid in todelete:
			self.delProt(apid)
		self.trimNonTransferedEnds(maxgapsize=maxgapsize)
		return newblock
			
	def checkGapProtSignal(self, reftree, bsthreshold, dirtreepickles, cursor):
		"""tests presence of signal against a transfer in gap proteins of the block.
		
		if a protein shows signal against the transfer event, splits the block at this protein.
		"""
		ncompat = 0
		nnotagainst = 0
		lagainst = []
		# sorts the block by position on replicon
		self.sortByPos()			
		for prot in self:
			# searches 'gap' proteins in blocks, i.e. those with no transfer detected as congruent with the block event (might be caused by lack of signal for a transfer or presence of signal against a transfer)
			if prot.getStatus()=='gap':
				## checks in gene tree(s) if there is at least one significant bootstrap on the path from donor to receptor species, meaning the real absence of transfer between both species. 
				# loads tree of the subfamily of the 'gap' protein
				ngbrtree = tree2.load_pickle("%s/%s.pickle"%(dirtreepickles, prot.family()))
				# searches in gene tree the non-redundant nodes corresponding to donor and receptor clades in reference tree
				dambiguous_rec = ngbrtree.map_to_nr_nodes(reftree, lnodes=self.getRecSet(), seed=prot.getid(), returnLabels=True)#, silent=(prot.getid()=='RHIL3_1.RBSK2')	# dictionnary of reference tree receptor nodes (values) mappable to a same equivalent gene tree node (key)
				dambiguous_don = ngbrtree.map_to_nr_nodes(reftree, lnodes=self.getDonSet(), seed=prot.getid(), returnLabels=True)#, silent=(prot.getid()=='RHIL3_1.RBSK2')	# dictionnary of reference tree donor nodes (values) mappable to a same equivalent gene tree node (key)
				if not dambiguous_rec:
					raise IndexError, "cannot map any of receptors % on %s tree"%(str(seslf.getRecSet()), prot.family())
				if not dambiguous_don:
					raise IndexError, "cannot map any of donors %s on %s tree"%(str(self.getDonSet()), prot.family())
				#~ if not (dambiguous_rec and dambiguous_don):
					#~ # cannot reject transfer
					#~ prot.setStatus('not_against')
				# for all the pairs of non-redundant rec/don gene tree nodes
				#~ for ngbrreccla in dambiguous_rec:
					#~ for ngbrdoncla in dambiguous_don:
					
				for ngbrrecclab in dambiguous_rec:
					for ngbrdonclab in dambiguous_don:
						ngbrreccla = ngbrtree[ngbrrecclab]
						ngbrdoncla = ngbrtree[ngbrdonclab]
						#~ if (prot.getid()=='RHIL3_1.RBSK2'):
							#~ print "ngbrreccla.go_root()==ngbrdoncla.go_root()"
							#~ print (ngbrreccla.go_root()==ngbrdoncla.go_root())
							#~ # they can be different... skip
						# checks the bootstrap on path between both nodes
						try:	
							lbs = ngbrreccla.bootstraps_on_path_to(ngbrdoncla, excludeLeaves=True, excludeTransfers=True)[1:-1] # checks on branches of the path excluding the extremes (the latter just support the nodes' integrity)
						except IndexError, e:
							print e
							print "prot", prot.getid()
							print "dambiguous_rec", dambiguous_rec, [n for n in dambiguous_rec], ngbrrecclab
							print "dambiguous_don", dambiguous_don, [n for n in dambiguous_don], ngbrdonclab
							print "root (ngbrreccla):", [ngbrreccla.go_root()], ngbrreccla.go_root(), ngbrreccla.go_root().label()
							print "root (ngbrdoncla):", [ngbrdoncla.go_root()], ngbrdoncla.go_root(), ngbrdoncla.go_root().label()
							#~ raise IndexError, e
							continue
						for bs in lbs:
							if bs > bsthreshold:
								break # the for bs loop
						# if no supported bipartition was found for this pair of nodes, sets the protein as compatible ('not_against')
						else:
							prot.setStatus('not_against')
							prot.addNotAgainstSet([dambiguous_rec[ngbrrecclab], dambiguous_don[ngbrdonclab]])
				# if pair(s) of nodes were compatible ('not_against') with the transfer history, merge the information of compatible pairs into rec/don sets
				if prot.getStatus() == 'not_against': 
					prot.mergeNotAgainstSets()
				# if not a pair was compatible, sets the protein as 'against' the transfer history.
				else:
					prot.setStatus('against')	
					# adds the protein to the list of 'against' protein used as split separator in splitDiscordantBlock()
					lagainst.append(prot)							
		return lagainst	
				
	def splitDiscordantBlock(self, lagainst, maxgapsize, nblock, maxsplit=1):
		"""splits the block at 'against' proteins"""
		lnewblocks = []	
		nsplit = 0		
		while lagainst and nsplit < maxsplit:
			# splits from the right, as self block will be the remainnig left part
			splitprot = lagainst[-1]
			# if splitprot was not removed during a previous trimming step, splits the block
			if splitprot.getid() in self.__lapid: 
				nsplit += 1
				nblock += 1
				newblock = self.splitBlock(splitprot.getid(), nblock, maxgapsize)
				lnewblocks.append(newblock)
			del lagainst[-1]
		return lnewblocks, nblock
		
	def countStatus(self):
		ncompat = 0
		nnotagainst = 0
		nagainst = 0
		for prot in self:
			if prot.getStatus() == 'not_against': nnotagainst += 1	
			elif prot.getStatus() == 'against': nagainst += 1
			else: ncompat += 1	
		return (ncompat, nnotagainst, nagainst)
		
	def assignBlock(self, nances, nblock, dancblocks, dfam_blocks, dfam_ancs, commonsetoperation='intersection', silent=True):
		"""search to assign the focal block to an existing ancestral block sharing events with it, or create a new one.
		
		return block ids increments and a list potentially containing a new leaf block made of uncompatible part of the block.
		"""
		# updates the dictionary of blocks involving a family
		self.updateLeafBlockDict(dfam_blocks)
		# search for ancestral block sharing gene families
		ancblock, nblock, newleafblocks = self.searchAncBlockDict(dfam_ancs, dancblocks, nblock, silent=silent, setOp=commonsetoperation)
		if newleafblocks:
			if not silent: print self.getid(), "here 10"
			for block in newleafblocks:
				nances, nblock = block.assignBlock(nances, nblock, dancblocks, dfam_blocks, dfam_ancs, commonsetoperation=commonsetoperation, silent=silent)
		else:
			if not ancblock:
				if not silent: print self.getid(), "here 11"
				nances += 1
				if not silent: print "new ancestral block:"
				ancblock = AncestralBlock(nances, self)
				if not silent: print ancblock
				dancblocks[nances] = ancblock
				
			if not silent: print self.getid(), "here 12"
			ancblock.updateAncBlockDict(dfam_ancs)
		return nances, nblock
		
	def searchAncBlockDict(self, dfam_ancs, dancblocks, nblock, silent=True, setOp='intersection'):	#mergeAncBlocks=True
		"""searches d{event_id: [ancestral blocks]} for an ancestral block compatible with self block and returns it if found"""
		funlog = "!!! LeafBlock.searchAncBlockDict():"
		newleafblocks = []	# list of new leaf blocks created when disassembling uncoherent blocks
		# builds list of putative ancestral blocks sharing gene families
		lputanc = []	
		dputanc_eventid = {}
		for eventid in self.eventids(compute=True):
			if eventid in dfam_ancs[self.eventtype()]:
				lputanc += dfam_ancs[self.eventtype()][eventid]
				for putanc in dfam_ancs[self.eventtype()][eventid]:
					dputanc_eventid.setdefault(putanc, []).append(eventid)
		lputanc = list(set(lputanc))
		ancblock = None
		print "lputanc", lputanc 
		for putanc in lputanc:
			putancblock = dancblocks.get(putanc)
			# debug exception to visualize the bad entry
			if not putancblock:
				for eventid in self.eventids():
					if eventid in dfam_ancs[self.eventtype()]:
						print eventid, dfam_ancs[self.eventtype()][eventid]
				#~ raise KeyError, "Unknown entry %s in dancblocks"%putanc
				print "Unknown entry %s in dancblocks"%putanc
				continue
			cev = self.compatibleEvent(putancblock, excludeTypes=False)
			if cev>0:
				if not ancblock:
					if not silent: print self.getid(), "here 1"
					ancblock = putancblock	
					#~ if not mergeAncBlocks:
						#~ # stop iteration to avoid encountering another ancestral block which event set overlaps leaf block
						#~ print self.getid(), "here 2"
						#~ return ancblock, nblock, newleafblocks		
				else:	
					# leaf block overlaps several ancestral blocks ; they must be merged
					if not silent:
						print self.getid(), "here 3"
						print funlog, "leaf block [%s] compatible with several ancestral blocks:"%str(self.getid())
						print "ancblock: %s"%str(ancblock)
						print "putancblock: %s"%str(putancblock)
					if ancblock.compatibleEvent(putancblock):
						if not silent: print self.getid(), "here 4"
						ancblock.mergeAncBlocks(putancblock, dancblocks, dfam_ancs, silent=silent, setOp=setOp)
					else:
						if not silent: print self.getid(), "here 5"
						# if not compatible histories, do not merge the blocks ; two blocks corresponding to two alternatives hypothesis for one transfer event will be issued
						if not silent: print funlog, "Leaf block [%s] compatible with two uncompatible ancestral blocks: %s, %s"%(str(self.getid()), str(ancblock.getid()), str(putancblock.getid()))
						
						shids = ancblock.sharedeventids(putancblock)
						if shids:
							if not silent: print self.getid(), "here 51"
							# the two uncompatible ancestral blocks share events
							#~ if len(shids)/len(putancblock.eventids()) <= len(shids)/len(ancblock.eventids()):
								#~ if not silent: print self.getid(), "here 52"
								#~ a = ancblock
								#~ b = putancblock
							#~ else:
								#~ if not silent: print self.getid(), "here 53"
								#~ b = ancblock
								#~ a = putancblock
							#~ # a is the one for which this common part is the most important in proportion: will be kept as is; b will be split
							Iab = np.mean([lb.eventInternality(shids) for lb in ancblock])
							Ipab = np.mean([lb.eventInternality(shids) for lb in putancblock])
							if Iab >= Ipab:
								if not silent: print self.getid(), "here 52"
								a = ancblock
								b = putancblock
							else:
								if not silent: print self.getid(), "here 53"
								b = ancblock
								a = putancblock
							# a is the one for which this common part is the located the most internally in blocks: will be kept as is; b will be split

							uncompevents = set()
							for lb in b:
								luncprots = lb.getUncompatibleProts(a) #, silent=silent)
								if luncprots:
									if not silent:
										print self.getid(), "here 531"
										print "lb", lb.getid()
										print "luncprots", [p.getid() for p in luncprots]
									uncomplb, nblock = lb.resolveBlockUncompatibility(a, nblock, luncompatprots=luncprots, silent=silent)
									if uncomplb == lb:
										print self.getid(), "here 532"
										# the whole leaf block must be detached from ancestral block
										b.delLeafBlock(lb)
										# must remove obsolete links of events to this ancestral block
										for ei in lb.eventids():
											if (b.getid() in dfam_ancs[b.eventtype()][ei]) and (ei not in b.eventids()):
												dfam_ancs[b.eventtype()][ei].remove(b.getid())
									newleafblocks += [uncomplb]
									uncompevents |= uncomplb.eventids()
							b.commonRecDonSet()
							# the compatible part of b is merged in a; the other part is returned as component leaf blocks for further treatment
							if not silent:
								print "uncompevents", uncompevents
								print "newleafblocks", [nlb.getid() for nlb in newleafblocks]
							a.mergeAncBlocks(b, dancblocks, dfam_ancs, silent=silent, setOp=setOp)
							
							if self.eventids() & uncompevents:
								if not silent: print self.getid(), "here 54"
								# self has common eventids with the uncompatible part of b, self must be split accordingly
								luncompatprots = []
								for eventid in self.eventids():
									if eventid in uncompevents:
										luncompatprots += self.eventidgetmembers(outid)
								uncomplb, nblock = self.resolveBlockUncompatibility(a, nblock, luncompatprots=luncompatprots, silent=silent)
								newleafblocks += [self, uncomplb]
								break	# for putanc loop
							else:
								if not silent: print self.getid(), "here 55"
								# self has not been modified, continue the for putanc loop
								if not self.compatibleEvent(a):
									raise ValueError, "merging the blocks made leaf block %s uncompatible"%str(self.getid())
								ancblock = a
						else:
							if not silent: print self.getid(), "here 56"
							# the two uncompatible ancestral blocks do not share events
							# just split the leaf block to distribute events to each ancestral blocks
							ashids = self.sharedeventids(ancblock)
							bshids = self.sharedeventids(putancblock)
							if len(bshids) <= len(ashids):
								if not silent: print self.getid(), "here 57"
								a = ancblock
								b = putancblock
							else:
								if not silent: print self.getid(), "here 58"
								b = ancblock
								a = putancblock
							# a shares the most events with self, will get the most of proteins from self; b will get only prteins linked to shared events
							luncompatprots = []
							for eventid in self.eventids():
								if eventid in b.eventids():
									if not silent: print self.getid(), "here 561"
									luncompatprots += self.eventidgetmembers(eventid)
							uncomplb, nblock = self.resolveBlockUncompatibility(a, nblock, luncompatprots=luncompatprots, silent=silent)
							newleafblocks += [self, uncomplb]
							break	# for putanc loop
						
			elif cev<0:		
				# ignore blocks of different types
				continue
				if not silent: print self.getid(), "here 6"
			else:
				if not silent:
					e = "block %s is not compatible with ancestral block %s while they share the same event(s):\n"%(str(self.getid()), str(putancblock.getid()))
					for eventid in dputanc_eventid[putanc]:
						e += "%d "%(eventid)
				#~ raise IndexError, e
					print funlog, e
					print putancblock
					print self.description()
				if not silent: print self.getid(), "here 7"
				uncblock, nblock = self.resolveBlockUncompatibility(putancblock, nblock, silent=silent)
				newleafblocks += [self, uncblock]
				if not self.compatibleEvent(putancblock):
					raise ValueError, "splitting the block into two blocks did no good : block %s is still uncompatible with ancestral block %s"%(str(self.getid()), str(putancblock.getid()))
				# will stop iteration here without linking the block to any ancestral block
				# but return newleafblocks list filled with both parts of the block to try again to assign them to an ancestral block
				break
		else:
			if not silent: print self.getid(), "here 8"
			if ancblock:
				if not silent:
					print self.getid(), "here 9"
					print "assign to existing ancestral block:"
					print ancblock
				ancblock.addLeafBlock(leafblock=self, setOp=setOp)
				
		return ancblock, nblock, newleafblocks
		
	def resolveBlockUncompatibility(self, uncompancblock, nblock, luncompatprots=None, silent=True):
		"""search for part of the block inducing uncompatibility with a target ancestral block
		 and resolve it into a compatible (reduced self) and an uncompatible part (new leaf block)"""			
		if not luncompatprots:
			luncompprots = []
			for uncompbroblock in uncompancblock:
				luncompprots += self.getUncompatibleProts(uncompbroblock, silent=silent)
			if not luncompprots:
				raise IndexError, "could not find source of uncompatibility between blocks"
		else:
			luncompprots = list(luncompatprots)
		luncompprots = list(set(luncompprots))
		#~ if not silent: print "self.getProts()", [p.getid() for p in self.getProts()]
		if set(luncompprots) == set(self.getProts()):
			# cannot delete every protein of the block ; return the intact block as a new one
			if not silent: print "whole leaf block to detach", self.getid()
			return self, nblock
		if not silent: print "proteins to detach from block %s to recover compatibility with ancestral block %s:\n%s"%(str(self.getid()), str(uncompancblock.getid()), '\n'.join([str(uncompprot) for uncompprot in luncompprots]))
		# delete prots from block
		for uncompprot in luncompprots:
			self.delProt(uncompprot.getid())
		# create new block for deleted proteins
		seed = luncompprots[0]
		seed.setStatus('seed')
		nblock += 1 
		uncblock = LeafBlock(nblock, seed)
		for uncompprot in luncompprots[1:]:
			if uncblock.compatibleEvent(uncompprot):
				if uncompprot.getStatus()=='seed': uncompprot.setStatus('compatible_set')
				uncblock.addProt(uncompprot)
			else:
				raise IndexError, "protein %s should be congruent to new block %:\n%s\n%s"%(uncompprot.getid(), uncblock.getid(), str(uncompprot), uncblock.description())
		self.commonRecDonSet()
		self.trimNonTransferedEnds(silent=silent)
		uncblock.trimNonTransferedEnds(silent=silent)
		if not silent: print "new blocks:\nself\n%s\nuncblock\n%s"%(self.description(), uncblock.description())
		return uncblock, nblock
		
	def getUncompatibleProts(self, ucblock, sharedeventids=None, refrecdonset=None, silent=True):
		"""search in the block the source of uncompatibility of rec/don sets between two leaf blocks that share (a) common event(s)"""
		if not sharedeventids: sharedids = self.sharedeventids(ucblock)
		else:  sharedids = sharedeventids
		if not silent: print "sharedids", sharedids
		luncompprots = []
		#~ if sharedids:
		if not refrecdonset: refrds = ucblock.getCommonRecDonSet()	# self.getCompatiblePattern(ucblock, sharedeventids=sharedids)
		else: refrds = refrecdonset
		if not silent: print "refrds", refrds
		dprotrds = self.getDictProtToRecDonSets()
		# find outlier (non-shared) events that may cause uncompatibility
		outids = self.eventids() - sharedids
		if not silent: print "outids", outids
		for outid in outids:
			loutprot = self.eventidgetmembers(outid)
			for outprot in loutprot:
				if not compatibleRecDonSet(refrds, dprotrds[outprot]):
					luncompprots.append(outprot)		
		return luncompprots
		
	def getCompatiblePattern(self, ucblock, sharedeventids=None):
		"""compatible part of rec/don sets between two leaf blocks that share (a) common event(s)"""
		if not sharedeventids: sharedids = self.sharedeventids(ucblock)
		else:  sharedids = sharedeventids
		if not sharedids:
			return None
		# define the rec/don coordinates of shared events
		dselfprotrds = self.getDictProtToRecDonSets()
		lsharedids = list(sharedids)
		seed = self.eventidgetmembers(lsharedids[0])[0]
		refrecset, refdonset = dselfprotrds[seed]
		for sharedid in lsharedids[1:]:
			#~ for prot in self.eventidgetmembers(sharedid):
			rds = dselfprotrds[self.eventidgetmembers(sharedid)[0]]
			refrecset &= rds[0]
			if refdonset: refdonset &= rds[1]
		return (refrecset, refdonset)
		
		
	def updateLeafBlockDict(self, dfam_blocks):
		"""updates the dictionary of blocks involving an event at a node in a gene family tree"""
		for eventid in self.eventids(compute=True):
			dfam_blocks.setdefault(eventid, []).append(self.getid())
			
	def description(self):
		s = str(self)+'\n'
		for prot in self:
			s+= str(prot)
			# searches 'gap proteins' in blocks, i.e. those with no transfer detected (might be caused by lack of signal for a transfer or presence of signal against a transfer)
			if prot.getStatus() == 'against':
				s += '\t// %s\n'%(str(prot.getStatus()).upper())
			else:
				s += '\t// %s\t[%s -> %s]\n'%( str(prot.getStatus()).upper(), ', '.join(list(prot.getDonSet())), ', '.join(list(prot.getRecSet())) )
		return s
		
	def writeBlockToFile(self, fout):
		fout.write(self.description()+'\n')
		
class AncestralBlock(BlockEvent):
	def __init__(self, ancblockid, leafblock):
		if not isinstance(leafblock, LeafBlock):
			raise TypeError, 'expected a LeafBlock object'
		super(AncestralBlock, self).__init__(ancblockid, leafblock)
		self.__famset = set(leafblock.getFamList(filterStatus=['against']))
		seedid = leafblock.getid()
		self.__lblock = [seedid]
		self.__dblock = dict({seedid: leafblock})
		
	def __str__(self):
		d = self.getFamStatusDict()
		s = "[ancestral block %s]\t%s [%s -> %s] :\n"%(self._id, self.eventtype(), " ".join(list(self._donset)), " ".join(list(self._recset)) )
		fams = []
		for key in ['compatible', 'gap', 'against']:
			fams.append("%s : %s"%(key, " ".join(d[key]) ) )
		s += "\t// ".join(fams)+"\n"
		s += "events: "+" ".join([str(eventid) for eventid in self.eventids()])+"\n"
		for nblock in self.__dblock:
			if nblock == self._seedid: isSeed = 'seed'
			else :  isSeed = ''
			s += "%s\t%s\n"%(isSeed, str(self.__dblock[nblock]))
		s.rstrip('\n')
		return s
		
	def __getitem__(self, nblock):
		return self.__dblock[nblock]
		
	def __iter__(self):
		return self.generator()
		
	def generator(self):
		for nblock in self.__lblock:
			yield self.__dblock[nblock]
			
	def getFamSet(self):
		return self.__famset
			
	def getFamList(self, filterStatus=[]):
		famlist = []
		for leafblock in self:
			famlist += leafblock.getFamList(filterStatus=filterStatus)
		return list(set(famlist))
		
	def getFamStatusDict(self):
		# list of families with at least one member (in any leaf block) with detected transfer and receptor / donor sets compatible to ancestral block
		lcompat = self.getFamList(filterStatus=['against', 'not_against', 'gap'])
		# list of families which no member has a compatible detected transfer and with at least one member being 'against'
		lagainst = list( set(self.getFamList(filterStatus=['seed', 'compatible_set','compatible_nni', 'not_against'])) - set(lcompat) )
		# list of families which members being only 'not_against'
		lnotagainst = list( set(self.getFamList(filterStatus=['seed', 'compatible_set','compatible_nni', 'against'])) - set(lcompat) - set(lagainst) )
		# list of families which members being only 'gap'
		lgap = list( set(self.getFamList(filterStatus=['seed', 'compatible_set','compatible_nni', 'against', 'not_against'])) - set(lcompat) - set(lagainst) - set(lnotagainst) )
		d = dict({'compatible': lcompat, 'not_against': lnotagainst, 'against': lagainst, 'gap':lgap})
		return d
		
	def getBlockList(self):
		return self.__lblock
			
	def addLeafBlock(self, leafblock, setOp='intersection'): 
		"""adds a LeafBlock object"""
		if not isinstance(leafblock, LeafBlock):
			raise TypeError, 'expected a LeafBlock object'
		nblock = leafblock.getid()
		if nblock in self.__lblock:
			raise IndexError, "cannot add twice a block"
		if len(self.__lblock)==0:
			# add to an empty block -> initiate a new one
			self = AncestralBlock(self.getid(), leafblock)
		else:
			self.__lblock.append(nblock)
			self.__dblock[nblock] = leafblock	
			fams = leafblock.getFamList(filterStatus=['against'])
			self.__famset |= set(fams)
			self.add_eventids(leafblock.eventids())
			if setOp=='intersection':
				self._recset &= leafblock.getRecSet()
				self._donset &= leafblock.getDonSet()
			elif setOp=='union':
				self._recset |= leafblock.getRecSet()
				self._donset |= leafblock.getDonSet()
			else:
				ValueError, 'wrong set operation definition'
			
	def delLeafBlock(self, leafblock): 
		"""deletes a LeafBlock object"""
		if not isinstance(leafblock, LeafBlock):
			raise TypeError, 'expected a LeafBlock object'
		nblock = leafblock.getid()
		if not nblock in self.__lblock:
			raise IndexError, "cannot delete a block not in self"
		self.__lblock.remove(nblock)
		del self.__dblock[nblock]
		if len(self.__lblock)>0:
			if not (self._seedid in self.__dblock):
				# change the seed to a current leaf block
				print "!!! AncestralBlock.delLeafBlock(): seed absent from ancestral block", str(self)
				self._seedid = self.__lblock[0]
				print "!!! AncestralBlock.delLeafBlock(): changed seed in ancestral block", str(self)
			self.getCommonRecDonSet(compute=True)
			self.__famset = set(self.getFamList(filterStatus=['against']))
			self.eventids(compute=True)
		else:
			print "!!! AncestralBlock.delLeafBlock(): created empty ancestral block", str(self)
			# empties the block
			self._recset = set([])
			self._donset = set([])
			self.__famset = set([])
			self._eventids = set([])
			self._seedid = None
		
	def searchLeafBlockDict(self, dfam_blocks, dblocks, setOp='intersection'):
		"""searches d{(family, gene tree node id): [leaf blocks]} for leaf blocks compatible with (self) ancestral block ; returns True if found, False otherwise"""
		# builds list of putative brother leaf blocks sharing gene families
		lputbro = []
		for eventid in self.eventids(compute=True):
			lputbro += dfam_blocks[self.eventtype()][eventid]
		lputbro = list(set(lputbro) - set(self.__lblock))
		foundleaf = False
		for putbro in lputbro:
			putbblock = dblocks[putbro]
			if self.compatibleEvent(putbblock):
				self.addLeafBlock(putbblock, setOp=setOp)
				foundleaf = True
		return foundleaf
		
	def updateAncBlockDict(self, dfam_ancs):
		"""updates the dictionary of blocks involving an event at a node in a gene family tree"""
		for eventid in self.eventids(compute=True):
			dfam_ancs[self.eventtype()][eventid] = dfam_ancs.setdefault(eventid, []) + [self.getid()]
			
	def mergeAncBlocks(self, putancblock, dancblocks, dfam_ancs, silent=True, setOp='intersection', allowpartcompat=True):
		if self.compatibleEvent(putancblock) or allowpartcompat:
			# !!! allowpartcompat=True can lead to loss of leaf blocks
			putanc = putancblock.getid()
			if not silent: print "merge ancestral blocks [%s] and [%s]"%(str(self.getid()), str(putanc))
			for putbroblock in putancblock:
				if (not allowpartcompat) or self.compatibleEvent(putbroblock):
					if not putbroblock.getid() in self.getBlockList():
						self.addLeafBlock(leafblock=putbroblock, setOp=setOp)
			if not silent:print "into a new ancestral block: %s"%str(self)
			for eventid in putancblock.eventids():
				if eventid in dfam_ancs:
					if putanc in dfam_ancs[eventid]:
						dfam_ancs[eventid].remove(putanc)
			del dancblocks[putanc]
		else:
			raise IndexError, "ancestral blocks to be merged are not compatible: %s\n%s"%(str(self), str(putancblock))
				
class RepliconMap(object):
	"""Mapping object of blocks on a replicon, given their (begin, end) coordinates."""
	def __init__(self, rid):
		self.__id = rid
		self.__lblocks = []
		self.__dblocks = {}
		self.__dcoord_blocks = {}
		self.__dbeg_blocks = {}
		self.__dend_blocks = {}
		self.__lbeg = []
		self.__lend = []
		self.__lcoord = []
		
	def __getitem__(self, n):
		return self.__dblocks[n]	
		
	def __str__(self):
		l = []
		self.__lcoord.sort()
		for coord in self.__lcoord:
			l.append('%s\t%s'%(str(coord), str(self.__dcoord_blocks[coord])))
		return '\n'.join(l)
			
	def __iter__(self):
		return self.generator()
		
	def generator(self):
		self.__lblocks.sort()
		for nblock in self.__lblocks:
			yield self.__dblocks[nblock]
			
	def getid(self):
		return self.__id
			
	def getCoordList(self):
		self.__lcoord.sort()
		return self.__lcoord
		
	def getCoordBlockDict(self):
		return self.__dcoord_blocks
	
	def addBlock(self, leafblock):
		nblock = leafblock.getid()
		if nblock in self.__lblocks:
			raise IndexError, "block %s already in replicon:\n%s"%(nblock, leafblock.description())
		self.__lblocks.append(nblock)
		self.__dblocks[nblock] = leafblock
			
	def addBlockCoords(self, leafblock, silent=True):
		nblock = leafblock.getid()
		coord = leafblock.getCoords()
		if not silent:
			print '\tadd coords', coord, leafblock
		self.__dcoord_blocks[coord] = self.__dcoord_blocks.setdefault(coord, []) + [nblock]
		self.__dbeg_blocks[coord[0]] = self.__dbeg_blocks.setdefault(coord[0], []) + [nblock]
		self.__dend_blocks[coord[1]] = self.__dend_blocks.setdefault(coord[1], []) + [nblock]
		if not coord in self.__lcoord:
			self.__lcoord.append(coord)
		if not coord[0] in self.__lbeg:
			self.__lbeg.append(coord[0])
		if not coord[1] in self.__lend:
			self.__lend.append(coord[1])
		
	def delBlockCoords(self, leafblock, silent=True):
		nblock = leafblock.getid()
		coord = leafblock.getCoords()
		if not silent:
			print '\tdel coords', coord, leafblock
		self.__dcoord_blocks[coord].remove(nblock)
		if self.__dcoord_blocks[coord] == []:
			del self.__dcoord_blocks[coord]
			self.__lcoord.remove(coord)
		self.__dbeg_blocks[coord[0]].remove(nblock)
		if self.__dbeg_blocks[coord[0]] == []:
			del self.__dbeg_blocks[coord[0]]
			self.__lbeg.remove(coord[0])
		self.__dend_blocks[coord[1]].remove(nblock)
		if self.__dend_blocks[coord[1]] == []:
			del self.__dend_blocks[coord[1]]
			self.__lend.remove(coord[1])		
		
	def delBlock(self, leafblock, silent=True):
		nblock = leafblock.getid()
		coord = leafblock.getCoords()
		if not silent:
			print '\tdel', coord, leafblock
		self.__lblocks.remove(nblock)
		del self.__dblocks[nblock]
		self.delBlockCoords(leafblock)
			
	def getBlockList(self):
		return self.__lblocks
		
	def purgeBlockList(self, reftree, silent=True):
		"""deletes redundant blocks from the record.
		
		Iters on blocks. Compared to the reference block, a block with coordinates equal or contained in the reference block coordinates is evaluated.
		If its common recetor/donor sets is the same 
		or if its protein set is the same or contained in the reference block protein set, 
		the block is redundant and is deleted.
		
		NB: doesn't deal with overlapping blocks ; need a fusion method.
		"""
		status = 0
		lnblocks = self.__lblocks
		for nblock in lnblocks:
			if not nblock in self.__dblocks:
				# if block has been deleted at a previous step, skips
				continue
			# set the reference block variables
			block = self[nblock]
			coord = block.getCoords()
			recdonset = block.getCommonRecDonSet()
			# sorts the indexes
			self.__lbeg.sort()
			self.__lend.sort()
			# get the possible coordinates of block contained in reference
			localbegs = []
			localends = []
			i = self.__lbeg.index(coord[0])
			while i < len(self.__lbeg) and self.__lbeg[i] < coord[1]:
				localbegs.append(self.__lbeg[i])
				i += 1
			j = self.__lend.index(coord[1])
			while j >= 0 and self.__lend[j] > coord[0]:
				localends.append(self.__lend[j])
				j -= 1
			# get list of contained blocks to evaluate
			lnb = []
			for beg in localbegs:
				for end in localends:
					if (beg, end) in self.__dcoord_blocks:
						lnb += self.__dcoord_blocks[(beg, end)]
			
			try:
				lnb.remove(nblock)
			except ValueError, e:
				print nblock, localbegs, localends, lnb
				print block.description()
				raise ValueError, e
			if lnb:
				if not silent: print 'ref', coord, block
				# evaluate the block
				for n in lnb:
					b = self[n]
					if b.eventtype()!=block.eventtype():
						continue
					b.commonRecDonSet()
					refmacrec = block.getMACReceptor(reftree)
					altmacrec = b.getMACReceptor(reftree)
					refmacdon = block.getMACDonor(reftree)
					altmacdon = b.getMACDonor(reftree)
					# evaluate the protein set of alternative block; considered redundant if smaller than protein set of reference block (favoring large gene block histories)
					if set(b.getProtList()) <= set(block.getProtList()):
						# evaluate the receptor / donor sets ; considered redundant if larger than or equal to reference one (less focused)					
						if ( b.getCommonRecDonSet(compute=False) >= recdonset ):
							# if redundant, deletes block and regenerates indexes
							if not silent: print b.getCommonRecDonSet(compute=False), recdonset
							self.delBlock(b, silent=silent)
							status = 1
		return status	
		
### generic compatibility test function
def compatibleRecDonSet(recdonset1, recdonset2, returnRecDon=False):
	"""if Receptor and Donor sets in tuples show non-null intersection, returns those intersections or boolean stating the compatibility"""
	#~ rds = [recdonset1[i] & recdonset2[i] for i in range(len(recdonset1))]
	commonRec = set(recdonset1[0]) & set(recdonset2[0])
	if commonRec:
		if recdonset1[1]:
			# case of a transfer event, a second coordinate is needed
			commonDon = set(recdonset1[1]) & set(recdonset2[1])
			if commonDon:
				if not returnRecDon: return 1
				else: return (commonRec, commonDon)
		else:
			if not returnRecDon: return 1
			else: return (commonRec, set())	
	else:
		if not returnRecDon: return 0
		else: return set()
		
def differenceRecDonSet(refrecdonset, testrecdonset):
	drds = [refrecdonset[i] - testrecdonset[i] for i in range(len(refrecdonset))]
	return drds
		
# moved to rec_to_db		
#~ #### generic query functions from Phylariane-like PostgreSQL database.
#~ def getGeneInfo(apid, fields, cursor, returnDict=True):
	#~ """queries the gene informations from a database. Returns adictionary."""
#~ #	fields = ['gene_id', 'replicon', 'genomic_begin', 'genomic_end', 'hogenom_gene_id', 'locus_tag', 'family_accession']
	#~ query = "SELECT %s FROM genome.gene AS G "%(', '.join(fields))
	#~ query += "WHERE G.hogenom_gene_id='%s'"%(apid)
	#~ cursor.execute( query )
	#~ cols = rec_to_db.colnames(cursor)
	#~ tgene = cursor.fetchone()
	#~ if not returnDict:
		#~ return tgene
	#~ else:
		#~ dgene = dict(zip(cols, tgene))
		#~ return dgene
	#~ 
#~ def getGeneTreeInfo(fam, fields, cursor, returnDict=True):
	#~ """queries the gene tree informations from a database. Returns a dictionary."""
	#~ query = "SELECT %s FROM phylogeny.gene_tree INNER JOIN phylogeny.reconciled_gene_tree ON gene_tree_id=input_gene_tree_id "%(', '.join(fields))
	#~ query += "WHERE family_accession='%s';"%(fam)
	#~ cursor.execute( query )
	#~ cols = rec_to_db.colnames(cursor)
	#~ tgenetree = cursor.fetchone()
	#~ if not returnDict:
		#~ return tgenetree
	#~ else:
		#~ dgentree = dict(zip(cols, tgenetree))
		#~ return dgenetree
	#~ 
#~ def getFamilyInfo(fam, fields, cursor, returnDict=True):
	#~ """queries the gene tree informations from a database. Returns a dictionary."""
	#~ query = "SELECT %s FROM genome.gene INNER JOIN phylogeny.gene_tree USING (family_accession) LEFT JOIN phylogeny.reconciled_gene_tree ON gene_tree_id=input_gene_tree_id "%(', '.join(fields))
	#~ query += "WHERE family_accession='%s';"%(fam)
	#~ cursor.execute( query )
	#~ cols = rec_to_db.colnames(cursor)
	#~ if not returnDict:
		#~ ttfamily = cursor.fetchall()
		#~ return ttfamily
	#~ else:
		#~ tdfamily = ()
		#~ for tfamily in cursor:
			#~ tdfamily += (dict(zip(cols, tgenetree)),)
		#~ return tdfamily
	#~ 
#~ def getEventInfo(recgiid, startnodeid, fields, cursor, returnDict=True):
	#~ """queries the event informations from a database. Returns a tuple of dictionaries."""
	#~ query = "SELECT %s FROM phylogeny.event "%(', '.join(fields))
	#~ query += "WHERE rec_gi_id=%d AND rec_gi_start_node=%d;"%(recgiid, startnodeid)
	#~ cursor.execute( query )
	#~ cols = rec_to_db.colnames(cursor)
	#~ if not returnDict:
		#~ ttevents = cursor.fetchall()
		#~ return ttevents
	#~ else:
		#~ tdevent = ()
		#~ for tevent in cursor:
			#~ devent = dict(zip(cols, tevent))
			#~ if 'event_type' in devent:	devent['eventtype'] = rec_to_db.dabrev_event[devent['event_type']]
			#~ tdevent += (devent,)
		#~ return tdevent
	#~ 
#~ def getEventLocation(eventid, eventtype, cursor):
	#~ """perform query to get event location details in a format similar to the 'eventlocation' part of the return value of GeneTree.getdicevent()"""
	#~ fields = ['number', 'sp_node_id', 'characteristic', 'reference_node']
	#~ query = "SELECT %s "%', '.join(fields)
	#~ query += "FROM phylogeny.species_node "
	#~ query += "INNER JOIN phylogeny.event_possible_species_node USING (sp_node_id) "
	#~ query += "WHERE event_id=%d ;"%eventid
	#~ cursor.execute( query )
	#~ cols = rec_to_db.colnames(cursor)
	#~ ttlocs = cursor.fetchall()
	#~ if eventtype=='transfer':
		#~ lrec = []
		#~ ldon = []
		#~ recid = None
		#~ donid = None
		#~ # fetch transfer child id
		#~ query = "SELECT rec_gi_end_node from phylogeny.event WHERE event_id=%d ;"%eventid
		#~ cursor.execute( query )
		#~ childid = cursor.fetchone()[0]
	#~ else:
		#~ lloc = []
		#~ refid = None
	#~ for tloc in ttlocs:
		#~ dloc = dict(zip(cols, tloc))
		#~ c = dloc['characteristic']
		#~ if c=='location':
			#~ lloc.append(dloc['number'])
			#~ if dloc['reference_node']==True:
				#~ refid = dloc['sp_node_id']
		#~ elif c=='rec':
			#~ lrec.append(dloc['number'])
			#~ if dloc['reference_node']==True:
				#~ recid = dloc['sp_node_id']
		#~ elif c=='don':
			#~ ldon.append(dloc['number'])
			#~ if dloc['reference_node']==True:
				#~ donid = dloc['sp_node_id']
	#~ if eventtype=='transfer':
		#~ return (recid, lrec, donid, ldon, childid)
	#~ else:
		#~ return (refid, lloc)
	#~ 
		#~ 
#~ def getEventLocsInfo(eventid, fields, cursor, returnDict=True):
	#~ """queries the event informations from a database. Returns a tuple of dictionaries."""
	#~ query = "SELECT %s FROM phylogeny.event_locs "%(', '.join(fields))
	#~ query += "WHERE event_id=%d;"%(eventid)
	#~ cursor.execute( query )
	#~ cols = rec_to_db.colnames(cursor)
	#~ if not returnDict:
		#~ ttevent = cursor.fetchall()
		#~ return ttevent
	#~ else:
		#~ tdevent = ()
		#~ for tevent in cursor:
			#~ tdevent += (dict(zip(cols, tevent)),)
		#~ return tdevent
		#~ 
#~ def getTreeEventsInfo(fam, startnodeid, fields, cursor, returnDict=True):
	#~ """queries the event informations from a database. Returns a tuple of dictionaries."""
	#~ query = "SELECT %s FROM phylogeny.tree_events "%(', '.join(fields))
	#~ query += "WHERE family_accession='%s' AND rec_gi_start_node=%d;"%(fam, startnodeid)
	#~ cursor.execute( query )
	#~ cols = rec_to_db.colnames(cursor)
	#~ if not returnDict:
		#~ ttevent = cursor.fetchall()
		#~ return ttevent
	#~ else:
		#~ tdevent = ()
		#~ for tevent in cursor:
			#~ tdevent += (dict(zip(cols, tevent)),)
		#~ return tdevent
		#~ 
#~ def getTreeEventRefLocsInfo(fam, startnodeid, fields, cursor, returnDict=True):
	#~ """queries the event informations from a database. Returns a tuple of dictionaries."""
	#~ query = "SELECT %s FROM phylogeny.tree_event_reflocs "%(', '.join(fields))
	#~ query += "WHERE family_accession='%s' AND rec_gi_start_node=%d;"%(fam, startnodeid)
	#~ cursor.execute( query )
	#~ cols = rec_to_db.colnames(cursor)
	#~ if not returnDict:
		#~ ttevent = cursor.fetchall()
		#~ return ttevent
	#~ else:
		#~ tdevent = ()
		#~ for tevent in cursor:
			#~ tdevent += (dict(zip(cols, tevent)),)
		#~ return tdevent
		
### block output functions
def writeAncestralBlockEvents(dancblocks, reftree, nfancblocks, nfmatrix, nfancstats, mincommonfams=2):
	"""write ancestral transfered blocks to ouput file and transfered event summary matrix"""
	statFields = ['nblock','receptor','donor','compatible','not_against','against', 'gap']
	lancestral = dancblocks.values() # list of AncestralBlock objects
	## preparing summary matrix of transfer events
	lcla = reftree.get_children_labels()
	dcla = dict( zip( lcla, range(len(lcla)) ))
	mattrans = np.zeros( (len(lcla),len(lcla)) )			
	## preparing ancestral block descriptions output file
	fancblocks = open(nfancblocks, 'w')
	fancstats = open(nfancstats, 'w')
	fancstats.write('\t'.join(statFields)+'\n')
	for ancblock in lancestral:
		macrec = ancblock.getMACReceptor(reftree, compute=False, returnLabel=True)
		macdon = ancblock.getMACDonor(reftree, compute=False, returnLabel=True)
		if ancblock.eventtype() == 'transfer':
			## filling summary matrix of block transfer events
			mattrans[ dcla[macdon], dcla[macrec] ] += 1
		## writing ancestral block description in output file for non-leaf blocks of minimum size 'mincommonfams'
		if len(ancblock.getFamSet()) >= mincommonfams and (not reftree[macrec].is_leaf()):
			fancblocks.write('\n%s\t\n'%str(ancblock))
		## writing ancestral block statistics
		d = ancblock.getFamStatusDict()
		fancstats.write('%d\t%s\t%s\t%s\n'%(ancblock.getid(), macrec, macdon, '\t'.join([str(len(d[k])) for k in ['compatible', 'not_against', 'against', 'gap']]) ) )
	fancblocks.close()
	fancstats.close()
	# writing summary matrix of block transfer events
	fmatout = open(nfmatrix, "w")
	fmatout.write("\t%s\n"%("\t".join(lcla)))
	for cla in lcla:
		l = mattrans[dcla[cla],].tolist()
		fmatout.write(cla)
		for t in l:
			fmatout.write("\t%s"%(str(t)) )
		fmatout.write("\n")
	fmatout.close()
	
def writeSQLBlockTables(dancblocks, dirsqlout, reftree, cursor, outputmode='tabular'):
	"""write SQL dump file to update tables concerning leaf and ancestral transfered blocks"""
	
	if outputmode == 'tabular':
		nfab = "%s/ancestral_block.tab"%dirsqlout
		fab = open(nfab, 'w')
		nflb = "%s/leaf_block.tab"%dirsqlout
		flb = open(nflb, 'w')
		nfbpsn = "%s/block_possible_species_node.tab"%dirsqlout
		fbpsn = open(nfbpsn, 'w')
		nfg2b = "%s/gene2block.tab"%dirsqlout
		fg2b = open(nfg2b, 'w')
		nff2b = "%s/fam2block.tab"%dirsqlout
		ff2b = open(nff2b, 'w')	
		nftempevt = "%s/temp4update_event.tab"%dirsqlout
		ftempevt = open(nftempevt, 'w')	
		fscript = open("%s/allblocktables_tabdump.sql"%dirsqlout, 'w')
	elif outputmode == 'SQL':
		fsqlout = open("%s/allblocktables_sqldump.sql"%dirsqlout, 'w')
	else:
		raise ValueError, "not a valid mode %s"%outputmode
		
	def block_species_node_location(block, blocktype, reftree):
		blockid = block.getid()
		macreclab = block.getMACReceptor(reftree, compute=False, returnLabel=True)
		macdonlab = block.getMACDonor(reftree, compute=False, returnLabel=True)
		for reclab in block.getRecSet():
			if reclab == macreclab: refnodebool = True
			else: refnodebool = False
			recid = reftree[reclab].nodeid()
			if outputmode == 'SQL':
				rec_to_db.writeSQLinsert(fsqlout, 'blocks.block_possible_species_node', [blockid, recid, 'rec', blocktype, refnodebool] ) 
			elif outputmode == 'tabular':
				rec_to_db.fillDataTableFile(fbpsn, 'blocks.block_possible_species_node', [blockid, recid, 'rec', blocktype, refnodebool] ) 
		for donlab in block.getDonSet():
			if donlab == macdonlab: refnodebool = True
			else: refnodebool = False
			donid = reftree[donlab].nodeid()
			if outputmode == 'SQL':
				rec_to_db.writeSQLinsert(fsqlout, 'blocks.block_possible_species_node', [blockid, donid, 'don', blocktype, refnodebool] ) 
			elif outputmode == 'tabular':
				rec_to_db.fillDataTableFile(fbpsn, 'blocks.block_possible_species_node', [blockid, donid, 'don', blocktype, refnodebool] ) 
	
	if outputmode == 'SQL':
		fsqlout.write("TRUNCATE TABLE blocks.leaf_block;\nTRUNCATE TABLE blocks.ancestral_block;\n")
	
	assignedevents = []	# list of event entries in DB that are assigned to an ancestral block
	lancestral = dancblocks.values() # list of AncestralBlock objects
	for ancblock in lancestral:
		ancblockid = ancblock.getid()		
		ancEtype = ancblock.eventtype()[0].upper()
		# writing ancestral block entry in SQL file
		if outputmode == 'SQL':
			dacid = dict(anc_block_id=ancblockid)
			rec_to_db.writeSQLinsert(fsqlout, 'blocks.ancestral_block', [ancblockid, ancEtype])
		elif outputmode == 'tabular':
			rec_to_db.fillDataTableFile(fab, 'blocks.ancestral_block', [ancblockid, ancEtype])
		# writing leaf block species node location entries
		block_species_node_location(ancblock, 'ancestral', reftree)
		# writing family to ancestral block mapping information
		lfam = ancblock.getFamList()
		for fam in lfam:
			if outputmode == 'SQL':
				rec_to_db.writeSQLinsert(fsqlout, 'blocks.fam2block', [ancblockid, fam])
			elif outputmode == 'tabular':
				rec_to_db.fillDataTableFile(ff2b, 'blocks.fam2block', [ancblockid, fam])
		for block in ancblock:
			blockid = block.getid()
			# writing leaf block entries
			if outputmode == 'SQL':
				rec_to_db.writeSQLinsert(fsqlout, 'blocks.leaf_block', [blockid, ancblockid] ) 
			elif outputmode == 'tabular':
				rec_to_db.fillDataTableFile(flb, 'blocks.leaf_block', [blockid, ancblockid] ) 
			# writing leaf block species node location entries
			block_species_node_location(block, 'leaf', reftree)		
			for prot in block:
				# writing gene to leaf block mapping information
				if outputmode == 'SQL':
					rec_to_db.writeSQLinsert(fsqlout, 'blocks.gene2block', [blockid, prot['gene_id']])
				elif outputmode == 'tabular':
					rec_to_db.fillDataTableFile(fg2b, 'blocks.gene2block', [blockid, prot['gene_id']])
				# get the corresponding events present in DB
				startnodeid = int(prot.eventnode().nodeid())
				fam = prot.family()
				eventtype = prot.eventtype()
				tdevent = rec_to_db.getTreeEventsInfo(fam, startnodeid, '*', cursor)
				recgiid = int(tdevent[0]['rec_gi_id'])
				for devent in tdevent:
					eventid = int(devent['event_id'])
					if not eventid in assignedevents:
						tdeventloc = rec_to_db.getEventLocsInfo(eventid, '*', cursor)
						lrec = []
						ldon = []
						for deventloc in tdeventloc:
							if deventloc['characteristic']=='don':
								ldon.append(deventloc['number'])
							# elif deventloc['characteristic'] in ['rec', 'location']:
							else:
								lrec.append(deventloc['number'])
						if ancblock.compatibleRecDonSet(set(lrec), set(ldon)):
							if outputmode == 'SQL':
								rec_to_db.writeSQLupdate(fsqlout, 'phylogeny.event', "WHERE event_id=%d"%(eventid), dacid)
							elif outputmode == 'tabular':
								#
								#	!!! this will create an array to be used as a temporary table for updating the existing `event` table ; proceed as follow:
								#
								devent['anc_block_id'] = ancblockid
								rec_to_db.fillDataTableFile(ftempevt, 'phylogeny.event', devent)
								
							assignedevents.append(eventid)
			
	if outputmode == 'SQL':
		fsqlout.close()
	elif outputmode == 'tabular':
		fab.close()
		script  = "\copy blocks.ancestral_block FROM '%s'\n"%nfab
		flb.close()
		script += "\copy blocks.leaf_block FROM '%s'\n"%nflb
		fbpsn.close()
		script += "\copy blocks.block_possible_species_node FROM '%s'\n"%nfbpsn
		fg2b.close()
		script += "\copy blocks.gene2block FROM '%s'\n"%nfg2b
		ff2b.close()
		script += "\copy blocks.fam2block FROM '%s'\n"%nff2b
		ftempevt.close()
		script += "CREATE TEMPORARY TABLE temp4update_event (LIKE phylogeny.event) ;\n"
		script += "\copy temp4update_event (%s) FROM '%s'\n"%(', '.join(rec_to_db.dtable_col['phylogeny.event']), nftempevt)
		script += "-- beware of removing all foreing key constraints referring to phylogeny.event before the DELETE step;\n"
		#~ script += "ALTER TABLE phylogeny.event DISABLE TRIGGER USER;\n"
		script += "ALTER TABLE phylogeny.event_possible_species_node DROP CONSTRAINT event_possible_species_node_event_id_fkey;\n"
		script += "ALTER TABLE phylogeny.represented_event DROP CONSTRAINT represented_event_event_id_fkey;\n"
		script += "ALTER TABLE phylogeny.transfer_confidence DROP CONSTRAINT transfer_confidence_event_id_fkey;\n"	
		script += "\d phylogeny.event\n"	
		script += "DELETE FROM phylogeny.event WHERE (event_id) IN (SELECT event_id FROM temp4update_event);\n"
		script += "INSERT INTO phylogeny.event (SELECT * FROM temp4update_event);\n"
		#~ script += "ALTER TABLE phylogeny.event ENABLE TRIGGER USER;\n"
		script += "ALTER TABLE phylogeny.event_possible_species_node ADD FOREIGN KEY (event_id) REFERENCES phylogeny.event (event_id) ON DELETE CASCADE ;\n"
		script += "ALTER TABLE phylogeny.represented_event ADD FOREIGN KEY (event_id) REFERENCES phylogeny.event (event_id) ON DELETE CASCADE ;\n"
		script += "ALTER TABLE phylogeny.transfer_confidence ADD FOREIGN KEY (event_id) REFERENCES phylogeny.event (event_id) ON DELETE CASCADE ;\n"
		#~ script += "ALTER TABLE phylogeny.event ADD FOREIGN KEY (anc_block_id) REFERENCES blocks.ancestral_block (anc_block_id) ON DELETE RESTRICT ;\n"
		#~ script += "DROP TABLE temp4update_event;\n"
		fscript.write(script)
		fscript.close()
		
def writeXMLBlockTreeCollections(dancblocks, dirxmltrees, dirconcattrees):
	for ancblockid in dancblocks:
		ancblock = dancblocks[ancblockid]
		lnfxmltrees = ancblock.getXMLTreeCollection(dirxmltrees)
		nfouttree = "%s/%s.xml"%(dirconcattrees, ancblockid)
		tree2.concat_phyloxml(lnfxmltrees, nfouttree)

