import os, sys
import shutil

def read_fasta(nffasta):
	ffasta = open(nffasta, 'r')
	fasta = ffasta.readlines()
	ffasta.close()
	return fasta
	
def write_sub_aln(fasta, seqs, nfout, keepProtLabels=False):
	#~ fout = open("%s/%s.%s"%(outfastadir, nust, alnext), 'w')
	fout = open(nfout, 'w')
	i = 0
	while i < len(fasta):
		if fasta[i].startswith('>'):
			#print fasta[i]
			if ' ' in fasta[i]:
				[prot, nr, residue] = fasta[i].lstrip('>').split()
			else:
				prot = fasta[i].strip('>\n')
				nr, residue = (None, None)
			[code, rep_peid] = prot.split('_')
			#print prot
			if prot in seqs:
				if keepProtLabels: p = prot
				else: p = code
				# write only the code of the species in alignment, as Prunier take in input a gene tree with species labels
				if residue:
					fout.write('>%s	       %s %s\n'%(p, nr, residue))
				else:
					fout.write('>%s\n'%p)
				i += 1
				while i < len(fasta) and not fasta[i].startswith('>'):
					seql = ''.join(fasta[i].split()).rstrip('\n')	# transform inter-spaced multi-line sequence into non-spaced mono-line sequence
					fout.write(seql)
					i += 1
				else:
					fout.write('\n')	# write EOL at the end of one sequence
			else:
					i += 1
		else:
			i += 1
	fout.close()



if len(sys.argv) < 5:
	print "Usage: python extract_sublan.py { path/subfam_seqdict_dir | path/subfam_seqlist-table } path/fasta_dir aln_ext path/output_fasta_dir [keep_protein_labels=True]"
	sys.exit(2)

famlistdir = sys.argv[1]
fastadir = sys.argv[2]
alnext = sys.argv[3]
outfastadir = sys.argv[4]
if len(sys.argv) > 5:
	try:
		s = sys.argv[5]
		trybool = s[0].upper()+s[1:].lower()
		keepProtLabels = eval(trybool)
	except NameError:
		keepProtLabels = False
else:
	keepProtLabels = False

if os.access(outfastadir, os.F_OK):
	shutil.rmtree(outfastadir)
os.mkdir(outfastadir)

# sample unicopy sub-tree leaves sequences from full family fasta (alignment) files to fasta (alignment) files 
n = 0
if  os.access(famlistdir, os.X_OK):
	# case where 'famlistdir' is a ditrectory
	print "list directory", famlistdir
	lfam = os.listdir("%s"%(famlistdir))
	for nfamdict in lfam:
		nfam = nfamdict.split('.')[0]
		nffasta = "%s/%s.%s"%(fastadir, nfam, alnext)
		if not os.access(nffasta, os.F_OK): continue
		n += 1
		sys.stdout.write("\r%d\t\t%s\t\t"%(n, nfam))
		nffam = "%s/%s"%(famlistdir, nfamdict)
		ffam = open(nffam, 'r')
		fasta = read_fasta(nffasta)
		for line in ffam:
			lsp = line.rstrip('\n').split('\t')
			nust = lsp[0]
			seqs = lsp[1].split(',')
			write_sub_aln(fasta, seqs, "%s/%s.%s"%(outfastadir, nust, alnext), keepProtLabels=keepProtLabels)

		ffam.close()
		
else:
	# case where 'famlistdir' is (expected to be) a table file
	print "open", famlistdir
	subfamtab = open(famlistdir, 'r')
	nfam = None
	ffam = None
	nust = None
	fasta = None
	seqs = None
	for line in subfamtab:
		[prot, ust] = line.rstrip('\n').split('\t')
		fam = ust.split('.')[0]
		if fam!=nfam:
			if ffam: ffam.close()
			nffasta = "%s/%s.%s"%(fastadir, fam, alnext)
			if os.access(nffasta, os.F_OK):
				nfam = fam
				n += 1
				sys.stdout.write("\r%d\t\t%s\t\t"%(n, nfam))
				fasta = read_fasta(nffasta)
			else:
				continue
		if ust!=nust:
			if nust: write_sub_aln(fasta, seqs, "%s/%s.%s"%(outfastadir, nust, alnext), keepProtLabels=keepProtLabels)
			nust = ust
			seqs = [prot]
		else:
			seqs.append(prot)
	else:
		if ffam: ffam.close()
		if nust: write_sub_aln(fasta, seqs, "%s/%s.%s"%(outfastadir, nust, alnext), keepProtLabels=keepProtLabels)
		subfamtab.close()
	
sys.stdout.write("\n")
