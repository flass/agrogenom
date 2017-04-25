#!/bin/csh

echo "JOB NAME : $PBS_JOBNAME"
echo "JOB ID :   $PBS_JOBID"
echo "WORK DIR		: $PWD"
echo "call: phyml.csh $1 $2 $3 $4 $5 $6"

cd /pandata/lassalle/agrogenom/big_trees

set curdir = $1
set tache = $2
set gbresults = $3
set treeresults = $4

setenv maxh $5

# specify the path to task manager program getnextbig
setenv utils /path/to/utils


if ( $6 != "" ) then
	set alnformat = $6
else
	set alnformat = 'fasta'
endif
set taskdate = `date`
echo "Task starts at $taskdate"

set file = $tache:t:r
set ext = $tache:e
echo "Tache   -  $tache  - file $file.$ext"


echo "Copie des donnees."
echo copy $tache
\cp $tache .  || (echo "erreur durant la copie de l'alignement $alnformat" ; exit(1))


if ( $alnformat != 'fasta' ) then
	echo "Conversion $alnformat -> fasta"
	echo "convert_ali $alnformat ./$file.$ext fasta ./$file.aln"
	$utils/convert_ali $alnformat ./$file.$ext fasta ./$file.aln || (echo "erreur durant la conversion au format fasta" ; exit(1))
	set ext = 'aln'
endif


# Verif si le calcul d'arbre est necessaire cad si le nb de seq est > 2
set buf = `grep ">" ./$file.$ext | wc -l`

set nbseqs = $buf

echo "Nombre de sequences = $nbseqs"

if  ($nbseqs <= 2) then
	echo "Nombre de sequences insufisant "
	echo "Tache consideree comme executee avec  succees,saute a la prochaine tache "
	exit(2)
endif

if ( ! -d $gbresults ) then
	echo "Traitement Gblocks"
	$utils/get_options_gb "  -t=d -b1=50 -b2=50 -b5=a  " ./$file.$ext > ./options || (echo "erreur durant le calcul des options gblocks" ; exit(1))
	set options = `cat ./options`
	echo "Options Gblocks: $options"
	echo Gblocks $file.$ext $options
	Gblocks $file.$ext $options

	echo "Conversion format phylip"
	echo "convert_phylip $file.$ext-gb"
	$utils/convert_phylip $file.$ext-gb || (echo "erreur durant la conversion au format phylip" ; exit(1))
	set alnext = "$ext-gb.phyl"

	echo "Copie de l'alignement traite"
	ls
	echo "\cp ./$file.$alnext $gbresults/"
	\cp ./$file.$alnext $gbresults/ || (echo " erreur durant la copie de l'alignement trime" ; exit(1))
else
	echo "No GBlocks trimming"
	echo "Conversion format phylip"
	echo "convert_phylip $file.$ext-gb"
	$utils/convert_phylip $file.$ext || (echo "erreur durant la conversion au format phylip" ; exit(1))
	set alnext = "$ext.phyl"
endif

echo "Calcul de l'arbre"
set phymlopt = "-i $file.$alnext -d nt -q -b -3  -m GTR -v e -c 8 -a e -s BEST --print_trace"
echo "phyml $phymlopt"
/usr/bin/phyml $phymlopt || (echo " erreur durant le calcul de l'arbre" ; exit(1))

echo "Copie de l'arbre"
ls
echo "\cp ./$file.$ext-gb.phyl_phyml_tree.txt $treeresults/$file.phb"
\cp ./$file.$ext-gb.phyl_phyml_tree.txt $treeresults/$file.phb || (echo " erreur durant la copie de l'arbre" ; exit(1))


