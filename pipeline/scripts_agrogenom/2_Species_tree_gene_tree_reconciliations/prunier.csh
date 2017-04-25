#!/bin/csh
echo "JOB NAME : $PBS_JOBNAME"
echo "JOB ID :   $PBS_JOBID"

set i = 0
set ier = 0

set tachelist = $1
set results = $2
set reftree = $3
set moltype = $4
set genetreedir = $5
setenv maxh $6
set bsc = $7
set bs = $8

# specify the path to task manager program getnextbig
setenv utils /path/to/utils
# specify the path to Prunier binary
setenv prunier_bin /path/to/prunier


set maxm = 0
@ maxm = $maxh * 60

echo "LISTE DE TACHES (ALIGNEMENTS DE GENES)   : $tachelist"
echo "LISTE D ARBRES DE GENES  : $genetreedir"
echo "ARBRE DE REFERENCE   : $reftree"
echo "RESULTATS         : $results"

# Calcul de la date

set dminute = `date +%M`
set dhour = `date +%H`
set dday = `date +%d`

set run = 1


while ($run != 0)

@ i = $i + 1


# Traitement
# ----------
echo ""
echo "Le traitement commence..."

# - Recuperation de la tache a effectuer
#   -------------------------------------
set gettache = `$utils/getnextbig $tachelist find `
echo $gettache

if ( "$gettache" == '***fini' ) then
	echo "======= All tasks have been processed  ===="
	set run = 0
	@ i = $i - 1
	date
	echo "Number of tasks  processed          : $i"
	exit(0)
endif

set numtache = `python -c "s='$gettache' ; print s.split()[-2]"`
set tache = `python -c "s='$gettache' ; print s.split()[-1]"`
set control = `python -c "s='$tache' ; print len(s.split('/'))"`
if ( $control < 2 ) then
	echo "scheduller bug yield improper task: $tache ; skip it"
	continue
endif

echo $tache:t

set fam = `python -c "s='$tache' ; print s.rsplit('/',1)[1].rsplit('.',1)[0]"`
set alnext = `python -c "s='$tache' ; print s.rsplit('/',1)[1].rsplit('.',1)[1]"`

if ( $genetreedir != "none" ) then
	set treefile = $genetreedir/$fam.phb
	if (! -e $treefile ) then	
		echo "pas d'arbre de gene"
		$utils/getnextbig $tachelist error $tache
		@ ier = $ier + 1
		echo "Arrete le job ici"
		break
		#~ echo "Fichier d'arbre $treefile manquant"
		#~ set treefile = "none"
		#~ echo "Prunier (TreeFinder) will compute it"
		#~ set bs = 100
	endif
else
	set treefile = "none"
endif

echo "Tache numero $numtache  -  $tache  - aln.file $fam.$alnext"
echo "Copie des donnees."

if ( -e $tache ) then
	echo "copy $tache"
	\cp -f $tache . 
	set stexec1 = $status
endif

echo "copy $reftree"
\cp -f $reftree . 
set stexec2 = $status
endif

set ref = $reftree:t

if ( $treefile != "none" ) then
	echo "copy $treefile"
	\cp -f $treefile . 
	set stexec3 = $status
else
	set stexec3 = 0
endif

@ stexec = $stexec1 + $stexec2 + $stexec3
echo "Statut : $stexec = $stexec1 + $stexec2 + $stexec3"
if ($stexec > 0 ) then
	echo "Erreur pendant la copie des donnees"
	$utils/getnextbig $tachelist error $tache
	@ ier = $ier + 1
	echo "Arrete le job ici"
	break
#	echo "Saute a la prochaine tache"
#	continue
endif

set gt = $treefile:t



# - Execution
#   ---------
echo "ls $PWD"
ls
echo "Recherche des transfers par PrunierFWD."
set pruniercall = "$prunier_bin/PrunierFWD_linux64 input.tree.file=$ref aln.file=$fam.$alnext genetree.file=$gt sequence.type=$moltype fwd.depth=2 aln.type=FASTA boot.thresh.conflict=$bsc max.bp=$bs multi_root.bool=false"
echo "Call: $pruniercall > $fam.prout"
$pruniercall > $fam.prout
set stexec = $status
echo "Statut : $stexec"
if ($stexec > 0 ) then
	echo "Erreur pendant PrunierFWD"
	$utils/getnextbig $tachelist error $tache
	@ ier = $ier + 1
	echo "Arrete le job ici"
	break
#	echo "Saute a la prochaine tache"
#	continue
endif


echo "Copie des resultats"
\cp ./$fam.prout $results/$fam.prout
set stexec1 = $status 
if ( $treefile == "none") then
	\cp ./$fam.$alnext\_tf.stats $results/$fam.$alnext\_tf.stats
	set stexec2a = $status 
	\cp ./$fam.$alnext\_tf.tree $results/$fam.$alnext\_tf.tree
	set stexec2b = $status
	@ stexec2 = $stexec2a + $stexec2b
else 
	set stexec2 = 0 
endif
#\cp ./$fam.$alnext.ref_tree_all_roots $results/$fam.$alnext.ref_tree_all_roots

@ stexec = $stexec1 + $stexec2
echo "Statut : $stexec = $stexec1 + $stexec2"
if ($stexec > 0 ) then
	echo "Erreur pendant la copie des resultats"
	$utils/getnextbig $tachelist error $tache
#	echo "Saute a la prochaine tache"
#	continue
endif

echo "Tache executee avec  succes, saute a la prochaine tache "
$utils/getnextbig $tachelist finish $tache


# Calcul de la date

set cminute = `date +%M`
set chour = `date +%H`
set cday = `date +%d`

set sday = 0
set shour = 0
set sminute = 0
set elaps = 0
@ sday = $cday - $dday
@ shour = $chour - $dhour
@ sminute = $cminute - $dminute

@ elaps = $sminute + ($shour * 60 )  + ($sday * 60 * 24 )

echo 
echo "Elapsed time : $elaps minutes / $maxm "

if ($elaps > (($maxm * 80) / 100) ) then
	echo "Job reached 80% of  calculation time. Job will stop now"
	set run = 0
	date
	echo "Number of tasks  processed          : $i"
	echo "Number of tasks  exiting with error : $ier"
	exit(0)

endif
 
#getnextbig $task finish $tache

date
echo "------------- End Task Number $i --------------"

end
echo "fin de boucle"

