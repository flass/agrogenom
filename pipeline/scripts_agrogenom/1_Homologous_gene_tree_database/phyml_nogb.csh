#!/bin/csh

echo "JOB NAME : $PBS_JOBNAME"
echo "JOB ID :   $PBS_JOBID"
echo "WORK DIR		: $PWD"
echo "call: phyml.csh $1 $2 $3 $4 $5 $6"

set i = 0
set ier = 0
set seuil = 75		# % de cpu avant arret
set curdir = $1
set tachelist = $2
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



echo "REPERTOIRE COURANT    : $curdir"
echo "LISTE DE TACHES       : $tachelist"
echo "RESULTATS GBLOCKS     : $gbresults"
echo "RESULTATS ARBRE       : $treeresults"
echo "MAX. NUMBER OF HOURS  : $maxh"

# Calcul de la date

set dminute = `date +%M`
set dhour = `date +%H`
set dday = `date +%d`


set i = 1

@ maxm = $maxh * 60


#while ($i < 1000000)



#@ i = $i + 1


while ( $i >= 0 )


# Traitement
# ----------

set taskdate = `date`
echo "Task starts at $taskdate"

# - Recuperation de la tache a effectuer
#   -------------------------------------
set gettache = `$utils/getnextbig $curdir/$tachelist find `

if ( "$gettache" == '***fini' ) then
	echo "======= All tasks have been processed  ===="
	@ i = $i - 1
	date
	echo "Number of tasks  processed          : $i"
	echo "Number of tasks  exiting with error : $ier"
	exit(0)
endif


set numtache = `python -c "s='$gettache' ; print s.split()[-2]"`
set tache = `python -c "s='$gettache' ; print s.split()[-1]"`
set control = `python -c "s='$tache' ; print len(s.split('/'))"`
if ( $control < 2 ) then
	echo "scheduller bug yield improper task: $tache ; skip it"
	continue
endif

echo "------------- Start Task Number $i of job --------------"
echo DEBUG $gettache
echo $tache | awk -F/ '{print $NF}'

set file = $tache:t:r
set ext = $tache:e
echo "Tache numero $numtache  -  $tache  - file $file.$ext"
echo 

# check for pre-existence of gene tree
if ( -e $treeresults/$file.phb) then
	" tree $treeresults/$file.phb already computed ; skip task"
	$utils/getnextbig $curdir/$tachelist finish $tache
	continue
endif

echo "Copie des donnees."
if ( -e $tache ) then
	echo copy $tache
	\cp $tache . 
	set stexec = $status
else
	echo copy $curdir/$tache
	\cp $curdir/$tache . 
	set stexec = $status
endif

echo Statut : $stexec

if ($stexec > 0 ) then
	echo "Erreur pendant la copie des donnees"
	$utils/getnextbig $curdir/$tachelist error $tache
	echo "Saute a la prochaine tache"
	@ ier = $ier + 1
	break
endif

if ( $alnformat != 'fasta' ) then
	echo "Conversion $alnformat -> fasta"
	echo "convert_ali $alnformat ./$file.$ext fasta ./$file.aln"
	$utils/convert_ali $alnformat ./$file.$ext fasta ./$file.aln
	set stexec = $status
	echo Statut : $stexec
	if ($stexec > 0 ) then
		echo "Erreur pendant la conversion"
		$utils/getnextbig $curdir/$tachelist error $tache
		echo "Saute a la prochaine tache"
		@ ier = $ier + 1
		break
	endif
	set ext = 'aln'
endif


# Verif si le calcul d'arbre est necessaire cad si le nb de seq est > 2
set buf = `grep ">" ./$file.$ext | wc -l`

set nbseqs = $buf

echo "Nombre de sequences = $nbseqs"

if  ($nbseqs <= 2) then

echo "Nombre de sequences insufisant "
echo "Tache consideree comme executee avec  succees,saute a la prochaine tache "
$utils/getnextbig $curdir/$tachelist finish $tache
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

if ($elaps > (($maxm * 85) / 100) ) then
	echo "Job reached 85% of  calculation time. Job will stop now"
	@ i = $i - 1
	date
	echo "Number of tasks  processed          : $i"
	echo "Number of tasks  exiting with error : $ier"
	exit(0)

endif
 

date
echo "------------- End Task Number $i --------------"

@ i = $i + 1

#break
#endif

else


echo "Conversion format phylip"
echo "convert_phylip $file.$ext"
$utils/convert_phylip $file.$ext
set stexec = $status
echo Statut : $stexec

echo "Copie de l'alignement traite"
ls
echo "\cp ./$file.$ext.phyl $gbresults/"
\cp ./$file.$ext.phyl $gbresults/
set stexec = $status
echo Statut : $stexec
if ($stexec > 0 ) then
	echo "Erreur pendant la copie de l'alignement trime"
	$utils/getnextbig $curdir/$tachelist error $tache
	echo "Saute a la prochaine tache"
	@ ier = $ier + 1
	break
endif


set phymlopt = "-i $file.$ext.phyl -d nt -q -b -3  -m GTR -v e -c 8 -a e -s BEST"
echo "phyml $phymlopt"
/usr/bin/phyml $phymlopt

set stexec = $status
echo Statut : $stexec
if ($stexec > 0 ) then
	echo "Erreur pendant phyml"
	$utils/getnextbig $curdir/$tachelist error $tache
	echo "Saute a la prochaine tache"
	@ ier = $ier + 1
	break
endif


echo "Copie de l'arbre"
ls
echo "\cp ./$file.$ext.phyl_phyml_tree.txt $treeresults/$file.phb"
\cp ./$file.$ext.phyl_phyml_tree.txt $treeresults/$file.phb
set stexec = $status
echo Statut : $stexec
if ($stexec > 0 ) then
	echo "Erreur pendant la copie de l'arbre"
	$utils/getnextbig $curdir/$tachelist error $tache
	echo "Saute a la prochaine tache"
	@ ier = $ier + 1
	break
endif

echo "Tache executee avec  succees,saute a la prochaine tache "
$utils/getnextbig $curdir/$tachelist finish $tache

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

if ($elaps > (($maxm * 85) / 100) ) then
	echo "Job reached 85% of  calculation time. Job will stop now"
	@ i = $i - 1
	date
	echo "Number of tasks  processed          : $i"
	echo "Number of tasks  exiting with error : $ier"
	exit(0)

endif
 
#getnextbig $task finish $tache
endif
date
echo "------------- End Task Number $i --------------"

@ i = $i + 1

end
echo fin de boucle

