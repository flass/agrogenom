#!/bin/csh
#########obsolete/usr/local/bin/csh
echo "JOB NAME : $PBS_JOBNAME"
echo "JOB ID :   $PBS_JOBID"
set i = 0
set seuil = 75		# % de cpu avant arret

set tachelist = $1
set results = $2

echo "LISTE DE TACHES   : $tachelist"
echo "RESULTATS         : $results"
echo "WORK DIR		: $PWD"


# specify the path to task manager program getnextbig
setenv utils /path/to/utils


while ($i < 100000)



@ i = $i + 1



# Traitement
# ----------
echo ""
set taskdate = `date`
echo "Task initiated at $taskdate"

# - Recuperation de la tache a effectuer
#   -------------------------------------
set gettache = `$utils/getnextbig $tachelist find `


if ( "$gettache" == '***fini' ) then
	echo "======= All tasks have been processed  ===="
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

set file = $tache:t:r
set ext = $tache:e

echo "Tache numero $numtache  -  $tache  - file $file.$ext"

set nseq = `grep ">" $tache | wc -l`
if ( $nseq < 3 ) then
	echo "$nseq sequences, pas assez pour calculer un arbre"
	$utils/getnextbig $tachelist finish $tache
	echo "Saute a la prochaine tache."
else

# - Copie des donnees
#   -------------------------------------

echo "Copie des donnees."
if ( -e $tache ) then
echo copy $tache
\cp $tache . 
set stexec = $status
else
echo copy $tache
\cp $tache . 
set stexec = $status
endif

echo Statut : $stexec
if ($stexec > 0 ) then
	echo "Erreur pendant la copie des donnees"
	$utils/getnextbig $tachelist error $tache
	echo "Finit le job ici"
	break
endif



# - Execution
#   ---------

echo "Calcul de l'alignement Muscle."
/usr/bin/muscle -in $file.$ext -out $file.aln -quiet
set stexec = $status
echo Statut : $stexec
if ($stexec > 0 ) then
	echo "Erreur pendant muscle"
	$utils/getnextbig $tachelist error $tache
	echo "Finit le job ici"
	break
endif


echo "Copie des resultats"
\cp ./$file.aln $results/$file.aln 
set stexec = $status
echo Statut : $stexec
if ($stexec > 0 ) then
	echo "Erreur pendant la copie des resultats"
	$utils/getnextbig $tachelist error $tache
	echo "Finit le job ici"
	break
endif

echo "Tache executee avec  succees, saute a la prochaine tache."
$utils/getnextbig $tachelist finish $tache

endif

end

echo "fin de boucle"

