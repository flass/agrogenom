#!/bin/csh
#########obsolete/usr/local/bin/csh
echo "JOB NAME : $PBS_JOBNAME"
echo "JOB ID :   $PBS_JOBID"
set i = 0
set ier = 0
set seuil = 75		# % de cpu avant arret
set curdir = $1
set tachelist = $2
set results = $3

setenv maxh $4

# specify the path to task manager program getnextbig
setenv utils /path/to/utils


#set fasta = 0
#if ($5 == "fasta") then
#	set fasta = 1
#endif



echo "REPERTOIRE COURANT    : $curdir"
echo "LISTE DE TACHES       : $tachelist"
echo "RESULTATS             : $results"
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


set numtache = $gettache[1]
set tache = $gettache[2]


echo "------------- Start Task Number $i of job --------------"
echo DEBUG $gettache
echo $tache | awk -F/ '{print $NF}'

set file = $tache:t:r
set ext = $tache:e
echo "Tache numero $numtache  -  $tache  - file $file.$ext"
echo 



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


#if ( $fasta == 1) then
## - Execution
##   ---------
#
#echo "Conversion fasta -> clustal"
#echo convert_ali fasta ./$file.$ext clustal ./$file.clustal
#$utils/convert_ali fasta ./$file.$ext clustal ./$file.clustal
#set stexec = $status
#echo Statut : $stexec
#if ($stexec > 0 ) then
#	echo "Erreur pendant la conversion"
#	$utils/getnextbig $curdir/$tachelist error $tache
#	echo "Saute a la prochaine tache"
#	@ ier = $ier + 1
#	break
#endif
#
#
#echo "Copie de l'alignement au format clustalw"
#echo \cp ./$file.clustal $results/$file.aln 
#\cp ./$file.clustal $results/$file.aln 
#set stexec = $status
#echo Statut : $stexec
#if ($stexec > 0 ) then
#	echo "Erreur pendant la copie"
#	$utils/getnextbig $curdir/$tachelist error $tache
#	echo "Saute a la prochaine tache"
#	@ ier = $ier + 1
#	break
#endif
#
#
#else
#
#echo \cp ./$file.$ext $results/$file.aln 
#\cp ./$file.$ext $results/$file.aln 
#
#echo "Conversion phylip -> fasta"
#echo convert_ali phylip ./$file.aln fasta ./$file.fasta
#$utils/convert_ali phylip ./$file.aln fasta ./$file.fasta
#set stexec = $status
#echo Statut : $stexec
#if ($stexec > 0 ) then
#	echo "Erreur pendant la conversion"
#	$utils/getnextbig $curdir/$tachelist error $tache
#	echo "Saute a la prochaine tache"
#	@ ier = $ier + 1
#	break
#endif
#
#echo "Copie de l'alignement au format fasta"
#echo \cp ./$file.fasta ./$file.aln 
#\cp ./$file.fasta ./$file.aln 
#set stexec = $status
#echo Statut : $stexec
#if ($stexec > 0 ) then
#	echo "Erreur pendant la copie"
#	$utils/getnextbig $curdir/$tachelist error $tache
#	echo "Saute a la prochaine tache"
#	@ ier = $ier + 1
#	break
#endif
#
#
#endif



# Verif si le calcul d'arbre est necessaire cad si le nb de seq est > 2
set buf = `grep ">" ./$file.aln |wc`

set nbseqs = $buf[1]

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

echo "Traitement Gblocks"
$utils/get_options_gb "  -t=d -b1=50 -b2=50 -b5=a  " ./$file.aln > ./options
set stexec = $status
echo Statut : $stexec
if ($stexec > 0 ) then
	echo "Erreur pendant le calcul des options gblocks"
	$utils/getnextbig $curdir/$tachelist error $tache
	echo "Saute a la prochaine tache"
	@ ier = $ier + 1
	break
endif
set options = `cat ./options`
echo "Options Gblocks: $options"
echo Gblocks $file.aln $options
Gblocks $file.aln $options
set stexec = $status



echo "Conversion format phylip"
echo "convert_phylip $file.aln-gb"
$utils/convert_phylip $file.aln-gb
set stexec = $status
echo Statut : $stexec

echo "Calculate ML tree"
set raxmlopt = "-s $file.aln-gb.phyl -n $file.tmp_tree -m GTRGAMMA -f d"
echo "raxml $raxmlopt"
/panhome/lassalle/RAxML-7.2.8-ALPHA/raxmlHPC-SSE3 $raxmlopt

set stexec = $status
echo Statut : $stexec
if ($stexec > 0 ) then
	echo "Erreur pendant phyml"
	$utils/getnextbig $curdir/$tachelist error $tache
	echo "Saute a la prochaine tache"
	@ ier = $ier + 1
	break
endif

echo "Calculate SH-like support on ML tree"
set raxmlopt = "-s $file.aln-gb.phyl -n $file.shlike_tree -m GTRGAMMA -f J -t $file.phb"
echo "raxml $raxmlopt"
/panhome/lassalle/RAxML-7.2.8-ALPHA/raxmlHPC-SSE3 $raxmlopt

set stexec = $status
echo Statut : $stexec
if ($stexec > 0 ) then
	echo "Erreur pendant phyml"
	$utils/getnextbig $curdir/$tachelist error $tache
	echo "Saute a la prochaine tache"
	@ ier = $ier + 1
	break
endif

echo "Copie des resultats"
ls
echo "\cp ./$file.shlike_tree* $results/"
\cp ./$file.shlike_tree* $results/
set stexec = $status
echo Statut : $stexec
if ($stexec > 0 ) then
	echo "Erreur pendant la copie des resultats"
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

