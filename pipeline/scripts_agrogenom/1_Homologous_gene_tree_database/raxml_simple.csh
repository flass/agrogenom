#!/bin/csh
#########obsolete/usr/local/bin/csh
echo "JOB NAME : $PBS_JOBNAME"
echo "JOB ID :   $PBS_JOBID"

set tache = $1
set results = $2
set alnformat = $3
set jnum = $4
set raxmlopt = `cat $5`

# specify the path to task manager program getnextbig
setenv utils /path/to/utils

echo "TACHE     			: $tache"
echo "RESULTATS             : $results"
echo "ALIGNMENT FORMAT		: $alnformat"
echo "RAxML OPTIONS			: $raxmlopt"

# Traitement
# ----------

set taskdate = `date`
echo "Task starts at $taskdate"

set fil = $tache:t:r
set file = $fil.$jnum
set ext = $tache:e
echo "Tache $tache  - file $file.$ext"
echo 

echo "Copie des donnees."
echo "copy $tache ./$file.$ext"
\cp $tache ./$file.$ext || (echo "error when copying alignment" ; exit(1))

if ( $alnformat != 'phylip' )then
	echo "Conversion format phylip"
	echo "convert_phylip $file.$ext"
	$utils/convert_phylip $file.$ext || (echo "error when converting alignment to phylip" ; exit(1))
	set ext = "$ext.phyl"
endif


echo "Calculate ML tree (with rapid bootstrapping)"
set raxmlcmd = "/panhome/lassalle/RAxML-7.2.8-ALPHA/raxmlHPC-SSE3 -s $file.$ext -n $file $raxmlopt"
echo $raxmlcmd
$raxmlcmd || (echo "error during RAxML" ; exit(1))


echo "Copie des resultats"
ls
echo "\cp ./$file* $results/"
\cp ./RAxML* $results/ || (echo "error when copying result files" ; exit(1))
