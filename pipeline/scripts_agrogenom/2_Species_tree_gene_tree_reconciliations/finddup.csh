#!/bin/csh
setenv tachelist $1
setenv alns $2
setenv results $3
setenv reftree $4
setenv maxh $5

# specify the path to task manager program getnextbig
setenv utils /path/to/utils
# specify the path to find_ancestral_duplications.py script
setenv scripts /path/to/script
@ maxm = $maxh * 60

# Calcul de la date

set dminute = `date +%M`
set dhour = `date +%H`
set dday = `date +%d`

set i = 1

while ( $i > 0 )
echo ""

set gettache = `$utils/getnextbig $tachelist find `

if ( "$gettache" == '***fini' ) then
	echo "======= All tasks have been processed  ===="
	@ i = $i - 1
	date
	break
endif

echo $gettache

set numtache = `python -c "s='$gettache' ; print s.split()[-2]"`
set tache = `python -c "s='$gettache' ; print s.split()[-1]"`
set control = `python -c "s='$tache' ; print len(s.split('/'))"`
if ( $control < 2 ) then
	echo "scheduller bug yield improper task: $tache ; skip it"
	continue
endif

echo "python $scripts/find_ancestral_duplications.py u $tache $alns $results $reftree 0.9 0.2"
python $scripts/find_ancestral_duplications.py u $tache $alns $results $reftree 0.9 0.2

set stexec = $status
if ( $stexec > 0 ) then
	$utils/getnextbig $tachelist error $tache
	echo "erreur durant le script find_ancestral_duplications.py"
	break
else
	$utils/getnextbig $tachelist finish $tache
endif
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
if ($elaps > (($maxm * 85) / 100) ) then
	echo "Job reached 85% of  calculation time. Job will stop now"
	@ i = $i - 1
	date
	echo "Number of tasks  processed          : $i"
	echo "Number of tasks  exiting with error : $ier"
	exit(0)

endif

end
