#!/bin/csh
setenv tache $1
setenv alns $2
setenv results $3
setenv reftree $4
setenv maxh $5

setenv penel /panhome/penel/penel/svn/bin/Debian/
@ maxm = $maxh * 60

# Calcul de la date

set dminute = `date +%M`
set dhour = `date +%H`
set dday = `date +%d`



echo "python scripts/find_ancestral_duplications.py u $tache $alns $results $reftree 0.9 0.2"
python /pandata/lassalle/agrogenom/scripts/find_ancestral_duplications.py u $tache $alns $results $reftree 0.9 0.2

set stexec = $status
if ( $stexec > 0 ) then
	echo "erreur durant le script find_ancestral_duplications.py"
endif
