#! /bin/tcsh
if ( $2 == "" ) then
echo "Usage: ./concat.sh listjkfam_dir/ outjkaln_dir/"
else

set jkfam=$1
set jkaln=$2

foreach famlist ( $jkfam/* )

set jk=$famlist:t
echo $jk

python /pandata/lassalle/agrogenom/scripts/concat.py $famlist $jkaln/$jk.aln
end
endif
