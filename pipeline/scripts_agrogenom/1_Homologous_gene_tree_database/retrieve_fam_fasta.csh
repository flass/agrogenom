#!/bin/csh

if ( $4 == "" ) then
	echo "Usage: path/acnucdb path/outfamlist path/outfastadir datatype"
	exit(2)
endif
setenv database $1
setenv outfamlist $2
setenv outfastadir $3
setenv datatype $4

if ( $datatype == 'aa' ) then
	setenv dt 2
else if ( $datatype == 'nuc' ) then
	setenv dt 1
else
	echo "Bad data type $datatype; must be 'aa' or 'nuc'"
	exit(2)
endif

setenv gcgacnuc $database/flat_files
setenv acnuc $database/index

\rm $outfamlist

/panhome/banques/debian-bin/query << EOF
terminal
n
select
kd(nk=gene_family)
save
list1
$outfamlist
stop
EOF

\rm $outfastadir/*

foreach fam ( `cat  $outfamlist` )
echo $fam
if ( $fam != 'GENE_FAMILY' ) then
/panhome/banques/debian-bin/query << EOF
terminal
n
select
k=$fam
extract
list1
y
3
$outfastadir/$fam.fasta
$dt
stop
EOF
endif
end

