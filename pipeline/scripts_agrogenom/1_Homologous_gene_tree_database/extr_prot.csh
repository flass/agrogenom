#!/bin/csh

if ( $2 == "" ) then
	echo "Usage: extr_prot.csh path/acnuc_database path/output_fasta_file [datatype]"
	exit(2)
endif

setenv database $1
setenv outfasta $2
setenv datatype $3

setenv acnuc $database/index/
setenv gcgacnuc $database/flat_files/

setenv outdir $outfasta:h
if !( -d $outdir ) then
	mkdir $outdir
endif
cd $outdir

if ( -e $outfasta ) then
	\rm $outfasta
endif
setenv fasta $outfasta:r

if ( $datatype == 'aa' ) then
	setenv dt 2
else 
	if ( $datatype == 'nuc' ) then
		setenv dt 1
	else 
		if ( $datatype == "" ) then
			echo "Default data type is 'aa'"
			setenv dt 2
		else
			echo "Bad data type $datatype; must be 'aa' or 'nuc'"
			exit(2)
		endif
	endif
endif

/panhome/banques/debian-bin/query << EOF
terminal
n
select
t=cds
extract
list1
y
3
$fasta
$dt
stop
EOF
