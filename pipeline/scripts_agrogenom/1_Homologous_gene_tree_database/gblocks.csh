#!/bin/csh

if ( $3 == "" ) then
	echo "Usage: gblocks_codons.csh path/codon_aln_list path/results gbocks_mode(nuc|prot|codon)"
	exit(2)
endif


setenv penel /panhome/penel/penel/svn/bin/Debian/

set tasklist = $1
set results = $2
set mode = $3

if ( $mode == "nuc" ) then
	set t = "d"
else if ( $mode == "prot" ) then
		set t = "p"
	else if ( $mode == "codon" ) then
			set t = "c"
		endif
	endif
endif

foreach aln ( `cat $tasklist` )
	$penel/get_options_gb "  -t=$t -b1=50 -b2=50 -b5=a  " $aln > ./gboptions
	set options = `cat ./gboptions`
	Gblocks $aln $options
	\mv $aln-gb $results
	endif
end

\rm ./gboptions
