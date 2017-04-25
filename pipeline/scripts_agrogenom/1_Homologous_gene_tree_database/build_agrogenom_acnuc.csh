#!/bin/csh

#if ( $1=="" ) then
#echo "Usage: build_agrogenom_acnuc.csh path/to/database_dir"
#exit(2)
#endif
setenv database $1

#if ( ! -d $database/flat_files ) then
#echo "'$database/' should contain a directory 'flat_files/' with genome flat files in EMBL format"
#exit(2)
#endif
setenv gcgacnuc $database/flat_files

#if ( ! -d $database/flat_files ) then
#mkdir $database/flat_files
#endif
setenv acnuc $database/index

# path to ACNUC binaries
# can be dowloaded from http://doua.prabi.fr/databases/acnuc in `Local installation' section.
setenv bins /path/to/acnuc/bin/
setenv log $database/agrogenom_build.log
setenv date `date`
echo "building agrogenom ACNUC database on $date" >> $log
echo "flat files at : $gcgacnuc" >> $log
echo "indexes at : $acnuc" >> $log
echo "" >> $log

cd $acnuc
$bins/initf embl
echo "initf status: $status" >> $log

cd $database
foreach dat (`ls $gcgacnuc`)
set name=`echo $dat | cut -d'.' -f1`
$bins/acnucgener d $name
echo "acnucgener status for $dat : $status" >> $log
end

$bins/ordnet K
echo "ordnet K status : $status" >> $log
$bins/ordnet S
echo "ordnet S status : $status" >> $log

$bins/newordalphab
echo "newordalphab status : $status" >> $log
