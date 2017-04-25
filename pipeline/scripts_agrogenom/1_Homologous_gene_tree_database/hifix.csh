#!/bin/csh

if ( $3 == "" ) then
echo "Usage: hifix.csh path/workdir fasta.file blast.outfile"
exit(2)
endif

setenv workdir $1
setenv fasta $2
setenv blastout $3
setenv dataset $fasta:t:r

setenv PATH /panhome/miele/HOGENOMv6/Programs/installdir/bin:${PATH}
# set PYTHONPATH for third-party Python libraries
setenv PYTHONPATH ${PYTHONPATH}:/bge/miele/HOGENOMv6/Programs/argparse-1.2.1:/bge/miele/HOGENOMv6/Programs/installdir/lib/python2.7/site-packages
# set path to MAFFT sequence aligner
setenv MAFFT_BINARIES /bge/miele/HOGENOMv6/Programs/installdir/libexec/mafft/

cd $workdir
# silix clustering
silix $fasta $blastout --net > $dataset.SLX.fnodes
/panhome/miele/HOGENOMv6/Programs/silix-1.2.5-p1/utils/silix-fsize $dataset.SLX.fnodes > $dataset.SLX.fsizes
# hifix family refining
hifix -d /dev/shm $fasta $dataset.net $dataset.SLX.fnodes > $dataset.HFX.fnodes
/panhome/miele/HOGENOMv6/Programs/silix-1.2.5-p1/utils/silix-fsize $dataset.HFX.fnodes > $dataset.HFX.fsizes
