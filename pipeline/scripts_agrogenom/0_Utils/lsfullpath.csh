#!/bin/csh

if ( $1 != "" ) then
setenv dir $1
else
setenv dir $PWD
endif

foreach file ( $dir/*$2 )
echo $file
end
