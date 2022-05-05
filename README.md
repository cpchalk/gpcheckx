# gpcheckx
kbmag program to find word differences quickly using ideas from MAF
gpcheckx is a binary that extends the set of binaries making up the kbmag system 
from Derek Holt at the University of Warwick. 
It can be used as a replacement for the gpcheckmult component
of gpmakefsa and the part of gpgeowa which calculates the word 
diffferences needed to build the correct geodesic word acceptor 
(''filename'.geowa'. It's aim is to find word differences as efficiently
as possible so that the resultant word difference file builds 
the correct word acceptor ('file name'.wa) and 
and general multiplier ('file name.gm).
 It consists of one source file gpcheckx.c.
 It requires an initial word difference file ('file name'.diff2)
which has been built by the KBMAG binary kbprog using 
the -wd option.
Example of use.
./bin/kbprog -v -wd -me 65000 -t 1000 f29
./bin/gpcheckx -v -p f29 +rptz
# f29.wa and f29.diff2 will be correct
./bin/gpcheckx -v -geo f29 +rptz
# f29.geowa will be correct
