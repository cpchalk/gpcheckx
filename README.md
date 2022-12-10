# gpcheckx
kbmag program to find word differences quickly using ideas from MAF (https://sourceforge.net/projects/maffsa/).

gpcheckx is a binary that extends the set of binaries making up the kbmag system 
from Derek Holt at the University of Warwick. 
It provides the equivalent functionality of the gpcheckmult component
of gpmakefsa and the part of gpgeowa which calculates the word 
differences needed to build the correct geodesic word acceptor 
('filename'.geowa). It's aim is to find word differences as efficiently
as possible so that the resultant word difference file ('file name'.diff2) builds 
the correct word acceptor ('file name'.wa) and so
the correct general multiplier ('file name'.gm).
 It consists of one source file gpcheckx.c.

To build, edit the kbmag src/makefile to include
gpcheckx in the list of binaries to be built, using, for example,
the build specification for gpcheckmult as a template.  

 It requires the presence of the word difference file ('file name'.diff2)
which has been built by the KBMAG binary kbprog using 
the -wd option.

Example of use. The file f29 defines 
the Fibonacci group F(2,9) and can be
found in the KBMAG test library.

./bin/kbprog -v -wd -me 65000 -t 1000 f29
(creates an initial file f29.diff2)

./bin/gpcheckx -v -p f29 +rptz

(f29.wa and f29.diff2 will be correct)

./bin/gpcheckx -v -geo f29 +rptz

( f29.geowa will be correct)

The difficult examples 3572 and h93.

The main motivation to write gpcheckx was
to enable kbmag to compute the automaticity 
and hyperbolicity of the groups 3572 and h93 within a matter
of hours rather than days or even weeks!

These groups are defined by text files 3572, h93 containing

_RWS := rec(
  isRWS := true,
  ordering := “shortlex”,
  generatorOrder := [a,A,b,B],
  inverses := [A,a,B,b],
 equations := [[a^3,IdWord],[b^5,IdWord],[(a\*b)^7,IdWord],[(a\*b\*A\*B)^2,IdWord]]
);

and 

_RWS := rec(
  isRWS := true,
  ordering := “shortlex”,
  generatorOrder := [a,A,b,B,c,C,d,D,e,E,f,F,g,G,h,H,i,I],
  inverses := [A,a,B,b,C,c,D,d,E,e,F,f,G,g,H,h,I,i],
  equations := [[a\*d,b], [b\*e,c],[c\*f,d],[d\*g,e],[e\*h,f],[f\*i,g],[g\*a,h], 
[h\*b,i], [i\*c,a]]
);

respectively.

The correct word acceptor, diff2 and geowa files can
then be computed as follows.

3572

./bin/kbprog -wd -t 1000 - me 200000 3572

./bin/gpcheckx -v -p -s ‘60000;10000’ 3572

./bin/gpcheckx -p -v 3572 +rptz

./bin/gpcheckx -geo -v 3572 +rptz


h93

./bin/kbprog -wd -t 1500 -me 500000 h93

./bin/gpcheckx -nf  -v  h93

./bin/gpcheckx -p -s 40000  -to 600 -v  h93

./bin/gpcheckx -p -s 60000  -to 600 -v  h93

./bin/gpcheckx -m -s 110000  -to 600 -v  h93

./bin/gpcheckx -m -s 150000  -to 900 -v  h93

./bin/gpcheckx -m -s 170000  -to 600 -v  h93

./bin/gpcheckx -p  -v h93 +rptz

./bin/gpcheckx -geo -v h93 +rptz

Another example of an automatic group which kbmag 
finds difficult to compute is the Fibonacci group
F(4,6) defined by the file f46 with contents:-

_RWS := rec(
  isRWS := true,
  ordering := “shortlex”,
  generatorOrder := [a,A,b,B,c,C,d,D,e,E,f,F],
  inverses := [A,a,B,b,C,c,D,d,E,e,F,f],
  equations := 
[[a\*b\*c\*d,e], [b\*c\*d\*e,f],
 [c\*d\*e\*f,a],[d\*e\*f\*a,b],
 [e\*f\*a\*b,c],[f\*a\*b\*c,d]]
);

The correct diff2 and wa files can then be
calculated as follows.

./bin/kbprog -wd -t 1000 -me 200000 f46

./bin/gpcheckx -v -p -s ‘4000;1000’ f46

./bin/gpcheckx -v -p f46 +rptz

Proving Non-hyperbolicity.

If we suspect that a group is not hyperbolic (for example the group
F(4,6), defined above, turns out not to be hyperbolic), it can be useful
to spot patterns, that might indicate the presence 
of 'fat triangles', by tracing the geodesic equations which give rise
to new geodesic word differences. This is
achieved by specifying the -ve option. For example

./bin/gpcheckx -geo -ve f46
