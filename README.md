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

# Proving Non-hyperbolicity.

If we suspect that a group is not hyperbolic (for example the group
F(4,6), defined above, turns out not to be hyperbolic), it can be useful
to spot patterns, that might indicate the presence 
of 'fat triangles', by tracing the geodesic equations which give rise
to new geodesic word differences. This is
achieved by specifying the -ve option. For example

./bin/gpcheckx -geo -ve f46

# Using diagonals to build the correct word acceptor.  
It has been observed that in nearly every case the complete 
set of word differences consists of so called 'diagonal' word 
differences. If two word differences wd1 and wd2 satisfy 
the equation wd2=g1^-1wd1g2, for some generators g1, g2,
then the word g1^-1wd1 is a diagonal of wd1. This leads 
the following procedure for calculating the word acceptor
of such automatic groups.

S0. Run kbprog for a short time to improve the likelihood
that all calculated word differences contained in gpname.diff2
are 'non-spurious'. 

S1. Calculate all possible diagonal words of diff2 and 
add these to make a larger diff2'.

S2. Calculate the word acceptor wa1 based on diff2 GPWA(diff2) , 
and calculate the word acceptor wa2 based on the larger diff2'
GPWA(diff2').

S3 Perform the  fsa operatio wa1 ANDNOT wa2 to create the 
fsa gp.andnot. gp.andnot will recognise lhs words which
are reducible in wa2 using the wd set diff2'but which are
are not reducible in wa1.

S4. Create a list of  lhs words sampled from gpname.andnot.
For each such lhs, calculate rhs=reduced(lhs) using diff2' 
Then calculate the word differences lhs(i)^-1rhs for i ranging from 1 to 
the length(lhs)-1 and any new ones to diff2.

Repeat steps S1 to S4 until wa1 and wa2 are equal.

The options -diagonals s e l  and -diff2name 'diff2suffix' are 
provided to perform the  above procedure. 
-diagonal s e l indicates that diagonals are to be calculated
and added to gpname.diff2'diff2suffix' according to the filter 
s(tart word difference), e(nd word difference and l(imit or maximum number of 
diagonals. Specifying -diagonals 0 0 0 indicates that all
possible diagonals will added to diff2'.

3572 example

./bin/kbprog -wd -t  -me 50000 3572

then repeatedly execute the cycle defined by the 
followin script  


'calculate 3572.wa1 using 3572.diff2'

./bin/gpcheckx -execwa './dowa 3572 >mfile' -waonly  -v  3572
cp 3572.wa 3572.wa1

'calculate all the diagonals of 3572.diff2 to make word difference file' 
'3572.diff2diaggoody'

./bin/gpcheckx  -diagonals 0 0 0 -diff2name diaggoody -w  -v  3572

'calculate 3572.wa2 using 3572.diff2diaggoody (using temporary file 3752,diff2d)'

cp 3572.diff2diaggoody 3572.diff2d
./dowa 3572 >mfile
rm 3572.diff2d
cp 3572.wa 3572.wa2

'calculate 3752.andnot'   

../../maf/bin/fsaandnot 3572.wa1 3572.wa2 3572.andnot >mfile2

'extract new word differences from the reducible lhs words in 3572.andnot.' 
'Create the lhs=rhs equations using 3572.wa1 and 3572.diff2diagody as the' 
'current word acceptor and word difference set respectively '

./bin/gpcheckx -t -to 500 -diff2name diaggoody -v 3572

where dowa is

if test -f $1.diff2d; then
	cp $1.diff2d $1.diff1c
else
	cp $1.diff2 $1.diff1c
fi
../../maf/bin/gpwa $1
cp $1.pwa $1.wa
if test -f $1.diff2c; then
	cp $1.diff2c $1.diff1c
fi


