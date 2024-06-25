# gpcheckx
kbmag program to find word differences quickly using ideas from MAF (https://sourceforge.net/projects/maffsa/).

gpcheckx is a binary that extends the set of binaries making up the kbmag system 
from Derek Holt at the University of Warwick. 
It provides the equivalent functionality of the gpcheckmult component
of gpmakefsa (-p option) and the part of gpgeowa which calculates the word 
differences needed to build the correct geodesic word acceptor 
'filename'.geowa (-geo option). It's aim is to find word differences as efficiently
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

# Using diagonals to build a 'small' and possibly correct word acceptor.  
It has been observed that, in many cases, the complete 
set of word differences of an automatic group consists entirely 
or almost entirely of so called 'diagonal' word differences. 
If two word differences wd1 and wd2 satisfy the equation 
wd2=g1^-1wd1g2, for some generators g1, g2,
then the word g1^-1wd1 is called a diagonal of wd1. 
This observation leads to considering the following 
speculative procedure for extracting new word differences . 


S0. Run kbprog for a short time to improve the likelihood
that all calculated word differences contained in 
gpname.diff2 (or diff2 for short) are 'non-spurious'. 

S1. Calculate all possible diagonal words of diff2 and 
add these to make a larger word difference set diff2'.

S2. Use the executable GPWA to calculate the word acceptor, wa1, 
based on the word diference set diff2.
Then calculate the word acceptor, wa2, based on the larger word 
difference set diff2'

S3  'Compare' the word acceptors wa1 and wa2, by 
performing the fsa operation wa1 ANDNOT wa2 to create the 
fsa gpname.andnot. This fsa will recognise reducible lhs words 
which fail to be recognised as such in wa1. These words will be 
reducible using the word difference set diff2', but not be reducible 
using the word difference set diff2.

S4. Create a list of lhs words sampled from gpname.andnot.
For each  lhs in the list, calculate its reduction rhs using diff2', 
so that lhs=rhs. Then calculate the word differences 
lhs(i)^-1rhs(i) for i ranging from 1 to the length(lhs)-1 and 
add any new ones to the word difference set diff2.

Repeat steps S1 to S4 until wa1 and wa2 are equal and, hopefully,
have a small size.

The aim of this procedure is to produce a word acceptor with 
a sufficiently small number of states so that it can 
then be used in the more memory intensive processes shown 
in the examples in the previous section  (for example calling gpcheckx with 
the -p option) in order to extract more word differences.

The gpcheckx options -diagonals s e l  and -diff2name 'diff2suffix' are 
provided to implement steps S1 and S4 of the  above 
procedure. 
The option '-diagonals s e l' indicates that diagonals are to be calculated
and added to gpname.diff2'diff2suffix' according to the filter 
s(tart word difference), e(nd word difference and l(imit or maximimum 
total of diagonals).
Specifying -diagonals 0 0 0 indicates that all possible diagonals will 
added to diff2'.

Example: 3572 calculation using diagonals.

./bin/kbprog -wd -t -me 50000 3572

then repeatedly execute the cycle defined by the 
following script (comments contained in '').  

'calculate 3572.wa1 using 3572.diff2'

./bin/gpcheckx -execwa './dowa 3572 >mfile' -waonly  -v  3572
cp 3572.wa 3572.wa1

'calculate all the diagonals of 3572.diff2 to make the 
larger word difference file, 3572.diff2diaggoody'

./bin/gpcheckx -diagonals 0 0 0 -diff2name diaggoody -w  -v  3572

'calculate 3572.wa2 using 3572.diff2diaggoody 
(using temporary file 3752,diff2d)'

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

where script 'dowa' is 

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


Results: 
After 11 cycles of comparing pairs of large word acceptors, 
each with 230000+ plus states, a 'small' word acceptor 
with 47611 states is built.
The 'build and check multiplier' type process

./bin/gpcheckx -p -v -w 3572 

uses this word acceptor to then extract more word differences 
to add to diff2. 
But these extra word differences result in the building 
of a large word acceptor of 220000+ states.
So we resume the cycle of adding diagonals to the diff2
word difference set and compare word acceptors to 
try and build a smaller word acceptor again. 
This time the process soon finishes with the building 
of a 'small' word acceptor of 47613 states which, this time,  
happens to be the correct word acceptor.
