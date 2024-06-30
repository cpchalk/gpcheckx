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

 In order to run, the presence of the word difference file ('file name'.diff2)
is required which has been built by the KBMAG binary kbprog. 

gpcheckx adresses the problem of when the the provisional word acceptor, 
is simply too large for further orocessing in the kbmag system and 
proposes two solutions.

1. Truncate the Word Acceptor to a more acceptable size.

See the 'Truncating' section.

2. Speculatively add so-called diagonal word differences
to the current word difference set in the hope that the
word acceptor built from this larger set becomes 'small'. 

See the 'Adding Diagonals' section.

# switches (selection)
-p 

build a word acceptor that accepts prefixed reduced
words, costruct an fsa that recognises 'fellow travellor'
lhs=rhs equations , derive nee lhs=rhs equations from
which new word differences can be calculated.

-m

as above, but build a word acceptor that accepts minimally
reducible words.

-geo 

compute the geodesic word acceptor and extract new 
word differences

-s N

truncate word acceptor to N states

-s 'Base;-Increment'

truncate the word acceptor by Base states but repeat
the action and increase the truncation value by Increment
states at a time until the truncation value is larger
than the number of states of the word acceptor.

+rptz 

repeat the specified until no new word differences are
found

# Examples of use. 

The file f29 defines 
the Fibonacci group F(2,9) and can be
found in the KBMAG test library.

./bin/kbprog -v -wd -me 65000 -t 1000 f29
(creates an initial file f29.diff2)

./bin/gpcheckx -v -p f29 +rptz

(f29.wa and f29.diff2 will be correct)

./bin/gpcheckx -v -geo f29 +rptz

( f29.geowa will be correct)

The 'difficult' examples 3572 and h93.

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

# Truncating a large word acceptor

If a word acceptor, W, consists of M states, we can 
truncate it to smaller number, N,  of states by treating 
all states with number greater than N as fail states. 
We specify this by the -s N switch.

# Using diagonals to build a 'small' and possibly correct word acceptor.  
It has been observed that, in many cases, the complete 
set of word differences of an automatic group consists entirely, 
or almost entirely, of so called 'diagonal' word differences. 
If two word differences wd1 and wd2 satisfy the equation 
wd2=g1^-1wd1g2, for some generators g1, g2,
then the word diagwd=g1^-1wd1 is called a diagonal of wd1. 


This observation suggests that adding all possible 
diagonals to a the word difference set in diff2 might 
result in the correct word acceptor being built.

The options

 -diagonals 'diff2suffix'

 -diff2name 'diff2suffix'

indicates that all possible diagonals are to be 
added to gpname.diff2'diff2suffix', and
this enlarged set of word differences is to be used 
to build the word acceptor and reduce words in
subsequent processing.

Example: f29 calculation using diagonals.

./bin/kbprog -wd -t 1000 -me 9000 f29

./bin/gpcheckx -diagonals -diff2name diaggoody -waonly f29

'f29.wa now correct'

./bin/gpcheckx -p -diff2name diaggoody -w +rptz

'f29.diff2 now correct'



