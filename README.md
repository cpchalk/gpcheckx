# gpcheckx

A kbmag program that can be used to efficiently find the word differences  
needed to build the correct automatic structure of a group gpname. 
It uses some ideas from the MAF system for building automatic structures. 
For background material and concepts, see  

https://maffsa.sourceforge.net/manpages/background.html

gpcheckx extends the set of binaries making up the kbmag system 
authored by Derek Holt at the University of Warwick. See

https://github.com/gap-packages/kbmag/tree/master/standalone


It provides the equivalent functionality of the gpcheckmult component
of gpmakefsa (-p option) and the part of gpgeowa which calculates the 
word differences needed to build the correct geodesic word acceptor 
gpname.geowa (-geo option). It's aim is to find word 
differences efficiently so that the resultant word difference file 
(gpname.diff2) can be used to build the correct word acceptor (gpname.wa) 
and correct general multiplier (gpname.gm).

It consists of one source file gpcheckx.c.

To build, edit the kbmag src/makefile to include
gpcheckx in the list of binaries to be built, using, for example,
the build specification for gpcheckmult as a template.  

In order to run, the presence of the word difference file (gpname.diff2)
is required which has been built by the kbmag binary kbprog. 

gpcheckx adresses the problem of when the provisional word acceptor
is simply too large for further processing in the kbmag system and 
proposes two solutions.

1. Truncate the word acceptor to a more acceptable size.

See the 'Truncating' section below.

2. Speculatively add so-called diagonal word differences
to the current word difference set in the hope that the
word acceptor built from this larger set becomes 'small'. 

See the 'Using diagonals' section below.

gpcheckx is run iteratively using the accumulated
information from previous runs. It finally finishes 
with success if the process

gpcheckx -p -v gpname 

finishes with no new word differences being discovered.

In this case, the file gpname.diff2 is deemed correct.

# Switches (selection)
-p 

build a word acceptor that accepts prefix reduced
words, construct an fsa that recognises 'fellow traveller'
lhs=rhs equations, derive new lhs=rhs equations from
which new word differences can be calculated.

-m

as with -p, but instead build a word acceptor that accepts minimally
reducible words.

-geo 

compute the geodesic word acceptor and extract new 
word differences

-s N

truncate the word acceptor to N states

-s 'Base;-Increment'

truncate the word acceptor to Base states but repeat
the action and increase the truncation value by Increment
states at a time until the truncation value is larger
than the number of states of the word acceptor.

+rptz 

repeat the specified action until no new word differences are
found.

-to S

timout scanning process after S seconds, see the 
'minimising time and memory' section.

-w 

read the current word acceptor, don't build a new one.

-v 

verbose mode, display progress messages etc.


# Examples of use. 

The file f29 defines the Fibonacci group F(2,9) 
and can be found in the kbmag folder ag_data.

./bin/kbprog -v -wd -me 65000 -t 1000 f29

(create an initial file f29.diff2)

./bin/gpcheckx -v -p f29 +rptz

(f29.wa and f29.diff2 will be correct)

./bin/gpcheckx -v -geo f29 +rptz

(f29.geowa will be correct)

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
F(4,6), defined by the file f46 containing

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

The correct diff2 and wa files can be
calculated as follows.

./bin/kbprog -wd -t 1000 -me 200000 f46

./bin/gpcheckx -v -p -s ‘4000;1000’ f46

./bin/gpcheckx -v -p f46 +rptz

See the 'Add diagonals' section for a quicker way to do this
calculation.

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
diagonals to a the word difference set might 
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

# The 'Add Diagonals' Script. Another way of calculating F(4,6).

This script has two parts.

Part 1 Add diagonals to the current set of word differences,
store this larger set in gpname.diff2diaggood and build
a word acceptor based on this which is then stored in gpname.wa2.

./bin/gpcheckx  -diagonals  -diff2name diaggood -waonly  -v  $1;

cp $1.wa $1.wa2

Part 2. Compute the word acceptor based on 
gpname.diff2, place this in gpname.wa1, then extract 
word differences from the differences found between the
word acceptors gpname.wa1 and gpname.wa2.

./bin/gpcheckx  -waonly  -v  $1

cp $1.wa $1.wa1

cp $1.wa2 $1.wa

./bin/gpcheckx -v -h -andnot $1.wa1 $1.wa2  $1

./bin/gpcheckx -t -to 400 -diff2name diaggood -v $1

F(4,6) can now also be calculated using the following 'recipe':-

./bin/kbprog -wd -t 1000 -me 60000 f46

Run parts 1 and 2 of the 'add diagonals' script 3 times.

Run part 2 only of the 'add diagonals' script 2 times, then do

./bin/gpcheckx -v -w -diff2name diaggood -m -s 80000 f46

./bin/gpcheckx -v -m f46 +rptz

./bin/gpcheckx -v -p f46 +rptz

f46.diff2 and f46.wa will now be correct.

# Minimising time and memory requirements

The gpcheckx operation can be divided into three 'costly' parts:-

building a word acceptor

scanning lists of lhs words for new word differences

building a 'triples' fsa to recognise fellow travelling
lhs=rhs equations and calculate new ones.

There is currently no ability to reduce the time and 
memory reqirements to build a word acceptor. The only solution 
to improve this is to either have more patience or use a faster 
computer with more memory.

The -to S option stops the scanning lists process after S seconds.
In addition, a SINGLE Control & C from the keypad will also cause
the scanning process to stop.

For a given N, the '-m -s N' option does the  building a 'triples'
part with the least memory/processing time requirement. The smaller N is, 
the smaller the processing time and memory requirements will be.
But if N is too small then no new word differences will be found.

The options '-p -s N' require more time and memory.
The option '-p' with no -s option requires the most
time and memory of all.

With the exception of the h93 calculation, the examples 
given here will run with 2 GB of available memory.
