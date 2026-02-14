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

build a word acceptor that accepts reducible prefix reduced
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

-s 'Base;Increment'

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

(or, more efficiently)

./bin/gpcheckx -v -p -s '2000;2000' f29 


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

(this calculation takes several hours to complete but would
only take several minutes if 3572 were defined in a different 
way. See the 'Minimising time...' section below for more
details.)

./bin/kbprog -wd -me 200000 -t 1000  3572

./bin/gpcheckx -v -p -s ‘60000;10000’ 3572 

./bin/gpcheckx -geo -v 3572 +rptz


h93

 (this calculation requires 8GB of memory)

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


See the 'Add diagonals' section for quicker and 
less memory intensive ways to carry out the 
calculations for h93 and f46.

# Proving non-hyperbolicity.

If we suspect that a group is not hyperbolic 
(for example the group F(4,6), defined above, 
turns out not to be hyperbolic - the interested 
reader is invited to use kbmag tools /bin/wordreduce and f46.wa to
show that the words aCCb, aDDa generate the 
free abelian group ZxZ), it can be useful
to spot patterns, that might indicate the presence 
of 'fat triangles', by tracing the geodesic equations 
which give rise to new geodesic word differences. 
This is achieved by specifying the -ve option. 

For example, the group G(5,3) (provided by Martin Edjvet), defined 
by the file,G53

_RWS := rec(
  isRWS := true,
  ordering := “shortlex”,
  generatorOrder := [t,T,u,U],
  inverses := [T,t,U,u],
 equations := [[(t^2)^5,IdWord],[t^5\*u\*t,u^3]]
);

is automatic and its diff2 and wa files can be calculated 
by, for example,

./bin/kbprog -wd -me 80000 -t 1000 -v G53

./bin/gpcheckx -p -v -s '8000;3000' G53

Then, by examining, the output of

./bin/gpcheckx -geo -ve G35

for new geodesic word differences which are powers
of U\*T^2, it can be straightforwardly shown in kbmag 
that the words

U^3\*T\*U\*t\*u\*T^3\*U\*t^2\*U^2\*T 

and

t\*U\*T\*u\*t\*u^3\*t^2*U^2\*T\*u\*t^2

commute and generate ZxZ.

# Truncating a large word acceptor

gpcheckx produces a list of lhs_words from which we 
attempt to calculate new word differences. 

Each such lhs_word, w, traces a unique path of 
numbered states in the current word acceptor. 
We call this path of states, the state set of w.

'Truncation' by a value N means that gpcheckx is 
constrained to produce just those lhs_words whose 
state set only contain states which are <= N.

We specify this truncation by the -s N or 
-s 'Base;Increment' switches.

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

 -diagonals -diff2name 'suffix1' -prevdiff2 'suffix2'

indicate that all possible diagonals are to be 
added to gpname.diff2 and saved to gpname.diff2'diff2suffix1'
with the exception of any diagonals that also happen
to be present in gpname.diff2'suffix2'.

This enlarged set of word differences is to be used 
to build the word acceptor and reduce words in
subsequent processing.

Example: f29 calculation using diagonals.

./bin/kbprog -wd -t 1000 -me 9000 f29

./bin/gpcheckx -diagonals -diff2name diags f29

'f29.wa now correct'

./bin/gpcheckx -p -diff2name diags -w +rptz

'f29.diff2 now correct'

# The 'Add Diagonals' script. 

This script has two parts. We use $1 as the variable used
to contain the group name (gpname).

Part 1. Add diagonals to the current set of word differences,
and store this larger set in gpname.diff2diags - but dont
include any diagonals present in any existing file
gpname.diff2diags - then build a word acceptor based on this 
which is then stored in gpname.wa2.

./bin/gpcheckx -diagonals -diff2name diags -prevdiff2 diags -v $1

cp $1.wa $1.wa2

Part 2. Compute the word acceptor based on 
gpname.diff2, place this in gpname.wa1, then extract 
word differences from the differences found between the
word acceptors gpname.wa1 and gpname.wa2.

./bin/gpcheckx  -waonly  -v $1

cp $1.wa $1.wa1

cp $1.wa2 $1.wa

./bin/gpcheckx -v  -andnot $1.wa1 $1.wa2  $1

./bin/gpcheckx -t  -diff2name diags -v $1

EXAMPLES

1. F(4,6)

F(4,6) can be calculated using the following 'recipe':-

./bin/kbprog -wd -t 1000 -me 60000 f46

Run parts 1 and 2 of the 'add diagonals' script four times.

 then do

./bin/gpcheckx -v -w -diff2name diags -m -s 40000 f46

./bin/gpcheckx -v -m -s '60000;20000' f46 

./bin/gpcheckx -v -p f46 +wrptz

f46.diff2 and f46.wa will now be correct.


2. h93 

(this calculation uses < 2GB of memory)

./bin/kbprog -wd -t 1000 -me 50000 h93

Run parts 1 and 2 of the add diagonals script. Then 
run part 2 a second time. Run part 1 once more. 
The word acceptor h93.wa  will then be corect. 
However, the scripts need some adjustments for this 
to work, and should appear as

(part 1)

./bin/gpcheckx  -diagonals  -diff2name diags -v h93

cp h93.wa h93.wa2

(part 2)

./bin/gpcheckx -waonly -v h93

cp h93.wa h93.wa1

./bin/gpcheckx -v -andnot h93.wa1 h93.wa2 h93

cp h93.wa2 h93.wa

./bin/gpcheckx -tt 1 -lineitems 150 0 -diff2name diags h93

(the scan of the lhs words will be split into
several smaller scans, each processing a line 
of output at a time)

(repeat part 2)  

./bin/gpcheckx -waonly -v h93

cp h93.wa h93.wa1

./bin/gpcheckx -v -andnot h93.wa1 h93.wa2 h93

cp h93.wa2 h93.wa

(process just 1 line of output)

./bin/gpcheckx -tt 1 -lineitems 100 1 
                 -v -diff2name diags h93

(do part 1 again)

./bin/gpcheckx  -diagonals -prevdiff2 diags 
           -notbigger -diff2name diags -v h93

(The switch -notbigger indicates that only those diagonals 
whose length is the same as the length of the word 
differences that they are a diagonal of will be selected.)

h93.wa will now be correct. Follow up with,

./bin/gpcheckx -w -diff2name diags -v -m  h93 

./bin/gpcheckx -v -p  -to 200 h93 +wrptz

h93.wa and h93.diff2 will now be correct. 

The correct geodesic word acceptor and correct geodesic 
word difference set, which show that h93 is hyperbolic, 
can then be calculated by

./bin/gpcheckx -geo -slowmin -v h93 +rptz

(the -slowmin switch causes the transitions table of 
very large fsa's to be read in a line at a time rather 
than loaded into memory.)

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

For a given N, the '-m -s N' option does the  building of a 'triples'
part with the least memory/processing time requirement. The smaller N is, 
the smaller the processing time and memory requirements will be.
But if N is too small then no new word differences will be found.

The options '-p -s N' require more time and memory.
The option '-p' with no -s option requires the most
time and memory of all.

Also, it cannot be emphasised enough that even a trivial
looking change to a groups presentation can affect the 
difficulty of the automaticity and hyperbolic calculations 
drastically.  

For example if five redundant generators and their 
inverses are added to the group 3572, defined above,  
so that its definition now reads

_RWS := rec
(
  isRWS := true,

  ordering := "shortlex",

  generatorOrder := [a,A,b,B

,ab,BA,ba,AB,bb,BB,ab2,BA2,ab3,BA3
],
  
inverses := [A,a,B,b

,BA,ab,AB,ba,BB,bb,BA2,ab2,BA3,ab3
],
  
equations :=
  [
    
[a^3,IdWord],
   
[b^5,IdWord],
    
[(a\*b)^7,IdWord],
   
 [(a\*b\*A\*B)^2,IdWord]
    
,[bb,b\*b], [ab,a\*b], [ab2,(a\*b)^2], 
    
 [ab3,(a\*b)^3], [ba,b\*a] 
  ]
);

then the automatic and hyperbolic calculations,
using, for example,

./bin/kbprog -wd -me 20000  3572;
./bin/gpcheckx -p -s '2000;2000' 3572;
./bin/gpcheckx -geo 3572 +rptz

are dramatically eased and will now complete in a matter of 
minutes rather than hours!

With the exception of the first h93 calculation,
the examples given here will run with 2 GB of available memory.

# References

The interested reader is referred to the following
publications that discuss the groups used in the 
examples above.

3572

G. Havas and D.F. Holt. 
On Coxeter’s families of group presentations. 
J. Algebra 324(5) (2010), 1076-1082. 48

h93

Ihechukwu Chinyere, Gerald Williams,
Hyperbolic groups of Fibonacci type and T(5) cyclically presented groups

J. Algebra, 580:104-126, 2021.

(https://arxiv.org/pdf/2008.08986)

f46, f29

Christopher P Chalk,
Fibonacci Groups with Aspherical Presentations
Comm. Algebra, 26(5):1511-1546, 1998

f29  

D.F. Holt, B. Eick and E.A. O’Brien, 
Handbook of Computational Group Theory, 
CRCPress, Boca Raton, 2005. Section 13.4

