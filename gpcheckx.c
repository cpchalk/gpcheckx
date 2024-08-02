#define WADIAG 0
#define BASE 2000
#define INC 1000
#define STATE_LIMIT 0
#define CALC_GEODIFF 0 // in geo mode calculate the optimal size of diff2
		       // that produces correct geowa
#define USEALT 1 // try to reduce duplicates by alternative prestate routes
#define NOMINMZ 1 // dont minimize the minred fsa
#define NOEXMZ 1 // dont minimize the exists fsa'
#define MZBH   0 // use fsa_minimize_big_hash instead of fsa_minimize
#define TRACE5 0 // diagnostic for fsa_wa_x
#define TRACE4 0 // diagnostic for minimize big hash
#define TRACE3 0 //  diagnostic for problems in one_level_reducible
#define TRACE2 0 // si<1500 //diagnostic for crash in processing hein group
#define TRACE1 0 //   diagnostic for tracing actions in one_level_reducible
#define MAXWDLENGTH 100
#define SCALE 3
#define NR_OF_DOTS 1600 // 20 lines of dots ??
#define DEFAULT_SIZE_OF_DOT 1000  // words processed
#define MAX_WORD_LEN 10000
#define DIFF2_INCREASE_FACTOR 3
#define BIG_HASH_SIZE 24000000
#define MAXPAIRS 4096 // max size of wd,wa state pairs in calculation
                      // for one_level_reducible
#define MAX_EXTRA_MEMORY 100  // number of times an extension of states
                              // table can be requested in fsa_triples_big_hash
#define DISPLAY_HASH_DEPTH 0 // display 'depth' of big hash table
#define MAXV 65536 /* The maximum number of vertices allowed. */

#define Printf if (kbm_print_level>1) printf
#define PPrintf if (kbm_print_level>0) printf
#include <stdio.h>
#include <signal.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include "defs.h"
#include "fsa.h"
#include "definitions.h"
#include "rws.h"
#include "hash.h"
#include "externals.h"
extern int (*reduce_word)();

static FILE *rfile, *wfile;

//void fsa_print ();
//int  fsa_minimize ();
fsa * fsa_genmult2_int_x();
//int genstrlen ();
//void genstrcpy ();
//int genstrcmp ();
//void fsa_clear_rws ();
//void fsa_set_is_accessible ();
//void srec_clear ();
//void fsa_table_init ();
//int hash_rec_len ();
void  interrupt_gpcheckx();
//fsa  *fsa_diff();
fsa *fsa_pfxred();
int process_words();
void  badusage_gpcheckx();
int  (*reduce_word)();
int getprefix();
void write_fsa ();
int  ***copy_table_dptr();
boolean look_for_diffs();
boolean one_level_reducible ();
void display_eqn ();
void update_table ();
void do_andnot ();
fsa *fsa_triples_big_hash();
fsa *fsa_geopairs_big_hash();
int fsa_minimize_big_hash();
int get_image_big_hash();
int get_image2_big_hash();
fsa *fsa_minred1 ();
//fsa *fsa_minred ();
fsa *fsa_beginswith ();
int  diff_reducex();
/* Functions used in this file and defined elsewhere */
//int sparse_target();
//void fsa_init();
//int  fsa_table_dptr_init();
//void fsa_set_is_accepting();
//void srec_copy();
//boolean  srec_equal();
//void fsa_clear();
//void compressed_transitions_read();
//int  hash_locate();
//int* hash_rec();
//void hash_init();
//void hash_clear();
//int  diff_reduce();
int  calculate_inverses();
//int stringlen();
//void  fsa_read();
//void  fsa_copy();
//int   add_wd_fsa();
//int  make_full_wd_fsa();
int  make_full_wd_fsax();
int add_diagonals_to_wd_fsa();
fsa *fsa_wa_x(); 
fsa *fsa_extractwdsfromwa(); 
fsa *fsa_extractwdsfrom_recovery(); 
fsa *fsa_ex_min(); 
//fsa  *fsa_minred();
//fsa  *fsa_minkb();
//fsa  *fsa_geopairs();
//fsa  *fsa_exists();
//fsa  *fsa_and();
//fsa  *fsa_and_not();
fsa  *fsa_and_not_first();
fsa  *fsa_binop_x();
//fsa  *fsa_not ();
int verify_word ();
void free_fsa ();
int main1 ();
int main2 ();
void analyse_diff2();
fsa *find_geo_wds();

// fsa_triples_big_hash structures
typedef struct mstate{
    unsigned int wa1;
    unsigned int wa2;
    unsigned int diff;
} mstate;

typedef struct kts{
    int state;
    mstate value;
    struct kts *next_value;
} kts;

// fsa_geopairs_big_hash structures
typedef struct mstate3{
    unsigned int wa;
    unsigned int diff;
} mstate3;

typedef struct kts3{
    int state;
    mstate3 value;
    struct kts3 *next_value;
} kts3;
boolean gpcheckx_onintr;
/*
  This program takes as input the file gpname.diff2. It then aims to
  discover new word differences and use these to produce a larger diff2 file.
  It tries to achieve this by doing the following steps.
  1) build a word accepter wa from the diff2
  2) Use the wa to build an fsa minred recognising the minimally reducible words.
     These correspond to the lhs of equations produced by the Knuth Bendix
     process.
  3) Use minred, wa and diff2 to build an fsa, minkb, which recognises all known
     lhs->rhs equations where lhs is minimally reducible, rhs is reduced and
     lhs=rhs.
  4) Use minred to build an fsa, Exists, which recognises all 'known' minimally
      reducible words.
  5) Use minred and Exists to build an fsa, Andnot, which recognises  unknown minimally
     reducible words (those recognised by fsa minred but not recognised by fsa Exists.
     Andnot is written to the file gpname.andnot.  A set, S, of words recognised by
     Andnot is produced where the number of words in S is the number of states of
     Andnot. For each state S a word is produced which is traces the shortest
     path from 1 to S and then traces the shortest path from S to the accept
     state.
  6) Sample some or all of the words of S (depending on the -d switch or default).
     Reduce each word using diff2 to obtain an equation from which new word differences
     can be discovered with which to build a bigger diff2.
  7) Optionally (-p switch) build an fsa pfxred recognising reducible words which
     are prefix reduced. pfxred is simply the wa with the fail state replaced by
     an accept state, Use pfxred in place of minred in steps 3), 4), 5) and 6) above
     to produce further equations from which further new word differences may be
     extracted.

   Repeating steps 2) to 7) frequently enough will produce a diff2 file from which
   the correct wa and gm can be built - subject to axiom checking.

     Usage: gpcheckx [-s] [-v][-ve][-p][-rp][-m][-w][-min][-t][-d n]
                    [-mx n][-mm n][-x] [-nf] [-ht n] [-minkb] 
		    [-tt n] [-andnot file1 file2] [-minred1] 
		    [-scale n] [-counter] [-mzbh] [-exmz] 
                    [-noalt] [-analyse] [-geo][-to][-gw] gpname [+rptz|+wrptz]

		    
     -v  set trace level 2
     -ve set equation tracing - the new word differences and the equation
                               they are extracted from will be displayed.
     -p  look for new wds from equations with a prefix reduced lhs
     -rp  look for new wds from equations with a 'restricted' prefix reduced 
	 lhs wg, where wg=v, g a generator and w and v are reduced and do not
	 share a prefix. This class includes minred but is contained in pfxred.
     -m  look for new wds from equations with a minimally reducible lhs.
	 (Default)
     -w  dont build the wa but read from an existing wa file. The file
         gpname.wa is expected to exist.
     -min  dont build the minred fsa but read from an existing wa file.
          The file gpname.minred is expected to exist.
     -minkb  dont build the minkb fsa but read from an existing minkb file.
          The file gpname.minkb is expected to exist.
     -t  dont build any fsa's but attempt to discover new word differences
         from an existing Andnot file and the diff2 file it was built from.
         The file gpname.andnot is expected to exist. 
     -tt n  As t above but start n words into the andnot file
     -d n  the process of examining words in the Andnot set, S, 
           is monitored by printing a dot, '.', for every 1000 words processed.
           Any discovered word differences are indicated by an 'x' or traced 
           in the case that -ve is set.
           n determines the number of '.' s displayed and so the number of
           words examined. If n is negative, -m say, then the first m*1000 words
           will be processed. If n is positive then all words from S will
           be sampled at regular intervals to cause n*1000 words to be processed.
     -mx n  only process words with length <= n
     -mm n  only process words with length => n
      -x   dont update the existing gpname.diff2 file.
      -nf  no filter - search the set of minimally reducible or prefix reduced
           words directly for missing wds. Omit the minkb, exists and andnot
           'filtering' steps.
     -ht n use the fsa_triples_big_hash procedure with a hash table size of n
           million  entries. A table built with more than n million states
           will be truncated to n million states but will still produce
           lhs words that can be examined for missing wds. If n is 0, then
           the kbmag standard fsa_minkb procedure is used.
           The default is to use fsa_triples_big_hash with a hash table size
           of 24 million entries.
    -noalt dont look for an alternative prefix when we detect we are about
	   to process a word identical to one that has already been processed.
     -l/-h large/huge hash-tables (for big examples)

   -andnot file1 file2 gpame  do an and_not on the fsa's defined by
           file1 and file2 and place the result in file gpname.andnot
   -scale n  adjust the number of words represented by a '.' The number is
             10 to the power of n where n is in the range 1 to 6.
   -counter  during the minkb fsa building phase display a count of the 
	     states processed against the states created.
   -mzbh     use a big hash table for minimisation 
		(still in experimental mode).
   -exmz     minimize the exists fsa
   -noexmz   dont minimize the exists fsa (Default)
   -minmz    minimize the minred fsa when -nf switch is set.
   -nominmz  dont minimize the minred fsa when -nf switch is set. (Default)
   -minred1  build an fsa which recognises words of the form gr, g a 
	     generator, r a word, where  r is minimally reducible. 
	     Scan lhs words from this fsa for possible new wd's..
   +rptz     repeat the command before the +rptz characters until no new 
	     word differences are found.
   +wrptz    repeat the command before the +wrptz characters until no new 
	     word differences are found. Read the current wa on the first call.
   -analyse  display a list of the current wd's in gpname.diff2  and whether 
             they belong to diff1, diff2 only or to neither diff1 nor diff2.
	     The presence of files gpname.diff1c and gpname.diff2c produced by
             gpminkb at the end of a successful fsa build is assumed. 
   -geo      search for word differences associated with a geodesic wa
             If none can be found output the calulated group.geowa.
   -to n     Set timeout to n seconds. This timeouts the scan of andnot after
	     n seconds and updates the diff2 with the wds extracted so far. 
   -gw       dont build the geowa but read from an existing geowa file. The file
             gpname.geowa is expected to exist.
      
  By default 1600000 (NR_OF_DOTS * 1000) words are sampled from the 
  set S of words recognised by the fsa Andnot.
 */
void  interrupt_gpcheckx()
{
  gpcheckx_onintr = TRUE;
  signal(SIGINT, SIG_DFL);
  signal(SIGKILL, SIG_DFL);
  signal(SIGQUIT, SIG_DFL);
}
 
int main(argc, argv)
        int             argc;
        char           *argv[];
{
	boolean v_option = FALSE;
	boolean l_option = FALSE;
	boolean h_option = FALSE;
	int options=0;
	int arg = 1;
	char *argvx_data;
	char **argvx;
	int outcome = 0;
	int wa_size = 0;
	char gpname [100];
        boolean gp_present = FALSE;
	//char argvx [10][100];
	//printf("calling main1");
	while (argc > arg) {
		if (argv[arg][0] != '-') 
			gp_present = TRUE;
    		else if (strcmp(argv[arg], "-v") == 0) {
        		if (arg >= argc)
        			badusage_gpcheckx(FALSE);
			v_option = TRUE;
			options++;
    		}
    		else if (strcmp(argv[arg], "-l") == 0) {
      			if (arg >= argc)
        			badusage_gpcheckx(FALSE);
			l_option = TRUE;
			options++;
		}
    		else if (strcmp(argv[arg], "-h") == 0) {
      			if (arg >= argc)
        			badusage_gpcheckx(FALSE);
			h_option = TRUE;
			options++;
		}
		else if (!strcmp(argv[arg],"-s")) {
      			if (arg >= argc)
        			badusage_gpcheckx(FALSE);
			char *qual=argv[++arg];
			// does it contain a ';'
			boolean found=FALSE;
			int len =strlen(qual);
			int i=0;
		        int semicolon1=0, semicolon2=0;	
			while (i<len) {
				if (qual[i]==';')
					if (semicolon1>0)
						semicolon2=i;
					else
					    semicolon1=i;
				if (semicolon2>0)
					break;
				i++;
			}		
			if (semicolon1==0)
				continue;
			// call main1 repeatedly with various values of -s
	qual[semicolon1]='\0';
	int state_limit=atoi(qual);
	if (state_limit==0)
		printf("invalid statelimit for -s switch\n");
	int max = 0;
	if (semicolon2 > 0) {
		qual[semicolon2]='\0';
		max=atoi(qual+semicolon2+1);
		if (max==0)
			printf("invalid maximum for -s switch\n");
	}
	int inc=atoi(qual+semicolon1+1);
	if (inc==0)
		printf("invalid increment for -s switch\n");
       // copy all the arguments apart from the two -s ones into this area. 
	i=1;
	int ll;
	int argcx=0;
	char ** argvx;
      	tmalloc(argvx,char *,argc+1);
	argvx[0]=argv[0];
	if (argv[argc-1][0]=='+')
	  // ignore any postfixes
		argc--;
      	for (ll=1;ll<argc;ll++) {
		if (strcmp(argv[ll],"-s"))
			argvx[i++]=argv[ll];
		else
		// skip -s and its argument
			ll++;
	}
	// there should be room for 3 more items
	int state_limit_pos=i;
	argcx=i+2;	
	argvx[i-1]="-s";
	argvx[argc-1]=argv[argc-1]; //gpname
	argvx[argc]="+rptz";
	wa_size=1000000000;
	while (state_limit<wa_size ) {
		if (max>0 && state_limit>max)
		{
			printf("maximum state limit exceeded\n");
			exit(0);
		}
                char s_value [20];
		sprintf(s_value,"%d",state_limit);
		argvx[state_limit_pos]=s_value;
		Printf("\n**** -s %d ****\n\n",state_limit);
		outcome=main1(argcx,argvx,&wa_size);
		//printf("wa size is %d\n",wa_size);
		state_limit+=inc;
	}
	// do a +rptz to finish off
	outcome=main1(argcx+1,argvx,&wa_size);
			exit(0);

		}
      		arg++;
	}
	if (argc == options+2 && gp_present)
	{
      		strcpy(gpname,argv[argc-1]);
		//printf("general option requested for %s %s!\n",gpname,argv[1]);
	}
	else {
		outcome=main1(argc,argv,&wa_size);
		//printf("return from specialised call %d\n",outcome);
		exit(outcome);
	}
	clock_t start_t, end_t, total_t;
  	start_t=clock();
	int op, rpt;
	int ll;
	tmalloc(argvx_data,char,(argc+1)*100);
      	tmalloc(argvx,char *,argc+1);
	strcpy(argvx[0],"not in use");
	arg=1;
	if (h_option)
	   strcpy(argvx[arg++],"-h");
	else if (l_option)
	   strcpy(argvx[arg++],"-l");
	if (v_option)
	   strcpy(argvx[arg++],"-v");
	op=arg;
	strcpy(argvx[arg++], "-nf");
	strcpy(argvx[arg++],gpname);
	rpt=arg;
	strcpy(argvx[arg++],"+rptz");
	printf("doing -nf option - scan minimally reducible lhs's for possible wds\n");
	outcome=main1(arg,argvx,&wa_size);
	strcpy(argvx[op], "-rp");
	strcpy(argvx[rpt], "+wrptz");
	printf("doing -rp option - scan prefix reduced lhs's for possible wds\n");
	outcome=main1(arg,argvx,&wa_size);
	strcpy(argvx[op], "-m");
	printf("doing -m option - calculate missing wds for minimally reducible lhs\n");
	outcome=main1(arg,argvx,&wa_size);
	strcpy(argvx[op], "-p");
	printf("doing -p option - calculate missing wds for prefix reduced reducible lhs\n");
	outcome=main1(arg,argvx,&wa_size);
	printf(", %s.wa should now be correct\n",gpname);
	printf("(type gpcheckx -geo [-v] [-l] %s +rptz to check for possible hyperbolicity)\n",gpname);
   	end_t=clock();
   	total_t=(end_t-start_t)/CLOCKS_PER_SEC;
   	printf("elapsed time %lds\n",total_t);
	tfree(argvx_data);
      	tfree(argvx);
}
int main1(argc, argv,wa_size)
        int             argc;
        char           *argv[];
	int *wa_size;
{
	int outcome=0;
	boolean read_last_wa=FALSE;
	int base=BASE;
	int inc=INC;
    	if (strcmp(argv[argc-1], "+wrptz") == 0) 
		read_last_wa=TRUE;
    	if (read_last_wa || (strcmp(argv[argc-1], "+rptz") == 0)) {
		Printf("repeat till zero word diffs!\n");
		outcome=1;
		while (outcome>0) {
			outcome=main2(argc-1,argv,read_last_wa,wa_size);
			read_last_wa=FALSE;
		}
	}
        else if ( argc>5 && !strcmp(argv[1],"-tt") && !strcmp(argv[3],"-lineitems") ) {  
		outcome=atoi(argv[2]);
		int repeats=atoi(argv[5]);
	        int cycle=1;
		while (outcome > 0) {
			Printf("Cycle %d\n",++cycle);
			repeats--;
                        char buf1 [100];
			int x=sprintf(buf1,"%d",outcome);
			argv[2]=buf1;
			outcome=main2(argc,argv,FALSE,wa_size);
			if (repeats == 0)
				break;
		}
	}
	else {
		outcome=main2(argc,argv,FALSE,wa_size);
	}
	return(outcome);
}

int main2(argc, argv, read_last_wa,wa_size)
        int             argc;
        char           *argv[];
	boolean         read_last_wa;
	int 		*wa_size;	

/* This program uses some ideas and code from MAF by Alun Williams.
   See http://maffdsa.sourceforge.net
*/

{ int arg, i, *inv, old_ndiff, numeqns, ngens;
  fsa  *diff2, *new_diff2, *diag_diff2,*diff2diag;
  char gpname[100],inf1[100], inf2[100], inf2x[100],inf3[100],inf4[100],inf5[100],
       inf6[100],inf7[100], inf8[100], inf9[100], outf[100], fsaname[100], tempfilename[100];
  fsa *wd_fsa; /* This is for doing word-reductions in the case that we
              * correct the diff2 machine
              */
  fsa *gpwa,  *gpwa2, *gpwa_beginswith1, *gpwa_beginswith2,*minred, *minred1,*pfxred,  
	*minkb, *fsaexists,  *fsaandnot_ptr, *geowa;

  boolean look_for_diff2_wds, p_switch,  use_andnot, use_andnotx, read_wa, read_geowa,
    read_minred, read_minred1,read_minkb, trace_equations;
  boolean diff2c = FALSE;
  boolean update_diff2 = TRUE;
  boolean restricted_pfxred= FALSE;
  boolean no_filter = FALSE;
  boolean verify = FALSE;
  boolean recovery = FALSE;
  boolean diagnostics = FALSE;
  boolean overview = FALSE;
  boolean dump_files = FALSE;
  boolean dogeowds = FALSE;
  boolean nobigger = FALSE;
  int maxwordlen=MAX_WORD_LEN;
  int minwordlen=1;
  reduction_equation *eqnptr;
  reduction_struct rs_wd;
  reduction_struct rs_wd2;
  storage_type ip_store = DENSE;
  int dr = 0;
  boolean seengpname;
  boolean outputwords = FALSE;
  int separator=0;
  storage_type op_store = DENSE;
  storage_type op_store2 = SPARSE;
  int state_limit=STATE_LIMIT;
  int no_dots=NR_OF_DOTS;
  unsigned int scale=SCALE;
  gen null_char = 0x0;
  int hash_table_size;
  boolean use_big_hash = TRUE;
  boolean minimize_big_hash=MZBH;
  int start_scan_from = 1; 
  int timeout = 0;
  gpcheckx_onintr=FALSE;
  setbuf(stdout,(char*)0);
  setbuf(stderr,(char*)0);

  inf1[0] = '\0';
  inf2[0] = '\0';
  inf2x[0] = '\0';
  inf3[0] = '\0';
  inf4[0]=0x0;

  outf[0] = '\0';
  arg = 1;
  seengpname=FALSE;
  look_for_diff2_wds=FALSE;
  p_switch=FALSE;
  use_andnot=FALSE;
  read_wa=read_last_wa;
  read_geowa=FALSE;
  read_minred=FALSE;
  read_minred1=FALSE;
  read_minkb=FALSE;
  trace_equations=FALSE;
  update_diff2=TRUE;
  hash_table_size=BIG_HASH_SIZE;
  boolean counter=FALSE;
  boolean minimize_andnot=FALSE;
  boolean no_exists_minimize=NOEXMZ;
  boolean no_minred_minimize=NOMINMZ;
  boolean use_alt_prestate = USEALT;
  boolean no_scan=FALSE;
  boolean do_minred1=FALSE;
  boolean do_mult2=FALSE;
  boolean prefixby1=FALSE;
  boolean do_beginswith=FALSE;
  char *bwqual1;
  boolean diff2name_command = FALSE;
  char *diff2namestr;
  char *usediff2str=NULL;
  char *prevdiff2str=NULL;
  boolean hashlimit_command = FALSE;
  int hashlimit=0;
  boolean use_wa = FALSE;
  boolean exmin_command = FALSE;
  boolean extractwdsfromwa_command = FALSE;
  boolean exwa_command = FALSE;
  boolean exminex_command = FALSE;
  char * ex_str_wa, *ex_str_minex;
  boolean do_waonly=FALSE;
  boolean do_minredonly=FALSE;
  boolean do_geo=FALSE;
  boolean add_diagonals=FALSE;
  boolean all_diagonals=FALSE;
  boolean check_diagonals=FALSE;
  int start_diagonals=0;
  int end_diagonals=0;
  int limit_diagonals=0;
  int max_state=0;
  int max_level=0;

  int verify_qualifier;
  int recovery_num;
  int line_items=0;
  int max_diff2_size=0;
  int min_wd_size=0;
  kbm_huge = TRUE;
  use_andnotx = FALSE;
  use_andnot = FALSE;
  while (argc > arg) {
    if (strcmp(argv[arg],"-dogeowds") == 0) {
      dogeowds = TRUE;
    }
    else if (strcmp(argv[arg],"-notbigger") == 0) {
      nobigger = TRUE;
    }
    else if (strcmp(argv[arg],"-diagnostics") == 0) {
      diagnostics = TRUE;
    }
    else if (strcmp(argv[arg], "-to") == 0) {
      arg++;
      if (arg >= argc)
        badusage_gpcheckx(FALSE);
      timeout = atoi(argv[arg]);
    }
    else if (strcmp(argv[arg], "-ht") == 0) {
      arg++;
      if (arg >= argc)
        badusage_gpcheckx(FALSE);
      hash_table_size = atoi(argv[arg]) * 1000000;
    }
    else if (strcmp(argv[arg], "-analyse") == 0) {
      arg++;
      if (arg >= argc)
        badusage_gpcheckx(FALSE);
      analyse_diff2(argv[arg],max_diff2_size,usediff2str);
      return 0;
    }
    else if (strcmp(argv[arg],"-lineitems")==0)
    {
      arg+=2;
      if (arg >= argc)
        badusage_gpcheckx(FALSE);
      line_items=atoi(argv[arg-1]);
      // repeats is used in main1 
    }
    else if (strcmp(argv[arg],"-recovery")==0)
    {
      arg++;
      if (arg >= argc)
        badusage_gpcheckx(FALSE);
      recovery_num=atoi(argv[arg]);
      recovery=TRUE;
    }
    else if (strcmp(argv[arg],"-maxwdsize")==0)
    {
      arg++;
      if (arg >= argc)
        badusage_gpcheckx(FALSE);
      max_diff2_size=atoi(argv[arg]);
      recovery=TRUE;
    }
    else if (strcmp(argv[arg],"-minwdsize")==0)
    {
      arg++;
      if (arg >= argc)
        badusage_gpcheckx(FALSE);
      min_wd_size=atoi(argv[arg]);
    }
    else if (strcmp(argv[arg],"-verify")==0)
    {
      use_andnot = TRUE;
      verify=TRUE;
      //update_diff2 = FALSE;	
      arg++;
      if (arg >= argc)
        badusage_gpcheckx(FALSE);
      if (strcmp(argv[arg],"pfx")==0)
		verify_qualifier=9999;
      else if (strcmp(argv[arg],"min")==0)
		verify_qualifier=1;
      else if (strcmp(argv[arg],"min1")==0)
		verify_qualifier=2;
      else if (strcmp(argv[arg],"notmin")==0)
		verify_qualifier=0;
      else
      	verify_qualifier = atoi(argv[arg]);
      printf("qual=%d\n",verify_qualifier);
    }
    else if (strcmp(argv[arg],"-diagonals")==0)
    {
	add_diagonals=TRUE;
	all_diagonals=TRUE;
	do_waonly=TRUE;
    }
    else if (strcmp(argv[arg],"-diagonalsonly")==0)
    {
	read_wa=TRUE; // trick to exit after diagonals 
    }
    else if (strcmp(argv[arg],"-maxlevel")==0)
    {
	arg+=1;
	max_level=atoi(argv[arg]);
    }
    else if (strcmp(argv[arg],"-maxstate")==0)
    {
	arg+=1;
	max_state=atoi(argv[arg]);
    }
    else if (strcmp(argv[arg],"-diagonalsextra")==0)
    {
      arg+=2;
      if (arg >= argc)
      {
	printf("-diagonalsextra should be followed by start and end values\n");
        badusage_gpcheckx(FALSE);
      }
	start_diagonals=atoi(argv[arg-1]);
	end_diagonals=atoi(argv[arg]);
	limit_diagonals=0;
	add_diagonals=TRUE;
	all_diagonals=TRUE;
	do_waonly=TRUE;
    }
    else if (strcmp(argv[arg],"-diagonalsextrax")==0)
    {
      arg+=3;
      if (arg >= argc)
      {
	printf("-diagonalsextrax should be followed by start and end and limit values\n");
        badusage_gpcheckx(FALSE);
      }
	start_diagonals=atoi(argv[arg-2]);
	end_diagonals=atoi(argv[arg-1]);
	limit_diagonals=atoi(argv[arg]);
	add_diagonals=TRUE;
	all_diagonals=TRUE;
	do_waonly=TRUE;
    }
    else if (strcmp(argv[arg],"-checkdiagonals")==0)
    {
	check_diagonals=TRUE;
	all_diagonals=TRUE;
	add_diagonals=TRUE;
    }
    else if (strcmp(argv[arg],"-diagonals0")==0)
    {
	add_diagonals=TRUE;
	//all_diagonals=TRUE;
    }
    else if (strcmp(argv[arg],"-geo")==0)
    {
	do_geo=TRUE;
	read_wa=TRUE;
    }
    else if (strcmp(argv[arg],"-rp")==0)
    {
	restricted_pfxred=TRUE;
	look_for_diff2_wds = TRUE;
    }
    else if (strcmp(argv[arg],"-ovw")==0)
    {
	overview=TRUE;
        badusage_gpcheckx(TRUE);
    }
    else if (strcmp(argv[arg],"-mult2")==0)
	do_mult2=TRUE;
    else if (strcmp(argv[arg],"-diff2c")==0)
	diff2c=TRUE;
    else if (strcmp(argv[arg],"-minredonly")==0)
	do_minredonly=TRUE;
    else if (strcmp(argv[arg],"-waonly")==0)
	do_waonly=TRUE;
    else if (strcmp(argv[arg],"-counter")==0)
	counter=TRUE;
    else if (strcmp(argv[arg],"-wt")==0)
	dump_files=TRUE;
    else if (strcmp(argv[arg],"-prefixby1")==0)
   	 prefixby1=TRUE;
    else if (strcmp(argv[arg],"-hashlimit")==0)
    {
		arg+=1;
      	if (arg >= argc)
      	{
		printf("-hashlimit should be followed by an integer \n");
        	badusage_gpcheckx(FALSE);
      	}
	hashlimit_command=TRUE;
       hashlimit=atoi(argv[arg]);
    }
    else if (strcmp(argv[arg],"-prevdiff2")==0)
    {
		arg+=1;
      	if (arg >= argc)
      	{
		printf("-prevdiff2 should be followed by one string\n");
        	badusage_gpcheckx(FALSE);
      	}
	prevdiff2str=argv[arg];
    }
    else if (strcmp(argv[arg],"-usediff2")==0)
    {
		arg+=1;
      	if (arg >= argc)
      	{
		printf("-usediff2 should be followed by one string\n");
        	badusage_gpcheckx(FALSE);
      	}
	usediff2str=argv[arg];
    }
    else if (strcmp(argv[arg],"-diff2name")==0)
    {
		arg+=1;
      	if (arg >= argc)
      	{
		printf("-diff2name should be followed by one string\n");
        	badusage_gpcheckx(FALSE);
      	}
	diff2name_command=TRUE;
	diff2namestr=argv[arg];
    }
    else if (strcmp(argv[arg],"-extractwdsfrommin")==0)
    {
	exmin_command=TRUE;
    }
    else if (strcmp(argv[arg],"-usewa")==0)
    {
	use_wa=TRUE;
    }
    else if (strcmp(argv[arg],"-extractwdsfromwa")==0)
    {
	extractwdsfromwa_command=TRUE;
    }
    else if (strcmp(argv[arg],"-execwa")==0)
    {
		arg+=1;
      	if (arg >= argc)
      	{
		printf("-execwa should be followed by one string\n");
        	badusage_gpcheckx(FALSE);
      	}
	exwa_command=TRUE;
	ex_str_wa=argv[arg];
    }
    else if (strcmp(argv[arg],"-execminex")==0)
    {
		arg+=1;
      	if (arg >= argc)
      	{
		printf("-execminex should be followed by one string\n");
        	badusage_gpcheckx(FALSE);
      	}
	exminex_command=TRUE;
	ex_str_minex=argv[arg];
    }
    else if (strcmp(argv[arg],"-minred1")==0){
      do_minred1=TRUE;
      //read_minred=TRUE;
      read_wa=TRUE;
    }
    else if (strcmp(argv[arg],"-mzbh")==0)
      minimize_big_hash=TRUE;
    else if (strcmp(argv[arg],"-exmz")==0)
      no_exists_minimize=FALSE;
    else if (strcmp(argv[arg],"-noalt")==0)
      use_alt_prestate=FALSE;
    else if (strcmp(argv[arg],"-minmz")==0) 
      no_minred_minimize=FALSE;
    else if (strcmp(argv[arg],"-noexmz")==0)
      no_exists_minimize=TRUE;
    else if (strcmp(argv[arg],"-alt")==0)
      use_alt_prestate=TRUE;
    else if (strcmp(argv[arg],"-nominmz")==0) {
      no_minred_minimize=TRUE;
      no_filter=TRUE;
    }
    else if (strcmp(argv[arg],"-min")==0)
    {
      read_minred=TRUE;
      read_wa=TRUE;
    }
    else if (strcmp(argv[arg],"-min1")==0)
    {
      read_minred1=TRUE;
      read_wa=TRUE;
    }
    else if (strcmp(argv[arg],"-minkb")==0)
    {
      read_minkb=TRUE;
      read_minred=TRUE;
      read_wa=TRUE;
    }
    else if (strcmp(argv[arg],"-nf")==0)
      no_filter=TRUE;
    else if (strcmp(argv[arg],"-v")==0)
      kbm_print_level = 2;
    else if (strcmp(argv[arg],"-ve")==0) {
      kbm_print_level = 2;
      trace_equations=TRUE;
    }
    else if (strcmp(argv[arg],"-p")==0) {
	look_for_diff2_wds = TRUE;
	p_switch = TRUE;
    }
    else if (strcmp(argv[arg],"-m")==0)
	// switch added for symmetry with -p
      look_for_diff2_wds = FALSE;
    else if (strcmp(argv[arg],"-t")==0)
      use_andnot = TRUE;
    else if (strcmp(argv[arg],"-tx")==0) {
      use_andnotx = TRUE;
      use_andnot = TRUE;
    }
    else if (strcmp(argv[arg], "-tt") == 0) {
      arg++;
      if (arg >= argc)
        badusage_gpcheckx(FALSE);
      use_andnot = TRUE;	
      start_scan_from = atoi(argv[arg]);
    }
    else if (strcmp(argv[arg], "-scale") == 0) {
      arg++;
      if (arg >= argc)
        badusage_gpcheckx(FALSE);
      scale = atoi(argv[arg]);
    }
    else if (strcmp(argv[arg], "-andnot") == 0) {
      char *biggerfsa, *smallerfsa;	
      arg+=3;
      if (arg >= argc)
      {
	printf("-andnot should be followed by two fsa names and group name\n");
        badusage_gpcheckx(FALSE);
      }
      do_andnot(argv[arg-2],argv[arg-1],argv[arg]);
      return 0;		
    }
    else if (strcmp(argv[arg],"-w")==0)
      read_wa = TRUE;
    else if (strcmp(argv[arg],"-gw")==0)
      read_geowa = TRUE;
    else if (strcmp(argv[arg], "-d") == 0) {
      arg++;
      if (arg >= argc)
        badusage_gpcheckx(FALSE);
      no_dots = atoi(argv[arg]);
      if (no_dots==0)
      {
	no_scan=TRUE;
	update_diff2=FALSE;
      }
    }
    else if (strcmp(argv[arg], "-s") == 0) {
      arg++;
      if (arg >= argc)
        badusage_gpcheckx(FALSE);
      state_limit = atoi(argv[arg]);
    }
    else if (strcmp(argv[arg], "-mx") == 0) {
      arg++;
      if (arg >= argc)
        badusage_gpcheckx(FALSE);
      maxwordlen = atoi(argv[arg]);
    }
    else if (strcmp(argv[arg], "-mm") == 0) {
      arg++;
      if (arg >= argc)
        badusage_gpcheckx(FALSE);
      minwordlen = atoi(argv[arg]);
    }
    else if (strcmp(argv[arg],"-x")==0)
      update_diff2 = FALSE;
    else if (strcmp(argv[arg],"-l")==0)
    {
      kbm_large = TRUE;
      kbm_huge = FALSE;
    }
    else if (strcmp(argv[arg],"-h")==0)
      kbm_huge = TRUE;
    else if (strcmp(argv[arg],"-noh")==0)
      kbm_huge = FALSE;
    else if (argv[arg][0] == '-')
      badusage_gpcheckx(FALSE);
    else  if (!seengpname) {
      seengpname=TRUE;
      strcpy(gpname,argv[arg]);
    }
    else
      badusage_gpcheckx(FALSE);
    arg++;
  }
  if (state_limit==0)
	dogeowds = FALSE;
  if (!seengpname)
    badusage_gpcheckx(FALSE);
  geowa=NULL;
  gpwa=NULL;  
  minred=NULL;  
  minred1=NULL;  
  pfxred=NULL;  
  minkb=NULL;  
  fsaexists=NULL;  
  fsaandnot_ptr=NULL;  
  diff2=NULL;  
  new_diff2=NULL;  
  diag_diff2=NULL;  
  if (no_filter == FALSE)
	no_minred_minimize=FALSE;
  if (hash_table_size == 0)
  {
      use_big_hash = FALSE;
  }
  strcpy(inf1,gpname);
  strcpy(inf3,inf1);
  strcpy(inf2,inf1);
  strcpy(inf2x,inf1);
  strcpy(inf4,inf1);
  strcat(inf4,".wa");
  strcpy(inf5,inf1);
  strcat(inf5,".minred");
  strcpy(inf6,inf1);
  strcat(inf6,".minkb");
  strcpy(inf7,inf1);
  strcat(inf7,".minred1");
  strcpy(inf8,inf1);
  strcat(inf8,".geowa");
  strcat(inf2,".diff2");
  strcat(inf2x,".diff2");

  if ((rfile = fopen(inf2,"r")) == 0) {
        fprintf(stderr,"Cannot open file %s.\n",inf2);
          exit(1);
   }
  tmalloc(diff2,fsa,1);
  fsa_read(rfile,diff2,DENSE,0,0,TRUE,fsaname);
  fclose(rfile);
  if (fsa_table_dptr_init(diff2)== -1) return -1;
  rs_wd.wd_fsa=diff2;
  rs_wd2.wd_fsa=diff2;
  if ((rfile = fopen(inf2,"r")) == 0) {
        fprintf(stderr,"Cannot open file %s.\n",inf2);
          exit(1);
  }
  tmalloc(new_diff2,fsa,1);
  int new_diff2_max=diff2->states->size * (DIFF2_INCREASE_FACTOR + 1);
  fsa_read(rfile,new_diff2,DENSE,0,new_diff2_max,TRUE,fsaname);
  fclose(rfile);
  if (fsa_table_dptr_init(new_diff2)== -1) return -1;
  //if (!add_diagonals && use_andnot)  
  if (extractwdsfromwa_command && !diff2name_command) {
	printf("please specify  -diff2name \n");
	exit(0);
  }
  if (exmin_command && !diff2name_command) {
	printf("please specify  -diff2name \n");
	exit(0);
  }
  if (!add_diagonals)  {
	if (diff2name_command) {
  		strcat(inf2x,diff2namestr); // gpname.diff2+diff2namestr
  		if ((rfile = fopen(inf2x,"r")) != 0) {
  			tmalloc(diff2diag,fsa,1);
  			fsa_read(rfile,diff2diag,DENSE,0,0,TRUE,fsaname);
  			fclose(rfile);
  			if (fsa_table_dptr_init(diff2diag)!= -1) {
  				rs_wd2.wd_fsa=diff2diag;
				Printf("%s used for reducing lhs words\n",inf2x);
			}
		}
		else {
        		fprintf(stderr,"Cannot open file %s.\n",inf2x);
          		exit(1);
		}
	}
  }
  reduce_word=diff_reducex;
  ngens = diff2->alphabet->base->size;
  diag_diff2=NULL;
  if (calculate_inverses(&inv,ngens,&rs_wd)==-1) return -1;
  if (add_diagonals)
  {
  	int diag_diff2_max=diff2->states->size * ngens;
	//int end1;
	//int end2;
	char inf2xx [100];
	if ((end_diagonals==999999)&&diff2name_command) {
  		strcpy(inf2xx,gpname); 
  		strcat(inf2xx,".diff2"); 
  		strcat(inf2xx,diff2namestr); 
	}
	else
		strcpy(inf2xx,inf2);
  	if ((rfile = fopen(inf2xx,"r")) == 0) {
       	 fprintf(stderr,"Cannot open file %s.\n",inf2xx);
          exit(1);
  	}
  	tmalloc(diag_diff2,fsa,1);
  	fsa_read(rfile,diag_diff2,DENSE,0,diag_diff2_max,TRUE,fsaname);
  	fclose(rfile);
  	if (fsa_table_dptr_init(diag_diff2)== -1) return -1;
  	add_diagonals_to_wd_fsa(diag_diff2,inv,&rs_wd,all_diagonals,check_diagonals,
		start_diagonals,end_diagonals,limit_diagonals,gpname,
		nobigger,diagnostics,min_wd_size,prevdiff2str);
  	Printf("\ndiff2 with diagonals has %d states\n",diag_diff2->states->size);
	if (diff2name_command) {
  		strcpy(inf2xx,".diff2"); 
  		strcat(inf2xx,diff2namestr); 
      		write_fsa (diag_diff2,gpname,inf2xx,fsaname);
	}
	else
	{
      		write_fsa (diag_diff2,gpname,".diff2d",fsaname);
	}
	if (read_wa)
	{
    		free_fsa(diag_diff2);
		exit(0);
	}
  }
  strcat(inf3,".andnot");
  strcpy(outf,inf3);
  strcpy(tempfilename,inf1);
  strcat(tempfilename,"temp_XXX");

  if (exmin_command ||extractwdsfromwa_command)
	read_wa=FALSE;
  if (!use_andnot && !read_wa && !exwa_command) {
      if (add_diagonals) {
        Printf("calling fsa_wa_x on %s.diff2 with diagonals added\n",gpname,hashlimit);
      	gpwa=fsa_wa_x(diag_diff2,op_store,tempfilename,FALSE);
    	free_fsa(diag_diff2);
      }
      else {
      	//gpwa=fsa_wa_x(diff2,op_store,tempfilename,FALSE);
	if (exmin_command) {
		char * fsatype=inf5;
		if (use_wa)
			fsatype=inf4;
		if (recovery)
      			gpwa2=fsa_extractwdsfrom_recovery(diff2,op_store,tempfilename,recovery_num);
		else
      			gpwa2=fsa_ex_min(diff2,op_store,tempfilename,fsatype,max_state,max_level,0); // .minred
		use_andnot=TRUE;
	}
	else if (extractwdsfromwa_command) {
		if (recovery)
      			gpwa2=fsa_extractwdsfrom_recovery(diff2,op_store,tempfilename,recovery_num);
		else
      			gpwa2=fsa_extractwdsfromwa(diff2diag,op_store,tempfilename,inf4,max_state,max_level,0); //.minred
		use_andnotx=TRUE;
	}
	else {
        	Printf("calling fsa_wa_x on %s.diff2\n",gpname);
      		gpwa=fsa_wa_x(diff2,op_store,tempfilename,FALSE,hashlimit);
	     }
      }
      if (exmin_command || extractwdsfromwa_command) {
	if (gpwa2==NULL) {
		Printf("No list of lhs words produced\n");
		exit(0);
	}
        Printf("  #Number of states of extractwdsfrom fsa  = %d.\n",
            gpwa2->states->size);
        write_fsa (gpwa2,gpname,".andnot",fsaname);
	free_fsa(gpwa2);
	// read wa
    	if ((rfile = fopen(inf4,"r")) == 0) {
       	 	fprintf(stderr,"Cannot open file %s.\n",inf4);
       	 	exit(1);
    	}
    	Printf("reading %s\n",inf4);
    	tmalloc(gpwa,fsa,1);
    	fsa_read(rfile,gpwa,DENSE,0,0,TRUE,fsaname);
    	Printf("%s has size %d\n",inf4,gpwa->states->size);
    	*wa_size=gpwa->states->size;
    	fclose(rfile);               
	//if (exmin_command)
	//	exit(0);
      }
      else { 
     
      if (kbm_print_level>1)
        printf("  #Number of states of gpwa before minimisation = %d.\n",
            gpwa->states->size);
      fsa_minimize_big_hash(gpwa,hash_table_size,counter,minimize_big_hash);
      if (kbm_print_level>1)
        printf("  #Number of states of gpwa after minimisation = %d.\n",
            gpwa->states->size);
      write_fsa (gpwa,gpname,".wa",fsaname);
      *wa_size=gpwa->states->size;
      if (do_waonly)
	 exit(0);	
      }
  }
  else {
    //-t, -exec -w or -geo switches
    // read existing wa file
    if (exwa_command && !read_wa) {
    	//char command [80];
    	//sprintf(command,"%s %s>cpcd",ex_str_wa,gpname);
	//if (add_diagonals) {
      	//	write_fsa (diag_diff2,gpname,".diff2d",fsaname);
	//}
    	Printf("calling %s\n",ex_str_wa);
    	int res=system(ex_str_wa);
    }
    if ((rfile = fopen(inf4,"r")) == 0) {
        fprintf(stderr,"Cannot open file %s.\n",inf4);
        exit(1);
    }
    Printf("reading %s\n",inf4);
    tmalloc(gpwa,fsa,1);
    fsa_read(rfile,gpwa,DENSE,0,0,TRUE,fsaname);
    Printf("%s has size %d\n",inf4,gpwa->states->size);
    *wa_size=gpwa->states->size;
    fclose(rfile);               
      if (do_waonly)
	 exit(0);	
  }
  if (do_beginswith) {
	printf("calling fsa_beginswith\n");
	gpwa_beginswith1=fsa_beginswith(gpwa,bwqual1);
  }
	

  if (use_andnot||use_andnotx) {
      // -t or -tx  switch
      // scan existing andnot file
      if ((rfile = fopen(inf3,"r")) != 0) {
          Printf("reading %s\n",inf3);
    	  tmalloc(fsaandnot_ptr,fsa,1);
          fsa_read(rfile,fsaandnot_ptr,DENSE,0,0,TRUE,fsaname);
          fclose(rfile);
      }
      else {
	printf("failed to open file %s\n",inf3);
	exit(1);
      }	
  }
  else if (do_geo) {
 	tmalloc(geowa,fsa,1);
	if (read_geowa) {
   if ((rfile = fopen(inf8,"r")) == 0) {
        fprintf(stderr,"Cannot open file %s.\n",inf8);
        exit(1);
    }
    Printf("reading %s\n",inf8);
    fsa_read(rfile,geowa,DENSE,0,0,TRUE,fsaname);
    fclose(rfile);
  }
  else {
	Printf("calling fsa_wa_x on %s.diff2 (geo option)\n",gpname);
	geowa=fsa_wa_x(diff2,DENSE,tempfilename,TRUE,0);
        if (kbm_print_level>1)
                printf(
                "  #Number of states of geowa before minimisation = %d.\n",
                         geowa->states->size);
        fsa_minimize(geowa);
        if (kbm_print_level>1)
        printf("  #Number of states of geowa after minimisation = %d.\n",
                geowa->states->size);
  }
	fsaandnot_ptr=find_geo_wds(diff2,gpwa,gpname, tempfilename,
			geowa,hash_table_size,counter,state_limit);
        if (!no_exists_minimize) {
      if (kbm_print_level>1)
        printf("  #Number of states of fsaandnot before minimisation = %d.\n",
            fsaandnot_ptr->states->size);
      fsa_minimize_big_hash(fsaandnot_ptr,hash_table_size,counter,minimize_big_hash);
      if (kbm_print_level>1)
        printf("  #Number of states of fsaandnot after minimisation = %d.\n",
            fsaandnot_ptr->states->size);
	}
  }
  else {
      if (look_for_diff2_wds)
         // -p or -rp  switches
      {
	  if (!dogeowds)
            Printf("calling fsa_pfxred on %s\n",inf4);
	  if (do_beginswith) // never used!
	  {
		printf("using gpwa_beginswith1 as lhs \n");
          	pfxred = fsa_pfxred(gpwa_beginswith1,op_store,
				tempfilename,0);
	  }
	  else
	  {
		if (!dogeowds) 
          		pfxred = fsa_pfxred(gpwa,op_store,tempfilename,state_limit);
		if (state_limit > gpwa->states->size)
		{
			state_limit = 0;
			dogeowds=FALSE;
		}
		if (state_limit) {
			Printf("wa state limit is %d\n",state_limit);
			if (!dogeowds) {
              		if (kbm_print_level>1)
                		printf(
			"#Number of states of pfxred before minimisation = %d.\n",
                    		pfxred->states->size);
              		if (fsa_minimize_big_hash(pfxred,hash_table_size,
				counter,minimize_big_hash)==-1) exit(1);
              		if (kbm_print_level>1)
                		printf(
				"#Number of states of pfxred after minimisation = %d.\n",
                    		pfxred->states->size);
			}
		}
	  }
          //pfxred= wa with failstate converted to accept
	  if (!dogeowds)
             if (pfxred==0) exit(1);
          if (!no_filter)
          {
              if (use_big_hash)
              {
		 if (!dogeowds) {
                  Printf(
		   "calling fsa_triples_big_hash on pfxred, gpwa and %s\n",inf2);
                  minkb =fsa_triples_big_hash (pfxred,gpwa,diff2,op_store2, FALSE,
				tempfilename,hash_table_size,counter, restricted_pfxred,state_limit);
		 }
		  if (do_mult2) {
   fsa *genmult2ptr = fsa_genmult2_int_x(minkb,op_store);
      if (genmult2ptr==0) exit(1);
      if (kbm_print_level>1)
        printf(
      "  #Number of states of genmult2 = %d.\n",genmult2ptr->states->size);
	exit(99);
   }
		  if (diff2c) {
		/*strcpy(outf,gpname);
      strcat(outf,".pfxred");
      base_prefix(fsaname);
      strcat(fsaname,".pfxred");
      wfile = fopen(outf,"w");
      Printf("\nwriting %s\n",outf);
      fsa_print(wfile,pfxred,fsaname);
      fclose(wfile);*/
              if (fsa_minimize_big_hash(minkb,hash_table_size,
			counter,minimize_big_hash)==-1) exit(1);
      printf("converting triples to DENSE format\n");
      write_fsa (minkb,gpname,".triples",fsaname);
      strcpy(outf,gpname);
      strcat(outf,".triples");
      rfile = fopen(outf,"r");
      fsa_read(rfile,minkb,DENSE,0,0,TRUE,fsaname);

			printf("calculating diff2c\n");
			fsa *diff2c=fsa_diff(minkb, &rs_wd, SPARSE);
			printf("diff2 has actual size %d\n",diff2c->states->size);
	printf("make full wd fsax\n");	
  	make_full_wd_fsax(diff2c,inv,0,&rs_wd);
  fclose(rfile);
   write_fsa (diff2c,gpname,".diff2c",fsaname);
	exit(99);
		  }
              }
              else
              {
                  Printf("calling fsa_minkb on pfxred, gpwa and %s\n",inf2);
                  minkb =fsa_minkb (pfxred,gpwa,diff2,op_store2,
                                FALSE,tempfilename);
              }
	      if (!dogeowds) { 
              if (minkb==0) exit(1);
              if (kbm_print_level>1)
                printf(
		"#Number of states of triples before minimisation = %d.\n",
                    minkb->states->size);
	      }
// cpcdaft
if (exminex_command) {
        write_fsa(minkb,gpname,".triples",fsaname);
	free_fsa(minkb);
        minkb=NULL;
        Printf("calling %s\n",ex_str_minex);
        int res=system(ex_str_minex);
	//exit(99);
      strcpy(outf,gpname);
      strcat(outf,".triples.min.exists");
      rfile = fopen(outf,"r");
  tmalloc(fsaexists,fsa,1);
Printf("reading %s\n",outf);
fsa_read(rfile,fsaexists,DENSE,0,0,TRUE,fsaname);
  unlink(outf);
               if (kbm_print_level>1)
                 printf("  Number of states of fsaexists  = %d.\n",
                    fsaexists->states->size);
}
else {
	      if (!dogeowds)
	      {
              if (fsa_minimize_big_hash(minkb,hash_table_size,
			counter,minimize_big_hash)==-1) exit(1);
              //if (fsa_minimize(minkb)==-1) exit(1);
              if (kbm_print_level>1)
                printf(
		"#Number of states of triples after minimisation = %d.\n",
                    minkb->states->size);
	      }
	      if (dogeowds&&(state_limit>0 && state_limit<gpwa->states->size))
              {
		/*printf("create diff2T from minkb\n");
		fsa * diff2T;
      		printf("converting triples to DENSE format\n");
      		write_fsa (minkb,gpname,".triples",fsaname);
      		strcpy(outf,gpname);
      		strcat(outf,".triples");
      		rfile = fopen(outf,"r");
      		fsa_read(rfile,minkb,DENSE,0,0,TRUE,fsaname);
      		fclose(rfile);
		printf("deleting %s\n",outf);
  		unlink(outf);
		printf("calculating diff2T\n");
		diff2T=fsa_diff(minkb, &rs_wd, SPARSE);
		printf("diff2T has  size %d\n",diff2T->states->size);
		free_fsa(minkb);
		minkb=NULL;
		printf("make full wd fsax\n");	
  		make_full_wd_fsax(diff2T,inv,1,&rs_wd);
		*/
		//printf("creating geowa on diff2T");
		//geowa=fsa_wa_x(diff2T,DENSE,tempfilename,TRUE);
	if (read_geowa) {
 	tmalloc(geowa,fsa,1);
   if ((rfile = fopen(inf8,"r")) == 0) {
        fprintf(stderr,"Cannot open file %s.\n",inf8);
        exit(1);
    }
    Printf("reading %s\n",inf8);
    fsa_read(rfile,geowa,DENSE,0,0,TRUE,fsaname);
    fclose(rfile);
  }
  else {
		printf("creating geowa on diff2");
		geowa=fsa_wa_x(diff2,DENSE,tempfilename,TRUE,0);
        if (kbm_print_level>1)
                printf(
                "  #Number of states of geowa before minimisation = %d.\n",
                         geowa->states->size);
        fsa_minimize(geowa);
        if (kbm_print_level>1)
        printf("  #Number of states of geowa after minimisation = %d.\n",
                geowa->states->size);
      write_fsa (geowa,gpname,".geowa",fsaname);
  }
        //fsaandnot_ptr=find_geo_wds(diff2T,gpwa,gpname, tempfilename,
                        //geowa,hash_table_size,counter,state_limit);
        fsaandnot_ptr=find_geo_wds(diff2,gpwa,gpname, tempfilename,
                        geowa,hash_table_size,counter,state_limit);
		//free_fsa(diff2T);
		//diff2T=NULL;
		free_fsa(geowa);
		geowa=NULL;
		//exit(99);
	      }
	  else {
              Printf("calling fsa_exists on triples\n");
              fsaexists = fsa_exists(minkb,op_store,TRUE,tempfilename);
              if (fsaexists==0) exit(1);
	     if (!no_exists_minimize)
	     {
              if (kbm_print_level>1)
                printf("  Number of states of fsaexists before minimisation = %d.\n",
                    fsaexists->states->size);
               if (fsa_minimize_big_hash(fsaexists,hash_table_size,
			counter,minimize_big_hash)== -1) exit(1);
               if (kbm_print_level>1)
                 printf("  Number of states of fsaexists after minimisation = %d.\n",
                    fsaexists->states->size);
	     }
	    } // not dogeowds
	    } // not exec_minex
	      if (!dogeowds) {
              Printf("fsa_and_not on pfxred and fsaexists\n");
	      if (pfxred->states->size == fsaexists->states->size)
		 minimize_andnot=TRUE;
              //fsaandnot_ptr = fsa_and_not(pfxred,fsaexists,op_store,TRUE,tempfilename);
              fsaandnot_ptr = fsa_and_not_first(pfxred,fsaexists,op_store,TRUE,tempfilename,state_limit);
	     }
         }
      }
      else
      {
          // looking for diff1 wds via minred
          minred=NULL;
          if (!read_minred&&!read_minred1)
          {
              Printf("calling fsa_minred on gpwa\n");
              minred = fsa_minred(gpwa,op_store,FALSE,tempfilename);
              //minred = fsa_minredx(gpwa,op_store,FALSE,tempfilename,state_limit);
	      //state_limit = 0;
              if (minred==0) exit(1);
              Printf("  #Number of states of minred before minimisation = %d.\n",
                    minred->states->size);
	      if (!no_minred_minimize)
	      {
              if (fsa_minimize_big_hash(minred,hash_table_size,
			counter,minimize_big_hash)==-1) exit(1);
              Printf("  #Number of states of minred after minimisation = %d.\n",
                    minred->states->size);
             // write minred to file
	      write_fsa (minred,gpname,".minred",fsaname);
	      if (do_minredonly)
		exit(0);
	      if (state_limit && state_limit < gpwa->states->size) {
		Printf("calling fsa_pfxred with state_limit %d\n",state_limit);
		pfxred = fsa_pfxred(gpwa,op_store,tempfilename,state_limit);
		Printf("anding pfxred with minred, then minimise\n");
		minred=fsa_and(pfxred,minred,op_store, TRUE,tempfilename);
		pfxred = NULL;
                if (fsa_minimize_big_hash(minred,hash_table_size,
			counter,minimize_big_hash)==-1) exit(1);
                Printf("  #Number of states of minred after anding with pfxred and minimisation = %d.\n",
                    minred->states->size);
		
	        state_limit = 0;
	      }
	      }
          }
          else {
            // -min or min1 switch  
	    char *fname;
	    if (read_minred)
		fname=inf5;
	    else
		fname=inf7;
            if ((rfile = fopen(fname,"r")) == 0) {
                fprintf(stderr,"Cannot open file %s.\n",fname);
                exit(1);
            }
            Printf("reading %s\n",fname);
            tmalloc(minred,fsa,1);
            fsa_read(rfile,minred,DENSE,0,0,TRUE,fsaname);
            fclose(rfile);               
          }
	    if (do_minred1) {
		// attempt to find diff2 wds quicker that gpcheckmult can
		Printf("calling fsa_minred1 using minred and gpwa\n");
		minred1=fsa_minred1(minred,gpwa,DENSE, FALSE,tempfilename);
		fsa_clear(minred);
		free(minred);
		minred=minred1;
                if (kbm_print_level>1)
                 printf(
  		"  #Number of states of minred1 before minimisation = %d.\n",
			minred1->states->size);
	      if (!no_minred_minimize){
		if (fsa_minimize(minred)==-1) exit(1);	
                if (kbm_print_level>1)
                 printf(
  		"  #Number of states of minred1 after minimisation = %d.\n",
			minred1->states->size);
 		// write minred to file
 		 write_fsa (minred,gpname,".minred",fsaname);
		 minred1=NULL; // minred1 is now minred!
	      }
           }
          if (!no_filter)
          {
	      if (!read_minkb)
	      {
              if (use_big_hash)
              {
		if (state_limit > minred->states->size)
			state_limit=0;
		if (state_limit)
			Printf("minred state limit is %d\n",state_limit);
                  Printf(
		   "calling fsa_triples_big_hash on minred, gpwa and %s\n",inf2);
                  minkb =
		   fsa_triples_big_hash (minred,gpwa,diff2,op_store2,FALSE,
                            tempfilename,hash_table_size,
				counter,!do_minred1,state_limit);
           if (dump_files) {
		write_fsa (minred,gpname,".minred",fsaname);
		write_fsa (minkb,gpname,".triples",fsaname);
        	exit(99);
                  }
              }
              else
              {
                  // -ht 0
                  Printf("calling fsa_minkb on minred, gpwa and %s\n",inf2);
                  minkb =fsa_minkb (minred,gpwa,diff2,op_store2,
                                FALSE,tempfilename);
              }
              if (minkb==0) exit(1);
              if (kbm_print_level>1)
                printf("  #Number of states of minkb before minimisation = %d.\n",
                    minkb->states->size);
              //if (fsa_minimize(minkb)==-1) exit(1);
              if (fsa_minimize_big_hash(minkb,hash_table_size,
			counter,minimize_big_hash)==-1) exit(1);

              if (kbm_print_level>1)
                printf("  #Number of states of minkb after minimisation = %d.\n",
                    minkb->states->size);
              }  
	      else
	      {
		// read minkb
                // -minkb` switch
	            if ((rfile = fopen(inf6,"r")) == 0) {
       		         fprintf(stderr,"Cannot open file %s.\n",inf6);
       		         exit(1);
       		     }
       		     Printf("reading %s\n",inf6);
       		     tmalloc(minkb,fsa,1);
       		     fsa_read(rfile,minkb,ip_store,dr,0,TRUE,fsaname);
       		     fclose(rfile);
	      }
              Printf("calling fsa_exists on minkb\n");
              fsaexists = fsa_exists(minkb,op_store,TRUE,tempfilename);
              if (fsaexists==0) exit(1);

	      if (!no_exists_minimize)
	      {		
              if (kbm_print_level>1)
                printf("  Number of states of fsaexists before minimisation = %d.\n",
                    fsaexists->states->size);
              if (fsa_minimize_big_hash(fsaexists,hash_table_size,
			counter,minimize_big_hash)== -1) exit(1);
              if (kbm_print_level>1)
                printf("  Number of states of fsaexists after minimisation = %d.\n",
                    fsaexists->states->size);
	      }	
              Printf("fsa_and_not on minred and fsaexists\n");
              if (minred->states->size == fsaexists->states->size)
                 minimize_andnot=TRUE;
              //fsaandnot_ptr = fsa_and_not(minred,fsaexists,op_store,TRUE,tempfilename);
             fsaandnot_ptr = fsa_and_not_first(minred,fsaexists,op_store,TRUE,tempfilename,state_limit);
          }
      }
      if (!no_filter)
      {
 
          if (fsaandnot_ptr->states->size == 0)
	  {
		Printf("empty fsa\n");
		exit(0);
	  }
          if (kbm_print_level>1)
            printf(
		"  #Number of states of fsaandnot before minimisation = %d.\n",
                	fsaandnot_ptr->states->size);
	  if (minimize_andnot)
	  {
	          if (fsa_minimize(fsaandnot_ptr)== -1) exit(1);
	          if (kbm_print_level>1)
	     		  printf(
     "  #Number of states of fsaandnot after minimisation = %d.\n",
	            fsaandnot_ptr->states->size); 
          	  if (fsaandnot_ptr->states->size == 0) exit(0);
	  }
      }
      else
          // -nf switch
      {
	  // avoid timing out too soon if wa and diff2 are inconsistent
	  if (!timeout)
		timeout=120000;
          if (look_for_diff2_wds)
          {
              fsaandnot_ptr = pfxred;
	      pfxred=NULL;
          }
          else {
              fsaandnot_ptr = minred;
              minred=NULL;
	  }
      }

  }
  if (!use_andnot && !use_andnotx)
    write_fsa(fsaandnot_ptr,gpname,".andnot",fsaname);
  if (no_scan) {
    printf("no scan\n");
  return 0;
  }	
  old_ndiff = new_diff2->states->size;
  kbm_print_level--;
  //ngens = diff2->alphabet->base->size;
  //reduce_word=diff_reduce;
  //if (calculate_inverses(&inv,ngens,&rs_wd)==-1) return -1;
  int pos_max=0;
  signal(SIGINT, interrupt_gpcheckx);
  signal(SIGKILL, interrupt_gpcheckx);
  signal(SIGQUIT, interrupt_gpcheckx);
  int old_diff_size = new_diff2->states->size;
  int wordnum=process_words(fsaandnot_ptr,&rs_wd,&rs_wd2,start_scan_from,
                gpwa,new_diff2,inv,no_dots,scale,
		minwordlen,maxwordlen,trace_equations,
		&pos_max,use_alt_prestate,verify,verify_qualifier,
		prefixby1,timeout,use_andnotx,line_items);
  kbm_print_level++;
  if (update_diff2)
  {
  	make_full_wd_fsa(new_diff2,inv,old_ndiff+1,&rs_wd);
  	if (kbm_print_level>=1)
	      printf("  #Word-difference machine now has %d states.\n",
               new_diff2->states->size);
      base_prefix(fsaname);
      strcat(fsaname,".diff2");
      wfile = fopen(inf2,"w");
      fsa_print(wfile,new_diff2,fsaname);
      fclose(wfile);
  }
  else
      printf("\ndiff2 not updated\n");
  int new_diff_size =
               new_diff2->states->size;
  free(inv);
  free_fsa(gpwa);
  free_fsa(minred);
  free_fsa(minred1);
  free_fsa(pfxred);
  free_fsa(minkb);
  free_fsa(fsaexists);
  free_fsa(fsaandnot_ptr);
  if (do_geo&&!read_geowa)
  {
	write_fsa (geowa,gpname,".geowa",fsaname);
	//printf("call gpgeowa to compute geodiff and geopairs\n");
  	if ((new_diff_size - old_diff_size) == 0)
		printf ("%s.geowa should be correct\n",gpname);
		
  }
  if (do_geo && CALC_GEODIFF)
  {
  char outfx [100];
  char outfy [100];
  fsa *geopairsptr, *geodiffptr;
  tmalloc(geopairsptr,fsa,1);
  strcpy(outfx,gpname);
  strcat(outfx,".geopairs");
  strcpy(outfy,gpname);
  strcat(outfy,".geodiff");
  rfile=fopen(outfx,"r");
  fsa_read(rfile,geopairsptr,DENSE,0,0,TRUE,fsaname);
  fclose(rfile);
  geodiffptr=fsa_diff(geopairsptr,&rs_wd,op_store2);
  fsa_clear(geopairsptr);
  tfree(geopairsptr);

  base_prefix(fsaname);
  strcat(fsaname,".geodiff");
  wfile = fopen(outfy,"w");
  fsa_print(wfile,geodiffptr,fsaname);
  fclose(wfile);
  if (kbm_print_level>0)
    printf("#Geodesic difference machine with %d states computed.\n",
            geodiffptr->states->size);
  fsa_clear(geodiffptr);
  tfree(geodiffptr);

  }
  free_fsa(diff2);
  free_fsa(new_diff2);
  free_fsa(geowa);
  if (dogeowds && (state_limit==0))
	if (new_diff_size == old_diff_size)
		exit(99);
  if (line_items) {
        //Printf("reached %d\n",wordnum);
	return wordnum;
  }
  else
  	return new_diff_size - old_diff_size;;
}

void free_fsa (fsa * fsa_to_free)
{
	if (fsa_to_free) {
		fsa_clear(fsa_to_free);
		free(fsa_to_free);
	}
}

fsa *fsa_pfxred(waptr,op_table_type,tempfilename,state_limit)
        fsa *waptr;
        storage_type op_table_type;
        char *tempfilename;
	int state_limit;
{
 /* This routine uses a longwinded method to produce an fsa which is simply
    the wa with its fail state replaced by an accept state! */

  int **table, ne, nsi, nsi1, ns, dr, *fsarow, nt, cstate, csa, csb, im, i, g,
      len = 0, ct, *ht_ptr;
  boolean dense_ip, dense_op;
  fsa *pfxred;
  hash_table ht;
  FILE *tempfile;

  if (kbm_print_level >= 3)
    printf("    #Calling fsa_pfxred.\n");
  if (!waptr->flags[DFA]) {
    fprintf(stderr, "Error: fsa_pfxred only applies to DFA's.\n");
    return 0;
  }

  tmalloc(pfxred, fsa, 1);
  fsa_init(pfxred);
  srec_copy(pfxred->alphabet, waptr->alphabet);
  pfxred->flags[DFA] = TRUE;
  pfxred->flags[ACCESSIBLE] = TRUE;
  pfxred->flags[BFS] = TRUE;

  ne = waptr->alphabet->size;
  nsi = waptr->states->size;
  nsi1 = nsi + 1;
  if (state_limit)
	nsi1=state_limit+1;
	

  table = waptr->table->table_data_ptr;

  dense_ip = waptr->table->table_type == DENSE;
  dr = waptr->table->denserows;
  dense_op = op_table_type == DENSE;

  pfxred->states->type = SIMPLE;
  pfxred->num_initial = 1;
  tmalloc(pfxred->initial, int, 2);
  pfxred->initial[1] = 1;
  pfxred->table->table_type = op_table_type;
  pfxred->table->denserows = 0;
  pfxred->table->printing_format = op_table_type;

  hash_init(&ht, TRUE, 1, 0, 0);
  ht_ptr = ht.current_ptr;
  ht_ptr[0] = waptr->initial[1];
  im = hash_locate(&ht, 1);
  if (im != 1) {
    fprintf(stderr, "Hash-initialisation problem in fsa_pfxred.\n");
    return 0;
  }
  if ((tempfile = fopen(tempfilename, "w")) == 0) {
    fprintf(stderr, "Error: cannot open file %s\n", tempfilename);
    return 0;
  }
  if (dense_op)
    tmalloc(fsarow, int, ne) else tmalloc(fsarow, int, 2 * ne + 1)

        cstate = 0;
  if (dense_op)
    len = ne; /* The length of the fsarow output. */
  nt = 0;     /* Number of transitions in pfxred */

  /* A state of *pfxred will be essentially a  state of *waptr,
   * with the transitions coming from those of *waptr.
   * The difference is that the new state will accept words  w
   * which are not accepted by *waptr but whose maximal prefix is,
   * So transitions to 0 in *waptr will go to a new accept-
   * state nsi1 instead (with no transitions from nsi1).
   * The initial state itself is non-accept.
   */
  while (++cstate <= ht.num_recs) {
    if (kbm_print_level >= 3) {
      if ((cstate <= 1000 && cstate % 100 == 0) ||
          (cstate <= 10000 && cstate % 1000 == 0) ||
          (cstate <= 100000 && cstate % 5000 == 0) || cstate % 50000 == 0)
        printf("    #cstate = %d;  number of states = %d.\n", cstate,
               ht.num_recs);
    }
    ht_ptr = hash_rec(&ht, cstate);
    csa = ht_ptr[0];
    if (!dense_op)
      len = 0;
    for (g = 1; g <= ne; g++) {
      /* Calculate action of generator g on state cstate */
      ht_ptr = ht.current_ptr;
      if (csa == 0 || csa == nsi1)
        ht_ptr[0] = 0;
      else {
        ht_ptr[0] = target(dense_ip, table, g, csa, dr);
        if (ht_ptr[0] == 0)
          ht_ptr[0] = nsi1;
	else {
      	  if (state_limit)
		if (ht_ptr[0]>state_limit)
			ht_ptr[0]=0;
	}
      }
      if (ht_ptr[0] == 0 )
        im = 0;
      else
        im = hash_locate(&ht, 1);
      if (dense_op)
        fsarow[g - 1] = im;
      else if (im > 0) {
        fsarow[++len] = g;
        fsarow[++len] = im;
      }
      if (im > 0)
        nt++;
    }
    if (!dense_op)
      fsarow[0] = len++;
    fwrite((void *)fsarow, sizeof(int), (size_t)len, tempfile);
  }
  fclose(tempfile);

  ns = pfxred->states->size = ht.num_recs;
  pfxred->table->numTransitions = nt;
  /* Locate the accept states: first count them and then record them. */
  ct = 0;
  for (i = 1; i <= ns; i++) {
    ht_ptr = hash_rec(&ht, i);
    if (ht_ptr[0] == nsi1)
      ct++;
  }
  pfxred->num_accepting = ct;
  if (ct == 1 || ct != ns) {
    tmalloc(pfxred->accepting, int, ct + 1);
    ct = 0;
    for (i = 1; i <= ns; i++) {
      ht_ptr = hash_rec(&ht, i);
      if (ht_ptr[0] == nsi1)
        pfxred->accepting[++ct] = i;
    }
  }
  hash_clear(&ht);
  tfree(fsarow);
  tfree(waptr->is_accepting);

  /* Now read the transition table back in */
  tempfile = fopen(tempfilename, "r");
  compressed_transitions_read(pfxred, tempfile);
  fclose(tempfile);

  unlink(tempfilename);

  return pfxred;
}
 
void badusage_gpcheckx(overview)
	boolean overview;
{
    fprintf(stderr,
	"\nUSAGE: gpcheckx [-p][-m][-nf] [-s state-limit][-geo][-to timeout][-w][-gw][-v][-ve] [-diagonals] [-diff2name suffix] gpname [+rptz|+wrptz]\n\n");
    if (!overview) {
    	fprintf(stderr,
	"Type 'gpcheckx -ovw' for more information\n");
    }
    else {
      fprintf(stderr,
	"\nOVERVIEW\n"
	" gpcheckx builds, scans and processes a list of left hand side (lhs) words, checking these for the presence of new word differences.\n" 
	" The objective is to build a correct word acceptor, gpname.wa, where gpname is an automatic group.\n"
	" The lhs word list can be built as non-filtered (-nf), prefix reduced (-p) or  minimally reducible (-m). See DEFINITIONS below\n" 
	" The lhs word list can be further restricted by specifying a state number, 'state-limit', of gpname.wa (-s state-limit).\n"
	" This means that the maximal prefix of any lhs word traces a path of states in gpname.wa and every state in this path should be less\n"
	" than 'state-limit'.\n"
	" Format of 'state-limit' is s-value|base;increment|base;increment;maximum. If 'base' and 'increment', but not 'maximum' are  specified the command will repeat\n"
	" with the 's-value' initially set to 'base', then increased by increment each time. The command will repeat until the wa size is less than the -s value\n"
	" and no new word differences have been found. If 'maximum' is specified, then the command will repeat until 's-value' is greater than 'maximum'.\n" 
	" gpcheckx can be run in geodesic mode (-geo). The objective, then, is to build a correct geodesic word acceptor, gpname.geowa.\n"  
	" Scanning of the lhs word list ends when either, all words are processed, a control & C signal is received or is timed out after\n" 
	" 'timeout' seconds (specified by -to timeout).\n"
	" gpcheckx can read in a previously built gpname.wa (-w) or gpname.geowa (-gw) instead of constructing it from scratch.\n" 
	" The output can be verbose,'-v', or verbose with equations traced,'-ve'\n"  
	" The program can be repeated (+rptz or +wrptz) with the same parameters until no new word differences can be found.\n"
	" +wrptz means the first run of a repeated action will read the wa from gpname.wa instead of constructing it.\n"
        " -diagonals causes 'diagonal' word differences to be temporarily added to the word difference list. This increases the chances of building a correct\n"
	" word acceptor at the first attempt, but at the cost of a much longer processing time.\n\n"
	"DEFINITIONS\n"
	" Reduced word. A word that is accepted by the word acceptor, gpname.wa.\n"
	" Prefix reduced word. A word that is reducible but all of its prefixes are reduced.\n"
	" Minimally reducible. A word that is reducible but whose every subword is reduced.\n"
	" left hand side word (lhs). A word that is prefix reduced. It reduces to a right hand side (rhs).\n"
	" Non-filtered word list. A word list that is derived directly from the fsa that accepts minimally reducible words, gpname.minred \n"
	" Diagonal word difference. A word obtained from a 'proper' word difference by multiplying on the left by a generator\n");
     }
     exit (99);
}

int process_words (fsaptr,rs_wd,rs_wd2,start_scan_from,gpwa,new_diff2,inv,
                   no_dots,scale,minwordlen, maxwordlen,trace_equations,
		   pos_max,use_alt_prestate,verify,verify_qualifier,
			prefixby1,timeout,andnotx,line_items)
	fsa *fsaptr;
	reduction_struct *rs_wd;
	reduction_struct *rs_wd2;
	int start_scan_from;
	fsa *gpwa;
	fsa *new_diff2;
	int *inv;
	int no_dots;
	unsigned int scale;
	int minwordlen;
	int maxwordlen;
	boolean trace_equations;
	int *pos_max;
	boolean use_alt_prestate;
	boolean verify;
	int verify_qualifier;
	boolean prefixby1;
	int timeout;
	boolean andnotx; //use fsaptr as  wa for onelevelreduced
        int line_items;
{
 // partly converted from Alun Williams MAF code

  int i, g1, g2, ne, ngens, ngens1, **table, dr, pass, found, si,
       firste, clength, clength1, clength2, cstate, im, *statelist;
  
  boolean dense, done, backtrack, foundword,printnow;
  gen  *rword,bufi,*prev;
  //int aim_state = fsaptr->accepting[1];
  int nr_states = fsaptr->states->size,letter;
  gen workbuf1[MAX_WORD_LEN];
  gen workbuf2[MAX_WORD_LEN];
  reduction_equation eqn;
  int now=3,hit=0,sincewin=0;
  gen lhs_word[MAX_WORD_LEN];
  gen rhs_word[MAX_WORD_LEN];

  fsa *diff2=rs_wd->wd_fsa;
  fsa *diff2_copy, *diff2_copy2;
  fsa *wa_for_1_level;
  int ***table_dptr, *prestates, *prestates2, *poststates, *distances;
  gen  *preletters, *preletters2,*postletters;
  int diff_size=diff2->states->size;
  gen *new_states; // contains word difference values
  int new_states_number [diff_size * DIFF2_INCREASE_FACTOR];
  boolean  *succeeding,*accepting_states, *duplicates, *statesdone; 
  gen null_char = 0x0;
  int count_dots=0;
  boolean first_time_display=TRUE;
  int dotscrosses = 0;
  
  if (fsaptr->num_initial==0)
  {
	fsa_clear(fsaptr);
        return -1;
  }
  if (andnotx)
	wa_for_1_level=fsaptr;
  else
	wa_for_1_level=gpwa; 

  int sample;
  int newdiff;
  int no_x=0;
  unsigned int *block1;
  unsigned int *block2;
  //int accept_state=fsaptr->accepting[1];
  int size_of_dot=DEFAULT_SIZE_OF_DOT;
  int maximum_length=1;
  if (scale < 7 && scale > 0)
  {
	int ii=0;
	size_of_dot=1;
	while (ii<scale)
	{
		ii++;
		size_of_dot*=10;
	}
  } 

 // Building prefix and postfix structures
  PPrintf("andnot fsa has %d states.\n",nr_states);
  PPrintf("building postfix and prefix structures\n"); 
  ngens = fsaptr->alphabet->size;

  dr = fsaptr->table->denserows;
  dense = fsaptr->table->table_type==DENSE;
  table = fsaptr->table->table_data_ptr;
  tmalloc(prestates,int,nr_states+1);
  tmalloc(preletters,gen,nr_states+1);
  if (use_alt_prestate) {
  	tmalloc(prestates2,int,nr_states+1);
  	tmalloc(preletters2,gen,nr_states+1);
  	memset(prestates2,0,sizeof(int)*(nr_states+1));
  	memset(preletters2,0,sizeof(gen)*(nr_states+1));
  }
  tmalloc(poststates,int,nr_states+1);
  tmalloc(postletters,gen,nr_states+1);
  tmalloc(distances,int,nr_states+1);
  tmalloc(succeeding,boolean,nr_states+1);
  tmalloc(accepting_states,boolean,nr_states+1);
  memset(succeeding,FALSE,sizeof(boolean)*(nr_states+1));
  memset(distances,0,sizeof(int)*(nr_states+1));
  memset(prestates,0,sizeof(int)*(nr_states+1));
  memset(preletters,0,sizeof(gen)*(nr_states+1));
  memset(poststates,0,sizeof(int)*(nr_states+1));
  memset(postletters,0,sizeof(gen)*(nr_states+1));
  memset(succeeding,FALSE,sizeof(boolean)*(nr_states+1));
  memset(accepting_states,FALSE,sizeof(boolean)*(nr_states+1));
  int acci=1;
  while (acci <= fsaptr->num_accepting)
  {
	succeeding[fsaptr->accepting[acci]]=TRUE;
	accepting_states[fsaptr->accepting[acci]]=TRUE;
	acci++;
  }

 // for each state > 1 find the lowest valued state that can
 // transition to the state and place this in the prestates entry
 // for that state
  found=1;	
  for (si = 1; si <= nr_states;si++)
  {
      int ti;
      for (ti = 1; ti <= ngens; ti++)
      {
          int ns = table[ti][si];
	  if (ns>nr_states)
		ns=nr_states; //pointer to an accept state
          if (ns > 0) {
              if (prestates[ns] == 0)
              {
                  prestates[ns] = si;
		  preletters[ns] = ti;
                  found++;
                  if (found >= nr_states)
                    break;
              }
	      else if (use_alt_prestate){
		// make this an alternative prestate
		prestates2[ns] = si;
		preletters2[ns] = ti;
	      }
          }    
      }
  }
 // find the shortest path from each state to an accept state -
 // distances contains, for each state, the length of the shortest path 
 // to an accept state 
  pass = 0;
  found=1;
  done=TRUE;
  while (done && found < nr_states)
  {
    done = FALSE;
    pass++;
    for (si = 1; si <= nr_states;si++)
    {
      int ti;
      if (!succeeding[si])
      {
        for (ti = 1; ti <= ngens; ti++)
        {
          int ns = table[ti][si];
	  if (ns>nr_states)
		ns=nr_states; //pointer to an accept state
          if  (succeeding[ns] && distances[ns]+1==pass)
          {
            distances[si] = pass;
	    poststates[si]=ns;
	    postletters[si]=ti;
            succeeding[si] = TRUE;
            found++;
            done = TRUE;
            break;
          }
        }
      }
    }
  }
  tfree(succeeding);
  tfree(distances);
  if (!andnotx)
  	fsa_clear(fsaptr);
  table=NULL;

  tmalloc(statesdone,boolean,nr_states+1);
  tmalloc(duplicates,boolean,nr_states+1);
  memset(statesdone,FALSE,sizeof(boolean)*(nr_states+1));
  memset(duplicates,FALSE,sizeof(boolean)*(nr_states+1));
  tfree(statesdone);
  // make table_dptr a copy of diff2->table->table_data_dptr[ll][rr][s]
  // with  extra room for new states and transitions. This will be
  // updated as new wd's are extracted from the lhs words
  table_dptr=copy_table_dptr (diff2);

 // determine how many dots will be displayed
  if (start_scan_from > nr_states)
	start_scan_from=1;
  if (line_items > 0)
	sample=1;
  else if (no_dots < 0) {
      int total_states = no_dots * -size_of_dot;
      sample=1;
      if (total_states < (nr_states-start_scan_from))
        nr_states=total_states + start_scan_from; 
  }
  else
      sample=(nr_states-start_scan_from)/(no_dots*size_of_dot) + 1;

  if (sample>size_of_dot - 1)
     sample=size_of_dot - 1;
  if (sample <= 0)
	sample=1;

  if (verify)
	printf("verifying..\n");
else
	PPrintf("\nprocessing %d words from %d, sample=%d, .=%d words,x new wd",
         nr_states-start_scan_from,start_scan_from,sample,size_of_dot);
	PPrintf("\n(do control & C to interrupt)\n");

  block1=malloc(MAXPAIRS*2*sizeof(unsigned int));
  block2=malloc(MAXPAIRS*2*sizeof(unsigned int));
  // new_states [diff_size * DIFF2_INCREASE_FACTOR][MAXWDLENGTH];
  new_states = malloc(diff_size * DIFF2_INCREASE_FACTOR * MAXWDLENGTH *sizeof(gen));
  new_states[0]=null_char;
  int ignored_duplicates=0;
  int processed_words=0;
  unsigned int total_length=0;
  int** watable=gpwa->table->table_data_ptr;

  int fails=0;
  ne = gpwa->alphabet->size;

  clock_t start_t, end_t, total_t;
  start_t=clock();

  for (si=start_scan_from;si<nr_states;si+=sample) {

      int i,si2=si,ti,j,tix,si2x;
      if (line_items>0 && dotscrosses>line_items)
      {
	PPrintf(" interrupted before %d\n",si);
	break;
      }
      end_t=clock();
      total_t=(end_t - start_t) / CLOCKS_PER_SEC;
      if (timeout)
	      if (total_t>timeout) {
		 PPrintf("\ntimeout of scan at word %d\n",si);
		 si=nr_states;
		 break; 
	      }
      if (gpcheckx_onintr) {
	 printf("interrupting at word %d\n",si);
	 break; 
      }
      if (TRACE2)
	printf("\nstate %d\n",si);
      if (si%(1+((size_of_dot/sample)*sample))==0) {
	  count_dots++;
	  dotscrosses++;
          PPrintf(".");
	  if (count_dots%50==0)
		PPrintf("%d",count_dots);
      }
      if (duplicates[si]) {
	  ignored_duplicates++;
	  continue;  // already done this word
      }	
      int postsi = poststates[si];
      gen postti=postletters[si];
      boolean duplicate=(prestates[postsi] == si) &&
                        (preletters[postsi]==postti);
      if (duplicate) {
        duplicates[postsi]=TRUE;
	if (use_alt_prestate) {
 	       int alternative=prestates2[postsi];	
	       if (alternative && (alternative != prestates[postsi])) {
			prestates[postsi]=alternative; 
			preletters[postsi]=preletters2[postsi];
			duplicates[postsi]=FALSE;
		}
	}
      }
    
      if (si>=nr_states)
          break;
      ti=postletters[si];
      if (ti==0x0) continue;
      if (TRACE2)
	printf("calling getprefix\n");
      i=getprefix(si,workbuf2,
		prestates,preletters,maxwordlen);
      if (TRACE2) {
	printf("returning from  getprefix\n");
	printf("prefix length %d\n",i);
      }	
      if (i <= 0)
          continue;
      if (i>maxwordlen) {
          continue;
      }
      i=genstrlen(workbuf2);
      if (i> maxwordlen)
          continue;
      genstrcpy(lhs_word,workbuf2);
  
      // get postfix
      while (!accepting_states[si2])
      {
         if (ti>ngens) {
          break;
         }
        if (i>maxwordlen)
             break;
        lhs_word[i++]=ti;
	si2=poststates[si2];
        ti=postletters[si2];
     }
      if (i==0)
        continue;
      if (i>maxwordlen)
        continue;
      if (i<minwordlen) {
      	if (i>maximum_length) {
	  maximum_length = i; 
          *pos_max=si;
	}
        continue;
      }
      lhs_word[i]=null_char;
      if (verify)
      {
	if (verify_word(lhs_word,gpwa,si,0))
		continue;
	/*int fail_state=verify_word(lhs_word,gpwa,si,verify_qualifier);
	if (fail_state)
	{
	 	printf("verify failed for state %d\n",fail_state);	
		exit(1);
	}
	continue;*/
      }
      if (TRACE2) {
            int ii=0;
            while (ii<i)
            {
                printf("%02x ",lhs_word[ii++]);
            }
            printf("\n");
      }
      //look for target state corresponding to last letter of lhs
      j=1;
      while (j<=diff_size) {
          gen * wd = diff2->states->words[j];
          if (genstrlen(wd)==1) 
              if (wd[0]==lhs_word[i-1])
                    break;
          j++;
      }
      if (j>diff_size)
          continue;
      if (prefixby1)
      {
         // insert a generator at the start of the  lhs, if lhs[0::len-1] is   
         // still recognised by wa, then use this lhs to look for word_diffs 
	 // (cheaper alternative to creating a minred1 fsa?)
	gen g;
	int sw=0;
	int len=genstrlen(lhs_word);
	for (g=1;g<=ne;g++)
	{
		sw=watable[g][1];
		int i=0;
		while (i<len-1)
		{
			sw=watable[lhs_word[i]][sw];
			if (sw==0)
				break;
			i++;
		}
		if (sw!=0)
			break;
	}
	if (sw!=0)
	{	
		int i=len;
       	 	lhs_word[i+1]=null_char;
		while (i>0)
		{
			i--;
			lhs_word[i+1]=lhs_word[i];
		}
		lhs_word[0]=g;
	}
	else
		continue;
      }
 //         genstrcpy(rhs_word,lhs_word);
  //        diff_reducex(rhs_word,rs_wd2);
//	  if (rhs_word[0] != lhs_word[0]) {
  // awful!!i
 //               display_eqn(si,lhs_word,rhs_word,
//				diff2->alphabet->base->names);
//		  break;
//	  }
	
      if (TRACE2)
		printf("one_level_reducible called. target is %d\n",j);
      processed_words++;
      total_length+=genstrlen(lhs_word);
      if (!one_level_reducible(lhs_word,table_dptr,wa_for_1_level,j,block1,block2))
      {
    	  int accept_state=0; 
    	  if (wa_for_1_level->num_accepting == 1) // using andnot file as the word acceptor
  		accept_state = wa_for_1_level->accepting[1];
	  if (TRACE2)
		printf("calling genstrcpy\n");
          genstrcpy(rhs_word,lhs_word);
	  if (TRACE2) {
                 printf("calling diff_reduce\n");
          }

    	  int **watable=wa_for_1_level->table->table_data_ptr;
	  int len=genstrlen(rhs_word);
	  int i=0; 
	  boolean reduced=TRUE;
	  int state=1;
	  // check if rhs_word is accepted by wa
	  while (i<len) {
		state=watable[rhs_word[i]][state];
		i++;
		if (state==0 || state==accept_state){
			reduced=FALSE;
			break;
		}
	  }
	  if (!reduced)
          //if (!genstrcmp(rhs_word,lhs_word)) 
          // cpcdaft
	   // if (1) 
	  {
	  //  diff_reduce doesnt always reduce a word  which is minimally reducible and reduces to a shorter word
	  //  but fsa_wa_x can build a wa, where such words are marked as reducible
	  // form inverse of rhs_word, reduce it, then reduce the inverse of the reduction.
		//PPrintf("!");
  		gen rhs_word2[MAX_WORD_LEN];
          	genstrcpy(rhs_word2,rhs_word);
	  	//int len = genstrlen(rhs_word);
		//int i=0;
		i=0;
	  	len--;
          	while (len >= 0) 
		{
			rhs_word[len]=inv[rhs_word2[i]];
			i++;
			len--;
		}
		//len=i-1;
		i=0;
		diff_reducex(rhs_word,rs_wd2);
	  	len = genstrlen(rhs_word);
		len--;
          	while (len >= i) 
		{
			if (len>i) {
				gen save_gen=rhs_word[i];
				rhs_word[i]=inv[rhs_word[len]];
				rhs_word[len]=inv[save_gen];
			}
			else
				rhs_word[len]=inv[rhs_word[i]];
			i++;
			len--;
		}
		diff_reducex(rhs_word,rs_wd2);
          	if (!genstrcmp(rhs_word,rhs_word2)) 
		{
			// need a better diff_reduce!
			//PPrintf("*!");
		  if (first_time_display) {
                  	display_eqn(si,lhs_word,rhs_word,
				diff2->alphabet->base->names);
			first_time_display=FALSE;
                  }
			//break;
			continue;
			
		}
		
	  }
	  if (TRACE2)
	  {
		printf("returned from diff_reduce\n");
		printf("look_for_diffs in lhs_word\n");
	  }
	 int ignore=0;
	 if (genstrlen(lhs_word)>(genstrlen(rhs_word)+2)) {
				//PPrintf("?");
		//bug in gpwa (kbprog, maf or fsa_wa_x?)
		ignore=1;
	}
        else if (look_for_diffs(lhs_word,rhs_word,rs_wd,new_states,
                    new_states_number,si,table_dptr,trace_equations,inv,&dotscrosses)) {
              eqn.lhs=lhs_word;
              eqn.rhs=rhs_word;
	      if (TRACE2)
		printf("add_wd_fsa called\n");
              add_wd_fsa(new_diff2,&eqn,inv,TRUE,rs_wd);                   
              if (trace_equations) {
                  display_eqn(si,lhs_word,rhs_word,
				diff2->alphabet->base->names);
              }
          }
	  else {
		// this shouldn't happen normally  - the diff2 set is inconsistent with the wa -
		// to avoid being stuck in 1000's of unnecessary diff_reduce's timeout after  so many secs
                // or bug or limit hit in one_level_reduce
		fails++;
		if (timeout==0)
		{
			if (!andnotx) { // we seem to get this error when the andnot is truncated
				PPrintf("!t");
				timeout = total_t + 300;
                  		display_eqn(si,lhs_word,rhs_word,
					diff2->alphabet->base->names);
			}
	 /*
	  {int state=1;
	  PPrintf("%d>",state);
	  len = genstrlen(lhs_word);
          i=0;
	  while (i<len) {
		state=watable[lhs_word[i]][state];
	  	PPrintf("%d>",state);
		i++;
		if (state==0 || state==accept_state){
			reduced=FALSE;
			break;
		}
	  }}
         */
		}
	  }
	 } // if not one_level_reducible
   } 
   if (si==nr_states)
	si=0;
   end_t=clock();
   total_t=(end_t-start_t)/CLOCKS_PER_SEC;
   PPrintf("\nscan duration=%lds, efficiency fails=%d\n",total_t,fails);

   //printf("\nprocessed %d words, ignored %d duplicate words\n",
   //		processed_words,ignored_duplicates);
   PPrintf("\n ignored %d duplicate words\n", ignored_duplicates);
   if (!verify&&(processed_words>0))
    PPrintf("average word length = %d\n",total_length/processed_words);

// delete the areas creating in copy_table_dptr
   ngens1=diff2->alphabet->base->size+1; 
   free(table_dptr[0][0]); // data area;
   int ll;
   for (ll=0;ll<=ngens1;ll++)
	free(table_dptr[ll]);
   free(table_dptr);
   free(block1);
   free(block2);
   free(new_states);
   free(prestates);
   free(poststates);
   free(preletters);
   free(postletters);
   free(accepting_states);
   free(duplicates);
   if (use_alt_prestate) {
  	free(prestates2);
  	free(preletters2);
   }
   return si;
}  
int verify_word(gen * lhs, fsa *gpwa, int si, int start_reduced)
{
	//printf("\b\b\b\b\b\b\b\b\b%d",si);
	int ** watable=	gpwa->table->table_data_ptr;
	int state=1;
	int i;
	int start_unreduced=0;
	int n=genstrlen(lhs);
//	printf("start reduced=%d\n",start_reduced);
	if (start_reduced != 9999) {
		if (start_reduced > n-2)
			return 0;
	}
	if (start_reduced==0)
		start_unreduced=1;
	else if (start_reduced > 1)
		start_unreduced=start_reduced - 1;
	// verify that slice [0:n] reduces only on the last letter
	int prev_state=1;
	for (i=0;i<n;i++)
	{
		prev_state=state;
		state=watable[lhs[i]][state];
		if (prev_state==state)
			return 1; // dont trust - maybe be truncated wa
		if (state==0) break;
	}
	return 0;
	if ((state!=0)||(i!=n-1))
	{
		printf("bad 1 on state %d\n",si);
		return si;
	}
	if (start_reduced == 9999)
		return 0;
	if (start_unreduced > 0)  
	{
	// verify word reduces on slice (start_unreduced:n) 
	state=1;
	for (i=start_unreduced;i<n;i++)
	{
		state=watable[lhs[i]][state];
		if (state==0) break;
	}
	if (state != 0)	
	{
		printf("bad 2 on state %d, red=%d,unr=%d\n",si,
			start_reduced,start_unreduced);
		return si;
	}
	}
	if (start_reduced==0)
		return 0;
	// verify word is reduced on slice (start_reduced:n) 
	state=1;
	for (i=start_reduced;i<n;i++)
	{
		state=watable[lhs[i]][state];
		if (state==0) break;
	}
	if (state == 0)	
	{
		printf("bad 3 on state %d\n",si);
		return si;
	}
	return 0;
}

int getprefix(int state, gen * workbuf2,int *prestates, 
			gen *preletters,int maxwordlen)
{
      int aim_state=state;
      int bufi=0;
      gen workbuf1 [MAX_WORD_LEN];
      gen null_char = 0x0;
      while (aim_state > 1) {
          int cstate=prestates[aim_state];
          // cstate is smallest state pointing at aim_state
          int letter=1;
          if (bufi==MAX_WORD_LEN) {
          	workbuf2[0]=null_char;
                return 0;
          }
          if (bufi>maxwordlen) {
                return bufi;
          }
          workbuf1[bufi++]= preletters[aim_state];
          aim_state=cstate;
      }
      int i=bufi;
      // reverse the word
      while (i>0) {
        int ii=i-1;
        workbuf2[bufi-i]=workbuf1[ii];
        i--;
      }
      workbuf2[bufi]=0x0;
      return bufi;
}

int ***copy_table_dptr(diff2)
    fsa *diff2;
{
      int ***table_dptr, *data_area;
      int diff_size=diff2->states->size;
      int ngens= diff2->alphabet->base->size;
      int ngens1=ngens+1;
      int ll, rr, s;
      int size_data_area =
          (ngens1+1)*(ngens1+1)*(diff_size+(diff_size*DIFF2_INCREASE_FACTOR)); 
      tmalloc(data_area,int,size_data_area);
      memset(data_area,0,sizeof(int)*size_data_area);	
      tmalloc(table_dptr,int **,ngens1+1);
      for (ll=0;ll<=ngens1;ll++) {
          tmalloc(table_dptr[ll],int *,ngens1+1);
          for (rr=0;rr<=ngens1;rr++) {
              table_dptr[ll][rr] = data_area +
                                  (ll*(ngens1+1)*(diff_size+(diff_size*DIFF2_INCREASE_FACTOR)))
                                             + rr*(diff_size+(diff_size*DIFF2_INCREASE_FACTOR));
              for (s=1;s<=diff_size;s++)  {
                  if (ll>0 && rr>0)
                    if (!((ll==ngens1) && (rr==ngens1)))
                      table_dptr[ll][rr][s]=diff2->table->table_data_dptr[ll][rr][s];
              }
          }
       }

       table_dptr[0][0][0]=diff_size;
       table_dptr[0][0][1]=ngens;
       table_dptr[0][0][2]=diff2->initial[1];
       return table_dptr;
}

void display_eqn(state,cwordx,rword,gpletters)
 int state;
 gen * cwordx;
 gen * rword;
 char ** gpletters;

{
          int   j=0;
          int clength=genstrlen(cwordx);
          printf(" %d. ",state);
          while (j<clength-1) {
              printf("%c",gpletters[(int) cwordx[j++]][0]);
          }
          printf("%c",gpletters[(int) cwordx[j]][0]);
	  if (rword) {
          printf("\n -> ");
          j=0;
          clength=genstrlen(rword); 
          while (j<clength-1) {
              printf("%c",gpletters[(int) rword[j++]][0]);
          }
          printf("%c",gpletters[(int) rword[j]][0]);
	  }
          printf("\n\n");
} 

boolean look_for_diffs(lhs,rhs,rs_wd,new_states,new_states_number,
                       state,table_dptr,trace_equations,inv,dotscrosses)
 gen *lhs;
 gen *rhs;
 reduction_struct *rs_wd;
 gen new_states [][MAXWDLENGTH];
 int *new_states_number;
 int state;
 int ***table_dptr;
 boolean trace_equations;
 int *inv;
 int *dotscrosses;
{
  // procedure should be called when we know that are new wdiffs to find.
  // Each new wdiff causes table_dptr to be updated.
  boolean newone=FALSE;
  fsa *diff2=rs_wd->wd_fsa;
  int diff_size=diff2->states->size;
  int xlen=genstrlen(lhs);

  char **gpletters=diff2->alphabet->base->names;

  int xi=0;
  int s=1;
  int old_s;
  xi=0;
  char wdx [2*xlen+1];
  gen wdr [2*xlen+1];
  gen wd [2*xlen+1];
  int si=state;
  gen null_char = 0x0;
  if (TRACE2)
      printf("Entering look_for_diffs\n");
  while (xi<xlen) {
    int wdi=xi;
    int xj=0;
    wd[2*(xi+1)]=0;
    while (wdi >= 0) {
        gen xx=lhs[wdi],xx_inv;
        //xx_inv = (xx%2 == 0) ? xx-1 : xx+1;
        xx_inv = inv[xx];
        wd[xi-wdi]=xx_inv;
        wd[2*(xi+1)-(wdi+1)]=rhs[xi-wdi];
        wdi--;
    } 
   // if (TRACE2)
    //    printf("calling diff_reduce on wd..");
    diff_reducex(wd,rs_wd);
    //if (TRACE2)
     //   printf("done\n");

    if (trace_equations)
    {
 	   wdi=genstrlen(wd);
 	   while (xj<wdi) {
	        if (gpletters[wd[xj]] == 0)
 	           return FALSE;
 	       wdx[xj]=gpletters[wd[xj]][0];
 	       xj++;
	    }
	    wdx[wdi]=0x0;
    }

    xj=1;
    while (xj<=diff_size) {
        //search the existing diffs
        int cmp=genstrcmp(diff2->states->words[xj],wd);
        if (cmp==0) {
		if (TRACE2)
		   printf ("wd%d ",xj);
                break;
        }
        xj++;
    }
    if (xj>diff_size) {
        xj=0;
        while (xj<=(diff_size*DIFF2_INCREASE_FACTOR)) {
            // search the new diffs that have been found already
          if (new_states[xj][0] == null_char)
            break;
          if (!genstrcmp(new_states[xj++],wd)) {
            xj=-1;
            break;
          }
        }
    }
    else
        xj=-1;
    if (xj >= 0) // new state
    {
        int nsi=xj;
        if ((genstrlen(wd)<MAXWDLENGTH) && (xj < (diff_size*DIFF2_INCREASE_FACTOR)))
        {
          int i=0;
          PPrintf("x");
	  (*dotscrosses)++;
          if (trace_equations) {
              printf( " new word difference %s (%d) ",wdx,xi+1);
	      int ii=0;
	     printf("\n,");
	     while (wdx[ii]!=0x0) 
		printf("%c\\*",wdx[ii++]);
	     printf("\n");
          }
          // calculate inverse
          xj=genstrlen(wd);
          while (xj > 0) {
            xj--;
            //wdr[i++]=wd[xj]%2?wd[xj]+1:wd[xj]-1;
            wdr[i++]=inv[wd[xj]];
          }
          wdr[i]=null_char;
          if (TRACE2)
            printf("\ncalling diff_reduce on wdr..");
          diff_reducex(wdr,rs_wd);
         if (TRACE2)
            printf("done\n");
	 if (trace_equations) {
          while (xj<i) {
            wdx[xj]=gpletters[wdr[xj]][0];
            xj++;
          }
	 }
         // update new_states and diff table
          new_states_number[nsi]=table_dptr[0][0][0]+1;
          genstrcpy(new_states[nsi++],wd);
          new_states[nsi][0]=null_char;
          table_dptr[0][0][0]++;
          if (TRACE2)
              printf("\ncalling update_table");
          update_table(table_dptr,rs_wd,wd,
			new_states,new_states_number,inv);
          if (genstrcmp(wd,wdr))
          {
          new_states_number[nsi]=table_dptr[0][0][0]+1;
          genstrcpy(new_states[nsi++],wdr);
          new_states[nsi][0]=null_char;
          table_dptr[0][0][0]++;
          if (TRACE2)
              printf("\ncalling update_table");
          update_table(table_dptr,rs_wd,wdr,
			new_states,new_states_number,inv);
  	  } 
          if (trace_equations) 
              printf ("(%s)\n",wdx);
          newone=TRUE;
        } // if strlen < 
    } // if xj > 0
    xi++;
  } // xi < xlen
  if (TRACE2)
      printf("\nreturning %d from look_for_diffs\n",newone);
  return newone;
}

void update_table(table_dptr,rs_wd,wd,new_states,new_states_number,inv)
 int ***table_dptr;
 reduction_struct *rs_wd;
 gen *wd;
 gen new_states [][MAXWDLENGTH];
 int *new_states_number;
 int *inv;
{
    fsa *diff2=rs_wd->wd_fsa;
    int diff_size=diff2->states->size;
    int ngens = diff2->alphabet->base->size;
    int ngens1 = ngens+1;
    int xi;
    int xj;
    int xn;
    gen buf[MAXWDLENGTH*sizeof(gen)];
    int newstate=table_dptr[0][0][0];
    int revi,revj;
    int ns=0; 
    gen null_char = 0x0;
    for (xi=1;xi<=ngens;xi++) {
        for (xj=1;xj<=ngens;xj++) {
            int len=1+genstrlen(wd)+1;
            buf[0]=xi;
            genstrcpy(buf+1,wd);
            buf[len-1]=xj;
            buf[len]=null_char;
            diff_reducex(buf,rs_wd);
            xn=1;
            while (xn<=diff_size) {
                //search the existing diffs
                int cmp=genstrcmp(diff2->states->words[xn],buf);
                if (cmp==0) 
                        break;
                xn++;
            }
            //revi = (xi%2 == 0) ? xi-1 : xi+1;
	    //revj = (xj%2 == 0) ? xj-1 : xj+1;
	    revi = inv[xi];
            revj = inv[xj];
 
            if (xn<=diff_size) {
 
                table_dptr[revi][xj][newstate]=xn;
                //and reverse
                table_dptr[xi][revj][xn]=newstate;
            }
            else {
                // search the new diffs
                ns=0;
                xn=0;
                while (new_states[xn][0])
                {
                    int cmp=genstrcmp(new_states[xn],buf);
                    if (cmp==0) {
                        ns=new_states_number[xn];
                        break;
                    }
                    xn++;
                }
                if (ns) {
                    table_dptr[revi][xj][newstate]=ns;
                    //and reverse
                    table_dptr[xi][revj][ns]=newstate;
                }
                else
                    table_dptr[revi][xj][newstate]=0;
            }
        }
    }
}

boolean one_level_reducible (lhs,dtable,gpwa,target_state,block1,block2)
        gen  *lhs;
        int ***dtable;
        fsa *gpwa;
        int target_state;
	unsigned int *block1;
	unsigned int *block2;
{
    // attempt at a fast way to determine if new wds
    // can be extracted from the lhs.

    // we create 'state's which contain pairs of fsa states (sd,sw), where sd is a state
    // from diff2 and sw is state from wa.
    // Let n be length of lhs.
    // Up to n+1 states are created. State 1 consists of the pair (1,1).
    // Each state x (2 to n+1) is derived from state x-1.
    // For each (sd,sw) in x-1 apply lhs[x-2],g2 to sd and g2 to sw to make
    // (sd',sw'). If both sw' and sd' are not=0 then append (sd',sw') to state x.
    //  g2 ranges over all generators. When x is n, g2 may also be $.
    // If x = n+1, then for all (sd,sw) in state x apply $,g2 to sd and g2 to sw.
    // If sd'=target and sw'>0 and x is n or n+1 then return TRUE,
    //  else return FALSE after all n+1 states have been created.
    
    int n, x, xmin2, ngens,dollar, **watable;
    int g1,g2;
    unsigned int  *xminus1_state, *x_state;
    int no_xpairs, no_xminus1pairs;
    boolean block1to2;
    int diff_identity;
    int accept_state=0; 
    if (gpwa->num_accepting == 1) // using andnot file as the word acceptor
  	accept_state = gpwa->accepting[1];
    block1to2=TRUE;
    x_state=block1;
    no_xpairs=0;
    diff_identity=dtable[0][0][2]; // stored here when we copied table
    x_state[0]=diff_identity;
    x_state[1]=gpwa->initial[1];  // initial state of gpwa
    // x_state[0]=gpwa->initial[1];  // initial state of gpwa
    no_xpairs=1;
    watable=gpwa->table->table_data_ptr;

    ngens=dtable[0][0][1]; // stored here when we copied table
    dollar=ngens+1;
    x=2;
    xmin2=0;
    n=genstrlen(lhs);
    if (TRACE3 || TRACE1){
        printf("One_level_reducible. Target is %d\n",target_state);
	int ii=0;
	while (ii<n)
		printf("%02x ",lhs[ii++]);
	printf("\n");
    }
    while (x<=n+1)
    {
        int sd,sw, remaining;
        if (block1to2) {
            xminus1_state=block1;
            x_state=block2;
            block1to2=FALSE;
        }
        else
        {
            xminus1_state=block2;
            x_state=block1;
            block1to2=TRUE;
        }
        no_xminus1pairs=no_xpairs;
        no_xpairs=0;
        if (TRACE3||TRACE1)
            printf("**state %d **\n",x);
        if (no_xminus1pairs==0)
        {
          return FALSE;
        }
        while (no_xminus1pairs>0) {
            int disp;
            int g2_range=ngens;
            no_xminus1pairs--;
            if (x==n)
              g2_range=dollar;
            disp=no_xminus1pairs*2;
	    if (TRACE3)
		printf("read sd,sw from xminus1 at %d\n",disp);
            sd=xminus1_state[disp];
            sw=xminus1_state[disp+1];
	    if (TRACE3)
		printf("values are %d,%d \n",sd,sw);
            if (x==n+1)
	    {
                g1=dollar;
		remaining=0;
	    }
            else
	    { 
                g1= lhs[xmin2]; 
		remaining=xmin2+1;
		if (TRACE3)
			printf("g1 is %d from %02x\n",g1,lhs[xmin2]);
	    }
            g2=1;
            while (g2<=g2_range) {
                int sd_dash=0, sw_dash=0;
		if (TRACE3)
			printf("read from dtable with %d, %d, %d\n",g1,g2,sd);
                sd_dash=dtable[g1][g2][sd];
                if (sd_dash != 0) {
                    if (g2==dollar)
                        sw_dash=sw;
                    else
		    {
			if (TRACE3)
				printf("read from watable with %d,%d\n",g2,sw);
                        sw_dash=watable[g2][sw];
	            }	
                    if (sw_dash !=0 && sw_dash != accept_state)
                    {
                        if ((x>=n) && (sd_dash==target_state)) {
                                if (TRACE1) {
                                 printf("TRUE %d->%d via %d,%d\n",sd,sd_dash,g1,g2);
                                 exit(99);
                                 }
                                return TRUE;
                        }
                        else
                        {
                            if (x<n+1) {
     				  //check if remainder of word is inverse to value of sd_dash
				  int next=remaining;
				  int next_letter;
				  int next_sd; 
			          if (next>0 && ((n-next)<=5)) {
					next_letter=lhs[next];
					next_sd=sd_dash;
					while (next_letter!=0) {
						next_sd=dtable[next_letter][dollar][next_sd];
						if (next_sd==0)
							break;
						if (next_sd==1) 
							return TRUE;
						next_letter=lhs[++next];
					}
				  }
				
                                // add (sd_dash,sw_dash) to state x
                                if (no_xpairs<MAXPAIRS)
                                {
                                    disp=no_xpairs*2;
                                    no_xpairs++;
					if (TRACE3)
					printf("write values to x_state at %d and %d\n",
						disp,disp+1);
                                    x_state[disp]=sd_dash;
                                    x_state[disp+1]=sw_dash;
                                }
                                if (TRACE1 || TRACE3) {
                                    printf("%d->%d via %d,%d\n",sd,sd_dash,g1,g2);
				    /*if (TRACE3)
				    {
				    	printf("havent crashed this time - try running with TRACE3 0\n");
			            	exit(99);
				    }*/
				}
                           } 
                        }
                    }
                }
                g2++;
            } // g2<=g2_range
        } // while no_xminus1pairs>=0
        x++;
	xmin2++;
    } // while x<=n+1
    if (TRACE1)
    {	
	printf("FALSE\n");		
        exit(99);
    }
    return FALSE;
}


fsa * fsa_triples_big_hash(minredptr,waptr,diffptr,op_table_type,
           destroy,tempfilename,hash_table_size,counter,minred_type,state_limit)
	fsa *minredptr, *waptr, *diffptr;
	storage_type op_table_type;
	boolean destroy;
	char *tempfilename;
        unsigned int hash_table_size;
	boolean counter;
	boolean minred_type;
        int state_limit;
	
{
  /* this is a version of fsa_minkb which uses one large hash table whose
     size is specified by the parameter hash_table_size. This value is also
     equal to the maximum number of states that the fsa minkb can contain.
     Please consult further the documentation contained in fsaminkb.c
  */
    
  int **minredtable, **watable, ***difftable, identity, ngens, ngens1, nswa1,
      ne, ns, *fsarow, nt, cstate, cswa1, cswa2, csdiff, im, i, e, len, ct,
      hashval, collisions;
  unsigned int *ht_ptr;
  boolean  dense_op;
  fsa *minkbptr;
  hash_table ht;
  FILE *tempfile, *fopen();
  gen g1, g2;
  reduction_struct rs_wd;
  int ht_i=0; 
  mstate *statetokey [MAX_EXTRA_MEMORY];
  kts *keytostate;
  unsigned int next_state;

  unsigned int no_states = hash_table_size+1;

  if (kbm_print_level>=3)
    printf("    #Calling fsa_triples_big_hash\n");

  if (!minredptr->flags[DFA] || !waptr->flags[DFA] || !diffptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_minkbptr only applies to DFA's.\n");
    return 0;
  }
  if (minredptr->alphabet->type!=IDENTIFIERS ||
                                       waptr->alphabet->type!=IDENTIFIERS) {
    fprintf(stderr,"Error in fsa_minkbptr: an fsa has wrong type.\n");
    return 0;
  }
  if (waptr->num_accepting != waptr->states->size) {
    fprintf(stderr,
       "Error in fsa_minkbptr: second fsa should be a word-acceptor.\n");
    return 0;
  }
  if (diffptr->alphabet->type!=PRODUCT || diffptr->alphabet->arity!=2) {
    fprintf(stderr,"Error in fsa_minkbptr: third fsa must be 2-variable.\n");
    return 0;
  }
  if (diffptr->states->type!=WORDS) {
    fprintf(stderr,
       "Error in fsa_minkbptr: third fsa must be word-difference type.\n");
    return 0;
  }
  if (!srec_equal(diffptr->alphabet->base,waptr->alphabet)) {
    fprintf(stderr,"Error in fsa_minkbptr: fsa's alphabet's don't match.\n");
    return 0;
  }
  if (minredptr->states->size>=MAXINT || waptr->states->size>=MAXINT
                                       || diffptr->states->size>=MAXINT) {
    fprintf(stderr,
       "Error in fsa_minkbptr: one of the fsa's has too many states.\n");
    return 0;
  }
  int diff_partition=hash_table_size/diffptr->states->size;

  if (fsa_table_dptr_init(diffptr)==-1) return 0;
  collisions=0;

  tmalloc(minkbptr,fsa,1);
  fsa_init(minkbptr);
  srec_copy(minkbptr->alphabet,diffptr->alphabet);
  minkbptr->flags[DFA] = TRUE;
  minkbptr->flags[ACCESSIBLE] = TRUE;
  minkbptr->flags[BFS] = TRUE;

  ngens = waptr->alphabet->size;
  ngens1 = ngens+1;
  ne = diffptr->alphabet->size;
  nswa1 = waptr->states->size + 1;

  if (ne != ngens1*ngens1-1) {
   fprintf(stderr,
       "Error: in a 2-variable fsa, alphabet size should = ngens^2 - 1.\n");
    return 0;
  }
// states to receive target of dollar transition
  int end_of_lhs = minredptr->states->size + 1;
  int end_of_rhs = waptr->states->size + 1;
  int dollar = ngens1; 
  // allocate hash and data areas
  int state_i=0;
  while (state_i<MAX_EXTRA_MEMORY) {
	statetokey[state_i]=NULL;
	state_i++;
  }	

  tmalloc(statetokey[0],mstate,no_states);
  tmalloc(keytostate,kts,no_states);
  for (i=0;i<no_states;i++) {
	keytostate[i].state=0;
	keytostate[i].next_value=0;
   }

  minredtable = minredptr->table->table_data_ptr;
  watable = waptr->table->table_data_ptr;
  difftable = diffptr->table->table_data_dptr;

  dense_op = op_table_type==DENSE;

  minkbptr->states->type = SIMPLE;

  minkbptr->num_initial = 1;
  tmalloc(minkbptr->initial,int,2);
  minkbptr->initial[1] = 1;
  minkbptr->table->table_type = op_table_type;
  minkbptr->table->denserows = 0;
  minkbptr->table->printing_format = op_table_type;
  //  we use the kbmag hash area as a work area
  hash_init(&ht,TRUE,3,0,0);
  ht_ptr = ht.current_ptr;
  ht_ptr[0] = minredptr->initial[1];
  ht_ptr[1] = waptr->initial[1];
  ht_ptr[2] = identity = diffptr->initial[1];
  im = hash_locate(&ht,3);
  if (im!=1) {
    fprintf(stderr,"Hash-initialisation problem in fsa_minkbptr.\n");
    return 0;
  }
  statetokey[ht_i][1].wa1=ht_ptr[0];
  statetokey[ht_i][1].diff=diffptr->initial[1];
  statetokey[ht_i][1].wa2=waptr->initial[1];
  
  if ((tempfile=fopen(tempfilename,"w"))==0){
    fprintf(stderr,"Error: cannot open file %s\n",tempfilename);
    return 0;
  }
  if (dense_op)
    tmalloc(fsarow,int,ne)
  else
    tmalloc(fsarow,int,2*ne+1)
 
  cstate = 0;
  if (dense_op)
    len = ne; /* The length of the fsarow output. */
  nt = 0; /* Number of transitions in minkbptr */
  next_state=ht.num_recs;
  if (counter)
      printf("processed/total states\n");
  //boolean truncate_msg_first_time=TRUE;
  int accept_state=minredptr->accepting[1];
  int local_state=0;
  int tmin1=no_states-1;
  while (++cstate <= next_state) {
        if (cstate%500==0 && counter)
            printf(
            "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%d  %d",
                  cstate,next_state);
    local_state=cstate;
    state_i=0;
    if (cstate>tmin1) {
	state_i=cstate/tmin1;
	local_state=cstate%tmin1;
	if (local_state==0)
	{
		local_state=tmin1;
		state_i--;
	}
    }
    cswa1=statetokey[state_i][local_state].wa1; //minred or pfxred
    csdiff=statetokey[state_i][local_state].diff;
    cswa2=statetokey[state_i][local_state].wa2; // wa

    if (dense_op)
      for (i=0;i<ne;i++) fsarow[i] = 0;
    else
      len = 0;
    e = 0; /* e is the number of the edge corresponding to the pair (g1,g2) */
    for (g1=1;g1<ngens1;g1++) for (g2=1;g2<=ngens1;g2++) {
      e++;
// Calculate action of generator-pair (g1,g2) on state cstate - 
//  we insert the code from MAF routine product_intersection
//  fsa.cpp lines 5087-5119 to make the unminimised minkb as small as possible.
      if (csdiff==1)
    	 if (minred_type && g1==g2) 
	   // can't occur if lhs is minimally reducible
           // or restricted prefix_reduced
	  continue;
      ht_ptr = ht.current_ptr;
      ht_ptr[2] = dense_dtarget(difftable,g1,g2,csdiff);
      if (ht_ptr[2]==0)
        continue;
      ht_ptr[0] = dense_target(minredtable,g1,cswa1);
      if (state_limit)
       if (ht_ptr[0]>state_limit)
         ht_ptr[0]=0;
      if (ht_ptr[0]==0)
          continue;
      if (cswa2  == end_of_rhs) 
	// x,$ followed by y,$ = xy,$, x,$ followed by y,z can't happen
	ht_ptr[1] = g2 == dollar ? end_of_rhs : 0;
      else if (g2==dollar)
	// x,$
	ht_ptr[1] = end_of_rhs;   
      else	
	ht_ptr[1] = dense_target(watable,g2,cswa2);
      if (ht_ptr[1]==0)
          continue;
      //t1++;
      if (ht_ptr[0] == accept_state && ht_ptr[2]!=1)
	// last letter of lhs, must join up with last letter of rhs 
	   continue;
	  //{c1++;continue;}
      if (cswa2 == end_of_rhs && (ht_ptr[2]!=1 || ht_ptr[0]!=accept_state))
	// g1,$ only relevant if g1 is last letter  of lhs
	   continue;
	  //{c2++;continue;}
	
      im=get_image_big_hash (ht_ptr,next_state,
                        keytostate,&collisions,no_states,diff_partition);
	
      // im is > 0 at this point 
      if (im>next_state) {
	   local_state=im;
	    state_i=0;
	   if (im > tmin1) {
		   state_i=im/tmin1;
		   local_state=im%tmin1;
		   if (local_state==0) {
			state_i--;
			local_state=tmin1;
		   }
		   if (state_i > 99)
		   {
			printf("out of memory!\n");
			exit(99);
		   }
		   if (statetokey[state_i]==NULL) {
			Printf("exceeded limit of %d states\n",
				state_i*tmin1);
			Printf("allocating  memory for %d more states\n",
				tmin1);
			tmalloc(statetokey[state_i],mstate,no_states);
		   }
	   }
           next_state=im;
           statetokey[state_i][local_state].wa1=ht_ptr[0];
           statetokey[state_i][local_state].wa2=ht_ptr[1];
           statetokey[state_i][local_state].diff=ht_ptr[2];
      }
      fsarow[++len] = e;
      fsarow[++len] = im;
      nt++;
    }   /* for (g1=1;g1<=ngens1; ... */
    fsarow[0] = len++;
    fwrite((void *)fsarow,sizeof(int),(size_t)len,tempfile);
  }  /*while (++cstate <= ht.num_recs) */
  //printf("nt %d t1 %d c1 %d c2 %d \n",nt,t1,c1,c2);
  fclose(tempfile);

  ns=next_state;
  minkbptr->states->size=ns;
  minkbptr->table->numTransitions = nt;

/* Now locate the accept states */

  fsa_set_is_accepting(minredptr);
  ct = 0;
  for (i=1;i<=ns;i++) {
    // the states are spread across statetokey tables each holding 
    // no_states -1 states  starting from index 1
    state_i=0;
    int local_i=i;
    if (i>tmin1)
    {
	    state_i= i /tmin1;
	    local_i= i%tmin1;
	    if (local_i==0) {
		state_i--;
		local_i=tmin1;
	    }
    }
    mstate mm = statetokey[state_i][local_i];
    if (minredptr->is_accepting[mm.wa1] && mm.diff==identity)
              ct++;
  }
  minkbptr->num_accepting = ct;
  if (ct==1 || ct != ns) {
    tmalloc(minkbptr->accepting,int,ct+1);
    ct = 0;
    for (i=1;i<=ns;i++) {
    // the states are spread across statetokey tables each holding 
    // no_states -1 states  starting from index 1
	    state_i=0;
	    int local_i=i;
	    if (i>tmin1)
	    {
		    state_i= i /tmin1;
		    local_i= i%tmin1;
		    if (local_i==0) {
			state_i--;
			local_i=tmin1;
		    }
	    }
	    mstate mm = statetokey[state_i][local_i];
	    if (minredptr->is_accepting[mm.wa1] && mm.diff==identity)
             minkbptr->accepting[++ct] = i;
    }
  }
// free the big_hash areas
   state_i=0;
   while (state_i<MAX_EXTRA_MEMORY) {
	if (statetokey[state_i]!=NULL)
	   tfree(statetokey[state_i]);
	state_i++;
   }
   float no_hash_entries=0;
   float no_next_values=0; 
   struct kts *this_kts, *next_kts, *to_free_kts;;
   int kts_i=0;	
   while (kts_i<no_states)
   {	
	this_kts=keytostate+kts_i;
	if (this_kts->state)
		no_hash_entries++;
	to_free_kts=this_kts->next_value;
	while (to_free_kts !=0)
	{		
		no_next_values++;
		next_kts=to_free_kts->next_value;
		free(to_free_kts);
		to_free_kts=next_kts;
	}
	kts_i++;
   } 
  tfree(keytostate);
  double depth=
	(no_next_values+no_hash_entries)/no_hash_entries;
      if (DISPLAY_HASH_DEPTH)
  printf("average hash depth (1 is perfect) = (%0.0f + %0.0f)/%0.0f = %0.2f\n",
       no_hash_entries,no_next_values,no_hash_entries,depth);

  hash_clear(&ht);
  tfree(minredptr->is_accepting);
  tfree(fsarow);
  if (destroy) {
    fsa_clear(minredptr); fsa_clear(waptr); fsa_clear(diffptr);
  }
/* Now read the transition table back in */
  tempfile = fopen(tempfilename,"r");
  compressed_transitions_read(minkbptr,tempfile);
  fclose(tempfile);

  unlink(tempfilename);

  return minkbptr;
}

int get_image_big_hash(rec,next_state,keytostate,
			collisions,no_states,diff_partition)
    mstate *rec;
    int next_state;
    kts *keytostate;
    unsigned int *collisions;
    unsigned int no_states;
    int diff_partition;
{
    unsigned int hashval=0;
    unsigned int m, hv, i, k;
    struct kts *this_kts, *next_kts;
    hashval=((rec->wa1)*101+(rec->diff)*diff_partition+(rec->wa2))%no_states;
    this_kts=keytostate+hashval;
    while (1) {
        if (this_kts->state==0) {
            //new states
            this_kts->state=next_state+1;
            memcpy(&this_kts->value,rec,sizeof(mstate));
            this_kts->next_value=0;
            break;
        }
        if (memcmp(&this_kts->value,rec,sizeof(mstate))==0) 
            break;
        if (this_kts->next_value==0) {
            (*collisions)++;
            tmalloc(next_kts,kts,1);
            this_kts->next_value=next_kts;
            this_kts=next_kts;
            this_kts->state=0;
            continue;
        }
        this_kts=this_kts->next_value;
    }
    return this_kts->state;
}

void do_andnot (fsafile1,fsafile2,gp)
	char *fsafile1;
	char *fsafile2;
        char *gp;
{
        FILE *rfile;
        fsa *fsa1, *fsa2, *fsaandnot;
        char fsaname [100], tempfilename [100],filename [100];
        storage_type op_store = DENSE;
        if ((rfile = fopen(fsafile1,"r")) == 0) {
                fprintf(stderr,"Cannot open file %s.\n",fsafile1);
                exit(1);
        }
        tmalloc(fsa1,fsa,1);
        fsa_read(rfile,fsa1,DENSE,0,0,TRUE,fsaname);
        fclose(rfile);
        
        if ((rfile = fopen(fsafile2,"r")) == 0) {
                fprintf(stderr,"Cannot open file %s.\n",fsafile2);
                exit(1);
        }
        tmalloc(fsa2,fsa,1);
        fsa_read(rfile,fsa2,DENSE,0,0,TRUE,fsaname);
        fclose(rfile);

	printf("doing fsa_andnot on %s and %s\n",fsafile1,fsafile2);
        if (fsa1->states->size != fsa2->states->size) {
            strcpy(tempfilename,gp);
            strcat(tempfilename,"temp_XXX");
            strcpy(fsaname,"_RWS");
            strcat(fsaname,".andnot");
            strcpy(filename,gp);
            strcat(filename,".andnot");
            if (!fsa1->accepting)
                // convert fail state to accepting state
                fsa1 = fsa_pfxred(fsa1,op_store,tempfilename,0);
            if (!fsa2->accepting)
                // convert fail state to accepting state
                fsa2 = fsa_pfxred(fsa2,op_store,tempfilename,0);

            fsaandnot = fsa_and_not_first(fsa2,fsa1,op_store,TRUE,tempfilename,0);
            //fsaandnot = fsa_and_not(fsa1,fsa2,op_store,TRUE,tempfilename);
            printf("writing %s, states=%d\n",filename,fsaandnot->states->size);
            wfile = fopen(filename,"w");
            fsa_print(wfile,fsaandnot,fsaname);
            fclose(wfile);
	    printf("scan using 'gpcheckx -t %s\n",gp);	
         }
}

// Code for minimizing the fsa
// by using a big hash table and a suitable hash calculator.
// set by -mzbh switch. If this isnt set call the normal fsa_minimize.

int   get_image_x();
void clear_hash ();

typedef struct  kts2{
      int state;
      int *value;
      int len;
      struct kts2 *next_value;
  } kts2;

int fsa_minimize_big_hash(fsaptr,hash_table_size,
				counter,minimize_big_hash)
        fsa     *fsaptr;
        int     hash_table_size;
        boolean counter;
	boolean minimize_big_hash;
/* Minimize the fsa *fsaptr. */
{ int *block_numa, *block_numb, *block_swap, i, j, k, l, len,
       *ptr, *ptrh, *ptr2, *ptr2e, *ht_ptr,
       ne, ns_orig, **table, ns_final, ns_new, num_iterations;
  hash_table ht;
  boolean fixed, acc;
  if (!minimize_big_hash)
	return fsa_minimize(fsaptr);
  Printf("calling fsa_minimize_big_hash\n");

  int next_state,  im, imh;
  kts2 *keytostate;
  int *values;
  int next_value;

  if (fsaptr->table->table_type==SPARSE && fsaptr->table->denserows>0) {
    fprintf(stderr,
"Sorry: fsa_minimize unavailable for sparse storage with dense rows.\n");
    return -1;
  }
 if (kbm_print_level>=3)
    printf("    #Calling fsa_minimize.\n");
  if (!fsaptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_minimize only applies to DFA's.\n");
    return -1;
  }
  if (fsaptr->flags[MINIMIZED])
    return 0;

  if (fsaptr->flags[RWS])
    fsa_clear_rws(fsaptr);

  acc = fsaptr->flags[ACCESSIBLE]  || fsaptr->flags[TRIM];
  if (!acc)
    fsa_set_is_accessible(fsaptr);

  ns_orig = fsaptr->states->size;
  if (ns_orig==0) {
    fsaptr->flags[TRIM] = TRUE;
    fsaptr->flags[MINIMIZED] = TRUE;
    tfree(fsaptr->is_accessible);
    return 0;
  }
 /* First throw away any existing structure on the state-set. */
  srec_clear(fsaptr->states);
  fsaptr->states->type = SIMPLE;
  ne = fsaptr->alphabet->size;
  table = fsaptr->table->table_data_ptr;
// big_hash areas
  int no_states=hash_table_size;
  tmalloc(keytostate,kts2,no_states);
  memset(keytostate,0,sizeof(kts2)*(no_states));

  fixed = fsaptr->table->table_type==DENSE;
  // work out size of values -
  // size of current table + 2 * ns_orig + table[ns_orig] = table[1]
  // + size of table[ns_orig] (2000 integers should cover it!)
  int values_size;
  if (!fixed)
      values_size=2*ns_orig + table[ns_orig] - table[1] + 2000;
  else
      values_size=ns_orig*(1+ne);
  tmalloc(values,int,values_size); 
  tmalloc(block_numa,int,ns_orig+1);
  tmalloc(block_numb,int,ns_orig+1);
  for (i=0;i<=ns_orig;i++)
    block_numb[i]=0;
/* Start with block_numa equal to the accept/reject division
 * Remember that state/block number 0 is always failure with no hope.
 */
  if (fsaptr->num_accepting == ns_orig) {
    block_numa[0]=0;
    for (i=1;i<=ns_orig;i++) if (acc || fsaptr->is_accessible[i])
      block_numa[i] = 1;
    else
      block_numa[i] = 0;
  }
  else {
    for (i=0;i<=ns_orig;i++)
      block_numa[i] = 0;
    for (i=1;i<=fsaptr->num_accepting;i++)
      if (acc || fsaptr->is_accessible[fsaptr->accepting[i]])
        block_numa[fsaptr->accepting[i]] = 1;
  }
 
  ns_new = 1;
  num_iterations = 0;
/* The main refinement loop follows. */
  if (counter)
        printf("iterations states\n");
  do {
    num_iterations++;
    if (counter)
        printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b%d  %d",num_iterations,ns_new);

    ns_final = ns_new;
    next_state=0;
/* Turn off excessive printing at this point */
    j=kbm_print_level; kbm_print_level=1;
    if (TRACE4)
      hash_init(&ht,fixed,ne+1,0,0);
    kbm_print_level=j;
    if (kbm_print_level>=3)
     printf("    #Iterating - number of state categories = %d.\n",ns_new);
    block_numb[0] = 0;
    next_value=0;
    for (i=1;i<=ns_orig;i++) if (acc || fsaptr->is_accessible[i]){
  /* Insert the encoded form of the transitions from state i into the hashtable
   * preceded by the current block number of i.
   */
      len = fixed ? ne+1 : table[i+1]-table[i]+1;
	if (TRACE4)
          ptrh = ht.current_ptr;
       ptr=values+next_value;
      *ptr = block_numa[i];
	if (TRACE4)
         *ptrh = block_numa[i];

      if (fixed) {
        for (j=1;j<len;j++) {
          ptr[j] = block_numa[table[j][i]];
	  if (TRACE4)
           ptrh[j] = block_numa[table[j][i]];
	}
        l = len;
      } // fixed
      else { // not fixed
       l = 0;
       for (j=1;j<len;j+=2) {
          k = block_numa[table[i][j]];
          if (k > 0){
	  	  if (TRACE4)
		  {
       		     ptrh[++l] = table[i][j-1];  
       		     ptr[l] = table[i][j-1];
       		     ptrh[++l] = k;
       		     ptr[l] = k;
       		  }
		  else {
       		     ptr[++l] = table[i][j-1];  
       		     ptr[++l] = k;
		  }
          }
       }
       if (l>0 || *ptr>0)
          l++;
/* For technical reasons, we want the zero record to be empty */
     } // not fixed
     if (TRACE4)
       imh=hash_locate(&ht,l);
     if (l>0)
        im=get_image_x(ptr,l,next_state,keytostate,no_states);
     else
        im=0;
     if (TRACE4) {
      if (im != imh)
      {
         printf("different!!\n");
         //exit(1);
      }
     }
     if (im>next_state) {
        next_state=im;
        next_value=next_value+l;
        if (!fixed) {
          // write and end of record marker
             *(values+next_value)=MAXINT;
              next_value++;  
        }
     }

      //block_numb[i] = hash_locate(&ht,l);
      block_numb[i]=im;
      if (block_numb[i]== -1) return -1;
    } // for i<= ns_orig
    else block_numb[i]=0;
    //ns_new = ht.num_recs;
    block_swap=block_numa; block_numa=block_numb; block_numb=block_swap;
    ns_new=next_state;
    if (ns_new > ns_final) {
      clear_hash (keytostate, no_states);
      if (TRACE4)
       hash_clear(&ht);
    }
	
  } while (ns_new > ns_final);

  if (counter)
	printf("\n");
  if (kbm_print_level>=4) { /* print out old-new state correspondence */
     printf("Old State   NewState\n");
     for (i=1;i<=ns_orig;i++)
       printf("   %6d     %6d\n",i,block_numa[i]);
  }

/* At this stage, either ns_final = ns_new, or the fsa has empty accepted
 * language, ns_new=0 and ns_final=1.
 */

  fsaptr->flags[TRIM] = TRUE;
  fsaptr->flags[MINIMIZED] = TRUE;

  if (ns_new==0) {
/* This is the awkward case of no states - always causes problems! */
    fsaptr->states->size=0;
    fsaptr->num_initial=0;
    tfree(fsaptr->initial);
    fsaptr->num_accepting = 0;
    tfree(fsaptr->accepting);
    tfree(fsaptr->table->table_data_ptr[0]);
    tfree(fsaptr->table->table_data_ptr);
  }
  else if (ns_final<ns_orig) {
/* Re-define the fsa fields  */
    fsaptr->states->size = ns_final;

    fsaptr->initial[1] = block_numa[fsaptr->initial[1]];

    if (fsaptr->num_accepting == ns_orig) {
      fsaptr->num_accepting = ns_final;
      if (ns_final==1) {

        tmalloc(fsaptr->accepting,int,2);
        fsaptr->accepting[1] = 1;
      }
    }
    else {
      tmalloc(fsaptr->is_accepting,boolean,ns_final+1);
      for (i=1;i<=ns_final;i++)
        fsaptr->is_accepting[i] = FALSE;
      for (i=1;i<=fsaptr->num_accepting;i++)
        fsaptr->is_accepting[block_numa[fsaptr->accepting[i]]] = TRUE;
      fsaptr->num_accepting = 0;
      for (i=1;i<=ns_final;i++) if (fsaptr->is_accepting[i])
        fsaptr->num_accepting++;
      tfree(fsaptr->accepting);
      tmalloc(fsaptr->accepting,int,fsaptr->num_accepting+1);
      j = 0;
      for (i=1;i<=ns_final;i++) if (fsaptr->is_accepting[i])
        fsaptr->accepting[++j] = i;
      tfree(fsaptr->is_accepting);
    }

/* Finally copy the transition table data from the hash-table back to the fsa */
    tfree(fsaptr->table->table_data_ptr[0]);
    tfree(fsaptr->table->table_data_ptr);
    if (fixed){
      fsa_table_init(fsaptr->table,ns_final,ne);
      table = fsaptr->table->table_data_ptr;
      int k=0;
      for (i=1;i<=ns_final;i++) {
        //ht_ptr = hash_rec(&ht,i);
        k++; // skip state value
        for (j=1;j<=ne;j++)
          table[j][i] = values[k++];
      }
    }
    else{
      tmalloc(fsaptr->table->table_data_ptr,int *,ns_final+2);
//      tmalloc(fsaptr->table->table_data_ptr[0],int,ht.tot_space-ns_final);
      tmalloc(fsaptr->table->table_data_ptr[0],int,next_value-2*ns_final);
      table = fsaptr->table->table_data_ptr;
      table[1] = ptr = table[0];
      /*for (i=1;i<=ns_final;i++){
        ht_ptr = hash_rec(&ht,i);
        ptr2 = ht_ptr+1;
        ptr2e = ht_ptr + hash_rec_len(&ht,i) - 1;
        while (ptr2 <= ptr2e)
          *(ptr++) = *(ptr2++);
        table[i+1] = ptr;
      }*/
      j=1;
      for (i=1;i<=ns_final+1;i++){
        table[i]=ptr;
        while (j<next_value) {
          if (values[j]==MAXINT) {
            // end of record marker
            j+=2; //skip record marker and the state number that follows
            break;
          }
          *(ptr++)=values[j++];
        }
      }
    }
  }
  // hash_clear(&ht);
  clear_hash (keytostate, no_states);
  tfree(values);
  tfree(keytostate);
  tfree(block_numa);
  tfree(block_numb);
  tfree(fsaptr->is_accessible);
  if (kbm_print_level>=3)
    printf("    #Number of iterations = %d.\n",num_iterations);
  return 0;
}

int get_image_x(rec,length,next_state,keytostate,no_states)
int *rec;
int length;
int next_state;
kts2 *keytostate;
int no_states;
{
    #define MOD 2039
    int hashval=0;
    struct kts2 *this_kts, *next_kts;
    int i=0;
    while (i<length) {
      hashval=(hashval+(((i+1)*rec[i])%MOD))%no_states;
      i++;
    }
    if (hashval==0) return 0;		
    this_kts=keytostate+hashval;
    while (1) {
        if (this_kts->state==0) {
            //new states
            this_kts->state=next_state+1;
            this_kts->value=rec;
            this_kts->len=length;
            this_kts->next_value=0;
            break;
        }
        if (this_kts->len==length)
          if (memcmp(this_kts->value,rec,length*sizeof(int))==0) 
              break;
        if (this_kts->next_value==0) {
            tmalloc(next_kts,kts2,1);
            this_kts->next_value=next_kts;
            this_kts=next_kts;
            this_kts->state=0;
            this_kts->next_value=0;
            continue;
        }
        this_kts=this_kts->next_value;
    }
    return this_kts->state;
}

void clear_hash (kts2 *keytostate,int no_states)
{
   struct kts2 *this_kts, *next_kts, *to_free_kts;;
   int kts_i=0;	
   while (kts_i<no_states)
   {	
      	this_kts=keytostate+kts_i;
        this_kts->state=0;
	this_kts->value=0;
	to_free_kts=this_kts->next_value;
        this_kts->next_value=0;
	while (to_free_kts !=0)
	{		
		next_kts=to_free_kts->next_value;
		free(to_free_kts);
		to_free_kts=next_kts;
	}
      	kts_i++;
   } 
}
fsa *
fsa_minred1(minred,waptr,op_table_type,destroy,tempfilename)
	fsa *minred;
	fsa *waptr;
	storage_type op_table_type;
	boolean destroy;
	char *tempfilename;
{
  int **watable, **table,  **mintable,ne, nsi, nsi1, ns, dr, *fsarow,
      nt, cstate, csa, csb, im, i, g, len, ct, *ht_ptr;
  boolean dense_ip, dense_op;
  fsa *minred1;
  hash_table ht;
  FILE *tempfile, *fopen();

  if (kbm_print_level>=3)
    printf("    #Calling fsa_minred1.\n");
  if (!waptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_minred only applies to DFA's.\n");
    return 0;
  }

  tmalloc(minred1,fsa,1);
  fsa_init(minred1);
  srec_copy(minred1->alphabet,waptr->alphabet);
  minred1->flags[DFA] = TRUE;
  minred1->flags[ACCESSIBLE] = TRUE;
  minred1->flags[BFS] = TRUE;

  ne = waptr->alphabet->size;
  nsi = waptr->states->size;
  nsi1 = nsi+1;

  watable = waptr->table->table_data_ptr;
  table=watable;
  mintable = minred->table->table_data_ptr;
  //return NULL;
  dense_ip = waptr->table->table_type==DENSE;
  dr = waptr->table->denserows;
  dense_op = op_table_type==DENSE;
  minred1->states->type = SIMPLE;
  minred1->num_initial = 1;
  tmalloc(minred1->initial,int,2);
  minred1->initial[1] = 1;
  minred1->table->table_type = op_table_type;
  minred1->table->denserows = 0;
  minred1->table->printing_format = op_table_type;
  
  hash_init(&ht,TRUE,2,0,0);
  ht_ptr = ht.current_ptr;
  ht_ptr[0] = waptr->initial[1];
  ht_ptr[1] = waptr->initial[1];
  im = hash_locate(&ht,2);
  if (im!=1) {
    fprintf(stderr,"Hash-initialisation problem in fsa_binop.\n");
    return 0;
  }
  if ((tempfile=fopen(tempfilename,"w"))==0){
    fprintf(stderr,"Error: cannot open file %s\n",tempfilename);
    return 0;
  }
  if (dense_op)
    tmalloc(fsarow,int,ne)
  else
    tmalloc(fsarow,int,2*ne+1)
 
  cstate = 0;
  if (dense_op)
    len = ne; /* The length of the fsarow output. */
  nt = 0; /* Number of transitions in minred */

/* A state of *minred1 will be essentially a pair of states of *waptr,
 * with the transitions coming from those of *waptr.
 * which are not accepted by *waptr but whose maximal prefix is,
 * whereas the right hand side will accept words  w which are not accepted
 * by *waptr but whose maximal suffix is.
 * Thus, on the lhs, transitions to 0 in *waptr will go to a new accept-
 * state nsi1 instead (with no transitions from nsi1) whereas on the rhs
 * the first transition will be back to the intiial state.
 * The initial state itself is non-accept.
 */
  while (++cstate <= ht.num_recs) {
    if (kbm_print_level>=3) {
      if ((cstate<=1000 && cstate%100==0)||(cstate<=10000 && cstate%1000==0)||
          (cstate<=100000 && cstate%5000==0) || cstate%50000==0)
       printf("    #cstate = %d;  number of states = %d.\n",cstate,ht.num_recs);
    }
    ht_ptr = hash_rec(&ht,cstate);
    csa = ht_ptr[0]; csb = ht_ptr[1];
    if (!dense_op)
      len = 0;
    for (g=1;g<=ne;g++) {
/* Calculate action of generator g on state cstate */
      ht_ptr = ht.current_ptr;
      if (csa==0 || csa==nsi1)
         ht_ptr[0] = 0;
      else {
        ht_ptr[0] = target(dense_ip,watable,g,csa,dr);
 //       if (ht_ptr[0]==0)
 //	        ht_ptr[0] = nsi1;
      }
      if (cstate==1)
        ht_ptr[1] = csb;
      else {
        ht_ptr[1] = csb ? target(dense_ip,mintable,g,csb,dr) : 0;
// check if ht_ptr[1] is an accept state
	int i_accepting = minred->num_accepting;
	int image=ht_ptr[1];
	while (i_accepting > 0) {
		if (image == minred->accepting[i_accepting]) {
			ht_ptr[0] = nsi1;
			break;
		}
		i_accepting--;
	}
      }
      if (ht_ptr[0]==0 || ht_ptr[1]==0)
        im = 0;
      else
        im = hash_locate(&ht,2);
      if (dense_op)
         fsarow[g-1] = im;
      else if (im>0) {
         fsarow[++len] = g;
         fsarow[++len] = im;
      }
      if (im>0)
        nt++;
    }
    if (!dense_op)
      fsarow[0] = len++;
    fwrite((void *)fsarow,sizeof(int),(size_t)len,tempfile);
  }
  fclose(tempfile);

  ns = minred1->states->size = ht.num_recs;
  minred1->table->numTransitions = nt;
/* Locate the accept states: first count them and then record them. */
  ct = 0;
  for (i=1;i<=ns;i++) {
    ht_ptr = hash_rec(&ht,i);
    if (ht_ptr[0]==nsi1)
        ct++;
  }
  minred1->num_accepting = ct;
  if (ct==1 || ct != ns) {
    tmalloc(minred1->accepting,int,ct+1);
    ct = 0;
    for (i=1;i<=ns;i++) {
      ht_ptr = hash_rec(&ht,i);
      if (ht_ptr[0]==nsi1)
        minred1->accepting[++ct] = i;
    }
  }
  hash_clear(&ht);
  tfree(fsarow);
  tfree(waptr->is_accepting);

  if (destroy) 
    fsa_clear(waptr);

/* Now read the transition table back in */
  tempfile = fopen(tempfilename,"r");
  compressed_transitions_read(minred1,tempfile);
  fclose(tempfile);

  unlink(tempfilename);

  return minred1;
}

fsa * fsa_beginswith(waptr,beginstring)
	fsa *waptr;
	char *beginstring;
{
	int **watable;
  	char **gpletters=waptr->alphabet->names;
  	int ne = waptr->alphabet->size;
	gen g;
	int ws=1;
	int lhsi=0;
	int len=strlen(beginstring);
	fsa *wa_beginswith;
 	tmalloc(wa_beginswith,fsa,1);
    	fsa_copy(wa_beginswith,waptr);
  	watable = wa_beginswith->table->table_data_ptr;
	while (lhsi<len)
	{
		int target=0;
		int i=1;

		while (i<=ne) {
			if (beginstring[lhsi]==gpletters[i][0]) {
				target=i;	
				break;
			}
			i++;
		}
		if (target==0)
		{
			printf("beginswith target invalid\n");
			fsa_clear(wa_beginswith);
			return waptr;
		}
		for (g=1;g<=ne;g++) { 
		    if (g !=target)
			watable[g][ws]=0;
		}
		ws=watable[target][ws];
		if (target==0)
		{
			printf("beginswith target invalid\n");
			fsa_clear(wa_beginswith);
			return waptr;
		}
		lhsi++;
	}	
	return wa_beginswith;
}

// MAF originated  simplification of fsa_wa which causes smaller intermediate
// wa files to be created. For further detail see kbmag/lib/fsawa.c and
// maf/mafauto.cpp from line 274 
fsa * fsa_wa_x (fsaptr,op_table_type,tempfilename,geodesic,hashlimit)
	fsa *fsaptr;
	storage_type op_table_type;
	char *tempfilename;
	boolean geodesic;
	int hashlimit;
{ int  ***dtable, ne, ngens, ndiff, ns, *fsarow, nt, cstate, cs, csdiff, csi,
       im, i, k, g1, g2, len, identity; 
  unsigned short int *ht_ptr, *ht_ptrb, *ht_ptre, *cs_ptr, *cs_ptre, *ptr;
  boolean dense_op, no_trans, good;
  char *cf;
  short_hash_table ht;
  fsa *wa;
  FILE *tempfile, *fopen();
  int SEEN_LHS_BETTER = 1;
  int SEEN_RHS_BETTER = 2;
  int SEEN_EQUAL = 3;

//void short_hash_init();
//int short_hash_locate();
//void short_hash_clear();
//unsigned short* short_hash_rec();

  if (!fsaptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_wa only applies to DFA's.\n");
    return 0;
  }

  if (fsaptr->alphabet->type != PRODUCT || fsaptr->alphabet->arity != 2) {
    fprintf(stderr,"Error in fsa_wa: fsa must be 2-variable.\n");
    return 0;
  }
  if (fsaptr->states->type != WORDS) {
    fprintf(stderr,"Error in fsa_wa: fsa must be word-difference type.\n");
    return 0;
  }

  tmalloc(wa,fsa,1);
  fsa_init(wa);
  srec_copy(wa->alphabet,fsaptr->alphabet->base);
  wa->flags[DFA] = TRUE;
  wa->flags[ACCESSIBLE] = TRUE;
  wa->flags[BFS] = TRUE;

  ne = fsaptr->alphabet->size;
  ngens = wa->alphabet->size;
  ndiff = fsaptr->states->size;

  if (ne != (ngens+1)*(ngens+1)-1) {
   fprintf(stderr,
       "Error: in a 2-variable fsa, alphabet size should = ngens^2 - 1.\n");
    return 0;
  }

  identity = fsaptr->accepting[1]; /* assumed to be unique */
  if (fsaptr->num_accepting!=1 || (identity != fsaptr->initial[1])) {
    fprintf(stderr,"Error: Input to fsa_wa not a word-difference machine.\n");
    return 0;
  }

  if (fsaptr->table->table_type!=DENSE) {
     fprintf(stderr,
      "Error: function fsa_wa can only be called with a densely-stored fsa.\n");
     return 0;
  }
  dense_op = op_table_type==DENSE;
  if (fsa_table_dptr_init(fsaptr)== -1) return 0;
  dtable = fsaptr->table->table_data_dptr;

  wa->states->type = SIMPLE;
  wa->num_initial = 1;
  tmalloc(wa->initial,int,2);
  wa->initial[1] = 1;
  wa->table->table_type = op_table_type;
  wa->table->denserows = 0;
  wa->table->printing_format = op_table_type;
  
  short_hash_init(&ht,FALSE,0,0,0);
  ht_ptr = ht.current_ptr;
  ht_ptr[0] = identity;
  im = short_hash_locate(&ht,1);
  if (im== -1) return 0;
  if (im!=1) {
    fprintf(stderr,"Hash-initialisation problem in fsa_wa.\n");
     return 0;
  }
  if ((tempfile=fopen(tempfilename,"w"))==0){
    fprintf(stderr,"Error: cannot open file %s\n",tempfilename);
     return 0;
  }
  if (dense_op)
    tmalloc(fsarow,int,ngens)
  else
    tmalloc(fsarow,int,2*ngens+1)
 
  cstate = 0;
  if (dense_op)
    len = ngens; /* The length of the fsarow output. */
  nt = 0; /* Number of transitions in exists */

  //tmalloc(cfcount,int,ndiff+1);
  //for (i=1;i<=ndiff;i++)
   //   cfcount[i] = 0;
  tmalloc(cf,char,ndiff+1);
  int dollarcount=0;
  int total_hashx=sizeof(short int);
  while (++cstate <= ht.num_recs) {
    if (kbm_print_level>=3) {
      if ((cstate<=1000 && cstate%100==0)||(cstate<=10000 && cstate%1000==0)||
          (cstate<=100000 && cstate%5000==0) || cstate%50000==0)
       printf("    #cstate = %d;  number of states = %d.\n",cstate,ht.num_recs);
    } 
    if (hashlimit > 0)
     if ( (cstate<=1000000 && cstate%5000==0) || cstate%50000==0)
       Printf("    #cstate = %d;  number of states = %d, hash space size =%d\n",cstate,ht.num_recs,total_hashx);
    cs_ptr = short_hash_rec(&ht,cstate);
    cs_ptre = short_hash_rec(&ht,cstate) + short_hash_rec_len(&ht,cstate) - 1;
    if (!dense_op)
      len = 0;
    for (g1=1;g1<=ngens;g1++) {
/* Calculate action of generator g1 on state cstate  - to get the image, we
 * have to apply (g1,g2) to each element in the subset corresponding to cstate,
 * and this for each generator g2 of the base-alphabet (including the padding
 * symbol).
 * Since we are excluding words that contain subwords w_1 s.t. (w_1,w_2) is
 * accepted by *fsaptr, we also have to apply (g1,g2) to the initial state
 * of *fsaptr.
 */
      int local_hashx=0;
      for (i=1;i<=ndiff;i++)
        cf[i] = 0;
      ptr = cs_ptr-1;
	if (TRACE5 && g1==1 && cstate<80) {
	printf ("old\n");
      while (ptr <= cs_ptre) {
	int gt_state, old_gt_state;;
        cs = ptr<cs_ptr ? identity : *ptr;
        csdiff = cs%ndiff;
	gt_state=1;
	if (cs>ndiff) gt_state=2;
        if (csdiff==0) csdiff = ndiff;
	printf(" %d-%d\n",csdiff,gt_state);
        ptr++;
	} 
      	ptr = cs_ptr-1; }
   /* csdiff is the state of *fsaptr corresponding to cs */
      no_trans = FALSE;
/* We will set no_trans to be true if we find that the transition leads to
 * failure.
 */
      while (ptr <= cs_ptre) {
/* We add the initial state of *fsaptr to the subset representing cstate */
	int gt_state, old_gt_state;;
        cs = ptr<cs_ptr ? identity : *ptr;
        csdiff = cs%ndiff;
        if (csdiff==0) csdiff = ndiff;
   /* csdiff is the state of *fsaptr corresponding to cs */
        ptr++;
	// MAF algorithm inserted here
	old_gt_state=SEEN_LHS_BETTER;
	if (cs>ndiff)
		old_gt_state=SEEN_RHS_BETTER;
	if (cs>2*ndiff)
		old_gt_state=SEEN_EQUAL;
        if (csdiff == identity) { 
		old_gt_state=SEEN_EQUAL;
	        for (g2=1;g2<ngens+1;g2++){
 			csi =  dense_dtarget(dtable,g1,g2,csdiff);
			if (csi==0)
       		     		continue;
 			gt_state = g2 < g1 ? SEEN_RHS_BETTER : SEEN_LHS_BETTER;
			if (geodesic)
				gt_state = SEEN_EQUAL;

          		if (csi==identity)
			{
			 if (gt_state == SEEN_RHS_BETTER)
			 {
            			no_trans = TRUE;
            			break;
			 }
			 continue;
			}
			if (g1 == g2)
				gt_state = SEEN_EQUAL;
            		if ( gt_state > cf[csi])
              			cf[csi] = gt_state;
		}
	}
	else {
		for (g2=1;g2<ngens+1;g2++){
		  csi =  dense_dtarget(dtable,g1,g2,csdiff);
		  if (csi==0)
		    continue;
		  gt_state=old_gt_state;	
		  if (csi==identity) {
           		if (gt_state == SEEN_RHS_BETTER)
			{
			      no_trans = TRUE;
			      break;
		    	}
			continue;
		  }
            	  if ( gt_state > cf[csi])
              		cf[csi] = gt_state;
		}
	}
        if (no_trans)
          break;
	/* now deal with pad on the right */
	  csi =  dense_dtarget(dtable,g1,g2,csdiff);
          if (csi==identity)
	  {
            	no_trans = TRUE;
            	break;
	  }
	  else {if (csi)
		    cf[csi] = SEEN_RHS_BETTER;
	       else if (WADIAG)
		    //printf("%d-%d,%d>?\n",csdiff,g1,g2);
		    dollarcount++;
	  }
      } // for ptr
      if (hashlimit && total_hashx>hashlimit) {
        // crude truncation!!!
	no_trans=FALSE;
      }
      if (no_trans) {
        if (dense_op)
          fsarow[g1-1] = 0;
        continue;
      }
/* Now we have the image stored in the array cf, and we translate it to a list
 * and insert it into the hash-table.
 */
	if (TRACE5 && cstate< 80) printf("%d****\n",cstate);
      ht_ptrb = ht.current_ptr;
      ht_ptre = ht_ptrb-1;
      for (i=1;i<=ndiff;i++) {
        k = cf[i];
        if (k>0) {
                local_hashx++;
        }
	if (TRACE5 && cstate < 80 && cf[i])
		printf(" %d-%d\n",i,cf[i]);
        if (k==SEEN_LHS_BETTER)
          *(++ht_ptre) = i;
        else if (k==SEEN_RHS_BETTER)
	{	             
	  //cfcount[i]++;
          *(++ht_ptre) = ndiff+i;
	}
        else if (k==SEEN_EQUAL)
          *(++ht_ptre) = ndiff+ndiff+i;
      }
      if (hashlimit && total_hashx>hashlimit) {
        // crude truncation!!!
	im=cstate;;
      }
      else
      	im = short_hash_locate(&ht,ht_ptre-ht_ptrb+1);
      if (im== -1) return 0;
	if (TRACE5 && cstate<80) printf("%d->%d\n",g1,im);
      if (dense_op)
         fsarow[g1-1] = im;
      else if (im>0) {
         fsarow[++len] = g1;
         fsarow[++len] = im;
      }
      if (im>0) {
        nt++;
        if (im==ht.num_recs) {
       		 if (!(hashlimit && total_hashx>hashlimit)) {
                        total_hashx+=(local_hashx*sizeof(short int));
               	} 
        }
      }
    }
    if (!dense_op)
      fsarow[0] = len++;
    fwrite((void *)fsarow,sizeof(int),(size_t)len,tempfile);
  }
  fclose(tempfile);

  short_hash_clear(&ht);
  tfree(fsarow);
  tfree(cf);
  if (WADIAG)
	printf("%d dollar calculations\n",dollarcount);

  ns = wa->states->size = ht.num_recs;
  wa->table->numTransitions = nt;

/* All states of wa will be accept states. */
  wa->num_accepting = ns;
  if (ns==1) {
    tmalloc(wa->accepting,int,2);
    wa->accepting[1] = 1;
  }
  tfree(fsaptr->is_accepting);

/* Now read the transition table back in */
  tempfile = fopen(tempfilename,"r");
  compressed_transitions_read(wa,tempfile);
  fclose(tempfile);

  unlink(tempfilename);

  //for (i=1;i<=ndiff;i++)
  //{
   //   Printf("%d.  %d\n",i,cfcount[i]);
  //}
  return wa;
}
void analyse_diff2(char * gpname, int max_size,char *diff2str)
{
  fsa *diff2, *diff1c, *diff2c;
  char inf1 [100];
  strcpy(inf1,gpname);
  strcat(inf1,".diff2");
  char inf2 [100];
  strcpy(inf2,gpname);
  strcat(inf2,".diff1c");
  char inf3 [100];
  strcpy(inf3,gpname);
  strcat(inf3,".diff2c");
  char fsaname [100];
  if (diff2str) {
	strcpy(inf2,diff2str);
	strcpy(inf3,diff2str);
  }

  if ((rfile = fopen(inf1,"r")) == 0) {
        fprintf(stderr,"Cannot open file %s.\n",inf1);
          exit(1);
   }
  tmalloc(diff2,fsa,1);
  fsa_read(rfile,diff2,DENSE,0,0,TRUE,fsaname);
  fclose(rfile);

  if ((rfile = fopen(inf2,"r")) == 0) {
        fprintf(stderr,"Cannot open file %s.\n",inf2);
          exit(1);
   }
  tmalloc(diff1c,fsa,1);
  fsa_read(rfile,diff1c,DENSE,0,0,TRUE,fsaname);
  fclose(rfile);

  if ((rfile = fopen(inf3,"r")) == 0) {
        fprintf(stderr,"Cannot open file %s.\n",inf3);
          exit(1);
   }
  tmalloc(diff2c,fsa,1);
  fsa_read(rfile,diff2c,DENSE,0,0,TRUE,fsaname);
  fclose(rfile);
  printf("diff2, diff1c and diff2c opened for %s \n",gpname);
  int diff2_size = diff2->states->size;
  int diff1c_size = diff1c->states->size;
  int diff2c_size = diff2c->states->size;
  printf("sizes %d, %d, %d\n",diff2_size,diff1c_size,diff2c_size);
  printf(
	"diff2 states are in diff1 (1), in diff2 only (2) or in neither (0)\n");
  printf("diff2 state, 2 or 0, 0 count - all other states are 1\n"); 
  gen **diff2_list=diff2->states->words;
  gen **diff1c_list=diff1c->states->words;
  gen **diff2c_list=diff2c->states->words;
  int i,j;
  int count0=0;
  int size_limit=2;
  int diff_count;
  int diffc_count;
  while (size_limit <= max_size)
  {
	diff_count=0;
	diffc_count=0;
  	for (i=2;i<= diff2c_size;i++) {
		if (genstrlen(diff2c_list[i])==size_limit)
			diffc_count++;
	}
		if (diffc_count > 0)
		{
  		for (i=2;i<= diff2_size;i++) {
			if (genstrlen(diff2_list[i])==size_limit)
				diff_count++;
		}
		printf("size %d wds,  %d, %d\n",size_limit,diffc_count, diff_count);
	}
	size_limit++;
  }
	
  
  for (i=2;i<= diff2_size;i++) {
	gen * wd=diff2_list[i];
        boolean found=FALSE;
	for (j=2;j<=diff1c_size;j++)
	   if (!genstrcmp(wd,diff1c_list[j])) {
		found=TRUE;
		break;
	   }
	if (found)
	  continue;
	for (j=2;j<=diff2c_size;j++)
	   if (!genstrcmp(wd,diff2c_list[j])) {
		printf("%d 2\n",i);
		found=TRUE;
		break;
	   }
	if (found)
	   continue;
	printf("%d 0 %d\n",i,++count0);
  }
  return;
}

fsa* find_geo_wds(fsa *diff2, fsa *gpwa, char * gpname, char *tempfilename,
			fsa *geowa,int hash_table_size,int counter,int state_limit)
{
	fsa *geopairs, *fsaexists;
	/*printf("calling fsa_wa_x on %s.diff2\n",gpname);
	fsa *geowa=fsa_wa_x(diff2,DENSE,tempfilename,TRUE);
	if (kbm_print_level>1)
		printf(
		"  #Number of states of geowa before minimisation = %d.\n",
	   		 geowa->states->size);
	fsa_minimize(geowa);
	if (kbm_print_level>1)
	printf("  #Number of states of geowa after minimisation = %d.\n",
    		geowa->states->size);
	fsa_copy(*geowa_ptr,geowa); */
	if (hash_table_size > 0) {
                Printf(
		   "calling fsa_geopairs_big_hash on  gpwa and %s.diff2\n",
			gpname);
		geopairs =
           	   fsa_geopairs_big_hash (gpwa,diff2,FALSE,
                            tempfilename,hash_table_size,
                                counter,state_limit);
	}
	else {
                printf(
		   "calling fsa_geopairs on gpwa and %s.diff2\n",
			gpname);
		geopairs =
              		fsa_geopairs (gpwa,diff2,SPARSE,FALSE,
                            tempfilename,TRUE);
	}
	if (kbm_print_level>1)
		printf(
		"  #Number of states of geopairs before minimisation = %d.\n",
	   		 geopairs->states->size);
	fsa_minimize(geopairs);
	if (kbm_print_level>1)
		printf(
		"  #Number of states of geopairs after minimisation = %d.\n",
    			geopairs->states->size);
if (CALC_GEODIFF) {
      char outfx[100];
      strcpy(outfx,gpname);
      strcat(outfx,".geopairs");
      wfile = fopen(outfx,"w");
      Printf("writing %s\n",outfx);
      fsa_print(wfile,geopairs,outfx);
      fclose(wfile);
}
	Printf("calling fsa_exists on geopairs\n");
        fsaexists = fsa_exists(geopairs,DENSE,TRUE,tempfilename);
	if (kbm_print_level>1)
                printf(
		"  Number of states of fsaexists before minimisation = %d.\n",
                    fsaexists->states->size);
        fsa *and_not;
	if (hash_table_size) {
		Printf("fsa_and_not_first on geowa and fsaexists\n");
		and_not = fsa_and_not_first(
			geowa,fsaexists,DENSE,FALSE,tempfilename,0);
	}
	else {
		Printf("fsa_and_not on geowa and fsaexists\n");
		and_not = fsa_and_not(
			geowa,fsaexists,DENSE,FALSE,tempfilename);
	}
	fsa_clear(fsaexists);
	return and_not;
}

fsa * fsa_geopairs_big_hash(waptr,diffptr,
           destroy,tempfilename,hash_table_size,counter,state_limit)
	fsa  *waptr, *diffptr;
	boolean destroy;
	char *tempfilename;
        unsigned int hash_table_size;
	boolean counter;
	int state_limit;
	
{
  /* this is a version of fsa_geopairs which uses one large hash table whose
     size is specified by the parameter hash_table_size. This value is also
     equal to the maximum number of states that the fsa geopairs can contain.
     Please consult further the documentation contained in fsageopairs.c
  */
    
  int  **watable, ***difftable, identity, ngens, ngens1,
      ne, ns, *fsarow, nt, cstate, cswa,  csdiff, im, i, e, len, ct,
      hashval, collisions;
  unsigned int *ht_ptr;
  boolean  dense_op;
  fsa *geopairsptr;
  hash_table ht;
  FILE *tempfile, *fopen();
  gen g1, g2;
  reduction_struct rs_wd;
  int ht_i=0; 
  mstate3 *statetokey [MAX_EXTRA_MEMORY];
  kts3 *keytostate;
  unsigned int next_state;

  unsigned int no_states = hash_table_size+1;

  if (kbm_print_level>=3)
    printf("    #Calling fsa_geopairs_big_hash\n");

  if (!waptr->flags[DFA] || !diffptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_geopairs only applies to DFA's.\n");

    return 0;
  }
  if (waptr->alphabet->type!=IDENTIFIERS) {
    fprintf(stderr,"Error in fsa_geopairs: an fsa has wrong type.\n");
    return 0;
  }
  if (waptr->num_accepting != waptr->states->size) {
    fprintf(stderr,
       "Error : first fsa should be a word-acceptor.\n");
    return 0;
  }
  if (diffptr->alphabet->type!=PRODUCT || diffptr->alphabet->arity!=2) {
    fprintf(stderr,"Error : second fsa must be 2-variable.\n");
    return 0;
  }
  if (diffptr->states->type!=WORDS) {
    fprintf(stderr,
       "Error : second fsa must be word-difference type.\n");
    return 0;
  }
  if (!srec_equal(diffptr->alphabet->base,waptr->alphabet)) {
    fprintf(stderr,"Error in fsa_geopairs: fsa's alphabet's don't match.\n");
    return 0;
  }
  if (waptr->states->size>=MAXINT || diffptr->states->size>=MAXINT) {
    fprintf(stderr,
       "Error: one of the fsa's has too many states.\n");
    return 0;
  }
  int diff_partition=hash_table_size/diffptr->states->size;

  if (fsa_table_dptr_init(diffptr)==-1) return 0;
  fsa_set_is_accepting(diffptr);
  collisions=0;

  tmalloc(geopairsptr,fsa,1);
  fsa_init(geopairsptr);
  srec_copy(geopairsptr->alphabet,diffptr->alphabet);
  geopairsptr->states->type = SIMPLE;
  geopairsptr->flags[DFA] = TRUE;
  geopairsptr->flags[ACCESSIBLE] = TRUE;
  geopairsptr->flags[BFS] = TRUE;

  ngens = waptr->alphabet->size;
  ngens1 = ngens+1;
  ne = diffptr->alphabet->size;

  if (ne != ngens1*ngens1-1) {
   fprintf(stderr,
       "Error: in a 2-variable fsa, alphabet size should = ngens^2 - 1.\n");
    return 0;
  }
  // allocate hash and data areas
  int state_i=0;
  while (state_i<MAX_EXTRA_MEMORY) {
	statetokey[state_i]=NULL;
	state_i++;
  }	

  tmalloc(statetokey[0],mstate3,no_states);
  tmalloc(keytostate,kts3,no_states);
  for (i=0;i<no_states;i++) {
	keytostate[i].state=0;
	keytostate[i].next_value=0;
   }

  watable = waptr->table->table_data_ptr;
  difftable = diffptr->table->table_data_dptr;

  dense_op = FALSE;

  geopairsptr->states->type = SIMPLE;

  geopairsptr->num_initial = 1;
  tmalloc(geopairsptr->initial,int,2);
  geopairsptr->initial[1] = 1;
  geopairsptr->table->table_type = SPARSE;
  geopairsptr->table->denserows = 0;
  geopairsptr->table->printing_format = SPARSE;
  //  we use the kbmag hash area as a work area
  hash_init(&ht,TRUE,2,0,0);
  ht_ptr = ht.current_ptr;
  ht_ptr[0] = waptr->initial[1];
  ht_ptr[1] = identity = diffptr->initial[1];
  im = hash_locate(&ht,2);
  if (im!=1) {
    fprintf(stderr,"Hash-initialisation problem in fsa_geopairs.\n");
    return 0;
  }
  statetokey[ht_i][1].wa=ht_ptr[0];
  statetokey[ht_i][1].diff=diffptr->initial[1];
  
  if ((tempfile=fopen(tempfilename,"w"))==0){
    fprintf(stderr,"Error: cannot open file %s\n",tempfilename);
    return 0;
  }
  if (dense_op)
    tmalloc(fsarow,int,ne)
  else
    tmalloc(fsarow,int,2*ne+1)
 
  cstate = 0;
  if (dense_op)
    len = ne; /* The length of the fsarow output. */
  nt = 0; /* Number of transitions in geopairsptr */
  next_state=ht.num_recs;
  if (counter)
      printf("processed/total states\n");
  int local_state=0;
  int tmin1=no_states-1;
  while (++cstate <= next_state) {
    if (cstate%500==0 && counter)
            printf(
            "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%d  %d",
                  cstate,next_state);
    local_state=cstate;
    state_i=0;
    if (cstate>tmin1) {
	state_i=cstate/tmin1;
	local_state=cstate%tmin1;
	if (local_state==0)
	{
		local_state=tmin1;
		state_i--;
	}
    }
    cswa=statetokey[state_i][local_state].wa; 
    csdiff=statetokey[state_i][local_state].diff;

    if (dense_op)
      for (i=0;i<ne;i++) fsarow[i] = 0;
    else
      len = 0;
    e = 0; /* e is the number of the edge corresponding to the pair (g1,g2) */
    for (g1=1;g1<=ngens1;g1++) for (g2=1;g2<=ngens1;g2++) {
	e++;
	if (g1==ngens1||g2==ngens1)
		continue;
	ht_ptr = ht.current_ptr;
	ht_ptr[1] = dense_dtarget(difftable,g1,g2,csdiff);
	if (ht_ptr[1]==0)
		im=0;
	else {
		ht_ptr[0] = dense_target(watable,g2,cswa);
		if (state_limit > 0)
			if (ht_ptr[0]>state_limit)
				ht_ptr[0]=0;
		if (ht_ptr[0]==0)
			im=0;
		else	
			im=get_image2_big_hash (ht_ptr,next_state,
			keytostate,&collisions,no_states,diff_partition);
	}   	
	
	if (im>next_state) {
	   local_state=im;
	    state_i=0;
	   if (im > tmin1) {
		   state_i=im/tmin1;
		   local_state=im%tmin1;
		   if (local_state==0) {
			state_i--;
			local_state=tmin1;
		   }
		   if (state_i > 99)
		   {
			printf("out of memory!\n");
			exit(99);
		   }
		   if (statetokey[state_i]==NULL) {
			printf("exceeded limit of %d states\n",
				state_i*tmin1);
			printf("allocating  memory for %d more states\n",
				tmin1);
			tmalloc(statetokey[state_i],mstate3,no_states);
		   }
	   }
	   next_state=im;
	   statetokey[state_i][local_state].wa=ht_ptr[0];
	   statetokey[state_i][local_state].diff=ht_ptr[1];
	}
	if (im>0) {
	      fsarow[++len] = e;
	      fsarow[++len] = im;
	      nt++;
	}
    }   /* for (g1=1;g1<=ngens1; ... */
    fsarow[0] = len++;
    fwrite((void *)fsarow,sizeof(int),(size_t)len,tempfile);
  }  /*while (++cstate <= ht.num_recs) */
  fclose(tempfile);

  ns=next_state;
  geopairsptr->states->size=ns;
  geopairsptr->table->numTransitions = nt;

/* Now locate the accept states */

  ct = 0;
  for (i=1;i<=ns;i++) {
    // the states are spread across statetokey tables each holding 
    // no_states -1 states  starting from index 1
    state_i=0;
    int local_i=i;
    if (i>tmin1)
    {
	    state_i= i /tmin1;
	    local_i= i%tmin1;
	    if (local_i==0) {
		state_i--;
		local_i=tmin1;
	    }
    }
    mstate3 mm = statetokey[state_i][local_i];
    if (diffptr->is_accepting[mm.diff])
              ct++;
  }
  geopairsptr->num_accepting = ct;
  if (ct==1 || ct != ns) {
    tmalloc(geopairsptr->accepting,int,ct+1);
    ct = 0;
    for (i=1;i<=ns;i++) {
    // the states are spread across statetokey tables each holding 
    // no_states -1 states  starting from index 1
	    state_i=0;
	    int local_i=i;
	    if (i>tmin1)
	    {
		    state_i= i /tmin1;
		    local_i= i%tmin1;
		    if (local_i==0) {
			state_i--;
			local_i=tmin1;
		    }
	    }
	    mstate3 mm = statetokey[state_i][local_i];
	    if (diffptr->is_accepting[mm.diff])
             geopairsptr->accepting[++ct] = i;
    }
  }
// free the big_hash areas
   state_i=0;
   while (state_i<MAX_EXTRA_MEMORY) {
	if (statetokey[state_i]!=NULL)
	   tfree(statetokey[state_i]);
	state_i++;
   }
   float no_hash_entries=0;
   float no_next_values=0; 
   struct kts3 *this_kts, *next_kts, *to_free_kts;;
   int kts_i=0;	
   while (kts_i<no_states)
   {	
	this_kts=keytostate+kts_i;
	if (this_kts->state)
		no_hash_entries++;
	to_free_kts=this_kts->next_value;
	while (to_free_kts !=0)
	{		
		no_next_values++;
		next_kts=to_free_kts->next_value;
		free(to_free_kts);
		to_free_kts=next_kts;
	}
	kts_i++;
   } 
   tfree(keytostate);
   double depth=
	(no_next_values+no_hash_entries)/no_hash_entries;
      if (DISPLAY_HASH_DEPTH)
  printf("average hash depth (1 is perfect) = (%0.0f + %0.0f)/%0.0f = %0.2f\n",
       no_hash_entries,no_next_values,no_hash_entries,depth);

  hash_clear(&ht);
  tfree(fsarow);
  if (destroy) {
     fsa_clear(waptr); fsa_clear(diffptr);
  }
/* Now read the transition table back in */
  tempfile = fopen(tempfilename,"r");
  compressed_transitions_read(geopairsptr,tempfile);
  fclose(tempfile);

  unlink(tempfilename);

  return geopairsptr;
}

int get_image2_big_hash(rec,next_state,keytostate,
			collisions,no_states,diff_partition)
    mstate3 *rec;
    int next_state;
    kts3 *keytostate;
    unsigned int *collisions;
    unsigned int no_states;
    int diff_partition;
{
    unsigned int hashval=0;
    unsigned int m, hv, i, k;
    struct kts3 *this_kts, *next_kts;
    hashval=((rec->wa)*101+(rec->diff)*diff_partition)%no_states;
    this_kts=keytostate+hashval;
    while (1) {
        if (this_kts->state==0) {
            //new states
            this_kts->state=next_state+1;
            memcpy(&this_kts->value,rec,sizeof(mstate3));
            this_kts->next_value=0;
            break;
        }
        if (memcmp(&this_kts->value,rec,sizeof(mstate3))==0) 
            break;
        if (this_kts->next_value==0) {
            (*collisions)++;
            tmalloc(next_kts,kts3,1);
            this_kts->next_value=next_kts;
            this_kts=next_kts;
            this_kts->state=0;
            continue;
        }
        this_kts=this_kts->next_value;
    }
    return this_kts->state;
}

// version of fsa_binop where and_not causes words to be accepted only if   
// prefixes of accepted words are not accepted. This implements the 
// and_not_first functionality found in MAF.

typedef enum {AND, OR, AND_NOT} kbm_binop;

fsa * fsa_binop_x(fsaptr1,fsaptr2,op_table_type,destroy,tempfilename,op,labels,state_limit)
        fsa *fsaptr1, *fsaptr2;
        storage_type op_table_type;
        boolean destroy;
        char *tempfilename;
        kbm_binop op;
        boolean labels;
        int state_limit;
{
  int **table1, **table2, ne, ns, dr1, dr2, *fsarow,
      nt, cstate, csa, csb, im, i, g, len, ct, *ht_ptr;
  boolean dense_ip1, dense_ip2, dense_op, accept;
  fsa *and_or_not;
  hash_table ht;
  setToLabelsType *lab;
  FILE *tempfile, *fopen();

  if (kbm_print_level>=3)
    printf("    #Calling fsa_binop_x.\n");
  if (!fsaptr1->flags[DFA] || !fsaptr2->flags[DFA]){
    fprintf(stderr,"Error: fsa_binop only applies to DFA's.\n");
    return 0;
  }
 if (!srec_equal(fsaptr1->alphabet,fsaptr2->alphabet)) {
    fprintf(stderr,"Error in fsa_binop: fsa's have different alphabets.\n");
    return 0;
  }

  if (fsaptr1->flags[RWS])
    fsa_clear_rws(fsaptr1);

  if (fsaptr2->flags[RWS])
    fsa_clear_rws(fsaptr2);

  tmalloc(and_or_not,fsa,1);
  if (fsaptr2->num_initial==0 && (op==AND_NOT || op==OR)) {
    fsa_copy(and_or_not,fsaptr1);
    and_or_not->table->printing_format = op_table_type;
    if (destroy) {
      fsa_clear(fsaptr1); fsa_clear(fsaptr2);
    }
    return and_or_not;
  }
  if (fsaptr1->num_initial==0 && op==OR) {
    fsa_copy(and_or_not,fsaptr2);
    and_or_not->table->printing_format = op_table_type;
   if (fsaptr1->states->type == LABELED) {
      srec_clear(and_or_not->states);
      srec_copy(and_or_not->states,fsaptr1->states);
      tfree(and_or_not->states->setToLabels);
      ns=fsaptr2->states->size;
      tmalloc(and_or_not->states->setToLabels,setToLabelsType,ns+1);
      for (i=0;i<=ns;i++)
        and_or_not->states->setToLabels[i]=0;
    }
    if (destroy) {
      fsa_clear(fsaptr1); fsa_clear(fsaptr2);
    }
    return and_or_not;
  }

  fsa_init(and_or_not);
  srec_copy(and_or_not->alphabet,fsaptr1->alphabet);
  and_or_not->flags[DFA] = TRUE;
  and_or_not->flags[ACCESSIBLE] = TRUE;
  and_or_not->flags[BFS] = TRUE;

  if (labels) {
    if (fsaptr1->states->type != LABELED) {
 fprintf(stderr,"Error in fsa_binop: first fsa should have labeled states.\n");
      return 0;
    }
    srec_copy(and_or_not->states,fsaptr1->states);
    tfree(and_or_not->states->setToLabels);
    lab = fsaptr1->states->setToLabels;
  }
  else and_or_not->states->type = SIMPLE;

  and_or_not->table->table_type = op_table_type;
  and_or_not->table->denserows = 0;
  and_or_not->table->printing_format = op_table_type;

  if (fsaptr1->num_initial==0 || fsaptr2->num_initial==0) {
    and_or_not->num_initial = 0;
    and_or_not->num_accepting = 0;
    and_or_not->states->size = 0;
    if (destroy) {
      fsa_clear(fsaptr1); fsa_clear(fsaptr2);
    }
    return and_or_not;
  }
 ne = fsaptr1->alphabet->size;

  fsa_set_is_accepting(fsaptr1);
  fsa_set_is_accepting(fsaptr2);
  table1 = fsaptr1->table->table_data_ptr;
  table2 = fsaptr2->table->table_data_ptr;

  dense_ip1 = fsaptr1->table->table_type==DENSE;
  dense_ip2 = fsaptr2->table->table_type==DENSE;
  dr1 = fsaptr1->table->denserows;
  dr2 = fsaptr2->table->denserows;
  dense_op = op_table_type==DENSE;

  and_or_not->num_initial = 1;
  tmalloc(and_or_not->initial,int,2);
  and_or_not->initial[1] = 1;

  hash_init(&ht,TRUE,2,0,0);
  ht_ptr = ht.current_ptr;
  ht_ptr[0] = fsaptr1->initial[1];
  ht_ptr[1] = fsaptr2->initial[1];
  im = hash_locate(&ht,2);
  if (im!=1) {
  fprintf(stderr,"Hash-initialisation problem in fsa_binop.\n");
    return 0;
  }
  if ((tempfile=fopen(tempfilename,"w"))==0){
    fprintf(stderr,"Error: cannot open file %s\n",tempfilename);
    return 0;
  }
  if (dense_op)
    tmalloc(fsarow,int,ne)
  else
    tmalloc(fsarow,int,2*ne+1)

  cstate = 0;
  if (dense_op)
    len = ne; /* The length of the fsarow output. */
  nt = 0; /* Number of transitions in and_or_not */

  while (++cstate <= ht.num_recs) {
    if (kbm_print_level>=3) {
      if ((cstate<=1000 && cstate%100==0)||(cstate<=10000 && cstate%1000==0)||
          (cstate<=100000 && cstate%5000==0) || cstate%50000==0)
       printf("    #cstate = %d;  number of states = %d.\n",cstate,ht.num_recs);
    }
   ht_ptr = hash_rec(&ht,cstate);
    csa = ht_ptr[0]; csb = ht_ptr[1];
    if (op==AND_NOT)
	// AND_NOT_FIRST functionality - prevent  accepted words
	// with  accepted prefixes
	if (fsaptr1->is_accepting[csa] && !fsaptr2->is_accepting[csb])
		csa=0;
    if (!dense_op)
      len = 0;
    for (g=1;g<=ne;g++) {
/* Calculate action of generator g on state cstate */
      ht_ptr = ht.current_ptr;
      ht_ptr[0] = csa ? target(dense_ip1,table1,g,csa,dr1) : 0;
      if (state_limit>0)
       if (ht_ptr[0] > state_limit)
        ht_ptr[0]=0;
      ht_ptr[1] = csb ? target(dense_ip2,table2,g,csb,dr2) : 0;
      if ( (op==AND && (ht_ptr[0]==0 || ht_ptr[1]==0)) ||
           (op==AND_NOT && ht_ptr[0]==0) )
        im = 0;
      else
        im = hash_locate(&ht,2);
      if (im== -1) return 0;
      if (dense_op)
         fsarow[g-1] = im;
      else if (im>0) {
         fsarow[++len] = g;
         fsarow[++len] = im;
      }
      if (im>0)
        nt++;
   }
    if (!dense_op)
      fsarow[0] = len++;
    fwrite((void *)fsarow,sizeof(int),(size_t)len,tempfile);
  }
  fclose(tempfile);

  ns = and_or_not->states->size = ht.num_recs;
  and_or_not->table->numTransitions = nt;
  if (labels)
    tmalloc(and_or_not->states->setToLabels,setToLabelsType,ns+1);
/* Locate the accept states: first count them and then record them. */
  ct = 0;
  for (i=1;i<=ns;i++) {
    ht_ptr = hash_rec(&ht,i);
    csa = ht_ptr[0]; csb = ht_ptr[1];
    accept =
      op==AND     ? fsaptr1->is_accepting[csa] && fsaptr2->is_accepting[csb] :
      op==OR      ? fsaptr1->is_accepting[csa] || fsaptr2->is_accepting[csb] :
      op==AND_NOT ? fsaptr1->is_accepting[csa] && !fsaptr2->is_accepting[csb] :
          FALSE; /* default cannot occur */
    if (accept)
        ct++;
   if (labels)
      and_or_not->states->setToLabels[i] = csa==0 ? 0 : lab[csa];
  }
  and_or_not->num_accepting = ct;
  if (ct==1 || ct != ns) {
    tmalloc(and_or_not->accepting,int,ct+1);
    ct = 0;
    for (i=1;i<=ns;i++) {
      ht_ptr = hash_rec(&ht,i);
      csa = ht_ptr[0]; csb = ht_ptr[1];
      accept =
      op==AND     ? fsaptr1->is_accepting[csa] && fsaptr2->is_accepting[csb] :
      op==OR      ? fsaptr1->is_accepting[csa] || fsaptr2->is_accepting[csb] :
      op==AND_NOT ? fsaptr1->is_accepting[csa] && !fsaptr2->is_accepting[csb] :
          FALSE; /* default cannot occur */
      if (accept)
        and_or_not->accepting[++ct] = i;
    }
  }
  hash_clear(&ht);
  tfree(fsarow);
  tfree(fsaptr1->is_accepting);
  tfree(fsaptr2->is_accepting);
if (destroy) {
    fsa_clear(fsaptr1); fsa_clear(fsaptr2);
  }
/* Now read the transition table back in */
  tempfile = fopen(tempfilename,"r");
  compressed_transitions_read(and_or_not,tempfile);
  fclose(tempfile);

  unlink(tempfilename);

  return and_or_not;
}

fsa *fsa_and_not_first (fsaptr1,fsaptr2, op_table_type, destroy,tempfilename,state_limit)
        fsa *fsaptr1, *fsaptr2;
        storage_type op_table_type;
        boolean destroy;
        char *tempfilename;
	int state_limit;
{
	return fsa_binop_x(fsaptr1,fsaptr2,op_table_type,
			destroy, tempfilename,AND_NOT,FALSE,state_limit);
}


fsa * fsa_genmult2_int_x(genmultptr,op_table_type)
	fsa *genmultptr;
	storage_type op_table_type;
{
#define LABHTSIZE 8192
  int **table, ne, ngens, ngens1, ns, *fsarow, e, e1, e2, es1, ef1, dr,
      nt, cstate, im, csa, csb, csima, csimb,
      i, j, j1, j2, g1, g2, g3, len, reclen, nlab, ct;
  int *ht_ptr, *ht_chptr, *ht_ptrb, *ht_ptre, *cs_ptr, *cs_ptre, *ptr;
  boolean dense_ip, dense_op, got, leftpad, rightpad, keeptable;
  setToLabelsType *label, l1, l2, *newlabel;
  gen **labellist1, **labellist2;
  hash_table ht, labelht;
  fsa *genmult2ptr;
  srec *labelset;
  FILE *tablefile, *fopen();

  if (kbm_print_level>=3)
    printf("    #Calling fsa_genmult2_int.\n");
  if (!genmultptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_genmult2 only applies to DFA's.\n");
    return 0;
  }

  if (genmultptr->alphabet->type!=PRODUCT || genmultptr->alphabet->arity!=2) {
    fprintf(stderr,"Error in fsa_genmult2: fsa must be 2-variable.\n");
    return 0;
  }

  if (genmultptr->states->type!=LABELED) {
    fprintf(stderr,"Error in fsa_genmult2: states of fsa must be labeled.\n");
    return 0;
  }


  tmalloc(genmult2ptr,fsa,1);
  fsa_init(genmult2ptr);
  srec_copy(genmult2ptr->alphabet,genmultptr->alphabet);
  genmult2ptr->num_initial = 1;
  tmalloc(genmult2ptr->initial,int,2);
  genmult2ptr->initial[1] = 1;
  genmult2ptr->flags[DFA] = TRUE;
  genmult2ptr->flags[ACCESSIBLE] = TRUE;
  genmult2ptr->flags[BFS] = TRUE;

  genmult2ptr->table->table_type = op_table_type;
  genmult2ptr->table->denserows = 0;
  genmult2ptr->table->printing_format = op_table_type;

  ne = genmultptr->alphabet->size;
  ngens = genmultptr->alphabet->base->size;
  ngens1 = ngens+1;

  if (ne != ngens1*ngens1-1) {
   fprintf(stderr,
     "Error: in a 2-variable fsa, alphabet size should = (ngens+1)^2 - 1.\n");
    return 0;
  }

  fsa_set_is_accepting(genmultptr);

  dense_ip = genmultptr->table->table_type==DENSE;
  dr = genmultptr->table->denserows;
  dense_op = op_table_type==DENSE;
  table = genmultptr->table->table_data_ptr;

  hash_init(&ht,FALSE,0,0,0);
  ht_ptr = ht.current_ptr;
  ht_ptr[0] = genmultptr->initial[1];
  ht_ptr[1] = genmultptr->initial[1];
  im = hash_locate(&ht,2);
/* Each state in the composite transition table will be represented as a subset
 * of the set of ordered pairs of states of *genmultptr.
 * The initial state contains just one such pair, of which both components are
 * the initial state of *genmultptr.
 * The subsets will be stored as variable-length records in the hash-table,
 * as a list of pairs in increasing order.
 * If a state is reached from a transition ($,x) or (x,$) (with $ the padding
 * symbol), then it needs to be marked as such, since we do not allow a $
 * to be followed by any other generator.
 * We do this by adding a 1 or a 2 to the end of the statelist - this is
 * recognisable by the fact that the length of the statelist then becomes odd.
 */
  if (im!=1) {
    fprintf(stderr,"Hash-initialisation problem in fsa_genmult2.\n");
    return 0;
  }
  if (dense_op)
    tmalloc(fsarow,int,ne)
  else
    tmalloc(fsarow,int,2*ne+1)
 
  label = genmultptr->states->setToLabels;
  cstate = 0;
  if (dense_op)
    len = ne; /* The length of the fsarow output. */
  nt = 0; /* Number of transitions in genmult2ptr */

  while (++cstate <= ht.num_recs) {
    if (kbm_print_level>=3) {
      if ((cstate<=1000 && cstate%100==0)||(cstate<=10000 && cstate%1000==0)||
          (cstate<=100000 && cstate%5000==0) || cstate%50000==0)
       printf("    #cstate = %d;  number of states = %d.\n",cstate,ht.num_recs);
    }
    reclen = hash_rec_len(&ht,cstate);
    cs_ptr = hash_rec(&ht,cstate);
    cs_ptre = hash_rec(&ht,cstate) + reclen - 1;
    if (reclen % 2 == 1) {
      if (*cs_ptre==1) {
        leftpad=TRUE;  rightpad=FALSE;
      } else {
        leftpad=FALSE;  rightpad=TRUE;
      }
      cs_ptre--;
    }
    else {
      leftpad=FALSE; rightpad=FALSE;
    }
      
    if (dense_op)
      for (i=0;i<ne;i++) fsarow[i] = 0;
    else
      len = 0;

    e = 0; /* e is the edge number of generator pair (g1,g3) */
    for (g1=1;g1<=ngens1;g1++)  for (g3=1;g3<=ngens1;g3++) {
/* Calculate action of pair (g1,g3) on state cstate  - to get the image, we
 * have to apply ( (g1,g2), (g2,g3) ) to each ordered pair in the subset
 * corresponding to cstate, * and this for each generator g2 of the
 * base-alphabet (including the padding symbol).
 */
       
      e++;
      if (g1==ngens1 && g3==ngens1)
        continue;

      if ((leftpad && g1<=ngens) || (rightpad && g3<=ngens))
        continue;
      ht_ptrb = ht.current_ptr;
      ht_ptre = ht_ptrb-1;

      es1 = (g1-1)*ngens1 + 1;
      ef1 = g1*ngens1;
/* As g2 ranges from 1 to ngens+1 in the pair (g1,g2), for fixed g1, the
 * corresponding edge number in the fsa ranges from es1 to ef1.
 *
 * As g2 ranges from 1 to ngens+1 in the pair (g2,g3), for fixed g3, the
 * corresponding edge number in the fsa ranges from g3 upwards, with
 * increments of size ngens+1.
 */

      ptr = cs_ptr;
      while (ptr <= cs_ptre) {
        csa = *(ptr++); csb = *(ptr++);
        for (g2=1,e1=es1,e2=g3; e1<=ef1; g2++,e1++,e2+=ngens1){
          csima = g1==ngens1 && g2==ngens1 ?
		   (genmultptr->is_accepting[csa] ? csa : 0) :
                   target(dense_ip,table,e1,csa,dr);
          if (csima==0)
            continue;

          csimb = g2==ngens1 && g3==ngens1 ?
                   (genmultptr->is_accepting[csb] ? csb : 0) :
                   target(dense_ip,table,e2,csb,dr);
          if (csimb==0)
            continue;
          if (ht_ptrb>ht_ptre || csima> *(ht_ptre-1) ||
                               (csima== *(ht_ptre-1) &&  csimb > *ht_ptre) ){
/* We have a new state for the image subset to be added to the end */
            *(++ht_ptre) = csima; *(++ht_ptre) = csimb;
          }
          else {
            ht_chptr = ht_ptrb;
            while (*ht_chptr < csima)
              ht_chptr+=2;
            while (*ht_chptr==csima && *(ht_chptr+1)<csimb)
                ht_chptr+=2;
            if (csima < *ht_chptr || csimb < *(ht_chptr+1)) {
/* we have a new state for the image subset to be added in the middle */
              ht_ptr = ht_ptre;
              ht_ptre += 2;
              while (ht_ptr >= ht_chptr) {
                *(ht_ptr+2) = *ht_ptr;
                ht_ptr--;
              }
              *ht_chptr = csima; *(ht_chptr+1) = csimb;
            }
          }
        }
      }
      if (ht_ptre > ht_ptrb) {
        if (g1==ngens1)
          *(++ht_ptre) = 1;
        else if (g3==ngens1)
          *(++ht_ptre) = 2;
      }
      im = hash_locate(&ht,ht_ptre-ht_ptrb+1);
      if (im > 0) {
        if (dense_op)
          fsarow[e-1] = im;
        else  {
         fsarow[++len] = e;
         fsarow[++len] = im;
        }
        nt++;
      }
    }
    if (!dense_op)
      fsarow[0] = len++;
    if (keeptable)
      fwrite((void *)fsarow,sizeof(int),(size_t)len,tablefile);
  }

  if (keeptable)
    fclose(tablefile);
  tfree(fsarow);

  ns =
  genmult2ptr->states->size = ht.num_recs;
  genmult2ptr->table->numTransitions = nt;

  if (kbm_print_level>=3) {
    printf("    #Calculated transitions - %d states, %d transitions.\n",ns,nt);
    printf("    #Now calculating state labels.\n");
  }

/* Locate the accept states for  Mult_(a*b)  each generator pair (a,b).
 * This is slightly cumbersome, since a state S
 * is an accept state if either the  subset representing S contains a
 * pair (s1,s2), where s1 is accept state for Mult_a and s2 for Mult_b,
 * or we can get from to such a state by applying ( ($,g), (g,$) ) to one
 * of the pairs in S, for some generator g.
 * It is most convenient to work out this information taking the states
 * S in reverse order.
 * The information on the accept-states will be stored as labels, which
 * (13/9/95 - change labels from words to lists of words)
 * will be lists of words in the generators. Each such list will have the form
 * [a1*b1,a2*b2,...,ar*br], where the (ai,bi) are the generator pairs for
 * which that particular state is an accept state. The labels will be
 * collected initially in a new hash-table.
 * The identity generator will be numbered ngens+1 and given name '_'
 * rather than represented by the empty word. This is so that we can
 * distinguish between multipliers for "_*a" and "a*_" if necessary.
 */

  genmult2ptr->states->type = LABELED;
  tmalloc(genmult2ptr->states->labels,srec,1);
  labelset = genmult2ptr->states->labels;
  labelset->type = LISTOFWORDS;
  labelset->alphabet_size = ngens1;
  for (i=1;i<=ngens;i++) {
    tmalloc(labelset->alphabet[i],char,
                        stringlen(genmultptr->states->labels->alphabet[i])+1);
    strcpy(labelset->alphabet[i],genmultptr->states->labels->alphabet[i]);
  }
  tmalloc(labelset->alphabet[ngens1],char,2);
  labelset->alphabet[ngens1][0]='_';
  labelset->alphabet[ngens1][1]=0;
  tmalloc(genmult2ptr->states->setToLabels,setToLabelsType,ns+1);
  newlabel = genmult2ptr->states->setToLabels;
  tmalloc(genmult2ptr->is_accepting,boolean,ns+1);
  for (i=1;i<=ns;i++)
    genmult2ptr->is_accepting[i] = FALSE;
  hash_init(&labelht,FALSE,0,LABHTSIZE,2*LABHTSIZE);

  es1 = ngens*ngens1 + 1;
  ef1 = ngens1*ngens1-1;
  for (cstate=ns;cstate>=1;cstate--) {
/* We work backwards through the states, since we wish to add on additional
 * elements at the end of the list in the hash-table - this destroys
 * later entries, but that doesn't matter, since we have already dealt
 * with them.
 */
    cs_ptr = hash_rec(&ht,cstate);
    reclen = hash_rec_len(&ht,cstate);
    if (reclen % 2 ==1)
      reclen--;
       /* The last entry only marks the fact that this is a "padded state"*/
    cs_ptre = hash_rec(&ht,cstate) + reclen - 1;
/* Apply generators ( ($,g2), (g2,$) ) and see if we get anything new.
 * We won't bother about having them in increasing order this time.
 */
    ptr = cs_ptr;
    while (ptr <= cs_ptre) {
      csa = *(ptr++); csb = *(ptr++);
      for (e1=es1,e2=ngens1; e1<=ef1; e1++,e2+=ngens1){
        csima = target(dense_ip,table,e1,csa,dr);
        if (csima==0)
          continue;
        csimb = target(dense_ip,table,e2,csb,dr);
        if (csimb==0)
          continue;
  /* see if (csima,csimb) is new */
        ht_chptr = cs_ptr;
        got = FALSE;
        while (ht_chptr < cs_ptre) {
          if (csima == ht_chptr[0] && csimb == ht_chptr[1]) {
            got = TRUE;
            break;
          }
          ht_chptr+=2;
        }
        if (!got) {
  /* add (csima,csimb) to the end */
          *(++cs_ptre) = csima; *(++cs_ptre) = csimb;
        }
      }
    }
  /* Now we see which pairs in the subset are of form (s,t), where s is
   * an accept state for a generator a, and t for a generator b.
   * The list of all such pairs (a,b) will form the label of the state, which
   * will be the list of words [a1*b1,a2*b2,..,ar*br], with the (ai,bi) coming
   * in lex-order.
   */

    ht_ptrb = labelht.current_ptr;
    ht_ptre = ht_ptrb-1;
    ptr = cs_ptr;
    while (ptr <= cs_ptre) {
      csa = *(ptr++); csb = *(ptr++);
      if (((l1=label[csa]) != 0) && ((l2=label[csb]) != 0)) {
        labellist1 = genmultptr->states->labels->wordslist[l1];
        labellist2 = genmultptr->states->labels->wordslist[l2];
        j1=0;
        while (labellist1[j1]) {
          g1 = labellist1[j1][0];
          if (g1==0) g1=ngens1;
          j2=0;
          while (labellist2[j2]) {
            genmult2ptr->is_accepting[cstate]=TRUE;
            g2 = labellist2[j2][0];
            if (g2==0) g2=ngens1;
            if (ht_ptrb>ht_ptre || g1> *(ht_ptre-1) ||
                                 (g1== *(ht_ptre-1) &&  g2 > *ht_ptre) ){
    /* We have a new generator pair to be added to the end */
              *(++ht_ptre) = g1; *(++ht_ptre) = g2;
            }
            else {
              ht_chptr = ht_ptrb;
              while (*ht_chptr < g1)
                ht_chptr+=2;
              while (*ht_chptr==g1 && *(ht_chptr+1)<g2)
                  ht_chptr+=2;
              if (g1 < *ht_chptr || g2 < *(ht_chptr+1)) {
    /* we have a new generator pair to be added in the middle */
                ht_ptr = ht_ptre;
                ht_ptre += 2;
                while (ht_ptr >= ht_chptr) {
                  *(ht_ptr+2) = *ht_ptr;
                  ht_ptr--;
                }
                *ht_chptr = g1; *(ht_chptr+1) = g2;
              }
            }
            j2++;
          }
          j1++;
        }
      }
    }
/* That completes the calculation of the label for cstate */
    newlabel[cstate] = hash_locate(&labelht,ht_ptre-ht_ptrb+1);
  } 
  hash_clear(&ht);

  ct=0;
  for (i=1;i<=ns;i++) if (genmult2ptr->is_accepting[i])
    ct++;
  genmult2ptr->num_accepting = ct;
  if (ct==1 || ct != ns) {
    tmalloc(genmult2ptr->accepting,int,ct+1);
    ct = 0;
    for (i=1;i<=ns;i++) if (genmult2ptr->is_accepting[i])
      genmult2ptr->accepting[++ct] = i;
  }
  tfree(genmult2ptr->is_accepting);
  tfree(genmultptr->is_accepting);

/* Finally copy the records from the label hash-table into the set of labels */
  nlab =
  labelset->size = labelht.num_recs;
  if (kbm_print_level>=3)
    printf("    #There are %d distinct labels.\n",nlab);
  tmalloc(labelset->wordslist,gen **,nlab+1);
  for (i=1;i<=nlab;i++) {
    len = hash_rec_len(&labelht,i)/2;
    tmalloc(labelset->wordslist[i],gen *,len+1);
    ht_ptr = hash_rec(&labelht,i);
    for (j=0;j<len;j++) {
      tmalloc(labelset->wordslist[i][j],gen,3);
      labelset->wordslist[i][j][0] = ht_ptr[2*j];
      labelset->wordslist[i][j][1] = ht_ptr[2*j+1];
      labelset->wordslist[i][j][2] = 0;
    }
    labelset->wordslist[i][len] = 0;
  }

  hash_clear(&labelht);


  return genmult2ptr;
}
void write_fsa (fsa_to_print,gpname,postfix,fsaname)
	fsa *fsa_to_print;
	char *gpname;
	char *postfix;
	char *fsaname;
{
      char outf [100];
      strcpy(outf,gpname);
      strcat(outf,postfix);
      base_prefix(fsaname);
      strcat(fsaname,postfix);
      wfile = fopen(outf,"w");
      Printf("writing %s\n",outf);
      fsa_print(wfile,fsa_to_print,fsaname);
      fclose(wfile);
}


//based on  make_full_wd_fsa(wd_fsaptr,inv,start_no,rsptr)
int add_diagonals_to_wd_fsa(wd_fsaptr,inv,rsptr,all_diagonals,
		check_diagonals,start,end,limit, gpname,nobigger,diagnostics,min_wd_size,prevdiff2str)
	fsa *wd_fsaptr;
	int *inv;
        reduction_struct *rsptr;
	boolean all_diagonals;
	boolean check_diagonals;
	int start;
	int end;
	int limit;
	char *gpname;
	boolean nobigger;
	boolean diagnostics;
	int min_wd_size;
        char *prevdiff2str;
/* Close the set of word-differences under action g,$, then add all possible
 * new transitions.
 */
#define xxx 0
{ int i, j, gen1, l, ns, ns_save, n, g1, g2, size_pba, **table, ***wd_table;
  gen **wdn, *stw, *ptr, *ptre, *ptr2;
  static boolean hadwarning=FALSE;
  gen testword[4096];

Printf("start=%d,end=%d\n",start,end);

  size_pba = 1 + wd_fsaptr->alphabet->base->size;
  ns = wd_fsaptr->states->size;
  ns_save=ns; // it could change!
  wdn = wd_fsaptr->states->words;
  table = wd_fsaptr->table->table_data_ptr;
  wd_table = wd_fsaptr->table->table_data_dptr;
  char * cf;
  int * cfn;
  int * letters;
  fsa *diff2c=NULL;
  fsa *diff2=NULL;
  fsa *prevdiff2=NULL;
  int total = 0;
  int total_diagonals=0;
  char fsaname [100];
  if (prevdiff2str)
   if ((rfile = fopen(prevdiff2str,"r")) != 0) {
	Printf("reading %s\n",prevdiff2str);
  	tmalloc(prevdiff2,fsa,1);
  	fsa_read(rfile,prevdiff2,DENSE,0,0,TRUE,fsaname);
  	fclose(rfile);
   }
  if (check_diagonals || (limit==999999)) {
  tmalloc(cf,char,ns+1);
  for (i=1;i<=ns;i++)
        cf[i] = 0;
  tmalloc(letters,int,size_pba+1);
  for (i=1;i<=size_pba;i++)
        letters[i] = 0;
  tmalloc(cfn,int,ns+1);
  for (i=1;i<=ns;i++)
        cfn[i] = 0;
  char inf3 [100];
  strcpy(inf3,gpname);
  strcat(inf3,".diff2c");
  if ((rfile = fopen(inf3,"r")) != 0) {
	Printf("reading %s.diff2c\n",gpname);
  	tmalloc(diff2c,fsa,1);
  	fsa_read(rfile,diff2c,DENSE,0,0,TRUE,fsaname);
  	fclose(rfile);
   }
   if (end==999999) {
  	strcpy(inf3,gpname);
  	strcat(inf3,".diff2");

  	if ((rfile = fopen(inf3,"r")) != 0) {
		Printf("reading %s.diff2\n",gpname);
  		tmalloc(diff2,fsa,1);
  		fsa_read(rfile,diff2,DENSE,0,0,TRUE,fsaname);
  		fclose(rfile);
		wdn=diff2->states->words; //read words from diff2 - the biggy
		ns=diff2->states->size;
   	}
  }
	
  }
  int start_no=wd_fsaptr->states->size+1;
  int trace=0;
  int startstate=2;
  int endstate=ns;
  if (start>2 && start<=ns)
	startstate=start;
  if (end>0 && end<ns && end>=start)
	endstate=end;
  ns=ns_save;
  for (gen1=1;gen1<size_pba;gen1++)
  //for (gen1=start;gen1<10;gen1++)
  {
	if (limit>0 && (total_diagonals>limit)){
		Printf("%d diagonals exceeds limit\n",total_diagonals);
		break;
	}
  //printf("gen1=%d,ns=%d\n",gen1,ns);
	// multiply every wd by a generator on the left
	// to get all possible 'diagonal word differences'
	// where a 'square' exists
  if (diff2 && diff2c) {
	startstate = start_no;
	// the assumption is that diff2 (the biggy) starts with same states
	// as wd_fsaptr. 
  }
  for (i=startstate;i<=endstate;i++){
    int gen2=1;
    int jj=0;
	if (limit>0 && (total_diagonals>limit))
	   break;	
    if (diff2c) {
    	jj=diff_no(diff2c,wdn[i]);
        if (diff2 && diff2c) {
		if (jj > 0)
	 	 // add to diff2 diag
		{
			total_diagonals++;
      			n = (++wd_fsaptr->states->size);
      			if (n > wd_fsaptr->table->maxstates){
       				fprintf(stderr,"Too many word-differences. Increase maxwdiffs.\n");
       				return -1;
			}
      			tmalloc(wd_fsaptr->states->words[n],gen,genstrlen(wdn[i])+1);
      			genstrcpy(wd_fsaptr->states->words[n],wdn[i]);
      			for (j=1;j<=wd_fsaptr->alphabet->size;j++)
       				set_dense_target(table,j,n,0);
      		}
		continue; //next state of diff2
	}
    	if (jj==0) {
		cf[i]=-1;
		cfn[i]=-1;
		continue;
    	}
   }

    for (gen2=1;gen2<size_pba;gen2++) {
	if (limit>0 && (total_diagonals>limit))
	   break;	
	if (dense_dtarget(wd_table,gen1,gen2,i)){
    		testword[0]=inv[gen1];
    		genstrcpy(testword+1,wdn[i]);
    		if ((*reduce_word)(testword,rsptr)==-1)
        		return -1;
   		if (!all_diagonals && (testword[0]==wdn[i][0]))
   		//if (!all_diagonals && (testword[0]==inv[gen1]))
			continue; 
		if (nobigger)
	        	if (genstrlen(testword)>genstrlen(wdn[i]))
				continue;
		if (min_wd_size > 0)
	        	if (genstrlen(testword)<min_wd_size)
				continue;
		j=diff_no(wd_fsaptr,testword);
		if (prevdiff2 && j==0)
			j=diff_no(prevdiff2,testword);
		jj=0;
		if (diff2c) {
			jj=diff_no(diff2c,testword);
			if (jj==0){
			    if (limit != 999999)
				j=-1;
			    else
			  	j=99; //pretend weve already got  a spurious wd
			}
		}
	        if (check_diagonals && j>0 && j<=ns) {
			if (xxx==0 || (j>xxx && i<=xxx)) {
				if (cf[j]==0)
					letters[gen1]+=1;
				cf[j]+=1;
				cfn[j]=i;
				total++;
			}
		}
		else {
			
    		//if (diff_no(wd_fsaptr,testword)==0)bracket /* new state */
     			//printf("new wd formed from gen %d and wd %i\n",gen1,i); 
			if (j==0)
			//if (j==0 && (((++total_diagonals)%2)))
			//if (j==0 && (!((++total_diagonals)%2)))
			{
				total_diagonals++;
      				n = (++wd_fsaptr->states->size);
      				if (n > wd_fsaptr->table->maxstates){
        				fprintf(stderr,"Too many word-differences. Increase maxwdiffs.\n");
        				return -1;
				}
      				tmalloc(wd_fsaptr->states->words[n],gen,genstrlen(testword)+1);
      				genstrcpy(wd_fsaptr->states->words[n],testword);
				if (diagnostics)
				{
                  			display_eqn(i,testword,NULL,
						wd_fsaptr->alphabet->base->names);
				}
      				for (j=1;j<=wd_fsaptr->alphabet->size;j++)
         				set_dense_target(table,j,n,0);
      			}
		}
		break;
	}
    } //for (gen2=1;gen2<size_pba;gen2++)
  } //for (i=2;i<=ns;i++
  if (diff2 && diff2c)
    break; // only interested in producing a diff2 with no spurious wd's
  } //for (gen1=1;gen1<size_pba;gen1++)
Printf("Total diagonals=%d\n",total_diagonals);
if (check_diagonals) {
	printf("checking..\n");
        printf("letters\n");
	for (i=1;i<=size_pba;i++){
		printf("%d ",letters[i]);
	}
	printf("\n");
	for (i=1;i<=ns;i++){
		if (cfn[i]==-1)
			printf("%d *\n",i);
		else
			printf("%d,%d <-%d\n",i,cf[i],cfn[i]);
	}
	printf("checked..%d\n",total);
	exit (99);
}
/* Now fill out table */
  Printf("%d diagonals added to diff2\n",wd_fsaptr->states->size - ns);
  ns = wd_fsaptr->states->size;
  if (start_no<1)
    start_no = 1;
  for (i=start_no;i<=ns;i++){
    if ((i-start_no)>0&&(i-start_no)%100 == 0)
	Printf("%d ",i-start_no);
    for (g1=1;g1<=size_pba;g1++) for (g2=1;g2<=size_pba;g2++){
      if (g1==size_pba && g2==size_pba)
        continue; /* Don't want padding-symbol on both sides. */
      if (dense_dtarget(wd_table,g1,g2,i) == 0) {
        stw=wd_fsaptr->states->words[i];
        l=genstrlen(stw);
        if (g1==size_pba) {
            genstrcpy(testword,stw);
            testword[l]=g2; testword[l+1]=0;
        }
        else if (g2==size_pba) {
            testword[0]=inv[g1];
            genstrcpy(testword+1,stw);
            testword[l+1]=0;
        }
        else {
            testword[0]=inv[g1];
            genstrcpy(testword+1,stw);
            testword[l+1]=g2;
            testword[l+2]=0;
        }
        if ((*reduce_word)(testword,rsptr)==-1)
          return -1;
        if (n=diff_no(wd_fsaptr,testword))
          set_dense_dtarget(wd_table,g1,g2,i,n);
        if (n>0 && n<start_no)
          set_dense_dtarget(wd_table,inv[g1],inv[g2],n,i);
      }
    }
  }
  return 0;
}

int
make_full_wd_fsax(wd_fsaptr,inv,start_no,rsptr)
	fsa *wd_fsaptr;
	int *inv;
        int start_no;
        reduction_struct *rsptr;
/* Close the set of word-differences under inversion, and add all possible
 * transitions, starting at state number start_no.
 */
{ int i, j, l, ns, n, g1, g2, size_pba, **table, ***wd_table;
  gen **wdn, *stw, *ptr, *ptre, *ptri;
  static boolean hadwarning=FALSE;
  gen testword[4096];
/* Now fill out table */
  ns = wd_fsaptr->states->size;
size_pba = 1 + wd_fsaptr->alphabet->base->size;
  ns = wd_fsaptr->states->size;
  wdn = wd_fsaptr->states->words;
  table = wd_fsaptr->table->table_data_ptr;
  wd_table = wd_fsaptr->table->table_data_dptr;
  if (start_no<1)
    start_no = 1;
  for (i=start_no;i<=ns;i++){
    for (g1=1;g1<=size_pba;g1++) for (g2=1;g2<=size_pba;g2++){
      if (g1==size_pba && g2==size_pba)
        continue; /* Don't want padding-symbol on both sides. */
      if (dense_dtarget(wd_table,g1,g2,i) == 0) {
        stw=wd_fsaptr->states->words[i];
        l=genstrlen(stw);
        if (g1==size_pba) {
            genstrcpy(testword,stw);
            testword[l]=g2; testword[l+1]=0;
        }
        else if (g2==size_pba) {
            testword[0]=inv[g1];
            genstrcpy(testword+1,stw);
            testword[l+1]=0;
        }
        else {
            testword[0]=inv[g1];
            genstrcpy(testword+1,stw);
            testword[l+1]=g2;
            testword[l+2]=0;
        }
        if ((*reduce_word)(testword,rsptr)==-1)
          return -1;
        if (n=diff_no(wd_fsaptr,testword)) {
          set_dense_dtarget(wd_table,g1,g2,i,n);
        //if (n>0 && n<start_no)
          set_dense_dtarget(wd_table,inv[g1],inv[g2],n,i);
	}
      }
    }
  }
  return 0;
}

int diff_reducex(w,rs_wd)
        gen              *w;
        reduction_struct *rs_wd;
/* w is the word to be reduced using the word-difference machine  *wd_fsa.
 * It is assumed that wd_fsa->table->table_data_dptr is set up.
 * This function allocates its own space.
 * NOTE: No checks on the validity of the word are carried out.
 n*/
{ int ndiff, ngens, identity, padsymbol, wordlen, ***difftab,
      gct, *gpref, level, gen1, gen2, diff, diffct, newdiff,
      olen, nlen, i, j;
  boolean deqi, donesub, *cf;
  fsa *wd_fsa = rs_wd->wd_fsa;
  int maxv = MAXV;
  struct vertexd {
       gen genno;
       int diffno;
       int sublen;
       struct vertexd *backptr;
  }
       *gptr, *ngptr, *substruc;

/* vertexd is the structure used to store a vertex in the graph of strings
   for possible substitution. The components are as follows.
   backptr - points back to another vertexd, or to zero.
   genno  - the number of the generator at the end of the string.
   diffno - the word difference number of the string defined by following
            backptr back to zero (using genno), relative to the corresponding
            part of the word being reduced.
   sublen - plus or minus the length of this string. sublen is positive if and
            only if the string lexicographically precedes the
            corresponding part of the word being reduced.
   sublen is put in to save time.
    Another essential component of a vertexd is its level (i.e. the length of
    the string got by chasing back to the beginning of the word)
    but we always calculate this, using
    the integers defined by gpref. (See below))
*/
  if (wd_fsa->alphabet->type != PRODUCT || wd_fsa->alphabet->arity != 2) {
    fprintf(stderr,
        "Error: diff_reduce must be called with a word-difference machine.\n");
    return -1;
  }
  gen * wcopy;
  ndiff = wd_fsa->states->size;
  ngens = wd_fsa->alphabet->base->size;
  identity = wd_fsa->initial[1];
  padsymbol = ngens+1;
  wordlen= genstrlen(w);
  tmalloc(wcopy,gen,wordlen+10); 
  genstrcpy(wcopy,w);
  int wordlencopy=wordlen;
  if (wordlen<=0)
    return 0;

  difftab = wd_fsa->table->table_data_dptr;
  w = w-1; /* since code below assumes word is from w[1] .. w[wordlen]. */

  tmalloc(cf,boolean,ndiff+1);
/* cf is used as a characteristic function, when constructing a subset of the
  set  D  of word differences.
*/

  tmalloc(gpref,int,wordlen+1);
  gct= -1;
  gpref[0]= -1;
/* gpref[n]+1 is the number of vertices that have been defined after reading the
  first n elements of the word. These vertices are gptr[0],...,gptr[gpref[n]].
  We start by allocating space for maxv vertices.
*/

  tmalloc(gptr,struct vertexd,maxv);

/* Now we start reading the word. */
  level=0;
  int startg2=0;
  while (++level<=wordlen) {
    for (i=1;i<=ndiff;i++)
      cf[i]=FALSE;
/* Read the element of the word at position level. */
    gen1= w[level];
    startg2=1;
   if (gen1==0) 
	printf("bad gen1!\n");
    //printf("%d at level %d\n",gen1,level);

/* The next loop is over the identity and the subset of D defined at the
   previous level, level-1.
*/
    diff = identity;
    while (1) {
      deqi= diff==identity;
/* First look for a possible substitution of a shorter string */
      newdiff = dense_dtarget(difftab,gen1,padsymbol,diff);
      if (newdiff==identity) {
/* Make substitution and reduce length of word by 1. */
        i = level-1;
        if (!deqi) {
          substruc = gptr+diffct;
          do {
            w[i] = substruc->genno;
	      //if (w[i]==0) printf("bad2!!\n");
            substruc = substruc->backptr;
            i--;
          } while (substruc);
        }
        for (j=level;j<wordlen;j++)
          w[j] = w[j+1];
	      //if (w[j]==0) printf("bad4!!\n");
        w[wordlen] = 0;
        wordlen--;
        level= i>0 ? i-1 : i;
    //remove any 0's from the new word!!!
		//printf("remove any zeros after sub and reduce by 1!!\n");
		int jj,kk=0;
        	for (jj=level+1;jj<=wordlen;jj++){
			if (w[jj]==0) {
				kk++;
				wordlen--;
			}
			if (kk>0) {
				if (w[jj+kk]==0)
          				w[jj] = w[jj+kk+1];
					//printf("NO1!!\n");
				else
          				w[jj] = w[jj+kk];
			}
		}
		w[wordlen+1]=0;

/* Whenever we make a substitution, we have to go back one level more than
   expected, because of our policy of looking ahead for substitutions
   that reduce the length by 2.
*/
        gct = gpref[level];
        break;
      }
      else if (newdiff && level<wordlen) {
        j = dense_dtarget(difftab,w[level+1],padsymbol,newdiff);
        if (j==identity)
/* Make substitution and reduce length of word by 2. */
        { i = level-1;
          if (!deqi) {
            substruc=gptr+diffct;
            do {
              w[i] = substruc->genno;
	      //if (w[i]==0) printf("bad1!!\n");
              substruc = substruc->backptr;
              i--;
            } while (substruc);
          }
          for (j=level;j<wordlen-1;j++)
            w[j] = w[j+2];
	      //if (w[j]==0) printf("bad3!!\n");
          w[wordlen-1] = 0;
          wordlen -= 2;
    //remove any 0's from the new word!!!
          level = i>0 ? i-1 : i;
		//printf("remove any zeros after sub and red by 2!!\n");
		int jj,kk=0;
        	for (jj=level+1;jj<=wordlen;jj++){
			if (w[jj]==0) {
				kk++;
				wordlen--;
			}
			if (kk>0) {
				if (w[jj+kk]==0)
          				w[jj] = w[jj+kk+1];
					//printf("NO2!!\n");
				else
          				w[jj] = w[jj+kk];
			}
		}
		w[wordlen+1]=0;
          gct = gpref[level];
          break;
        }
        else  {
		if (!deqi) {
			startg2=0;
			//printf("shorter word %d -> %d\n",diff,newdiff);
		}
        }
      }

      donesub = FALSE;
/* Now we loop over the generator that is a candidate for substitution
   at this point.
*/
      for (gen2=startg2;gen2<=ngens;gen2++)
       if ((gen2==0)||(newdiff = dense_dtarget(difftab,gen1,gen2,diff))) {
	//printf("\n%d (%d,%d)-> %d",diff,gen1,gen2,newdiff); 
        if (newdiff==identity) {
          if (deqi) {
            if (gen2<gen1) {
              w[level] = gen2;
	      //printf("\n** %d = %d",level, gen2);
              level = level>1 ? level-2 : level-1;
              gct = gpref[level];
              donesub = TRUE;
              break;
            }
          }
          else if (gptr[diffct].sublen > 0) {
/* Make a substitution (by a string of equal length). */
	    int shorterby=0; 
            w[level] = gen2;
	      //printf("\n** %d = %d",level, gen2);
            i = level-1;
            substruc = gptr+diffct;
            do {
	      if (substruc->genno==0)
	      {
		//printf("shorter at %d!\n",i);
		shorterby++;
		//w[i]=1;
	      }
              w[i] = substruc->genno;
	      //printf("\n* %d = %d",i,  w[i]);
              substruc = substruc->backptr;
              i--;
            } while (substruc);
            level = i>0 ? i-1 : i;
            gct = gpref[level];
	    if (shorterby) {
		//printf("shorterby!!\n");
		int j,k=0;
        	for (j=level+1;j<=wordlen;j++){
			if (w[j]==0) {
				k++;
				wordlen--;
			}
			if (k>0) {
				if (w[j+k]==0) {
					w[j]=w[j+k+1];
				}
				else
          				w[j] = w[j+k];
			}
		}
        	w[wordlen+1] = 0;
	    }
            donesub = TRUE;
            break;
          }
        }
        else {
          if (gen2>0 && cf[newdiff])
/* We have this word difference stored already, but we will check to see if
   the current string precedes the existing one.
*/
            for (i=gpref[level-1]+1;;i++) {
              substruc=gptr+i;
              if (substruc->diffno == newdiff) {
                olen = substruc->sublen;
                nlen = deqi ? (gen2 < gen1 ? 1 : -1) : 
                       (j = (gptr[diffct].sublen))>0 ? j+1 : j-1;
                if (nlen > olen) {
/* The new string is better than the existing one */
                  substruc->genno = gen2;
                  substruc->sublen = nlen;
                  substruc->backptr = deqi ? 0 : gptr+diffct;
                }
                break;
              }
            }
          else
/* This is a new word difference at this level, so we define a new vertexd in
   graph.
*/
          { gct++;
            if (gct >= maxv) {
/* We need more space for vertices. Allocate twice the preceding space and
   copy existing data.
*/
              tmalloc(ngptr,struct vertexd,2*maxv);
              if (kbm_print_level>=3)
                printf("    #Allocating more space in diff_reduce.\n");
              for (i=0;i<maxv;i++) {
                ngptr[i].genno = gptr[i].genno;
                ngptr[i].diffno = gptr[i].diffno;
                ngptr[i].sublen = gptr[i].sublen;
                substruc = gptr[i].backptr;
                if (substruc==0)
                  ngptr[i].backptr = 0;
                else
                  for (j=i-1;;j--) if (substruc==gptr+j) {
                    ngptr[i].backptr = ngptr+j;
                    break;
                  }
              }
              tfree(gptr);
              gptr=ngptr;
              maxv *= 2;
            }
/* Define the new vertexd. */
            substruc = gptr+gct;
            nlen = deqi ? (gen2<gen1 ? 1 : -1) : 
                       (j = (gptr[diffct].sublen))>0 ? j+1 : j-1;
            substruc->genno = gen2;
            substruc->diffno = newdiff;
            substruc->sublen = nlen;
            substruc->backptr = deqi ? 0 : gptr+diffct;
	    if (gen2!=0)
               cf[newdiff] = TRUE;
	    if (gen2==0) {
		if ((gptr+diffct)->genno==0)
			// dont want two $s in succession - forget about it!!
			gct--;	
		else {
                 cf[newdiff] = TRUE;
		 if (nlen<0)
			substruc->sublen=nlen*-1;
		}
		//gct--;
		//printf("shorter word %d -> %d\n",diff,newdiff);
		startg2=1;
	    }
          }
        }
      } /*End of loop over gen2 */

      if (donesub)
        break;

/* Go on to next word difference from previous level. */
      if (diff==identity) {
        if (level==1) break;
        diffct = gpref[level-2]+1;
      }
      else
        diffct++;
      if (diffct > gpref[level-1])
         break;
      diff = gptr[diffct].diffno;
    } /* end of loop over word differences at previous level */

    gpref[level] = gct;
  }
  tfree(gptr);
  tfree(cf);
  tfree(gpref);
  tfree(wcopy);
  return 0;
}

fsa * fsa_ex_min (fsaptr,  op_table_type,tempfilename,minstr,max_state,max_level)
	fsa *fsaptr;
	storage_type op_table_type;
	char *tempfilename;
	char *minstr;
	int max_state;
	int max_level;
{ 

// extract wds using minred
// create list of lhs's suitable for scanning using gpcheckx -t option
// invoke by -extractwdsfrommin 
int  ***dtable, ne, ngens, ndiff, ns, *fsarow,  nt, cstate, cs, csdiff, csi,
       im, i, k, g1, g2, len, identity;
 unsigned short int *ht_ptr, *ht_ptr_save, *ht_ptrb, *ht_ptre, *cs_ptr, *cs_ptre, *ptr;
  boolean dense_op,  no_trans_by_wa,good;
  char *cf;
  int *wa1states;
  int *historys;
  int *historyl;
  short_hash_table ht;
  fsa *wa;
  fsa *gpmin;
  FILE *tempfile;
  int **watable;
  int **mintable;
  char fsaname [100];
  FILE  *fopen();
  unsigned int	total_elements=0;
  unsigned int tot_space_save;
  unsigned int block_space_save;
  unsigned int fail_state;
unsigned int no_fails=0;
unsigned int no_accepts=0;
unsigned int no_prev_accepts=0;
boolean watype=FALSE;
Printf("Extracting new word differences from minred?? and diff2\n");	
    if ((rfile = fopen(minstr,"r")) == 0) {
        fprintf(stderr,"Cannot open file %s.\n",minstr);
        exit(1);
    }
    Printf("reading %s\n",minstr);
    tmalloc(gpmin,fsa,1);
    fsa_read(rfile,gpmin,DENSE,0,0,TRUE,fsaname);
    Printf("%s has size %d\n",minstr,gpmin->states->size);
    fclose(rfile);               
  watable=gpmin->table->table_data_ptr;
  if (gpmin->num_accepting==1)
	fail_state=gpmin->accepting[1];
  else {
	fail_state=0;
 	watype=TRUE; 
  }
  if (!fsaptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_wa only applies to DFA's.\n");
    return 0;
  }

  if (fsaptr->alphabet->type != PRODUCT || fsaptr->alphabet->arity != 2) {
    fprintf(stderr,"Error in fsa_wa: fsa must be 2-variable.\n");
    return 0;
  }
  if (fsaptr->states->type != WORDS) {
    fprintf(stderr,"Error in fsa_wa: fsa must be word-difference type.\n");
    return 0;
  }

  tmalloc(wa,fsa,1);
  fsa_init(wa);
  srec_copy(wa->alphabet,fsaptr->alphabet->base);
  wa->flags[DFA] = TRUE;
  wa->flags[ACCESSIBLE] = TRUE;
  wa->flags[BFS] = TRUE;

  ne = fsaptr->alphabet->size;
  ngens = wa->alphabet->size;
  ndiff = fsaptr->states->size;

  if (ne != (ngens+1)*(ngens+1)-1) {
   fprintf(stderr,
       "Error: in a 2-variable fsa, alphabet size should = ngens^2 - 1.\n");
    return 0;
  }

  identity = fsaptr->accepting[1]; /* assumed to be unique */
  if (fsaptr->num_accepting!=1 || (identity != fsaptr->initial[1])) {
    fprintf(stderr,"Error: Input to fsa_wa not a word-difference machine.\n");
    return 0;
  }

  if (fsaptr->table->table_type!=DENSE) {
     fprintf(stderr,
      "Error: function fsa_wa can only be called with a densely-stored fsa.\n");
     return 0;
  }
  dense_op = op_table_type==DENSE;
  if (fsa_table_dptr_init(fsaptr)== -1) return 0;
  dtable = fsaptr->table->table_data_dptr;

  wa->states->type = SIMPLE;
  wa->num_initial = 1;
  tmalloc(wa->initial,int,2);
  wa->initial[1] = 1;
  wa->table->table_type = op_table_type;
  wa->table->denserows = 0;
  wa->table->printing_format = op_table_type;
  
  short_hash_init(&ht,FALSE,0,0,0);
  ht_ptr = ht.current_ptr;
  *(unsigned int *) ht_ptr=identity;
  //ht_ptr[0] = identity; // wa id
  *(ht_ptr+2)=identity;
  //ht_ptr[1] = identity;
  im = short_hash_locate(&ht,3);
  if (im== -1) return 0;
  if (im!=1) {
    fprintf(stderr,"Hash-initialisation problem in fsa_wa.\n");
     return 0;
  }

  if ((tempfile=fopen(tempfilename,"w"))==0){
    fprintf(stderr,"Error: cannot open file %s\n",tempfilename);
     return 0;
  }
  tmalloc(fsarow,int,ngens);
  len=ngens;
 
  cstate = 0;
  nt = 0; /* Number of transitions in exists */
  tmalloc(cf,char,ndiff+1);
  int incx=0; 
  int incx2=0; 
  int kk=0;
  int no_statesx=1;
  int total_hashx=sizeof(int);
  int total_rhs=0;
  int total_equal=0;
  int currentlevel=0;
  int maxthislevel=1;
  //boolean first_coincidence=TRUE;
  tmalloc(wa->accepting,int,1);
  wa->num_accepting = 1;
  wa->accepting[1] = 1000000000;
  while (++cstate <= ht.num_recs) {
    if (max_state>0)
     if (cstate > max_state)
	break;
    if (max_level>0)
     if (currentlevel > max_level)
	break;
    if (cstate>maxthislevel) {
	maxthislevel=ht.num_recs;
	currentlevel++;
	if (no_accepts > 0){
		Printf("level %d at state %d,accepts %d\n",currentlevel,cstate,no_accepts);
		no_prev_accepts = no_accepts;
	}
    }
    if (kbm_print_level>1 && no_prev_accepts==no_accepts) {
      if (
          //(cstate<=1000000 && cstate%5000==0) || cstate%50000==0)
          cstate%50000==0)
       Printf("    #cstate = %d;  number of states = %d, level = %d\n",cstate,ht.num_recs, currentlevel);
    }
  ht_ptr_save=ht.current_ptr;
  block_space_save = ht.block_space;
  tot_space_save=ht.tot_space;
    cs_ptr = short_hash_rec(&ht,cstate);
    cs_ptre = short_hash_rec(&ht,cstate) + short_hash_rec_len(&ht,cstate) - 1;
    if (!dense_op)
      len = 0;
    for (g1=1;g1<=ngens;g1++) {
/* Calculate action of generator g1 on state cstate  - to get the image, we
 * have to apply (g1,g2) to each element in the subset corresponding to cstate,
 * and this for each generator g2 of the base-alphabet (including the padding
 * symbol).
 */
      int wa1state=watable[g1][*(unsigned int *)cs_ptr];
      no_trans_by_wa=FALSE;
      boolean looking_for_1s=FALSE;
      boolean found_1s=FALSE;
      boolean accept_state = FALSE;
      if (wa1state==fail_state) { 
	looking_for_1s=TRUE;
	//no_fails++;
	no_trans_by_wa=TRUE;
/* check if any of the wd's of this state have a transition by g1/x to 1 */
/* If not, then this must be a lhs which doesnt fellow travel yet => new wd(s) */
      }
      if (wa1state==0)
      {
	no_trans_by_wa=TRUE;
      }
      if (!looking_for_1s) {
      for (i=1;i<=ndiff;i++) {
        cf[i] = 0;
      }
      }
      ptr = cs_ptr+2; // first wd
   /* csdiff is the state of *fsaptr corresponding to cs */
      while (ptr <= cs_ptre) {
	int gt_state, old_gt_state;;
        csdiff =  *ptr; // dont add initial state
        //csdiff = cs%ndiff;
        //if (csdiff==0) csdiff = ndiff;

    	// following code only applicable for .wa - where is always satisfied 
	// if we allow 1-g1/g1>1 which we must!
   /* .wa code */
	if (csdiff<5000&&watype && looking_for_1s) {
		int target=g1; 
		if (g1%2)
			target+=2;
        	if (csdiff==target) 
	   		found_1s=TRUE; 
	        else for (g2=1;g2<ngens+1;g2++){
			 if (dense_dtarget(dtable,g1,g2,csdiff)==identity){
				found_1s=TRUE;
				break;
			}
		}
       }
       if (csdiff>5000&&watype && looking_for_1s) {
	   		found_1s=TRUE; 
	}
   /*.wa code */

   /* csdiff is the state of *fsaptr corresponding to cs */
        ptr++;
	old_gt_state=1;
        if (csdiff == identity || csdiff == identity+5000) { 
	        for (g2=1;g2<ngens+1;g2++){
			if (!watype && g1==g2)
		        	continue;
                         int csdiffx=csdiff;
			 if (csdiff>5000)
			{
				csdiffx=csdiff-5000;
                        }
 			csi =  dense_dtarget(dtable,g1,g2,csdiffx);
			if (csi==0)
       		     		continue;
			if (watype && (csi==identity||csdiff>5000))
              			cf[csi] = 2;
			else
				cf[csi] = 1;
		}
	} 
	else { //not identity
		for (g2=1;g2<ngens+1;g2++){
                         int csdiffx=csdiff;
			 if (csdiff>5000)
			{
				csdiffx=csdiff-5000;
                        }
		  csi =  dense_dtarget(dtable,g1,g2,csdiffx);
		  if (csi==0)
		    continue;
		  gt_state=old_gt_state;	
		  if (looking_for_1s && (csdiff < 5000)) {
		  	if (csi==identity) {
				found_1s = TRUE;
				break;
		  	}
		  	continue;
		  }
		  if (looking_for_1s && (csdiff > 5000)) {
			found_1s = TRUE;
			break;
		  }
		  if (csdiff>5000) 
              		cf[csi] = 2;
		  else
			cf[csi] = 1;
		} // for g2

	  if (looking_for_1s && (csdiff < 5000)) {
		/* now deal with pad on the right */
	  	csi =  dense_dtarget(dtable,g1,g2,csdiff);
          	if (csi==identity)
	  	{
			found_1s = TRUE;
			break;
	  	}
	  }
	  if (looking_for_1s && (csdiff > 5000)) {
		found_1s = TRUE;
		break;
	  }
	} // not identity
      } // for ptr
      if (looking_for_1s && !found_1s&& csdiff<5000) {
		accept_state=TRUE;
      }
      else
      {
	accept_state=FALSE;
      }
      if (no_trans_by_wa) {
	if (accept_state) {
		fsarow[g1-1] = 1000000000; // make the last state + 1 the accept state
		//if (no_accepts<20)
		//	Printf("%d-%d>accept\n",cstate,g1);
		no_accepts++;
	}
	else {
		fsarow[g1-1] = 0;
		no_fails++;
	}
        continue;
      }
      
      ht_ptrb = ht.current_ptr;
      ht_ptre = ht_ptrb-1;
      *(unsigned int *)(++ht_ptre)=wa1state;
      ht_ptre++;
      int local_hashx=0;
      int local_rhs=0;
      int local_equal=0;
      for (i=1;i<=ndiff;i++) {
        k = cf[i];
        if (k>0) {
		local_hashx++;
	}
	if (k==1)
         *(++ht_ptre) = i;
	if (k==2)
         *(++ht_ptre) = i+5000;
      }
      im = short_hash_locate(&ht,ht_ptre-ht_ptrb+1);
      if (im== -1) return 0;
      fsarow[g1-1] = im; 
	
      if (max_state > 0)
     	 if (im>max_state) {
      		fsarow[g1-1] = 0; 
	}
    } // for gens
    fwrite((void *)fsarow,sizeof(int),(size_t)len,tempfile);
  } // while states
  // write an accept state
  int lenx=len; 
  while (lenx>0) {
	fsarow[--lenx]=0;
  }
  fwrite((void *)fsarow,sizeof(int),(size_t)len,tempfile);
  
  fclose(tempfile);
  Printf("last state=%d\n",cstate);
  Printf(" Accepts=%d\n",no_accepts);
  short_hash_clear(&ht);
  tfree(cf);
  ns = wa->states->size = ht.num_recs+1;
  wa->accepting[1] = ns;
  wa->table->numTransitions = nt;
/* Now read the transition table back in */
  tempfile = fopen(tempfilename,"r");
  if (no_accepts>0)
  //compressed_transitions_read(wa,tempfile);
  {
  int ns, ne, nt, **table, *tab_ptr, len, i, j;

  ns = wa->states->size;
  ne = wa->alphabet->size;
  nt = wa->table->numTransitions;

  if (wa->table->table_type == DENSE) {
    /* it is a pity that the table was output in transposed form! */
    fsa_table_init(wa->table, ns, ne);
    table = wa->table->table_data_ptr;
    for (i = 1; i <= ns; i++)
      for (j = 1; j <= ne; j++) {
        size_t s = fread((void *)(table[j] + i), sizeof(int), (size_t)1, tempfile);
	if (*(table[j]+i) > ns) 
		*(table[j]+i)=ns;
        (void)s; // HACK to silence compiler warning
      }
  }
  } // compressed_transitions
  fclose(tempfile);

  unlink(tempfilename);
  if (no_accepts>0)
  	return wa; // scan using -t switch
  return NULL;

}

fsa * fsa_extractwdsfromwa (fsaptr, op_table_type,tempfilename,wastr,max_state,max_level)
	fsa *fsaptr;
	storage_type op_table_type;
	char *tempfilename;
	char * wastr; 
	int max_state;
	int max_level;
{ int  ***dtable, ne, ngens, ndiff, ns, *fsarow, nt, cstate, cs, csdiff, csi,
       im, i, k, g1, g2, len, identity;
  unsigned short int *ht_ptr, *ht_ptrb, *ht_ptre, *cs_ptr, *cs_ptre, *ptr;
  boolean dense_op, no_trans, no_trans_by_wa, good;
  char *cf;
  short_hash_table ht;
  fsa *wa, *gpwa;
  FILE *tempfile, *fopen();
  char fsaname [100];
  int **watable;
  int fail_state=0;

  int SEEN_LHS_BETTER = 1;
  int SEEN_RHS_BETTER = 2;
  int SEEN_EQUAL = 3;

  unsigned int no_accepts=0;
  unsigned int no_prev_accepts=0;

    Printf("Extracting new word differences from wa and diff2 with diagonals added\n");
    if ((rfile = fopen(wastr,"r")) == 0) {
        fprintf(stderr,"Cannot open file %s.\n",wastr);
        exit(1);
    }
    Printf("reading %s\n",wastr);
    tmalloc(gpwa,fsa,1);
    fsa_read(rfile,gpwa,DENSE,0,0,TRUE,fsaname);
    Printf("%s has size %d\n",wastr,gpwa->states->size);
    fclose(rfile);               
  watable=gpwa->table->table_data_ptr;
  if (!fsaptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_wa only applies to DFA's.\n");
    return 0;
  }

  if (fsaptr->alphabet->type != PRODUCT || fsaptr->alphabet->arity != 2) {
    fprintf(stderr,"Error in fsa_wa: fsa must be 2-variable.\n");
    return 0;
  }
  if (fsaptr->states->type != WORDS) {
    fprintf(stderr,"Error in fsa_wa: fsa must be word-difference type.\n");
    return 0;
  }

  tmalloc(wa,fsa,1);
  fsa_init(wa);
  srec_copy(wa->alphabet,fsaptr->alphabet->base);
  wa->flags[DFA] = TRUE;
  wa->flags[ACCESSIBLE] = TRUE;
  wa->flags[BFS] = TRUE;

  ne = fsaptr->alphabet->size;
  ngens = wa->alphabet->size;
  ndiff = fsaptr->states->size;

  if (ne != (ngens+1)*(ngens+1)-1) {
   fprintf(stderr,
       "Error: in a 2-variable fsa, alphabet size should = ngens^2 - 1.\n");
    return 0;
  }

  identity = fsaptr->accepting[1]; /* assumed to be unique */
  if (fsaptr->num_accepting!=1 || (identity != fsaptr->initial[1])) {
    fprintf(stderr,"Error: Input to fsa_wa not a word-difference machine.\n");
    return 0;
  }

  if (fsaptr->table->table_type!=DENSE) {
     fprintf(stderr,
      "Error: function fsa_wa can only be called with a densely-stored fsa.\n");
     return 0;
  }
  dense_op = op_table_type==DENSE;
  if (fsa_table_dptr_init(fsaptr)== -1) return 0;
  dtable = fsaptr->table->table_data_dptr;

  wa->states->type = SIMPLE;
  wa->num_initial = 1;
  tmalloc(wa->initial,int,2);
  wa->initial[1] = 1;
  wa->table->table_type = op_table_type;
  wa->table->denserows = 0;
  wa->table->printing_format = op_table_type;
  
  short_hash_init(&ht,FALSE,0,0,0);
  ht_ptr = ht.current_ptr;
  *(unsigned int *) ht_ptr=identity;
  //*(unsigned int *) (ht_ptr+2)=identity; //wa2state
  *(ht_ptr+2)=identity;
  //*(ht_ptr+4)=identity; //wa2state
  im = short_hash_locate(&ht,3);
  //im = short_hash_locate(&ht,5);//wa2state
  //ht_ptr[0] = identity;
  //im = short_hash_locate(&ht,1);
  if (im== -1) return 0;
  if (im!=1) {
    fprintf(stderr,"Hash-initialisation problem in fsa_wa.\n");
     return 0;
  }
  if ((tempfile=fopen(tempfilename,"w"))==0){
    fprintf(stderr,"Error: cannot open file %s\n",tempfilename);
     return 0;
  }
  if (dense_op)
    tmalloc(fsarow,int,ngens)
  else
    tmalloc(fsarow,int,2*ngens+1)
 
  cstate = 0;
  if (dense_op)
    len = ngens; /* The length of the fsarow output. */
  nt = 0; /* Number of transitions in exists */
  tmalloc(cf,char,ndiff+1);
  int dollarcount=0;
  int total_hashx=sizeof(short int);
  tmalloc(wa->accepting,int,1);
  wa->num_accepting = 1;
  wa->accepting[1] = 1000000000;
  int currentlevel=0;
  int maxthislevel=1;
  while (++cstate <= ht.num_recs) {
    if (max_state>0)
     if (cstate > max_state)
        break;
    if (max_level>0)
     if (currentlevel > max_level)
	break;
    if (cstate>maxthislevel) {
	maxthislevel=ht.num_recs;
	currentlevel++;
	if (no_accepts > 0){
		Printf("level %d at state %d,accepts %d\n",currentlevel,cstate,no_accepts);
		no_prev_accepts = no_accepts;
	}
    }
    if (kbm_print_level>1 && no_prev_accepts==no_accepts) {
      if (
          //(cstate<=1000000 && cstate%5000==0) || cstate%50000==0)
          cstate%50000==0)
       Printf("    #cstate = %d;  number of states = %d, level %d\n",cstate,ht.num_recs,currentlevel);
    }
    cs_ptr = short_hash_rec(&ht,cstate);
    cs_ptre = short_hash_rec(&ht,cstate) + short_hash_rec_len(&ht,cstate) - 1;
    if (!dense_op)
      len = 0;
    for (g1=1;g1<=ngens;g1++) {
/* Calculate action of generator g1 on state cstate  - to get the image, we
 * have to apply (g1,g2) to each element in the subset corresponding to cstate,
 * and this for each generator g2 of the base-alphabet (including the padding
 * symbol).
 * Since we are excluding words that contain subwords w_1 s.t. (w_1,w_2) is
 * accepted by *fsaptr, we also have to apply (g1,g2) to the initial state
 * of *fsaptr.
 */
      int local_hashx=0;
      int wa1state=watable[g1][*(unsigned int *)cs_ptr];
      //int wa2state=watable[g1][*(unsigned int *)cs_ptr+2];//wa2state
      //if (cstate==1) //wa2state
      //		wa2state=1; //wa2state
      // wa2state ensures that only minimally reducible words are accepted
      no_trans_by_wa=FALSE;
      boolean accept_state = FALSE;

      for (i=1;i<=ndiff;i++)
        cf[i] = 0;
      ptr = cs_ptr-1;
      ptr+=2; // skip the wa1state
      //ptr+=2; // skip the wa2state
   /* csdiff is the state of *fsaptr corresponding to cs */
      no_trans = FALSE;
/* We will set no_trans to be true if we find that the transition leads to
 * failure.
 */
      while (ptr <= cs_ptre) {
	int gt_state, old_gt_state;;
/* We add the initial state of *fsaptr to the subset representing cstate */
        cs = ptr<(cs_ptr+2) ? identity : *ptr;
        //cs = ptr<(cs_ptr+4) ? identity : *ptr; //wa2state
        csdiff = cs%ndiff;
        if (csdiff==0) csdiff = ndiff;
   /* csdiff is the state of *fsaptr corresponding to cs */
        ptr++;
	// MAF algorithm inserted here
	old_gt_state=SEEN_LHS_BETTER;
	if (cs>ndiff)
		old_gt_state=SEEN_RHS_BETTER;
	if (cs>2*ndiff)
		old_gt_state=SEEN_EQUAL;
        if (csdiff == identity) { 
		old_gt_state=SEEN_EQUAL;
	        for (g2=1;g2<ngens+1;g2++){
 			csi =  dense_dtarget(dtable,g1,g2,csdiff);
			if (csi==0)
       		     		continue;
 			gt_state = g2 < g1 ? SEEN_RHS_BETTER : SEEN_LHS_BETTER;

          		if (csi==identity)
			{
			 if (gt_state == SEEN_RHS_BETTER)
			 {
            			no_trans = TRUE;
            			break;
			 }
			 continue;
			}
			if (g1 == g2)
				gt_state = SEEN_EQUAL;
            		if ( gt_state > cf[csi])
              			cf[csi] = gt_state;
		}
	}
	else {
		for (g2=1;g2<ngens+1;g2++){
		  csi =  dense_dtarget(dtable,g1,g2,csdiff);
		  if (csi==0)
		    continue;
		  gt_state=old_gt_state;	
		  if (csi==identity) {
           		if (gt_state == SEEN_RHS_BETTER)
			{
			      no_trans = TRUE;
			      break;
		    	}
			continue;
		  }
            	  if ( gt_state > cf[csi])
              		cf[csi] = gt_state;
		}
	}
        if (no_trans)
          break;
	/* now deal with pad on the right */
	  csi =  dense_dtarget(dtable,g1,g2,csdiff);
          if (csi==identity)
	  {
            	no_trans = TRUE;
            	break;
	  }
	  else {if (csi)
		    cf[csi] = SEEN_RHS_BETTER;
	       else if (WADIAG)
		    //printf("%d-%d,%d>?\n",csdiff,g1,g2);
		    dollarcount++;
	  }
      } // for ptr
      if (no_trans) {
        if (dense_op)
          fsarow[g1-1] = 0;
        //if (wa1state > 0 && wa2state > 0)  //wa2state
        if (wa1state > 0 ) {
               fsarow[g1-1] = 1000000000; // make the last state + 1 the accept state
                no_accepts++;	
	}
        continue;
      }
/* Now we have the image stored in the array cf, and we translate it to a list
 * and insert it into the hash-table.
 */
	if (TRACE5 && cstate< 80) printf("%d****\n",cstate);
      ht_ptrb = ht.current_ptr;
      ht_ptre = ht_ptrb-1;
      *(unsigned int *)(++ht_ptre)=wa1state;
      ht_ptre++;
      //*(unsigned int *)(++ht_ptre)=wa2state; //wa2state
      //ht_ptre++;			     //wa2state
      for (i=1;i<=ndiff;i++) {
        k = cf[i];
        if (k>0) {
                local_hashx++;
        }
	if (TRACE5 && cstate < 80 && cf[i])
		printf(" %d-%d\n",i,cf[i]);
        if (k==SEEN_LHS_BETTER)
          *(++ht_ptre) = i;
        else if (k==SEEN_RHS_BETTER)
          *(++ht_ptre) = ndiff+i;
        else if (k==SEEN_EQUAL)
          *(++ht_ptre) = ndiff+ndiff+i;
      }
      im = short_hash_locate(&ht,ht_ptre-ht_ptrb+1);
      if (im== -1) return 0;
	if (TRACE5 && cstate<80) printf("%d->%d\n",g1,im);
      if (dense_op)
         fsarow[g1-1] = im;
      else if (im>0) {
         fsarow[++len] = g1;
         fsarow[++len] = im;
      }
      if (im>0) 
        nt++;
        
      if (max_state > 0)
     	 if (im>max_state) {
      		fsarow[g1-1] = 0; 
	}
    } //for g1 
    if (!dense_op)
      fsarow[0] = len++;
    fwrite((void *)fsarow,sizeof(int),(size_t)len,tempfile);
  } // while states
  // write an accept state
  int lenx=len;
  while (lenx>0) {
        fsarow[--lenx]=0;
  }
  fwrite((void *)fsarow,sizeof(int),(size_t)len,tempfile);
  fclose(tempfile);

  Printf("last state=%d\n",cstate);
  Printf("Accepts=%d\n",no_accepts);
  short_hash_clear(&ht);
  tfree(fsarow);
  tfree(cf);
  if (WADIAG)
	printf("%d dollar calculations\n",dollarcount);

  wa->table->numTransitions = nt;
  ns = wa->states->size = ht.num_recs+1;
  wa->accepting[1] = ns;
  wa->table->numTransitions = nt;

/* Now read the transition table back in */
  tempfile = fopen(tempfilename,"r");
  if (no_accepts>0)
  //compressed_transitions_read(wa,tempfile);
  {
  int ns, ne, nt, **table, *tab_ptr, len, i, j;

  ns = wa->states->size;
  ne = wa->alphabet->size;
  nt = wa->table->numTransitions;

  if (wa->table->table_type == DENSE) {
    /* it is a pity that the table was output in transposed form! */
    fsa_table_init(wa->table, ns, ne);
    table = wa->table->table_data_ptr;
    for (i = 1; i <= ns; i++)
      for (j = 1; j <= ne; j++) {
        size_t s = fread((void *)(table[j] + i), sizeof(int), (size_t)1, tempfile);
        if (*(table[j]+i) > ns)
                *(table[j]+i)=ns;
        (void)s; // HACK to silence compiler warning
      }
  }
  } // compressed_transitions
  fclose(tempfile);

  unlink(tempfilename);
  if (no_accepts>0)
  	return wa;
  return NULL;
}
fsa * fsa_extractwdsfrom_recovery (fsaptr, op_table_type,tempfilename,recovery_num)
	fsa *fsaptr;
	storage_type op_table_type;
	char *tempfilename;
        int recovery_num;
{ int  ***dtable, ne, ngens, ndiff, ns, *fsarow, nt, cstate, cs, csdiff, csi,
       im, i, k, g1, g2, len, identity;
  unsigned short int *ht_ptr, *ht_ptrb, *ht_ptre, *cs_ptr, *cs_ptre, *ptr;
  boolean dense_op, no_trans, no_trans_by_wa, good;
  char *cf;
  short_hash_table ht;
  fsa *wa, *gpwa;
  FILE *tempfile, *fopen();
  char fsaname [100];
  int **watable;
  int fail_state=0;

  int SEEN_LHS_BETTER = 1;
  int SEEN_RHS_BETTER = 2;
  int SEEN_EQUAL = 3;

  unsigned int no_accepts=0;
  unsigned int no_prev_accepts=0;

    Printf("Recovery!!!!\n");

  if (!fsaptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_wa only applies to DFA's.\n");
    return 0;
  }

  if (fsaptr->alphabet->type != PRODUCT || fsaptr->alphabet->arity != 2) {
    fprintf(stderr,"Error in fsa_wa: fsa must be 2-variable.\n");
    return 0;
  }
  if (fsaptr->states->type != WORDS) {
    fprintf(stderr,"Error in fsa_wa: fsa must be word-difference type.\n");
    return 0;
  }

  tmalloc(wa,fsa,1);
  fsa_init(wa);
  srec_copy(wa->alphabet,fsaptr->alphabet->base);
  wa->flags[DFA] = TRUE;
  wa->flags[ACCESSIBLE] = TRUE;
  wa->flags[BFS] = TRUE;

  ne = fsaptr->alphabet->size;
  ngens = wa->alphabet->size;
  ndiff = fsaptr->states->size;

  if (ne != (ngens+1)*(ngens+1)-1) {
   fprintf(stderr,
       "Error: in a 2-variable fsa, alphabet size should = ngens^2 - 1.\n");
    return 0;
  }

  identity = fsaptr->accepting[1]; /* assumed to be unique */
  if (fsaptr->num_accepting!=1 || (identity != fsaptr->initial[1])) {
    fprintf(stderr,"Error: Input to fsa_wa not a word-difference machine.\n");
    return 0;
  }

  if (fsaptr->table->table_type!=DENSE) {
     fprintf(stderr,
      "Error: function fsa_wa can only be called with a densely-stored fsa.\n");
     return 0;
  }
  dense_op = op_table_type==DENSE;
  if (fsa_table_dptr_init(fsaptr)== -1) return 0;
  dtable = fsaptr->table->table_data_dptr;

  wa->states->type = SIMPLE;
  wa->num_initial = 1;
  tmalloc(wa->initial,int,2);
  wa->initial[1] = 1;
  wa->table->table_type = op_table_type;
  wa->table->denserows = 0;
  wa->table->printing_format = op_table_type;
  
  if ((tempfile=fopen(tempfilename,"r"))==0){
    fprintf(stderr,"Error: cannot open file %s\n",tempfilename);
     return 0;
  }
  if (dense_op)
    tmalloc(fsarow,int,ngens)
  else
    tmalloc(fsarow,int,2*ngens+1)
 
  cstate = 0;
  if (dense_op)
    len = ngens; /* The length of the fsarow output. */
  nt = 0; /* Number of transitions in exists */
  tmalloc(cf,char,ndiff+1);
  int dollarcount=0;
  int total_hashx=sizeof(short int);
  tmalloc(wa->accepting,int,1);
  wa->num_accepting = 1;
  wa->accepting[1] = 1000000000;
  int currentlevel=0;
  int maxthislevel=1;
  recovery_num++;

  tfree(fsarow);
  tfree(cf);
  if (WADIAG)
	printf("%d dollar calculations\n",dollarcount);

  wa->table->numTransitions = 9999;
  wa->accepting[1] = recovery_num;

/* Now read the transition table back in */
  no_accepts=1;
  if (no_accepts>0)
  //compressed_transitions_read(wa,tempfile);
  {
  int ns, ne, nt, **table, *tab_ptr, len, i, j;

  wa->states->size=recovery_num;
  ne = wa->alphabet->size;
  nt = wa->table->numTransitions;
  if (wa->table->table_type == DENSE) {
    size_t s;
    /* it is a pity that the table was output in transposed form! */
    fsa_table_init(wa->table, recovery_num, ne);
    table = wa->table->table_data_ptr;
    for (i = 1; i <= recovery_num; i++)
      for (j = 1; j <= ne; j++) {
	if (i==recovery_num)
		table[j][i]=0;
	else
        	 s = fread((void *)(table[j] + i), sizeof(int), (size_t)1, tempfile);
        if (*(table[j]+i) == 1000000000)
                *(table[j]+i)=recovery_num;
        if (*(table[j]+i) > recovery_num)
                *(table[j]+i)=0;
        (void)s; // HACK to silence compiler warning
      }
  }
  } // compressed_transitions
  fclose(tempfile);

  //unlink(tempfilename);
  if (no_accepts>0)
  	return wa;
  return NULL;
}