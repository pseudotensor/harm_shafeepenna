

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// MNEMONICS or other things that should rarely change, or things that on depend on the above items (e.g. if statements, loops, etc.)
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//equals to unity if the interpolation of gamma is performed and if requested to use prim. reduction
#define STORE_GAMMA_PRIM_REDUCTION_FRACTION (WENO_USE_PRIM_REDUCTION && (VARTOINTERP == PRIMTOINTERP_3VEL_GAMMA || VARTOINTERP == PRIMTOINTERP_RHOV_GAMMA || VARTOINTERP == PRIMTOINTERP_3VELREL_GAMMAREL) )



//
// below for para2() and para3()
//

#define PMC 100
// PMC: originally used in PPM paper -- Xiaoyue says fails on nonlinear tests since has no limiter -- but they used it, right?

#define MINM_STEEPENER 101

// whether to allow extremum in para that has monotonicity check at end
#define PARALINE2EXTREME 1
#define PARALINE2LIM MC

#if(WHICHPARA==PARA2)
#define PARA2LIM VANL
// jon tests show for PARA2:
// PARA2LIM = MC completely fails
// PARA2LIM = VANL best ok
// PARA2LIM = PMC kinda ok
#elif(WHICHPARA==PARA3)
// MC recommended by Xiaoyue
#define PARA2LIM MC
#elif(WHICHPARA==PARA4)
// MC recommended by Xiaoyue
#define PARA2LIM MC
#endif


#if( MERGEDC2EA2CMETHOD )
#define NUM_CVT_TYPES 8  //4 usual recon types + 4 derivatives
#else
#define NUM_CVT_TYPES 4
#endif



// whether to check if rho<=0 or u<=0 when restarting.
#if( (EOMTYPE==EOMGRMHD)||(EOMTYPE==EOMCOLDGRMHD) )
#define CHECKRHONEGZERORESTART 1
#else
#define CHECKRHONEGZERORESTART 0
#endif


#if((WHICHCURRENTCALC==0)||(WHICHCURRENTCALC==2))
#define NUMCURRENTSLOTS 5
#elif(WHICHCURRENTCALC==1)
#define NUMCURRENTSLOTS 6
#endif

#if(SPLITPRESSURETERMINFLUXMA)
// where to put gas pressure part of flux that is T^i_i
#define FLUXSPLITPMA(dir) (B1+dir-1) // put in unused magnetic field part
#else
#define FLUXSPLITPMA(dir) (-100) // so FLUXSPLITPMA(dir)==pl is always 0
#endif

#if(SPLITPRESSURETERMINFLUXEM)
// where to put magnetic pressure part of flux that is T^i_i
#define FLUXSPLITPEM(dir) (RHO) // put in unused mass advection part
#else
#define FLUXSPLITPEM(dir) (-100) // so FLUXSPLITPEM(dir)==pl is always 0
#endif


#if((USESTOREDSPEEDSFORFLUX==0)||(STOREWAVESPEEDS==0))
#define USEGLOBALWAVE 0
#else
#define USEGLOBALWAVE 1
#endif


#define COMPUTE4FIELDforFLUX (MAXWELL==GENMAXWELL || USEGLOBALWAVE==0)



#if(MCOORD==CARTMINKMETRIC)
#define ISTORE(i) (0)
#define JSTORE(j) (0)
#define KSTORE(k) (0)
#else
#define ISTORE(i) (i)
#define JSTORE(j) (j)
#define KSTORE(k) (k)
#endif



// setup for various boundary situations
// so doesn't produce differences in irrelevant directions, whether boundary zones or not
// mac(?) macros are for use in definitions of other macros since macro with args needs to be directly a function of what's hardcoded, not some "replacement" since nothing is replaced, not to be used inside code for USING a macro.

#if(N1>1)
#define im1 (i-1)
#define im1mac(i) (i-1)
#define ip1 (i+1)
#define ip1mac(i) (i+1)

#define irefshift (2*N1-i-1)

#else

#define im1 (i)
#define im1mac(i) (i)
#define ip1 (i)
#define ip1mac(i) (i)

#define irefshift (i)

#endif

#if(N2>1)
#define jm1 (j-1)
#define jm1mac(j) (j-1)
#define jp1 (j+1)
#define jp1mac(j) (j+1)
#define jrefshift (2*N2-j-1)
#else
#define jm1 (j)
#define jm1mac(j) (j)
#define jp1 (j)
#define jp1mac(j) (j)
#define jrefshift (j)
#endif

#if(N3>1)
#define km1 (k-1)
#define km1mac(k) (k-1)
#define kp1 (k+1)
#define kp1mac(k) (k+1)
#define krefshift (2*N3-k-1)
#else
#define km1 (k)
#define km1mac(k) (k)
#define kp1 (k)
#define kp1mac(k) (k)
#define krefshift (k)
#endif


// used to tell if N>1 or not (can just ask directly)
#define N1NOT1 ((N1>1) ? 1 : 0)
#define N2NOT1 ((N2>1) ? 1 : 0)
#define N3NOT1 ((N3>1) ? 1 : 0)

// GODMARK: looks like I can set just above and rest are set for me for any case

// GODMARK: check these new conditions

// restrict loops only over relevant domain in reduced dimension case



// 3 maximum boundary zones needed if doing Parabolic interpolation
// maximum number of boundary zones needed for all calculations

// have to set manually if going to set DOENOFLUX at runtime.



// x1
#define SHIFT1 N1NOT1
#define N1BND MAXBND*N1NOT1

#define INFULL1 (-N1BND)
#define INFULLP11 (-N1BND+SHIFT1)
#define OUTFULL1 (N1-1+N1BND) 
#define OUTFULLM11 (N1-1+N1BND-SHIFT1) 
#define OUTFULLP11 (N1-1+N1BND+SHIFT1)

#define INHALF1 (-N1BND/2)
#define OUTHALF1 (N1-1+N1BND/2)
#define INP11 (-N1BND+SHIFT1)
#define OUTP11 (N1-1+N1BND-SHIFT1)
#define INM1 -SHIFT1
#define OUTM1 N1-1+SHIFT1

// unlike other loops limits that should reduce to 0 when the N=1 to as if like dimension didn't exist,
// this one should force loop to not happen at all when N=1 since acts on boundary zones don't exist
#define INBOUNDLO1 (-N1BND)
#define INBOUNDHI1 (-1)

#define OUTBOUNDLO1 (N1)
#define OUTBOUNDHI1 (N1+N1BND-1)

//#define INFACEBOUNDLO1 (-N1BND) // (-N1BND+1) // GODMARK: large domain used for easy checking of fluxes after bound_flux().
//#define INFACEBOUNDHI1 (-1+SHIFT1)

// up to -1 since 0 is actually defined with original primitives
// generalize (expand) a bit:
#define INFACEBOUNDLO1 (-N1BND)
#define INFACEBOUNDHI1 (-1)


// from N1+1 since N1 is actually defined with original primitives (only true if not FLUXCTSTAG)
// generalize (expand) a bit:
//#define OUTFACEBOUNDLO1 (N1+1)
#define OUTFACEBOUNDLO1 (N1)
#define OUTFACEBOUNDHI1 (N1+N1BND-1)

// x2
#define SHIFT2 N2NOT1
#define N2BND MAXBND*N2NOT1

#define INFULL2 (-N2BND)
#define INFULLP12 (-N2BND+SHIFT2)
#define OUTFULL2 (N2-1+N2BND)
#define OUTFULLM12 (N2-1+N2BND-SHIFT2) 
#define OUTFULLP12 (N2-1+N2BND+SHIFT2)

#define INHALF2 (-N2BND/2)
#define OUTHALF2 (N2-1+N2BND/2)
#define INP12 (-N2BND+SHIFT2)
#define OUTP12 (N2-1+N2BND-SHIFT2)
#define INM2 -SHIFT2
#define OUTM2 N2-1+SHIFT2

// unlike other loops limits that should reduce to 0 when the N=1 to as if like dimension didn't exist,
// this one should force loop to not happen at all when N=1 since acts on boundary zones don't exist
#define INBOUNDLO2 (-N2BND)
#define INBOUNDHI2 (-1)

#define OUTBOUNDLO2 (N2)
#define OUTBOUNDHI2 (N2+N2BND-1)

//#define INFACEBOUNDLO2 (-N2BND)
//#define INFACEBOUNDHI2 (-1+SHIFT2)

// generalize (expand) a bit:
//#define INFACEBOUNDLO2 (-N2BND+1)
#define INFACEBOUNDLO2 (-N2BND)
#define INFACEBOUNDHI2 (-1)


// generalize (expand) a bit:
//#define OUTFACEBOUNDLO2 (N2+1)
#define OUTFACEBOUNDLO2 (N2)
#define OUTFACEBOUNDHI2 (N2+N2BND-1)


// x3
#define SHIFT3 N3NOT1
#define N3BND MAXBND*N3NOT1

#define INFULL3 (-N3BND)
#define INFULLP13 (-N3BND+SHIFT3)
#define OUTFULL3 (N3-1+N3BND)
#define OUTFULLM13 (N3-1+N3BND-SHIFT3)
#define OUTFULLP13 (N3-1+N3BND+SHIFT3)

#define INHALF3 (-N3BND/2)
#define OUTHALF3 (N3-1+N3BND/2)
#define INP13 (-N3BND+SHIFT3)
#define OUTP13 (N3-1+N3BND-SHIFT3)
#define INM3 -SHIFT3
#define OUTM3 N3-1+SHIFT3

// unlike other loops limits that should reduce to 0 when the N=1 to as if like dimension didn't exist,
// this one should force loop to not happen at all when N=1 since acts on boundary zones don't exist
#define INBOUNDLO3 (-N3BND)
#define INBOUNDHI3 (-1)

#define OUTBOUNDLO3 (N3)
#define OUTBOUNDHI3 (N3+N3BND-1)

//#define INFACEBOUNDLO3 (-N3BND)
//#define INFACEBOUNDHI3 (-1+SHIFT3)

// generalize (expand) a bit:
//#define INFACEBOUNDLO3 (-N3BND+1)
#define INFACEBOUNDLO3 (-N3BND)
#define INFACEBOUNDHI3 (-1)


// generalize (expand) a bit:
//#define OUTFACEBOUNDLO3 (N3+1)
#define OUTFACEBOUNDLO3 (N3)
#define OUTFACEBOUNDHI3 (N3+N3BND-1)




/* allocated memory uses this for active zones 0-N1-1 and bc beyond that */
#define N1M (N1+N1BND*2)
#define N2M (N2+N2BND*2)
#define N3M (N3+N3BND*2)


/* NBIG is bigger of N1 and N2 and N3 */
#define NBIG1 ((N1>N2) ? N1 : N2)
#define NBIG  ((NBIG1>N3) ? NBIG1 : N3)

#define NBIGBND1 ((N1BND>N2BND) ? N1BND : N2BND)
#define NBIGBND  ((NBIGBND1>N3BND) ? NBIGBND1 : N3BND)

// N?OFF and N?NOT1 are a bit redundant
//#define N1OFF (((N1BND>0)&&(N1>1)) ? 1 : 0)
//#define N2OFF (((N2BND>0)&&(N2>1)) ? 1 : 0)
//#define N3OFF (((N3BND>0)&&(N3>1)) ? 1 : 0)



/* NBIGM is bigger of N1M and N2M and N3M */
// not currently used
#define NBIG1M ((N1M>N2M) ? N1M : N2M)
#define NBIGM  ((NBIG1M>N3M) ? NBIG1M : N3M)

// surface areas of sides WITH boundary zones
// reduces to 0 if that dimension not used
#define MAXSURFA1 (N2M*N3M*N1NOT1)
#define MAXSURFA2 (N1M*N3M*N2NOT1)
#define MAXSURFA3 (N1M*N2M*N3NOT1)

// maximal surface of boundary exchange
// notice that this works in any number of dimensions and any N1,N2,N3
// that is, it reduces correctly when a dimension is degenerate (N=1)
#define NBIGS1M ((MAXSURFA1>MAXSURFA2) ? MAXSURFA1 : MAXSURFA2)
#define NBIGSM ((NBIGS1M>MAXSURFA3) ? NBIGS1M : MAXSURFA3)

#if(LIMADJUST!=LIMITERFIXED && FLUXADJUST!=FLUXFIXED)
#define NUMPFLAGSBOUND (NUMPFLAGS) // all
#elif(LIMADJUST!=LIMITERFIXED)
#define NUMPFLAGSBOUND (1+FLAGREALLIM) // only 0,1
#else
#define NUMPFLAGSBOUND (1+FLAGUTOPRIMFAIL) // only 0
#endif

// GODMARK: Could make a volume that is not NBIGBND*NBIGSM but may be smaller?
// used in init_mpi.c for workbc and workbc_int


// number of interpolation variables for staggered field method per B and v each
#define NUMCORNINTERP 4

// for wavespeeds.c and fluxcompute.c
#define NUMCS 2
#define CMIN 0
#define CMAX 1

#define MINMAX(q,a,b) ( ((q)==CMIN) ? MIN(a,b) : MAX(a,b) )





#if((DODISS||DODISSVSR)&&(DOENTROPY==DONOENTROPY))
#error Turn on entropy evolution if want dissipation
#endif




// processed version of npr definition since no immediate evaluation in standard c preprocessor
// this "global.defnprs.h" is generated by the program generatenprs.c
#include "global.defnprs.h"





#define NMAXBOUND ((NPRBOUND>NFLUXBOUND) ? NPRBOUND : NFLUXBOUND)





// cent,face1,face2,face3,corn


// NPR*4 = 1 dUgeom and 3 dUriemanns ; 3 directions for F1,F2,F3 and pl pr and F(pl) and F(pr)
#define NUMFLUXDUMP (NPR*4 + NPR*3*(1+2+2))


// 4 space-time directions with only spatial parts used for now
#define NUMVPOTDUMP 4



// size of certain dumped tavg quantities

// was 29 =>
#define NUMNORMDUMP (NPR+1+4*4+6) // number of "normal" dump variables
// for above see diag.c and set_varstavg()
#define NUMFARADAY 6
#define NUMOTHER 1
#define NUMSTRESSTERMS (NUMFLUXTERMS*NDIM*NDIM)

/** GLOBAL ARRAY SECTION **/






// size of data type used for all floats
#define FLOATTYPE 0
#define DOUBLETYPE 1
#define LONGDOUBLETYPE 2
#define LONGLONGINTTYPE 3


#if(REALTYPE>SENSITIVE)
god=deathadflkjasdflkjasdlfkja242424
#endif

#ifndef FLT_EPSILON
#define FLT_EPSILON (1.19209290e-07F)
#endif

#ifndef DBL_EPSILON
#define DBL_EPSILON (2.2204460492503131e-16L)
#endif

#ifndef LDBL_EPSILON
#define LDBL_EPSILON (1.08420217248550443401e-19L)
#endif


// need not change below datatype stuff
#if(REALTYPE==FLOATTYPE)
#define MINNUMREPRESENT FLT_MIN
#define NUMEPSILON FLT_EPSILON
#elif(REALTYPE==DOUBLETYPE)
#define MINNUMREPRESENT DBL_MIN
#define NUMEPSILON DBL_EPSILON
#elif(REALTYPE==LONGDOUBLETYPE)
#define MINNUMREPRESENT DBL_MIN
#define NUMEPSILON LDBL_EPSILON
#endif

#if(SENSITIVE==FLOATTYPE) // for sensitive counters
#define SFTYPE float
#elif(SENSITIVE==DOUBLETYPE)
#define SFTYPE double
#elif(SENSITIVE==LONGDOUBLETYPE)
#define SFTYPE long double
#endif


// used for numerical differencing
#define NUMSQRTEPSILON (sqrt(NUMEPSILON))
#define NUMSQEPSILON (NUMEPSILON*NUMEPSILON)



#if(COUNTTYPE==DOUBLETYPE)
#define CTYPE double
#define CTYPEHEADERONEOUT "%21.15g"
#define CTYPEHEADERONEIN "%lf"
#elif(COUNTTYPE==LONGLONGINTTYPE)
#define CTYPE long long int
#define CTYPEHEADERONEOUT "%lld"
#define CTYPEHEADERONEIN "%lld"
#endif


#if(PFLAGTYPE==INTTYPE)
#define PFTYPE signed int
#define PFTYPEHEADERONEOUT "%d"
#define PFTYPEHEADERONEIN "%d"
#elif(PFLAGTYPE==CHARTYPE)
  // on Harvard's BlueGene/L char is by default unsigned!, so force signed as required by our code
#define PFTYPE signed char
#define PFTYPEHEADERONEOUT "%d"
#define PFTYPEHEADERONEIN "%d"
#endif


// GODMARK: NUMENERVAR outdated?
#define NUMENERVAR (6+NPR+NPR+3)


/* numerical convenience */
#if(REALTYPE==FLOATTYPE)
#define VERYBIG (1.e30)
#else
#define VERYBIG (1.e150)
#endif

#define BIG (1.e+50)
#define SMALL	(1.e-50)

#define SLEPSILON (1.e-6)


/* size of step in numerical derivative evaluations */
#define HSTEP	(1.e-5)




#if(SENSITIVE==LONGDOUBLETYPE)
#define SFTYPEHEADERONEIN "%Lf"
// assume sensitive>=realtype in precision
#if(REALTYPE==LONGDOUBLETYPE) // was FLOATTYPE==REALTYPE and SENS=DOUBLETYPE
#define HEADERONEIN "%Lf"
#define HEADER2IN "%Lf %Lf"
#define HEADER3IN "%Lf %Lf %Lf"
#define HEADER4IN "%Lf %Lf %Lf %Lf"
#define HEADER5IN "%Lf %Lf %Lf %Lf %Lf"
#define HEADER6IN "%Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER7IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER8IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER9IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define RESTARTHEADER "%d %d %d "\
	      "%Lf %Lf %ld %Lf %Lf %Lf %Lf %Lf "\
	      "%Lf %d %d %d %d %d %d %d %d "\
	      "%Lf %Lf %Lf %Lf %d "\
              "%d %d %d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#elif(REALTYPE==DOUBLETYPE)
#define HEADERONEIN "%lf"
#define HEADER2IN "%lf %lf"
#define HEADER3IN "%lf %lf %lf"
#define HEADER4IN "%lf %lf %lf %lf"
#define HEADER5IN "%lf %lf %lf %lf %lf"
#define HEADER6IN "%lf %lf %lf %lf %lf %lf"
#define HEADER7IN "%lf %lf %lf %lf %lf %lf %lf"
#define HEADER8IN "%lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER9IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define RESTARTHEADER "%d %d %d "\
	      "%Lf %Lf %ld %Lf %Lf %Lf %lf %lf "\
	      "%Lf %d %d %d %d %d %d %d %d "\
	      "%lf %lf %lf %lf %d "\
              "%d %d %d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#elif(REALTYPE==FLOATTYPE)
#define HEADERONEIN "%f"
#define HEADER2IN "%f %f"
#define HEADER3IN "%f %f %f"
#define HEADER4IN "%f %f %f %f"
#define HEADER5IN "%f %f %f %f %f"
#define HEADER6IN "%f %f %f %f %f %f"
#define HEADER7IN "%f %f %f %f %f %f %f"
#define HEADER8IN "%f %f %f %f %f %f %f %f"
#define HEADER9IN "%f %f %f %f %f %f %f %f %f"
#define RESTARTHEADER "%d %d %d "\
	      "%Lf %Lf %ld %Lf %Lf %Lf %f %f "\
	      "%Lf %d %d %d %d %d %d %d %d "\
	      "%f %f %f %f %d "\
              "%d %d %d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f"
#endif

#elif(SENSITIVE==DOUBLETYPE)
#define SFTYPEHEADERONEIN "%lf"
// assume sensitive>=realtype in precision
#if(REALTYPE==DOUBLETYPE)
#define HEADERONEIN "%lf"
#define HEADER2IN "%lf %lf"
#define HEADER3IN "%lf %lf %lf"
#define HEADER4IN "%lf %lf %lf %lf"
#define HEADER5IN "%lf %lf %lf %lf %lf"
#define HEADER6IN "%lf %lf %lf %lf %lf %lf"
#define HEADER7IN "%lf %lf %lf %lf %lf %lf %lf"
#define HEADER8IN "%lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER9IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define RESTARTHEADER "%d %d %d "\
	      "%lf %lf %ld %lf %lf %lf %lf %lf "\
	      "%lf %d %d %d %d %d %d %d %d "\
	      "%lf %lf %lf %lf %d "\
              "%d %d %d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#elif(REALTYPE==FLOATTYPE)
#define HEADERONEIN "%f"
#define HEADER2IN "%f %f"
#define HEADER3IN "%f %f %f"
#define HEADER4IN "%f %f %f %f"
#define HEADER5IN "%f %f %f %f %f"
#define HEADER6IN "%f %f %f %f %f %f"
#define HEADER7IN "%f %f %f %f %f %f %f"
#define HEADER8IN "%f %f %f %f %f %f %f %f"
#define HEADER9IN "%f %f %f %f %f %f %f %f %f"
#define RESTARTHEADER "%d %d %d "\
	      "%lf %lf %ld %lf %lf %lf %f %f "\
	      "%lf %d %d %d %d %d %d %d %d "\
	      "%f %f %f %f %d "\
              "%d %d %d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f"
#endif

#elif(SENSITIVE==FLOATTYPE)
#define SFTYPEHEADERONEIN "%f"

#if(REALTYPE==DOUBLETYPE)
#define RESTARTHEADER "" // dumb, so crash on compile
#elif(REALTYPE==FLOATTYPE)
#define HEADERONEIN "%f"
#define HEADER2IN "%f %f"
#define HEADER3IN "%f %f %f"
#define HEADER4IN "%f %f %f %f"
#define HEADER5IN "%f %f %f %f %f"
#define HEADER6IN "%f %f %f %f %f %f"
#define HEADER7IN "%f %f %f %f %f %f %f"
#define HEADER8IN "%f %f %f %f %f %f %f %f"
#define HEADER9IN "%f %f %f %f %f %f %f %f %f"
#define RESTARTHEADER "%d %d %d "\
	      "%f %f %ld %f %f %f %f %f "\
	      "%f %d %d %d %d %d %d %d %d "\
	      "%f %f %f %f %d "\
              "%d %d %d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f"
#endif

#endif


// 2
// 6
// 10
// 3
// 5
// 4
// SUM=30
// 23+10=33
// total=63
#define WRITERESTARTHEADER "%d %d %d " \
		 "%26.20g %26.20g %ld %26.20g %26.20g %26.20g %26.20g %26.20g " \
		 "%26.20g %d %d %d %d %d %d %d %d " \
		 "%26.20g %26.20g %26.20g %26.20g %d " \
                 "%d %d %d %d %d %d " \
                 "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g %26.20g "

#define HEADERONEOUT "%26.20g"
#define SFTYPEHEADERONEOUT "%31.26g"




