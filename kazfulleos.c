////////////////////////////
//
// EOS from KAZ FULL
//
////////////////////////////



// TODO:
//
// 0) Make sure icurr,jcurr,kcurr,pcurr set before EOS called!!!! GODMARK
// 1) choose case where checking if within EOS is not done (assume always there and truncate instead of extend)
//    TWO if's removed then: if(iswithin_eostable) and in each type of calculation the call: if(get_eos_fromtable()) are removed : might speed things up
// 2) Is my interpolation best/fastest/correct?  Thompson et al. (2003) use bivariate interpolation that looks alot more complicated


// TODO NEW:
// 1) New H calculation (MAKE MPI)
// 3) check additional code that's presently in debug mode
// 4) check MPI stuff here and for gravity
// 5) change bounds so loops over bounding vars rather than only up to B3
// 6) rescale() depends often on final quantity being B3

//////////////////////////////
//
// Some often changable variables
//
///////////////////////////////

// whether to allow use of full table (if 0, then others must be turned on)
#define ALLOWFULLTABLE 1

// whether to only use full table (0) for allow use of simple tables if can (1)
#define ALLOWSIMPLETABLE 1

// whether to use simplezoom table if can
// zoom not needed anymore with new degen offset method
// indeed, set to 0 unless make new simplezoom table
#define ALLOWSIMPLEZOOMTABLE 0

// how many dimensions to consider.  Other dimensions' values will be consider as the dimension's lowest value
#define WHICHEOSDIMEN 4


// whether to use degen offset (otherwise assume degen offset from file is 0 even if read-in differently)
#define ALLOWDEGENOFFSET 1


// whether to check if table returns a valid EOS value by using existence of stored inversion to temperature
// so far only setup for F(rho0,u)
// if some invalid, then don't use those data points.  If all surrounding points are invalid, then use them as if they were valid
#define CHECKIFVALIDEOSDATA 1


// using log interpolation results in much smoother results, and consistent with eos_extract.m for interpolation
// That is, using integer position is log-interp since all independents are log on the grid
// And those functions in eos_extract.m interpolated as log are here interpolated as log
// 0 or 1
#define DOLOGINTERP 1


// which EOS to reduce to if beyond table
// mignone doesn't make sense
#define REDUCE2WHICHEOS IDEALGAS
// ensure that gamideal is chosen
// GODMARK: Could choose nearest tabulated value of dp/du|rho0 and dp/dchi|rho0 for gamideal when indeps are rho0,u and rho0,chi


//////////////////////////////
//
// More permenant variables
//
///////////////////////////////

// assume degen table is always stored along with corresponding normal table
#define NUMTBLS 3
#define FULLTABLE 0
#define SIMPLETABLE 1
#define SIMPLEZOOMTABLE 2


// 5 true dimensions (rhob,u/p/chi/tk,hcm,tdynorye,tdynorynu)
#define NUMINDEPDIMENS 5
// what is contained within file (rhob,utotdiff,ptotdiff,chidiff,hcm,tdynorye,tdynorynu)
#define NUMEOSINDEPS 7


#define NUMEXTRAINDEP 3 // other besides rho and u-versions

// these are fixed in order and number as consistent with what's read-in from file
#define RHOEOS 0 // rest-mass density
#define UEOS   1    // internal energy density: used for u
#define PEOS   2    // pressure energy density: used for p
#define CHIEOS 3    // enthalpy energy density: used for \chi
#define TEOS   4    // tdynorye: dynamical time or Y_e
#define YNUEOS 5    // tdynorynu: dynamical time or Y_\nu
#define HEOS   6    // height (not used normally anymore)


// below ending # corresponds to whichrnpmethod
// monotonized is 22 + extras
// for input, total is 7 + extras
// totals for input are 8 23 16 is present original data for version1,2,3, then eos_extract adds:
// for normal table, eos_extract.m adds 3 (degens) + 3 (tk's) = 6 total normal added
// totals for input after eos_extract are: 7+extra+6
// eos_extract always makes 24 + extras = 10 iterators + 14 eosquantities + extras
// for degen table, eox_extract.m has 9 total

// number of base quantities to *store* from table made by eos_extract.m
#define NUMEOSQUANTITIESBASE (2 + 2 + 1 + 3 + 3 + 3)  // 14 : eosquantities

// for memory optimization, specifiy which datatype
#define WHICHDATATYPEGENERAL 4


#if(WHICHEOSDIMEN==4 && WHICHDATATYPEGENERAL!=4)
#error WHICHEOSDIMEN and WHICHDATATYPEGENERAL inconsistent
#endif

// for different datatypes have different extra things
#define MAXNUMDATATYPES 4
#define NUMEOSQUANTITIESTYPE1 (NUMEOSQUANTITIESBASE + 1)  // (Full EOS, rho,T,H)
#define NUMEOSQUANTITIESTYPE2 (NUMEOSQUANTITIESBASE + 16) // (non-neutrino EOS, rho,T,Y_e)
#define NUMEOSQUANTITIESTYPE3 (NUMEOSQUANTITIESBASE + 11) // (Full EOS, rho,T,Y_e,Y_nu but H-dependent)
#define NUMEOSQUANTITIESTYPE4 (NUMEOSQUANTITIESBASE + 24) // (non-neutrino EOS, rho,T,Y_e,Y_\nu)

#if(WHICHDATATYPEGENERAL==3)
#define NUMEOSQUANTITIESMEM NUMEOSQUANTITIESTYPE3
#else
#define NUMEOSQUANTITIESMEM (MAX(MAX(MAX(NUMEOSQUANTITIESTYPE1,NUMEOSQUANTITIESTYPE2),NUMEOSQUANTITIESTYPE3),NUMEOSQUANTITIESTYPE4))
#endif

int numeosquantities[NUMTBLS]; // true number
int numeosdegenquantities[NUMTBLS]; // true number

// these are fixed in order and number from read-in file
#define PofRHOU 0      // p(rho0,u)
#define UofRHOP 1      // u(rho0,p)

#define DPDRHOofRHOU 2 // dpdrho0 |u (rho0,u)
#define DPDUofRHOU 3   // dp/du |rho0 (rho0,u)

#define CS2ofRHOU 4    // cs^2(rho0,u)

#define SofRHOU 5      // S(rho0,u)
#define DSDRHOofRHOU 6 // dS/drho0 |u (rho0,u)
#define DSDUofRHOU 7   // dS/du |\rho0 (rho0,u)

#define PofRHOCHI 8    // P(rho0,\chi)  \chi = u+p
#define IDRHO0DP 9     // 1/(d\rho0/dp) |\chi
#define IDCHIDP 10      // 1/(d\chi/dp) |\rho0

#define TEMPU 11 // temperature in Kelvin (doesn't need to be function of H or TDYNORYE, but can change later)
#define TEMPP 12 // temperature in Kelvin (doesn't need to be function of H or TDYNORYE, but can change later)
#define TEMPCHI 13 // temperature in Kelvin (doesn't need to be function of H or TDYNORYE, but can change later)

// so-called "extras" in eos_extract.m: Those things that didn't require extra processing as independent variables or derivatives -- just interpolated from T -> U only
// extras:

#define EXTRA1  14
#define EXTRA2  15
#define EXTRA3  16
#define EXTRA4  17
#define EXTRA5  18
#define EXTRA6  19
#define EXTRA7  20
#define EXTRA8  21
#define EXTRA9  22
#define EXTRA10 23
#define EXTRA11 24
#define EXTRA12 25
#define EXTRA13 26
#define EXTRA14 27
#define EXTRA15 28
#define EXTRA16 29
#define EXTRA17 30
#define EXTRA18 31
#define EXTRA19 32
#define EXTRA20 33
#define EXTRA21 34
#define EXTRA22 35
#define EXTRA23 36
#define EXTRA24 37

// used to map request to correct EXTRA for given table that uses certain whichdatatype
#define LAMBDATOT -100
#define QDOTNU -101


// whichdatatype==1
// EXTRA1: Neutrino cooling rate (erg/s/cm^2)
#define DATATYPE1_EXTRAFINAL EXTRA1

// whichdatatype==2
// EXTRA1:  qtautelohcm
// EXTRA2:  qtauaelohcm
// EXTRA3:  qtautmuohcm
// EXTRA4:  qtauamuohcm
// EXTRA5:  qtauttauohcm
// EXTRA6:  qtauatauohcm
// EXTRA7:  ntautelohcm
// EXTRA8:  ntauaelohcm
// EXTRA9:  ntautmuohcm
// EXTRA10:  ntauamuohcm
// EXTRA11:  ntauttauohcm
// EXTRA12:  ntauatauohcm
// EXTRA13:  gammapeglobal+gammaAeglobal
// EXTRA14:  gammapnuglobal+gammapenuglobal
// EXTRA15:  gammanglobal + gammaneglobal
// EXTRA16: gammannuglobal
#define DATATYPE2_EXTRAFINAL EXTRA16

// whichdatatype==3
// for now this is opimal choice for simplicity, although big table unresolved in H and Ynu
// EXTRA1:  Qphoton
// EXTRA2:  Qm
// EXTRA3:  graddotrhouyl
// EXTRA4:  Tthermaltot
// EXTRA5:  Tdifftot
// EXTRA6:  lambdatot
// EXTRA7:  lambdaintot
// EXTRA8:  Enuglobal
// EXTRA9:  Enueglobal
// EXTRA10:  Enuebarglobal
// EXTRA11: Ynuthermal
#define DATATYPE3_EXTRAFINAL EXTRA11

// whichdatatype==4
// GODMARK: if this is going to work, need also the energy density parts as functions of \chi, but for now doesn't seem this method will be useful due to need to iterate within table to find du from u
// EXTRA1:  qtautnueohcm 5E7
// EXTRA2:  qtauanueohcm 1E-50
// EXTRA3:  qtautnuebarohcm 1E36
// EXTRA4:  qtauanuebarohcm 1E36
// EXTRA5:  qtautmuohcm 1E-9
// EXTRA6:  qtauamuohcm 1E-14
// EXTRA7:  ntautnueohcm 5E7
// EXTRA8:  ntauanueohcm 3E-48
// EXTRA9:  ntautnuebarohcm 2E35
// EXTRA10:  ntauanuebarohcm 2E35
// EXTRA11:  ntautmuohcm 1E-9
// EXTRA12:  ntauamuohcm 1E-14
// EXTRA13:  unue0 8E59
// EXTRA14:  unuebar0 3E-26
// EXTRA15:  unumu0 4E25
// EXTRA16:  nnue0 1E57
// EXTRA17:  nnuebar0 3E-20
// EXTRA18:  nnumu0 1E31
// EXTRA19:  lambdatot 3E-34
// EXTRA20:  lambdaintot 3E-34
// EXTRA21:  tauphotonohcm 2E9
// EXTRA22:  tauphotonabsohcm 2E9
// EXTRA23:  nnueth0
// EXTRA24:  nnuebarth0
#define DATATYPE4_EXTRAFINAL EXTRA24




// processed quantities
// MAXPROCESSEDEXTRAS is set in nondepnmenmics.h
#define QPHOTON 0
#define QNEUTRINO 1
#define GRADDOTRHOUYL 2
#define TTHERMAL 3
#define TDIFF 4
#define RHONU 5
#define PNU 6
#define SNU 7
#define YNULOCAL 8
#define YNUTHERMAL 9
#define ENUAVG 10
#define ENUE 11
#define ENUEBAR 12



// degen offset quantities for independent variables utot,ptot, chi
#define NUMEOSDEGENQUANTITIESMEM (3) // utotdiff, ptotdiff, chidiff

#define NUMEOSDEGENINDEPS (NUMINDEPDIMENS-1) // rho, tdynorye, tdynorynu, H

// independents for degen table
#define RHOEOSDEGEN 0
#define TEOSDEGEN 1
#define YNUEOSDEGEN 2
#define HEOSDEGEN 3


// these are fixed in order and number from read-in file
// so can access functional degentable of independent variables by whichindep-1
#define UTOTDIFF (UEOS-1)   // should always resolve to: 0
#define PTOTDIFF (PEOS-1)   // :1
#define CHIDIFF  (CHIEOS-1) // :2







//////////////
//
// here:
// [independent variables]
// rhob = rho0 -- rest-mass density in g/cc
// tk = Temperature in Kelvin
// hcm = height of medium in cm
// tdynorye = dynamical time in seconds (assumed NSE time) OR Y_e
//
// [dependent variables]
// p_tot = total pressure in cgs units
// u_tot = u = internal relativistic comoving energy (no rest-mass) in cgs units
// s_tot = total entropy density in comoving frame in cgs units






///////////////
//
// EOSN1=rho0
// EOSN2=u or p or chi
// EOSN3=tdynorye in seconds or Y_e
// EOSN4=tdynorynu in seconds or Y_\nu
// EOSN5=height in cm

// full EOS table w/ neutrino part to be computed during run-time
#define EOSN1 100
#define EOSN2 50
#define EOSN3 50
#define EOSN4 20
#define EOSN5 1         // H not tabulated

// EOS table with assumed Height (small so optically thin problem) and assumed Y_e(\rhob,T)
// degenerate table assumes to be same size except N2=1
#define EOSSIMPLEN1 50
#define EOSSIMPLEN2 50
#define EOSSIMPLEN3 50
#define EOSSIMPLEN4 1   // here Y_\nu~0 (optically thin)
#define EOSSIMPLEN5 1   // H not tabulated

// EOS table with assumed Height (small) and assumed tdynorye (large)
// degenerate table assumes to be same size except N2=1
// NOT USED RIGHT NOW
#define EOSSIMPLEZOOMN1 1
#define EOSSIMPLEZOOMN2 1
#define EOSSIMPLEZOOMN3 1
#define EOSSIMPLEZOOMN4 1
#define EOSSIMPLEZOOMN5 1

// GODMARK: could have a table for Ynu=thermalized and have an array that stores when source term forces Ynu to be perfectly thermal, and use that table in that case.
// generating it now




///////////////
//
// Table limits

#define UPDOWN 2 // 0=down 1=up

#define TBLITEMS (UPDOWN+2+2)
#define TBLLINEARITEMS (UPDOWN+2)

// first [4] : above 4 types of independent variables
// second [4] : 0 = lower log_base limit, 1 = upper log_base limit, 2=step, 3 = divisor of grid position 4=base of log, 5 = log10 value of offset for log_base stepping so can control how resolved
FTYPE inputtablelimits[NUMTBLS][NUMEOSINDEPS][TBLITEMS];
// below is same as input, but converted for easy use
FTYPE tablelimits[NUMTBLS][NUMEOSINDEPS][TBLITEMS];
// first [4] as above, second [2] : 0=lower linear limit, 1=upper linear limit, 2 = log(-base) 3 = linear offset
FTYPE lineartablelimits[NUMTBLS][NUMEOSINDEPS][TBLLINEARITEMS];
// first [4] : as first [4] above
int tablesize[NUMTBLS][NUMEOSINDEPS];

int vartypearray[NUMINDEPDIMENS+1];
int varnormalcompare2degentypearray[NUMEOSDEGENINDEPS+1]; // to be used to compare against degen table
int vardegentypearray[NUMEOSDEGENINDEPS+1];



////////////////////////////////////
//
// Allocate space for table
//
///////////////////////////////////////////////
//
// This entire file assumes all things are FTYPE except eostable itself that is explicitly double
//
///////////////////////////////////////////////

// presumed to be several functions as functions of 4 other quantities that are the independent variables (UEOS,PEOS,CHIEOS share same index)
// Notice that eostable has density as EOSN1 so fastest indicies are related to density and pressure-related quantities since this will result in fastest lookup for fixed H and Tdynorye
// due to float limit of 1E38, need doubles
// the 1 on N2 is where the DEGEN table is stored for this EOS table

#if(WHICHEOS==KAZFULL)
double a_eostable[NUMEOSQUANTITIESMEM][EOSN5][EOSN4][EOSN3][EOSN2][EOSN1];
#endif
static double (*eostable)[EOSN5][EOSN4][EOSN3][EOSN2][EOSN1];

#if(WHICHEOS==KAZFULL)
double a_eosdegentable[NUMEOSDEGENQUANTITIESMEM][EOSN5][EOSN4][EOSN3][1][EOSN1];
#endif
static double (*eosdegentable)[EOSN5][EOSN4][EOSN3][1][EOSN1];


// simple density-internal energy dependent EOS table
#if(WHICHEOS==KAZFULL)
double a_eossimpletable[NUMEOSQUANTITIESMEM][EOSSIMPLEN5][EOSSIMPLEN4][EOSSIMPLEN3][EOSSIMPLEN2][EOSSIMPLEN1];
#endif
static double (*eossimpletable)[EOSSIMPLEN5][EOSSIMPLEN4][EOSSIMPLEN3][EOSSIMPLEN2][EOSSIMPLEN1];

#if(WHICHEOS==KAZFULL)
double a_eosdegensimpletable[NUMEOSDEGENQUANTITIESMEM][EOSSIMPLEN5][EOSSIMPLEN4][EOSSIMPLEN3][1][EOSSIMPLEN1];
#endif
static double (*eosdegensimpletable)[EOSSIMPLEN5][EOSSIMPLEN4][EOSSIMPLEN3][1][EOSSIMPLEN1];


// simplezoom density-internal energy dependent EOS table
#if(WHICHEOS==KAZFULL)
double a_eossimplezoomtable[NUMEOSQUANTITIESMEM][EOSSIMPLEZOOMN5][EOSSIMPLEZOOMN4][EOSSIMPLEZOOMN3][EOSSIMPLEZOOMN2][EOSSIMPLEZOOMN1];
#endif
static double (*eossimplezoomtable)[EOSSIMPLEZOOMN5][EOSSIMPLEZOOMN4][EOSSIMPLEZOOMN3][EOSSIMPLEZOOMN2][EOSSIMPLEZOOMN1];

#if(WHICHEOS==KAZFULL)
double a_eosdegensimplezoomtable[NUMEOSDEGENQUANTITIESMEM][EOSSIMPLEZOOMN5][EOSSIMPLEZOOMN4][EOSSIMPLEZOOMN3][1][EOSSIMPLEZOOMN1];
#endif
static double (*eosdegensimplezoomtable)[EOSSIMPLEZOOMN5][EOSSIMPLEZOOMN4][EOSSIMPLEZOOMN3][1][EOSSIMPLEZOOMN1];




////////////////////////////////////
//
// Some global position variables used to determine EOS
//

// should be 4
#define NUMHDIRECTIONS 4


// NOTE: must be in same order and number as EOS independent vars
// GODMARK: must also change MAXPARLIST in nondepnmemonics.h
#define NUMNONSTANDARD 6
#define NUMEOSGLOBALS (NUMEXTRAINDEP+NUMNONSTANDARD)   // number of per CPU position-based data for EOS

// these should be ordered and numbered such that correspond to EOS table independent variables
// rho, u/p/chi, Y_e, Y_\nu, H
// 1    2        3    4      5 == 5D
// do NOT correspond to expanded independent variables list from EOS as read-in (i.e. not RHOEOS, UEOS, PEOS, CHIEOS, TEOS, YNUEOS,  HEOS)
#define FIRSTEOSGLOBAL 3 // 1=rho, 2=u/p/chi, 3=Y_e, 4=Y_\nu 5=H
#define TDYNORYEGLOBAL (FIRSTEOSGLOBAL)             // Tdyn or Y_e depending upon whichrnpmethod
#define YNUGLOBAL (TDYNORYEGLOBAL+1) // Tdyn or Y_\nu depending upon whichynumethod
#define HGLOBAL (YNUGLOBAL+1)        // scale-height (for method that uses this for EOS, some averaged version of H
#define H2GLOBAL (HGLOBAL+1)         // 2,3,4 are other directions for axisymmetric emission
#define H3GLOBAL (H2GLOBAL+1) 
#define H4GLOBAL (H3GLOBAL+1) 
#define UNUGLOBAL (H4GLOBAL+1)       // extra non-standard variable used to speed up iterative process when doing whichdatatype==4
#define PNUGLOBAL (UNUGLOBAL+1)      // extra non-standard variable used to speed up iterative process when doing whichdatatype==4
#define SNUGLOBAL (PNUGLOBAL+1)      // extra non-standard variable used to speed up iterative process when doing whichdatatype==4


// EOSglobal is per CPU data for EOS
#if(WHICHEOS==KAZFULL)
FTYPE a_EOSglobal[NUMEOSGLOBALS][N1M][N2M][N3M];
#endif
FTYPE (*EOSglobal)[N1M][N2M][N3M];


///////////////
//
// some globally accessed functions
//
static int get_eos_fromtable(int whichfun, FTYPE quant1, FTYPE quant2, FTYPE *answer);


// tolerance to check whether repeated case for i,j,k,rho0,u
#define OLDTOLERANCE (1E-14)

// tolerance for checks on input values of table
#define TABLETOL (1E-14)

// value of read-in temperature such that below this is treated as indicating an invalid (rho0,u) EOS pair
// actual read-in value is 1E-20, but using 5E-20 guarantees no machine-error choices and works with floats too
// also, generally is more accurate as temperature since problems with inversion are near T~0
#define INVALIDTEMP (5E-20)
// code value of invalid temperature and log10 version
FTYPE invalidtempcode,invalidlogtempcode;

int whichrnpmethod[NUMTBLS], whichynumethod[NUMTBLS], whichhcmmethod[NUMTBLS];
int whichdatatype[NUMTBLS],numc[NUMTBLS],numextras[NUMTBLS];

int primarytable=-1; // -1 indicates no EOS setup/read-in yet

int extralimits[MAXNUMDATATYPES][2];

#define EOSHEADNAME "eosnew.head"
#define EOSTABLENAME "eosnew.dat"
#define EOSDEGENTABLENAME "eosdegennew.dat"

#define EOSSIMPLEHEADNAME "eossimplenew.head"
#define EOSSIMPLETABLENAME "eossimplenew.dat"
#define EOSDEGENSIMPLETABLENAME "eosdegensimplenew.dat"

#define EOSSIMPLEZOOMHEADNAME "eossimplezoomnew.head"
#define EOSSIMPLEZOOMTABLENAME "eossimplezoomnew.dat"
#define EOSDEGENSIMPLEZOOMTABLENAME "eosdegensimplezoomnew.dat"



FTYPE NUCLEAROFFSET;
int didsetupkazeos;







// reads in table for EOS and sets up the table parameters
void read_setup_eostable(void)
{
  FILE *inhead;
  FILE *intable;
  FILE *indegentable;
  int	ii,jj;
  int m,n,o,p,q; // 5 dimension labels
  int iii,jjj,kkk,lll,mmm; // for 5 dimensions of table
  int ppp,qqq;
  int totalindex[NUMEOSINDEPS];
  FTYPE indep[NUMEOSINDEPS];
  FTYPE indepplusdegen[NUMEOSDEGENINDEPS];
  FTYPE indepdegen[NUMEOSDEGENINDEPS];
  FTYPE lstepdep,diff,lindep,lindeptry;
  char headername[NUMTBLS][MAXFILENAME];
  char tablename[NUMTBLS][MAXFILENAME];
  char degentablename[NUMTBLS][MAXFILENAME];
  int tableiter;
  int tablesizeexpected[NUMTBLS][NUMEOSINDEPS];
  double tabletemp[NUMEOSQUANTITIESMEM]; // same type as tables
  double degentabletemp[NUMEOSDEGENQUANTITIESMEM]; // same type as tables
  int sizematch; // 0 = table sizes don't match  1 = they match
  double valuetemp[NUMEOSQUANTITIESMEM],degenvaluetemp[NUMEOSDEGENQUANTITIESMEM];
  FTYPE errordegen,errordegen2;
  int i;
  int numeosquantitiestype[MAXNUMDATATYPES];






  trifprintf("Setting up Kaz EOS table\n");





  if(ALLOWFULLTABLE) primarytable=FULLTABLE;
  else if(ALLOWSIMPLETABLE) primarytable=SIMPLETABLE;
  else if(ALLOWSIMPLEZOOMTABLE) primarytable=SIMPLEZOOMTABLE;



  ///////////////////////////////////////////////
  //
  // Set expected sizes of full tables (not degen tables, for which N2=1)
  //
  ///////////////////////////////////////////////



  i=0;
  numeosquantitiestype[i]=NUMEOSQUANTITIESTYPE1; extralimits[i][0]=EXTRA1; extralimits[i][1]=DATATYPE1_EXTRAFINAL; i++;
  numeosquantitiestype[i]=NUMEOSQUANTITIESTYPE2; extralimits[i][0]=EXTRA1; extralimits[i][1]=DATATYPE2_EXTRAFINAL; i++;
  numeosquantitiestype[i]=NUMEOSQUANTITIESTYPE3; extralimits[i][0]=EXTRA1; extralimits[i][1]=DATATYPE3_EXTRAFINAL; i++;
  numeosquantitiestype[i]=NUMEOSQUANTITIESTYPE4; extralimits[i][0]=EXTRA1; extralimits[i][1]=DATATYPE4_EXTRAFINAL; i++;
  if(i!=MAXNUMDATATYPES){
    dualfprintf(fail_file,"numeosquantitiestype not setup for that many data types: i=%d\n",i);
    myexit(3206881);
  }


  // number of things corresponds to read-in number of quantities, not final table sizes
  i=0;
  tablesizeexpected[FULLTABLE][i]=EOSN1; i++;
  tablesizeexpected[FULLTABLE][i]=EOSN2; i++;
  tablesizeexpected[FULLTABLE][i]=EOSN2; i++;
  tablesizeexpected[FULLTABLE][i]=EOSN2; i++;
  tablesizeexpected[FULLTABLE][i]=EOSN3; i++;
  tablesizeexpected[FULLTABLE][i]=EOSN4; i++;
  tablesizeexpected[FULLTABLE][i]=EOSN5; i++;
  if(i!=NUMEOSINDEPS){
    dualfprintf(fail_file,"tablesizeexpected(full) not setup for that many indepdimens: i=%d\n",i);
    myexit(3206882);
  }

  i=0;
  tablesizeexpected[SIMPLETABLE][i]=EOSSIMPLEN1; i++;
  tablesizeexpected[SIMPLETABLE][i]=EOSSIMPLEN2; i++;
  tablesizeexpected[SIMPLETABLE][i]=EOSSIMPLEN2; i++;
  tablesizeexpected[SIMPLETABLE][i]=EOSSIMPLEN2; i++;
  tablesizeexpected[SIMPLETABLE][i]=EOSSIMPLEN3; i++;
  tablesizeexpected[SIMPLETABLE][i]=EOSSIMPLEN4; i++;
  tablesizeexpected[SIMPLETABLE][i]=EOSSIMPLEN5; i++;
  if(i!=NUMEOSINDEPS){
    dualfprintf(fail_file,"tablesizeexpected(simple) not setup for that many indepdimens: i=%d\n",i);
    myexit(3206883);
  }

  i=0;
  tablesizeexpected[SIMPLEZOOMTABLE][i]=EOSSIMPLEZOOMN1; i++;
  tablesizeexpected[SIMPLEZOOMTABLE][i]=EOSSIMPLEZOOMN2; i++;
  tablesizeexpected[SIMPLEZOOMTABLE][i]=EOSSIMPLEZOOMN2; i++;
  tablesizeexpected[SIMPLEZOOMTABLE][i]=EOSSIMPLEZOOMN2; i++;
  tablesizeexpected[SIMPLEZOOMTABLE][i]=EOSSIMPLEZOOMN3; i++;
  tablesizeexpected[SIMPLEZOOMTABLE][i]=EOSSIMPLEZOOMN4; i++;
  tablesizeexpected[SIMPLEZOOMTABLE][i]=EOSSIMPLEZOOMN5; i++;
  if(i!=NUMEOSINDEPS){
    dualfprintf(fail_file,"tablesizeexpected(simplezoom) not setup for that many indepdimens: i=%d\n",i);
    myexit(3206884);
  }


  // setup what q1-q5 mean associated with the 5 dimensions of the arrays
  i=1;
  vartypearray[i]=RHOEOS;  i++;
  vartypearray[i]=UEOS;    i++; // UEOS used for reading table, but later changed for each whichindep
  vartypearray[i]=TEOS;    i++;
  vartypearray[i]=YNUEOS;  i++;
  vartypearray[i]=HEOS;    i++;
  if(i!=NUMINDEPDIMENS+1){
    dualfprintf(fail_file,"vartypearray not setup for that many indepdimens: i=%d\n",i);
    myexit(3206885);
  }




  // setup what degen type q1-q"5" mean associated with the 5 dimensions of the arrays
  i=1;
  vardegentypearray[i]=RHOEOSDEGEN;  varnormalcompare2degentypearray[i] = RHOEOS; i++;
  vardegentypearray[i]=TEOSDEGEN;    varnormalcompare2degentypearray[i] = TEOS;   i++;
  vardegentypearray[i]=YNUEOSDEGEN;  varnormalcompare2degentypearray[i] = YNUEOS; i++;
  vardegentypearray[i]=HEOSDEGEN;    varnormalcompare2degentypearray[i] = HEOS;   i++;
  if(i!=NUMEOSDEGENINDEPS+1){
    dualfprintf(fail_file,"vardegentypearray and varnormalcompare2degentypearray not setup for that many indepdimens: i=%d\n",i);
    myexit(3206886);
  }

  


  ///////////////////////////////////////////////
  //
  // assign global pointer to eostable (static assignment)
  //
  //////////////////////////////////////////////
#if(WHICHEOS==KAZFULL)
  eostable = (double (*)[EOSN5][EOSN4][EOSN3][EOSN2][EOSN1])(&(a_eostable[0][0][0][0][0][0]));
  eosdegentable = (double (*)[EOSN5][EOSN4][EOSN3][1][EOSN1])(&(a_eosdegentable[0][0][0][0][0][0]));

  eossimpletable = (double (*)[EOSSIMPLEN5][EOSSIMPLEN4][EOSSIMPLEN3][EOSSIMPLEN2][EOSSIMPLEN1])(&(a_eossimpletable[0][0][0][0][0][0]));
  eosdegensimpletable = (double (*)[EOSSIMPLEN5][EOSSIMPLEN4][EOSSIMPLEN3][1][EOSSIMPLEN1])(&(a_eosdegensimpletable[0][0][0][0][0][0]));

  eossimplezoomtable = (double (*)[EOSSIMPLEZOOMN5][EOSSIMPLEZOOMN4][EOSSIMPLEZOOMN3][EOSSIMPLEZOOMN2][EOSSIMPLEZOOMN1])(&(a_eossimplezoomtable[0][0][0][0][0][0]));
  eosdegensimplezoomtable = (double (*)[EOSSIMPLEZOOMN5][EOSSIMPLEZOOMN4][EOSSIMPLEZOOMN3][1][EOSSIMPLEZOOMN1])(&(a_eosdegensimplezoomtable[0][0][0][0][0][0]));

  EOSglobal = (FTYPE (*)[N1M][N2M][N3M])(&(a_EOSglobal[-FIRSTEOSGLOBAL][N1BND][N2BND][N3BND])); // -FIRSTEOSGLOBAL so EOSglobal[FIRSTEOSGLOBAL] is a_EOSglobal[0]
#else
  eostable = (double (*)[EOSN5][EOSN4][EOSN3][EOSN2][EOSN1])(0);
  eosdegentable = (double (*)[EOSN5][EOSN4][EOSN3][1][EOSN1])(0);

  eossimpletable = (double (*)[EOSSIMPLEN5][EOSSIMPLEN4][EOSSIMPLEN3][EOSSIMPLEN2][EOSSIMPLEN1])(0);
  eosdegensimpletable = (double (*)[EOSSIMPLEN5][EOSSIMPLEN4][EOSSIMPLEN3][1][EOSSIMPLEN1])(0);

  eossimplezoomtable = (double (*)[EOSSIMPLEZOOMN5][EOSSIMPLEZOOMN4][EOSSIMPLEZOOMN3][EOSSIMPLEZOOMN2][EOSSIMPLEZOOMN1])(0);
  eosdegensimplezoomtable = (double (*)[EOSSIMPLEZOOMN5][EOSSIMPLEZOOMN4][EOSSIMPLEZOOMN3][1][EOSSIMPLEZOOMN1])(0);

  EOSglobal = (FTYPE (*)[N1M][N2M][N3M])(0);
#endif




  // GODMARK: should really read-in from table, although expected tabular value is 9.14Mev and this gives 7.108 for smooth connection for initial conditions
  // GODMARK: The below is for Shen EOS only
  NUCLEAROFFSET=7.57E-3;




  ///////////////////////////////////////////////
  //
  // Get tables from files
  //
  //////////////////////////////////////////////

  // do only over first CPU (myid==0) and then send data to rest of CPUs
  if(myid==0){


    ///////////////////////////////////////////////
    //
    // Set file names
    //
    ///////////////////////////////////////////////

    // header and table names
    strcpy(headername[FULLTABLE],EOSHEADNAME);
    strcpy(tablename[FULLTABLE],EOSTABLENAME);
    strcpy(degentablename[FULLTABLE],EOSDEGENTABLENAME);


    strcpy(headername[SIMPLETABLE],EOSSIMPLEHEADNAME);
    strcpy(tablename[SIMPLETABLE],EOSSIMPLETABLENAME);
    strcpy(degentablename[SIMPLETABLE],EOSDEGENSIMPLETABLENAME);

    strcpy(headername[SIMPLEZOOMTABLE],EOSSIMPLEZOOMHEADNAME);
    strcpy(tablename[SIMPLEZOOMTABLE],EOSSIMPLEZOOMTABLENAME);
    strcpy(degentablename[SIMPLEZOOMTABLE],EOSDEGENSIMPLEZOOMTABLENAME);





    ///////////////
    //
    // loop over normal tables (degen table read-in with its normal version)
    //
    //////////////

    for(tableiter=0;tableiter<NUMTBLS;tableiter++){

      // avoid tables not needed since should not expect user has copied or created such tables if not using them
      if(tableiter==0 && ALLOWFULLTABLE==0) continue;
      if(tableiter==1 && ALLOWSIMPLETABLE==0) continue;
      if(tableiter==2 && ALLOWSIMPLEZOOMTABLE==0) continue;



      ///////////////////////////////////
      //
      // first read-in the header file (only normal table has header file since degen table is same except one dimension (u/p/chi))
      //
      //////////////////////////////////
      
      // open header file
      if( (inhead = fopen(headername[tableiter],"rb"))==NULL){
	dualfprintf(fail_file,"No such file %s\n",headername[tableiter]);
	myexit(16622);
      }
      
      // get number of output colums and check that code and data file agree
      // also get whichmethod using
      fscanf(inhead,"%d %d %d",&whichrnpmethod[tableiter],&whichynumethod[tableiter],&whichhcmmethod[tableiter]);
      fscanf(inhead,"%d %d %d",&whichdatatype[tableiter],&numc[tableiter],&numextras[tableiter]);




      ////////////////////////////
      //
      // perform some checks
      //
      ////////////////////////////

      if(whichrnpmethod[tableiter]==0){
	fprintf(fail_file,"This method is not setup\n");
	myexit(3966738);
      }

      if(numextras[tableiter]>MAXNUMEXTRAS){
	dualfprintf(fail_file,"Increase MAXNUMEXTRAS to %d\n",numextras[tableiter]);
	myexit(626236);
      }

      if(whichdatatype[tableiter]>MAXNUMDATATYPES){
	dualfprintf(fail_file,"Increase MAXNUMDATATYPES to %d\n",whichdatatype[tableiter]);
	myexit(626237);
      }




      // check number of quantities
      numeosquantities[tableiter]=numeosquantitiestype[whichdatatype[tableiter]-1];
      numeosdegenquantities[tableiter]=NUMEOSDEGENQUANTITIESMEM;
      if(numc[tableiter]!=NUMINDEPDIMENS + numeosdegenquantities[tableiter] + NUMEOSINDEPS + numeosquantities[tableiter]){
	dualfprintf(fail_file,"numcolumns=%d but code expected %d\n",numc[tableiter],NUMINDEPDIMENS + numeosdegenquantities[tableiter] + NUMEOSINDEPS + numeosquantities[tableiter]);
	dualfprintf(fail_file,"tableiter=%d whichdatatype=%d :: %d %d %d %d\n",tableiter,whichdatatype[tableiter],NUMINDEPDIMENS,numeosdegenquantities[tableiter],NUMEOSINDEPS,numeosquantities[tableiter]);
	myexit(16623);
      }
      


      for(ii=0;ii<NUMEOSINDEPS;ii++){
	fscanf(inhead,"%d",&tablesize[tableiter][ii]);
	// second [4] : 0 = lower log_base limit, 1 = upper log_base limit, 2=step, 3 = divisor of grid position 4=base of log, 5 = linear value of offset for log_base stepping so can control how resolved
 	fscanf(inhead,HEADERONEIN,&inputtablelimits[tableiter][ii][0]);
	fscanf(inhead,HEADERONEIN,&inputtablelimits[tableiter][ii][1]);
	fscanf(inhead,HEADERONEIN,&inputtablelimits[tableiter][ii][2]);
	fscanf(inhead,HEADERONEIN,&inputtablelimits[tableiter][ii][4]); // base of log offset
	fscanf(inhead,HEADERONEIN,&inputtablelimits[tableiter][ii][5]); // linear offset


	// convert base 10 to actual base with actual offset so can be easily used for mapping
	// x_in:=
	tablelimits[tableiter][ii][0] = log10(pow(10.0,inputtablelimits[tableiter][ii][0])-inputtablelimits[tableiter][ii][5])/log10(inputtablelimits[tableiter][ii][4]);
	// x_out:=
	tablelimits[tableiter][ii][1] = log10(pow(10.0,inputtablelimits[tableiter][ii][1])-inputtablelimits[tableiter][ii][5])/log10(inputtablelimits[tableiter][ii][4]);
	
	tablelimits[tableiter][ii][4]=inputtablelimits[tableiter][ii][4]; // just copy base, assumed to always be log_base() where the base itself is in base 10 format
	tablelimits[tableiter][ii][5]=inputtablelimits[tableiter][ii][5]; // shouldn't need to use linear version of offset, but just copy for now
	


	if(fabs(inputtablelimits[tableiter][ii][2])<SMALL || tablesize[tableiter][ii]==1 ){
	  // then assume that direction has no dimensionality
	  inputtablelimits[tableiter][ii][2] = 1.0;

	  tablelimits[tableiter][ii][2] = 0.0; // forced
      
	  tablelimits[tableiter][ii][3] = 0.0; // forced so index=0 always

	}
	else{

	  // compute log_base step to be used in formula:
	  // r = r_0 + base^(x_in+dx*i) for i=0..N-1 where x = x_in+dx*i
	  // dx = (x_out-x_in)/(N-1)
	  // dx=:
	  tablelimits[tableiter][ii][2] = (tablelimits[tableiter][ii][1]-tablelimits[tableiter][ii][0])/((FTYPE)tablesize[tableiter][ii]-1.0);
      
	  // the below definition is consistent with Kaz's code, matlab eos_extract.m and elsewhere in this code
	  //	tablelimits[tableiter][ii][3] = ((FTYPE)tablesize[tableiter][ii]-1.0)/(tablelimits[tableiter][ii][1] - tablelimits[tableiter][ii][0]);


	  // i = (log_b (r-r_0) - x_in)/dx
	  // 1/dx:=
	  tablelimits[tableiter][ii][3] = 1.0/tablelimits[tableiter][ii][2];

	}
	
	// below used to truncate(limit) input values so lookup doesn't have as many conditionals
	// lineartablelimits does NOT include degen offset, so these are combined with offset when used since offset is function of rhob,hcm,tdynorye
	lineartablelimits[tableiter][ii][0]=pow(10.0,inputtablelimits[tableiter][ii][0]); // linear inner
	lineartablelimits[tableiter][ii][1]=pow(10.0,inputtablelimits[tableiter][ii][1]); // linear outer
	lineartablelimits[tableiter][ii][2]=1.0/log10(inputtablelimits[tableiter][ii][4]); // 1.0 divided by "normal log(base)" (used for mapping)
	lineartablelimits[tableiter][ii][3]=inputtablelimits[tableiter][ii][5]; // linear offset
      }

      // Temperature table limits not read in, just for user use


      // check size of tables
      sizematch=1; // assume match
      for(ii=0;ii<NUMEOSINDEPS;ii++){
	sizematch *= tablesizeexpected[tableiter][ii]==tablesize[tableiter][ii];
      }

      if(!sizematch){
	dualfprintf(fail_file,"Size of table does not match\n");
	for(ii=0;ii<NUMEOSINDEPS;ii++){
	  dualfprintf(fail_file,"tablesize[%d]=%d\n",ii,tablesize[tableiter][ii]);
	}
	dualfprintf(fail_file,"EOSN's = ");
	for(ii=0;ii<NUMEOSINDEPS;ii++){
	  dualfprintf(fail_file," %d",tablesizeexpected[tableiter][ii]);
	}
	dualfprintf(fail_file,"\n");
	myexit(16624);
      }
      
      fclose(inhead);
      






      ////////////////////////////
      //
      // at this point header is consistent with expectations, so read in the table
      //
      // now read-in table as written
      //
      //////////////////////////////////
      
      // open eostable file
      if( (intable = fopen(tablename[tableiter],"rb"))==NULL){
	dualfprintf(fail_file,"No such file %s\n",tablename[tableiter]);
	myexit(16625);
      }


      // open eosdegentable file
      if( (indegentable = fopen(degentablename[tableiter],"rb"))==NULL){
	dualfprintf(fail_file,"No such file %s\n",degentablename[tableiter]);
	myexit(166255);
      }


      // file has rhob as fastest index
      // here assumes tablesize of UEOS, PEOS, and CHIEOS are same
      // assume jjj=0 to start since degen checks below depend on that
      // notice that vartypearray has same size as dimension of arrays and loops.
      for(mmm=0;mmm<tablesize[tableiter][vartypearray[5]];mmm++)for(lll=0;lll<tablesize[tableiter][vartypearray[4]];lll++)for(kkk=0;kkk<tablesize[tableiter][vartypearray[3]];kkk++)for(jjj=0;jjj<tablesize[tableiter][vartypearray[2]];jjj++)for(iii=0;iii<tablesize[tableiter][vartypearray[1]];iii++){

	// first, get positions to make sure everything is consistent
	fscanf(intable,"%d %d %d %d %d",&m,&n,&o,&p,&q);

	if(m!=iii || n!=jjj || o!=kkk || p!=lll || q!=mmm){
	  dualfprintf(fail_file,"Read-in table (%d) indicies inconsistent with expected indicies: m=%d iii=%d n=%d jjj=%d o=%d kkk=%d p=%d lll=%d q=%d mmm=%d\n",tableiter,m,iii,n,jjj,o,kkk,p,lll,q,mmm);
	  dualfprintf(fail_file,"whichrnpmethod=%d whichynumethod=%d whichhcmmethod=%d\n",whichrnpmethod[tableiter],whichynumethod[tableiter],whichhcmmethod[tableiter]);
	  dualfprintf(fail_file,"whichdatatype=%d\n",whichdatatype[tableiter]);
	  for(jj=0;jj<NUMEOSINDEPS;jj++){
	    dualfprintf(fail_file,"tablesize[%d][%d]=%d\n",tableiter,jj,tablesize[tableiter][jj]);
	  }
	  myexit(16626);
	}

	if(jjj==0){ // degen table
	  n=jjj;
	  fscanf(indegentable,"%d %d %d %d",&m,&o,&p,&q);
	  if(m!=iii || n!=jjj || o!=kkk || p!=lll || q!=mmm){
	    dualfprintf(fail_file,"Read-in degentable indicies inconsistent with expected indicies: m=%d iii=%d n=%d jjj=%d o=%d kkk=%d p=%d lll=%d q=%d mmm=%d\n",tableiter,m,iii,n,jjj,o,kkk,p,lll,q,mmm);
	    myexit(166265);
	  }
	}





	// second, read in the independent variable values and compare with expected value
	for(ii=0;ii<NUMEOSINDEPS;ii++){
	  fscanf(intable,"%lf",&indep[ii]); // rhob, Udiff, Pdiff, CHIdiff, tdynorye, tdynorynu, H   for given grid value
	}
	// read-in UofUdiff, PofPdiff, CHIofCHIdiff -- used to check degen offset calculation in HARM
	for(ii=0;ii<numeosdegenquantities[tableiter];ii++){
	  fscanf(intable,"%lf",&indepplusdegen[ii]);
	}

	// true independent dimensions associated with free index
	totalindex[vartypearray[1]]   = iii;
	totalindex[CHIEOS] = totalindex[PEOS]=totalindex[UEOS]   = jjj; // notice that UEOS,PEOS,CHIEOS are actually same index
	totalindex[vartypearray[3]]   = kkk;
	totalindex[vartypearray[4]]   = lll;
	totalindex[vartypearray[5]]   = mmm;
	

	if(jjj==0){
	  for(ii=0;ii<NUMEOSDEGENINDEPS;ii++){
	    fscanf(indegentable,"%lf",&indepdegen[ii]); // rho, tdynorye, tdynorynu, H
	    // check consistency between normal and degen tablef or independent variables (assumes jjj!=0 in normal table is same for these quantities as jjj==0)
	  
	    if(fabs(indepdegen[vardegentypearray[ii]]-indep[varnormalcompare2degentypearray[ii]])>TABLETOL){
	      dualfprintf(fail_file,"degen table not consistent with normal table (tableiter=%d ii=%d) for %d %d: %21.15g %21.15g\n",tableiter,ii,vardegentypearray[ii],vartypearray[ii],indepdegen[vardegentypearray[ii]],indep[varnormalcompare2degentypearray[ii]]);
	    }
	  }
	}



	// third, check gridding of independent variables in normal table
	for(ii=0;ii<NUMEOSINDEPS;ii++){

	  if(tablesize[tableiter][ii]>1){ // only check actual table and assume degen table consistent

	    // get step (consistent with how step is computed in Kaz's code and in matlab script eos_extract.m)
	    // really only has to be consistent with eos_extract.m
	    //	    lstepdep = (-tablelimits[tableiter][ii][0])/((FTYPE)tablesize[tableiter][ii]-1.0);
	    lstepdep = inputtablelimits[tableiter][ii][2];
	    // compare step sizes to read-in step sizes
	    diff = fabs(lstepdep - tablelimits[tableiter][ii][2])/(fabs(lstepdep)+fabs(tablelimits[tableiter][ii][2]));
	    if(diff > TABLETOL){
	      dualfprintf(fail_file,"Grid step size is incorrect: mmm=%d lll=%d kkk=%d jjj=%d iii=%d :: ii=%d readin-value=%21.15g lstepdep=%21.15g\n",mmm,lll,kkk,jjj,iii,ii,tablelimits[tableiter][ii][2],lstepdep);
	      dualfprintf(fail_file,"tablelimits[%d][%d][0]=%21.15g tablelimits[%d][%d][1]=%21.15g\n",tableiter,ii,tablelimits[tableiter][ii][0],tableiter,ii,tablelimits[tableiter][ii][1]);
	      myexit(16627);
	    }


	    // grid used is:
	    //
	    // x = log(r-r_0)/log(base) such that r = r_0 + base^x
	    //
	    //
	    

	    // get read-in value of independent variable
	    // get x (here x is now such things as rhob, utotoffset, ptotoffset, chioffset, hcm, tdynorye)
	    //	    lindep=log10(indep[ii]);
	    lindep=log10(indep[ii]-lineartablelimits[tableiter][ii][3])*lineartablelimits[tableiter][ii][2];
	    // get computed value of independent variable (used later for lookup, so verifies lookup method)
	    // get x using lookup from tablular index (didn't need to change)
	    lindeptry=tablelimits[tableiter][ii][0] + totalindex[ii]*lstepdep;
	    // compare to be sure same
	    diff = fabs(lindep-lindeptry)/(fabs(lindep)+fabs(lindeptry));
	    if(diff>TABLETOL){
	      dualfprintf(fail_file,"Grid position data is incorrect: mmm=%d lll=%d kkk=%d jjj=%d iii=%d :: ii=%d readin-lindep=%21.15g lindeptry=%21.15g diff=%21.15g\n",mmm,lll,kkk,jjj,iii,ii,lindep,lindeptry,diff);
	      dualfprintf(fail_file,"tablelimits[%d][%d][0]=%21.15g totalindex[%d]=%d lstepdep=%21.15g\n",tableiter,ii,tablelimits[tableiter][ii][0],ii,totalindex[ii],lstepdep);
	      dualfprintf(fail_file,"%21.15g %21.15g %21.15g\n",log10(indep[ii]),lineartablelimits[tableiter][ii][3],lineartablelimits[tableiter][ii]);
	      myexit(16628);
	    }
	  }
	}



	// fourth, since everything is consistent, now read in columns of actual dependent variable values
	for(ppp=0;ppp<numeosquantities[tableiter];ppp++){
	  fscanf(intable,"%lf",&valuetemp[ppp]); // double values
	  if(tableiter==FULLTABLE){
	    eostable[ppp][mmm][lll][kkk][jjj][iii]=valuetemp[ppp];
	  }
	  else if(tableiter==SIMPLETABLE){
	    eossimpletable[ppp][mmm][lll][kkk][jjj][iii]=valuetemp[ppp];
	  }
	  else if(tableiter==SIMPLEZOOMTABLE){
	    eossimplezoomtable[ppp][mmm][lll][kkk][jjj][iii]=valuetemp[ppp];
	  }
	}



	if(jjj==0){ // degen table
	  for(ppp=0;ppp<numeosdegenquantities[tableiter];ppp++){
	    fscanf(indegentable,"%lf",&degenvaluetemp[ppp]); // double values
	    if(tableiter==FULLTABLE){
	      eosdegentable[ppp][mmm][lll][kkk][0][iii]=degenvaluetemp[ppp];
	    }
	    else if(tableiter==SIMPLETABLE){
	      eosdegensimpletable[ppp][mmm][lll][kkk][0][iii]=degenvaluetemp[ppp];
	    }
	    else if(tableiter==SIMPLEZOOMTABLE){
	      eosdegensimplezoomtable[ppp][mmm][lll][kkk][0][iii]=degenvaluetemp[ppp];
	    }
	  }
	}
	else{
	  // this is just for checks, store back into simple array the degen table
	  for(ppp=0;ppp<numeosdegenquantities[tableiter];ppp++){
	    if(tableiter==FULLTABLE){
	      degenvaluetemp[ppp]=eosdegentable[ppp][mmm][lll][kkk][0][iii];
	    }
	    else if(tableiter==SIMPLETABLE){
	      degenvaluetemp[ppp]=eosdegensimpletable[ppp][mmm][lll][kkk][0][iii];
	    }
	    else if(tableiter==SIMPLEZOOMTABLE){
	      degenvaluetemp[ppp]=eosdegensimplezoomtable[ppp][mmm][lll][kkk][0][iii];
	    }
	  }
	}



#if(0)
	// fifth, check degen table against normal table for proper use of independent variable offset
	// check for all jjj
	// now ensure that indep[normal table UEOS,PEOS,CHIEOS] + eosdegentable[UTOTDIFF,PTOTDIFF,CHIDIFF] = indepplusdegen[UTOTDIFF,PTOTDIFF,CHIDIFF] for degen'ed quantities
	// first and last are just read-in, while eosdegentable should be for all jjj, so use stored array for each tableiter
	for(ii=0;ii<numeosdegenquantities[tableiter];ii++){
	  // assume hit jjj==0 first time in loop
	  
	  errordegen=fabs(indep[ii+UEOS] + degenvaluetemp[ii] - indepplusdegen[ii])/fabs(fabs(indep[ii+UEOS]) + fabs(degenvaluetemp[ii]) + fabs(indepplusdegen[ii]));
	  errordegen2=fabs(indep[ii+UEOS])/(fabs(degenvaluetemp[ii] - indepplusdegen[ii])+SMALL);
	  if(errordegen>TABLETOL && errordegen2>0.5){
	    //	  if(errordegen2>1.0){
	    dualfprintf(fail_file,"Degen not correct: iii=%d jjj=%d kkk=%d lll=%d mmm=%d:: ii=%d :: error=%21.15g :: %21.15g %21.15g %21.15g :: %21.15g\n",iii,jjj,kkk,lll,mmm,ii,errordegen,indep[ii+UEOS],degenvaluetemp[ii],indepplusdegen[ii],errordegen2);
	  }

	}
#endif






	// continue onto next row
      }// end loop over all rows
      
      // done reading in table from file
      fclose(intable);
      fclose(indegentable);
      
      
      // report read-in
      trifprintf("Done reading in EOS table: tableiter=%d of %d\n",tableiter,NUMTBLS-1);







      //////////////////////////////////////////
      //
      // convert quantities to code units
      //
      //////////////////////////////////////////
    
      for(jj=0;jj<TBLLINEARITEMS;jj++){
	if(jj!=2){ // log of base has no units conversion, but rest do
	  if(rho0unittype==0) lineartablelimits[tableiter][RHOEOS][jj]/=rhounit;
	  else lineartablelimits[tableiter][RHOEOS][jj]/=rhomassunit;
	  lineartablelimits[tableiter][UEOS][jj]/=Pressureunit;
	  lineartablelimits[tableiter][PEOS][jj]/=Pressureunit;
	  lineartablelimits[tableiter][CHIEOS][jj]/=Pressureunit;
	  if(whichrnpmethod[tableiter]==0) lineartablelimits[tableiter][TEOS][jj]/=Tunit; // otherwise no conversion for dimensionless Y_e
	  if(whichynumethod[tableiter]==0) lineartablelimits[tableiter][YNUEOS][jj]/=Tunit; // otherwise no conversion for dimensionless Y_\nu
	  lineartablelimits[tableiter][HEOS][jj]/=Lunit;
	}
      }
      
      // recompute tablelimits (UPDOWN=[0,1])
      for(jj=0;jj<NUMEOSINDEPS;jj++) for(ii=0;ii<UPDOWN;ii++){
	// log_base(Rin-R0), log_base(Rout-R0), which gives same result as if shifted x_in and x_out by -log_base(r_units)
	// log(Rin-R0)/log(base) for both Rout and Rin
	tablelimits[tableiter][jj][ii]=log10(lineartablelimits[tableiter][jj][ii]-lineartablelimits[tableiter][jj][3])*lineartablelimits[tableiter][jj][2];
      }
      
      // convert log10 step to code units
      // unit division results in constant uniform offset in log-space
      // below ??[2] not ever used
      
      // for new generalized log grid, transformation is that \tilde{r} = r/r_units such that:
      // \tilde{x_in} = x_in - log_b(r_units)
      // new formula is just:
      // 
      // \tilde{r} = \tilde{r_0} + b^{\tilde{x_in} + i dx}
      // that is, dx and b are unchanged, there was only an offset of the starting position
      // shift x_in and x_out the same way is consistent with above recomputation of x_in and x_out

      // dx stays the same now [2]
      //      if(rho0unittype==0) tablelimits[tableiter][RHOEOS][2]=log10(pow(10.0,tablelimits[tableiter][RHOEOS][2])/rhounit);
      //      else tablelimits[tableiter][RHOEOS][2]=log10(pow(10.0,tablelimits[tableiter][RHOEOS][2])/rhomassunit);
      //      tablelimits[tableiter][UEOS][2]=log10(pow(10.0,tablelimits[tableiter][UEOS][2])/Pressureunit);
      //      tablelimits[tableiter][TEOS][2]=log10(pow(10.0,tablelimits[tableiter][TEOS][2])/Tunit);
      //      tablelimits[tableiter][YNUEOS][2]=log10(pow(10.0,tablelimits[tableiter][YNUEOS][2])/Tunit);
      //      tablelimits[tableiter][HEOS][2]=log10(pow(10.0,tablelimits[tableiter][HEOS][2])/Lunit);
      
      // just recompute the below
      for(ii=0;ii<NUMEOSINDEPS;ii++){
	tablelimits[tableiter][jj][3] = ((FTYPE)tablesize[tableiter][ii]-1.0)/(tablelimits[tableiter][ii][1] - tablelimits[tableiter][ii][0]);
      }

      // assume base stays same [4]
      // shouldn't need log10 version of offset, so assume not changed or needed anymore

      // now 0-3 (and 4-5) of tablelimits set, lineartablelimits set, tablesize is same
      
      
      // output to logfull_file so have information
      trifprintf("Table information in code units:\n");
      for(jj=0;jj<NUMEOSINDEPS;jj++){
	trifprintf("tablesize[%d][%d]=%d\n",tableiter,jj,tablesize[tableiter][jj]);
	for(ii=0;ii<TBLLINEARITEMS;ii++) trifprintf("lineartablelimits[%d][%d][%d]=%21.15g\n",tableiter,jj,ii,lineartablelimits[tableiter][jj][ii]);
	for(ii=0;ii<TBLITEMS;ii++) trifprintf("tablelimits[%d][%d][%d]=%21.15g\n",tableiter,jj,ii,tablelimits[tableiter][jj][ii]);
      }
      
      // report read-in
      trifprintf("Done converting units for header information: tableiter=%d of %d\n",tableiter,NUMTBLS-1);
      
          


      
      ///////////////////////////////
      //
      // convert eostable units
      //
      ///////////////////////////////


      // UEOS gives same size as PEOS and CHIEOS
      for(mmm=0;mmm<tablesize[tableiter][vartypearray[5]];mmm++)for(lll=0;lll<tablesize[tableiter][vartypearray[4]];lll++)for(kkk=0;kkk<tablesize[tableiter][vartypearray[3]];kkk++)for(jjj=0;jjj<tablesize[tableiter][vartypearray[2]];jjj++)for(iii=0;iii<tablesize[tableiter][vartypearray[1]];iii++){
	
	// temp store
	if(tableiter==FULLTABLE) for(jj=0;jj<numeosquantities[tableiter];jj++) tabletemp[jj]=eostable[jj][mmm][lll][kkk][jjj][iii];
	else if(tableiter==SIMPLETABLE) for(jj=0;jj<numeosquantities[tableiter];jj++) tabletemp[jj]=eossimpletable[jj][mmm][lll][kkk][jjj][iii];
	else if(tableiter==SIMPLEZOOMTABLE) for(jj=0;jj<numeosquantities[tableiter];jj++) tabletemp[jj]=eossimplezoomtable[jj][mmm][lll][kkk][jjj][iii];

	if(jjj==0){// degen table
	  if(tableiter==FULLTABLE) for(jj=0;jj<numeosdegenquantities[tableiter];jj++) degentabletemp[jj]=eosdegentable[jj][mmm][lll][kkk][jjj][iii];
	  else if(tableiter==SIMPLETABLE) for(jj=0;jj<numeosdegenquantities[tableiter];jj++) degentabletemp[jj]=eosdegensimpletable[jj][mmm][lll][kkk][jjj][iii];
	  else if(tableiter==SIMPLEZOOMTABLE) for(jj=0;jj<numeosdegenquantities[tableiter];jj++) degentabletemp[jj]=eosdegensimplezoomtable[jj][mmm][lll][kkk][jjj][iii];
	}


	tabletemp[PofRHOU]/=Pressureunit;
	tabletemp[UofRHOP]/=Pressureunit;
	// dPdRHO0ofRHOU dimensionless
	// dPdUofRHOU dimensionless
	tabletemp[CS2ofRHOU]/=(Vunit*Vunit);
	// entropy density (erg/K/cc)
	// kb doesn't convert units, but changes kb T to erg
	// presumes entropy is used with energy as in first-law: dQ = (kbT)dS where kbT is in ergs
	// previously had units of entropy density as erg/K/cc, but now units read-in are just 1/cc to avoid use of pointless kb
	tabletemp[SofRHOU]/=(1.0/pow(Lunit,3.0));
	// Note that often people will plot "entropy per baryon" that is "SofRHOU"/(\rho/m_b) that is dimensionless entropy per baryon
	// From HARM quantities, convert back to cgs and then compute above
	// below is (1/cc) / (erg/cc) \propto 1/erg since we have rho as rho c^2
	tabletemp[DSDRHOofRHOU]/=(1.0/energyunit);
	tabletemp[DSDUofRHOU]/=(1.0/energyunit);
	
	tabletemp[PofRHOCHI]/=Pressureunit;
	// IDRHO0DP is dimensionless
	// IDCHIDP is dimensionless
	
	// TEMP used for diagnostics, not used otherwise
	tabletemp[TEMPU]/=Tempunit;
	tabletemp[TEMPP]/=Tempunit;
	tabletemp[TEMPCHI]/=Tempunit;


	////////////////////////////////
	//
	// deal with extra quantities
	if(whichdatatype[tableiter]==1){
	  // Qm is in erg/s/cm^3 (Qvol, not Qsurf)
	  // this is divided by H when used as a volume rate
	  tabletemp[EXTRA1]/=(edotunit/(Lunit*Lunit*Lunit));
	}
	else if(whichdatatype[tableiter]==2){
	  // \tau/H
	  tabletemp[EXTRA1]/=(1.0/Lunit);
	  tabletemp[EXTRA2]/=(1.0/Lunit);
	  tabletemp[EXTRA3]/=(1.0/Lunit);
	  tabletemp[EXTRA4]/=(1.0/Lunit);
	  tabletemp[EXTRA5]/=(1.0/Lunit);
	  tabletemp[EXTRA6]/=(1.0/Lunit);
	  tabletemp[EXTRA7]/=(1.0/Lunit);
	  tabletemp[EXTRA8]/=(1.0/Lunit);
	  tabletemp[EXTRA9]/=(1.0/Lunit);
	  tabletemp[EXTRA10]/=(1.0/Lunit);
	  tabletemp[EXTRA11]/=(1.0/Lunit);
	  tabletemp[EXTRA12]/=(1.0/Lunit);
	  //	  \Gamma = 1/s
	  tabletemp[EXTRA13]/=(1.0/Tunit);
	  tabletemp[EXTRA14]/=(1.0/Tunit);
	  tabletemp[EXTRA15]/=(1.0/Tunit);
	  tabletemp[EXTRA16]/=(1.0/Tunit);
	}
	else if(whichdatatype[tableiter]==3){
	  // Qphoton=erg/s/cm^3
	  tabletemp[EXTRA1]/=(edotunit/(Lunit*Lunit*Lunit));

	  // Qm=erg/s/cm^3
	  tabletemp[EXTRA2]/=(edotunit/(Lunit*Lunit*Lunit));
	  
	  // graddotrhouyl=\rho/sec = m_b/s/cm^3
	  // \nabla_\mu (\rho_0 u^\mu Y_e) =  (m_b/H) (\dot{N}_{\bar{\nu}_e} - \dot{N}_{\bar{\nu}_e}) 
	  if(rho0unittype==0) tabletemp[EXTRA3]/=(rhounit/Tunit);
	  else tabletemp[EXTRA3]/=(rhomassunit/Tunit);
	  
	  // Tthermaltot
	  tabletemp[EXTRA4]/=Tunit;
	  tabletemp[EXTRA5]/=Tunit;

	  // lambdatot = mean free path such that H = [\int (dr/\lambda)] / [1/\lambda]
	  // and \tau = \int dr/\lambda
	  tabletemp[EXTRA6]/=Lunit;
	  tabletemp[EXTRA7]/=Lunit;

	  // Enuglobal : erg
	  tabletemp[EXTRA8]/=energyunit;
	  tabletemp[EXTRA9]/=energyunit;
	  tabletemp[EXTRA10]/=energyunit;

	  // no conversion for Ynuthermal that's dimensionless EXTRA11
	}
	else if(whichdatatype[tableiter]==4){
	  // \tau/H
	  tabletemp[EXTRA1]/=(1.0/Lunit);
	  tabletemp[EXTRA2]/=(1.0/Lunit);
	  tabletemp[EXTRA3]/=(1.0/Lunit);
	  tabletemp[EXTRA4]/=(1.0/Lunit);
	  tabletemp[EXTRA5]/=(1.0/Lunit);
	  tabletemp[EXTRA6]/=(1.0/Lunit);
	  tabletemp[EXTRA7]/=(1.0/Lunit);
	  tabletemp[EXTRA8]/=(1.0/Lunit);
	  tabletemp[EXTRA9]/=(1.0/Lunit);
	  tabletemp[EXTRA10]/=(1.0/Lunit);
	  tabletemp[EXTRA11]/=(1.0/Lunit);
	  tabletemp[EXTRA12]/=(1.0/Lunit);
	  // energy densities
	  tabletemp[EXTRA13]/=(Pressureunit);
	  tabletemp[EXTRA14]/=(Pressureunit);
	  tabletemp[EXTRA15]/=(Pressureunit);
	  // number densities
	  tabletemp[EXTRA16]/=(1.0/pow(Lunit,3.0));
	  tabletemp[EXTRA17]/=(1.0/pow(Lunit,3.0));
	  tabletemp[EXTRA18]/=(1.0/pow(Lunit,3.0));
	  // mean free paths (length)
	  tabletemp[EXTRA19]/=(Lunit);
	  tabletemp[EXTRA20]/=(Lunit);
	  // \tau/H for photons
	  tabletemp[EXTRA21]/=(1.0/Lunit);
	  tabletemp[EXTRA22]/=(1.0/Lunit);

	  // thermal number densities
	  tabletemp[EXTRA23]/=(1.0/pow(Lunit,3.0));
	  tabletemp[EXTRA24]/=(1.0/pow(Lunit,3.0));
	}


	
	//	for(jj=0;jj<numeosquantities[tableiter];jj++) dualfprintf(fail_file,"tableiter=%d tabletemp[%d][%d][%d][%d][%d][%d]=%21.15g\n",tableiter,jj,mmm,lll,kkk,jjj,iii,tabletemp[jj]);

	if(jjj==0){// degen table
	  for(jj=0;jj<numeosdegenquantities[tableiter];jj++){
	    degentabletemp[jj]/=Pressureunit; // all are pressure units
	  }
	}
	

	// temp restore
	if(tableiter==FULLTABLE) for(jj=0;jj<numeosquantities[tableiter];jj++) eostable[jj][mmm][lll][kkk][jjj][iii]=tabletemp[jj];
	else if(tableiter==SIMPLETABLE) for(jj=0;jj<numeosquantities[tableiter];jj++) eossimpletable[jj][mmm][lll][kkk][jjj][iii]=tabletemp[jj];
	else if(tableiter==SIMPLEZOOMTABLE) for(jj=0;jj<numeosquantities[tableiter];jj++) eossimplezoomtable[jj][mmm][lll][kkk][jjj][iii]=tabletemp[jj];

	// DEBUG:
	//	if(tableiter==FULLTABLE) for(jj=0;jj<numeosquantities[tableiter];jj++) dualfprintf(fail_file,"jj=%d eosvalue=%21.15g\n",jj,eostable[jj][mmm][lll][kkk][jjj][iii]);


	if(jjj==0){// degen table
	  // temp restore
	  if(tableiter==FULLTABLE) for(jj=0;jj<numeosdegenquantities[tableiter];jj++) eosdegentable[jj][mmm][lll][kkk][jjj][iii]=degentabletemp[jj];
	  else if(tableiter==SIMPLETABLE) for(jj=0;jj<numeosdegenquantities[tableiter];jj++) eosdegensimpletable[jj][mmm][lll][kkk][jjj][iii]=degentabletemp[jj];
	  else if(tableiter==SIMPLEZOOMTABLE) for(jj=0;jj<numeosdegenquantities[tableiter];jj++) eosdegensimplezoomtable[jj][mmm][lll][kkk][jjj][iii]=degentabletemp[jj];
	}


      } // end loop over table

            
      // report conversion
      trifprintf("Done converting units of EOS table: tableiter=%d of %d\n",tableiter,NUMTBLS-1);

    }// end loop over tables


    // compute code version of invalid temperature
    invalidtempcode=INVALIDTEMP/Tempunit;


    trifprintf("Sden_convertfactor=%21.15g\n",(1.0/pow(Lunit,3.0)));

  } // end if myid==0
  




  //////////////////////////////////////
  //
  // send data to all CPUs
  // so all CPUs have parameters and full eostable data
  //
  ///////////////////////////////////////
#if(USEMPI)

  MPI_Bcast(&numeosquantities,1,MPI_INT,NUMTBLS,MPI_COMM_GRMHD);
  MPI_Bcast(&numeosdegenquantities,1,MPI_INT,NUMTBLS,MPI_COMM_GRMHD);

  MPI_Bcast(&numeosquantitiestype,1,MPI_INT,MAXNUMDATATYPES,MPI_COMM_GRMHD);

  MPI_Bcast(&inputtablelimits[0][0][0],1,MPI_FTYPE,NUMTBLS*NUMEOSINDEPS*TBLITEMS,MPI_COMM_GRMHD);
  MPI_Bcast(&tablelimits[0][0][0],1,MPI_FTYPE,NUMTBLS*NUMEOSINDEPS*TBLITEMS,MPI_COMM_GRMHD);
  MPI_Bcast(&lineartablelimits[0][0][0],1,MPI_FTYPE,NUMTBLS*NUMEOSINDEPS*TBLLINEARITEMS,MPI_COMM_GRMHD);
  MPI_Bcast(&tablesize[0][0],1,MPI_INT,NUMTBLS*NUMEOSINDEPS,MPI_COMM_GRMHD);

  MPI_Bcast(&vartypearray[0],1,MPI_INT,NUMEOSINDEPS+1,MPI_COMM_GRMHD);
  MPI_Bcast(&vardegentypearray[0],1,MPI_INT,NUMEOSDEGENINDEPS+1,MPI_COMM_GRMHD);
  MPI_Bcast(&varnormalcompare2degentypearray[0],1,MPI_INT,NUMEOSDEGENINDEPS+1,MPI_COMM_GRMHD);

  
  // eostable is double, so MPI_DOUBLE
  MPI_Bcast(&(eostable[0][0][0][0][0][0]),numeosquantities[FULLTABLE]*EOSN5*EOSN4*EOSN3*EOSN2*EOSN1,MPI_DOUBLE,0,MPI_COMM_GRMHD);
  MPI_Bcast(&(eosdegentable[0][0][0][0][0][0]),numeosdegenquantities[FULLTABLE]*EOSN5*EOSN4*EOSN3*1*EOSN1,MPI_DOUBLE,0,MPI_COMM_GRMHD);
  
  // eossimpletable is double, so MPI_DOUBLE
  MPI_Bcast(&(eossimpletable[0][0][0][0][0][0]),numeosquantities[SIMPLETABLE]*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*EOSSIMPLEN2*EOSSIMPLEN1,MPI_DOUBLE,0,MPI_COMM_GRMHD);
  MPI_Bcast(&(eosdegensimpletable[0][0][0][0][0][0]),numeosdegenquantities[SIMPLETABLE]*EOSSIMPLEN5*EOSSIMPLEN4*EOSSIMPLEN3*1*EOSSIMPLEN1,MPI_DOUBLE,0,MPI_COMM_GRMHD);

  // eossimplezoomtable is double, so MPI_DOUBLE
  MPI_Bcast(&(eossimplezoomtable[0][0][0][0][0][0]),numeosquantities[SIMPLEZOOMTABLE]*EOSSIMPLEZOOMN5*EOSSIMPLEZOOMN4*EOSSIMPLEZOOMN3*EOSSIMPLEZOOMN2*EOSSIMPLEZOOMN1,MPI_DOUBLE,0,MPI_COMM_GRMHD);
  MPI_Bcast(&(eosdegensimplezoomtable[0][0][0][0][0][0]),numeosquantities[SIMPLEZOOMTABLE]*EOSSIMPLEZOOMN5*EOSSIMPLEZOOMN4*EOSSIMPLEZOOMN3*1*EOSSIMPLEZOOMN1,MPI_DOUBLE,0,MPI_COMM_GRMHD);

  MPI_Bcast(&invalidtempcode,1,MPI_FTYPE,1,MPI_COMM_GRMHD);
  MPI_Bcast(&invalidlogtempcode,1,MPI_FTYPE,1,MPI_COMM_GRMHD);

  MPI_Bcast(&whichrnpmethod[0],1,MPI_INT,NUMTBLS,MPI_COMM_GRMHD);
  MPI_Bcast(&whichynumethod[0],1,MPI_INT,NUMTBLS,MPI_COMM_GRMHD);
  MPI_Bcast(&whichhcmmethod[0],1,MPI_INT,NUMTBLS,MPI_COMM_GRMHD);

  MPI_Bcast(&whichdatatype[0],1,MPI_INT,NUMTBLS,MPI_COMM_GRMHD);
  MPI_Bcast(&numc[0],1,MPI_INT,NUMTBLS,MPI_COMM_GRMHD);
  MPI_Bcast(&numextras[0],1,MPI_INT,NUMTBLS,MPI_COMM_GRMHD);

  MPI_Bcast(&primarytable,1,MPI_INT,1,MPI_COMM_GRMHD);

  MPI_Bcast(&extralimits[0][0],1,MPI_INT,MAXNUMDATATYPES*2,MPI_COMM_GRMHD);

#endif
  

  didsetupkazeos=1;
  trifprintf("Done with reading in EOS tables\n");
  dualfprintf(log_file,"proc: %d: Done reading in EOS table\n",myid);
  
  
  
  
}



// iswithin_eostable() is written for each type of dataset(s) to used -- too complicated to make general

// check if density and pressure-like quantity are within table
// all quantities are real code units (i.e. linear mapping instead of log10 mapping)
// here we only presume to care about density and internal energy, while H and T are presumed to be truncated rather than extended with some alternative
// ifdegencheck: 0 = normal table check  1 = ignores q2 (u,p,chi) since generally q2 is largest range possible and later will restrict/check if within the full table
// whichindep = which independent variable (based upon which function looking up)
int iswithin_eostable(int ifdegencheck, int whichindep, int *vartypearray, FTYPE *qarray, int *whichtable)
{
  // don't assume tables are setup so q1 and q2 always have same range
  // assume all quantities depend on limits in same way, so "which" doesn't matter
  // assume HEOS and TEOS are always withing range for now



  /////////////////////////////
  //
  // Note that q3,q4 ((TDYN or YE),(TDYN or YNU), Hcm) are computed such that is forced to be within FULLTABLE limits so consistently used
  //
  /////////////////////////////

  if(ALLOWFULLTABLE
	  &&
	  qarray[1]>=lineartablelimits[FULLTABLE][vartypearray[1]][0] && qarray[1]<=lineartablelimits[FULLTABLE][vartypearray[1]][1]
	  &&
	  (ifdegencheck || qarray[2]>=lineartablelimits[FULLTABLE][vartypearray[2]][0] && qarray[2]<=lineartablelimits[FULLTABLE][vartypearray[2]][1])
	  ){
    *whichtable=FULLTABLE;
    return(1);
  }
  else if(ALLOWSIMPLETABLE
	  &&
	  qarray[1]>=lineartablelimits[SIMPLETABLE][vartypearray[1]][0] && qarray[1]<=lineartablelimits[SIMPLETABLE][vartypearray[1]][1]
	  &&
	  (ifdegencheck || qarray[2]>=lineartablelimits[SIMPLETABLE][vartypearray[2]][0] && qarray[2]<=lineartablelimits[SIMPLETABLE][vartypearray[2]][1])
     ){
    *whichtable=SIMPLETABLE;
    return(1);
  }
  else{

    if(debugfail>=2){ // DEBUG GODMARK: was turned on when debugging EOS
      dualfprintf(fail_file,"NOT IN LOOKUP: ifdegencheck=%d whichindep=%d qarray[1]=%21.15g qarray[2]=%21.15g :: i=%d j=%d k=%d p=%d\n",ifdegencheck,whichindep,qarray[1],qarray[2],icurr,jcurr,kcurr,pcurr);
      dualfprintf(fail_file,"lin0=%g lin1=%g\n",lineartablelimits[SIMPLETABLE][RHOEOS][0],lineartablelimits[SIMPLETABLE][RHOEOS][1]);
    }
    return(0);
  }
  // GODMARK:
  // alternative to this function is that we simply place a floor on the linear values of rho, u, H, and T so always within table values
  // Use: lineartablelimits[ii][0,1]
  //
  // problem is that for inversion, truncation still will be an if check so as slow as above anyways, so can just truncate here if wanted and then more general.  Only slows down simpler single-type calls to EOS.
}










// rho0,u,H,T have independent lookups for mapping to i,j,k,l and already done in eos_lookup_degen()
// u here stands for utotdiff,ptotdiff,and chidiff depending upon whichindep
// here only need to look up jeos associated with u
// notice that qarray and vartypearray start at 1 not 0
void eos_lookup_degen(int begin, int end, int skip, int whichtable, int whichindep, int *vartypearray, FTYPE *qarray, FTYPE *indexarray)
{
  FTYPE logq[NUMINDEPDIMENS+1];
  FTYPE prelogq[NUMINDEPDIMENS+1];
  int qi;

  // is this expensive?
  // log_base(R-R0)
  // avoid nan's by choosing minimal density

  // old simple log10 way:
  // i = (x - x0)*[(N-1)/(x1-x0)] such that if x=x0, then i=0, if x=x1, then i=N-1
  
  // new generalized log way:
  // i = (log_base(R-R0) - x_in)/dx

  for(qi=begin;qi<=end;qi++){
    if(qi==skip) continue;
    
    prelogq[qi] = MAX(qarray[qi]    -lineartablelimits[whichtable][vartypearray[qi]]  [3],SMALL);
    logq[qi]    = log10(prelogq[qi])*lineartablelimits[whichtable][vartypearray[qi]][2];

    indexarray[qi] = (logq[qi]   -tablelimits[whichtable][vartypearray[qi]]  [0])*tablelimits[whichtable][vartypearray[qi]]  [3];

    // DEBUG:
    dualfprintf(fail_file,"qi=%d qarray=%21.15g prelogq=%21.15g logq=%21.15g indexarray=%21.15g\n",qi,qarray[qi],prelogq[qi],logq[qi],indexarray[qi]);

  }


}


// whichdegen: 1 = do degen lookup  0 = do normal lookup without energy indepenent variable lookup
void eos_lookup_prepost_degen(int whichdegen, int whichtable, int whichindep, int *vartypearray, FTYPE *qarray, FTYPE *indexarray)
{

  void eos_lookup_degen(int begin, int end, int skip, int whichtable, int whichindep, int *vartypearray, FTYPE *qarray, FTYPE *indexarray);
  int skip;


  if(whichdegen==0){// nonsense value to avoid skipping
    // if normal table after degen table, then only need to get q2 location
    skip = -1000;
    eos_lookup_degen(2,2,skip, whichtable, whichindep, vartypearray, qarray, indexarray);
  }
  else{
    // if doing degentable, lookall up except q2=UEOS/PEOS/CHIEOS
    skip=2; 
    eos_lookup_degen(1,WHICHEOSDIMEN,skip,whichtable, whichindep, vartypearray, qarray, indexarray);
    indexarray[2]=0.0; // set 2 as 0 since degentable has this as 0 (may be floating -> integer issue in how this floating value is used)
  }

}


// quad-linear interpolation
FTYPE get_eos_fromlookup(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichfun, int whichindep, int *vartypearray, FTYPE *indexarray)
{
  FTYPE get_eos_fromlookup_nearest_dumb(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichfun, int whichindep, int *vartypearray, FTYPE *indexarray);
  FTYPE get_eos_fromlookup_nearest(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichfun, int whichindep, int *vartypearray, FTYPE *indexarray);
  FTYPE get_eos_fromlookup_linear(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichfun, int whichindep, int *vartypearray, FTYPE *indexarray);
  FTYPE get_eos_fromlookup_parabolic(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichfun, int whichindep, int *vartypearray, FTYPE *indexarray);
  FTYPE get_eos_fromlookup_parabolicfull(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichfun, int whichindep, int *vartypearray, FTYPE *indexarray);



#if(0)
  // neartest dumb for normal table
  if(degentable==0){
    return(get_eos_fromlookup_nearest_dumb(repeatedeos,tabledimen,degentable, whichtable, whichfun, whichindep, vartypearray, indexarray));
  }
  else{
    // linear
    return(get_eos_fromlookup_linear(repeatedeos,tabledimen,degentable, whichtable, whichfun, whichindep, vartypearray, indexarray));
  }
#elif(0)
  // nearest dumb for all
  return(get_eos_fromlookup_nearest_dumb(repeatedeos,tabledimen,degentable, whichtable, whichfun, whichindep, vartypearray, indexarray));
#elif(1)
  // linear
  return(get_eos_fromlookup_linear(repeatedeos,tabledimen,degentable, whichtable, whichfun, whichindep, vartypearray, indexarray));
#elif(0)
  // parabolic for density and tri-linear otherwise
  return(get_eos_fromlookup_parabolic(repeatedeos,tabledimen,degentable, whichtable, whichfun, whichindep, vartypearray, indexarray));
#elif(0)
  // parabolic for all quantities
  return(get_eos_fromlookup_parabolicfull(repeatedeos,tabledimen,degentable, whichtable, whichfun, whichindep, vartypearray, indexarray));
#endif

}







// tri-linear + parabolic (density) interpolation
FTYPE get_eos_fromlookup_parabolicfull(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichfun, int whichindep, int *vartypearray, FTYPE *indexarray)
{
  static int ii,jj,kk,ll,iio,jjo,kko,llo;
  static int startiii,startjjj,startkkk,startlll;
  static int endiii,endjjj,endkkk,endlll;
  FTYPE tempcheck;
  FTYPE totalf[3][3][3][3]; // 3 values for parabolic interpolation
  FTYPE totalffinal;
  FTYPE (*tfptr)[3][3][3];
  FTYPE *totalfptr;
  int iii,jjj,kkk,lll;
  int whichtemp,whichdegenfun;
  int whichinterp1,whichinterp2,loginterp;
  FTYPE xmx0,AA,BB;
  FTYPE ieos,jeos,keos,leos,meos;
  int EXTRASTART,EXTRAFINISH;


  ieos=indexarray[1];
  jeos=indexarray[2];
  keos=indexarray[3];
  leos=indexarray[4];




  tfptr=(FTYPE (*)[3][3][3]) (&(totalf[1][1][1][1])); // so tfptr[-1,0,1]

  // definition consistent with numerical assignments of indecies of arrays
  whichdegenfun = whichindep-1;


  // determine which temperature to use to check inversion
  // GODMARK: this could be stored in array where index is whichindep and output is whichtemp
  if(whichindep==UEOS) whichtemp=TEMPU;
  else if(whichindep==PEOS) whichtemp=TEMPP;
  else if(whichindep==CHIEOS) whichtemp=TEMPCHI;





#if(DOLOGINTERP)

  EXTRASTART=extralimits[whichdatatype[whichtable]-1][0];
  EXTRAFINISH=extralimits[whichdatatype[whichtable]-1][1];

  // GODMARK: Can make array that stores this info, looked up by whichfun as index
  // functions (F) F(rho0,u)
  whichinterp1=(whichfun==PofRHOCHI||whichfun==UofRHOP||whichfun==TEMPP||whichfun==PofRHOU||whichfun==CS2ofRHOU||whichfun==SofRHOU||(whichfun>=EXTRASTART && whichfun<=EXTRAFINISH)||whichfun==TEMPU||whichfun==TEMPCHI);
  // functions (F) F(rho0,p)
  whichinterp2=(whichfun==DPDRHOofRHOU||whichfun==DPDUofRHOU||whichfun==DSDRHOofRHOU||whichfun==DSDUofRHOU||whichfun==IDRHO0DP||whichfun==IDCHIDP);

  if(whichinterp1||degentable==1) loginterp=1;
  else if(whichinterp2) loginterp=0;
  else{
    dualfprintf(fail_file,"Undefined whichfun=%d in get_eos_fromlookup_linear()\n",whichfun);
    myexit(62662);
  }
#endif



  if(repeatedeos==0){
    // assume jeos, keos, leos set correctly for given degentable and tabledimen (i.e. are 0 or integers in those cases so d's computed correctly)
    ii=ROUND2INT(ieos); //ii=(int)ieos; // round instead when doing parabolic interpolation with 3 points
    // limit to within table
    if(ii<1) ii=1; // 0 is minimum
    if(ii>tablesize[whichtable][RHOEOS]-2) ii=tablesize[whichtable][RHOEOS]-2; // N-1 is maximum

    jj=ROUND2INT(jeos); //jj=(int)jeos; // round instead when doing parabolic interpolation with 3 points
    // limit to within table
    if(jj<1) jj=1; // 0 is minimum
    if(jj>tablesize[whichtable][whichindep]-2) jj=tablesize[whichtable][whichindep]-2; // N-1 is maximum

    kk=ROUND2INT(keos); //kk=(int)keos; // round instead when doing parabolic interpolation with 3 points
    // limit to within table
    if(kk<1) kk=1; // 0 is minimum
    if(kk>tablesize[whichtable][TEOS]-2) kk=tablesize[whichtable][TEOS]-2; // N-1 is maximum

    ll=ROUND2INT(leos); //ll=(int)leos; // round instead when doing parabolic interpolation with 3 points
    // limit to within table
    if(ll<1) ll=1; // 0 is minimum
    if(ll>tablesize[whichtable][YNUEOS]-2) ll=tablesize[whichtable][YNUEOS]-2; // N-1 is maximum

    // if only choosing 1 value, then choose rounded version
    iio=ROUND2INT(ieos);
    jjo=ROUND2INT(jeos);
    kko=ROUND2INT(keos);
    llo=ROUND2INT(leos);
    

    // set range of loops for different table types
    // tabledimen overrules table type (i.e. take section out of fuller table -- assumed to be 0 index)
    // GODMARK: these things couuld be stored as functions of whichtable/degentable/tabledimen
    if(tablesize[whichtable][RHOEOS]!=1) { startiii=-1; endiii=1;}
    else{ startiii=endiii=0; }

    if(degentable==0 || tablesize[whichtable][whichindep]!=1){ startjjj=-1; endjjj=1;}
    else {startjjj=endjjj=0;}

    if(tablesize[whichtable][TEOS]!=1){ startkkk=-1; endkkk=1;}
    else {startkkk=endkkk=0;}
     
    if(tablesize[whichtable][YNUEOS]!=1){ startlll=-1; endlll=1;}
    else {startlll=endlll=0;}
  }


    
  // Loop over nearby table values and determine bi-linearly interpolated value
  // 4-D means 2^4=16 positions
  // 2-D means 2^2=4 positions
  // get 3 values as function of density
  for(iii=startiii;iii<=endiii;iii++) for(jjj=startjjj;jjj<=endjjj;jjj++)for(kkk=startkkk;kkk<=endkkk;kkk++)for(lll=startlll;lll<=endlll;lll++){


#if(CHECKIFVALIDEOSDATA)
    if(degentable==0){
      // don't use values of table that have no inversion to temperature
      if(whichtable==SIMPLEZOOMTABLE) tempcheck=eossimplezoomtable[whichtemp][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
      else if(whichtable==SIMPLETABLE) tempcheck=eossimpletable[whichtemp][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
      else if(whichtable==FULLTABLE) tempcheck=eostable[whichtemp][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
    }
    //    dualfprintf(fail_file,"temp=%21.15g\n",tempcheck);
    if(degentable==1 || tempcheck>invalidtempcode) // Avoid invalid inversions if T>Tbad, but only deal with temperature if degentable==0
#else
    if(1)
#endif
    {
      if(whichtable==SIMPLEZOOMTABLE){
	if(degentable==0) tfptr[iii][jjj][kkk][lll] = eossimplezoomtable[whichfun][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
	else tfptr[iii][jjj][kkk][lll] = eosdegensimplezoomtable[whichdegenfun][0][ll+lll][kk+kkk][0][ii+iii];
      }
      else if(whichtable==SIMPLETABLE){
	if(degentable==0) tfptr[iii][jjj][kkk][lll] = eossimpletable[whichfun][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
	else tfptr[iii][jjj][kkk][lll] = eosdegensimpletable[whichdegenfun][0][ll+lll][kk+kkk][0][ii+iii];
      }
      else if(whichtable==FULLTABLE){
	if(degentable==0) tfptr[iii][jjj][kkk][lll] = eostable[whichfun][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
	else tfptr[iii][jjj][kkk][lll] = eosdegentable[whichdegenfun][0][ll+lll][kk+kkk][0][ii+iii];

      }


      // DEBUG
      //      dualfprintf(fail_file,"tabledimen=%d degentable=%d whichtable=%d whichfun=%d whichdegenfun=%d ii=%d iii=%d jj=%d jjj=%d kk=%d kkk=%d ll=%d lll=%d :: f=%21.15g dist=%21.15g totalf=%21.15g\n",tabledimen, degentable, whichtable,whichfun,whichdegenfun,ii,iii,jj,jjj,kk,kkk,ll,lll,f[iii][jjj][kkk][lll],dist[iii][jjj][kkk][lll],totalf);

    }// end if good temperature or doing degentable
    else{

      // just use nearest neighbor if no valid inversion
      if(whichtable==SIMPLEZOOMTABLE){
	if(degentable==0) tfptr[iii][jjj][kkk][lll] = eossimplezoomtable[whichfun][0][llo+lll][kko+kkk][jjo+jjj][iio+iii];
	else  tfptr[iii][jjj][kkk][lll] = eosdegensimplezoomtable[whichdegenfun][0][llo+lll][kko+kkk][0][iio+iii];
      }
      else if(whichtable==SIMPLETABLE){
	if(degentable==0) tfptr[iii][jjj][kkk][lll] = eossimpletable[whichfun][0][llo+lll][kko+kkk][jjo+jjj][iio+iii];
	else  tfptr[iii][jjj][kkk][lll] = eosdegensimpletable[whichdegenfun][0][llo+lll][kko+kkk][0][iio+iii];
      }
      else if(whichtable==FULLTABLE){
	if(degentable==0) tfptr[iii][jjj][kkk][lll] = eostable[whichfun][0][llo+lll][kko+kkk][jjo+jjj][iio+iii];
	else tfptr[iii][jjj][kkk][lll] = eosdegentable[whichdegenfun][0][llo+lll][kko+kkk][0][iio+iii];
      }

      if(debugfail>=2){ // DEBUG GODMARK: was turned on when debugging EOS
	// DEBUG
	dualfprintf(fail_file,"Out of Bounds based upon temperature: icurr=%d :: whichfun=%d ii=%d iii=%d jj=%d jjj=%d kk=%d kkk=%d ll=%d lll=%d temp=%21.15g invalidtempcode=%21.15g\n",icurr,whichfun,ii,iii,jj,jjj,kk,kkk,ll,lll,tempcheck,invalidtempcode);
      }

    }



#if(DOLOGINTERP)
      if(loginterp){
	//	if(tfptr[iii][jjj][kkk][lll]<=0.0){
	//	  dualfprintf(fail_file,"Negative of log10 icurr=%d\n",icurr);
	//	}
	tfptr[iii][jjj][kkk][lll] = log10(tfptr[iii][jjj][kkk][lll]);
      }
#endif


  }// end loop over dimensions





  // perform interpolation over values
  for(jjj=startjjj;jjj<=endjjj;jjj++)for(kkk=startkkk;kkk<=endkkk;kkk++)for(lll=startlll;lll<=endlll;lll++){


    //    totalfptr=(FTYPE (*)[3][3][3]) (&(totalf[1][1][1][1])); // so tfptr[-1,0,1]

    // now use 3 data points to get density-parabolic distribution and value at ieos
    // have tfptr[iii][jjj][kkk][lll] @ iii=-1,0,1 with i=0 meaning ii and offset being ROUND2INT(ieos)
    // Form parabolic answer
    xmx0 = (ieos-(FTYPE)ii);
    AA = 0.5*(tfptr[1][jjj][kkk][lll]-tfptr[-1][jjj][kkk][lll]);
    BB = 0.5*(tfptr[1][jjj][kkk][lll]+tfptr[-1][jjj][kkk][lll]-2.0*tfptr[0][jjj][kkk][lll]);
    tfptr[0][jjj][kkk][lll] = tfptr[0][jjj][kkk][lll] + AA*xmx0 + BB*xmx0*xmx0;

    if(!isfinite(tfptr[0][jjj][kkk][lll])){
      dualfprintf(fail_file,"1notfinite, icurr=%d :: %d %d %d\n",icurr,jjj,kkk,lll);
    }
    
  } // end over iii,jjj,kkk,lll

  // perform interpolation over values
  for(kkk=startkkk;kkk<=endkkk;kkk++)for(lll=startlll;lll<=endlll;lll++){

    xmx0 = (jeos-(FTYPE)jj);
    AA = 0.5*(tfptr[0][1][kkk][lll]-tfptr[0][-1][kkk][lll]);
    BB = 0.5*(tfptr[0][1][kkk][lll]+tfptr[0][-1][kkk][lll]-2.0*tfptr[0][0][kkk][lll]);
    tfptr[0][0][kkk][lll] = tfptr[0][0][kkk][lll] + AA*xmx0 + BB*xmx0*xmx0;

    if(!isfinite(tfptr[0][0][kkk][lll])){
      dualfprintf(fail_file,"2notfinite, icurr=%d :: %d %d %d\n",icurr,0,kkk,lll);
    }


  } // end over kkk,lll

  // perform interpolation over values
  for(lll=startlll;lll<=endlll;lll++){

    xmx0 = (keos-(FTYPE)kk);
    AA = 0.5*(tfptr[0][0][1][lll]-tfptr[0][0][-1][lll]);
    BB = 0.5*(tfptr[0][0][1][lll]+tfptr[0][0][-1][lll]-2.0*tfptr[0][0][0][lll]);
    tfptr[0][0][0][lll] = tfptr[0][0][0][lll] + AA*xmx0 + BB*xmx0*xmx0;

    if(!isfinite(tfptr[0][0][0][lll])){
      dualfprintf(fail_file,"3notfinite, icurr=%d :: %d %d %d\n",icurr,0,0,lll);
    }

    
  } // end over lll

  // final interpolation
  xmx0 = (leos-(FTYPE)ll);
  AA = 0.5*(tfptr[0][0][0][1]-tfptr[0][0][0][-1]);
  BB = 0.5*(tfptr[0][0][0][1]+tfptr[0][0][0][-1]-2.0*tfptr[0][0][0][0]);
  tfptr[0][0][0][0] = tfptr[0][0][0][0] + AA*xmx0 + BB*xmx0*xmx0;
  
  if(!isfinite(tfptr[0][0][0][0])){
    dualfprintf(fail_file,"4notfinite, icurr=%d :: %d %d %d\n",icurr,0,0,0);
  }
  dualfprintf(fail_file,"4finite, icurr=%d :: %d %d %d : %21.15g\n",icurr,0,0,0,tfptr[0][0][0][0]);


  totalffinal=tfptr[0][0][0][0];



#if(DOLOGINTERP)
  if(loginterp){
    totalffinal=pow(10.0,totalffinal);
  }
#endif



  return(totalffinal);

}











// tri-linear + parabolic (density) interpolation
FTYPE get_eos_fromlookup_parabolic(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichfun, int whichindep, int *vartypearray, FTYPE *indexarray)
{
  static int ii,jj,kk,ll;
  static FTYPE dj[2],dk[2],dl[2];
  static int startiii,startjjj,startkkk,startlll;
  static int endiii,endjjj,endkkk,endlll;
  FTYPE totaldist[3];
  FTYPE *tdist;
  FTYPE dist[2][2][2],f[2][2][2];
  FTYPE tempcheck;
  FTYPE totalf[3]; // 3 values for parabolic interpolation
  FTYPE totalffinal;
  FTYPE *tfptr;
  int iii,jjj,kkk,lll;
  int whichtemp,whichdegenfun;
  int whichinterp1,whichinterp2,loginterp;
  FTYPE xmx0,AA,BB;
  FTYPE ieos,jeos,keos,leos,meos;
  int EXTRASTART,EXTRAFINISH;


  ieos=indexarray[1];
  jeos=indexarray[2];
  keos=indexarray[3];
  leos=indexarray[4];




  tfptr=&totalf[1]; // so tfptr[-1,0,1]
  tdist=&totaldist[1];

  // definition consistent with numerical assignments of indecies of arrays
  whichdegenfun = whichindep-1;


  // determine which temperature to use to check inversion
  // GODMARK: this could be stored in array where index is whichindep and output is whichtemp
  if(whichindep==UEOS) whichtemp=TEMPU;
  else if(whichindep==PEOS) whichtemp=TEMPP;
  else if(whichindep==CHIEOS) whichtemp=TEMPCHI;





#if(DOLOGINTERP)

  EXTRASTART=extralimits[whichdatatype[whichtable]-1][0];
  EXTRAFINISH=extralimits[whichdatatype[whichtable]-1][1];

  // GODMARK: Can make array that stores this info, looked up by whichfun as index
  // functions (F) F(rho0,u)
  whichinterp1=(whichfun==PofRHOCHI||whichfun==UofRHOP||whichfun==TEMPP||whichfun==PofRHOU||whichfun==CS2ofRHOU||whichfun==SofRHOU||(whichfun>=EXTRASTART && whichfun<=EXTRAFINISH)||whichfun==TEMPU||whichfun==TEMPCHI);
  // functions (F) F(rho0,p)
  whichinterp2=(whichfun==DPDRHOofRHOU||whichfun==DPDUofRHOU||whichfun==DSDRHOofRHOU||whichfun==DSDUofRHOU||whichfun==IDRHO0DP||whichfun==IDCHIDP);

  if(whichinterp1||degentable==1) loginterp=1;
  else if(whichinterp2) loginterp=0;
  else{
    dualfprintf(fail_file,"Undefined whichfun=%d in get_eos_fromlookup_linear()\n",whichfun);
    myexit(62662);
  }
#endif



  if(repeatedeos==0){
    // assume jeos, keos, leos set correctly for given degentable and tabledimen (i.e. are 0 or integers in those cases so d's computed correctly)
    ii=ROUND2INT(ieos); //ii=(int)ieos; // round instead when doing parabolic interpolation with 3 points
    // limit to within table
    if(ii<1) ii=1; // 0 is minimum
    if(ii>tablesize[whichtable][RHOEOS]-2) ii=tablesize[whichtable][RHOEOS]-2; // N-1 is maximum

    // normal behavior for other independent variables
    jj=(int)jeos;
    kk=(int)keos;
    ll=(int)leos;
	
    dj[1]=jeos-(FTYPE)jj;
    dj[0]=1.0-dj[1];
    dk[1]=keos-(FTYPE)kk;
    dk[0]=1.0-dk[1];
    dl[1]=leos-(FTYPE)ll;
    dl[0]=1.0-dl[1];


    // set range of loops for different table types
    // tabledimen overrules table type (i.e. take section out of fuller table -- assumed to be 0 index)
    // GODMARK: these things couuld be stored as functions of whichtable/degentable/tabledimen
    if(tablesize[whichtable][vartypearray[1]]!=1) endiii=1;
    else endiii=0;

    if(degentable==0 || tablesize[whichtable][vartypearray[2]]!=1) endjjj=1;
    else endjjj=0;

    if(tablesize[whichtable][vartypearray[3]]!=1) endkkk=1;
    else endkkk=0;
    
    if(tablesize[whichtable][vartypearray[4]]!=1) endlll=1;
    else endlll=0;

  }


    
  // Loop over nearby table values and determine bi-linearly interpolated value
  // 4-D means 2^4=16 positions
  // 2-D means 2^2=4 positions
  // get 3 values as function of density
  for(iii=startiii;iii<=endiii;iii++){
    tfptr[iii]=0.0;
    tdist[iii]=0.0;
  }

  for(iii=startiii;iii<=endiii;iii++) for(jjj=0;jjj<=endjjj;jjj++)for(kkk=0;kkk<=endkkk;kkk++)for(lll=0;lll<=endlll;lll++){


#if(CHECKIFVALIDEOSDATA)
    if(degentable==0){
      // don't use values of table that have no inversion to temperature
      if(whichtable==SIMPLEZOOMTABLE) tempcheck=eossimplezoomtable[whichtemp][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
      else if(whichtable==SIMPLETABLE) tempcheck=eossimpletable[whichtemp][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
      else if(whichtable==FULLTABLE) tempcheck=eostable[whichtemp][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
    }
    //    dualfprintf(fail_file,"temp=%21.15g\n",tempcheck);
    if(degentable==1 || tempcheck>invalidtempcode) // Avoid invalid inversions if T>Tbad, but only deal with temperature if degentable==0
#else
    if(1)
#endif
    {
      tdist[iii] += dist[jjj][kkk][lll] = dj[jjj]*dk[kkk]*dl[lll];
      if(whichtable==SIMPLEZOOMTABLE){
	if(degentable==0) f[jjj][kkk][lll] = eossimplezoomtable[whichfun][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
	else f[jjj][kkk][lll] = eosdegensimplezoomtable[whichdegenfun][0][ll+lll][kk+kkk][0][ii+iii];
      }
      else if(whichtable==SIMPLETABLE){
	if(degentable==0) f[jjj][kkk][lll] = eossimpletable[whichfun][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
	else f[jjj][kkk][lll] = eosdegensimpletable[whichdegenfun][0][ll+lll][kk+kkk][0][ii+iii];
      }
      else if(whichtable==FULLTABLE){
	if(degentable==0) f[jjj][kkk][lll] = eostable[whichfun][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
	else f[jjj][kkk][lll] = eosdegentable[whichdegenfun][0][ll+lll][kk+kkk][0][ii+iii];

      }

#if(DOLOGINTERP)
      if(loginterp){
	f[jjj][kkk][lll] = log10(f[jjj][kkk][lll]);
      }
#endif

      //dualfprintf(fail_file,"f[%d][%d][%d][%d]=%21.15g\n",iii,jjj,kkk,lll,f[jjj][kkk][lll]);
      tfptr[iii] +=f[jjj][kkk][lll]*dist[jjj][kkk][lll];

      // DEBUG
      //      dualfprintf(fail_file,"tabledimen=%d degentable=%d whichtable=%d whichfun=%d whichdegenfun=%d ii=%d iii=%d jj=%d jjj=%d kk=%d kkk=%d ll=%d lll=%d :: f=%21.15g dist=%21.15g totalf=%21.15g\n",tabledimen, degentable, whichtable,whichfun,whichdegenfun,ii,iii,jj,jjj,kk,kkk,ll,lll,f[iii][jjj][kkk][lll],dist[iii][jjj][kkk][lll],totalf);

    }// end if good temperature or doing degentable
    else{

      if(debugfail>=2){ // DEBUG GODMARK: was turned on when debugging EOS
	// DEBUG
	dualfprintf(fail_file,"Out of Bounds based upon temperature: icurr=%d :: whichfun=%d ii=%d iii=%d jj=%d jjj=%d kk=%d kkk=%d ll=%d lll=%d temp=%21.15g invalidtempcode=%21.15g\n",icurr,whichfun,ii,iii,jj,jjj,kk,kkk,ll,lll,tempcheck,invalidtempcode);
      }


    }
  }// end loop over dimensions


  for(iii=startiii;iii<=endiii;iii++) {
    ////////////////////////////
    //
    // finally normalize
    if(tdist[iii]==0.0){
      // Good to know if not even close
      // GODMARK: Should probably store mapping to closest EOS data in table and use that mapping in this case (like extrapolation, but since no rho0,u pair, assume not using valid "u" and so should really change "u" if we want to be consistent.

      // if only choosing 1 value, then choose rounded version
      ii=ROUND2INT(ieos);
      jj=ROUND2INT(jeos);
      kk=ROUND2INT(keos);
      ll=ROUND2INT(leos);


      // just use nearest neighbor if no valid inversion
      if(whichtable==SIMPLEZOOMTABLE){
	if(degentable==0) tfptr[iii] = eossimplezoomtable[whichfun][0][ll][kk][jj][ii+iii];
	else  tfptr[iii] = eosdegensimplezoomtable[whichdegenfun][0][ll][kk][0][ii+iii];
      }
      else if(whichtable==SIMPLETABLE){
	if(degentable==0) tfptr[iii] = eossimpletable[whichfun][0][ll][kk][jj][ii+iii];
	else  tfptr[iii] = eosdegensimpletable[whichdegenfun][0][ll][kk][0][ii+iii];
      }
      else if(whichtable==FULLTABLE){
	if(degentable==0) tfptr[iii] = eostable[whichfun][0][ll][kk][jj][ii+iii];
	else tfptr[iii] = eosdegentable[whichdegenfun][0][ll][kk][0][ii+iii];
      }

#if(DOLOGINTERP)
      if(loginterp){
	tfptr[iii] = log10(tfptr[iii]);
      }
#endif


      dualfprintf(fail_file,"Never found temperature: degentable=%d whichtable=%d whichfun=%d whichindep=%d whichdegenfun=%d ieos=%21.15g jeos=%21.15g keos=%21.15g leos=%21.15g tfptr[iii]=%21.15g\n",degentable,whichtable,whichfun,whichindep,whichdegenfun,ieos,jeos,keos,leos,tfptr[iii]);

    }
    else{
      tfptr[iii] /=tdist[iii];

    }

    // GODMARK DEBUG:
    if(whichfun==TEMPU){
      dualfprintf(fail_file,"ii=%d jj=%d kk=%d ll=%d :: degen=%d :: %g\n",ii,jj,kk,ll,degentable,tfptr[iii]);
    }

  } // end over iii


  // now use 3 data points to get density-parabolic distribution and value at ieos
  // have tfptr[iii] @ iii=-1,0,1 with i=0 meaning ii and offset being ROUND2INT(ieos)
  // Form parabolic answer
  xmx0 = (ieos-(FTYPE)ii);
  AA = 0.5*(tfptr[1]-tfptr[-1]);
  BB = 0.5*(tfptr[1]+tfptr[-1]-2.0*tfptr[0]);
  totalffinal = tfptr[0] + AA*xmx0 + BB*xmx0*xmx0;

  // DEBUG:
  //  if(whichfun==TEMPU){
  //    dualfprintf(fail_file,"ieos=%g ii=%d :: tfptr: %21.15g %21.15g %21.15g\n",ieos,ii,tfptr[-1],tfptr[0],tfptr[1]);
  //  }



#if(DOLOGINTERP)
  if(loginterp){
    totalffinal=pow(10.0,totalffinal);
  }
#endif



  return(totalffinal);

}












// full 5D linear interpolation
FTYPE get_eos_fromlookup_linear(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichfun, int whichindep, int *vartypearray, FTYPE *indexarray)
{
  static int ii,jj,kk,ll,mm;
  static FTYPE di[2],dj[2],dk[2],dl[2],dm[2];
  static int endiii,endjjj,endkkk,endlll,endmmm;
  FTYPE totaldist;
  FTYPE dist[2][2][2][2][2],f[2][2][2][2][2];
  FTYPE tempcheck;
  FTYPE totalf;
  int iii,jjj,kkk,lll,mmm;
  int whichtemp,whichdegenfun;
  int whichinterp1,whichinterp2,loginterp;
  FTYPE ieos,jeos,keos,leos,meos;
  int EXTRASTART,EXTRAFINISH;
  int qi;




  ieos=indexarray[1];
  jeos=indexarray[2];
  keos=indexarray[3];
  leos=indexarray[4];
  meos=indexarray[5];
  



  // definition consistent with numerical assignments of indecies of arrays
  whichdegenfun = whichindep-1;


  // determine which temperature to use to check inversion
  // GODMARK: this could be stored in array where index is whichindep and output is whichtemp
  if(whichindep==UEOS) whichtemp=TEMPU;
  else if(whichindep==PEOS) whichtemp=TEMPP;
  else if(whichindep==CHIEOS) whichtemp=TEMPCHI;





#if(DOLOGINTERP)

  EXTRASTART=extralimits[whichdatatype[whichtable]-1][0];
  EXTRAFINISH=extralimits[whichdatatype[whichtable]-1][1];

  // GODMARK: Can make array that stores this info, looked up by whichfun as index
  // functions (F) F(rho0,u)
  whichinterp1=(whichfun==PofRHOCHI||whichfun==UofRHOP||whichfun==TEMPP||whichfun==PofRHOU||whichfun==CS2ofRHOU||whichfun==SofRHOU||(whichfun>=EXTRASTART && whichfun<=EXTRAFINISH)||whichfun==TEMPU||whichfun==TEMPCHI);
  // functions (F) F(rho0,p)
  whichinterp2=(whichfun==DPDRHOofRHOU||whichfun==DPDUofRHOU||whichfun==DSDRHOofRHOU||whichfun==DSDUofRHOU||whichfun==IDRHO0DP||whichfun==IDCHIDP);

  if(whichinterp1||degentable==1) loginterp=1;
  else if(whichinterp2) loginterp=0;
  else{
    dualfprintf(fail_file,"Undefined whichfun=%d in get_eos_fromlookup_linear(): %d %d %d %d %d %d\n",whichfun, repeatedeos, tabledimen, degentable, whichtable, whichfun, whichindep);
    for(qi=1;qi<=NUMINDEPDIMENS+1;qi++) dualfprintf(fail_file,"%d : vartypearray=%d indexarray=%d\n",qi,vartypearray[qi],indexarray[qi]);
    myexit(62662);
  }
#endif



  if(repeatedeos==0){
    // assume jeos, keos, leos set correctly for given degentable and tabledimen (i.e. are 0 or integers in those cases so d's computed correctly)
    ii=(int)ieos;
    //if(degentable) jj=0; else jj=(int)jeos; // force to be 0
    jj=(int)jeos; // (int)(0.0) = 0, so ok as is
    kk=(int)keos;
    ll=(int)leos;
    mm=(int)meos;
	
    di[1]=ieos-(FTYPE)ii;
    di[0]=1.0-di[1];
    dj[1]=jeos-(FTYPE)jj;
    dj[0]=1.0-dj[1];
    dk[1]=keos-(FTYPE)kk;
    dk[0]=1.0-dk[1];
    dl[1]=leos-(FTYPE)ll;
    dl[0]=1.0-dl[1];
    dm[1]=meos-(FTYPE)mm;
    dm[0]=1.0-dm[1];


    // set range of loops for different table types
    // tabledimen overrules table type (i.e. take section out of fuller table -- assumed to be 0 index)
    // GODMARK: these things couuld be stored as functions of whichtable/degentable/tabledimen
    if(tablesize[whichtable][vartypearray[1]]!=1) endiii=1;
    else endiii=0;

    if(degentable==0 || tablesize[whichtable][vartypearray[2]]!=1) endjjj=1;
    else endjjj=0;

    if(tablesize[whichtable][vartypearray[3]]!=1) endkkk=1;
    else endkkk=0;
    
    if(tablesize[whichtable][vartypearray[4]]!=1) endlll=1;
    else endlll=0;

    if(tablesize[whichtable][vartypearray[5]]!=1) endmmm=1;
    else endmmm=0;


  }


    
  // Loop over nearby table values and determine bi-linearly interpolated value
  // 4-D means 2^4=16 positions
  // 2-D means 2^2=4 positions
  totaldist=0.0;
  totalf=0.0;
  for(mmm=0;mmm<=endmmm;mmm++)for(lll=0;lll<=endlll;lll++)for(kkk=0;kkk<=endkkk;kkk++)for(jjj=0;jjj<=endjjj;jjj++)for(iii=0;iii<=endiii;iii++){


#if(CHECKIFVALIDEOSDATA)
    if(degentable==0){
      // don't use values of table that have no inversion to temperature
      if(whichtable==SIMPLEZOOMTABLE) tempcheck=eossimplezoomtable[whichtemp][mm+mmm][ll+lll][kk+kkk][jj+jjj][ii+iii];
      else if(whichtable==SIMPLETABLE) tempcheck=eossimpletable[whichtemp][mm+mmm][ll+lll][kk+kkk][jj+jjj][ii+iii];
      else if(whichtable==FULLTABLE) tempcheck=eostable[whichtemp][mm+mmm][ll+lll][kk+kkk][jj+jjj][ii+iii];
    }
    //    dualfprintf(fail_file,"temp=%21.15g\n",tempcheck);
    if(degentable==1 || tempcheck>invalidtempcode) // Avoid invalid inversions if T>Tbad, but only deal with temperature if degentable==0
#else
    if(1)
#endif
    {
      totaldist += dist[iii][jjj][kkk][lll][mmm] = di[iii]*dj[jjj]*dk[kkk]*dl[lll]*dl[mmm];
      if(whichtable==SIMPLEZOOMTABLE){
	if(degentable==0) f[iii][jjj][kkk][lll][mmm] = eossimplezoomtable[whichfun][mm+mmm][ll+lll][kk+kkk][jj+jjj][ii+iii];
	else f[iii][jjj][kkk][lll][mmm] = eosdegensimplezoomtable[whichdegenfun][mm+mmm][ll+lll][kk+kkk][0][ii+iii];
      }
      else if(whichtable==SIMPLETABLE){
	if(degentable==0) f[iii][jjj][kkk][lll][mmm] = eossimpletable[whichfun][mm+mmm][ll+lll][kk+kkk][jj+jjj][ii+iii];
	else f[iii][jjj][kkk][lll][mmm] = eosdegensimpletable[whichdegenfun][mm+mmm][ll+lll][kk+kkk][0][ii+iii];
      }
      else if(whichtable==FULLTABLE){
	if(degentable==0) f[iii][jjj][kkk][lll][mmm] = eostable[whichfun][mm+mmm][ll+lll][kk+kkk][jj+jjj][ii+iii];
	else f[iii][jjj][kkk][lll][mmm] = eosdegentable[whichdegenfun][mm+mmm][ll+lll][kk+kkk][0][ii+iii];

      }

#if(DOLOGINTERP)
      if(loginterp){
	f[iii][jjj][kkk][lll][mmm] = log10(f[iii][jjj][kkk][lll][mmm]);
      }
#endif

      //dualfprintf(fail_file,"f[%d][%d][%d][%d]=%21.15g\n",iii,jjj,kkk,lll,f[iii][jjj][kkk][lll]);
      totalf +=f[iii][jjj][kkk][lll][mmm]*dist[iii][jjj][kkk][lll][mmm];

      // DEBUG
      //      dualfprintf(fail_file,"tabledimen=%d degentable=%d whichtable=%d whichfun=%d whichdegenfun=%d ii=%d iii=%d jj=%d jjj=%d kk=%d kkk=%d ll=%d lll=%d :: f=%21.15g dist=%21.15g totalf=%21.15g\n",tabledimen, degentable, whichtable,whichfun,whichdegenfun,ii,iii,jj,jjj,kk,kkk,ll,lll,f[iii][jjj][kkk][lll],dist[iii][jjj][kkk][lll],totalf);

    }// end if good temperature or doing degentable
    else{


      if(1|| debugfail>=2){ // DEBUG GODMARK: was turned on when debugging EOS
	// DEBUG
	dualfprintf(fail_file,"Out of Bounds based upon temperature: icurr=%d :: whichfun=%d ii=%d iii=%d jj=%d jjj=%d kk=%d kkk=%d ll=%d lll=%d mm=%d mmm=%d :: temp=%21.15g invalidtempcode=%21.15g\n",icurr,whichfun,ii,iii,jj,jjj,kk,kkk,ll,lll,mm,mmm,tempcheck,invalidtempcode);
      }

      
    }
  }// end loop over dimensions


  ////////////////////////////
  //
  // finally normalize
  if(totaldist==0.0){
    // Good to know if not even close
    // GODMARK: Should probably store mapping to closest EOS data in table and use that mapping in this case (like extrapolation, but since no rho0,u pair, assume not using valid "u" and so should really change "u" if we want to be consistent.

    // if only choosing 1 value, then choose rounded version
    ii=ROUND2INT(ieos);
    jj=ROUND2INT(jeos);
    kk=ROUND2INT(keos);
    ll=ROUND2INT(leos);
    mm=ROUND2INT(meos);


    // just use nearest neighbor if no valid inversion
    if(whichtable==SIMPLEZOOMTABLE){
      if(degentable==0) totalf = eossimplezoomtable[whichfun][mm][ll][kk][jj][ii];
      else  totalf = eosdegensimplezoomtable[whichdegenfun][mm][ll][kk][0][ii];
    }
    else if(whichtable==SIMPLETABLE){
      if(degentable==0) totalf = eossimpletable[whichfun][mm][ll][kk][jj][ii];
      else  totalf = eosdegensimpletable[whichdegenfun][mm][ll][kk][0][ii];
    }
    else if(whichtable==FULLTABLE){
      if(degentable==0) totalf = eostable[whichfun][mm][ll][kk][jj][ii];
      else totalf = eosdegentable[whichdegenfun][mm][ll][kk][0][ii];
    }

    dualfprintf(fail_file,"Never found temperature: degentable=%d whichtable=%d whichfun=%d whichindep=%d whichdegenfun=%d ieos=%21.15g jeos=%21.15g keos=%21.15g leos=%21.15g totalf=%21.15g\n",degentable,whichtable,whichfun,whichindep,whichdegenfun,ieos,jeos,keos,leos,totalf);

  }
  else{
    totalf /=totaldist;

#if(DOLOGINTERP)
    if(loginterp){
      totalf=pow(10.0,totalf);
    }
#endif

  }

  //  dualfprintf(fail_file,"ii=%d jj=%d kk=%d ll=%d\n",ii,jj,kk,ll);
  //  dualfprintf(fail_file,"totaldist=%21.15g totalf=%21.15g\n",totaldist,totalf);

  // GODMARK DEBUG:
  if(whichfun==TEMPU){
    dualfprintf(fail_file,"ii=%d jj=%d kk=%d ll=%d mm=%d :: degen=%d :: %g\n",ii,jj,kk,ll,mm,degentable,totalf);
  }


  return(totalf);

}










// nearest neighbor interpolation but caution to temperature defined or not
FTYPE get_eos_fromlookup_nearest(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichfun, int whichindep, int *vartypearray, FTYPE *indexarray)
{
  static int ii,jj,kk,ll;
  static int endiii,endjjj,endkkk,endlll;
  FTYPE totaldist;
  FTYPE tempcheck;
  FTYPE totalf;
  int iii,jjj,kkk,lll;
  int whichtemp,whichdegenfun;
  FTYPE ieos,jeos,keos,leos,meos;


  ieos=indexarray[1];
  jeos=indexarray[2];
  keos=indexarray[3];
  leos=indexarray[4];




  // definition consistent with numerical assignments of indecies of arrays
  whichdegenfun = whichindep-1;


  // determine which temperature to use to check inversion
  // GODMARK: this could be stored in array where index is whichindep and output is whichtemp
  if(whichindep==UEOS) whichtemp=TEMPU;
  else if(whichindep==PEOS) whichtemp=TEMPP;
  else if(whichindep==CHIEOS) whichtemp=TEMPCHI;




  if(repeatedeos==0){
    // assume jeos, keos, leos set correctly for given degentable and tabledimen (i.e. are 0 or integers in those cases so d's computed correctly)
    ii=(int)ieos;
    jj=(int)jeos;
    kk=(int)keos;
    ll=(int)leos;
	

    // set range of loops for different table types
    // tabledimen overrules table type (i.e. take section out of fuller table -- assumed to be 0 index)
    // GODMARK: these things couuld be stored as functions of whichtable/degentable/tabledimen
    if(tablesize[whichtable][vartypearray[1]]!=1) endiii=1;
    else endiii=0;

    if(degentable==0 || tablesize[whichtable][vartypearray[2]]!=1) endjjj=1;
    else endjjj=0;

    if(tablesize[whichtable][vartypearray[3]]!=1) endkkk=1;
    else endkkk=0;
    
    if(tablesize[whichtable][vartypearray[4]]!=1) endlll=1;
    else endlll=0;
 
  }



  totaldist=0.0;
  // loop over to ensure obtained best nearest neighbor    
  for(iii=0;iii<=endiii;iii++)for(jjj=0;jjj<=endjjj;jjj++)for(kkk=0;kkk<=endkkk;kkk++)for(lll=0;lll<=endlll;lll++){
    
    
#if(CHECKIFVALIDEOSDATA)
    if(degentable==0){
      // don't use values of table that have no inversion to temperature
      if(whichtable==SIMPLEZOOMTABLE) tempcheck=eossimplezoomtable[whichtemp][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
      else if(whichtable==SIMPLETABLE) tempcheck=eossimpletable[whichtemp][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
      else if(whichtable==FULLTABLE) tempcheck=eostable[whichtemp][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
    }
    //    dualfprintf(fail_file,"temp=%21.15g\n",tempcheck);
    if(degentable==1 || tempcheck>invalidtempcode) // Avoid invalid inversions if T>Tbad, but only deal with temperature if degentable==0
#else
    if(1)
#endif
    {
      if(whichtable==SIMPLEZOOMTABLE){
	if(degentable==0) totalf = eossimplezoomtable[whichfun][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
	else totalf = eosdegensimplezoomtable[whichdegenfun][0][ll+lll][kk+kkk][0][ii+iii];
      }
      else if(whichtable==SIMPLETABLE){
	if(degentable==0) totalf = eossimpletable[whichfun][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
	else totalf = eosdegensimpletable[whichdegenfun][0][ll+lll][kk+kkk][0][ii+iii];
      }
      else if(whichtable==FULLTABLE){
	if(degentable==0) totalf = eostable[whichfun][0][ll+lll][kk+kkk][jj+jjj][ii+iii];
	else totalf = eosdegentable[whichdegenfun][0][ll+lll][kk+kkk][0][ii+iii];
      }
      // terminate all 4 loops once have a single value
      iii=endiii;
      jjj=endjjj;
      kkk=endkkk;
      lll=endlll;
      totaldist=1.0;

    }// end if good temperature or doing degentable
    else{


      if(debugfail>=2){ // DEBUG GODMARK: was turned on when debugging EOS
	// DEBUG
	dualfprintf(fail_file,"Out of Bounds based upon temperature: icurr=%d :: whichfun=%d ii=%d iii=%d jj=%d jjj=%d kk=%d kkk=%d ll=%d lll=%d temp=%21.15g invalidtempcode=%21.15g\n",icurr,whichfun,ii,iii,jj,jjj,kk,kkk,ll,lll,tempcheck,invalidtempcode);
      }


    }
  }// end loop over dimensions


  ////////////////////////////
  //
  // finally normalize
  if(totaldist==0.0){
    dualfprintf(fail_file,"Never found temperature: degentable=%d whichtable=%d whichfun=%d whichindep=%d whichdegenfun=%d ieos=%21.15g jeos=%21.15g keos=%21.15g leos=%21.15g totalf=%21.15g\n",degentable,whichtable,whichfun,whichindep,whichdegenfun,ieos,jeos,keos,leos,totalf);

  }

  return(totalf);

}


//#define thiseostable(name,fun,index5,index4,index3,index2,index1) name[fun][index5][index4][index3][index2][index1]

//#define index2array(name,fun) name[fun][indexarray[5]][indexarray[4]][indexarray[3]][indexarray[2]][indexarray[1]]

// nearest neighbor interpolation with no temperature check
FTYPE get_eos_fromlookup_nearest_dumb(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichfun, int whichindep, int *vartypearray, FTYPE *indexarray)
{
  static int ii,jj,kk,ll;
  FTYPE totalf;
  int whichdegenfun;
  FTYPE ieos,jeos,keos,leos,meos;


  ieos=indexarray[1];
  jeos=indexarray[2];
  keos=indexarray[3];
  leos=indexarray[4];



  // definition consistent with numerical assignments of indecies of arrays
  whichdegenfun = whichindep-1;


#if(0)
  ii=(int)ieos;
  jj=(int)jeos;
  kk=(int)keos;
  ll=(int)leos;
#else
  // if only choosing 1 value, then choose rounded version
  ii=ROUND2INT(ieos);
  jj=ROUND2INT(jeos);
  kk=ROUND2INT(keos);
  ll=ROUND2INT(leos);
#endif

  
  if(whichtable==SIMPLEZOOMTABLE){
    if(degentable==0) totalf = eossimplezoomtable[whichfun][0][ll][kk][jj][ii];
    else totalf = eosdegensimplezoomtable[whichdegenfun][0][ll][kk][0][ii];
  }
  else if(whichtable==SIMPLETABLE){
    if(degentable==0) totalf = eossimpletable[whichfun][0][ll][kk][jj][ii];
    else totalf = eosdegensimpletable[whichdegenfun][0][ll][kk][0][ii];
  }
  else if(whichtable==FULLTABLE){
    if(degentable==0) totalf = eostable[whichfun][0][ll][kk][jj][ii];
    else totalf = eosdegentable[whichdegenfun][0][ll][kk][0][ii];
  }

  if(whichfun==TEMPU){
    dualfprintf(fail_file,"TEMPUCHECK1 :: icurr=%d :: ieos=%g jeos=%g keos=%g leos=%g :: ii=%d jj=%d kk=%d ll=%d :: degen=%d :: %g\n",icurr,ieos,jeos,keos,leos,ii,jj,kk,ll,degentable,totalf);
  }


  return(totalf);

}










////////////////////////////
//
// Choose behavior of loop depending upon whether including degen offset
//
////////////////////////////
// GODMARK: Could set this 
#if(ALLOWDEGENOFFSET)
// then need to lookup degen part first
#define whichdegenstart 1
#define whichdegenend 0
#else
#define whichdegenstart 0
#define whichdegenend 0
#endif

// for now assume TDYNORYE is locked to timestep that is order 1E-6 seconds, so can use lowest leos
// which = which function looking up
// q1 = rhob
// q2 = u or P or \chi
// notice that q1-q2 are pass by value, so internally changed but externally remains unchanged
static int get_eos_fromtable(int whichfun, FTYPE quant1, FTYPE quant2, FTYPE *answer)
{
  FTYPE get_eos_fromlookup(int repeatedeos, int tabledimen, int degentable, int whichtable, int whichfun, int whichindep, int *vartypearray, FTYPE *indexarray);
  int iswithin_eostable(int ifdegencheck, int whichindep, int *vartypearray, FTYPE *qarray, int *whichtable);

  void eos_lookup(int whichtable, int whichindep, int *vartypearray, FTYPE *qarray, FTYPE *indexarray);
  void eos_lookup_prepost_degen(int whichdegen, int whichtable, int whichindep, int *vartypearray, FTYPE *qarray, FTYPE *indexarray);

  int whichcheck1,whichcheck2,whichcheck3;
  int whichindep;
  FTYPE myanswer[2]; // 0 = non-degen answer  1=degen answer
  int whichdegen;
  FTYPE qarray[NUMINDEPDIMENS+1];
  FTYPE qfloor[NUMINDEPDIMENS+1];
  int i;
  int qi;
  int getwhichtable;
  int repeatedeos;
  // variables need to keep so can quickly get EOS values if repeated input q1-q5
  static int gottable[2]={0,0}; // for degen and non-degen case : should end up same!
  static int whichtable[2];
  static FTYPE indexarray[NUMINDEPDIMENS+1];
  static FTYPE qoldarray[NUMINDEPDIMENS+1];
  static FTYPE resultold[NUMEOSQUANTITIESMEM];
  static int repeatedfun[NUMEOSQUANTITIESMEM];
  int get_whichindep(int whichfun);




  if(didsetupkazeos==0){
    dualfprintf(fail_file,"Requested WHICHEOS=KAZEOS but didn't initialize kaz table by calling read_setup_eostable()\n");
    myexit(151068); 
  }

  //////////////////////
  //
  // Determine which independent variable to use for "temperature": utot, ptot, or chi (or utotdiff, ptotdiff, chidiff for degen case)
  //
  //////////////////////

  whichindep = get_whichindep(whichfun); // gets UEOS, PEOS, or CHIEOS
  vartypearray[2]=whichindep; // used to translate q1-q5 to INDEP0-INDEP6, and used with limits and sizes of tables




  ////////////////////////////
  //
  // Set independent EOS quantities
  //
  ////////////////////////////


  // setup array
  qarray[1]=quant1;
  qarray[2]=quant2;
  // GODMARK: assumes anytime EOS is accessed that icurr/jcurr/kcurr set!
  // qarray[FIRSTEOSGLOBAL+] are always stored in EOSglobal
  for(qi=FIRSTEOSGLOBAL;qi<=NUMINDEPDIMENS;qi++){
    qarray[qi] = EOSglobal[qi][icurr][jcurr][kcurr];
  }


  //////////////////////////////////////////////////////
  //
  // check if already got here with these q values (done after assignment of q3-q5, so distinct -- done before floors on q since they make less distinct anyways)
  //
  //////////////////////////////////////////////////////

  repeatedeos=1;
  //  for(qi=1;qi<=NUMINDEPDIMENS;qi++){
  for(qi=1;qi<=WHICHEOSDIMEN;qi++){ // only need to repeat used independent variables, not all
    repeatedeos*=(fabs(qarray[qi]-qoldarray[qi])<OLDTOLERANCE);
  }

  // setup old values
  if(repeatedeos==0){
    for(qi=1;qi<=NUMINDEPDIMENS;qi++) qoldarray[qi]=qarray[qi];
    // once value of independent variables changes, reset all functions as not found before
    for(qi=0;qi<NUMEOSQUANTITIESMEM;qi++) repeatedfun[qi]=0;
    gottable[0]=gottable[1]=0;
    whichtable[0]=whichtable[0]=0;
  }



  //////////////////////////////////////////////////////
  //
  // Compute result or retrieve already computed result
  //
  //////////////////////////////////////////////////////


  if(repeatedeos){

    // check if already got the function wanted
    if(repeatedfun[whichfun]){
      // choice of table irrelevant if already found solution
      *answer = resultold[whichfun];
      dualfprintf(fail_file,"REPEATALL\n"); // DEBUG
      // done!
    }
    else{
      // get result
      whichdegen=whichdegenend; // whatever is done after degen table lookup and normal table lookup
      
      if(gottable[whichdegen]==0 || gottable[1]!=gottable[0] || whichtable[1]!=whichtable[0]) return(1);// indicates "failure" and no answer within table
      else{
	dualfprintf(fail_file,"REPEATFROMLOOKUP\n"); // DEBUG
	*answer=get_eos_fromlookup(repeatedeos,WHICHEOSDIMEN,whichdegen, whichtable[0], whichfun, whichindep, vartypearray, indexarray);
	// now save result if got result
	resultold[whichfun]=*answer;
	repeatedfun[whichfun]=1;
	// done!
      }
    }

  }
  else{







    /////////////////
    //
    // assume table has low density and internal energy and don't want to go much below this since then becomes negative
    // recall that lineartablelimits is based upon utotdiff/ptotdiff/chidiff when dealing with q2
    //
    ///////////////////
    
#if(ALLOWDEGENOFFSET) // only allow floor on u if utotdiff==0 implies T=0, and assume T<0 means really T=0

#if(ALLOWSIMPLETABLE && ALLOWFULLTABLE)
    qfloor[2]   = MIN(lineartablelimits[FULLTABLE][whichindep][0],lineartablelimits[SIMPLETABLE][whichindep][0]);
#else
    qfloor[2]   = lineartablelimits[primarytable][whichindep][0];
#endif
#else
    qfloor[2] = qarray[2];// no floor on density -- if outside table (low density) will use Mignone EOS
#endif

    if(qarray[2]<qfloor[2]) qarray[2]=qfloor[2];





    ////////////////////////////
    //
    // Obtain interpolated function
    //
    ////////////////////////////
    // first get whichtable for degen table -- which determines ieos,keos,leos
    // then determine jeos after have degen offset from interpolation of ieos,keos,leos
    for(whichdegen=whichdegenstart;whichdegen>=whichdegenend;whichdegen--){


      // Once looked-up ieos,keos,leos, don't need to do full lookup again (but do need to do full interpolation again)
      // that is, only jeos changes
      // see if within table, and if so lookup interpolated result
    
#if(ALLOWDEGENOFFSET==0)
      getwhichtable=1;
#else
      //    if(whichdegen==1) getwhichtable=1;
      //    else getwhichtable=0;
      getwhichtable=1; // can't predict that when taking subtraction on utot -> utotdiff that will end up within same table
#endif
    
      if(getwhichtable) gottable[whichdegen]=iswithin_eostable(whichdegen,whichindep, vartypearray, qarray, &whichtable[whichdegen]);
      // otherwise assume same condition for gottable

      // check if got a table, and if so then use it, otherwise return(1) and assume using off-table values
      if(gottable[whichdegen]){

	// check that degen and normal table resolved to same table, since can't use degen value from one table to get values from another table
	if(whichdegen==whichdegenend){

	  if(whichtable[0]!=whichtable[1]){
	    dualfprintf(fail_file,"Degen and normal table selections different: %d %d\n",whichtable[0],whichtable[1]);
	    return(1); // revert to non-tabulated EOS
	  }
	}

	// Lookup the positions: ieos,jeos,keos,leos for a given table, independent variable definition, and indep values "qarray"
#if(ALLOWDEGENOFFSET==0)
	eos_lookup(whichtable[whichdegen], whichindep, vartypearray, qarray, indexarray);
#else
	// below also sets indexarray[2]=0 for accessing degen table
	eos_lookup_prepost_degen(whichdegen, whichtable[whichdegen], whichindep, vartypearray, qarray, indexarray);
#endif
	for(qi=WHICHEOSDIMEN+1;qi<=NUMINDEPDIMENS;qi++){
	  indexarray[qi]=0.0; // assume bottom of table is default unless changed here
	}





	// DEBUG:
	//dualfprintf(fail_file,"LOOKUPANDGET: whichdegen=%d\n",whichdegen); // DEBUG


	for(qi=1;qi<=WHICHEOSDIMEN;qi++){
	  if(indexarray[qi]<0.0){
	    dualfprintf(fail_file,"out of bounds (whichdegen=%d icurr=%d jcurr=%d kcurr=%d) for indexarray[%d]=%21.15g qarray=%21.15g\n",whichdegen,icurr,jcurr,kcurr,qi,indexarray[qi],qarray[qi]);
	  }
	}

	// now compute result
	myanswer[whichdegen]=get_eos_fromlookup(repeatedeos,WHICHEOSDIMEN,whichdegen, whichtable[whichdegen], whichfun, whichindep, vartypearray, indexarray);

	// BEGIN DEBUG
	//dualfprintf(fail_file,"q1=%21.15g q2=%21.15g q3=%21.15g q4=%21.15g\n",q1,q2,q3,q4);
	//dualfprintf(fail_file,"ieos=%21.15g jeos=%21.15g keos=%21.15g leos=%21.15g :: whichdegen=%d whichtable=%d whichfun=%d whichindep=%d :: myanswer=%21.15g\n",ieos,jeos,keos,leos,whichdegen,whichtable[whichdegen],whichfun,whichindep,myanswer[whichdegen]);
	// END DEBUG
      
	/////////////////////////
	//      
	// subtract offset from actual code value to see where within table we are in terms of the offset
	// when doing normal table this change won't matter (i.e. q2 not used again, and meaning of q2 is undefined when mixing function with independent variable, but avoid if statement)
	// if whichdegen==1, after below line then q2 contains utotoffset, ptotoffset, or chioffset consistent with independent variables used for EOSQUANTITIES




	// GODMARK: DEBUG
	if(whichfun==TEMPU){
	  dualfprintf(fail_file,"i=%d :: whichtable=%d :: ieos=%21.15g jeos=%21.15g :: answer=%21.15g\n",icurr,whichtable[whichdegen],indexarray[1],indexarray[2],myanswer[whichdegen]);
	}

	if(whichdegen==1){

	  // GODMARK: DEBUG
	  if(whichfun==TEMPU){
	    dualfprintf(fail_file,"i=%d :: ieos=%21.15g jeos=%21.15g :: q2orig=%21.15g q2new=%21.15g degen=%21.15g\n",icurr,indexarray[1],indexarray[2],qarray[2],qarray[2]-myanswer[whichdegen],myanswer[whichdegen]);
	  }
			 

	  qarray[2] -= myanswer[whichdegen];
	  if(qarray[2]<0.0){
	    // DEBUG:
	    dualfprintf(fail_file,"Got negative q2=%21.15g, forcing to be %21.15g: icurr=%d :: myanswer[%d]=%21.15g :: whichdegen=%d whichtable=%d whichfun=%d whichindep=%d\n",qarray[2],SMALL,icurr,whichdegen,myanswer[whichdegen],whichdegen,whichtable[whichdegen],whichfun,whichindep);
	    for(qi=1;qi<=NUMINDEPDIMENS;qi++){
	      dualfprintf(fail_file,"indexarray[%d]=%d\n",qi,indexarray[qi]);
	    }

#if(ALLOWDEGENOFFSET)
	    qarray[2]=1.001*lineartablelimits[primarytable][whichindep][0]; // can't be negative and use log10 interpolation, so assume very small.  This works well when using degen tables since q2~0 corresponds to T~0
	    // this way results in T~0 instead of going outside table.
#else
	    qarray[2]=SMALL; // no scale in table, so just assume bad and will likely not be within table now
#endif
	  }
	  else{
	    // DEBUG:
	    //dualfprintf(fail_file,"Got positive q2=%21.15g\n",qarray[2]);
	  
	  }
	}
      
      
	//    if(which==NUCOOL){
	//      dualfprintf(fail_file,"whichtable=%d whichindep=%d q1=%21.15g q2=%21.15g :: *answer=%21.15g\n",whichtable[whichdegen], whichindep,qarray[1],qarray[2],*answer);
	//    }



      } // end if within some table
      else return(1);// indicates "failure" and no answer within table
    }// end loop over degenerate and then normal table (or if ALLOWDEGENOFFSET==0, then only for normal table)
    
    // finally, normal table lookup gives answer
    *answer=myanswer[0];


  }






  return(0);


}





int get_whichindep(int whichfun)
{
  int whichcheck1,whichcheck2,whichcheck3;
  int whichindep;


  // GODMARK: Can make array that stores this info, looked up by whichfun as index
  // functions (F) F(rho0,u)
  // GODMARK: >=EXTRA1 assumes all extras are F(rho0,u), although for whichdatatype==4 some must be functions of \chi -- need to change if going to use that method
  whichcheck1=(whichfun==PofRHOU||whichfun==DPDRHOofRHOU||whichfun==DPDUofRHOU||whichfun==CS2ofRHOU||whichfun==SofRHOU||whichfun==DSDRHOofRHOU||whichfun==DSDUofRHOU||(whichfun>=EXTRA1)||whichfun==TEMPU);

  // functions (F) F(rho0,p)
  whichcheck2=(whichfun==UofRHOP||whichfun==TEMPP);

  // functions (F) F(rho0,\chi=u+p)
  whichcheck3=(whichfun==PofRHOCHI||whichfun==IDRHO0DP||whichfun==IDCHIDP||whichfun==TEMPCHI);


  // determine which temperature to use to check inversion
  if(whichcheck1) whichindep=UEOS;
  else if(whichcheck2) whichindep=PEOS;
  else if(whichcheck3) whichindep=CHIEOS;
  else{
    dualfprintf(fail_file,"Undefined whichfun=%d in get_whichindep()\n",whichfun);
    myexit(26867);
  }

  return(whichindep);

}








// assumes this is computed every timestep (or substep) or at least on some timescale that H changes
void compute_Hglobal(FTYPE (*prim)[N2M][N3M][NPR])
{
  int i,j,k;
  int ii,jj,kk;
  int pl;
  FTYPE r,th,phi;
  FTYPE X[NDIM],V[NDIM],dxdxp[NDIM][NDIM];
  FTYPE dr,rdth,rsinthdphi,lambdatot;
  FTYPE rho0,u;
  FTYPE Htest1,Htest2;
  int lambdatotextra;
  int hi;
  int whichfun;


  //  return; //assume initial H is good -- use to compare input and output for EOS

  // this calculation makes sense only if using spherical polar coordinates
  if(!ISSPCMCOORD(MCOORD)){
    dualfprintf(fail_file,"Hglobal not setup for anything except spherical polar coordinates: %d\n",ISSPCMCOORD(MCOORD));
    myexit(2673);
  }


  // first store lookups and computations per point
  COMPFULLLOOP{

    coord_ijk(i,j,k,CENT,X); // doesn't matter if CENT or anyother thing is used since just an estimate
    bl_coord_ijk(i,j,k,CENT,V);
    dxdxprim_ijk(i,j,k,CENT,dxdxp);
    r=V[1];
    th=V[2];
    phi=V[3];
    // assumes original coordinates are r,\theta,\phi
    dr = MAX(fabs(dxdxp[1][1]*dx[1]),fabs(dxdxp[1][2]*dx[2])); // just an estimate -- assume dxdxp[1][1]~1
    rdth = r*MAX(fabs(dxdxp[2][1]*dx[1]),fabs(dxdxp[2][2]*dx[2])); // just an estimate -- assume dxdxp[2][2]~1
#if(0)
    rsinthdphi = r*sin(th)*fabs(dxdxp[3][3]*dx[3]); // just an estimate -- assume dxdxp[3][3]~1
#endif

    rho0=prim[i][j][k][RHO];
    u=prim[i][j][k][UU];



    /////////////////////
    //
    // determine what whichfun0 means for this table
    // could put this in lookup with primarytable->whichtable, but expensive to put there and assume all tables same whichdatatype, so this is ok
    //
    /////////////////////
    if(whichdatatype[primarytable]==3){
      whichfun=EXTRA5;
    }
    else if(whichdatatype[primarytable]==4){
      whichfun=EXTRA19;
    }
    else{
      dualfprintf(fail_file,"Shouldn't request whichfun0=%d if primarytable=%d\n",LAMBDATOT,primarytable);
    }

    if(get_eos_fromtable(whichfun,rho0,u,&lambdatot)){
      lambdatot=1.0E30; // then assume optically thin (worry if outside when rho>>the limit in the table?) GODMARK
    }
    ptemparray[i][j][k][0] = lambdatot;
    ptemparray[i][j][k][1] = dr/(lambdatot+SMALL); // dtau for radial integral
    ptemparray[i][j][k][2] = rdth/(lambdatot+SMALL); // dtau for angular integral
#if(0)
    ptemparray[i][j][k][3] = rsinthddphi/(lambdatot+SMALL); // dtau for radial integral
#endif

    
  }



  //////////////////////////////
  //
  // get +- r direction H
  // GODMARK : NOT YET FOR MPI
  COMPFULLLOOP{ 
    Htest1 = Htest2 = 0.0;
    jj=j;
    kk=k;
    for(ii=i;ii<N1+N1BND;ii++){// outward pointing integral (in MPI, should include other integrals that go to r=infinity)
      Htest1 +=ptemparray[ii][jj][kk][1];
    }
#if(0)
    // convergence of neutrinos to r=0 makes no sense, so only use outgoing rays as estimate
    Htest2 = ptemparray[0][jj][kk][1]; // double count ii=0
    for(ii=0;ii<i;ii++){
      Htest2 +=ptemparray[ii][jj][kk][1];
    } 
    EOSglobal[H2GLOBAL][i][j][k] = (Htest1+Htest2*2.0)*ptemparray[i][j][k][0];
#endif
    // This is scale-height used by Kaz's EOS
    EOSglobal[HGLOBAL][i][j][k]  = Htest1*ptemparray[i][j][k][0];
    EOSglobal[H2GLOBAL][i][j][k] = Htest1*ptemparray[i][j][k][0];
  }


  //////////////////////////////
  //
  // get +- \theta direction H (assume spherical polar)
  // GODMARK : NOT YET FOR MPI
  COMPFULLLOOP{ 
    Htest1 = Htest2 = 0.0;
    ii=i;
    kk=k;
    for(jj=j;jj<N2;jj++){
      Htest1 +=ptemparray[ii][jj][kk][1];
    }
    for(jj=j;jj>=0;jj--){
      Htest2 +=ptemparray[ii][jj][kk][1];
    }
    // This is scale-height used by Kaz's EOS
    // add outgoing radial part to hack photon trajectory for +-z for disk near BH or NS
    EOSglobal[H3GLOBAL][i][j][k] = (Htest1 + EOSglobal[HGLOBAL][i][1][k]/(ptemparray[i][1][k][0]+SMALL))*ptemparray[i][j][k][0];
    EOSglobal[H4GLOBAL][i][j][k] = (Htest2 + EOSglobal[HGLOBAL][i][N2-2][k]/(ptemparray[i][N2-2][k][0]+SMALL))*ptemparray[i][j][k][0];
  }


#if(0)
  //////////////////////////////
  //
  // get +- \phi direction H (assume spherical polar)
  // GODMARK : NOT YET FOR MPI
  COMPFULLLOOP{ 
    Htest1 = Htest2 = 0.0;
    ii=i;
    jj=j;
    for(kk=k;kk<N3;kk++){
      Htest1 +=ptemparray[ii][jj][kk][1];
    }
    for(kk=k;kk>=0;kk--){
      Htest2 +=ptemparray[ii][jj][kk][1];
    }
    // This is scale-height used by Kaz's EOS
    EOSglobal[H5GLOBAL][i][j][k] = Htest1*ptemparray[i][j][k][0];
    EOSglobal[H6GLOBAL][i][j][k] = Htest2*ptemparray[i][j][k][0];
  }
#endif


  if(whichdatatype[primarytable]!=4){ // otherwise no need to limit since not using table lookup for H
    for(hi=0;hi<NUMHDIRECTIONS;hi++){
      COMPFULLLOOP{ 
	
	
	// since there exists no information beyond bounds of table, limit H to table values so return value is consistent
	// This is also needed since H is used with other EOS table values and needs to be consistent
	// Note that when using simpletable that these restrictions won't matter
	// offset from strict linear table limits because log-linear conversion leads to error
	if(EOSglobal[HGLOBAL+hi][i][j][k]>0.999*lineartablelimits[primarytable][HEOS][1]) EOSglobal[HGLOBAL+hi][i][j][k] = 0.999*lineartablelimits[primarytable][HEOS][1];
	if(EOSglobal[HGLOBAL+hi][i][j][k]<1.001*lineartablelimits[primarytable][HEOS][0]) EOSglobal[HGLOBAL+hi][i][j][k] = 1.001*lineartablelimits[primarytable][HEOS][0];
      }
    }
  }




}











// assumes this is computed every timestep (or substep) before anything else computed
void compute_TDYNORYE_YNU_global(FTYPE (*prim)[N2M][N3M][NPR])
{
  int i,j,k;
  //  int jj;
  //  int pl;
  //  FTYPE r,th;
  //  FTYPE X[NDIM],V[NDIM],dxdxp[NDIM][NDIM];
  FTYPE TDYNORYEtouse, YNUtouse;


#if(DOYL==DONOYL || DOYNU==DONOYNU)
  return;
#endif


  COMPFULLLOOP{
    
    //    coord_ijk(i,j,k,CENT,X); // doesn't matter if CENT or anyother thing is used since just an estimate
    //    bl_coord_ijk(i,j,k,CENT,V);
    //    dxdxprim_ijk(i,j,k,CENT,dxdxp);
    //    r=V[1];

    
#if(! (DOYL==DONOYL && DOYNU==DONOYNU) )
    TDYNORYEtouse = prim[i][j][k][YL] - prim[i][j][k][YNU];

    if(whichdatatype[primarytable]==4){
      // Below really YNU0
      YNUtouse = EOSglobal[YNUGLOBAL][i][j][k];
    }
    else{
      YNUtouse = prim[i][j][k][YNU];
    }
#endif  

    //    dualfprintf(fail_file,"COMPUTEYE: %21.15g %21.15g %21.15g\n",prim[i][j][k][YL],prim[i][j][k][YNU],TDYNORYEtouse);


    // since there exists no information beyond bounds of table, limit TDYNORYE to table values so return value is consistent
    // This is also needed since TDYNORYE is used with other EOS table values and needs to be consistent
    // Note that when using simpletable that these restrictions won't matter
    // offset from strict linear table limits because log-linear conversion leads to error
    if(TDYNORYEtouse>0.999*lineartablelimits[primarytable][TEOS][1]) TDYNORYEtouse = 0.999*lineartablelimits[primarytable][TEOS][1];
    if(TDYNORYEtouse<1.001*lineartablelimits[primarytable][TEOS][0]) TDYNORYEtouse = 1.001*lineartablelimits[primarytable][TEOS][0];

    if(YNUtouse>0.999*lineartablelimits[primarytable][YNUEOS][1]) YNUtouse = 0.999*lineartablelimits[primarytable][YNUEOS][1];
    if(YNUtouse<1.001*lineartablelimits[primarytable][YNUEOS][0]) YNUtouse = 1.001*lineartablelimits[primarytable][YNUEOS][0];

    // make final assignment
    EOSglobal[TDYNORYEGLOBAL][i][j][k]=TDYNORYEtouse;
    EOSglobal[YNUGLOBAL][i][j][k]=YNUtouse;        // Below really YNU0 if whichdatatype==4

  }


}




// assumes this is computed every timestep (or substep) or at least on some timescale that H changes
void compute_ups_global(FTYPE (*prim)[N2M][N3M][NPR])
{
  int get_rhops_nu(int whichd, int i, int j, int k, FTYPE *pr, FTYPE *rho_nu, FTYPE *p_nu, FTYPE *s_nu);
  FTYPE rho0,u;
  FTYPE rho_nu, p_nu, s_nu;
  int i,j,k;



  if(whichdatatype[primarytable]!=4) return; // doesn't need to be done if !=4


  COMPFULLLOOP{

    get_rhops_nu(UTOTDIFF, i, j, k, prim[i][j][k], &rho_nu, &p_nu, &s_nu);

    EOSglobal[UNUGLOBAL][i][j][k] =  rho_nu;
    EOSglobal[PNUGLOBAL][i][j][k] =  p_nu;
    EOSglobal[SNUGLOBAL][i][j][k] =  s_nu;

  }

}





// tabulated dfun[du] and need to get u=u+unu and then fun= dfun+fun_nu
FTYPE dfun2fun_kazfull(int whichfun, int whichd, FTYPE quant1, FTYPE quant2, FTYPE *dfinalreturn)
{
  FTYPE final;
  FTYPE unu,snu,pnu;
  FTYPE dfinal,dquant2;
  FTYPE quant2mod;


  if(whichdatatype[primarytable]==4){
    unu = EOSglobal[UNUGLOBAL][icurr][jcurr][kcurr];
    pnu = EOSglobal[PNUGLOBAL][icurr][jcurr][kcurr];
    snu = EOSglobal[SNUGLOBAL][icurr][jcurr][kcurr];
    if(whichd==UTOTDIFF)      dquant2  = quant2 - unu;
    else if(whichd==PTOTDIFF) dquant2  = quant2 - pnu;
    else if(whichd==CHIDIFF)  dquant2  = quant2 - (unu+pnu);
  }
  else{
    unu=snu=pnu=0.0;
    dquant2 = quant2;
  }

  if(get_eos_fromtable(whichfun,quant1,dquant2,&dfinal)){ // input quant1,dquant2 and get dfinal
    quant2mod = (quant2/quant1 - NUCLEAROFFSET)*quant1; // set nuclear per baryon offset so can smoothly connect to ideal gas EOS
    // otherwise use TM EOS
#if(REDUCE2WHICHEOS==MIGNONE)
    if(whichfun==PofRHOU)        final = pressure_rho0_u_mignone(quant1, quant2mod); // use total quant2mod
    else if(whichfun==SofRHOU)   final = compute_entropy_mignone(quant1, quant2mod);
    else if(whichfun==UofRHOP)   final = u_rho0_p_mignone(quant1, quant2mod);
    else if(whichfun==PofRHOCHI) final = pressure_wmrho0_mignone(quant1, quant2mod);
#elif(REDUCE2WHICHEOS==IDEALGAS)
    // use ideal EOS
    if(whichfun==PofRHOU)        final = pressure_rho0_u_idealgas(quant1, quant2mod); // use total quant2mod
    else if(whichfun==SofRHOU)   final = compute_entropy_idealgas(quant1, quant2mod);
    else if(whichfun==UofRHOP)   final = u_rho0_p_idealgas(quant1, quant2mod);
    else if(whichfun==PofRHOCHI) final = pressure_wmrho0_idealgas(quant1, quant2mod);    
#endif
    dfinal=final;
  }
  else{
    // relevance of neutrinos can be estimated from size of dp vs. p_nu
    if(whichfun==PofRHOU){
      final = dfinal + pnu;
      // DEBUG: GODMARK -- was on when debugging EOS
      //      dualfprintf(fail_file,"PofRHOU final=%21.15g dfinal=%21.15g pnu=%21.15g\n",final,dfinal,pnu);
    }
    else if(whichfun==SofRHOU)   final = dfinal + snu;
    else if(whichfun==UofRHOP)   final = dfinal + unu;
    else if(whichfun==PofRHOCHI) final = dfinal + pnu;
  }

  // sometimes need pure "gas" (non-neutrino) part of answer
  *dfinalreturn=dfinal;
  
  return(final);
}


// general function to interpolate between non-neutrino and neutrino values
// used for those quantities not yet setup for exactly correct answer
FTYPE fudgefrac_kazfull(int whichfun, int whichd, FTYPE quant1, FTYPE quant2)
{
  FTYPE final;
  FTYPE dquant2,quant2nu,pnu,rhonu,chinu;
  FTYPE pressure_dp_rho0_u_kazfull(FTYPE rho0, FTYPE u, FTYPE *dp);
  FTYPE pressure_dp_wmrho0_kazfull(FTYPE rho0, FTYPE wmrho0, FTYPE *dp);
  FTYPE nonneutrino,neutrino;
  FTYPE total;
  FTYPE ptot;
  FTYPE frac;
  FTYPE quant2mod,dquant2mod;
  FTYPE pgas, chigas, chitot;



  // first get total pressure in order to frac-fudge the answer
  if(whichd==UTOTDIFF){
    ptot = pressure_dp_rho0_u_kazfull(quant1, quant2, &pgas); // need pgas to form chi since otherwise don't have it
    chigas = quant2 + pgas;
  }
  else if(whichd==CHIDIFF){
    ptot =pressure_dp_wmrho0_kazfull(quant1, quant2, &pgas); // don't need separate pgas
    chigas = quant2;
  }
  else{
    dualfprintf(fail_file,"fudgefrac_kazfull() not setup for whichd=%d\n",whichd);
    myexit(19672606);
  }


  /////////////////////////////////////
  //
  // now get neutrino quantities alone
  //
  /////////////////////////////////////
  if(whichdatatype[primarytable]==4){

    pnu = EOSglobal[PNUGLOBAL][icurr][jcurr][kcurr];
    rhonu = EOSglobal[UNUGLOBAL][icurr][jcurr][kcurr];
    chinu = rhonu + pnu;

    if(whichd==UTOTDIFF)      quant2nu = rhonu;
    else if(whichd==PTOTDIFF) quant2nu = pnu;
    else if(whichd==CHIDIFF)  quant2nu = chinu;

    dquant2 = quant2 - quant2nu;

  }
  else{
    dquant2 = quant2;
    quant2nu = chinu = pnu = rhonu = 0.0;
  }


  ///////////////
  //
  // set nuclear per baryon offset so can smoothly connect to ideal gas form of EOS
  //
  ///////////////
  dquant2mod = (quant2/quant1 - NUCLEAROFFSET)*quant1 - quant2nu;

  // set nuclear per baryon offset so can smoothly connect to ideal gas EOS
  quant2mod = (quant2/quant1 - NUCLEAROFFSET)*quant1;



  ////////////////////////////////////////
  //
  // determine fake-neutrino part
  //
  ////////////////////////////////////////
#if(REDUCE2WHICHEOS==MIGNONE)
  // GODMARK: Assume at least that if neutrino pressure is dominant then treat like mixture of gas+neutrinos
  if(whichfun==CS2ofRHOU)         neutrino = cs2_compute_mignone(quant1, quant2nu);
  else if(whichfun==DPDUofRHOU)   neutrino = dpdu_rho0_u_mignone(quant1, quant2nu);
  else if(whichfun==DPDRHOofRHOU) neutrino = dpdrho0_rho0_u_mignone(quant1, quant2nu);
  else if(whichfun==IDRHO0DP)     neutrino = compute_idrho0dp_mignone(quant1, quant2nu);
  else if(whichfun==IDCHIDP)      neutrino = compute_idwmrho0dp_mignone(quant1, quant2nu);
  else if(whichfun==DSDRHOofRHOU) neutrino = compute_dSdrho_mignone(quant1, quant2nu);
  else if(whichfun==DSDUofRHOU)   neutrino = compute_dSdu_mignone(quant1, quant2nu);
#elif(REDUCE2WHICHEOS==IDEALGAS)
  if(whichfun==CS2ofRHOU)         neutrino = cs2_compute_idealgas(quant1, quant2nu);
  else if(whichfun==DPDUofRHOU)   neutrino = dpdu_rho0_u_idealgas(quant1, quant2nu);
  else if(whichfun==DPDRHOofRHOU) neutrino = dpdrho0_rho0_u_idealgas(quant1, quant2nu);
  else if(whichfun==IDRHO0DP)     neutrino = compute_idrho0dp_idealgas(quant1, quant2nu);
  else if(whichfun==IDCHIDP)      neutrino = compute_idwmrho0dp_idealgas(quant1, quant2nu);
  else if(whichfun==DSDRHOofRHOU) neutrino = compute_dSdrho_idealgas(quant1, quant2nu);
  else if(whichfun==DSDUofRHOU)   neutrino = compute_dSdu_idealgas(quant1, quant2nu);
#endif

  // reduce to using totals if not within table
#if(REDUCE2WHICHEOS==MIGNONE)
  // GODMARK: Assume at least that if neutrino pressure is dominant then treat like mixture of gas+neutrinos
  if(whichfun==CS2ofRHOU)         total = cs2_compute_mignone(quant1, quant2mod);
  else if(whichfun==DPDUofRHOU)   total = dpdu_rho0_u_mignone(quant1, quant2mod);
  else if(whichfun==DPDRHOofRHOU) total = dpdrho0_rho0_u_mignone(quant1, quant2mod);
  else if(whichfun==IDRHO0DP)     total = compute_idrho0dp_mignone(quant1, quant2mod);
  else if(whichfun==IDCHIDP)      total = compute_idwmrho0dp_mignone(quant1, quant2mod);
  else if(whichfun==DSDRHOofRHOU) total = compute_dSdrho_mignone(quant1, quant2mod);
  else if(whichfun==DSDUofRHOU)   total = compute_dSdu_mignone(quant1, quant2mod);
#elif(REDUCE2WHICHEOS==IDEALGAS)
  if(whichfun==CS2ofRHOU)         total = cs2_compute_idealgas(quant1, quant2mod);
  else if(whichfun==DPDUofRHOU)   total = dpdu_rho0_u_idealgas(quant1, quant2mod);
  else if(whichfun==DPDRHOofRHOU) total = dpdrho0_rho0_u_idealgas(quant1, quant2mod);
  else if(whichfun==IDRHO0DP)     total = compute_idrho0dp_idealgas(quant1, quant2mod);
  else if(whichfun==IDCHIDP)      total = compute_idwmrho0dp_idealgas(quant1, quant2mod);
  else if(whichfun==DSDRHOofRHOU) total = compute_dSdrho_idealgas(quant1, quant2mod);
  else if(whichfun==DSDUofRHOU)   total = compute_dSdu_idealgas(quant1, quant2mod);
#endif


  /////////////////////////////////
  //
  // now get non-neutrino value
  // nonneutrino will then be only "gas" part without neutrinos
  //
  /////////////////////////////////
  if(get_eos_fromtable(whichfun,quant1,dquant2,&nonneutrino)){
    final=total;
  }
  else{
    // then estimate importance of neutrinos and its control over sound speed
    if(whichfun==CS2ofRHOU){
      // for sound speed more accurate to do:
      // c_s^2[total] = c_s^2[gas]*(hgas/htotal) + c_s^2[neutrino]*(hneutrino/htotal)
      // hgas = (rho_0 + ugas + pgas)/rho_0
      // hneutrino = (rho_0 + rhonu + pnu)/rho_0
      chitot = chigas + chinu;
      final = (nonneutrino*(quant1+chigas) + neutrino*(quant1+chinu))/(quant1+chitot);
    }
    else{
      //    frac = (pnu/(ptot-pnu));
      frac = pnu/ptot; // -> 1 if ptot=pnu.  Suppose photons+electrons+neutrinos equally dominate.  Then cs2 still correct.
      // choose cs2total if neutrino-dominated, otherwise choose non-neutrino cs2
      final = total*frac + nonneutrino*(1.0-frac);
    }

  }

  return(final);

}






// p(rho0, u) (needed to get initial guess for W)
FTYPE pressure_rho0_u_kazfull(FTYPE rho0, FTYPE u)
{
  FTYPE dp;
  FTYPE dfun2fun_kazfull(int whichfun, int whichd, FTYPE quant1, FTYPE quant2, FTYPE *dfinalreturn);
  return(dfun2fun_kazfull(PofRHOU, UTOTDIFF, rho0, u, &dp));
}

// p(rho0, u) (used internally)
FTYPE pressure_dp_rho0_u_kazfull(FTYPE rho0, FTYPE u, FTYPE *dp)
{
  FTYPE dfun2fun_kazfull(int whichfun, int whichd, FTYPE quant1, FTYPE quant2, FTYPE *dfinalreturn);
  return(dfun2fun_kazfull(PofRHOU, UTOTDIFF, rho0, u, dp));
}

// exactly correct answer (not frac-fudged)
// entropy as function of rho0 and internal energy (u)
// S(rho0,u)
// tabulated ds(du), so first compute du and then ds and then s
FTYPE compute_entropy_kazfull(FTYPE rho0, FTYPE u)
{
  FTYPE ds;
  FTYPE dfun2fun_kazfull(int whichfun, int whichd, FTYPE quant1, FTYPE quant2, FTYPE *dfinalreturn);
  return(dfun2fun_kazfull(SofRHOU, UTOTDIFF, rho0, u, &ds));
}


// u(rho0, p) (used for initial conditions)
FTYPE u_rho0_p_kazfull(FTYPE rho0, FTYPE p)
{
  FTYPE du;
  FTYPE dfun2fun_kazfull(int whichfun, int whichd, FTYPE quant1, FTYPE quant2, FTYPE *dfinalreturn);
  return(dfun2fun_kazfull(UofRHOP, PTOTDIFF, rho0, p, &du));
}


// p(rho0, w-rho0 = u+p)
// Notice that using EOSglobal bypasses need to have other quantities as functions of wmrho0 unless want more direct (by iteration) result compared to what EOSglobal gives
FTYPE pressure_wmrho0_kazfull(FTYPE rho0, FTYPE wmrho0)
{
  FTYPE dp;
  FTYPE dfun2fun_kazfull(int whichfun, int whichd, FTYPE quant1, FTYPE quant2, FTYPE *dfinalreturn);
  return(dfun2fun_kazfull(PofRHOCHI, CHIDIFF, rho0, wmrho0, &dp));
}

// used internally
FTYPE pressure_dp_wmrho0_kazfull(FTYPE rho0, FTYPE wmrho0, FTYPE *dp)
{
  FTYPE dfun2fun_kazfull(int whichfun, int whichd, FTYPE quant1, FTYPE quant2, FTYPE *dfinalreturn);
  return(dfun2fun_kazfull(PofRHOCHI, CHIDIFF, rho0, wmrho0, dp));
}






// frac-fudged because dpdu is complicated with p=dp+pnu and u=du+unu
FTYPE dpdu_rho0_u_kazfull(FTYPE rho0, FTYPE u)
{
  FTYPE fudgefrac_kazfull(int whichfun, int whichd, FTYPE quant1, FTYPE quant2);
  return(fudgefrac_kazfull(DPDRHOofRHOU,UTOTDIFF,rho0, u));
}


// for this could store dp_\nu/drho0, which for whichdatatype==4 is 0 given other independent variables are rho0,T,Y_e,Y_\nu,H
// for whichdatatype==4 check that dpnu/drho0 = 0 for whichdatatype==4 if ever want to use that method (not now)
// dp(rho0, u)/drho0
// GODMARK: frac-fudged (but perhaps shouldn't)
FTYPE dpdrho0_rho0_u_kazfull(FTYPE rho0, FTYPE u)
{
  FTYPE fudgefrac_kazfull(int whichfun, int whichd, FTYPE quant1, FTYPE quant2);
  return(fudgefrac_kazfull(DPDUofRHOU,UTOTDIFF,rho0, u));
}


// sound speed squared (for vchar.c) -- important for treatment of shocks
// frac-fudged
FTYPE cs2_compute_kazfull(FTYPE rho0, FTYPE u)
{
  FTYPE fudgefrac_kazfull(int whichfun, int whichd, FTYPE quant1, FTYPE quant2);
  return(fudgefrac_kazfull(CS2ofRHOU,UTOTDIFF,rho0, u));
}


// used for dudp_calc
// frac-fudged
FTYPE compute_dSdrho_kazfull(FTYPE rho0, FTYPE u)
{
  FTYPE fudgefrac_kazfull(int whichfun, int whichd, FTYPE rho0, FTYPE u);
  return(fudgefrac_kazfull(DSDRHOofRHOU,UTOTDIFF,rho0, u));
}


// used for dudp_calc
// frac-fudged
FTYPE compute_dSdu_kazfull(FTYPE rho0, FTYPE u)
{
  FTYPE fudgefrac_kazfull(int whichfun, int whichd, FTYPE rho0, FTYPE u);
  return(fudgefrac_kazfull(DSDUofRHOU,UTOTDIFF, rho0, u));
}


// u(rho0,S)
// here entropy is entropy density?
// not needed for now
// GODMARK: not corrected for ds -- well, doesn't even use table
FTYPE compute_u_from_entropy_kazfull(FTYPE rho0, FTYPE entropy)
{
  FTYPE u;
  FTYPE compute_u_from_entropy_mignone(FTYPE rho0, FTYPE entropy);

  // GODMARK: no kaz solution yet since not needed yet
  u=compute_u_from_entropy_mignone(rho0, entropy);

  // if(get_eos_fromtable(UofRHOS,rho0,entropy,Hglobal,TDYNORYEglobal, &u)){
  //  u=0.0; // GODMARK: not set yet
  //}

  return(u);

}


// 1 / (drho0/dp) holding wmrho0 fixed
// frac-fudged
FTYPE compute_idrho0dp_kazfull(FTYPE rho0, FTYPE wmrho0)
{
  FTYPE fudgefrac_kazfull(int whichfun, int whichd, FTYPE rho0, FTYPE u);
  return(fudgefrac_kazfull(IDRHO0DP,CHIDIFF, rho0, wmrho0));
}



// 1 / (d(u+p)/dp) holding rho0 fixed
// frac-fudged
FTYPE compute_idwmrho0dp_kazfull(FTYPE rho0, FTYPE wmrho0)
{
  FTYPE fudgefrac_kazfull(int whichfun, int whichd, FTYPE rho0, FTYPE u);
  return(fudgefrac_kazfull(IDRHO0DP,CHIDIFF, rho0, wmrho0));
}


// used for fully tabulated quantities that are functions of du/dp/dchi with whichdatatype==4
FTYPE compute_tabulated_kazfull(int whichfun, int whichd, FTYPE rho0, FTYPE quant2)
{
  FTYPE final;
  FTYPE dquant2;


  if(whichdatatype[primarytable]==4){
    if(whichd==UTOTDIFF) dquant2 = quant2 - EOSglobal[UNUGLOBAL][icurr][jcurr][kcurr];
    else if(whichd==PTOTDIFF) dquant2 = quant2 - EOSglobal[PNUGLOBAL][icurr][jcurr][kcurr];
    else if(whichd==CHIDIFF) dquant2 = quant2 - (EOSglobal[UNUGLOBAL][icurr][jcurr][kcurr]+EOSglobal[PNUGLOBAL][icurr][jcurr][kcurr]);
  }
  else{
    dquant2 = quant2;
  }

  if(get_eos_fromtable(whichfun,rho0,dquant2,&final)){ // uses dquant2
    if(whichfun==TEMPU || whichfun==TEMPP || whichfun==TEMPCHI) final=0.0;
    else if(whichfun==QDOTNU) final=0.0;
    else final=0.0; // use mignone?
  }  

  return(final);

}


// volume heating rate(rho0,u)
FTYPE compute_qdot_kazfull(FTYPE rho0, FTYPE u)
{
  FTYPE compute_tabulated_kazfull(int whichfun, int whichd, FTYPE rho0, FTYPE u);
  int whichfun;

  if(whichdatatype[primarytable]==1){
    whichfun=EXTRA1;
  }
  else if(whichdatatype[primarytable]==3){
    whichfun=EXTRA1;
  }
  else{
    dualfprintf(fail_file,"Shouldn't request whichfun0=%d if primarytable=%d\n",QDOTNU,primarytable);
  }

  return(compute_tabulated_kazfull(whichfun, UTOTDIFF, rho0, u));
}


// temperature(rho0,u)
FTYPE compute_temp_kazfull(FTYPE rho0, FTYPE u)
{
  FTYPE compute_tabulated_kazfull(int whichfun, int whichd, FTYPE rho0, FTYPE u);
  return(compute_tabulated_kazfull(TEMPU, UTOTDIFF, rho0, u));
}


// temperature (direct lookup from differential quant2: dquant2)
// whichfun = TEMPU + 0,1,2 for TEMPU,TEMPP,TEMPCHI
FTYPE compute_temp_whichd_kazfull(int whichd, FTYPE rho0, FTYPE dquant2)
{
  FTYPE temp;

  if(get_eos_fromtable(TEMPU+whichd,rho0,dquant2,&temp)){
    temp=0.0;
  }

  return(temp);



}


// uses global variable "numextras" that should be set before this function is called
// f(rho0,du)
// direct lookup
void compute_allextras_du_kazfull(int justnum, FTYPE rho0, FTYPE du, int *numextrasreturn, FTYPE *extras)
{
  int i;

  if(justnum==0){

    *numextrasreturn=numextras[primarytable]; // by default assume tablulated data exists


    // assume all tables have same number of extras or else this doesn't make sense in general
    for(i=0;i<numextras[primarytable];i++){
      if(get_eos_fromtable(EXTRA1+i,rho0,du,&(extras[i]))){
	extras[i]=0.0;
	*numextrasreturn=0; // tells calling function that no tabulated extras
      }
    }
    // set rest to 0
    for(i=numextras[primarytable];i<MAXNUMEXTRAS;i++){
      extras[i] = 0.0;
    }
  }

}


// uses global variable "numextras" that should be set before this function is called
// f(rho0,u)
void compute_allextras_kazfull(int justnum, FTYPE rho0, FTYPE u, int *numextrasreturn, FTYPE *extras)
{
  void compute_allextras_du_kazfull(int justnum, FTYPE rho0, FTYPE du, int *numextrasreturn, FTYPE *extras);
  int i;
  FTYPE du;

  if(whichdatatype[primarytable]==4){
    du = u - EOSglobal[UNUGLOBAL][icurr][jcurr][kcurr];
  }
  else{
    du = u;
  }

  compute_allextras_du_kazfull(justnum, rho0, du, numextrasreturn, extras);

}










// f2c prototype
#include "f2c.h"
#include "tau_neededbyharm.P"
// not linking with libf2c since don't want that dependence and conversion doesn't need it since the original code was simple

// get neutrino terms (and return extras also)
// here function is (rho0,u) -- never used as function of p or chi
int get_extrasprocessed_kazfull(int doall, int i, int j, int k, FTYPE *pr, FTYPE *extras, FTYPE *processed)
{
  FTYPE quant1,quant2;
  // assume all tables have same number of extras or else this doesn't make sense in general
  FTYPE tempk,kbtk,rho0,u,H,ye,ynu;
  int numextrasreturn,ei;
  void compute_allextras_du_kazfull(int justnum, FTYPE rho0, FTYPE du, int *numextrasreturn, FTYPE *extras);
  FTYPE compute_temp_kazfull(FTYPE rho0, FTYPE u);
  FTYPE unue0,unuebar0,unumu0,
    qtautnueohcm,qtautnuebarohcm,qtautmuohcm,
    qtauanueohcm,qtauanuebarohcm,qtauamuohcm,
    nnue0,nnuebar0,nnumu0,
    ntautnueohcm,ntautnuebarohcm,ntautmuohcm,
    ntauanueohcm,ntauanuebarohcm,ntauamuohcm,
    lambdatot,lambdaintot,
    tauphotonohcm,tauphotonabsohcm,
    nnueth0,nnuebarth0;
  FTYPE qphoton,qm,graddotrhouyl,tthermaltot,tdifftot,rho_nu,p_nu,s_nu,ynulocal,Ynuthermal,enu,enue,enuebar;
  FTYPE qphoton_a[NUMHDIRECTIONS],qm_a[NUMHDIRECTIONS],graddotrhouyl_a[NUMHDIRECTIONS],tthermaltot_a[NUMHDIRECTIONS],tdifftot_a[NUMHDIRECTIONS],rho_nu_a[NUMHDIRECTIONS],p_nu_a[NUMHDIRECTIONS],s_nu_a[NUMHDIRECTIONS],ynulocal_a[NUMHDIRECTIONS],Ynuthermal_a[NUMHDIRECTIONS],enu_a[NUMHDIRECTIONS],enue_a[NUMHDIRECTIONS],enuebar_a[NUMHDIRECTIONS];
  FTYPE dquant2;
  FTYPE qarray[NUMINDEPDIMENS+1];
  static FTYPE extrasold[MAXNUMEXTRAS];
  static FTYPE processedold[MAXPROCESSEDEXTRAS];
  static FTYPE qoldarray[NUMINDEPDIMENS+1];
  int repeatedeos;
  int qi;
  int hi;
  FTYPE frac;
  int notintable;
  int whichd;



  
  if(whichdatatype[primarytable]!=4){
    dualfprintf(fail_file,"Neutrino calculation not setup for other datatypes=%d primarytable=%d\n",whichdatatype[primarytable],primarytable);
    myexit(2954654);
  }

  whichd=UTOTDIFF; // only one allowed for now
  // for EOS functions
  icurr=i;
  jcurr=j;
  kcurr=k;

  quant1=pr[RHO];
  quant2=pr[UU];

  // first get dquant2 from quant2
  // avoid getting densities effectively twice, so just get UNUGLOBAL directly
  // Get dquant2 = u - u_\nu
  dquant2 = quant2 - EOSglobal[UNUGLOBAL][i][j][k];


  // setup array
  qarray[1]=quant1;
  qarray[2]=dquant2;
  // qarray[3+] are always stored in EOSglobal
  for(qi=3;qi<=NUMINDEPDIMENS;qi++){
    qarray[qi] = EOSglobal[qi][i][j][k];
  }


  // check if repeated case (doesn't matter if i,j,k same since result only depends on all qarray values)
  repeatedeos=1;
  for(qi=1;qi<=WHICHEOSDIMEN;qi++){ // only need to repeat used independent variables, not all
    repeatedeos*=(fabs(qarray[qi]-qoldarray[qi])<OLDTOLERANCE);
  }



  if(repeatedeos){
    if(doall){
      // then repeated case, so just return old result
      for(ei=0;ei<MAXNUMEXTRAS;ei++) extras[ei]=extrasold[ei];
      for(ei=0;ei<MAXPROCESSEDEXTRAS;ei++) processed[ei]=processedold[ei];
    }
    else{
      processed[RHONU]=processedold[RHONU];
      processed[PNU]=processedold[PNU];
      processed[SNU]=processedold[SNU];
    }    
  }
  else{

    // set extra parameters
    rho0=qarray[1];
    u=qarray[2]; // "u" not used
    ye=qarray[3]; // "ye" not used
    ynu=pr[YNU]; // Ynu[orig] = Ynu
    H=qarray[5];

    // not normal temperature(rho0,u), but temp(rho0,d), so direct lookup as tabulated
    // if used temp(rho0,u) then would get densities and this is already being done here, so avoid extra iteration
    // tempk is always stored in dimensionless form of: T[k] k_b/(m_b c^2), so since need kbtk in energy units, then always multiply tempk by m_b c^2 no matter what rho0unittype is
    tempk=compute_temp_whichd_kazfull(whichd, quant1, dquant2);
    kbtk=tempk*mbcsq; // in code energy units


    if(doall){
      // get extras (function of rho0,du as tabulated)
      compute_allextras_du_kazfull(0, quant1, dquant2, &numextrasreturn, extras);


      if(numextrasreturn!=0){
	notintable=0;
	// MAXNUMEXTRAS entries
	// not same order as passed to function: computefinal_fromhcm()
	ei=0;
	qtautnueohcm=extras[ei]; ei++;
	qtauanueohcm=extras[ei]; ei++;
	qtautnuebarohcm=extras[ei]; ei++;
	qtauanuebarohcm=extras[ei]; ei++;
	qtautmuohcm=extras[ei]; ei++;
	qtauamuohcm=extras[ei]; ei++;

	ntautnueohcm=extras[ei]; ei++;
	ntauanueohcm=extras[ei]; ei++;
	ntautnuebarohcm=extras[ei]; ei++;
	ntauanuebarohcm=extras[ei]; ei++;
	ntautmuohcm=extras[ei]; ei++;
	ntauamuohcm=extras[ei]; ei++;

	unue0=extras[ei]; ei++;
	unuebar0=extras[ei]; ei++;
	unumu0=extras[ei]; ei++;

	nnue0=extras[ei]; ei++;
	nnuebar0=extras[ei]; ei++;
	nnumu0=extras[ei]; ei++;

	lambdatot=extras[ei]; ei++; // EXTRA19
	lambdaintot=extras[ei]; ei++; // EXTRA20

	tauphotonohcm=extras[ei]; ei++;
	tauphotonabsohcm=extras[ei]; ei++;

	nnueth0=extras[ei]; ei++;
	nnuebarth0=extras[ei]; ei++;
      }
      else{
	notintable=1;
	// then not in table, so revert to estimated non-tabulated values (i.e. optically thin)
	// extra quantities were assumed set to 0 if not in table, so only need to change non-zero things
	// optical depths and densities are complicated and if out of table then they aren't determined easily
	lambdatot=lambdaintot=1.0E30; // optically thin
      }

    }
    else{
      notintable=0;
      // get extras (only those needed)
      if(get_eos_fromtable(EXTRA1,quant1,dquant2,&qtautnueohcm)){
	qtautnueohcm=0.0; // not in tables, assume optically thin
	notintable=1;
      }
      if(!notintable){
	// EXTRA1:  qtautnueohcm
	// EXTRA2:  qtauanueohcm
	// EXTRA3:  qtautnuebarohcm
	// EXTRA4:  qtauanuebarohcm
	// EXTRA5:  qtautmuohcm
	// EXTRA6:  qtauamuohcm
	// EXTRA13:  unue0
	// EXTRA14:  unuebar0
	// EXTRA15:  unumu0
	get_eos_fromtable(EXTRA2,quant1,dquant2,&qtauanueohcm);
	get_eos_fromtable(EXTRA3,quant1,dquant2,&qtautnuebarohcm);
	get_eos_fromtable(EXTRA4,quant1,dquant2,&qtauanuebarohcm);
	get_eos_fromtable(EXTRA5,quant1,dquant2,&qtautmuohcm);
	get_eos_fromtable(EXTRA6,quant1,dquant2,&qtauamuohcm);
	get_eos_fromtable(EXTRA13,quant1,dquant2,&unue0);
	get_eos_fromtable(EXTRA14,quant1,dquant2,&unuebar0);
	get_eos_fromtable(EXTRA15,quant1,dquant2,&unumu0);
      }
      else{
	// not in table, then non-zero quantities need to be set to something else
	// so far choosing all to be 0.0 is correct for optically thin approximation to off-table quantities
      }
    }




    for(hi=0;hi<NUMHDIRECTIONS;hi++){
      // now process neutrino variables into final form
      H = EOSglobal[HGLOBAL+hi][i][j][k];

      if(doall){

	if(!notintable){
	  // calling f2c generated code, which assumes a certain dimensionless form for inputs
	  // rhob/mb should be a number density, and mb*rate/cc should be in rhounit*rate
	  // this means if rho0 in c^2 units, then mb should be too
	  computefinal_fromhcm__(&Ccode,&mbwithrhounit,&rho0,&kbtk,&H,
				 &unue0,&unuebar0,&unumu0,
				 &qtautnueohcm,&qtautnuebarohcm,&qtautmuohcm,
				 &qtauanueohcm,&qtauanuebarohcm,&qtauamuohcm,
				 &nnue0,&nnuebar0,&nnumu0,
				 &ntautnueohcm,&ntautnuebarohcm,&ntautmuohcm,
				 &ntauanueohcm,&ntauanuebarohcm,&ntauamuohcm,
				 &lambdatot,&lambdaintot,
				 &tauphotonohcm,&tauphotonabsohcm,
				 &nnueth0,&nnuebarth0,
				 // outputs below
				 &qphoton_a[hi],&qm_a[hi],&graddotrhouyl_a[hi],&tthermaltot_a[hi],&tdifftot_a[hi],&rho_nu_a[hi],&p_nu_a[hi],&s_nu_a[hi],&ynulocal_a[hi],&Ynuthermal_a[hi],&enu_a[hi],&enue_a[hi],&enuebar_a[hi]);
	}
	else{
	  // then set to approximate optically thin values (i.e. ~0)
	  qphoton_a[hi]=qm_a[hi]=graddotrhouyl_a[hi]=rho_nu_a[hi]=p_nu_a[hi]=s_nu_a[hi]=ynulocal_a[hi]=Ynuthermal_a[hi]=enu_a[hi]=enue_a[hi]=enuebar_a[hi]=0.0;
	  tthermaltot_a[hi]=tdifftot_a[hi]=1E50; // force non-evolution of Ynu
	}
      }
      else{
	// DEBUG:
	dualfprintf(fail_file,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",qtautnueohcm,qtauanueohcm,qtautnuebarohcm,qtauanuebarohcm,qtautmuohcm,qtauamuohcm,unue0,unuebar0,unumu0);
	dualfprintf(fail_file,"CCODE: icurr=%d :: %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",icurr,Ccode,mbwithrhounit,rho0,kbtk,H,unue0,unuebar0,unumu0,qtautnueohcm,qtautnuebarohcm,qtautmuohcm,qtauanueohcm,qtauanuebarohcm,qtauamuohcm);
	dualfprintf(fail_file,"mb[cgs]=%21.15g T[K] = %21.15g\n",mbwithrhounit*Munit/Vunit/Vunit,kbtk*energyunit/kb);

	if(!notintable){
	  // get rho_nu, p_nu, s_nu from extras
	  computefinal_justdensities_fromhcm__(&Ccode,&mbwithrhounit,
					       &rho0,&kbtk,&H,
					       &unue0,&unuebar0,&unumu0,
					       &qtautnueohcm,&qtautnuebarohcm,&qtautmuohcm,
					       &qtauanueohcm,&qtauanuebarohcm,&qtauamuohcm,
					       // outputs below
					       &rho_nu_a[hi],&p_nu_a[hi],&s_nu_a[hi]); // notice no & since already pointer
	}
	else{
	  // not in table
	  rho_nu_a[hi]=p_nu_a[hi]=s_nu_a[hi]=0.0;
	}
      }
    }
    

    // now get appropriate sum
    // if all h's same then get normal result, otherwise optical depth supression occurs per direction
    // hack for real method
    // GODMARK: a heating method might want access to these different directional quantities
    qphoton=qm=graddotrhouyl=tthermaltot=tdifftot=rho_nu=p_nu=s_nu=ynulocal=Ynuthermal=enu=enue=enuebar=0.0;
    frac=1.0/((FTYPE)NUMHDIRECTIONS);
    for(hi=0;hi<NUMHDIRECTIONS;hi++){

      if(doall){
	qphoton += qphoton_a[hi]*frac;
	qm += qm_a[hi]*frac;
	graddotrhouyl += graddotrhouyl_a[hi]*frac;
	tthermaltot += tthermaltot_a[hi]*frac;
	tdifftot += tdifftot_a[hi]*frac;
	ynulocal += ynulocal_a[hi]*frac;
	Ynuthermal += Ynuthermal_a[hi]*frac;
	enu += enu_a[hi]*frac;
	enue += enue_a[hi]*frac;
	enuebar += enuebar_a[hi]*frac;
      }
      
      rho_nu += rho_nu_a[hi]*frac;
      p_nu += p_nu_a[hi]*frac;
      s_nu += s_nu_a[hi]*frac;

    }


    if(doall){
      // MAXPROCESSEDEXTRAS entries
      ei=-1; // just so easy copy/paste form
      ei++; processed[ei]=qphoton;
      ei++; processed[ei]=qm;
      ei++; processed[ei]=graddotrhouyl;
      ei++; processed[ei]=tthermaltot;
      ei++; processed[ei]=tdifftot;
      ei++; processed[ei]=rho_nu;
      ei++; processed[ei]=p_nu;
      ei++; processed[ei]=s_nu;
      ei++; processed[ei]=ynulocal;
      ei++; processed[ei]=Ynuthermal;
      // below are energies of *escaping* neutrinos
      ei++; processed[ei]=enu;
      ei++; processed[ei]=enue;
      ei++; processed[ei]=enuebar;
    }
    else{
      processed[RHONU]=rho_nu;
      processed[PNU]=p_nu;
      processed[SNU]=s_nu;
    }

    //  MAXPROCESSEDEXTRAS entries
    // below not used right now


    // setup old values
    for(qi=1;qi<=NUMINDEPDIMENS;qi++) qoldarray[qi]=qarray[qi];
    if(doall){
      for(ei=0;ei<MAXNUMEXTRAS;ei++) extrasold[ei]=extras[ei];
      for(ei=0;ei<MAXPROCESSEDEXTRAS;ei++) processedold[ei]=processed[ei];
    }
    else{
      processedold[RHONU]=processed[RHONU];
      processedold[PNU]=processed[PNU];
      processedold[SNU]=processed[SNU];
    }


    // assume DO WANT to update EOSglobal[UNUGLOBAL,PNUGLOBAL,SNUGLOBAL][i][j][k] since input is actual value on grid
    EOSglobal[UNUGLOBAL][i][j][k] = rho_nu;
    EOSglobal[PNUGLOBAL][i][j][k] = p_nu;
    EOSglobal[SNUGLOBAL][i][j][k] = s_nu;
    // Will use below to compute Ynu0 = Ynu[orig] + dYnu, where dYnu = -Ynu[Ynu0,rho,du,Ye,H] + Ynu0
    //    EOSglobal[YNUGLOBAL][i][j][k] = ynu + (-ynulocal+EOSglobal[YNUGLOBAL][i][j][k]);
    // iterate Ynu0 for table lookup
    EOSglobal[YNUGLOBAL][i][j][k] += (ynu - ynulocal);

    // need to limit Ynu0 to table (assumes reasonable behavior for above iteration!)
    if(EOSglobal[YNUGLOBAL][i][j][k]>0.999*lineartablelimits[primarytable][YNUEOS][1]) EOSglobal[YNUGLOBAL][i][j][k] = 0.999*lineartablelimits[primarytable][YNUEOS][1];
    if(EOSglobal[YNUGLOBAL][i][j][k]<1.001*lineartablelimits[primarytable][YNUEOS][0]) EOSglobal[YNUGLOBAL][i][j][k] = 1.001*lineartablelimits[primarytable][YNUEOS][0];


  }



  return(0);


}






// whichd = UTOTDIFF,PTOTDIFF,CHIDIFF
// get neutrino rho_nu, p_nu, s_nu
// function is (rho0,u/p/chi)
// for now assume only needed as function of quant2=u.  Assumes this function called infrequently outside iterative routines
// Point is that anyways we are doing an approximation, and save memory and gain simplicity by avoiding tabulating other f(rho0,chi) quantities
int get_rhops_nu(int whichd, int i, int j, int k, FTYPE *pr, FTYPE *rho_nu, FTYPE *p_nu, FTYPE *s_nu)
{
  int get_extrasprocessed_kazfull(int doall, int i, int j, int k, FTYPE *pr, FTYPE *extras, FTYPE *processed);
  int toreturn;
  FTYPE extras[MAXNUMEXTRAS];
  FTYPE processed[MAXPROCESSEDEXTRAS];


  icurr=i;
  jcurr=j;
  kcurr=k;

  // assume WHICHD==UTOTDIFF only
  toreturn=get_extrasprocessed_kazfull(0, i, j, k, pr, extras, processed);
  *rho_nu=processed[RHONU];
  *p_nu=processed[PNU];
  *s_nu=processed[SNU];

  return(toreturn);
}









// compute sources to equations of motion for HARM
// duplicate of compute_neutrino() in sources.c
// Ui and dUother in UNOTHING form
int compute_sources_EOS_kazfull(FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *Ui, FTYPE *dUother, FTYPE(*dUcomp)[NPR])
{
  int ii,jj,kk,pp;
  int j;
  FTYPE X[NDIM],V[NDIM],r,th,R;
  FTYPE du;
  FTYPE rho,u,yl,ynu;
  FTYPE cofactor;
  FTYPE dUdtau;
  FTYPE rph,photoncapture;
  FTYPE extras[MAXNUMEXTRAS];
  void compute_allextras_kazfull(int justnum, FTYPE rho0, FTYPE u, int *numextrasreturn, FTYPE *extras);
  int get_extrasprocessed_kazfull(int doall, int i, int j, int k, FTYPE *pr, FTYPE *extras, FTYPE *processed);
  FTYPE processed[MAXPROCESSEDEXTRAS];
  int ei;
  FTYPE qphoton,qm,graddotrhouyl,tthermaltot,tdifftot,rho_nu,p_nu,s_nu,ynulocal,enu,enue,enuebar;
  FTYPE Ynuthermal,lambdatot,lambdaintot;
  FTYPE graddotrhouynu;
  FTYPE dtau,dtlimit;
  FTYPE gradUU[NDIM];
  FTYPE Unewlocal[NPR];
  FTYPE Uylupdate;
  FTYPE Urhoupdate;
  FTYPE graddotrhouynulimit;
  int numextrasreturn;




  ii=geom->i;
  jj=geom->j;
  kk=geom->k;
  pp=geom->p;

  coord_ijk(ii,jj,kk,pp,X) ;
  bl_coord_ijk(ii,jj,kk,pp,V) ;

  r=V[1];
  th=V[2];
  R = r*sin(th) ;
  rho = pr[RHO];
  u = pr[UU];
#if(DOYNU!=DONOYNU)
  ynu = pr[YNU];
#else
  ynu = 0.0;
#endif
#if(DOYL!=DONOYL)
  yl = pr[YL];
#else
  yl = 0.0;
#endif


  // approximately account for photon capture by black hole
  // 50% of masslesss particles from **stationary** isotropic emitters are captured by black hole
  // See Shapiro & Teukolsky (1983) and equation 2.81 and 2.82 in Shibata, Sekiguchi, Takahashi (2007)
  // Note that if MBH=0 then rph=0, so works for NS

  // For now, treat capture process as suppression of cooling/heating rates or any changes in fluid
  // In this way, capture by BH is treated as assuming fluid advected them in since assume if trapping photons then trapping fluid for sure.
  // In reality have to follow those photons/neutrinos to check if this is true, but good assumption.
  rph = 2.0*MBH*(1.0 + cos(2.0/3.0*(acos(-a/MBH))));
  photoncapture = (r>rph) ? 1.0 : 0.0;


  if(whichdatatype[primarytable]==4){
    // get neutrino and photon quantities
    get_extrasprocessed_kazfull(1, ii, jj, kk, pr, extras, processed); // (rho0,u)


    // MAXPROCESSEDEXTRAS entries
    qphoton=processed[QPHOTON];
    qm = processed[QNEUTRINO]; 
    graddotrhouyl = processed[GRADDOTRHOUYL]; 
    tthermaltot = processed[TTHERMAL]; 
    tdifftot = processed[TDIFF]; 
    rho_nu = processed[RHONU]; 
    p_nu = processed[PNU]; 
    s_nu = processed[SNU]; 
    ynulocal = processed[YNULOCAL];  // don't actually need this
    Ynuthermal = processed[YNUTHERMAL];
    // below are energies of *escaping* neutrinos
    enu = processed[ENUAVG]; 
    enue = processed[ENUE]; 
    enuebar = processed[ENUEBAR]; 

    // also store below additional as directly wanted quantities from extras
    lambdatot = extras[EXTRA19-EXTRA1]; // extras[0] is start
    lambdaintot = extras[EXTRA20-EXTRA1];
  }
  else if(whichdatatype[primarytable]==3){

    compute_allextras_kazfull(0, rho, u, &numextrasreturn, extras);

    qphoton = extras[EXTRA1-EXTRA1];
    qm = extras[EXTRA2-EXTRA1];
    graddotrhouyl = extras[EXTRA3-EXTRA1];
    tthermaltot = extras[EXTRA4-EXTRA1];
    tdifftot = extras[EXTRA5-EXTRA1];
    lambdatot = extras[EXTRA6-EXTRA1];
    lambdaintot = extras[EXTRA7-EXTRA1];
    enu = extras[EXTRA8-EXTRA1];
    enue = extras[EXTRA9-EXTRA1];
    enuebar = extras[EXTRA10-EXTRA1];
    Ynuthermal = extras[EXTRA11-EXTRA1];
  }


  // compute some other things
  //  dtau = dt/(q->ucon[TT]);

  /////////////
  //
  // Neutrino thermalization
  // Neutrino thermalization source term (due to neutrinos thermalizing fluid)
  //
  /////////////
  // regardless of the proper Lorentz transformation, this is the right dtau to use to limit ynu to Ynuthermal in the fluid frame
  // well, not really quite true due to advection terms
  // since U[RHO] will get updated, no point in trying to strictly limit using present U[RHO] due to advection terms
  // limit graddotrhouynu so that ynu can't go beyond Ynuthermal (haven't yet accounted for flux, so is not quite right constraint)
  // GODMARK: could add constraints later once flux known -- new function
  //  graddotrhouynu=rho*(Ynuthermal-ynu)/MAX(dtau,tthermaltot);
  //  if(rho*ucon[TT]*ynu+graddotrhouynu*dt>rho*ucon[TT]*Ynuthermal){
  //    graddotrhouynu = rho*(Ynuthermal - ynu)/dtau;
  //  }

  // equation is \nablda_\mu (\rho u^\mu Y_\nu) = \rho dY_\nu/d\tau = \rho (Y_{\nu,thermal} - Y_\nu)/\tau_{thermal} = -graddotrhouynu
  //  graddotrhouynu=-rho*(Ynuthermal-ynu)/MAX(dtau,tthermaltot);

  Urhoupdate = Ui[RHO] + (dUother[RHO])*dt; // assumes no physics source here for RHO (should add it if there is)
  // then based upon this criterion that doesn't necessarily imply exact dtau, strictly limit
  // Y_nu[new] = Unew[Ynu]/Unew[RHO] = (Ui[Ynu] + dUother[Ynu]*dt + dUsource[Ynu]*dt)/Unew[RHO] and solved for dUsource[Ynu]
  // this guarantees that post-advection an overshoot is truncated and leads to Ynu=Ynuthermal
  graddotrhouynulimit = (Ynuthermal*Urhoupdate - Ui[YNU])/dt - dUother[YNU];
  
  // dtlimit is chosen so the expression for graddotrhouynu is consistent
  // that is graddotrhoynu(original) = graddotrhouynulimit and solve for dtlimit
  dtlimit = -rho*(Ynuthermal-ynu)/graddotrhouynulimit;
  // GODMARK: should compare dtlimit with estimate of dtau -- for example can dtlimit<0?
  graddotrhouynu = -rho*(Ynuthermal-ynu)/MAX(dtlimit,tthermaltot);


  /////////////
  //
  // Lepton losses
  // Lepton source term (due to neutrinos escaping fluid)
  //
  // equation is:
  //
  // \nablda_\mu (\rho u^\mu Y_l) = \rho dY_l/d\tau = graddotrhouyl
  //
  // Equation is Unew = Ui + dUother*dt + gradU*dt
  //
  /////////////
  // predicted update
  Uylupdate = Ui[YL] + (dUother[YL] + graddotrhouyl)*dt;
  // limit update so Uylupdate<=0
  // limit Y_l update (Am I sure Y_l can't go negative? GODMARK)
  // here limit is perfect since any change in U[RHO] doesn't matter since U[YL]=0 -> Y_l=0 for all U[RHO]
  if(Uylupdate<0.0){
    graddotrhouyl = -Ui[YL]/dt - dUother[YL]; // so yl=0 after full Riemann + geometry + this source update
  }



  /////////////
  //
  // Lepton+Photon Cooling
  //
  /////////////
  // No obvious limit on gradUU since energy can be negative in ergosphere? and momentum and be + or - . GODMARK
  DLOOPA(j){
    gradUU[j]=-(qphoton+qm)*(q->ucov[j]);
  }





  /////////////////////////////
  //
  // now apply source terms
  //
  ////////////////////////////
  DLOOPA(j){
    // provide energy density per second cooled
    // source term has: \detg dU/d\tau u_\mu
    // dUcomp[][j+1] because j=0 implies UU as required (i.e. dUcomp in NPR form, not NDIM form)
    if(t > 0.){
      dUcomp[RADSOURCE][j+1]=photoncapture*gradUU[j]; // photons+neutrino energy lost in fluid frame.  Some may be captured by black hole, rest reaches infinity
    }
    else{
      dUcomp[RADSOURCE][j+1]=0.0;
    }
  }

  // GODMARK: need a way to track photons that reach infinity vs. black hole
  //*photoncapture


#if(DOYL!=DONOYL)
  j=YL; // accessing NPR-type quantity, so acces dUcomp[][j]
  if(t > 0.){
    dUcomp[RADSOURCE][j]=photoncapture*graddotrhouyl; // GODMARK: check sign of graddotrhouyl
  }
  else{
    dUcomp[RADSOURCE][j]=0;
  }
#endif

#if(DOYNU!=DONOYNU)
  j=YNU; // accessing NPR-type quantity, so acces dUcomp[][j]
  if(t > 0.){
    dUcomp[RADSOURCE][j]=photoncapture*graddotrhouynu; // GODMARK: check sign
  }
  else{
    dUcomp[RADSOURCE][j]=0;
  }
#endif
  


  return(0) ;




}












// compute things beyond simple EOS independent variables
void compute_EOS_parms_kazfull(FTYPE (*prim)[N2M][N3M][NPR])
{
  void compute_Hglobal(FTYPE (*prim)[N2M][N3M][NPR]);
  void compute_TDYNORYE_YNU_global(FTYPE (*prim)[N2M][N3M][NPR]);
  void compute_ups_global(FTYPE (*prim)[N2M][N3M][NPR]);
  

  //return; // see if HARM H is so different from stellar estimate of H that cuases bad temperature
  // didn't matter

  // compute H
  compute_Hglobal(prim);
  // compute Tdyn or Ye depending upon whichrnpmethod
  compute_TDYNORYE_YNU_global(prim);
  // compute neutriono u_\nu, p_\nu, and s_\nu
  compute_ups_global(prim);


}




void store_EOS_parms_kazfull(int i, int j, int k, int numparms, FTYPE *parlist)
{
  FTYPE TDYNORYEtouse;
  FTYPE YNUtouse;
  FTYPE Htouse[NUMHDIRECTIONS];
  FTYPE unu,pnu,snu;
  int hi;


  if(numparms!=NUMEOSGLOBALS){
    dualfprintf(fail_file,"numparms=%d NUMEOSGLOBALS=%d\n",numparms,NUMEOSGLOBALS);
    myexit(917692);
  }

  // must be right order when this function is used
  TDYNORYEtouse=parlist[TDYNORYEGLOBAL-FIRSTEOSGLOBAL];
  YNUtouse=parlist[YNUGLOBAL-FIRSTEOSGLOBAL];
  for(hi=0;hi<NUMHDIRECTIONS;hi++){
    Htouse[hi]=parlist[HGLOBAL-FIRSTEOSGLOBAL+hi]; // 5-3+hi and so for hi=0 gives 2 as required
  }
 

  if(TDYNORYEtouse>lineartablelimits[primarytable][TEOS][1]) TDYNORYEtouse = lineartablelimits[primarytable][TEOS][1];
  if(TDYNORYEtouse<lineartablelimits[primarytable][TEOS][0]) TDYNORYEtouse = lineartablelimits[primarytable][TEOS][0];

  if(YNUtouse>lineartablelimits[primarytable][YNUEOS][1]) YNUtouse = lineartablelimits[primarytable][YNUEOS][1];
  if(YNUtouse<lineartablelimits[primarytable][YNUEOS][0]) YNUtouse = lineartablelimits[primarytable][YNUEOS][0];

  if(whichdatatype[primarytable]!=4){
    for(hi=0;hi<NUMHDIRECTIONS;hi++){
      if(Htouse[hi]>lineartablelimits[primarytable][HEOS][1]) Htouse[hi] = lineartablelimits[primarytable][HEOS][1];
      if(Htouse[hi]<lineartablelimits[primarytable][HEOS][0]) Htouse[hi] = lineartablelimits[primarytable][HEOS][0];
    }
  }

  
  EOSglobal[TDYNORYEGLOBAL][i][j][k] = TDYNORYEtouse;
  EOSglobal[YNUGLOBAL][i][j][k] = YNUtouse;

  for(hi=0;hi<NUMHDIRECTIONS;hi++){
    EOSglobal[HGLOBAL+hi][i][j][k]=Htouse[hi];
  }

  EOSglobal[UNUGLOBAL][i][j][k]=parlist[UNUGLOBAL-FIRSTEOSGLOBAL];
  EOSglobal[PNUGLOBAL][i][j][k]=parlist[PNUGLOBAL-FIRSTEOSGLOBAL];
  EOSglobal[SNUGLOBAL][i][j][k]=parlist[SNUGLOBAL-FIRSTEOSGLOBAL];


}

void get_EOS_parms_kazfull(int i, int j, int k, int*numparms, FTYPE *parlist)
{
  int hi;

  // expected use of output should be know this ordering
  parlist[0]=EOSglobal[TDYNORYEGLOBAL][i][j][k];
  parlist[1]=EOSglobal[YNUGLOBAL][i][j][k];
  for(hi=0;hi<NUMHDIRECTIONS;hi++){
    parlist[HGLOBAL-FIRSTEOSGLOBAL+hi]=EOSglobal[HGLOBAL+hi][i][j][k];
  }
  parlist[UNUGLOBAL-FIRSTEOSGLOBAL]=EOSglobal[UNUGLOBAL][i][j][k];
  parlist[PNUGLOBAL-FIRSTEOSGLOBAL]=EOSglobal[PNUGLOBAL][i][j][k];
  parlist[SNUGLOBAL-FIRSTEOSGLOBAL]=EOSglobal[SNUGLOBAL][i][j][k];

  *numparms=NUMEOSGLOBALS;

}


