/////////////////////////
//
//
// HYBRID PPM SETTINGS (only used when lim=WENO5*)
//
//////////////////////////


// whether to use para with WENO
// CHANGINGMARK: Not correct with MPI -- e.g. produces monopoles for any FLUXB
#define PARAMODWENO 1   // default: 1

// whether to use para for only contacts (0) or full hybrid method based upon indicators (e.g. stiffindicator)
#define FULLHYBRID 1   // default: 1

// see interppoint.para.c on the meaning of these
// contact steepener is invalid procedure for general problems with gravity and such things (e.g. can steepend density profiles due to force balance! -- like disks!!!)
#define DOPPMCONTACTSTEEPMODWENO 0    // default: 0
// not sure if field method works well enough to use generally // not used? GODMARK
#define DOPPMSTEEPVARTYPEMODWENO 0    // default: 0
#define DOPPMREDUCEMODWENO 1    // default: 0 (stay speculative) -- well, causes problems.  In general strong shock will be too oscillatory so need to keep and hope that shock indicator only activates in strong shocks

// whether to use MC when WENO reduces (only replaces WENO part, not PARA part)
// Point is that WENO3 is much more diffusion than MC
// may introduce slight switching problem when MC goes sharp
#define USEMCFORLOWERORDERWENO 0   // default: 0










/////////////////////////
//
// PPM settings for para4gen() used by *all* para methods (lim=PARA,PARAFLAT,PARALINE)
//
//////////////////////////

#define JONPARAREDUCE 0  // default: 0

// actually helps keep para flat, but also stepeends in perfect limit of infinite drop in density
// Helps #2 (double rarefaction) reach low u/rho! (not significantly better)
#define JONPARASTEEP 0  // default: 0

// whether to assume input y's are cell average or point values
#define AVGINPUT 0  // default: 0

// don't have to use PARA to get original l,r
// 0 = para   1 = PARA2LIM limiter directly (then acts like MCSTEEP)
#define WHICHLIMITERTOUSEFORLR 0  // default: 0

// Noticed with MCSTEEP that allowing extremum causes catastrophic problems with torus problem with field
 // assume monocheck on, or else results in crazy behavior
// Even with PARA this (monotonic condition check) causes problems (e.g. horizon failures and polar blow-ups)
// With PARAFLAT or PARALINE, enough ddq's (second derivatives) to ensure really a parabola at extrema
#define PARAGENDQALLOWEXTREMUM 1  // default: 1



/////////////////////////
//
// PPM SETTINGS for checkparamonotonicity()
//
//////////////////////////
 
// whether to reduce to DONOR or MONO method when parabola is not monotonic
// DONOR: 0
// PARA2LIM: 1 (not robust)
#define NONMONOLIM 0  // default: 0




///////////////////////////////////////
//
// FOR PARAPL (i.e. lim=PARAFLAT) only
//
///////////////////////////////////////

// DOPPMREDUCE probably needed if using contact steepener and other such things, but does result in more diffusion result that makes benefit of general PPM not good
#define DOPPMREDUCE 1 // default: 1
// whether to use PPM contact steepener (often results in bad behvior in general even if works very well for 1D LW2003 tests)
#define DOPPMCONTACTSTEEP 0 // default: 0


/////////////////////////
//
// PPM SETTINGS for ftilde() (lim=PARAFLAT and PARALINE only)
//
//////////////////////////

//#define SP0 0.65 // for strong shocks
#define SP0 0.75 // FLASH

#define SPSMOOTH 0.5  // default: 0.5


/////////////////////////
//
//
// PPM SETTINGS for parasteep() (only for lim=PARAFLAT and PARALINE)
//
//////////////////////////

 // 0 == steepend only density   1 == steepend density and cross-fields (i.e. current sheets)
#define DOPPMSTEEPVARTYPE 0   // default: 0





