
#define PRIMECOORDS -1 // whatever the prime coordinate/metric is, used in transforms.c
#define CARTMINKMETRIC 0 // cartesian that is
#define BLCOORDS 1 // just set Rin oustide horizon for normal star
#define KSCOORDS 2 // as with BLCOORDS
#define HTMETRIC 3 // this defines exterior only, need specific Rin[\theta]
#define CYLMINKMETRIC 4 // cylindrical minkowski
#define HTMETRICACCURATE 5 // consistent expansion form of HT metric
#define SPCMINKMETRIC 6
#define UNIGRAVITY 7 // for RT test problem
#define KS_BH_TOV_COORDS 8 // KS BH+TOV mixed metric
#define KS_TOV_COORDS 9 // KS TOV metric
#define BL_TOV_COORDS 10 // BL TOV metric

// 0 : Boyer-Lindquist (based on r theta)
// 1 : Kerr-Schild (based on bl coords)
// contains metric definitions
// Kerr-Schild: future regularized, ep=-1, k=1

// check for metric.c on whether using spherical polar grid so can identify polar axis and r=0 singularities
#define ISSPCMCOORDNATIVE(whichcoord) (whichcoord==BLCOORDS || whichcoord==KSCOORDS || whichcoord==HTMETRIC || whichcoord==HTMETRICACCURATE || whichcoord==SPCMINKMETRIC || whichcoord==KS_BH_TOV_COORDS || whichcoord==KS_TOV_COORDS || whichcoord==BL_TOV_COORDS )

#define ISSPCMCOORD(whichcoord) (whichcoord==PRIMECOORDS && ISSPCMCOORDNATIVE(MCOORD) || ISSPCMCOORDNATIVE(whichcoord))


// black hole metrics
#define ISBLACKHOLEMCOORD(whichcoord) (whichcoord==BLCOORDS ||whichcoord==KSCOORDS ||whichcoord==KS_BH_TOV_COORDS ||whichcoord==KS_TOV_COORDS ||whichcoord==BL_TOV_COORDS)


