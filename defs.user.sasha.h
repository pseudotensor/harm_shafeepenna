//GODMARK: don't know how to use FTYPE here -- FTYPE defined in global.h 
// need not static so other files can access settings
struct Ccoordparams coordparams;

FTYPE interpn( int order, FTYPE x_eval,  FTYPE x1, FTYPE f1, FTYPE x2, FTYPE f2, FTYPE x3, FTYPE f3, FTYPE x4, FTYPE f4, FTYPE x5, FTYPE f5, FTYPE x6, FTYPE f6 );
int bounds_generate( int i, int j, int k, FTYPE *prim );
