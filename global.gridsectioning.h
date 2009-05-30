// Sasha related enerregion's to grid sectioning, e.g. NUMENERREGIONS



#if( DOGRIDSECTIONING )

//corrections to be applied to loops
#define SHIFTX1DN (enerposreg[ACTIVEREGION][X1DN]-0)
#define SHIFTX1UP (enerposreg[ACTIVEREGION][X1UP]-(N1-1))
#define SHIFTX2DN (enerposreg[ACTIVEREGION][X2DN]-0)
#define SHIFTX2UP (enerposreg[ACTIVEREGION][X2UP]-(N2-1))
#define SHIFTX3DN (enerposreg[ACTIVEREGION][X3DN]-0)
#define SHIFTX3UP (enerposreg[ACTIVEREGION][X3UP]-(N3-1))

#else

//no sectioning -- zero corrections
#define SHIFTX1DN (0)
#define SHIFTX1UP (0)
#define SHIFTX2DN (0)
#define SHIFTX2UP (0)
#define SHIFTX3DN (0)
#define SHIFTX3UP (0)

#endif

#if( DOGRIDSECTIONING )
#define WITHINACTIVESECTION(ri,rj,rk) ( ri >= enerposreg[ACTIVEREGION][X1DN] && ri <= enerposreg[ACTIVEREGION][X1UP] \
                                     && rj >= enerposreg[ACTIVEREGION][X2DN] && rj <= enerposreg[ACTIVEREGION][X2UP] \
                                     && rk >= enerposreg[ACTIVEREGION][X3DN] && rk <= enerposreg[ACTIVEREGION][X3UP] )
#else
#define WITHINACTIVESECTION(ri,rj,rk) (1)  //always within active section since no sectioning
#endif



