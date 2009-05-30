// full loop + 1 (shift away from boundary) on inner edge for comptuing emf or corner quantities
//#define LOOPINFP13 for(k=INFULLP13;k<=OUTFULL3;k++)  //SASMARK SECTIONMARK these three loops not used anywhere anymore?
//#define LOOPINFP12 for(j=INFULLP12;j<=OUTFULL2;j++)
//#define LOOPINFP11 for(i=INFULLP11;i<=OUTFULL1;i++)

#define COMPLOOPINFP13 for(k=INFULLP13+SHIFTX3DN;k<=OUTFULL3+SHIFTX3UP;k++)
#define COMPLOOPINFP12 for(j=INFULLP12+SHIFTX2DN;j<=OUTFULL2+SHIFTX2UP;j++)
#define COMPLOOPINFP11 for(i=INFULLP11+SHIFTX1DN;i<=OUTFULL1+SHIFTX1UP;i++)



//#define LOOPOUTFM13 for(k=INFULL3;k<=OUTFULLM13;k++)
//#define LOOPOUTFM12 for(j=INFULL2;j<=OUTFULLM12;j++)
//#define LOOPOUTFM11 for(i=INFULL1;i<=OUTFULLM11;i++)

#define COMPLOOPOUTFM13 for(k=INFULL3+SHIFTX3DN;k<=OUTFULLM13+SHIFTX3UP;k++)
#define COMPLOOPOUTFM12 for(j=INFULL2+SHIFTX2DN;j<=OUTFULLM12+SHIFTX2UP;j++)
#define COMPLOOPOUTFM11 for(i=INFULL1+SHIFTX1DN;i<=OUTFULLM11+SHIFTX1UP;i++)


// SECTIONMARK: needs to be limited
#define COMPLOOPF COMPLOOPF3 COMPLOOPF2 COMPLOOPF1
#define COMPLOOPH CPMPLOOPH3 COMPLOOPH2 COMPLOOPH1
#define COMPLOOPP1 COMPLOOPP13 COMPLOOPP12 COMPLOOPP11
#define COMPLOOPFMH CMPLOOPFMH3 CMPLOOPFMH2 CMPLOOPFMH1   //SASMARK where is COMPLOOPFMH used?
#define COMPLOOPHMFP CMPLOOPHMFP3 CMPLOOPHMFP2 CMPLOOPHMFP1
#define COMPLOOPHP CMPLOOPHP3 CMPLOOPHP2 CMPLOOPHP1

// these loops used for computational purposes
#define COMPLOOPF3 for(k=INFULL3+SHIFTX3DN;k<=OUTFULL3+SHIFTX3UP;k++)
#define COMPLOOPF2 for(j=INFULL2+SHIFTX2DN;j<=OUTFULL2+SHIFTX2UP;j++)
#define COMPLOOPF1 for(i=INFULL1+SHIFTX1DN;i<=OUTFULL1+SHIFTX1UP;i++)

#define NDIMENNOT1(dimen) (NDIMEN(dimen)!=1)
#define COMPLOOPFDIMEN(dimen,ind)  for( ind =  0                 - NDIMENBND(dimen) + SHIFTDIR(DIR(dimen,-1)); \
                                        ind <= NDIMEN(dimen) - 1 + NDIMENBND(dimen) + SHIFTDIR(DIR(dimen,+1)); \
                                        ind++ )


#define COMPLOOPH3 for(k=INHALF3+SHIFTX3DN;k<=OUTHALF3+SHIFTX3UP;k++)
#define COMPLOOPH2 for(j=INHALF2+SHIFTX2DN;j<=OUTHALF2+SHIFTX2UP;j++)
#define COMPLOOPH1 for(i=INHALF1+SHIFTX1DN;i<=OUTHALF1+SHIFTX1UP;i++)

#define COMPLOOPP13 for(k=INP13+SHIFTX3DN;k<=OUTP13+SHIFTX3UP;k++)
#define COMPLOOPP12 for(j=INP12+SHIFTX2DN;j<=OUTP12+SHIFTX2UP;j++)
#define COMPLOOPP11 for(i=INP11+SHIFTX1DN;i<=OUTP11+SHIFTX1UP;i++)

#define COMPLOOPFMHP3 for(k=INFULL3+SHIFTX3DN;k<=OUTHALF3+SHIFTX3UP;k++)
#define COMPLOOPFMHP2 for(j=INFULL2+SHIFTX2DN;j<=OUTHALF2+SHIFTX2UP;j++)
#define COMPLOOPFMHP1 for(i=INFULL1+SHIFTX1DN;i<=OUTHALF1+SHIFTX1UP;i++)

#define COMPLOOPHP3 for(k=0+SHIFTX3DN;k<=OUTHALF3+SHIFTX3UP;k++)
#define COMPLOOPHP2 for(j=0+SHIFTX2DN;j<=OUTHALF2+SHIFTX2UP;j++)
#define COMPLOOPHP1 for(i=0+SHIFTX1DN;i<=OUTHALF1+SHIFTX1UP;i++)



// SECTIONMARK: needs to be limited
#define COMPLOOPN3 for(k=0+SHIFTX3DN;k<=N3-1+SHIFTX3UP;k++)
#define COMPLOOPN2 for(j=0+SHIFTX2DN;j<=N2-1+SHIFTX2UP;j++)
#define COMPLOOPN1 for(i=0+SHIFTX1DN;i<=N1-1+SHIFTX1UP;i++)

#define COMPLOOPN COMPLOOPN3 COMPLOOPN2 COMPLOOPN1

// SECTIONMARK: should be limited
#define COMPFULLLOOP COMPLOOPF1 COMPLOOPF2 COMPLOOPF3

#define COMPLOOPF3 for(k=INFULL3+SHIFTX3DN;k<=OUTFULL3+SHIFTX3UP;k++)
#define COMPLOOPF2 for(j=INFULL2+SHIFTX2DN;j<=OUTFULL2+SHIFTX2UP;j++)
#define COMPLOOPF1 for(i=INFULL1+SHIFTX1DN;i<=OUTFULL1+SHIFTX1UP;i++)


// SECTIONMARK: leave LOOPF? alone
//#define FULLLOOP LOOPF1 LOOPF2 LOOPF3

// computing emf for FLUXCT
// SECTIONMARK: needs to be limited
// change names too so looks "computational"
#define COMPLOOPINFP1 COMPLOOPINFP11 COMPLOOPINFP12 COMPLOOPINFP13

#define COMPLOOPINFP1dir1full COMPLOOPF1 COMPLOOPINFP12 COMPLOOPINFP13

#define COMPLOOPINFP1dir2full COMPLOOPINFP11 COMPLOOPF2 COMPLOOPINFP13

#define COMPLOOPINFP1dir3full COMPLOOPINFP11 COMPLOOPINFP12 COMPLOOPF3

#define COMPLOOPINFP1dir23full COMPLOOPINFP11 COMPLOOPF2 COMPLOOPF3

#define COMPLOOPINFP1dir13full COMPLOOPF1 COMPLOOPINFP12 COMPLOOPF3

#define COMPLOOPINFP1dir12full COMPLOOPF1 COMPLOOPF2 COMPLOOPINFP13


// computing emf for FLUXCD
#define COMPLOOPOUTFM1 COMPLOOPOUTFM11 COMPLOOPOUTFM12 COMPLOOPOUTFM13

#define COMPLOOPOUTFM1dir1full COMPLOOPF1 COMPLOOPOUTFM12 COMPLOOPOUTFM13

#define COMPLOOPOUTFM1dir2full COMPLOOPOUTFM11 COMPLOOPF2 COMPLOOPOUTFM13

#define COMPLOOPOUTFM1dir3full COMPLOOPOUTFM11 COMPLOOPOUTFM12 COMPLOOPF3


// larger loop than full for cornered quantities such as emf defined on corners that need to be initialized for boundary condition reasons
//#define FULLLOOPP1 LOOPFP11 LOOPFP12 LOOPFP13
#define COMPFULLLOOPP1 COMPLOOPFP11 COMPLOOPFP12 COMPLOOPFP13

// full loop + 1 on outer edge for emf or corner quantities
#define COMPLOOPFP13 for(k=INFULL3+SHIFTX3DN;k<=OUTFULLP13+SHIFTX3UP;k++)
#define COMPLOOPFP12 for(j=INFULL2+SHIFTX2DN;j<=OUTFULLP12+SHIFTX2UP;j++)
#define COMPLOOPFP11 for(i=INFULL1+SHIFTX1DN;i<=OUTFULLP11+SHIFTX1UP;i++)


// divb loop
//#define LOOPDIVB LOOPP11 LOOPP12 LOOPP13
// boundary zones may not require divb=0 since proxy for flux
//#define LOOPDIVB LOOPC1 LOOPC2 LOOPC3
#define COMPLOOPDIVB COMPLOOPN1 COMPLOOPN2 COMPLOOPN3
//ZSLOOP(-N1BND+1,N1+1,-1,N2+1)






//SASMARK SECTIONMARK PUT MACRODEFS HERE


////////////////////////////////////////////
//
// normal non-MPI or standard MPI  (was after #include "mympi.h")
//
////////////////////////////////////////////

#if( SIMULBCCALC <= 0 )  //SASMARK not sure if need to require = -1 or <= 0 

// SECTIONMARK: needs to be limited
#define COMPZLOOP ZSLOOP(0+SHIFTX1DN,N1-1+SHIFTX1UP,0+SHIFTX2DN,N2-1+SHIFTX2UP,0+SHIFTX3DN,N3-1+SHIFTX3UP)

// can use COMPZSLOOP if don't care about extra calculations 
//SASMARK SECTIONMARK: modified is/ie, etc. globally assuming everything depends on the computational box linearly
#define COMPZSLOOP(istart,istop,jstart,jstop,kstart,kstop) ZSLOOP(istart+SHIFTX1DN,istop+SHIFTX1UP,jstart+SHIFTX2DN,jstop+SHIFTX2UP,kstart+SHIFTX3DN,kstop+SHIFTX3UP)

// these are special loops with very careful calculation ranges to avoid extra calculations
// but assumes are global variables which are assigned, so remains intact under whatever circumstances needed next
#define COMPFZLOOP(istart,jstart,kstart) ZSLOOP(istart,OUTM1,jstart,OUTM2,kstart,OUTM3)  //not used anywhere for SIMULBCCALC == -1
#define COMPEMFZLOOP ZSLOOP(0,OUTM1,0,OUTM2,0,OUTM3)             //not used anywhere at all (only in "//" comments)
#define COMPPREEMFZLOOP ZSLOOP(INM1,OUTM1,INM2,OUTM2,INM3,OUTM3) //not used anywhere at all (only in "//" comments)
#define COMPF1CTZLOOP ZSLOOP(0,OUTM1,0,N2-1,0,N3-1)              //not used anywhere at all (only in "//" comments)
#define COMPF2CTZLOOP ZSLOOP(0,N1-1,0,OUTM2,0,N3-1)              //not used anywhere at all (only in "//" comments)
#define COMPF3CTZLOOP ZSLOOP(0,N1-1,0,N2-1,0,OUTM3)              //not used anywhere at all (only in "//" comments)
#define COMPPREDQZLOOP FULLLOOP                                  //not used anywhere at all (only in "//" comments)
#define COMPDQZLOOP COMPZSLOOP(INM1,OUTM1,INM2,OUTM2,INM3,OUTM3)

////////////////////////////////////////////
//
// normal non-MPI or standard MPI
//
////////////////////////////////////////////



//#define CZLOOP ZLOOP
//// can use CZSLOOP if don't care about extra calculations
//#define CZSLOOP(istart,istop,jstart,jstop,kstart,kstop) ZSLOOP(istart,istop,jstart,jstop,kstart,kstop)
//// these are special loops with very careful calculation ranges to avoid extra calculations
//// but assumes are global variables which are assigned, so remains intact under whatever circumstances needed next
//#define FZLOOP(istart,jstart,kstart) ZSLOOP(istart,OUTM1,jstart,OUTM2,kstart,OUTM3)
//#define EMFZLOOP ZSLOOP(0,OUTM1,0,OUTM2,0,OUTM3)
//#define PREEMFZLOOP ZSLOOP(INM1,OUTM1,INM2,OUTM2,INM3,OUTM3)
//#define COMPF1CTZLOOP ZSLOOP(0,OUTM1,0,N2-1,0,N3-1)
//#define COMPF2CTZLOOP ZSLOOP(0,N1-1,0,OUTM2,0,N3-1)
//#define COMPF3CTZLOOP ZSLOOP(0,N1-1,0,N2-1,0,OUTM3)
//#define COMPPREDQZLOOP FULLLOOP
//#define COMPDQZLOOP ZSLOOP(INM1,OUTM1,INM2,OUTM2,INM3,OUTM3)


#endif
