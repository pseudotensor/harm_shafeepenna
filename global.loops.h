

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// LOOPS
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







/////////////////////////////////////
/////////////////////////////////////
//
// MOST GENERAL LOOPS
//
/////////////////////////////////////
/////////////////////////////////////


// Sasha's dad and google say memcopy() and memmove() can be much faster than loop if compiler doesn't recognize this optimization.  So best to force it in simple cases
#define USE_MEMCPY 1

// applies to continugous memory regions (i.e. feed single pointer for source, size of data to copy in its own dimensions)
#if(USE_MEMCPY)
// iter not used
#define GENFORALL(iter,src,dest,numelem) memcpy(dest, src, numelem*sizeof(src[0]))
#define GENFORALLOVERLAP(iter,src,dest,numelem) memmove(dest, src, numelem*sizeof(src[0]))
#else
// assumes iteration is single (i.e. all data from start to finish)
#define GENFORALL(iter,src,dest,numelem) \
for (iter = 0; iter < numelem; iter++) \
{\
dest[iter] = src[iter];\
}
#define GENFORALLOVERLAP(iter,src,dest,numelem) GENFORALL(iter,src,dest,numelem)
#endif


// these loops used for general purposes
#define LOOPF3 for(k=INFULL3;k<=OUTFULL3;k++)
#define LOOPF2 for(j=INFULL2;j<=OUTFULL2;j++)
#define LOOPF1 for(i=INFULL1;i<=OUTFULL1;i++)


// full loop + 1 on outer edge for emf or corner quantities
#define LOOPFP13 for(k=INFULL3;k<=OUTFULLP13;k++)
#define LOOPFP12 for(j=INFULL2;j<=OUTFULLP12;j++)
#define LOOPFP11 for(i=INFULL1;i<=OUTFULLP11;i++)


// full loop + 1 (shift away from boundary) on inner edge for comptuing emf or corner quantities
//#define LOOPINFP13 for(k=INFULLP13;k<=OUTFULL3;k++)
//#define LOOPINFP12 for(j=INFULLP12;j<=OUTFULL2;j++)
//#define LOOPINFP11 for(i=INFULLP11;i<=OUTFULL1;i++)


//#define LOOPOUTFM13 for(k=INFULL3;k<=OUTFULLM13;k++)
//#define LOOPOUTFM12 for(j=INFULL2;j<=OUTFULLM12;j++)
//#define LOOPOUTFM11 for(i=INFULL1;i<=OUTFULLM11;i++)

#define LOOPH3 for(k=INHALF3;k<=OUTHALF3;k++)
#define LOOPH2 for(j=INHALF2;j<=OUTHALF2;j++)
#define LOOPH1 for(i=INHALF1;i<=OUTHALF1;i++)

#define LOOPP13 for(k=INP13;k<=OUTP13;k++)
#define LOOPP12 for(j=INP12;j<=OUTP12;j++)
#define LOOPP11 for(i=INP11;i<=OUTP11;i++)

#define LOOPN3 for(k=0;k<=N3-1;k++)
#define LOOPN2 for(j=0;j<=N2-1;j++)
#define LOOPN1 for(i=0;i<=N1-1;i++)

#define LOOPFMHP3 for(k=INFULL3;k<=OUTHALF3;k++)
#define LOOPFMHP2 for(j=INFULL2;j<=OUTHALF2;j++)
#define LOOPFMHP1 for(i=INFULL1;i<=OUTHALF1;i++)

#define LOOPHMFP3 for(k=INHALF3;k<=OUTFULL3;k++)
#define LOOPHMFP2 for(j=INHALF2;j<=OUTFULL2;j++)
#define LOOPHMFP1 for(i=INHALF1;i<=OUTFULL1;i++)

#define LOOPHP3 for(k=0;k<=OUTHALF3;k++)
#define LOOPHP2 for(j=0;j<=OUTHALF2;j++)
#define LOOPHP1 for(i=0;i<=OUTHALF1;i++)

// below used for initialization and such, not a computational issue
#define LOOPF LOOPF3 LOOPF2 LOOPF1
#define LOOPH LOOPH3 LOOPH2 LOOPH1
#define LOOPP1 LOOPP13 LOOPP12 LOOPP11
#define LOOP LOOPN3 LOOPN2 LOOPN1
#define LOOPFMHP LOOPFMHP3 LOOPFMHP2 LOOPFMHP1
#define LOOPHMFP LOOPHMFP3 LOOPHMFP2 LOOPHMFP1
#define LOOPHP LOOPHP3 LOOPHP2 LOOPHP1

#define LOOPINT3 for(k=intix3;k<intox3;k++)
#define LOOPINT2 for(j=intix2;j<intox2;j++)
#define LOOPINT1 for(i=intix1;i<intox1;i++)



#define LOOPFC LOOPF
#define LOOPHC LOOPH
#define LOOPFMHPC LOOPFMHP
#define LOOPHMFPC LOOPHMFP
#define LOOPHPC LOOPHP


#define LOOPC3 LOOPN3
#define LOOPC2 LOOPN2
#define LOOPC1 LOOPN1

#define LOOPC LOOPC3 LOOPC2 LOOPC1

// general loop, but assumes i,j,k used
#define ZSLOOP(istart,istop,jstart,jstop,kstart,kstop) \
	for(i=istart;i<=istop;i++)\
	for(j=jstart;j<=jstop;j++)\
	for(k=kstart;k<=kstop;k++)

// general loop for any indicies
#define GENLOOP(i,j,k,istart,istop,jstart,jstop,kstart,kstop) \
        for((i)=(istart);(i)<=(istop);(i)++)\
        for((j)=(jstart);(j)<=(jstop);(j)++)\
        for((k)=(kstart);(k)<=(kstop);(k)++)

// general loop for any indicies
#define SUPERGENLOOP(i,j,k,istart,istop,jstart,jstop,kstart,kstop,di,dj,dk) \
        for((i)=(istart);(i)<=(istop);(i)+=(di))\
        for((j)=(jstart);(j)<=(jstop);(j)+=(dj))\
        for((k)=(kstart);(k)<=(kstop);(k)+=(dk))

/* loop over all active zones */
//#define ZLOOP for(i=0;i<=N1-1;i++)for(j=0;j<=N2-1;j++)for(k=0;k<=N3-1;k++)
#define ZLOOP ZSLOOP(0,N1-1,0,N2-1,0,N3-1)


//#define FULLLOOP ZSLOOP(-N1BND, N1 -1 + N1BND, -N2BND, N2 -1 + N2BND, -N3BND, N3 -1 + N3BND)
#define FULLLOOP LOOPF1 LOOPF2 LOOPF3

//#if(BOUNDPLPR&&NOFLUXCTONX1DN)
//#define FULLLOOPflux1 LOOPF1 LOOPF2 LOOPF3 if(!(startpos[1]+i==0))
//#else
//#define FULLLOOPflux1 LOOPF1 LOOPF2 LOOPF3
//#endif




// computing emf for FLUXCT
		    //#define LOOPINFP1 LOOPINFP11 LOOPINFP12 LOOPINFP13

//#define LOOPINFP1dir1full LOOPF1 LOOPINFP12 LOOPINFP13

//#define LOOPINFP1dir2full LOOPINFP11 LOOPF2 LOOPINFP13

//#define LOOPINFP1dir3full LOOPINFP11 LOOPINFP12 LOOPF3

//#define LOOPINFP1dir23full LOOPINFP11 LOOPF2 LOOPF3

//#define LOOPINFP1dir13full LOOPF1 LOOPINFP12 LOOPF3

//#define LOOPINFP1dir12full LOOPF1 LOOPF2 LOOPINFP13


// computing emf for FLUXCD
		    //#define LOOPOUTFM1 LOOPOUTFM11 LOOPOUTFM12 LOOPOUTFM13

//#if(BOUNDPLPR&&NOFLUXCTONX1DN)
//#define LOOPOUTFM1dir1fullflux LOOPF1 LOOPOUTFM12 LOOPOUTFM13 if(!(startpos[1]+i==0))
//#else
//#define LOOPOUTFM1dir1fullflux LOOPF1 LOOPOUTFM12 LOOPOUTFM13
//#endif

		    //#define LOOPOUTFM1dir1full LOOPF1 LOOPOUTFM12 LOOPOUTFM13

//#define LOOPOUTFM1dir2full LOOPOUTFM11 LOOPF2 LOOPOUTFM13

//#define LOOPOUTFM1dir3full LOOPOUTFM11 LOOPOUTFM12 LOOPF3





// larger loop than full for cornered quantities such as emf defined on corners that need to be initialized for boundary condition reasons
#define FULLLOOPP1 LOOPFP11 LOOPFP12 LOOPFP13

//#define WSPEEDLOOP ZSLOOP(-SHIFT1,N1-1+2*SHIFT1,-SHIFT2,N2-1+2*SHIFT2,-SHIFT3,N3-1+2*SHIFT3)


//#define PLUSLOOP ZSLOOP(-1, N1, -1, N2, -1, N3)

// below same as FULLLOOP if NBND=2
//#define PLUSPLUSLOOP ZSLOOP(-2, N1+1, -2, N2+1, -2, N3+1)




// divb loop (for diagnostics only)
//#define LOOPDIVB LOOPP11 LOOPP12 LOOPP13
// boundary zones may not require divb=0 since proxy for flux
#define LOOPDIVB LOOPC1 LOOPC2 LOOPC3
//ZSLOOP(-N1BND+1,N1+1,-1,N2+1)








/////////////////////////////////////
/////////////////////////////////////
//
// BOUNDARY CONDITION RELATED LOOPS
//
/////////////////////////////////////
/////////////////////////////////////

// general loop for boundary conditions
#define BOUNDLOOPF(i,in,out) for(i=in;i<=out;i++)

// for X1-dir
#define LOOPX1dir BOUNDLOOPF(j,innormalloop[2],outnormalloop[2]) BOUNDLOOPF(k,innormalloop[3],outnormalloop[3])

// for X2-dir
#define LOOPX2dir BOUNDLOOPF(i,inboundloop[1],outboundloop[1]) BOUNDLOOPF(k,innormalloop[3],outnormalloop[3])

// for X3-dir
#define LOOPX3dir BOUNDLOOPF(i,inboundloop[1],outboundloop[1]) BOUNDLOOPF(j,inboundloop[2],outboundloop[2])


// loop over all directions (block of data in i,j,k only going out as far as the boundary conditions would go out)
#define LOOPXalldir BOUNDLOOPF(i,inboundloop[1],outboundloop[1]) BOUNDLOOPF(j,inboundloop[2],outboundloop[2])  BOUNDLOOPF(k,inboundloop[3],outboundloop[3])

/////////////////////////////////////
//
// boundary loops for CENT quantities
//
/////////////////////////////////////
//
// normally start is a lower index than end and then loops take into account how to deal with this
//
// most general loop ignorant of start begin before end
#define LOOPBOUNDGENMORE(i,iloopstart,iloopend,iloopstep) for(i=iloopstart;i<=iloopend;i+=iloopstep)
// make start near boundary so periodic works for MAXBND up to 2N
#define LOOPBOUNDINGEN(i,start,end) for(i=end;i>=start;i--)
// bound entire region inside non-evolved portion of grid
// special horizon code still needed for multiple CPUs -- GODMARK -- unless chose range of CPUs inside horizon to apply inner BCs
#define LOOPBOUNDINMOREGEN(i,start,end,ri) for(i=ri-1;i>=start;i--)
#define LOOPBOUNDOUTGEN(i,start,end) for(i=start;i<=end;i++)

// user should use set_boundloop() to control shifts
#define LOOPBOUND1IN  LOOPBOUNDINMOREGEN(i,inoutlohi[POINTDOWN][POINTDOWN][1],inoutlohi[POINTDOWN][POINTUP][1],ri)
#define LOOPBOUND1OUT LOOPBOUNDOUTGEN   (i,inoutlohi[POINTUP][POINTDOWN][1],inoutlohi[POINTUP][POINTUP][1])

#define LOOPBOUND2IN  LOOPBOUNDINGEN (j,inoutlohi[POINTDOWN][POINTDOWN][2],inoutlohi[POINTDOWN][POINTUP][2])
#define LOOPBOUND2OUT LOOPBOUNDOUTGEN(j,inoutlohi[POINTUP][POINTDOWN][2],inoutlohi[POINTUP][POINTUP][2])

#define LOOPBOUND3IN  LOOPBOUNDINGEN (k,inoutlohi[POINTDOWN][POINTDOWN][3],inoutlohi[POINTDOWN][POINTUP][3])
#define LOOPBOUND3OUT LOOPBOUNDOUTGEN(k,inoutlohi[POINTUP][POINTDOWN][3],inoutlohi[POINTUP][POINTUP][3])

#define LOWERBOUND1 (inoutlohi[POINTDOWN][POINTDOWN][1])
#define UPPERBOUND1 (inoutlohi[POINTUP][POINTUP][1])










/////////////////////////////////////
/////////////////////////////////////
//
// DIAGNOSTIC RELATED LOOPS
//
/////////////////////////////////////
/////////////////////////////////////



/* want dump output to be ordered in radius first!! */
#define DUMPLOOP(istart,istop,jstart,jstop,kstart,kstop) \
	for(k=kstart;k<=kstop;k++)\
	for(j=jstart;j<=jstop;j++)\
	for(i=istart;i<=istop;i++)

#if(FULLOUTPUT==0)
#define EXTRADUMP1 0
#define EXTRADUMP2 0
#define EXTRADUMP3 0
#else
#define EXTRADUMP1T FULLOUTPUT*N1NOT1
#define EXTRADUMP2T FULLOUTPUT*N2NOT1
#define EXTRADUMP3T FULLOUTPUT*N3NOT1

#define EXTRADUMP1 ((EXTRADUMP1T>N1BND) ? N1BND : EXTRADUMP1T)
#define EXTRADUMP2 ((EXTRADUMP2T>N2BND) ? N2BND : EXTRADUMP2T)
#define EXTRADUMP3 ((EXTRADUMP3T>N3BND) ? N3BND : EXTRADUMP3T)

#endif

#if(FULLOUTPUT==0)
#define DUMPGENLOOP DUMPLOOP(0,N1-1,0,N2-1,0,N3-1)
#else
#define DUMPGENLOOP DUMPLOOP(-EXTRADUMP1,N1-1+EXTRADUMP1,-EXTRADUMP2,N2-1+EXTRADUMP2,-EXTRADUMP3,N3-1+EXTRADUMP3)
#endif


// defines whether within the enerregion
// considered loop-related macro
#define WITHINENERREGION(theenerpos,i,j,k) (i>=theenerpos[X1DN])&&(i<=theenerpos[X1UP])&&(j>=theenerpos[X2DN])&&(j<=theenerpos[X2UP])&&(k>=theenerpos[X3DN])&&(k<=theenerpos[X3UP]) 





  //#define IMAGELOOP(istart,istop,jstart,jstop,kstart,kstop) \
  //	for(k=kstart;k<=kstop;k++)\
  //	for(j=jstart;j<=jstop;j++)\
  //	for(i=istart;i<=istop;i++)

//#define OLDIMAGELOOP for(j=N2-1;j>=0;j--) for(i=0;i<N1;i++)	// nasty 
								// to
								// deal 
								// with



/////////////////////////////////////
/////////////////////////////////////
//
// full CPU 1-D loops
//
/////////////////////////////////////
/////////////////////////////////////

#define NUMGRAVPOS (ncpux1*N1+N1BND*2+1)
#define GRAVLOOP(ii) for(ii=-N1BND;ii<ncpux1*N1+N1BND+1;ii++)
#define GRAVLOOPACTIVE(ii) for(ii=0;ii<ncpux1*N1;ii++)
#define DUMPGRAVLOOP(ii) GRAVLOOPACTIVE(ii) // normal
//#define DUMPGRAVLOOP(ii) GRAVLOOP(ii) // abnormal, debugging stuff in SM






/////////////////////////////////////
/////////////////////////////////////
//
// PER-POINT LOOPS
//
/////////////////////////////////////
/////////////////////////////////////


  // check for existence in bad form using:
// grep "PLOOP" *.c | grep --invert-match "PLOOP("

// PLOOP controls looping over conserved or primitive -type quantities except during interpolation or other listed below
#if(1) // for now always control interpolated quantities (used to be only for SLPITNPR)

#define PLOOP(pl) for(plglobal=nprstart,pl=nprlist[plglobal];plglobal<=nprend;plglobal++,pl=nprlist[plglobal])
#define PLOOPINTERP(pl) for(pl2global=npr2interpstart,pl=npr2interplist[pl2global];pl2global<=npr2interpend;pl2global++,pl=npr2interplist[pl2global])
#define PLOOPNOTINTERP(pl) for(pl3global=npr2notinterpstart,pl=npr2notinterplist[pl3global];pl3global<=npr2notinterpend;pl3global++,pl=npr2notinterplist[pl3global])

// loop over all bounding Primitive variables
#define PBOUNDLOOP(pl) for(pl4global=nprboundstart,pl=nprboundlist[pl4global];pl4global<=nprboundend;pl4global++,pl=nprboundlist[pl4global])

// loop over all bounding flux variables
#define PFLUXBOUNDLOOP(pl) for(pl5global=nprfluxboundstart,pl=nprfluxboundlist[pl5global];pl5global<=nprfluxboundend;pl5global++,pl=nprfluxboundlist[pl5global])

// loop over all dumped Primitive variables
#define PDUMPLOOP(pd) for(pl6global=nprdumpstart,pd=nprdumplist[pl6global];pl6global<=nprdumpend;pl6global++,pd=nprdumplist[pl6global])

// loop over all inversion Primitive variables -- only over 5 quantities (not field)
#define PINVERTLOOP(pi) for(pl7global=nprinvertstart,pl=nprinvertlist[pl7global];pl7global<=nprinvertend;pl7global++,pl=nprinvertlist[pl7global])


// to be used locally:
#define NUMPRIMLOOP(pl) for(pllocal=nprlocalstart,pl=nprlocallist[pllocal];pllocal<=nprlocalend;pllocal++,pl=nprlocallist[pllocal])

#define NUMPRIMLOOPGEN(pllocal,pl) for(pllocal=nprlocalstart,pl=nprlocallist[pllocal];pllocal<=nprlocalend;pllocal++,pl=nprlocallist[pllocal])




#else


/* loop over all Primitive variables */
#define PLOOP(pl) for(pl=0;pl<NPR;pl++) // original
/* loop over all center to edge variables */
#define PLOOPINTERP(pl) for(pl=0;pl<NPR2INTERP;pl++)
#define PLOOPNOTINTERP(pl) for(pl=0;pl<0;pl++) // do nothing
/* loop over all bounding Primitive variables */
#define PBOUNDLOOP(pb) for(pb=0;pb<NPRBOUND;pb++)

/* loop over all dumped Primitive variables */
#define PDUMPLOOP(pd) for(pd=0;pd<NPRDUMP;pd++)
// loop over all inversion Primitive variables -- only over 5 quantities (not field)
#define PINVERTLOOP(pi) for(pi=0;pi<NPRINVERT;pi++)

#define NUMPRIMLOOP(pl) for(pl=0;pl<numprims;pl++)
#define NUMPRIMLOOPGEN(pl) NUMPRIMLOOP(pl,pl)



#endif




// always goes over all conserved
#define PALLLOOP(pl) for(pl=0;pl<NPR;pl++)


// always goes over all conserved
#define PALLREALLOOP(pl) for(pl=0;pl<NPRREALSET;pl++)


// goes over all/any npr lists for copying one list to another
#define PMAXNPRLOOP(pl) for(pl=0;pl<MAXNPR;pl++)

// always goes over all standard invertable quantities for inversion to operate normally
#define PLOOPALLINVERT(pl) for(pl=0;pl<=B3;pl++)

// always goes over all interpolatable primitives
#define PLOOPALLINTERP(pl) for(pl=0;pl<NPR2INTERP;pl++)


// always goes over all primitives
#define PDIAGLOOP(pl) PALLLOOP(pl)





#define PLOOPNOB1(pl) for(pl=0;pl<B1;pl++)
#define PLOOPBONLY(pl) for(pl=B1;pl<=B3;pl++)
#define PLOOPNOB2(pl) for(pl=B3+1;pl<NPR;pl++)



/* loop over all Dimensions; second rank loop */
#define DLOOP(j,k) for(j=0;j<NDIM;j++)for(k=0;k<NDIM;k++)
/* loop over all Dimensions; first rank loop */
#define DLOOPA(j) for(j=0;j<NDIM;j++)
/* loop over all Space dimensions; second rank loop */
#define SLOOP(j,k) for(j=1;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* loop over all Space dimensions; first rank loop */
#define SLOOPA(j) for(j=1;j<NDIM;j++)
/* loop over all for j and Space for k; second rank loop */
#define DSLOOP(j,k) for(j=0;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* space-space */
#define SSLOOP(j,k) for(j=1;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* loop over all for k and Space for j; second rank loop */
#define SDLOOP(j,k) for(j=1;j<NDIM;j++)for(k=0;k<NDIM;k++)

// loop over directions
#define DIRLOOP(dir) for(dir=0;dir<COMPDIM*2;dir++)

#define DIRSIGNLOOP(dirsign) for(dirsign=-1;dirsign<=1;dirsign+=2)


#define DIMENLOOP(dir) for(dir=1;dir<=COMPDIM;dir++)





/////////////////////////////////////
/////////////////////////////////////
//
// PER-POINT LOOPS (MORE RELATED TO NOT BEING PRIMITIVE OR CONSERVED)
//
/////////////////////////////////////
/////////////////////////////////////


// loop over fail flag in boundary code
#define FBOUNDLOOP(ff) for(ff=0;ff<NUMPFLAGSBOUND;ff++)

// loop over jet regions
#define JETLOOP(jetio) for(jetio=0;jetio<NUMJETS;jetio++)

// loop over ener/flux regions
#define ENERREGIONLOOP(enerregion) for(enerregion=0;enerregion<NUMENERREGIONS;enerregion++) if(dothisenerreg[enerregion])

// loop over ALL ener/flux regions
#define ENERREGIONALLLOOP(enerregion) for(enerregion=0;enerregion<NUMENERREGIONS;enerregion++)

// loop over fair/floor types
#define FLOORLOOP(floor) for(floor=0;floor<NUMFAILFLOORFLAGS;floor++)

// loop over debug time scales
#define TSCALELOOP(tscale) for(tscale=0;tscale<NUMTSCALES;tscale++)


// loop over ALL sources (including geometry)
#define SCLOOP(sc) for(sc=0;sc<NUMSOURCES;sc++)

// loop over all sources EXCEPT geometry (assumes GEOMSOURCE==0 or at least nothing before GEOMSOURCE matters)
#define SCPHYSICSLOOP(sc) for(sc=GEOMSOURCE+1;sc<NUMSOURCES;sc++)

// loop over fluxterms
#define FLLOOP(fl) for(fl=0;fl<NUMFLUXTERMS;fl++)


// loop over pflag flags
#define PFLAGLOOP(pf) for(pf=0;pf<NUMPFLAGS;pf++)

// for USEMPI&&USEROMIO==1
#define ROMIOCOLLOOP(romiocoliter) for(romiocoliter=0;romiocoliter<romiocloopend;romiocoliter++)

#define BUFFERINIT nextbuf=0
#define COLINIT nextcol=0

// for mpicombie==0
#define COLLOOP(coliter) for(coliter=0;coliter<numfiles;coliter++)


#define DTSTAGELOOP(dtstage) for(dtstage=0;dtstage<MAXDTSTAGES;dtstage++)

#define INTERPLOOP(interpi) for(interpi=0;interpi<NUMINTERPTYPES;interpi++)

#define ENODEBUGLOOP(enodebugi) for(enodebugi=0;enodebugi<NUMENODEBUGS;enodebugi++)



