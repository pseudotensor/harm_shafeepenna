// whether to check symmetry along diagonal (for test==151 in Sasha's tests)
#define ASYMDIAGCHECK 0

// whether to use accurate (bit more expensive) diagnostics.  Only necessary for timeorder>2
#define ACCURATEDIAG 1

// whether to use sin/cos exactly symmetrized around pi/2
// if one doesn't do this, then \theta=\pi pole can lead to \detg=small number instead of 0.0
#define ACCURATESINCOS 1

// whether to try cold inversion if hot fails
#define HOT2COLD 0


// *NUMBER* OF DIMENSIONS FOR COMPUTATION
// Can choose 3 and still do 1D optimized problems
// Choosing 2 or even 1 reduces some memory load for those things that are accessed by using 1,2,3 instead of arbitrarily accessed
// not sure if 1 works.
#define COMPDIM 3

// whether to have get_geometry() point to global gcov,gcon,gcovpert -- faster and no limitations
#define GETGEOMUSEPOINTER 1

// whether to split NPR PLOOPs
#define SPLITNPR 0



// whether need higher order memory stuff
#define HIGHERORDERMEM 1

// whether to split off MA and EM parts of flux for higher order
// for now always do if doing higher order
// GODMARK: If SPLITMAEM==0 and pointfield method, then only EMF terms have fixed stencil
//#define SPLITMAEMMEM (HIGHERORDERMEM)
#define SPLITMAEMMEM (0) // CHANGINGMARK

//#define DOMONOINTERP NOMONOINTERP
//#define DOMONOINTERP JMONOINTERP
#define DOMONOINTERP SMONOINTERP

// whether to split isotropic pressure off of flux for separate interpolation
// this splitting only used ultimately if splitmaem==1 (even if computed separately)
#define SPLITPRESSURETERMINFLUXMA 0 // CHANGINGMARK
// below not setup yet:
#define SPLITPRESSURETERMINFLUXEM 0 // CHANGINGMARK


// whether to do a2c during c2e
// see interpline.c
#define MERGEDC2EA2CMETHOD 0


#define HIGHERORDERROUGH 0
#define HIGHERORDERSMOOTH 1
#define HIGHERORDERVERYSMOOTH 2

//////////////////////////////////////////////////////////////////
//
// CHOOSE SMOOTHNESS EXPECTED OF SOLUTION for HIGHERORDER TERMS:
//
//////////////////////////////////////////////////////////////////
#define HIGHERORDERTYPE HIGHERORDERSMOOTH



#if(HIGHERORDERTYPE==HIGHERORDERROUGH)// non-smooth version:
// for pointfield method doesn't keep divb=0 properly



#define DO_SPLITC2E NOSPLITA2C
#define DO_SPLITA2C ENERGY_CONTROLS_ALL_WEIGHTS // c2a and a2c for conserved quants
#define DO_SPLITA2CSMOOTH DO_SPLITA2C
#define DO_SPLITA2C4FLUX DO_SPLITA2C // more generally
#define DO_SPLITA2C4MAFLUX MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS // couples weights together
#define DO_SPLITA2C4SMOOTHFLUX DO_SPLITA2C
#define DO_SPLITSOURCE NOSPLITA2C // c2a for ENOSOURCETERM

#define EMFLUXFIELDTERMTYPE ENOFLUX
#define FIELDTERMTYPE ENOFLUX
#define FIELDFVTYPE ENOCONSERVED



#elif(HIGHERORDERTYPE==HIGHERORDERSMOOTH)  // smooth version


// should be used with pointfield method for general problems
// SMOOTH-named variables used for magnetic field operations

#define DO_SPLITC2E NOSPLITA2C
#define DO_SPLITA2C ENERGY_CONTROLS_ALL_WEIGHTS // c2a and a2c for conserved quants
#define DO_SPLITA2CSMOOTH CONSTANT_ALL_WEIGHTS
//#define DO_SPLITA2C4FLUX CONSTANT_ALL_WEIGHTS // needed for TESTNUMBER==49 if doing FLUXRECON -- GODMARK: perhaps good for smooth problems
//#define DO_SPLITA2C4FLUX NOSPLITA2C // so pressure changes don't control other terms (for example pressure can be constant but field or velocity can vary) -- causes problems for test=106 since E-M not coupled enough when correction occurs
//#define DO_SPLITA2C4FLUX ENERGY_IS_ALL_WEIGHTS
//#define DO_SPLITA2C4FLUX MINIMIZE_ALL_WEIGHTS
//#define DO_SPLITA2C4FLUX ENERGY_CONTROLS_ALL_WEIGHTS
// CHANGINGMARK
#define DO_SPLITA2C4FLUX MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS
#define DO_SPLITA2C4MAFLUX DO_SPLITA2C4FLUX
//#define DO_SPLITA2C4FLUX MINIMIZE_ALL_WEIGHTS
#define DO_SPLITA2C4SMOOTHFLUX CONSTANT_ALL_WEIGHTS
#define DO_SPLITSOURCE NOSPLITA2C // c2a for ENOSOURCETERM

#define EMFLUXFIELDTERMTYPE ENOSMOOTHFLUX
#define FIELDTERMTYPE ENOSMOOTHFLUX
#define FIELDFVTYPE ENOSMOOTHCONSERVED



#elif(HIGHERORDERTYPE==HIGHERORDERVERYSMOOTH) // supersmooth version
// for problems with no very strong discontinuities (often good even for strong discontinuities)



#define DO_SPLITC2E NOSPLITA2C
#define DO_SPLITA2C CONSTANT_ALL_WEIGHTS
#define DO_SPLITA2CSMOOTH CONSTANT_ALL_WEIGHTS
#define DO_SPLITA2C4FLUX CONSTANT_ALL_WEIGHTS
#define DO_SPLITA2C4MAFLUX CONSTANT_ALL_WEIGHTS
#define DO_SPLITA2C4SMOOTHFLUX CONSTANT_ALL_WEIGHTS
#define DO_SPLITSOURCE NOSPLITA2C

#define EMFLUXFIELDTERMTYPE ENOSMOOTHFLUX
#define FIELDTERMTYPE ENOSMOOTHFLUX
#define FIELDFVTYPE ENOSMOOTHCONSERVED

#endif





///////////////
// for MASSENERGYMOMENTUM_IS_COUPLED_WEIGHTS
// used in reconstruct.c and interpline.mono.c
#define VELTERMSMINIMIZE(pl) (pl==RHO || pl==UU || pl==UU+dir)
#define ORTHOVEL1TERMSMINIMIZE(pl) (pl==UU+odir1)
#define ORTHOVEL2TERMSMINIMIZE(pl) (pl==UU+odir2)
#define PRESSUREMINIMIZE(pl) (pl==FLUXSPLITPMA(dir)) // pressure term when splitmaem==1 and doing MA term
#define EMFTERMS(pl) (pl>=B1 && pl<=B3) // only used if splitmaem==0

// ALLOTHERS corresponds to all those pl's not specially handled by the minimization routine.  EMFTERMS only handled specially if emffixedstencil==1
// below is always true unless hit certain terms that are true
// 
#define ALLOTHERSMINIMIZE(pl) (! (VELTERMSMINIMIZE(pl) || ORTHOVEL1TERMSMINIMIZE(pl) || ORTHOVEL2TERMSMINIMIZE(pl) || PRESSUREMINIMIZE(pl) || ((emffixedstencil==1)&&EMFTERMS(pl))  ) )
      //#define OTHERADVECTMINIMIZE(pl) (pl==YL && DOYL!=DONOYL || pl==YNU && DOYNU!=DONOYNU || pl==ENTROPY && DOENTROPY!=DONOENTROPY)



// whether to use staggered fields (which requires memory for it)
#define FIELDSTAGMEM 1

// whether to store get_state() data into array for later use since expensive to compute get_geometry() and ucon_calc() type things
#define STOREFLUXSTATE 1

// whether to define Toth memory stuff
//#define FIELDTOTHMEM (!FIELDSTAGMEM) // CHANGINGMARK
#define FIELDTOTHMEM 1 // CHANGINGMARK


// whether to track and dump A_i the vector potential
#define TRACKVPOT 0

// whether to redefine fields from A_i during evolution
// currently only done after each full timestep to avoid being expensive
// only machine error different from evolution of field itself
#define EVOLVEWITHVPOT TRACKVPOT // choice

// whether to specify gdet at end when setting EMF or to have internal to variables before averaging.
// point is that gdet at end is probably better, esp. at coordinate singularities.
#define CORNGDETVERSION 1
//#define CORNGDETVERSION 0


// size of global arrays
// ALSO CONTROLS WHETHER MEMORY IS ALLOCATED FOR SPECIAL 1D/2D/3D directional arrays
// e.g. if N3==1, then second dimension is neglected in all calculations

#define N1	64		// number of zones in 1-direction
#define N2	64       	// number of zones in 2-direction
#define N3      1               // number of zones in 3-direction

//#if(DOENOFLUX==ENOFINITEVOLUME)
#define MAXBND 5
//#else
//#define MAXBND 3
//#endif


// choose metric
#define MCOORD KSCOORDS


// defines how one forms the EM stress-energy tensor
// GENMAXWELL is perfectly fine for relativistic problems that are not near machine error for highly relativistic flows that usually is only for tests
//#define MAXWELL GENMAXWELL
// prim version avoids catastrophic cancellation in non-rel limit
#define MAXWELL PRIMMAXWELL




#define DO_VORTICITY_IMAGE 0


#define PRODUCTION 0
// 0: full images, dumps, etc., few more key failure stderr/fail_file messages
// 1: only log density image since too many images (takes alot of time), no utoprim failure messages -- assume debug.out and debug???? will have info if needed


// 0: normal computational zones outputted on diagnostics
// else, FULLOUTPUT # of extra boundary zones on each side (if allowed by dimensionality)
// If # of requested boundary zones is larger than real #, then real # used

// ONLY CAN BE USED WITH USEMPI==0
#define FULLOUTPUT 0


#define MAILWHENDONE 1
#define MAILFROMREMOTE 0
#define REMOTEHOST "bh.astro.uiuc.edu"
#define EMAILADDRESS "jmckinney@cfa.harvard.edu"
#define EMAILMESSAGE "Done with GRMHD run DEFAULT"



// whether doing super long doubles
#define SUPERLONGDOUBLE 0


#define PERFTEST 0
// 0: don't perform performance test
// 1: do

#define DOAVG 0
// 0: don't do time average dumps, so don't allocate memory
// 1: do

// GODMARK: only setup for full Pi theta grid, not Pi/2 or any other theta grid.
// e.g. will fail for defcoord=3 Pi/2 grid.  Causes code to crash
#define DOJETDIAG 1
// 0: don't do jet diagnostics (ener, not flener that is always done)
// 1: do

#define DOAVG2 0 // only make 1 if DOAVG 1 above
// 0: don't split AVG file
// 1: do

#define DODEBUG 1
// 0: don't output debug dump file or ener file(ener is based on dump counts)
// 1: do


// whether to dump vector potential
#define DOVPOTDUMP (1 && TRACKVPOT)

#define DOENODEBUG 0
// whether to do ENO debug




#define DOGRIDSECTIONING 0
// 0: don't do grid sectioning
// 1: do (only part of the grid in the active section is evolved)



#define WENO_USE_PRIM_REDUCTION 1
// 0: do not perform the c2a limiting on the fluxes
// 1: perform the limiting of c2a flux reconstruction based on the limiting of the a2c primitive correction:
//    if the a2c reconstruction of the conserved quantities leads to very different primitives, then no a2c reconstruction is done; in this case
//    no reconstruction is done on the fluxes either.

#define WENO_EXTRA_A2C_MINIMIZATION 0 // CHANGINGMARK

#define WENO_REDUCE_A2C_LOOK_OTHER_DIRECTIONS 1
// 0: do not perform additional reduction
// 1: perform additional reduction:
//    for NDIM > 1, subsequent 1d a2c reconstructions reduce to lower order (=no reconstruction) if any of the previous reconstructions reduced

#define LIMIT_FLUXC2A_PRIM_CHANGE 0
// 0: do not limit primtive correction due to c2a done on fluxes
// 1: limit c2a correction
 
#define DO_WENO_DEBUG 0
// whether to debug WENO

#define DOSUPERDEBUG 0
// 0: don't output super debug output
// 1: do

// whether to evolve metric value of M and a and Q
#define DOEVOLVEMETRIC 0

// 0 = evolve only infrequently after many steps
// 1 = evolve on every substep in full RK form
// 2 = use accdt and ndt1/2/3 to determine whether one should evolve metric on substeps or longsteps
#define EVOLVEMETRICSUBSTEP 0

// whether to limit timestep using source term at t=0 and full flux term for future times
#define LIMITDTWITHSOURCETERM 0

// whether to limit source update by some constraint (gives dUriemann to source())
#define LIMITSOURCES 0


// whether to use gravitydt to limit dt
#define USEGRAVITYDTINDTLIMIT 0


// 0: no restriction beyond normal ghost+active bounds
// 1: don't limit via accdt or gravitydt inside horizon
// 2: in addition to #1, don't limit due to velocity inside horizon (can't be used if boundary condition drives solution and would limit dt more than active region)
#define RESTRICTDTSETTINGINSIDEHORIZON 1


#define DODISS 0

// see diag_source()
#define DOLUMVSR 0

// see diag_source()
#define DODISSVSR 0


// see metric.c
#define DOSELFGRAVVSR 0

#define DOFIELDLINE 1
// 0: don't output energy@infinity and field line stuff
// 1: do

// whether to use Roe-averaged state in determining diffusive term in flux
#define ROEAVERAGEDWAVESPEED 0
// whether to use Athena's version of Roe averaging and estimating wave speeds
#define ATHENAROE 0

// whether to store wave speeds over whole grid before computing flux
// useful to avoid extra calculations if doing "local" LAXF or "global" LAXF.
// default HARM was using VERY local LAXF (only wavespeeds from primitives interpolated to the edge).
#define STOREWAVESPEEDS 0

// whether to use stored wave speeds for flux calculation (allows one to store wave speeds for interp.c but use true VERYLOCALVCHAR that is vchar's estimated from boundary as in standard HARM -- rather than maximum from center zones as done by STORED version of VERYLOCALVCHAR)
// silly choice would be 0 if VCHARTYPE=GLOBALVCHAR since interp.c doesn't use the stored wave speeds if this is the case.  So shouldn't store in this case since no point.
#define USESTOREDSPEEDSFORFLUX (STOREWAVESPEEDS) // choice really independent of STOREWAVESPEEDS, but generall normally want to couple them

#define VCHARTYPE VERYLOCALVCHAR

// whether to check on inversion and report problem
#define CHECKONINVERSION 1

#define PRECISEINVERSION 1
// whether we use ultimately precise inversion or "workable" inversion precision (see utoprim's for how defined)

#define WHICHVEL VELREL4
// which velocity to compute in (init can be anything (see init.c's for transformations)


// if WHICHEOM==WITHNOGDET, then below determines which EOMs get what geometric prefactor.  Notice (as described in phys.c's source_conn() ) that geometry issue applies AFTER additions/subtractions of EOMs (as done by REMOVERESTMASSFROMUU).
#define REMOVERESTMASSFROMUU 1
// whether to subtract rest-mass from energy equation for divT=0 equation of motion
// 0: use MHD stress tensor with explicit rest-mass included
// 1: as 0, but subtract out rest-mass from conservation and flux terms for evolution
// 2: use MHD stress tensor withOUT rest-mass
// this changes mhd_calc() in phys.c and assumes rest of code uses mhd stress-energy tensor without restmass also!!
// DIAGNOSTICS also without rest-mass for UU terms.

#define RELTYPE RELEOM

#define EOMTYPE EOMGRMHD
//#define EOMTYPE EOMCOLDGRMHD

// whether to try other methods for the inversion if primary choices fails
// created because utoprim_2d_final fails for large b^2/rho when other methods (even utoprim_2d) do not fail.
// seems that in some cases causes code to stall while 1D method takes forever on many iterations
#define UTOPRIMTRYAGAIN 0


// which EOS (EoS) (equation of state) to use
#define WHICHEOS IDEALGAS
//#define WHICHEOS MIGNONE
//#define WHICHEOS GRBPWF99
//#define WHICHEOS KAZFULL



// whether to coevolve the entropy to check for shock+reconnection dissipation
#define DOENTROPY DONOENTROPY // normal total energy equation
//#define DOENTROPY DOEVOLVECOMPAREENTROPY // coevolve and compare
//#define DOENTROPY DOEVOLVEDIRECTENTROPY // directly evolve entropy equation instead of total energy equation


// generically which type of entropy evolution to do, except not used when doing DOENTROPY==DOEVOLVEDIRECTENTROPY
#define WHICHENTROPYEVOLVE EVOLVESIMPLEENTROPY

// whether to call fixup() after initialization
#define FIXUPAFTERINIT 1

// whether to call fixup() after restart
#define FIXUPAFTERRESTART 1

// checks if solution is reasonble and fails a point if not
// only checks if b^2/\rho>>BSQORHOMAX, etc.
#define CHECKSOLUTION 1



// factor by which zone-zone values can be different for \gamma and internal energy, used of CHECKSOLUTION==1
#define GAMMAPERCDIFFMAX 2.0
#define UPERCDIFFMAX 10.0


// whether to interpolate extra quantity (e.g. v^2) and use to renormalize velocities after they are interpolated
#define DOEXTRAINTERP 0
// must also set RESCALEINTERP=1


// whether to evolve Y_l (see kazeosfull.c)
#define DOYL 0
// whether to evolve Y_\nu (see kazeosfull.c)
#define DOYNU DOYL


// whether to control the limiter
#define LIMADJUST LIMITERFIXED
//#define LIMADJUST LIMITERBSQORHOANDU

// whether to limit all variables or just hydro variables
#define HYDROLIMADJUSTONLY 1

// determine how flux is computed per zone
#define FLUXADJUST FLUXFIXED
//#define FLUXADJUST FLUXBSQORHOANDU

// whether to change flux calculation for all variables or just hydro variables
#define HYDROFLUXADJUSTONLY 1

// whether to allow negative internal energy in substep
// UTOPRIMFAILRETURNTYPE==UTOPRIMRETURNADJUSTED should be set in global.h
#define STEPOVERNEGU NEGDENSITY_FIXONFULLSTEP
// seems to make things worse (more failures, but also worse failures)

#define STEPOVERNEGRHO NEGDENSITY_FIXONFULLSTEP

#define STEPOVERNEGURHO NEGDENSITY_FIXONFULLSTEP




//#define UTOPRIMADJUST UTOPRIMSTATIC
#define UTOPRIMADJUST UTOPRIMAVG

// 0: return non-updated quantities
// 1: if u<0 or rho<0, then return updated quantities, else return non-updated quantities
#define UTOPRIMFAILRETURNTYPE UTOPRIMRETURNADJUSTED

// whether to smooth singularity of black hole over grid scale
#define SMOOTHSING 1

#define COORDSINGFIX 1
// whether to move polar axis to a bit larger theta
// theta value where singularity is displaced to

#define COORDSINGFIXCYL 0 //whether perform the same fix for CYLMINKMETRIC's z axis

//#define SINGSMALL (1E-3)
#define SINGSMALL (1E-20)
// Hawley uses 0.06283 (0.02Pi)

#define VOLUMEDIFF 0
// whether to use volume regularization or not

#define DOSTOREPOSITIONDATA 1
// whether to store X,V,dxdxp in case are expensive to compute

// which derivative type to use when computing connection coefficients
#define CONNDERTYPE DIFFGAMMIE

// whether to correct connection so body forces are 0 when p=constant
// only makes sense for non-higher order scheme
#define CONNMACHINEBODY 0

#define WHICHEOM WITHGDET
//#define WHICHEOM WITHNOGDET


// whether to include memory and fill with gdetvol_func()
#define GDETVOLDIFF 0

// whether to make gdet finite volume-like when reverting to 1D spherical polar coordinates
#define FIXGDETSPC_WHEN_1DRADIAL 1


#define MINDT 1.e-20 // minimum dt

#define JONCHECKS 1 // for vel=3 extra checks (standard things)
#define JONCHECKS2 1 // check_pr() crazy thing for vel3 (crazy u^t check and fix -- 2D only)


#define FLOORDIAGS 1 // whether to compute floor diagnostics

// whether to use analytic g^{\mu\nu} if available
#define ANALYTICGCON 0

#define ANALYTICCONNECTION 0 // whether to use analytic connection
// only applies to certain metric's and coordinates, see set_grid.c

#define ANALYTICSOURCE 0 // whether to use analytic source function
// very slow with gcc, some extend pgcc, not a bit problem with icc
// only applies to certain metric's and coordaintes, see phys.c

// whether to (with ENO method) avoid boundary zones if doing outflow on that boundary
#define OUTFLOWAVOIDBC 0

#define FLUXDIMENSPLIT PERFECTUNSPLIT
#define A2CDIMENSPLIT PERFECTUNSPLIT

// whether to bound flux if doing fluxrecon
#define BOUNDFLUXRECON 0

// whether to bound unew (GODMARK: no separate algorithm yet, so only works for non-outflow boundary conditions) -- only needed if unewisavg==1
#define BOUNDUNEW 0

// whether to have dq's : only used by second order methods
#define DODQMEMORY 1

// whether to adjust interpolation based upon boundary conditions
// defunct.  Not used anymore by ENO scheme and shouldn't be activated for other schemes.
// LEAVE 0!
#define BOUNDARYINTERPADJUST 0

// whether to compute \dot{F_r} as in diag.c
//disable the analysis if N1 <= 1 because leads to segfault if F1 == NULL (as is in the case of N1 == 1)
#define COMPUTEFRDOT 0

// whether to compute faraday and currents and output them to dumps
#define CALCFARADAYANDCURRENTS 1

#define WHICHCURRENTCALC CURRENTCALC1


// whether want faraday at t=0 for dump
#define FARADAYT0 1

// whether want partial currents for t=0 dump
#define CURRENTST0 1

// whether to flip gdet sign over coordinate singularities
#define FLIPGDETAXIS 0

#define BOUNDPLPR 0
#define NOFLUXCTONX1DN 0


// boundary condition mnemonics
// can be reset/added to by user init.h
#define OUTFLOW	0
#define SYMM	1
#define ASYMM	2
#define FIXED	3
#define POLARAXIS 4
#define FIXEDOUTFLOW 5 // means fixed inflow but allows outflow -- basically outflow if no inflow, but if inflow then set values to fixed quantities
#define NSSURFACE 6 // whatever in bounds.c for NS surface
#define PERIODIC 7 // periodic boundary conditions
#define OUTFLOWNOINFLOW 8  //copy if velocity directed outward; force velocity to zero and copy if directed inward (into the grid from the boundary)
#define RAMESHOUTFLOW 9 //OUTFLOW quantities according the Ramesh's power-law solution
#define BCEXTRAP 100  //atch: extrapolation into the boundary with high order (can be specified)
#define CYLAXIS 101   //atch: cylindrical axis BC
#define BCEXTRAP_VEL3 102  //atch: the same as BCEXTRAP but extrapolates 3-velocity rather than 4-velocity; important for Hubble flow
#define JETINJECTION 103 //Fixed boundary condition for jet injection surrounded by an outflow condition
#define BCU1EXTRAPOTHERFIXED 104
#define BCEXTRAPCONSTRAINED 105
#define RESCALEOUTFLOW 106
#define RESCALEFIXEDALLOUTFLOWU1 107
#define FIXED_RESCALEOUTFLOWU1 108
#define BONDIMDOTOUTFLOW 109
#define BONDIINTOUTFLOW 110
#define DISKSURFACE 111
#define FREEOUTFLOW 112
#define R0SING 113


#define EVOLVECHECKS 0
// whether to check boundary conditions and limit gamma during advance()

// whether and which type of fixups to be used
#define FIXUPZONES FIXUP1ZONE

/** algorithmic choices **/

#define HLLBOUNDARY 0
// use HLL on boundary fluxes


#define FIXUPFLUX 0
// fix up the flux using fix_flux() in fixup.c
// caution...
// GODMARK: Not sure how metric and how emf works around axes.  Not unlike \Omega and such things that are "symmetric" while v^\phi might be antisymmetric?


#define ZEROOUTFLOWFLUX 0
// 0: don't do anything special, let boundary conditions control fux
// 1: zero the x1-velocity at the outflow boundary if inflow in flux routine so interpolation is overridden
// seems to cause problems at outer edge in internal energy

// seems to cause problems with EOMFFDE at boudaries.
// GODMARK



// assumes inner and outer theta boundaries are at 0 and Pi.
#define ZEROPOLEFLUX 0
// 0: don't do anything special
// 1: zero polar theta flux since should be 0

// seems to cause problems at 0 and pi boundaries with FFDE, or even just one boundary???
// GODMARK

// REASON: Cannot just set F=0.  The flux of the kinetic energy still has the pressure term even if the velocity into the wall is 0!



// whether to rescale interpolation
#define RESCALEINTERP 0
// 0: don't rescale
// 1: do rescale

#define BDIRCONT 1
// 0: don't
// 1: make field along flux direction continuous

// whether to use hyperbolic function for cmin and cmax instead of sharp max(0,l,r)
#define HYPERHLL 0

// whether to shift stencil inside horizon
#define HORIZONSUPERFAST 0



//////////////////////////////////
//
// which variable to interpolate
//
/////////////////////////////////


//#define VARTOINTERP PRIMTOINTERP_VSQ
#define VARTOINTERP PRIMTOINTERP



/////////////////////////////////////////////////////////
//
// STUFF FOR INTERPLINE.C
//
/////////////////////////////////////////////////////////

// \detg rho u^i used or not (then just u^i or primvel^i)
#define VLINEWITHGDETRHO 1

// whether Pline in interpline.c for PARALINE without weno has b^2/2 term in pressure (should!)
#define PLINEWITHFIELD 1

// whether PARALINE uses MONO or not
#define PARALINEUSESMONO 0

// whether to use appropriately centered primitive
// GODMARK: averaging is diffusive, so diffuses shock -- may make indicator weak
#define USEAVGPRIMITIVEFORWENOFLAT 1


// whether to invert average conserved quantity to get a primitive to use as shock indicator
// then correctly positioned in time, but more diffused value
// may be a bad idea since inversion is expensive, may find no solution, or may lead to negative internal energy in shock.
#define USEPRIMITIVEFROMAVGCONSERVED 0


// whether to use contact discontinuity indicator to steepen contacts in interpline.c
#define CONTACTINDICATOR 0

// send Sasha dP
#define COMPUTEDRHODP 1


// whether to reduce to 1-point stencil in case of superfast divergence at that point
#define SUPERFASTDIVREDUCE 0 //atch

// minimum preferred order
#define MINPREFORDER 3

// whether to compute shock indicator
#define SHOCKINDICATOR 1

// interp general stuff
#define WHICHPARA PARA4


#define BONDI_BOUNDARY_SET_PL_PR 0  //do not analytically set p_l & p_r at the outer boundary for the Bondi problem


#define NOSPECIALFIELD 0
#define PULSARFIELD 1
#define PULSARFIELD2 2
#define PULSARFIELD3 3

#define VARTOINTERPFIELD NOSPECIALFIELD

#define NUMPANALYTICOTHER 0
#define DODUMPOTHER 0 // whether to dump other stuff

// whether to do flux dumps and have memory for it
#define FLUXDUMP 0

// default number of extra things needed
#define NUMPOTHER 0


