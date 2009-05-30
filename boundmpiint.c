#include "decs.h"


int bound_mpi_int(int boundstage, int boundvartype, PFTYPE prim[][N2M][N3M][NUMPFLAGS])
{
  int dir;

#if(USEMPI)
  /* These arrays contain designations that identify 
   * each recv and send */
  static MPI_Request requests[COMPDIM * 2 * 2];
  // format of map for requests[dir*2+recv/send(0/1)]
#endif

#if(USEMPI)

  ///////////////
  // dir=1  
  if((boundstage==STAGE0)||(boundstage==STAGEM1)){
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[boundvartype][dir][DIRIF]) pack_int(dir,boundvartype,prim,workbc_int);
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[boundvartype][dir][DIRIF]) sendrecv_int(dir,boundvartype,workbc_int,requests);
  }
  if((boundstage==STAGE1)||(boundstage==STAGEM1)){
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[boundvartype][dir][DIRIF]) recvwait(dir,requests);
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[boundvartype][dir][DIRIF]) unpack_int(dir,boundvartype,workbc_int,prim);
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[boundvartype][dir][DIRIF]) sendwait(dir,requests);
  }

  ///////////////
  // dir=2
  if((boundstage==STAGE2)||(boundstage==STAGEM1)){
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[boundvartype][dir][DIRIF]) pack_int(dir,boundvartype,prim,workbc_int);
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[boundvartype][dir][DIRIF]) sendrecv_int(dir,boundvartype,workbc_int,requests);
  }
  if((boundstage==STAGE3)||(boundstage==STAGEM1)){
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[boundvartype][dir][DIRIF]) recvwait(dir,requests);
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[boundvartype][dir][DIRIF]) unpack_int(dir,boundvartype,workbc_int,prim);
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[boundvartype][dir][DIRIF]) sendwait(dir,requests);
  }

  ///////////////
  // dir=3
  if((boundstage==STAGE4)||(boundstage==STAGEM1)){
    for(dir=X3UP;dir<=X3DN;dir++) if(dirset[boundvartype][dir][DIRIF]) pack_int(dir,boundvartype,prim,workbc_int);
    for(dir=X3UP;dir<=X3DN;dir++) if(dirset[boundvartype][dir][DIRIF]) sendrecv_int(dir,boundvartype,workbc_int,requests);
  }
  if((boundstage==STAGE5)||(boundstage==STAGEM1)){
    for(dir=X3UP;dir<=X3DN;dir++) if(dirset[boundvartype][dir][DIRIF]) recvwait(dir,requests);
    for(dir=X3UP;dir<=X3DN;dir++) if(dirset[boundvartype][dir][DIRIF]) unpack_int(dir,boundvartype,workbc_int,prim);
    for(dir=X3UP;dir<=X3DN;dir++) if(dirset[boundvartype][dir][DIRIF]) sendwait(dir,requests);
  }


  // now corner zones will be filled correctly
  // GODMARK: If made fixup_utoprim() and check_solution() not use corner zones could bound all directions at once -- probably not important performance hit


  // end if mpi
#endif

  return(0);

}	
// end function


// PACKLOOP allows one to alter which i or j is faster iterated
#define PACKLOOP_INT(i,j,k,istart,istop,jstart,jstop,kstart,kstop) GENLOOP(i,j,k,istart,istop,jstart,jstop,kstart,kstop) FBOUNDLOOP(pl)

// packs data for shipment
void pack_int(int dir, int boundvartype,PFTYPE prim[][N2M][N3M][NUMPFLAGS],PFTYPE workbc_int[][COMPDIM * 2][NUMPFLAGSBOUND * NBIGBND * NBIGSM])
{
  // dir=direction sending
  int i,j,k;
  int pl;
  int bci;


  bci=0;
  PACKLOOP_INT(i,j,k
	   ,dirset[boundvartype][dir][DIRPSTART1]
	   ,dirset[boundvartype][dir][DIRPSTOP1]
	   ,dirset[boundvartype][dir][DIRPSTART2]
	   ,dirset[boundvartype][dir][DIRPSTOP2]
	   ,dirset[boundvartype][dir][DIRPSTART3]
	   ,dirset[boundvartype][dir][DIRPSTOP3]
	       ){
    /*
    if(bci>=dirset[boundvartype][dir][DIRSIZE]){
      dualfprintf(fail_file,"pack memory leak: bci: %d dirset[%d][DIRSIZE]: %d\n",bci,dirset[boundvartype][dir][DIRSIZE]);
      myexit(10);
    }
    */
    workbc_int[PACK][dir][bci++] = prim[i][j][k][pl];
  }
}

#if(USEMPI)
void sendrecv_int(int dir, int boundvartype,PFTYPE workbc_int[][COMPDIM * 2][NUMPFLAGSBOUND * NBIGBND * NBIGSM],MPI_Request *requests)
{
  MPI_Irecv(workbc_int[UNPACK][dir],
	    dirset[boundvartype][dir][DIRSIZE],
	    MPI_PFTYPE,
	    dirset[boundvartype][dir][DIROTHER],
	    dirset[boundvartype][dir][DIRTAGR],
	    MPI_COMM_GRMHD,
	    &requests[dir*2+REQRECV]);

  MPI_Isend(workbc_int[PACK][dir],
	    dirset[boundvartype][dir][DIRSIZE],
	    MPI_PFTYPE,
	    dirset[boundvartype][dir][DIROTHER],
	    dirset[boundvartype][dir][DIRTAGS],
	    MPI_COMM_GRMHD,
	    &requests[dir*2+REQSEND]);

}
#endif


void unpack_int(int dir, int boundvartype,PFTYPE workbc_int[][COMPDIM * 2][NUMPFLAGSBOUND * NBIGBND * NBIGSM],PFTYPE prim[][N2M][N3M][NUMPFLAGS])
{
  // dir is direction receiving from
  int i,j,k;
  int pl;
  int bci;

  bci=0;
  PACKLOOP_INT(i,j,k
	   ,dirset[boundvartype][dir][DIRUSTART1]
	   ,dirset[boundvartype][dir][DIRUSTOP1]
	   ,dirset[boundvartype][dir][DIRUSTART2]
	   ,dirset[boundvartype][dir][DIRUSTOP2]
	   ,dirset[boundvartype][dir][DIRUSTART3]
	   ,dirset[boundvartype][dir][DIRUSTOP3]
	       ){
    prim[i][j][k][pl]=workbc_int[UNPACK][dir][bci++];
  }
}

