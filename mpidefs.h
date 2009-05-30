


#ifndef USINGLIAISON
// for only grmhd code
#include "supermpidefs.h"
#endif


// for both grmhd and liaison codes whether or not DOINGGRMHDLIAISON==1
#include "mpidefs.mpi_grmhd_grray_liaison.h"




int romiocoliter;
int periodicx1, periodicx2, periodicx3;
int mpiperiodicx1, mpiperiodicx2, mpiperiodicx3;
int skipix1, reflectix1, reflectox1;
int skipix2, reflectix2, reflectox2;
int skipix3, reflectix3, reflectox3;
int intix1, intox1, intix2, intox2, intix3, intox3;
int skipintix1, skipintix2, skipintix3;
int ncpux1, ncpux2, ncpux3;
int truenumprocs;
int myid, myid_world, numprocs;
char myidtxt[MAXFILENAME];
int totalzones, realtotalzones;
int rtotalzones;
int itotalzones;
int sizes[COMPDIM + 1][MAXCPUS];
int isizes[COMPDIM + 1][MAXCPUS];
int totalsize[COMPDIM + 1];
int itotalsize[COMPDIM + 1];
int mycpupos[COMPDIM + 1];		// my position amongst the cpus
int dirset[NUMBOUNDTYPES][COMPDIM*2][DIRNUMVARS];
int srdir[3*2];			// which direction this cpu
				// sends/receives normal interior data
int startpos[COMPDIM + 1];
int endpos[COMPDIM + 1];		// startj and endj are where this CPU
				// located on full grid 
int *startpos0[COMPDIM+1];
int *endpos0[COMPDIM+1];
int *mycpupos0[COMPDIM+1];

int plmpiglobal;


int procnamelen;
#if(USEMPI)
MPI_Group MPI_GROUP_WORLD;
char processor_name[MPI_MAX_PROCESSOR_NAME];
MPI_Status mpichstatus;
#endif

// MPI transmit vars, so minimum local code changes
FTYPE ndtsend, bsq_maxsend;

// for data output
int nextbuf,numcolumns;
int bufferoffset;
int joniosize,writebufsize;


