// large memory items for mpidefs.h
#if(USEMPI)
FTYPE workbca[2][COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM];	// [1=out/2=in][0=right,2=up,1=left,3=down,4=out,5=in][datawidth]
PFTYPE workbc_inta[2][COMPDIM * 2][NUMPFLAGSBOUND * NBIGBND * NBIGSM];	// [1=out/2=in][0=right,2=up,1=left,3=down,4=out,5=in][datawidth]
#endif
FTYPE (*workbc)[COMPDIM * 2][NMAXBOUND * NBIGBND * NBIGSM];
PFTYPE (*workbc_int)[COMPDIM * 2][NUMPFLAGSBOUND * NBIGBND * NBIGSM];
