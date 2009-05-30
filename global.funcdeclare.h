//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// function declarations
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


extern int main(int argc, char *argv[]);

extern int init(int *argc, char **argv[]);
extern void parainitchecks(void);
extern void myargs(int argc, char *argv[]);


// SECTIONMARK:
extern int init_gridsectioning(void);
extern int bound_gridsectioning(int ispstag, FTYPE (*prim)[N2M][N3M][NPR]);
extern int findandsetactivesection(int timeorder, int numtimeorders, long int thenstep, FTYPE thetime );
extern int findindexfromradius(FTYPE xr, int *xcpupos1, int *xi);
extern int setactivesection(int (*abs)[NDIM], int doprintout);


// stepping
extern int step_ch_full(void);
extern void get_truetime_fluxdt(int numtimeorders, SFTYPE localdt, FTYPE (*CUf)[4], FTYPE (*Cunew)[4], SFTYPE *fluxdt, SFTYPE *boundtime, SFTYPE *tstepparti, SFTYPE *tsteppartf);

extern void set_normal_realisinterp(int *realisinterp);

extern int fluxcalc(int stage, FTYPE pr[][N2M][N3M][NPR],
		    FTYPE F1[][N2M][N3M][NPR], 
		    FTYPE F2[][N2M][N3M][NPR], 
		    FTYPE F3[][N2M][N3M][NPR], 
		    FTYPE CUf,
		    FTYPE fluxdt,
		    FTYPE *ndt1,
		    FTYPE *ndt2,
		    FTYPE *ndt3
		    );
extern void diag_source_comp(struct of_geom *ptrgeom, FTYPE (*dUcomp)[NPR],SFTYPE Dt);
extern void diag_source_all(struct of_geom *ptrgeom, FTYPE *dU,SFTYPE Dt);
extern int diag_flux(FTYPE prim[][N2M][N3M][NPR], FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR],FTYPE F3[][N2M][N3M][NPR],SFTYPE Dt);
extern int compute_new_metric_substep(FTYPE pb[][N2M][N3M][NPR],FTYPE *CUf, FTYPE *Cunew);  
extern int diss_compute(int evolvetype, int inputtype, FTYPE *U, struct of_geom *ptrgeom, FTYPE *prbefore, FTYPE *pr);
extern int advance(int stage, FTYPE pi[][N2M][N3M][NPR],FTYPE pb[][N2M][N3M][NPR], FTYPE pf[][N2M][N3M][NPR],
		   FTYPE ui[][N2M][N3M][NPR],FTYPE uf[][N2M][N3M][NPR], FTYPE ucum[][N2M][N3M][NPR],
		   FTYPE *CUf,FTYPE *Cunew,SFTYPE fluxdt, SFTYPE boundtime, int timeorder, int numtimeorders, FTYPE *ndt);


extern int avg2cen_interp(int *locpl, int *whichpltoavg,  int *ifnotavgthencopy, int whichquantity, int whichavg2cen, FTYPE (*prims_from_avg_cons)[N2M][N3M][NPR], FTYPE (*in)[N2M][N3M][NPR], FTYPE (*out)[N2M][N3M][NPR]);



extern void set_defaults_performance_checks(void);
extern void set_file_versionnumbers(void);

extern int timecheck(int whichlocation, SFTYPE comptstart);
extern int gocheck(int whichlocation);
extern int output_steptimedt_info(SFTYPE comptstart);


extern int error_check(int wherefrom);
extern int find_horizon(int fromwhere);

// initialize DUMP stuff
extern int init_dumps(void);
extern int init_linklists(void);
int setuplinklist(int numcolumns,int which);
extern struct blink * addlink(struct blink * clinkptr);

// ENER file stuff
extern int dump_ener(int doener, int dordump, int call_code);

extern int diag(int call_code, FTYPE time, long localnstep, long localrealnstep);

extern void frdotout(void);

extern void report_systeminfo(FILE * fileout);
extern int IsLittleEndian(void);
extern void *SwapEndian(void* Addr, const int Nb);

extern void makedirs(void);

extern void appendener(FILE* ener_file,SFTYPE pdot_tot[][NPR],SFTYPE*fladd_tot,SFTYPE*sourceadd_tot);

extern void divbmaxavg(FTYPE p[][N2M][N3M][NPR],FTYPE*ptrdivbmax,FTYPE*ptrdivbavg);
extern void gettotal(int numvars, SFTYPE* vars[],int*sizes,SFTYPE*vars_tot[]);
extern void gettotali(int numvars, int* vars[],int*sizes,int*vars_tot[]);
extern int constotal(int enerregion, SFTYPE *vars_tot);
extern int integrate(int numelements, SFTYPE * var,SFTYPE *var_tot,int type, int enerregion);

extern int counttotal(int enerregion, CTYPE *vars_tot, int num);
extern int integratel(int numelements, CTYPE * var,CTYPE *var_tot,int type, int enerregion);


// DUMP file stuff
extern int isenoughfreespace(unsigned long long need);

extern int dump_gen(int readwrite, long dump_cnt, int bintxt, int whichdump,MPI_Datatype datatype, char *fileprefix, char *fileformat, char *filesuffix, int (*headerfun) (int bintxt, FILE*headerptr),int (*content) (int i, int j, int k, MPI_Datatype datatype, void*setbuf));

extern int header1_gen(int readwrite, int bintxt, int bcasthead, void *ptr, size_t size, char *format, size_t nmemb, MPI_Datatype datatype, FILE *stream);


extern int dump(long dump_cnt);
extern int dump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);
extern int dump_header(int bintxt, FILE *headerptr);

extern int avgdump(long avg_cnt);
extern int avg_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int avgdump2(long avg_cnt);
extern int avg2_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int debugdump(long debug_cnt);
extern int debug_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int gdump(long gdump_cnt);
extern int gdump_content(int i, int j, int k, MPI_Datatype datatype, void *writebuf);

extern int fieldlinedump(long fieldline_cnt);
extern int fieldline_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int dissdump(long dump_cnt);
extern int dissdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int dumpother(long dump_cnt);
extern int dumpother_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);


extern int fluxdumpdump(long dump_cnt);
extern int fluxdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int eosdump(long dump_cnt);
extern int eosdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int vpotdump(long dump_cnt);
extern int vpotdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern int fakedump(void);
extern int fakedump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);
extern int fakedump_header(int bintxt, FILE *headerptr);


extern int image_dump(long image_cnt);
extern int imagedefs(int whichk, int scale, int limits, int vartype);
extern int image(long dump_cnt, int whichk, int scale, int limits, int vartype);
extern int image_header(int bintxt, FILE *headerptr);
extern int image_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);
extern void prminmaxsum(FTYPE p[][N2M][N3M][NPR], int start,int nmemb, FTYPE *max, FTYPE*min,FTYPE*sum);

extern int restart_init(int which);
extern int restart_init_checks(int which);


// restart dump
extern int restart_read(long which);
extern int check_fileformat(int readwrite, int bintxt, int whichdump, int numcolumns, int docolsplit, int mpicombine, int sizeofdatatype, FILE *stream);
extern int read_restart_header(int bintxt, FILE* headerptr);
extern int restart_read_defs(void);
extern int rdump_read_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);
extern int restart_write(long dump_cnt);
extern int write_restart_header(int bintxt, FILE* headerptr);
extern int rdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

// old metric restart dump
extern int restartmetric_read(long which);
extern int read_restartmetric_header(int bintxt, FILE* headerptr);
extern int rmetricdump_read_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);
extern int restartmetric_write(long dump_cnt);
extern int write_restartmetric_header(int bintxt, FILE* headerptr);
extern int rmetricdump_content(int i, int j, int k, MPI_Datatype datatype,void *writebuf);

extern void myfopen(char*fname, char*fmt, char*message, FILE ** fileptr);
extern void myfclose(FILE ** fileptr,char*message);

extern void myset(MPI_Datatype datatype, void *ptr, int start, int nmemb, void*writebuf);
extern void myget(MPI_Datatype datatype, void *ptr, int start, int nmemb, void*writebuf);

extern void myfwrite(int bintxt, MPI_Datatype datatype, void *ptr, int start, int nmemb, int i, int j, int k, FILE**stream,void*writebuf);

extern void myfread(int bintxt, MPI_Datatype datatype, void *ptr, int start, int nmemb, int i, int j, int k, FILE**stream,void*writebuf);



// initialize stuff
// specific to init.c's and used in initbase.c and init.c, so leave global
extern int post_init_specific_init(void);
extern int pre_init_specific_init(void);
extern int prepre_init_specific_init(void);
extern int init_consts(void);
extern int init_grid(void);
extern int init_global(void);

extern int init_primitives(FTYPE p[][N2M][N3M][NPR]);

extern int init_vpot(FTYPE p[][N2M][N3M][NPR]);
extern int vpot2field(FTYPE A[][N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3],FTYPE p[][N2M][N3M][NPR]);
extern int update_vpot(int stage, FTYPE pr[][N2M][N3M][NPR],FTYPE (*ptrfluxvec[NDIM])[N2M][N3M][NPR], FTYPE CUf);
extern int normalize_field_withnorm(FTYPE norm);
extern int assign_fieldconservatives_pointvalues(FTYPE U[][N2M][N3M][NPR]);
extern void setfdivb(FTYPE *divb, FTYPE (*p)[N2M][N3M][NPR],FTYPE (*U)[N2M][N3M][NPR], int i, int j, int k);

extern int copy_vpot2flux(FTYPE (*A)[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3]);
extern int evolve_withvpot(FTYPE (*prim)[N2M][N3M][NPR]);

extern int init_vpot_user(int *whichcoord, int l, int i, int j, int k, FTYPE p[][N2M][N3M][NPR], FTYPE *V, FTYPE *A);
extern int init_vpot2field_user(FTYPE A[][N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3],FTYPE pr[][N2M][N3M][NPR]);

extern int transform_primitive_vB(int whichvel, int whichcoord, int i,int j, int k, FTYPE p[][N2M][N3M][NPR], FTYPE pstag[][N2M][N3M][NPR]);
extern int transform_primitive_pstag(int whichvel, int whichcoord, int i,int j, int k, FTYPE p[][N2M][N3M][NPR], FTYPE pstag[][N2M][N3M][NPR]);
extern int init_zero_field(FTYPE p[][N2M][N3M][NPR]);

extern int pi2Uavg(int *fieldfrompotential, FTYPE (*prim)[N2M][N3M][NPR], FTYPE (*Upoint)[N2M][N3M][NPR], FTYPE (*Uavg)[N2M][N3M][NPR]);

extern void set_default_nprlists(void);

extern int addremovefieldfromnpr(int doadd, int *whichpltoavg, int *ifnotavgthencopy, int interptype, int dir, int *nprlocalstart, int *nprlocalend, int *nprlocallist, FTYPE (*current_in)[N2M][N3M][NPR], FTYPE (*current_out)[N2M][N3M][NPR]);
extern int addremovefromnpr(int doadd, int *whichpltoavg, int *ifnotavgthencopy, int *nprlocalstart, int *nprlocalend, int *nprlocallist, FTYPE (*in)[N2M][N3M][NPR], FTYPE (*out)[N2M][N3M][NPR]);

extern int addremovefromanynpr(int doadd, int *whichpltoavg, int *ifnotavgthencopy, int *anynprstart, int *anynprend, int *anynprlist, int *nprlocalstart, int *nprlocalend, int *nprlocallist, FTYPE (*in)[N2M][N3M][NPR], FTYPE (*out)[N2M][N3M][NPR]);

// called in restart.c and initbase.c
extern void set_grid(int whichtime,FTYPE *CUf, FTYPE *Cunew);

extern int higherorder_set(int whichquantity, int recontype, int*weightsplittype);

extern int get_fluxpldirs(int *Nvec, int dir, int *fluxdir, int* pldir, int *plforflux, FTYPE *signflux);
extern void get_odirs(int dir,int *odir1,int *odir2);
extern int set_location_fluxasemforvpot(int dir, int *numdirs, int *odir1, int *odir2, int *loc);
extern int get_numdirs_fluxasemforvpot(int *numdirs, int *fieldloc);

extern int plstart_set(int whichquantity, int dir, int recontype, int *plstart);


// some physics

extern int sourcephysics(FTYPE *ph, struct of_geom *geom, struct of_state *q, FTYPE *Ugeomfree, FTYPE *dUother, FTYPE (*dUcomp)[NPR]);

extern void postdt(void);
extern int primtoU(int returntype, FTYPE *p, struct of_state *q, struct of_geom *geom,
		   FTYPE *U);

extern int ucon_calc_3vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon);
extern int ucon_calc_rel4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon);
extern int ucon_calc_4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon);
extern int ucon_calc_4vel_bothut(FTYPE *pr, struct of_geom *geom, FTYPE *ucon, FTYPE *ucon2);

#if(RELTYPE==RELEOM)

#if(WHICHVEL==VEL4)
#define ucon_calc ucon_calc_4vel
#define dudp_calc dudp_calc_gen
#elif(WHICHVEL==VEL3)
#define ucon_calc ucon_calc_3vel
#define dudp_calc dudp_calc_3vel
#elif(WHICHVEL==VELREL4)
#define ucon_calc ucon_calc_rel4vel
#define dudp_calc dudp_calc_gen

#elif(RELTYPE==NONRELEOM) // not really right
#define ucon_calc ucon_calc_nonrel
#define dudp_calc dudp_calc_nonrel
#endif

#endif

extern int ucon_calcother(FTYPE *pr, FTYPE *ucon);
extern void ucon_precalc(FTYPE *ucon, FTYPE *AA, FTYPE *BB,
			 FTYPE *CC, FTYPE *discr);



extern FTYPE ranc(int seed);

extern FTYPE interpn( int order, FTYPE x_eval,  FTYPE x1, FTYPE f1, FTYPE x2, FTYPE f2, FTYPE x3, FTYPE f3, FTYPE x4, FTYPE f4, FTYPE x5, FTYPE f5, FTYPE x6, FTYPE f6 );


// fixup stuff

extern int check_pr(FTYPE *pr, FTYPE *prmodel, struct of_geom *geom, int modelpos,int finalstep);
extern int ucon_fix(FTYPE disc, FTYPE AA, FTYPE BB, FTYPE CC,
		    FTYPE *ucon);

/* // dudp stuff */

/* extern void dutdui_calc(FTYPE *ucon, FTYPE *dutdui); */
/* extern void duiduj_calc(FTYPE *ucon, FTYPE *dutdui); */
/* extern void dbtdui_calc(FTYPE *dutdui, FTYPE *pr, FTYPE *dbtdui); */
/* extern void dbiduj_calc(FTYPE *dbtdui, FTYPE *dutdui, FTYPE *ucon, */
/* 			FTYPE *b, FTYPE dbiduj[][NDIM]); */
/* extern void db2dui_calc(FTYPE dbiduj[][NDIM], FTYPE *b, */
/* 			FTYPE *db2dui); */
/* extern void duudud_calc(FTYPE *ucon, FTYPE duudud[][NDIM]); */

/* extern void dbsqdui_calc(FTYPE dbiduj[][NDIM], FTYPE *b, */
/* 			 FTYPE *dbsqdui); */
/* extern void dgdvi_calc(FTYPE *pr,FTYPE *dgdvi); */
/* extern void duidvj_calc(FTYPE *dgdv,FTYPE duidvj[][NDIM]); */
/* extern void dudduu_calc(FTYPE*dutdui, FTYPE dudduu[][NDIM]); */
/* extern void dbdiduj_calc(FTYPE dbiduj[][NDIM],FTYPE dbdiduj[][NDIM]); */
/* extern void ducon_dv3_calc(struct of_state *q,FTYPE ducon_dv[][NDIM]); */
extern int sp_stress_calc(FTYPE *pr, FTYPE tens_matt[][NDIM],
			  FTYPE tens_em[][NDIM], FTYPE *b,
			  FTYPE *ucon);



// log file stuff
extern void myfprintf(FILE* fileptr, char *format, ...);
extern void dualfprintf(FILE* fileptr,char *format, ...);
extern void logsfprintf(char *format, ...);
extern void trifprintf(char *format, ...);

// boundary stuff
extern void set_boundloop(int boundvartype, int *inboundloop, int*outboundloop, int*innormalloop, int*outnormalloop, int (*inoutlohi)[NUMUPDOWN][NDIM], int *riin, int *riout, int *rjin, int *rjout, int *rkin, int *rkout);
extern int report_bound_loop(void);
extern void set_numbnd(int boundvartype, int *numbnd, int *numnpr);

// below are for particular purposes
extern int bound_allprim(int boundstage, SFTYPE boundtime, FTYPE prim[][N2M][N3M][NPR]);
extern int bound_evolveprim(int boundstage, SFTYPE boundtime, FTYPE prim[][N2M][N3M][NPR]);
extern int bound_prim(int boundstage, SFTYPE boundtime, FTYPE prim[][N2M][N3M][NPR]);
extern int bound_pstag(int boundstage, SFTYPE boundtime, FTYPE prim[][N2M][N3M][NPR]);
extern int bound_beforeevolveprim(int boundstage, SFTYPE boundtime, FTYPE prim[][N2M][N3M][NPR]);

// below can choose boundvartype
extern int bound_anyallprim(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE prim[][N2M][N3M][NPR]);
extern int bound_anyprim(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE prim[][N2M][N3M][NPR]);
extern int bound_anypstag(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE prim[][N2M][N3M][NPR]);
extern int bound_uavg(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE uavg[][N2M][N3M][NPR]);

// only pflag doesn't have boundvartype
extern int bound_pflag(int boundstage, SFTYPE boundtime, PFTYPE primbase[][N2M][N3M][NUMPFLAGS]);


extern int bound_flux(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR], FTYPE F3[][N2M][N3M][NPR]);

// user bounds:
extern int bound_prim_user_dir(int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE prim[][N2M][N3M][NPR]);
extern int bound_pstag_user_dir(int boundstage, SFTYPE boundtime, int whichdir, int boundvartype, FTYPE pstag[][N2M][N3M][NPR]);
extern int bound_prim_user_after_mpi_dir(int boundstage, SFTYPE boundtime, int whichdir, FTYPE prim[][N2M][N3M][NPR]);
extern int bound_flux_user(int boundstage, SFTYPE boundtime, int boundvartype, FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR], FTYPE F3[][N2M][N3M][NPR]);
extern int bound_pflag_user(int boundstage, SFTYPE boundtime, int boundvartype, PFTYPE prim[][N2M][N3M][NUMPFLAGS]);



extern int inflow_check_4vel(int dir, FTYPE *pr, struct of_geom *ptrgeom, int finalstep);
extern int inflow_check_3vel(int dir, FTYPE *pr, struct of_geom *ptrgeom, int finalstep);
extern int inflow_check_rel4vel(int dir, FTYPE *pr, struct of_geom *ptrgeom, int finalstep);


// transform stuff
extern int bl2met2metp2v(int whichvel, int whichcoord, FTYPE *pr, int ii, int jj, int kk);
extern int bl2met2metp2v_genloc(int whichvel, int whichcoord, FTYPE *pr, int ii, int jj, int kk, int loc);
extern int bl2met2metp2v_gen(int whichvel, int whichcoord, int newwhichvel, int newwhichcoord, FTYPE *pr, int ii, int jj, int kk);

extern int ucov_whichcoord2primecoords(int whichcoord, int ii, int jj, int kk, int loc, FTYPE *ucov);

extern int metp2met2bl(int whichvel, int whichcoord, FTYPE *pr, int ii, int jj, int kk);
extern int pr2ucon(int whichvel, FTYPE *pr, struct of_geom *geom, FTYPE*ucon);
extern int coordtrans(int whichcoordin, int whichcoordout, int ii, int jj, int kk, int loc, FTYPE*ucon);
extern void bltoks(int ii, int jj, int kk, int loc, FTYPE*ucon);
extern void kstobl(int ii, int jj, int kk, int loc, FTYPE*ucon);
extern void mettometp(int ii, int jj, int kk, FTYPE*ucon);
extern void metptomet(int ii, int jj, int kk, FTYPE*ucon);
extern void mettometp_genloc(int ii, int jj, int kk, int loc, FTYPE*ucon);
extern void metptomet_genloc(int ii, int jj, int kk, int loc, FTYPE*ucon);
extern void mettometp_simple(FTYPE idxdxp[][NDIM], FTYPE*ucon);
extern void metptomet_simple(FTYPE dxdxp[][NDIM], FTYPE*ucon);
extern void metptomet_ucov_simple(FTYPE idxdxp[][NDIM], FTYPE*ucon);
extern void mettometp_ucov_simple(FTYPE dxdxp[][NDIM], FTYPE*ucon);
extern void metptomet_Tud(int ii, int jj, int kk, FTYPE Tud[][NDIM]);
extern void metptomet_simple_Tud(FTYPE dxdxp[][NDIM], FTYPE idxdxp[][NDIM], FTYPE Tud[][NDIM]);
extern void ucon2pr(int whichvel, FTYPE *ucon, struct of_geom *geom, FTYPE *pr);
extern int vcon2pr(int whichvel, FTYPE *vcon, struct of_geom *geom, FTYPE *pr);



// metric stuff
extern void gset_genloc(int getprim, int whichcoord, int i, int j, int k, int loc, struct of_geom *geom);
extern void gset(int getprim, int whichcoord, int i, int j, int k, struct of_geom *geom);
extern FTYPE gdet_func(int whichcoord, FTYPE gcov[][NDIM]);
extern FTYPE gdet_func_singcheck(int whichcoord, FTYPE *V, FTYPE gcov[][NDIM]);
extern void gdetvol_func(struct of_geom *ptrgeom, FTYPE (*gdet)[N2M+SHIFT2][N3M+SHIFT3][NPG], FTYPE *gdetvol);
//extern FTYPE bl_gdet_func(FTYPE r, FTYPE th);
//extern void bl_gcov_func(FTYPE r, FTYPE th, FTYPE gcov[][NDIM]);
//extern void bl_gcon_func(FTYPE r, FTYPE th, FTYPE gcon[][NDIM]);
extern void conn_func(int whichcoord, FTYPE *X, struct of_geom *geom,
		      FTYPE lconn[][NDIM][NDIM],FTYPE *conn2);
extern void mks_unitheta_idxvol_func(int i, int j, int k, FTYPE *idxvol);

//extern void gcov_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM]);
//extern void gcon_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE gcon[][NDIM]);
extern void matrix_inverse(int whichcoord, FTYPE gcov[][NDIM], FTYPE gcon[][NDIM]);
extern void alphalapse_func(struct of_geom *ptrgeom, int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE gcon[][NDIM], FTYPE *alphalapse);



// coordinate stuff
extern void set_coord_parms(void);
extern void write_coord_parms(void);
extern void read_coord_parms(void);
extern void coord(int i, int j, int k, int loc, FTYPE *X);
extern void coord_ijk(int i, int j, int k, int loc, FTYPE *X);
extern void coord_free(int i, int j, int k, int loc, FTYPE *X);

extern void bl_coord(FTYPE *X, FTYPE *V);
extern void bl_coord_ijk(int i, int j, int k, int loc, FTYPE *V);
extern void bl_coord_ijk_2(int i, int j, int k, int loc, FTYPE *X, FTYPE *V);

extern void dxdxprim(FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM]);
extern void dxdxprim_ijk(int i, int j, int k, int loc, FTYPE (*dxdxp)[NDIM]);
extern void dxdxprim_ijk_2(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE (*dxdxp)[NDIM]);

extern void idxdxprim(FTYPE (*dxdxp)[NDIM], FTYPE (*idxdxp)[NDIM]);
extern void idxdxprim_ijk(int i, int j, int k, int loc, FTYPE (*idxdxp)[NDIM]);
extern void idxdxprim_ijk_2(struct of_geom *ptrgeom, FTYPE *X, FTYPE *V, FTYPE (*idxdxp)[NDIM]);


extern int setihor(void);
extern FTYPE setRin(int ihor);

extern int is_inside_surface(int dir, int ii, int jj, int kk, int pp);
extern int is_on_surface(int dir, int ii, int jj, int kk, int pp);

extern FTYPE mysign(FTYPE x);
extern FTYPE myfabs(FTYPE x);

extern FTYPE mysin(FTYPE th);
extern FTYPE mycos(FTYPE th);


// eos stuff
extern int pickeos_eomtype(int whicheos, int whicheom);

extern FTYPE pressure_rho0_u(FTYPE rho0, FTYPE u);
extern FTYPE u_rho0_p(FTYPE rho0, FTYPE p);
extern FTYPE pressure_rho0_w(FTYPE rho0, FTYPE w);
extern FTYPE dpdu_rho0_u(FTYPE rho0, FTYPE u);
extern FTYPE dpdrho0_rho0_u(FTYPE rho0, FTYPE u);
extern FTYPE cs2_compute(FTYPE rho0, FTYPE u);
extern FTYPE compute_entropy(FTYPE rho0, FTYPE u);
extern FTYPE compute_u_from_entropy(FTYPE rho0, FTYPE entropy);
extern FTYPE pressure_wmrho0(FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_idwmrho0dp(FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_idrho0dp(FTYPE rho0, FTYPE wmrho0);
extern FTYPE compute_qdot(FTYPE rho0, FTYPE u);
extern int compute_sources_EOS(FTYPE *pr, struct of_geom *geom, struct of_state *q, FTYPE *Ui, FTYPE *dUother, FTYPE(*dUcomp)[NPR]);
extern void compute_allextras(int justnum, FTYPE rho0, FTYPE u,int *numextrasreturn,FTYPE*extras);
extern int get_extrasprocessed(int doall, int i, int j, int k, FTYPE *pr, FTYPE *extras, FTYPE *processed);
extern FTYPE compute_temp(FTYPE rho0, FTYPE u);
extern void compute_EOS_parms(FTYPE (*prim)[N2M][N3M][NPR]);
extern void store_EOS_parms(int i, int j, int k, int numparms, FTYPE *parlist);
extern void get_EOS_parms(int i, int j, int k, int*numparms, FTYPE *parlist);



// physics stuff
extern int set_zamo_velocity(int whichvel, struct of_geom *ptrgeom, FTYPE *pr);
extern int set_zamo_ucovuconplus1ud0(struct of_geom *ptrgeom, FTYPE *ucov, FTYPE *ucon, FTYPE *plus1ud0);
extern int set_zamo_ucon(struct of_geom *ptrgeom, FTYPE *ucon);
extern FTYPE contract(FTYPE *vcon, FTYPE *wcon);

extern int bsq_calc(FTYPE *pr, struct of_geom *geom, FTYPE *b2);
extern void b_calc(FTYPE *pr, FTYPE *ucon, FTYPE *b);


extern int gamma_calc(FTYPE *pr, struct of_geom *geom,FTYPE *gamma);


extern int dudp_calc_gen(int whichcons,FTYPE *pr, struct of_state *q,
		     struct of_geom *geom, FTYPE **alpha);

extern int dudp_calc_3vel(int whichcons,FTYPE *pr, struct of_state *q,
		     struct of_geom *geom, FTYPE **alpha);


extern int sol(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmax, FTYPE *vmin);

extern void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
extern void UtoU_evolve2diag(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);



extern int Utoprimgen(int finalstep, int evolvetype, int inputtype, FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr);
extern int Utoprimloop(FTYPE unew[][N2M][N3M][NPR],FTYPE pf[][N2M][N3M][NPR]);
extern int primtoUloop(FTYPE pi[][N2M][N3M][NPR],FTYPE unew[][N2M][N3M][NPR]);

extern int Utoprim(int entropyeom, FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr);
extern int Utoprim_ldz(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr);
extern int Utoprim_1d(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr);
extern int Utoprim_1d_opt(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr);
extern int Utoprim_2d(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr);
extern int Utoprim_1d_final(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr);
extern int Utoprim_2d_final(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr);
//extern int Utoprim_2d_final_nonrelcompat_inputnorestmass(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr);  //wrong function name, corrected by atch, see below
extern int Utoprim_jon_nonrelcompat_inputnorestmass(int whicheom, FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr);
extern int Utoprim_5d2_final(FTYPE *U, struct of_geom *geom, PFTYPE *lpflag, FTYPE *pr);




extern void tetr_func(FTYPE tetr_cov[][NDIM], FTYPE tetr_con[][NDIM]);
extern void get_geometry(int i, int j, int k, int loc, struct of_geom *geom);
extern void get_geometry_gdetonly(int ii, int jj, int kk, int pp, struct of_geom *geom);
extern void get_geometry_geomeonly(int ii, int jj, int kk, int pp, struct of_geom *geom);

extern void get_allgeometry(int i, int j, int k, int loc, struct of_allgeom *allgeom, struct of_geom *geom);
extern void set_igeom(struct of_geom *geom);
extern void set_igeomsimple(struct of_geom *geom);



extern int get_state(FTYPE *pr, struct of_geom *geom,struct of_state *q);
extern int get_stateforsource(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);
extern int get_stateforfluxcalc(int dimen, int isleftright, FTYPE *pr, struct of_geom *ptrgeom, struct of_state **qptr);
extern int get_stateforUdiss(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q);


extern int primtoflux(int returntype, FTYPE *pa, struct of_state *q, int dir,
	       struct of_geom *geom, FTYPE *fl);
extern int primtoflux_splitmaem(int returntype, FTYPE *pa, struct of_state *q, int fluxdir, int fundir, struct of_geom *geom, FTYPE *flma, FTYPE *flem);


extern int flux_compute_general(int i, int j, int k, int dir, struct of_geom *geom, FTYPE CUf, FTYPE *p_c, FTYPE *p_l, FTYPE *p_r, FTYPE *F, FTYPE *ctop);
extern int flux_compute_splitmaem(int i, int j, int k, int dir, struct of_geom *geom, FTYPE CUf, FTYPE *p_c, FTYPE *p_l, FTYPE *p_r, FTYPE *F, FTYPE *FEM, FTYPE *ctop);
extern int get_global_wavespeeds(int dir,FTYPE *pr,struct of_geom *ptrgeom, FTYPE *wspeedtemp);



extern int deconvolve_flux(int dir, struct of_geom *ptrgeom, FTYPE *p_c, FTYPE *p_l, FTYPE *p_r, struct of_state *state_c, struct of_state *state_l, struct of_state *state_r, FTYPE *F_c, FTYPE *F_l, FTYPE *F_r, FTYPE *U_c, FTYPE *U_l, FTYPE *U_r);
extern int deconvolve_flux_ma(int dir, struct of_geom *ptrgeom, FTYPE *p_c, FTYPE *p_l, FTYPE *p_r, struct of_state *state_c, struct of_state *state_l, struct of_state *state_r, FTYPE *F_c, FTYPE *F_l, FTYPE *F_r, FTYPE *U_c, FTYPE *U_l, FTYPE *U_r);
extern int deconvolve_flux_em(int dir, struct of_geom *ptrgeom, FTYPE *p_c, FTYPE *p_l, FTYPE *p_r, struct of_state *state_c, struct of_state *state_l, struct of_state *state_r, FTYPE *F_c, FTYPE *F_l, FTYPE *F_r, FTYPE *U_c, FTYPE *U_l, FTYPE *U_r);

extern int get_wavespeeds(int dir, struct of_geom *ptrgeom, FTYPE *p_l, FTYPE *p_r, FTYPE *U_l, FTYPE *U_r, FTYPE *F_l, FTYPE *F_r, struct of_state *state_l, struct of_state * state_r, FTYPE *cminmax_l, FTYPE *cminmax_r, FTYPE *cminmax, FTYPE *ctop);



extern void mks_source_conn(FTYPE *ph, struct of_geom *ptrgeom,
		     struct of_state *q,FTYPE *dU);
extern int source(FTYPE *pa, struct of_geom *geom, struct of_state *q, FTYPE *Ui, FTYPE *dUriemann,
		  FTYPE (*Uacomp)[NPR], FTYPE *Ua);

extern FTYPE taper_func(FTYPE R,FTYPE rin) ;
extern FTYPE rhor_calc(int which);
extern FTYPE rmso_calc(int which) ;
extern FTYPE uphi_isco_calc(int which,FTYPE r);

extern int set_atmosphere(int whichcond, int whichvel, struct of_geom *geom, FTYPE *pr);
extern int set_density_floors_default(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *scaler);
extern int set_density_floors(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *scaler);
extern int fixup(int stage,FTYPE (*var)[N2M][N3M][NPR],int finalstep);
extern int fixup1zone(FTYPE *pr,struct of_geom *ptrlgeom, int finalstep);
extern int diag_fixup(FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, int finalstep,int whocalled);

extern int superdebug(FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, int whocalled);

extern int limit_gamma(FTYPE gammamax, FTYPE*pr,struct of_geom *geom, int finalstep);

extern int fixup_checksolution(int stage, FTYPE (*pv)[N2M][N3M][NPR],int finalstep);
extern int fixup_utoprim(int stage, FTYPE (*pv)[N2M][N3M][NPR],FTYPE (*pbackup)[N2M][N3M][NPR], int finalstep);
extern int post_fixup(int stage, SFTYPE boundtime, FTYPE (*pv)[N2M][N3M][NPR],FTYPE (*pbackup)[N2M][N3M][NPR],int finalstep);
extern int pre_fixup(int stage, FTYPE (*pv)[N2M][N3M][NPR]);
extern int get_bsqflags(int stage, FTYPE (*pv)[N2M][N3M][NPR]);


extern void tet_func(FTYPE metr[][NDIM], FTYPE tetr[][NDIM]);
extern int dsyev_(char *jobz, char *uplo, int *n, FTYPE *a, int *lda,
		  FTYPE *w, FTYPE *work, int *lwork, int *iwork,
		  int *liwork, int *info);
extern int fail(int fail_type);
extern void setfailresponse(int restartonfail);
extern int setflux(void);
extern int sethorizonflux(void);
extern int settrueglobalregion(void);
extern int settrueglobalwithbndregion(void);
extern int setjetflux(void);
extern void setrestart(int*appendold);


extern int vchar(FTYPE *pr, struct of_state *q, int dir,
		 struct of_geom *geom, FTYPE *cmax, FTYPE *cmin,int *ignorecourant);
extern FTYPE chk_disp(FTYPE v);
extern void make_co_to_comov(FTYPE *ucon, FTYPE ecov[][NDIM],
			     FTYPE econ[][NDIM]);
extern void transform(FTYPE *vec, FTYPE t[][NDIM]);
extern void coeff_set(FTYPE rho, FTYPE u);
extern void transform(FTYPE *ucon, FTYPE t[][NDIM]);

extern void mhd_calc(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd);
extern void mhd_calc_0(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd);

extern void mhd_calc_em(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd);
extern void mhd_calc_ma(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd, FTYPE *mhddiagpress);

extern int area_map(int call_code, int type, int size, int i, int j, int k, FTYPE prim[][N2M][N3M][NPR]);
extern void bcon_calc(FTYPE *pr, FTYPE *ucon, FTYPE *ucov,
		      FTYPE *bcon);

extern FTYPE lc4(int updown, FTYPE detg, int mu,int nu,int kappa,int lambda);
extern void faraday_calc(int which, FTYPE *b, FTYPE *u, struct of_geom *geom, FTYPE faraday[][NDIM]);
extern void current_precalc(int which, struct of_geom *geom, struct of_state *q, SFTYPE Dt,FTYPE faraday[][3]);
extern void init_varstavg(void);
extern void final_varstavg(FTYPE IDT);
extern int set_varstavg(FTYPE tfrac);
extern void current_calc(FTYPE cfaraday[][N2M][N3M][NUMCURRENTSLOTS][3]);
extern int current_doprecalc(int which, FTYPE p[][N2M][N3M][NPR]);
extern int average_calc(int doavg);



// interpolation stuff
extern int get_loop(int pointorlinetype, int interporflux, int dir, struct of_loop *loop);
extern int set_interpalltypes_loop_ranges(int pointorlinetype, int interporflux, int dir, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be);


// line types:
extern void set_interp_loop_gen(int withshifts, int interporflux, int dir, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be);
extern void set_interp_loop(int withshifts, int interporflux, int dir, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be);
extern void set_interp_loop_expanded(int withshifts, int interporflux, int dir, int *intdir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk, int *bs, int *ps, int *pe, int *be);

// point types:
extern int set_interppoint_loop_ranges(int interporflux, int dir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk);
extern void set_interppoint_loop(int interporflux, int dir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk);
extern void set_interppoint_loop_expanded(int interporflux, int dir, int *is, int *ie, int *js, int *je, int *ks, int *ke, int *di, int *dj, int *dk);


extern int vpot2field_useflux(int *fieldloc,FTYPE pfield[][N2M][N3M][NPR],FTYPE ufield[][N2M][N3M][NPR]);
extern int vpot2field_centeredfield(FTYPE (*A)[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3],FTYPE pfield[][N2M][N3M][NPR],FTYPE ufield[][N2M][N3M][NPR]);
extern int vpot2field_staggeredfield(FTYPE (*A)[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3],FTYPE pfield[][N2M][N3M][NPR],FTYPE ufield[][N2M][N3M][NPR]);
extern int interpolate_ustag2fieldcent(int stage, SFTYPE boundtime, int timeorder, int numtimeorders, FTYPE preal[][N2M][N3M][NPR],FTYPE pstag[][N2M][N3M][NPR],FTYPE ucent[][N2M][N3M][NPR],FTYPE pcent[][N2M][N3M][NPR]);
extern int vectorpot_fluxreconorfvavg(int stage, FTYPE pr[][N2M][N3M][NPR], FTYPE (*A)[N1M+SHIFT1][N2M+SHIFT2][N3M+SHIFT3]);
extern int deaverage_fields_fv(FTYPE (*primreal)[N2M][N3M][NPR], FTYPE (*in)[N2M][N3M][NPR], FTYPE (*out)[N2M][N3M][NPR]);
extern int field_integrate_fluxrecon(int stage, FTYPE pr[][N2M][N3M][NPR], FTYPE quasifield[][N2M][N3M][NPR], FTYPE pointfield[][N2M][N3M][NPR]);
extern int vectorpot_useflux(int stage, FTYPE pr[][N2M][N3M][NPR]);
extern int field_Bhat_fluxrecon(FTYPE pr[][N2M][N3M][NPR], FTYPE pointfield[][N2M][N3M][NPR], FTYPE quasifield[][N2M][N3M][NPR]);
extern int ucons2upointppoint(SFTYPE boundtime, FTYPE pfield[][N2M][N3M][NPR],FTYPE unew[][N2M][N3M][NPR],FTYPE ulast[][N2M][N3M][NPR],FTYPE pcent[][N2M][N3M][NPR]);

extern int deaverage_ustag2pstag(FTYPE preal[][N2M][N3M][NPR], FTYPE ustag[][N2M][N3M][NPR], FTYPE pstag[][N2M][N3M][NPR]);


// functions for loop stuff
extern void  setup_nprlocalist(int whichprimtype, int *nprlocalstart, int *nprlocalend,int *nprlocallist, int *numprims);

/////////////////////////////////////
//
// NR STUFF
//
/////////////////////////////////////
extern int ludcmp(FTYPE **a, int n, int *indx, FTYPE *d);
extern void lubksb(FTYPE **a, int n, int *indx, FTYPE *d);
//extern FTYPE zbrent(FTYPE (*func) (FTYPE), FTYPE v1, FTYPE v2,
//		     FTYPE tol);


/* NR routines from nrutil.h */
extern int *ivector(long nl, long nh);
extern void free_ivector(int *v, long nl, long nh);
extern FTYPE *dvector(long nl, long nh);
extern void free_dvector(FTYPE *v, long nl, long nh);
extern FTYPE **dmatrix(long nrl, long nrh, long ncl, long nch);
extern void free_dmatrix(FTYPE **m, long nrl, long nrh, long ncl,
			 long nch);
extern FTYPE ***dtensor(long nrl, long nrh, long ncl, long nch,
			 long ndl, long ndh);
extern void free_dtensor(FTYPE ***t, long nrl, long nrh, long ncl,
			 long nch, long ndl, long ndh);
extern void nrerror(char error_text[]);

//////////////////////////////
//
// specialty functions
//
//////////////////////////////
extern void bondi_solve(FTYPE K, FTYPE gam, FTYPE *Rs, FTYPE *Urs,
			FTYPE *Edot);
extern FTYPE bondi_trace(FTYPE K, FTYPE gam, FTYPE edotf, FTYPE r,
			  FTYPE rs, FTYPE urs);
extern void timestep(FTYPE ndtr, FTYPE ndth);
extern FTYPE dtset(FTYPE ndtr, FTYPE ndth);

extern FTYPE bondi_trace(FTYPE K, FTYPE gam, FTYPE edotf,
			  FTYPE r, FTYPE rs, FTYPE urs);
extern void bondi_solve(FTYPE K, FTYPE gam, FTYPE *Rs,
			FTYPE *Urs, FTYPE *Edot);
extern FTYPE edot_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE dedr_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE dedur_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE d2edr2_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE d2edur2_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE d2edrdur_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);

extern void lower_vec(FTYPE *a, struct of_geom *geom, FTYPE *b);
extern void lowerf(FTYPE *a, struct of_geom *geom, FTYPE *b);
extern void raise_vec(FTYPE *v1, struct of_geom *geom, FTYPE *v2);
extern int gaussj(FTYPE **tmp, int n, FTYPE **b, int m);
extern void set_points(void);
// extern FTYPE delta(int j, int k) ;
// extern FTYPE mink(int j, int k) ;
extern void make_tetr(FTYPE *ucon, FTYPE econ[][NDIM]);


extern FTYPE sign_bad(FTYPE a);
extern FTYPE sign_func(FTYPE a);
#if(SUPERLONGDOUBLE)
#define sign(a) (sign_bad(a))
#else
#define sign(a) (copysign(1.0,a))
#endif


extern FTYPE signavoidzero(FTYPE a);

#ifndef WIN32
extern FTYPE max(FTYPE a, FTYPE b);

extern FTYPE min(FTYPE a, FTYPE b);
#endif
///////////////////////////////////
//
// SUPERLONGDOUBLE declarations
//
///////////////////////////////////

#if(SUPERLONGDOUBLE)
#include "mconf.h"
extern long double fabsl ( long double );
extern long double sqrtl ( long double );
extern long double cbrtl ( long double );
extern long double expl ( long double );
extern long double logl ( long double );
extern long double tanl ( long double );
extern long double atanl ( long double );
extern long double sinl ( long double );
extern long double asinl ( long double );
extern long double cosl ( long double );
extern long double acosl ( long double );
extern long double powl ( long double, long double );
extern long double tanhl ( long double );
extern long double atanhl ( long double );
extern long double sinhl ( long double );
extern long double asinhl ( long double );
extern long double coshl ( long double );
extern long double acoshl ( long double );
extern long double exp2l ( long double );
extern long double log2l ( long double );
extern long double exp10l ( long double );
extern long double log10l ( long double );
extern long double gammal ( long double );
extern long double lgaml ( long double );
extern long double jnl ( int, long double );
extern long double ynl ( int, long double );
extern long double ndtrl ( long double );
extern long double ndtril ( long double );
extern long double stdtrl ( int, long double );
extern long double stdtril ( int, long double );
extern long double ellpel ( long double );
extern long double ellpkl ( long double );
long double lgammal(long double);
extern int isfinitel ( long double );
#define finite(arg) isfinitel(arg)
#define copysign( a, b ) ( fabsl(a) * sign(b) ) 
extern int merror;
#else



#include <math.h>

#ifdef WIN32
#define finite(arg) _finite(arg)
#define isfinite(arg) _finite(arg)
#endif

#ifndef WIN32
//#if USINGICC==0
#if( !defined(isfinite))
// needed for Sauron
#define isfinite(arg) finite(arg)  //atch -- on mako, in force-free it would complain about multiply-defined __finite() if not include this line
#endif
#endif // end if not defined WIN32

#endif


#ifdef WIN32
#define copysign( a, b ) ( fabs(a) * sign(b) ) 
#endif




#if(!DO_ASSERTS)
#define assert assert_func_empty
#else
#define assert assert_func
#endif


extern int assert_func( int is_bad_val, char *s, ... );
extern int assert_func_empty( int is_bad_val, char *s, ... );


