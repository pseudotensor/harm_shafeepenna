
#include "decs.h"

/* diagnostics subroutine */

#define DIAGREPORT {trifprintf("t=%21.15g to do: tener=%21.15g (dt=%21.15g): dump_cnt=%ld @ t=%21.15g (dt=%21.15g) : avg_cnt=%ld @ t=%21.15g (dt=%21.15g) : debug_cnt=%ld @ t=%21.15g (dt=%21.15g) : image_cnt=%ld @ t=%21.15g (dt=%21.15g): restart=%d @ nstep=%ld (dt=%ld)\n",t,tdumpgen[DTENER],DTdumpgen[DTENER],dumpcntgen[DTDUMP],tdumpgen[DTDUMP],DTdumpgen[DTDUMP],dumpcntgen[DTAVG],tdumpgen[DTAVG],DTdumpgen[DTAVG],dumpcntgen[DTDEBUG],tdumpgen[DTDEBUG],DTdumpgen[DTDEBUG],dumpcntgen[DTIMAGE],tdumpgen[DTIMAGE],DTdumpgen[DTIMAGE],whichrestart,nrestart,DTr);}


static int get_dodumps(int call_code, int firsttime, SFTYPE localt, long localnstep, long localrealnstep, FTYPE *tdumpgen, FTYPE *tlastgen, FTYPE tlastareamap, long long int nlastrestart, long long int nrestart, int *dogdump, int *dordump, int *doareamap, int *dodumpgen);


int diag(int call_code, FTYPE localt, long localnstep, long localrealnstep)
{
  //////////////////
  //
  // BEGIN DIAG VARS
  //
  static int firsttime = 1;
  int pl;
  FTYPE asym[NPR], norm[NPR], maxasym[NPR];
  
  FILE *dumpcnt_filegen[NUMDTDS];
  static SFTYPE tlastgen[NUMDTDS];
  static SFTYPE tlastareamap;
  static long long int dumpcgen[NUMDTDS];
  static long long int restartc;
  static SFTYPE tdumpgen[NUMDTDS];
  static long long int nlastrestart,nrestart;
  int dogdump,dordump,doareamap;
  //  int doavg, dordump,dodump,doener,doimagedump,doareamap,dodebug,dofieldlinedump;
  int dodumpgen[NUMDTDS];
  int dtloop;

  int floor,i,j,k;
  int dir,interpi,enodebugi;
  int dissloop;
  int enodebugdump(long dump_cnt);
  int asym_compute_1(FTYPE (*prim)[N2M][N3M][NPR]);
  int whichDT;





  ///////////////////////////
  //
  // setup timing of writes to files
  //

  if ((call_code == INIT_OUT) || (firsttime == 1)) {
    
    // no need to repeat if restarting
    nlastrestart = (long) (DTr*(SFTYPE)((localrealnstep/DTr))-1);
    tlastareamap = localt-SMALL;

    for(dtloop=0;dtloop<NUMDTDS;dtloop++){
      tlastgen[dtloop] = DTdumpgen[dtloop]*(SFTYPE)((long long int)(t/DTdumpgen[dtloop]))-SMALL;
    }

    //    tlastener  = DTener*(SFTYPE)((int)(t/DTener))-SMALL;
    //tlastdump  = DTd*(SFTYPE)((int)(t/DTd))-SMALL;
    //tlastavg  = DTavg*(SFTYPE)((int)(t/DTavg))-SMALL;
    //tlastimage = DTi*(SFTYPE)((int)(t/DTi))-SMALL;
    //tlastdebug  = DTdebug*(SFTYPE)((int)(t/DTdebug))-SMALL;

    for(dtloop=0;dtloop<NUMDTDS;dtloop++){
      dumpcgen[dtloop]=0;
    }
    //    dumpc = avgc = imagec = restartc = enerc = debugc = 0;

    if (RESTARTMODE == 0) {

      for(dtloop=0;dtloop<NUMDTDS;dtloop++){
	tdumpgen[dtloop]=localt;
      }


      nrestart = localnstep;
      // override defaults:
      tdumpgen[DTAVG]=localt+DTdumpgen[DTAVG]; // do next time

      //      tdump = timage = tener = t;
      //      tavg=t+DTavg; // do next time

      for(dtloop=0;dtloop<NUMDTDS;dtloop++){
	dumpcntgen[dtloop]=0;
      }

      //      dump_cnt = 0;
      //      image_cnt = 0;
      //      rdump_cnt = 0;
      //      avg_cnt = 0;
      //      debug_cnt = 0;

      appendold = 0;
    } else {
      setrestart(&appendold);

      // assuming started at t=0 and localnstep=0 for original run
      // time to dump NEXT
      nrestart = (long)(DTr*(SFTYPE)((localrealnstep/DTr)+1));

      for(dtloop=0;dtloop<NUMDTDS;dtloop++){
	tdumpgen[dtloop]=DTdumpgen[dtloop]*(SFTYPE)((long long int)(localt/DTdumpgen[dtloop])+1);
      }

      DIAGREPORT;
    }
  }// end if firsttime or diag(INIT_OUT) called




  
  ///////////////////////
  //
  // setup what we will dump this call to diag()
  //
  get_dodumps(call_code, firsttime, localt, localnstep, localrealnstep, tdumpgen, tlastgen, tlastareamap, nlastrestart, nrestart, &dogdump, &dordump, &doareamap, dodumpgen);



  /// check if only wanted to know if next time would lead to making a dump file
  // used so only compute expensive things during step for diagnostics if outputting diagnostics
  if(call_code==FUTURE_OUT){
    if(dodumpgen[DTDUMP]) return(DOINGFUTUREOUT);
    else return(0);
  }






  //////////////////////////////////////////////////////////////////////
  //
  ////////////////// NOW WRITE TO FILES and update tdumpgen and tlastgen
  //
  // For non-blocking ROMIO to be effective as written, need to have largest dumps last
  //
  //////////////////////////////////////////////////////////////////////


#if(ASYMDIAGCHECK)
  asym_compute_1(pdump);
#endif







  // extra bounding for diagnostics
  if(
     // if doing simulbccalc type loop in step_ch.c then need to bound when doing diagnostics since not done yet
     (SIMULBCCALC>=1 && (dodumpgen[DTDUMP]||dordump||dodumpgen[DTENER]))
     // assume if PRODUCTION==1 then user knows divB won't be computed at MPI boundaries and that's ok.  Don't want to compute so avoid bounding that slows things down.  Assume ok to still compute dump file version, just not ener version that may be too often
     || (PRODUCTION==0 && (dodumpgen[DTDUMP]||dodumpgen[DTENER]))
     || (PRODUCTION==1 && (dodumpgen[DTDUMP]))
     ){

    // for dump, rdump, and divb in ener
    bound_allprim(STAGEM1,localt,pdump);
    if(DOENOFLUX != NOENOFLUX){
      // bound udump (unew) so divb can be computed at MPI boundaries (still real boundaries won't be computed correctly for OUTFLOW types)
      // Notice only need to bound 1 cell layer (BOUNDPRIMSIMPLETYPE) since divb computation only will need 1 extra cell
      bound_mpi(STAGEM1,BOUNDPRIMSIMPLETYPE,udump,NULL,NULL,NULL);
    }
  }
  



#if(PRODUCTION==0)
  // DEBUG:
  // GODMARK: need to disable to see if hurts internal solution since bounding this may help algorithm itself and we don't want that
  // so far don't see problem with stag + FV or FLUXRECON

  // bound conserved quantities so divb is computed correctly when doing higher-order methods
  // note that divb will always be in error at boundary since uses unew
  // don't need ustag in boundary cells, but do this to check divb on boundary when testing

  // not generally a good idea
  //  bound_uavg(STAGEM1,localt,udump);
#endif

  //////////////////////
  //
  // ener dump (integrated quantities: integrate and possibly  dump them too (if doener==1))
  //
  /////////////////////

  if(dodumpgen[DTENER]||dordump||(call_code==FINAL_OUT)||(call_code==INIT_OUT)){ // need integratd quantities for restart dump, but don't write them to ener file.

    dump_ener(dodumpgen[DTENER],dordump,call_code);
    if(COMPUTEFRDOT){
      if(dodumpgen[DTENER]||(call_code==FINAL_OUT)||(call_code==INIT_OUT)){
	frdotout(); // need to include all terms and theta fluxes on horizon/outer edge at some point GODMARK
      }
    }
    whichDT=DTENER;
    // below is really floor to nearest integer plus 1
    dumpcgen[whichDT] = 1 + MAX(0,(long long int)((localt-tdumpgen[whichDT])/DTdumpgen[whichDT]));
    tdumpgen[whichDT] = (ROUND2LONGLONGINT(tdumpgen[whichDT]/DTdumpgen[whichDT]) + dumpcgen[whichDT])*DTdumpgen[whichDT];
    tlastgen[DTENER]=localt;
  }




  ///////////////////////
  //
  // RESTART DUMP
  //
  ///////////////////////
  

  if(dordump){
    DIAGREPORT;
    trifprintf("dumping: restart: %d localnstep: %ld nlastrestart: %ld nrestart: %ld restartc: %d\n", whichrestart,localrealnstep,nlastrestart,nrestart,restartc);

    restart_write((long)whichrestart);	// 0 1 0 1 0 1 ...
    if(DOEVOLVEMETRIC) restartmetric_write((long)whichrestart);

    restartsteps[whichrestart] = localrealnstep;
    whichrestart = !whichrestart;

    restartc = 1 + MAX(0,(long long int)(((FTYPE)localrealnstep-(FTYPE)nrestart)/((FTYPE)DTr)));
    nrestart = (ROUND2LONGLONGINT((FTYPE)localrealnstep/((FTYPE)DTr)) + restartc) * DTr;
    nlastrestart=localrealnstep;
  }
      
  



  ///////////////////////
  //
  // DEBUG DUMP
  //
  ///////////////////////
  
  if(dodumpgen[DTDEBUG]){
    DIAGREPORT;
    trifprintf("debug dumping: debug_cnt=%ld : t=%21.15g tlastdebug=%21.15g tdebug=%21.15g debugc=%d\n", dumpcntgen[DTDEBUG],localt,tlastgen[DTDEBUG],tdumpgen[DTDEBUG],dumpcgen[DTDEBUG]);
    
    /* make regular dump file */
    if (debugdump(dumpcntgen[DTDEBUG]) >= 1){
      dualfprintf(fail_file,"unable to print debug dump file\n");
      return (1);
    }

    if(DOENODEBUG){
      /* make regular dump file */
      if (enodebugdump(dumpcntgen[DTDEBUG]) >= 1){
	dualfprintf(fail_file,"unable to print enodebug dump file\n");
	return (1);
      }
    }

    // iterate counter
    dumpcntgen[DTDEBUG]++;
    whichDT=DTDEBUG;
    // below is really floor to nearest integer plus 1
    dumpcgen[whichDT] = 1 + MAX(0,(long long int)((localt-tdumpgen[whichDT])/DTdumpgen[whichDT]));
    tdumpgen[whichDT] = (ROUND2LONGLONGINT(tdumpgen[whichDT]/DTdumpgen[whichDT]) + dumpcgen[whichDT])*DTdumpgen[whichDT];
    // output number of dumps
    myfopen("dumps/0_numdebug.dat","w","error opening debug dump count file\n",&dumpcnt_filegen[DTDEBUG]);      
    myfprintf(dumpcnt_filegen[DTDEBUG], "# Number of debug dumps\n%ld\n", dumpcntgen[DTDEBUG]);
    myfclose(&dumpcnt_filegen[DTDEBUG],"Couldn't close debugcnt_file");
    tlastgen[DTDEBUG]=localt;
  }



  ///////////////////////
  //
  // AREA MAP
  //
  ///////////////////////
  
  if(doareamap){
    if(area_map(call_code, TIMESERIESAREAMAP, 20, ifail, jfail, kfail, pdump)>=1) return(1);
    tlastareamap=t;
  }

  ///////////////////////
  //
  // IMAGE
  //
  ///////////////////////

  if(dodumpgen[DTIMAGE]){
    DIAGREPORT;
    trifprintf("image dump %ld : t=%21.15g tlastimage=%21.15g timage=%21.15g imagec=%d\n", dumpcntgen[DTIMAGE], localt,tlastgen[DTIMAGE],tdumpgen[DTIMAGE],dumpcgen[DTIMAGE]);
    
    
    /* make regular image file */
    if(image_dump(dumpcntgen[DTIMAGE])>=1) return(1);

    // iterate counter
    dumpcntgen[DTIMAGE]++;

    whichDT=DTIMAGE;
    // below is really floor to nearest integer plus 1
    dumpcgen[whichDT] = 1 + MAX(0,(long long int)((localt-tdumpgen[whichDT])/DTdumpgen[whichDT]));
    tdumpgen[whichDT] = (ROUND2LONGLONGINT(tdumpgen[whichDT]/DTdumpgen[whichDT]) + dumpcgen[whichDT])*DTdumpgen[whichDT];
    // output number of images
    myfopen("images/0_numimages.dat","w","error opening image count file\n",&dumpcnt_filegen[DTIMAGE]);      
    myfprintf(dumpcnt_filegen[DTIMAGE], "# Number of images\n%ld\n", dumpcntgen[DTIMAGE]);
    myfclose(&dumpcnt_filegen[DTIMAGE],"Couldn't close imagecnt_file");
    tlastgen[DTIMAGE]=localt;
  }


  ///////////////////////
  //
  // FIELDLINE
  //
  ///////////////////////

  if(dodumpgen[DTFIELDLINE]){
    DIAGREPORT;
    trifprintf("fieldline dump %ld : t=%21.15g tlastfieldline=%21.15g tfieldline=%21.15g fieldlinec=%d\n", dumpcntgen[DTFIELDLINE], localt,tlastgen[DTFIELDLINE],tdumpgen[DTFIELDLINE],dumpcgen[DTFIELDLINE]);
    
    // (after processing) equivalent to image in interest
    if(fieldlinedump(dumpcntgen[DTFIELDLINE])>=1) return(1);

    // iterate counter
    dumpcntgen[DTFIELDLINE]++;

    whichDT=DTFIELDLINE;
    // below is really floor to nearest integer plus 1
    dumpcgen[whichDT] = 1 + MAX(0,(long long int)((localt-tdumpgen[whichDT])/DTdumpgen[whichDT]));
    tdumpgen[whichDT] = (ROUND2LONGLONGINT(tdumpgen[whichDT]/DTdumpgen[whichDT]) + dumpcgen[whichDT])*DTdumpgen[whichDT];
    // output number of fieldlines
    myfopen("dumps/0_numfieldlines.dat","w","error opening fieldline count file\n",&dumpcnt_filegen[DTFIELDLINE]);
    myfprintf(dumpcnt_filegen[DTFIELDLINE], "# Number of fieldlines\n%ld\n", dumpcntgen[DTFIELDLINE]);
    myfclose(&dumpcnt_filegen[DTFIELDLINE],"Couldn't close fieldlinecnt_file");
    tlastgen[DTFIELDLINE]=localt;
  }






  ///////////////////////
  //
  // DUMP
  //
  ///////////////////////
  
  if(dodumpgen[DTDUMP]){
    DIAGREPORT;
    trifprintf("dumping: dump_cnt=%ld : t=%21.15g tlastdump=%21.15g tdump=%21.15g dumpc=%d\n", dumpcntgen[DTDUMP],localt,tlastgen[DTDUMP],tdumpgen[DTDUMP],dumpcgen[DTDUMP]);


    if(dnumcolumns[EOSDUMPCOL]>0){ // otherwise no point
      if (eosdump(dumpcntgen[DTDUMP]) >= 1){
	dualfprintf(fail_file,"unable to print eosdump file\n");
	return (1);
      }
    }

    if(DOVPOTDUMP){
      if (vpotdump(dumpcntgen[DTDUMP]) >= 1){
	dualfprintf(fail_file,"unable to print vpotdump file\n");
	return (1);
      }
    }


    // so can restart at a dump without reconstructing the rdump from a dump.
    // Also, if run out of disk space then normal rdump's can be corrupted
    restart_write(-(long)dumpcntgen[DTDUMP]-1);
    if(DOEVOLVEMETRIC) restartmetric_write(-(long)dumpcntgen[DTDUMP]-1);


    if(FLUXDUMP){
      if (fluxdumpdump(dumpcntgen[DTDUMP]) >= 1){
	dualfprintf(fail_file,"unable to print fluxdump file\n");
	return (1);
      }
    }

    if(DODUMPOTHER){
      if (dumpother(dumpcntgen[DTDUMP]) >= 1){
	dualfprintf(fail_file,"unable to print dumpother file\n");
	return (1);
      }
    }
    


    if(DODISS){
      /* make dissdump file */
      if (dissdump(dumpcntgen[DTDUMP]) >= 1){
	dualfprintf(fail_file,"unable to print dissdump file\n");
	return (1);
      }
    }


    /* make regular dump file */
    if (dump(dumpcntgen[DTDUMP]) >= 1){
      dualfprintf(fail_file,"unable to print dump file\n");
      return (1);
    }


    if(DOEVOLVEMETRIC&&(N2==1)&&(N3==1)&&(DOGDUMPDIAG)&&(!GAMMIEDUMP)&&((RESTARTMODE==0))){
      // only reasonable to do if 1-D
      gdump(dumpcntgen[DTDUMP]);
    }

    // iterate counter
    dumpcntgen[DTDUMP]++;
    whichDT=DTDUMP;
    // below is really floor to nearest integer plus 1
    dumpcgen[whichDT] = 1 + MAX(0,(long long int)((localt-tdumpgen[whichDT])/DTdumpgen[whichDT]));
    tdumpgen[whichDT] = (ROUND2LONGLONGINT(tdumpgen[whichDT]/DTdumpgen[whichDT]) + dumpcgen[whichDT])*DTdumpgen[whichDT];
    // output number of dumps
    myfopen("dumps/0_numdumps.dat","w","error opening dump count file\n",&dumpcnt_filegen[DTDUMP]);      
    myfprintf(dumpcnt_filegen[DTDUMP], "# Number of dumps\n%ld\n", dumpcntgen[DTDUMP]);
    myfclose(&dumpcnt_filegen[DTDUMP],"Couldn't close dumpcnt_file");
    tlastgen[DTDUMP]=localt;
  }




  ///////////////////////
  //
  // AVG
  //
  ///////////////////////

  if(DOAVGDIAG){
    // do every time step
    // assume can't fail, but can apparently
    if(average_calc(dodumpgen[DTAVG])>=1) return(1);
  }

  
  if(dodumpgen[DTAVG]){
    DIAGREPORT
    trifprintf("avging dump: avg_cnt=%ld : t=%21.15g tlastavg=%21.15g tavg=%21.15g avgc=%d\n", dumpcntgen[DTAVG],localt,tlastgen[DTAVG],tdumpgen[DTAVG],dumpcgen[DTAVG]);
    
    /* make avg dump file */
    if (avgdump(dumpcntgen[DTAVG]) >= 1){
      dualfprintf(fail_file,"unable to print avg file\n");
      return (1);
    }
#if(DOAVG2)
    /* make avg dump file */
    if (avgdump2(dumpcntgen[DTAVG]) >= 1){
      dualfprintf(fail_file,"unable to print avg2 file\n");
      return (1);
    }
#endif

    // iterate counter
    dumpcntgen[DTAVG]++;

    whichDT=DTAVG;
    // below is really floor to nearest integer plus 1
    dumpcgen[whichDT] = 1 + MAX(0,(long long int)((localt-tdumpgen[whichDT])/DTdumpgen[whichDT]));
    tdumpgen[whichDT] = (ROUND2LONGLONGINT(tdumpgen[whichDT]/DTdumpgen[whichDT]) + dumpcgen[whichDT])*DTdumpgen[whichDT];
    // output number of avgs
    myfopen("dumps/0_numavgs.dat","w","error opening avg count file\n",&dumpcnt_filegen[DTAVG]);      
    myfprintf(dumpcnt_filegen[DTAVG], "# Number of avgs\n%ld\n", dumpcntgen[DTAVG]);
    myfclose(&dumpcnt_filegen[DTAVG],"Couldn't close avgcnt_file");
    tlastgen[DTAVG]=localt;
  }


  //////////////////////
  //
  // Grid dump
  //
  /////////////////////

  if(dogdump) gdump(-1); // -1 means no file number on filename











  ///////////////////////
  //
  // DIAG clearings
  //
  // finally some clearing based upon what called
  //
  ////////////////////////

  if(DODEBUG){ // shouldn't clear these till after ener and debug dump done so both have all timescale data.
    if(dodumpgen[DTENER]){
      // cleanse the ener time scale for the failure diag
      ZLOOP FLOORLOOP(floor) failfloorcount[i][j][k][ENERTS][floor]=0;
    }
    if(dodumpgen[DTDEBUG]){
      // clense failure diag
      ZLOOP FLOORLOOP(floor) failfloorcount[i][j][k][DEBUGTS][floor]=0;
    }
    if(dodumpgen[DTIMAGE]){
      // clense the failure counts for this time scale
      ZLOOP FLOORLOOP(floor) failfloorcount[i][j][k][IMAGETS][floor]=0;
    }
  }

  if(DOENODEBUG){
    if(dodumpgen[DTDEBUG]){
      FULLLOOP DIMENLOOP(dir) INTERPLOOP(interpi) PDIAGLOOP(pl) ENODEBUGLOOP(enodebugi){
	if(dir<=2 && pl<=U2){
	  enodebugarray[i][j][k][dir-1][interpi][pl][enodebugi]=0;
	}
      }
    }
  }

  if(0&&DODISS){
    // clear diss (I don't see why one should do this and tie results to dump frequency)
    if(dodumpgen[DTDUMP]){
      FULLLOOP{
	for(dissloop=0;dissloop<NUMDISSVERSIONS;dissloop++){
	  // clear failures too
	  dissfunpos[i][j][k][dissloop]=0.0;
	}
      }
    }
  }


/*	//symmetry check
	j = 0;
	k = 0;
	PDIAGLOOP(pl) maxasym[pl] = 0.0;

	PDIAGLOOP(pl) {
		if( pl >= U1 && pl <= U3 ) {
			norm[pl] = coordparams.timescalefactor;
		}
		else if( pl == UU ) {
			norm[pl] = coordparams.timescalefactor * coordparams.timescalefactor;
		}
		else {
			norm[pl] = 1.;
		}

		norm[pl] = 1 / norm[pl];
	}
	
	for( i = 0; i < N1; i++ ) {
		PDIAGLOOP( pl ) {
			if( pl == U1 ) {
				asym[pl] = fabs( p[i][j][k][pl] + p[N1 - 1 - i][j][k][pl] );
			}
			else {
				asym[pl] = fabs( p[i][j][k][pl] - p[N1 - 1 - i][j][k][pl] );
			}


			asym[pl] /= norm[pl];

			maxasym[pl] = MAX( asym[pl], maxasym[pl] );
		}
	}

	dualfprintf( fail_file, "asym %ld ", localrealnstep );
		
	PDIAGLOOP(pl) if( pl <= U1 ) dualfprintf( fail_file, " %22.16g", maxasym[pl] );

	dualfprintf( fail_file, "\n" );
*/
  firsttime = 0;
  return (0);
}







int get_dodumps(int call_code, int firsttime, SFTYPE localt, long localnstep, long localrealnstep, FTYPE *tdumpgen, FTYPE *tlastgen, FTYPE tlastareamap, long long int nlastrestart, long long int nrestart, int *dogdump, int *dordump, int *doareamap, int *dodumpgen)
{

  // output grid (probaly want both fullgrid (to make sure ok) and compute grid to compare with data dumps
  if((DOGDUMPDIAG)&&(!GAMMIEDUMP)&&(firsttime&&(RESTARTMODE==0))){
    *dogdump=1;
  }
  else *dogdump=0;

  if((DORDUMPDIAG)&&( ((nlastrestart!=nrestart)&&(failed == 0) && (localrealnstep >= nrestart ))||(call_code==FINAL_OUT) ) ){
    *dordump=1;
  }
  else *dordump=0;


  if((DOAREAMAPDIAG)&&(localt != tlastareamap)&&(dofailmap)&&(localnstep>=steptofailmap)){
    *doareamap=1;
  }
  else *doareamap=0;


  ////////////////////////////
  //
  // rest are on NUMDTDS list
  //
  ////////////////////////////


  if( ((DODUMPDIAG)&&(DODIAGEVERYSUBSTEP||((localt!=tlastgen[DTDUMP])&&(localt >= tdumpgen[DTDUMP] || (RESTARTMODE&&dofaildump&&(localnstep>=steptofaildump)) || call_code==FINAL_OUT ))) )  ){
    dodumpgen[DTDUMP]=1;
  }
  else dodumpgen[DTDUMP]=0;

  if( ((DODEBUG)&&(DOENODEBUGEVERYSUBSTEP||DODIAGEVERYSUBSTEP||((localt!=tlastgen[DTDEBUG])&&(localt >= tdumpgen[DTDEBUG] || call_code==FINAL_OUT))) )  ){
    dodumpgen[DTDEBUG]=1;
  }
  else dodumpgen[DTDEBUG]=0;

  
  if((DOAVGDIAG)&&((localt!=tlastgen[DTAVG])&&(localt >= tdumpgen[DTAVG] || call_code==FINAL_OUT))){
    dodumpgen[DTAVG]=1;
  }
  else dodumpgen[DTAVG]=0;
    

  // t!=tlast avoids duplicate entries
  if(DOENERDIAG&&(DODIAGEVERYSUBSTEP||((localt!=tlastgen[DTENER])&&( (localt >= tdumpgen[DTENER])||(call_code==INIT_OUT)||(call_code==FINAL_OUT)||firsttime)))){
    dodumpgen[DTENER]=1;
  }
  else dodumpgen[DTENER]=0;


  /* image dump at regular intervals */
  if((DOIMAGEDIAG)&&(DODIAGEVERYSUBSTEP||((localt!=tlastgen[DTIMAGE])&&(localt >= tdumpgen[DTIMAGE] || call_code==FINAL_OUT))) ){
    dodumpgen[DTIMAGE]=1;
  }
  else dodumpgen[DTIMAGE]=0;


  /* fieldline dump at regular intervals */
  // use image time period
  if((DOFIELDLINE)&&(DODIAGEVERYSUBSTEP||((localt!=tlastgen[DTFIELDLINE])&&(localt >= tdumpgen[DTFIELDLINE] || call_code==FINAL_OUT))) ){
    dodumpgen[DTFIELDLINE]=1;
    dumpcntgen[DTFIELDLINE]=dumpcntgen[DTIMAGE]; // force to go with image dump (needed also so restart() knows the count)
  }
  else dodumpgen[DTFIELDLINE]=0;




  return(0);
}












/** some diagnostic routines **/





int asym_compute_1(FTYPE (*prim)[N2M][N3M][NPR])
{
  int i,j,k;
  int pl;

#if(0)
// for implosion problem
  FULLLOOP{

    if(prim[i][j][k][RHO]!=prim[j][i][k][RHO]){
      dualfprintf(fail_file,"ASYM in RHO %d %d %d : %23.16g %23.16g\n",i,j,k,prim[i][j][k][RHO],prim[j][i][k][RHO]);
    }

    if(prim[i][j][k][UU]!=prim[j][i][k][UU]){
      dualfprintf(fail_file,"ASYM in UU %d %d %d : %23.16g %23.16g\n",i,j,k,prim[i][j][k][UU],prim[j][i][k][UU]);
    }


    if(prim[i][j][k][U1]!=prim[j][i][k][U2]){
      dualfprintf(fail_file,"ASYM in U1 %d %d %d : %23.16g %23.16g\n",i,j,k,prim[i][j][k][U1],prim[j][i][k][U2]);
    }
     
  }
#endif
#if(1)
  // for any periodic boundary conditions
  FULLLOOP{
    if(i<N1/2 && j<N2/2 && k<N3/2){
      PLOOP(pl){
	if(prim[i][j][k][pl]!=prim[i+N1][j][k][pl]){
	  dualfprintf(fail_file,"ASYM nstep=%ld steppart=%d in pl=%d dir=1 :: %d %d %d : %23.16g %23.16g\n",nstep,steppart,pl,i,j,k,prim[i][j][k][pl],prim[i+N1][j][k][pl]);
	}
	if(prim[i][j][k][pl]!=prim[i][j+N2][k][pl]){
	  dualfprintf(fail_file,"ASYM nstep=%ld steppart=%d in pl=%d dir=2 :: %d %d %d : %23.16g %23.16g\n",nstep,steppart,pl,i,j,k,prim[i][j][k][pl],prim[i][j+N2][k][pl]);
	}
      }
    }
  }

  // for any periodic boundary conditions
  FULLLOOP{
    if(i<N1/2 && j<N2/2 && k<N3/2){
      PLOOP(pl){
	if(unew[i][j][k][pl]!=unew[i+N1][j][k][pl]){
	  dualfprintf(fail_file,"ASYMUNEW nstep=%ld steppart=%d in pl=%d dir=1 :: %d %d %d : %23.16g %23.16g\n",nstep,steppart,pl,i,j,k,unew[i][j][k][pl],unew[i+N1][j][k][pl]);
	}
	if(unew[i][j][k][pl]!=unew[i][j+N2][k][pl]){
	  dualfprintf(fail_file,"ASYMUNEW nstep=%ld steppart=%d in pl=%d dir=2 :: %d %d %d : %23.16g %23.16g\n",nstep,steppart,pl,i,j,k,unew[i][j][k][pl],unew[i][j+N2][k][pl]);
	}
      }
    }
  }

#endif


  return(0);
}

// for implosion problem
int asym_compute_2(FTYPE (*prim)[N2M][N3M][NPR])
{
  int i,j,k;
  int pl;


#if(0)
  ZLOOP{

    if(prim[i][j][k][RHO]!=prim[j][i][k][RHO]){
      dualfprintf(fail_file,"ASYM in RHO %d %d %d : %23.16g %23.16g\n",i,j,k,prim[i][j][k][RHO],prim[j][i][k][RHO]);
    }

    if(prim[i][j][k][UU]!=prim[j][i][k][UU]){
      dualfprintf(fail_file,"ASYM in UU %d %d %d : %23.16g %23.16g\n",i,j,k,prim[i][j][k][UU],prim[j][i][k][UU]);
    }


    if(prim[i][j][k][U1]!=prim[j][i][k][U2]){
      dualfprintf(fail_file,"ASYM in U1 %d %d %d : %23.16g %23.16g\n",i,j,k,prim[i][j][k][U1],prim[j][i][k][U2]);
    }
    
  }
#endif
#if(1)
  // for any periodic boundary conditions
  ZLOOP{
    if(i<N1/2 && j<N2/2 && k<N3/2){
      PLOOP(pl){
	if(prim[i][j][k][pl]!=prim[i+N1][j][k][pl]){
	  dualfprintf(fail_file,"ASYM nstep=%ld steppart=%d in pl=%d dir=1 :: %d %d %d : %23.16g %23.16g\n",nstep,steppart,pl,i,j,k,prim[i][j][k][pl],prim[i+N1][j][k][pl]);
	}
	if(prim[i][j][k][pl]!=prim[i][j+N2][k][pl]){
	  dualfprintf(fail_file,"ASYM nstep=%ld steppart=%d in pl=%d dir=2 :: %d %d %d : %23.16g %23.16g\n",nstep,steppart,pl,i,j,k,prim[i][j][k][pl],prim[i][j+N2][k][pl]);
	}
      }
    }
  }
#endif



  return(0);
}



// 2D only for now since really only useful for 2D imaging

/* map out region around failure point */
int area_map(int call_code, int type, int size, int i, int j, int k, FTYPE prim[][N2M][N3M][NPR])
{
  int pl;
  int l,m,ll,mm;
  FTYPE vmin1, vmax1, vmin2, vmax2;
  int ignorecourant;
  struct of_geom geom;
  struct of_state q;
  FTYPE X[NDIM],V[NDIM];
  FTYPE divb;
  FTYPE tens_em[NDIM][NDIM], tens_matt[NDIM][NDIM], b[NDIM],
    ucon[NDIM];
  FTYPE U[NPR];
  int lowersizex1,uppersizex1;
  int lowersizex2,uppersizex2;

  static FILE* fileptr;
  static int firsttime=1;
  static int domap=0;
  static int doclose=0;



  trifprintf("\nStart area_map function ... ");


  k=N3/2+SHIFT3; // middle cell // GODMARK : 2D only below


  if(i-(-N1BND)<size/2) lowersizex1=i-(-N1BND);
  else lowersizex1=size/2;
  if((N1-1+N1BND)-i<size/2) uppersizex1=(N1-1+N1BND)-i;
  else uppersizex1=size/2;
  
  if(j-(-N2BND)<size/2) lowersizex2=j-(-N2BND);
  else lowersizex2=size/2;
  if((N2-1+N2BND)-j<size/2) uppersizex2=(N2-1+N2BND)-j;
  else uppersizex2=size/2;


  if(firsttime){
    if((type==TIMESERIESAREAMAP)&&(dofailmap)){
      if((fileptr=fopen("areamap","wt"))==NULL){
	dualfprintf(fail_file,"Cannot open ./areamap on proc=%d\n",myid);
	domap=0;
      }
      else domap=1;
    }
  }

  if((type==TIMESERIESAREAMAP)&&domap&&(call_code==2)){
    doclose=1;
  }
  else doclose=0;

  if(type==FINALTDUMPAREAMAP){
    dualfprintf(fail_file, "area map\n");
    dualfprintf(fail_file, "failure at: i=%d j=%d k=%d\n",i+startpos[1],j+startpos[2],k+startpos[3]);
    coord(i,j,k,CENT,X);
    dualfprintf(fail_file, "failure at: i=%d j=%d k=%d\n",i+startpos[1],j+startpos[2],k+startpos[3]);



    PDIAGLOOP(pl) {// all vars
      
      dualfprintf(fail_file, "variable %d \n", pl);
      
      dualfprintf(fail_file, "i = \t ");
      for(l=i-lowersizex1;l<=i+uppersizex1;l++){
	ll=l+startpos[1];
	dualfprintf(fail_file, "%21d", ll);
      }
      dualfprintf(fail_file, "\n");
      for(m=j-lowersizex2;m<=j+uppersizex2;m++){
	mm=m+startpos[2];
	dualfprintf(fail_file, "j = %d \t ",mm);
	for(l=i-lowersizex1;l<=i+lowersizex1;l++){
	  ll=l+startpos[1];
	  dualfprintf(fail_file, "%21.15g ",prim[l][m][k][pl]);
	}
	dualfprintf(fail_file, "\n");
      }
    }
  }
  else if((type==TIMESERIESAREAMAP)&&(domap)){
    if(firsttime){
      // GODMARK: 2D only, not outputting z-related stuff so function remains consistent with SM macro
      fprintf(fileptr,"%21.15g %d %d %21.15g %21.15g %21.15g %21.15g %d %d %d %d %d %d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",
	      t,totalsize[1],totalsize[2],startx[1],startx[2],dx[1],dx[2],lowersizex1,uppersizex1,lowersizex2,uppersizex2,startpos[1]+i,startpos[2]+j,gam,a,R0,Rin,Rout,hslope);
      fflush(fileptr);
    }
    for(m=j-size/2;m<=j+size/2;m++){
      if((m<-N2BND)||(m>N2-1+N2BND)) continue;
      mm=m+startpos[2];
      for(l=i-size/2;l<=i+size/2;l++){
	if((l<-N1BND)||(l>N1-1+N1BND)) continue;
	ll=l+startpos[1];
	
	
	coord(l, m, k, CENT, X);
	bl_coord(X, V); 
	get_geometry(l, m, k, CENT, &geom);
	if (!failed) {
	  if (get_state(prim[l][m][k], &geom, &q) >= 1)
	    FAILSTATEMENT("diag.c:areamap()", "get_state() dir=0", 1);
	  if (vchar(prim[l][m][k], &q, 1, &geom, &vmax1, &vmin1,&ignorecourant) >= 1)
	    FAILSTATEMENT("diag.c:areamap()", "vchar() dir=1or2", 1);
	  if (vchar(prim[l][m][k], &q, 2, &geom, &vmax2, &vmin2,&ignorecourant) >= 1)
	    FAILSTATEMENT("diag.c:areamap()", "vchar() dir=1or2", 2);
	  // GODMARK: no 3-direction char.
	}
	if((l>=-1)&&(l<=N1+1)&&(m>=-1)&&(m<=N2+1)&&(k>=-1)&&(k<=N3+1) ){ setfdivb(&divb, prim, udump, l, m, k);} // udump here must be set externally GODMARK
	else divb=0.0;

	// same order as dump.c for first columns (easy sm read)
	fprintf(fileptr,
		"%d %d "
		"%21.15g %21.15g "
		"%21.15g %21.15g "
		"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g "
		"%21.15g "
		"%21.15g %21.15g %21.15g %21.15g "
		"%21.15g %21.15g %21.15g %21.15g "
		"%21.15g %21.15g %21.15g %21.15g "
		"%21.15g %21.15g %21.15g %21.15g "
		"%21.15g %21.15g %21.15g %21.15g "
		"%21.15g "
		"%21.15g %ld\n",
		ll,mm,
		X[1],X[2],
		V[1],V[2],
		prim[l][m][k][0],
		prim[l][m][k][1],
		prim[l][m][k][2],
		prim[l][m][k][3],
		prim[l][m][k][4],
		prim[l][m][k][5],
		prim[l][m][k][6],
		prim[l][m][k][7],
		divb,
		q.ucon[0],q.ucon[1],q.ucon[2],q.ucon[3],
		q.ucov[0],q.ucov[1],q.ucov[2],q.ucov[3],
		q.bcon[0],q.bcon[1],q.bcon[2],q.bcon[3],
		q.bcov[0],q.bcov[1],q.bcov[2],q.bcov[3],
		vmin1,vmax1,vmin2,vmax2,
		geom.g,
		t,realnstep);
      }
    }
    fflush(fileptr);
  }

  if(doclose) if(fileptr!=NULL) fclose(fileptr);


  /* print out other diagnostics here */

  firsttime=0;
  trifprintf("end area_map function.\n");  
  return(0);
}



/* evaluate fluxed based diagnostics; put results in global variables */

#define JETBSQORHO (3.162)


// notice that F1 and F2 have arbitrary eomfunc that must be divided out to get real flux
int diag_flux(FTYPE prim[][N2M][N3M][NPR], FTYPE F1[][N2M][N3M][NPR], FTYPE F2[][N2M][N3M][NPR], FTYPE F3[][N2M][N3M][NPR], SFTYPE Dt)
{
  int fluxdir;
  int i, j, k, pl, l, dir,fl,enerregion;
  SFTYPE surface,mysurface,surf2;
  SFTYPE surgdet;
  FTYPE (*flux)[N2M][N3M][NPR];
  int start1,start2,start3,stop1,stop2,stop3;
  int gpos;
  int ii;
  FTYPE ftemp;
  FTYPE ftemp0,ftemp1,ftemp2,ftemp3,ftemp4,ftemp5,ftemp6;
  FTYPE pgas,bsq,bsqorho;
  struct of_geom geom;
  struct of_state q;
  FTYPE X[NDIM],V[NDIM];
  FTYPE b[NDIM],ucon[NDIM];
  FTYPE U[NPR];
  int condjet2;
  FTYPE Ftemp[NPR],Ftempdiag[NPR];
  int firstinloop;
  FTYPE pr[NPR];
  int *localdoflux, *localenerpos;
  SFTYPE (*localpcum)[NPR];
  SFTYPE (*localpdot)[NPR];
  SFTYPE (*localpdotterms)[NUMFLUXTERMS][NPR];



  // initialize
  ENERREGIONLOOP(enerregion){
    localdoflux=dofluxreg[enerregion];
    localenerpos=enerposreg[enerregion];
    localpcum=pcumreg[enerregion];
    localpdot=pdotreg[enerregion];
    localpdotterms=pdottermsreg[enerregion];
    

    DIRLOOP(dir) PDIAGLOOP(pl){
      localpdot[dir][pl]=0;
      FLLOOP(fl) localpdotterms[dir][fl][pl]=0;
      if(enerregion==0) FLLOOP(fl) pdottermsjet2[dir][fl][pl]=0;
    }
    // true outer boundary surface fluxes (per direction, per conserved variable)
    DIRLOOP(dir) {
      if (localdoflux[dir] >= 0) {// then this cpu is involved in flux BOUNDARY calculation
	// otherwise don't add to it


	///////
	// assumes rectangular region
	//
	if((dir==X1UP)||(dir==X1DN)){

	  surface = dx[2]*dx[3];

	  start1=stop1=localdoflux[dir];

	  start2=localenerpos[X2DN];
	  stop2=localenerpos[X2UP];

	  start3=localenerpos[X3DN];
	  stop3=localenerpos[X3UP];

	  flux=F1;
	  gpos=FACE1;
	  fluxdir=1;
	}
	else if((dir==X2UP)||(dir==X2DN)){

	  surface = dx[1]*dx[3];

	  start2=stop2=localdoflux[dir];

	  start1=localenerpos[X1DN];
	  stop1=localenerpos[X1UP];

	  start3=localenerpos[X3DN];
	  stop3=localenerpos[X3UP];

	  flux=F2;
	  gpos=FACE2;
	  fluxdir=2;
	}
	else if((dir==X3UP)||(dir==X3DN)){

	  surface = dx[1]*dx[2];

	  start3=stop3=localdoflux[dir];

	  start1=localenerpos[X1DN];
	  stop1=localenerpos[X1UP];

	  start2=localenerpos[X2DN];
	  stop2=localenerpos[X2UP];

	  flux=F3;
	  gpos=FACE3;
	  fluxdir=3;
	}
	else{
	  dualfprintf(fail_file,"no such direction: %d\n",dir);
	  myexit(1);
	}

	// zero out summation quantities
	PDIAGLOOP(pl){
	  localpdot[dir][pl]=0;
	  FLLOOP(fl) localpdotterms[dir][fl][pl]=0;
	  if(enerregion==0) FLLOOP(fl) pdottermsjet2[dir][fl][pl]=0;
	}
	// GENLOOP set by enerregion so set correctly for GRID SECTION
	GENLOOP(i,j,k,start1,stop1,start2,stop2,start3,stop3){
	  // now add up all zones for each conserved quantity (PDIAGLOOP)



	  ////////////
	  //
	  // GET STATE
	  //
	  ////////////



	  ////////////
	  get_geometry(i, j, k, CENT, &geom);
	  PDIAGLOOP(pl) pr[pl]=prim[i][j][k][pl];

	  //coord(i, j, k, CENT, X);
	  //bl_coord(X, V);
	    // if failed, then data output for below invalid, but columns still must exist    
	    // CENT since p at center
	  if (!failed) {
	    if (get_state(pr, &geom, &q) >= 1)
	      FAILSTATEMENT("diag.c:diag_flux()", "get_state() dir=0", 1);
	  }
	  else {// do a per zone check, otherwise set to 0
	    whocalleducon=1; // force no failure mode, just return like failure, and don't return if failure, just set to 0 and continue
	    if (get_state(pr, &geom, &q) >= 1){
	      for (l = 0; l < NDIM; l++)
		q.ucon[l]=0;
	      for (l = 0; l < NDIM; l++)
		q.ucov[l]=0;
	      for (l = 0; l < NDIM; l++)
		q.bcon[l]=0;
	      for (l = 0; l < NDIM; l++)
		q.bcov[l]=0;
	    }
	    whocalleducon=0; // return to normal state
	  }
	  
	  // somewhat like function which averages stress terms
	  pgas = pressure_rho0_u(pr[RHO],pr[UU]);
	  bsq=0; for(l=0;l<NDIM;l++) bsq+=(q.bcon[l])*(q.bcov[l]);




	  // modify effective position of integral in case adding energy/angular momentum to horizon
	  // apparently not correct
	  if(0 && enerregion==OUTSIDEHORIZONENERREGION){
	    mysurface = surface /(q.ucov[TT]);
	  }
	  else mysurface = surface;




	  ////////////
	  //
	  // do standard flux for accounting
	  //
	  ////////////

	  get_geometry(i, j, k, gpos, &geom);
	  PDIAGLOOP(pl) Ftemp[pl]=flux[i][j][k][pl]*mysurface; // in UEVOLVE form
	  // GODMARK: for finite volume method, below doesn't change the result if eomfunc=gdet.
	  // Otherwise flux would have to be completely recomputed for gdet case JUST for diagnostic to be consistent at higher order
	  UtoU_evolve2diag(UEVOLVE,UDIAG,&geom,Ftemp,Ftempdiag); // convert to diag form
	  PDIAGLOOP(pl){
	    if(!isfinite(Ftempdiag[pl])){
	      //	      dualfprintf(fail_file,"pl=%d i=%d j=%d k=%d enerregion=%d\n",pl,i,j,k,enerregion);
	      Ftempdiag[pl]=0.0; // assume if nan that just box is beyond where flux defined
	    }
	    localpdot[dir][pl]  += Ftempdiag[pl];
	  }
	  if(REMOVERESTMASSFROMUU==2){
	    //	    localpdot[dir][UU] += -Ftempdiag[RHO]; // add rest-mass back in (GODMARK: Problem is that non-relativistic problems will have energy term swamped by mass term)
	  }
	  else{
	    // then already correct
	    // so fix this instead!
	    localpdot[dir][UU] += Ftempdiag[RHO]; // remove rest-mass term.
	  }

	  /*
	  // DEBUG
	  PDIAGLOOP(pl) if(!finite(localpdot[dir][pl])){
	  dualfprintf(fail_file,"not finite: i=%d j=%d k=%d :: dir=%d pl=%d %g : %g : %g %g :: %g %g\n",i,j,k,dir,pl,localpdot[dir][pl],flux[i][j][k][pl],geom.g,mysurface,Ftemp[pl],Ftempdiag[pl]);
	  }
	  */

	



	  ////////////
	  //
	  // do term-by-term flux for physics accounting
	  // somewhat like dumps and somewhat like avg of terms
	  // GODMARK: For finite volume method this does not differentiate between point and average values!
	  //







	  if(enerregion==0){
	    bsqorho=bsq/pr[RHO]; // b^2/\rho
	    // we assume user will check if this condition makes sense for a particular simulation.
	    // condition answer for recording for jet2 region
	    //
	    // a real analysis would trace a field line from the horizon to the outer edge in the funnel-type region and only include the region inside, where eventually we have at the outer edge (or will have) unbound/outbound flow.
	    if(dir==X1DN){
	      condjet2=(bsqorho>JETBSQORHO); // assumes that plunging region never develops such large values.  Can occur, but generally not so.  Can raise if problems.
	    }
	    else if(dir==X1UP){ // assumes in jet2 region, but non-jet2 region could have this property.
	      condjet2=((q.ucon[1]>0.0)&&(-q.ucov[0]-1.0>0.0)); // outgoing and unbound at outer edge
	    }
	    else{
	      condjet2=0;
	    }
	  }
	  else condjet2=0; // never touches pdottermsjet2 then
	  
	  
	  surgdet=(geom.g)*mysurface;

	  // loop and if's since some calculations are redundantly simple for similar pl
	  PDIAGLOOP(pl){
	    if(pl==RHO){
	      ftemp0=pr[pl]*(q.ucon[fluxdir])*surgdet;
	      localpdotterms[dir][0][pl]+=ftemp0; // only one part
	      if(condjet2) pdottermsjet2[dir][0][pl]+=ftemp0; // only one part
	      //	  localpdot[dir][pl]  += ftemp0 * surgdet;
	    }
	    // part0-6
	    else if((pl>=UU)&&(pl<=U3)){
	      l=pl-UU;
	      // we currently DO NOT add flux[RHO] to flux[UU], just assume reader knows this is from the native stress
	      ftemp0=pgas*(q.ucon[fluxdir])*(q.ucov[l])*surgdet;
	      localpdotterms[dir][0][pl]+=ftemp0;
	      if(condjet2)	    pdottermsjet2[dir][0][pl]+=ftemp0;

	      ftemp1=p[i][j][k][RHO]*(q.ucon[fluxdir])*(q.ucov[l])*surgdet;
	      localpdotterms[dir][1][pl]+=ftemp1;
	      if(condjet2)	    pdottermsjet2[dir][1][pl]+=ftemp1;


	      ftemp2=p[i][j][k][UU]*(q.ucon[fluxdir])*(q.ucov[l])*surgdet;
	      localpdotterms[dir][2][pl]+=ftemp2;
	      if(condjet2)	    pdottermsjet2[dir][2][pl]+=ftemp2;


	      ftemp3=bsq*(q.ucon[fluxdir])*(q.ucov[l])*surgdet;
	      localpdotterms[dir][3][pl]+=ftemp3;
	      if(condjet2)	    pdottermsjet2[dir][3][pl]+=ftemp3;


	      ftemp4=pgas*delta(fluxdir,pl-UU)*surgdet;
	      localpdotterms[dir][4][pl]+=ftemp4;
	      if(condjet2)	    pdottermsjet2[dir][4][pl]+=ftemp4;


	      ftemp5=0.5*bsq*delta(fluxdir,pl-UU)*surgdet;
	      localpdotterms[dir][5][pl]+=ftemp5;
	      if(condjet2)	    pdottermsjet2[dir][5][pl]+=ftemp5;


	      ftemp6=-(q.bcon[fluxdir])*(q.bcov[l])*surgdet;
	      localpdotterms[dir][6][pl]+=ftemp6;
	      if(condjet2)	    pdottermsjet2[dir][6][pl]+=ftemp6;

	    }
	    else if(pl==B1){
	      ftemp0=(q.bcon[1])*(q.ucon[fluxdir])*surgdet; // flux_b1 term1
	      localpdotterms[dir][0][pl]+=ftemp0;
	      if(condjet2)	    pdottermsjet2[dir][0][pl]+=ftemp0;


	      ftemp1=-(q.bcon[fluxdir])*(q.ucon[1])*surgdet; // flux_b1 term2
	      localpdotterms[dir][1][pl]+=ftemp1;
	      if(condjet2)	    pdottermsjet2[dir][1][pl]+=ftemp1;

	    }
	    else if(pl==B2){
	      ftemp0=(q.bcon[2])*(q.ucon[fluxdir])*surgdet; // flux_b2 term1
	      localpdotterms[dir][0][pl]+=ftemp0;
	      if(condjet2)	    pdottermsjet2[dir][0][pl]+=ftemp0;

	    
	      ftemp1=-(q.bcon[fluxdir])*(q.ucon[2])*surgdet; // flux_b2 term2
	      localpdotterms[dir][1][pl]+=ftemp1;
	      if(condjet2)	    pdottermsjet2[dir][1][pl]+=ftemp1;

	    }
	    else if(pl==B3){
	      ftemp0=(q.bcon[3])*(q.ucon[fluxdir])*surgdet; // flux_b3 term1
	      localpdotterms[dir][0][pl]+=ftemp0;
	      if(condjet2)	    pdottermsjet2[dir][0][pl]+=ftemp0;

	      ftemp1=-(q.bcon[fluxdir])*(q.ucon[3])*surgdet; // flux_b3 term2
	      localpdotterms[dir][1][pl]+=ftemp1;
	      if(condjet2)	    pdottermsjet2[dir][1][pl]+=ftemp1;

	    }

	  }// end PDIAGLOOP over term-by-term fluxes

	}// end GENLOOP

	// cumulative only based upon localpdot
	// localpdot contains entire sum over grid of relevant surface integral value for this direction and ener-region.
	PDIAGLOOP(pl) localpcum[dir][pl]+=localpdot[dir][pl]*Dt; // localpdot is already corrected for REMOVERESTMASSFROMUU

      }// end if localdoflux
    }// end DIRLOOP
  }// end ENERloop
#if(COMPUTEFRDOT)
  // radial flux vs. radius
  flux=F1;
  surface = dx[2]*dx[3];
  for(i=0;i<N1;i++){

    PDIAGLOOP(pl) frdot[i][pl]=0;

    for(j=0;j<N2;j++) for(k=0;k<N3;k++){
      get_geometry(i, j, k, FACE1, &geom);
      PDIAGLOOP(pl) Ftemp[pl]=flux[i][j][k][pl]*surface; // UEVOLVE form
      // GODMARK: for finite volume method, below doesn't change the result if eomfunc=gdet.
      // Otherwise flux would have to be completely recomputed for gdet case JUST for diagnostic to be consistent at higher order
      UtoU_evolve2diag(UEVOLVE,UDIAG,&geom,Ftemp,Ftempdiag); // convert to diag form
      PDIAGLOOP(pl) frdot[i][pl]+=Ftempdiag[pl];
      if(REMOVERESTMASSFROMUU==2){
	//	frdot[i][UU]+= -Ftempdiag[RHO]; // GODMARK: Problem is rest-mass term will dominate and destroy ability to recover energy flux
      }
      else{
	// already correct
	// so fix this instead!
	frdot[i][UU] += Ftempdiag[RHO]; // remove rest-mass term.

      }
    }

  }
#endif
  // GODMARK
  // want all fluxes vs theta on horizon


  return(0);

}


#define DEBUGFRLOOP 1

// write the flux vs. radius
void frdotout(void)
{
  int i,j,k,pl,l;
  SFTYPE ftemp;
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
  SFTYPE frdottemp[N1][NPR];
  SFTYPE *frtot;
  int ospos1;
  FILE*frout;
  static int firsttime=1;

  if(numprocs==1){
    frtot=(SFTYPE (*))(&frdot[0][0]);
  }
  else{
#if(USEMPI)
    if(myid==0){
      frtot=(SFTYPE*) malloc(sizeof(SFTYPE)*totalsize[1]*NPR);
      if(frtot==NULL){
	dualfprintf(fail_file,"Cannot get frtot memory\n");
	myexit(1);
      }
      else{
	for(i=0;i<totalsize[1];i++) PDIAGLOOP(pl){
	  frtot[i*NPR+pl]=0;
	}
      }
      for(l=0;l<numprocs;l++){ // just go over all cpus and assume only added to frdot per cpu for correct cpus.
	ospos1=(l%ncpux1)*N1;
	if(l==0){ // assumes cpu=0 is main cpu and is on horizon
	  for(i=0;i<N1;i++) PDIAGLOOP(pl){
	    frdottemp[i][pl]=frdot[i][pl];
	  }
	}
	else{
	  MPI_Irecv(frdottemp,N1*NPR,MPI_SFTYPE,l,l,MPI_COMM_GRMHD,&rrequest);
	  MPI_Wait(&rrequest,&mpichstatus);
	}
	for(i=0;i<N1;i++) PDIAGLOOP(pl){
	  frtot[(ospos1+i)*NPR+pl]+=frdottemp[i][pl];
#if(DEBUGFRLOOP)
	  if((ospos1+i)*NPR+pl>=totalsize[1]*NPR){
	    dualfprintf(fail_file,"outside bounds: %d\n",(ospos1+i)*NPR+pl);
	    myexit(1);
	  }
#endif
	}
      }
    }
    else{
      MPI_Isend(frdot,N1*NPR,MPI_SFTYPE,0,myid,MPI_COMM_GRMHD,&srequest);
      MPI_Wait(&srequest,&mpichstatus);
    }
#endif
  }
  // now we have frtot with full fluxes vs. radius (totalsize[1]), so output

  if(myid==0){
    frout=fopen("frdot.out","at");
    if(frout==NULL){
      dualfprintf(fail_file,"Cannot open frdot.out\n");
      myexit(1);
    }
    if(firsttime){
      fprintf(frout,"%21.15g %ld %d %d %d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",
	      t,realnstep,totalsize[1],totalsize[2],totalsize[3],startx[1],startx[2],startx[3],dx[1],dx[2],dx[3]);
      fflush(frout);
    }

    for(i=0;i<totalsize[1];i++){
      fprintf(frout,"%21.15g %d ",t,i);
      PDUMPLOOP(pl){// dump only dump prims
	fprintf(frout,"%21.15g ",frtot[i*NPR+pl]);
      }
      fprintf(frout,"\n");
    }
    
    if(frout!=NULL) fclose(frout);
#if(USEMPI)
    if( numprocs != 1 ) free(frtot);  //atch corrected: without the "if" on 1 proc. in MPI this leads to freeing the stack memory
#endif
  }
  firsttime=0;
}





void init_varstavg(void)
{
  int i,j,k,ii;

  ZLOOP{
    for(ii=0;ii<NUMNORMDUMP;ii++){
      normalvarstavg[i][j][k][ii]=0.0;
      anormalvarstavg[i][j][k][ii]=0.0;
    }      
    for(ii=0;ii<NDIM;ii++){
#if(CALCFARADAYANDCURRENTS)
      jcontavg[i][j][k][ii]=0.0;
      jcovtavg[i][j][k][ii]=0.0;
      ajcontavg[i][j][k][ii]=0.0;
      ajcovtavg[i][j][k][ii]=0.0;
#endif
      massfluxtavg[i][j][k][ii]=0.0;
      amassfluxtavg[i][j][k][ii]=0.0;
    }
    for(ii=0;ii<NUMOTHER;ii++){
      othertavg[i][j][k][ii]=0.0;      
      aothertavg[i][j][k][ii]=0.0;      
    }
#if(CALCFARADAYANDCURRENTS)
    for(ii=0;ii<NUMFARADAY;ii++){
      fcontavg[i][j][k][ii]=0.0;
      fcovtavg[i][j][k][ii]=0.0;
      afcontavg[i][j][k][ii]=0.0;
      afcovtavg[i][j][k][ii]=0.0;
    }
#endif
    for(ii=0;ii<NUMSTRESSTERMS;ii++){
      tudtavg[i][j][k][ii]=0.0;
      atudtavg[i][j][k][ii]=0.0;
    }      
    
  }
}

void final_varstavg(FTYPE IDT)
{
  int i,j,k,ii;

  ZLOOP{
    for(ii=0;ii<NUMNORMDUMP;ii++){
      normalvarstavg[i][j][k][ii]=normalvarstavg[i][j][k][ii]*IDT;
      anormalvarstavg[i][j][k][ii]=anormalvarstavg[i][j][k][ii]*IDT;
    }      
    for(ii=0;ii<NDIM;ii++){
#if(CALCFARADAYANDCURRENTS)
      jcontavg[i][j][k][ii]=jcontavg[i][j][k][ii]*IDT;
      jcovtavg[i][j][k][ii]=jcovtavg[i][j][k][ii]*IDT;
      ajcontavg[i][j][k][ii]=ajcontavg[i][j][k][ii]*IDT;
      ajcovtavg[i][j][k][ii]=ajcovtavg[i][j][k][ii]*IDT;
#endif
      massfluxtavg[i][j][k][ii]*=IDT;
      amassfluxtavg[i][j][k][ii]*=IDT;
    }
    for(ii=0;ii<NUMOTHER;ii++){
      othertavg[i][j][k][ii]*=IDT;
      aothertavg[i][j][k][ii]*=IDT;
    }
#if(CALCFARADAYANDCURRENTS)
    for(ii=0;ii<NUMFARADAY;ii++){
      fcontavg[i][j][k][ii]=fcontavg[i][j][k][ii]*IDT;
      fcovtavg[i][j][k][ii]=fcovtavg[i][j][k][ii]*IDT;
      afcontavg[i][j][k][ii]=afcontavg[i][j][k][ii]*IDT;
      afcovtavg[i][j][k][ii]=afcovtavg[i][j][k][ii]*IDT;
    }
#endif
    for(ii=0;ii<NUMSTRESSTERMS;ii++){
      tudtavg[i][j][k][ii]=tudtavg[i][j][k][ii]*IDT;
      atudtavg[i][j][k][ii]=atudtavg[i][j][k][ii]*IDT;
    }      
  }
}


int set_varstavg(FTYPE tfrac)
{
  int i,j,k;
  int iii;
  int ll;
  int l,ii,aii;
  FTYPE ftemp;
  FTYPE ftemp0,ftemp1,ftemp2,ftemp3,ftemp4,ftemp5,ftemp6;
  FTYPE pgas,bsq;
  FTYPE jcov[NDIM];
  FTYPE fcov[NUMFARADAY];
  FTYPE V[NDIM], vmin[NDIM], vmax[NDIM];
  int ignorecourant;
  struct of_geom geom;
  struct of_state q;
  FTYPE X[NDIM];
  FTYPE divb;
  FTYPE b[NDIM],ucon[NDIM];
  FTYPE U[NPR];


  ZLOOP{

    // just like dumps
    coord(i, j, k, CENT, X);
    bl_coord(X, V);
    // if failed, then data output for below invalid, but columns still must exist    
    get_geometry(i, j, k, CENT, &geom);
    if (!failed) {
      if (get_state(pdump[i][j][k], &geom, &q) >= 1)
	FAILSTATEMENT("diag.c:set_varstavg()", "get_state() dir=0", 1);

      if (vchar(pdump[i][j][k], &q, 1, &geom, &vmax[1], &vmin[1],&ignorecourant) >= 1)
	FAILSTATEMENT("diag.c:set_varstavg()", "vchar() dir=1or2or3", 1);
      if (vchar(pdump[i][j][k], &q, 2, &geom, &vmax[2], &vmin[2],&ignorecourant) >= 1)
	FAILSTATEMENT("diag.c:set_varstavg()", "vchar() dir=1or2or3", 2);
      if (vchar(pdump[i][j][k], &q, 3, &geom, &vmax[3], &vmin[3],&ignorecourant) >= 1)
	FAILSTATEMENT("diag.c:set_varstavg()", "vchar() dir=1or2or3", 3);
    }
    else {// do a per zone check, otherwise set to 0
      whocalleducon=1; // force no failure mode, just return like failure, and don't return if failure, just set to 0 and continue
      if (get_state(pdump[i][j][k], &geom, &q) >= 1){
	for (iii = 0; iii < NDIM; iii++)
	  q.ucon[iii]=0;
	for (iii = 0; iii < NDIM; iii++)
	  q.ucov[iii]=0;
	for (iii = 0; iii < NDIM; iii++)
	  q.bcon[iii]=0;
	for (iii = 0; iii < NDIM; iii++)
	  q.bcov[iii]=0;
      }
      if (vchar(pdump[i][j][k], &q, 1, &geom, &vmax[1], &vmin[1],&ignorecourant) >= 1){
	vmax[1]=vmin[1]=0;
      }
	
      if (vchar(pdump[i][j][k], &q, 2, &geom, &vmax[2], &vmin[2],&ignorecourant) >= 1){
	vmax[2]=vmin[2]=0;
      }

      if (vchar(pdump[i][j][k], &q, 2, &geom, &vmax[3], &vmin[3],&ignorecourant) >= 1){
	vmax[3]=vmin[3]=0;
      }


      whocalleducon=0; // return to normal state

    }

    setfdivb(&divb, pdump, udump, i, j, k); // pdump,udump set externally GODMARK


    ii=0;
    aii=0;

    for(iii=0;iii<NPR;iii++){ // always NPR here
      normalvarstavg[i][j][k][ii++]+=pdump[i][j][k][iii]*tfrac;
      anormalvarstavg[i][j][k][aii++]+=fabs(pdump[i][j][k][iii])*tfrac;
    }
    normalvarstavg[i][j][k][ii++]+=divb*tfrac;
    anormalvarstavg[i][j][k][aii++]+=fabs(divb)*tfrac;

    for (iii = 0; iii < NDIM; iii++) normalvarstavg[i][j][k][ii++]+=q.ucon[iii]*tfrac;
    for (iii = 0; iii < NDIM; iii++) anormalvarstavg[i][j][k][aii++]+=fabs(q.ucon[iii])*tfrac;

    for (iii = 0; iii < NDIM; iii++) normalvarstavg[i][j][k][ii++]+=q.ucov[iii]*tfrac;
    for (iii = 0; iii < NDIM; iii++) anormalvarstavg[i][j][k][aii++]+=fabs(q.ucov[iii])*tfrac;

    for (iii = 0; iii < NDIM; iii++) normalvarstavg[i][j][k][ii++]+=q.bcon[iii]*tfrac;
    for (iii = 0; iii < NDIM; iii++) anormalvarstavg[i][j][k][aii++]+=fabs(q.bcon[iii])*tfrac;

    for (iii = 0; iii < NDIM; iii++) normalvarstavg[i][j][k][ii++]+=q.bcov[iii]*tfrac;
    for (iii = 0; iii < NDIM; iii++) anormalvarstavg[i][j][k][aii++]+=fabs(q.bcov[iii])*tfrac;

    normalvarstavg[i][j][k][ii++]+=vmin[1]*tfrac;
    anormalvarstavg[i][j][k][aii++]+=fabs(vmin[1])*tfrac;

    normalvarstavg[i][j][k][ii++]+=vmax[1]*tfrac;
    anormalvarstavg[i][j][k][aii++]+=fabs(vmax[1])*tfrac;

    normalvarstavg[i][j][k][ii++]+=vmin[2]*tfrac;
    anormalvarstavg[i][j][k][aii++]+=fabs(vmin[2])*tfrac;

    normalvarstavg[i][j][k][ii++]+=vmax[2]*tfrac;
    anormalvarstavg[i][j][k][aii++]+=fabs(vmax[2])*tfrac;

    normalvarstavg[i][j][k][ii++]+=vmin[3]*tfrac;
    anormalvarstavg[i][j][k][aii++]+=fabs(vmin[3])*tfrac;

    normalvarstavg[i][j][k][ii++]+=vmax[3]*tfrac;
    anormalvarstavg[i][j][k][aii++]+=fabs(vmax[3])*tfrac;


#if(CALCFARADAYANDCURRENTS)
    lower_vec(jcon[i][j][k],&geom,jcov);
    for(ii=0;ii<NDIM;ii++){
      jcontavg[i][j][k][ii]+=jcon[i][j][k][ii]*tfrac;
      jcovtavg[i][j][k][ii]+=jcov[ii]*tfrac;
      ajcontavg[i][j][k][ii]+=fabs(jcon[i][j][k][ii])*tfrac;
      ajcovtavg[i][j][k][ii]+=fabs(jcov[ii])*tfrac;
    }
#endif
    
    for(ii=0;ii<NDIM;ii++){
      ftemp=(geom.g)*pdump[i][j][k][RHO]*(q.ucon[ii]);
      massfluxtavg[i][j][k][ii]+=ftemp*tfrac;
      amassfluxtavg[i][j][k][ii]+=fabs(ftemp)*tfrac;
    }

    ii=0;
    aii=0;
    ftemp=(q.ucon[3])/(q.ucon[0]);
    othertavg[i][j][k][ii++]=ftemp*tfrac;
    aothertavg[i][j][k][aii++]=fabs(ftemp)*tfrac;

#if(CALCFARADAYANDCURRENTS)
    lowerf(fcon[i][j][k],&geom,fcov);
    for(ii=0;ii<NUMFARADAY;ii++){
      fcontavg[i][j][k][ii]+=fcon[i][j][k][ii]*tfrac;
      fcovtavg[i][j][k][ii]+=fcov[ii]*tfrac;
      afcontavg[i][j][k][ii]+=fabs(fcon[i][j][k][ii])*tfrac;
      afcovtavg[i][j][k][ii]+=fabs(fcov[ii])*tfrac;
    }
#endif

    pgas = pressure_rho0_u(pdump[i][j][k][RHO],pdump[i][j][k][UU]);
    bsq=0; for(iii=0;iii<NDIM;iii++) bsq+=(q.bcon[iii])*(q.bcov[iii]);

    // part0
    ii=0;
    for(iii=0;iii<NDIM;iii++) for(l=0;l<NDIM;l++){
      ftemp0=pgas*(q.ucon[iii])*(q.ucov[l]);
      tudtavg[i][j][k][ii]+=ftemp0*tfrac;
      atudtavg[i][j][k][ii]+=fabs(ftemp0)*tfrac;
      ii++;
    }
    // part1
    for(iii=0;iii<NDIM;iii++) for(l=0;l<NDIM;l++){
      ftemp1=pdump[i][j][k][RHO]*(q.ucon[iii])*(q.ucov[l]);
      tudtavg[i][j][k][ii]+=ftemp1*tfrac;
      atudtavg[i][j][k][ii]+=fabs(ftemp1)*tfrac;
      ii++;
    }
    // part2
    for(iii=0;iii<NDIM;iii++) for(l=0;l<NDIM;l++){
      ftemp2=pdump[i][j][k][UU]*(q.ucon[iii])*(q.ucov[l]);
      tudtavg[i][j][k][ii]+=ftemp2*tfrac;
      atudtavg[i][j][k][ii]+=fabs(ftemp2)*tfrac;
      ii++;
    }
    // part3
    for(iii=0;iii<NDIM;iii++) for(l=0;l<NDIM;l++){
      ftemp3=bsq*(q.ucon[iii])*(q.ucov[l]);
      tudtavg[i][j][k][ii]+=ftemp3*tfrac;
      atudtavg[i][j][k][ii]+=fabs(ftemp3)*tfrac;
      ii++;
    }
    // part4
    for(iii=0;iii<NDIM;iii++) for(l=0;l<NDIM;l++){
      ftemp4=pgas*delta(iii,l);
      tudtavg[i][j][k][ii]+=ftemp4*tfrac;
      atudtavg[i][j][k][ii]+=fabs(ftemp4)*tfrac;
      ii++;
    }
    // part5
    for(iii=0;iii<NDIM;iii++) for(l=0;l<NDIM;l++){
      ftemp5=0.5*bsq*delta(iii,l);
      tudtavg[i][j][k][ii]+=ftemp5*tfrac;
      atudtavg[i][j][k][ii]+=fabs(ftemp5)*tfrac;
      ii++;
    }
    // part6
    for(iii=0;iii<NDIM;iii++) for(l=0;l<NDIM;l++){
      ftemp6=-(q.bcon[iii])*(q.bcov[l]);
      tudtavg[i][j][k][ii]+=ftemp6*tfrac;
      atudtavg[i][j][k][ii]+=fabs(ftemp6)*tfrac;
      ii++;
    }

  }

  return(0);

}


// if doavg==1, then assume this is call before dumping
int average_calc(int doavg)
{
  static FTYPE lastdt;
  static int calls=0;
  static FTYPE tavgi,tavgf;
  static int tavgflag=1;
 
  if(calls>0){ // since need 2 times

    if(tavgflag){
      // gets reached on next call after dump call or first time
      init_varstavg();
      tavgflag=0;
      tavgi=t;
    }

    // always do
    if(set_varstavg(0.5*(lastdt+dt))>=1) return(1);

    if(doavg==1){
      tavgflag=1;
      tavgf=t;
      final_varstavg(1.0/(tavgf-tavgi));
      // expect to dump after this function ends and before next call to this function
    }
  }
  calls++;
  lastdt=dt;

  return(0);
}


void diag_source_all(struct of_geom *ptrgeom, FTYPE *dU,SFTYPE Dt)
{
  int pl,enerregion;
  FTYPE ftemp[NPR];
  FTYPE ftempdiag[NPR];
  SFTYPE *localsourceadd;
  int *localenerpos;
  SFTYPE (*localsourceaddterms)[NPR];


  // does not matter what stage
  if(Dt>0.0){

    //    dualfprintf(fail_file,"got here: i=%d j=%d t=%21.15g\n",ptrgeom->i,ptrgeom->j,t);


    ENERREGIONLOOP(enerregion){
      localenerpos=enerposreg[enerregion];
      localsourceaddterms=sourceaddtermsreg[enerregion];
      localsourceadd=sourceaddreg[enerregion];

      if(WITHINENERREGION(localenerpos,ptrgeom->i,ptrgeom->j,ptrgeom->k)){
	PDIAGLOOP(pl) ftemp[pl]=Dt*dVF*dU[pl]; // in UEVOLVE form
	// GODMARK: for finite volume method, below doesn't change the result if eomfunc=gdet.
	// Otherwise source would have to be completely recomputed for gdet case JUST for diagnostic to be consistent at higher order
	UtoU_evolve2diag(UEVOLVE,UDIAG,ptrgeom,ftemp,ftempdiag); // convert to diag form

	// now assign diagnostic form of source
	PDIAGLOOP(pl){
	  localsourceadd[pl]+=ftempdiag[pl];
#if(DOLUMVSR)
	  // GODMARK: only correct for diagonal coordinate Jacobian in which each i is same radius for all j
	  if(pl==UU) if(enerregion==0){
	      lumvsr[startpos[1]+ptrgeom->i]+=ftempdiag[pl];
	      //	      if(!isfinite(ftempdiag[pl])){
	      //		dualfprintf(fail_file,"ii=%d ftempdiag=%21.15g\n",startpos[1]+ptrgeom->i,ftempdiag[pl]);
		//	      }
	    }
#endif
	} // end PDIAGLOOP on diag
      }
    }
  }

}



void diag_source_comp(struct of_geom *ptrgeom, FTYPE (*dUcomp)[NPR],SFTYPE Dt)
{
  int sc,pl,enerregion;
  FTYPE ftemp[NPR];
  FTYPE ftempdiag[NPR];
  SFTYPE *localsourceadd;
  int *localenerpos;
  SFTYPE (*localsourceaddterms)[NPR];


  // does not matter what stage
  if(Dt>0.0){

    ENERREGIONLOOP(enerregion){
      localenerpos=enerposreg[enerregion];
      localsourceaddterms=sourceaddtermsreg[enerregion];
      localsourceadd=sourceaddreg[enerregion];

      if(WITHINENERREGION(localenerpos,ptrgeom->i,ptrgeom->j,ptrgeom->k)){
	SCLOOP(sc){
	  PDIAGLOOP(pl) ftemp[pl]=Dt*dVF*dUcomp[sc][pl]; // in UEVOLVE form
	  // GODMARK: for finite volume method, below doesn't change the result if eomfunc=gdet.
	  // Otherwise source would have to be completely recomputed for gdet case JUST for diagnostic to be consistent at higher order
	  UtoU_evolve2diag(UEVOLVE,UDIAG,ptrgeom,ftemp,ftempdiag); // convert to diag form

	  // now assign diagnostic form of source
	  PDIAGLOOP(pl){
	    localsourceaddterms[sc][pl]+=ftempdiag[pl];
	  } // end PDIAGLOOP on diag
	} // end SCLOOP
      }
    }
  }

}





// whether to not really do full inversion since can be expensive
#define AVOIDFULLINVERSION 1

// compute dissipated energy due to (e.g.) shocks and reconnection
// If not evolving entropy for full&direct evolution, then this function computes entropy update
// If doing "comparison" then also do dissipation stuff
// In case if doing dissipation stuff, then can compute dissipation in 2 ways (from inversion or from trivialized inversion)
// Output both for now so can make comparison
int diss_compute(int evolvetype, int inputtype, FTYPE *U, struct of_geom *ptrgeom, FTYPE *prbefore, FTYPE *pr)
{
  extern int Utoprimdiss(int evolvetype, int inputtype, FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr, int *otherfail);
  FTYPE prother[NUMDISSVERSIONS][NPR];
  FTYPE Uother[NUMDISSVERSIONS][NPR];
  FTYPE Uenergy[NPR];
  int otherfail;
  int enerregion;
  FTYPE dissenergy[NUMDISSVERSIONS];
  int pl;
  struct of_state q,qsimple,qfull,qenergy;
  FTYPE Unew[NPR];
  FTYPE primtoUcons;
  extern int invertentropyflux_calc(FTYPE entropyflux,int dir, struct of_state *q, FTYPE*pr);
  FTYPE Ugeomfree[NPR];
  extern int ufromentropy_calc(FTYPE entropy, FTYPE *pr);
  extern int entropy_calc(FTYPE *pr, FTYPE *entropy);
  int loopinv;
  int choseninv;
  FTYPE entropy,entropydiss;



#if(DOENTROPY==DONOENTROPY || ENTROPY==-100)
  dualfprintf(fail_file,"Should not be doing diss_compute() unless DOENTROPY!=DONOENTROPY\n");
  myexit(7626762);
#endif



  if(DOENTROPY==DOEVOLVECOMPAREENTROPY){



    ///////////////////////////////////
    //
    // Idea is following:
    //
    // 1) We follow some fluid deformed fluid element that at the new time is the exact shape of the cell
    // 2) Now within that cell the conserved entropy and the entropy deduced from the energy equation can be compared
    //    One compares S=s u^t tracked with S deduced from energy equation either by:
    //    a) s = \rho_0 S/D where D=\rho_0 u^t where \rho_0 is from energy evolution (assumes fluid elements match)
    //    b) Or do full inversion of conserved quantities comparing energy version with entropy version (full inversion seems unreliable)
    // 3) To obtain lab-frame (at infinity) values do:
    //    a) dU/d\tau u_\mu where u is internal energy added extra due to dissipation (assumes what reaches inf is all dissipated energy)
    //    b) One can also compare T^t_t[energy evolution] and T^t_t[entropy version of $u$ otherwise same]
    //    This assumes emission is instantaneous (no cooling timescale) so cooling time doesn't enter into issue of d\tau vs. dt that otherwise would be an issue (i.e. one would have to have physical model of cooling)
    //
    // 4) Entropy in comoving frame can be understood as useful as related to entropy source term to entropy equation
    //    d/d\tau(entropy/rho)=0 -> \nabla_\mu(entropy u^\mu)=0
    //    or \nabla_\mu(entropy u^\mu) = \rho_0 d/d\tau(entropy/rho)
    //    Like with internal energy, assume all entropy goes away on each timestep (so no issue with d\tau vs. dt)
    //    Then source is just \rho_0 \delta (entropy/rho) for the timestep being studied
    //    So \rho_o d(s/\rho_0/dt = ds/dt - s/\rho_0 d\rho_0/dt so if ds is comoving source then term we need is:
    //    ds - s/\rho_0 d\rho_0  , so we need how \rho_0 changed in fluid element from previous time
    //    However, we are not comparing previous time with current time, we are instead comparing the same
    //      fluid element evolved differently.  So the only source of variance between the two is assumed to be
    //      (for SIMPLEINV) the entropy, not density changes. GODMARK: Is this right?
    //
    // I noticed that the shocks can be entropy-violating in general (i.e. total entropy decreases)
    // So also output versions without MAX that match better true conserved entropy
    //
    ////////////////////////////////////



    


    ////////////////////////////////
    //
    // Obtain simple inversion (all other components (velocity, etc.) assumed same in final state -- this is consistent with fluid element being same)
    //
    // This is \partial U/\partial\tau or internal energy changed in comoving frame
    // Source term to energy-momentum equation is (\nabla T)_\mu = dU/d\tau u_\mu
    // so energy-momentum transferred to infinity is: {L_\inf}_\mu = dU/d\tau u_\mu
    // So the conserved energy at infinity is {L_\inf}_t = dU/d\tau u_t
    // 
    //
    ////////////////////////////////

    // set rho,v,B
    // ufromentropy_calc needs prother to at least have RHO set for ideal gas EOS
    PALLLOOP(pl) prother[DISSSIMPLEINVCO][pl]=pr[pl];

    UtoU(inputtype,UNOTHING,ptrgeom,U,Ugeomfree);
    // Ugeomfree[ENTROPY] is \rho_0 u^t (entropy)
    //    invertentropyflux_calc(Ugeomfree[ENTROPY],TT, *q, FTYPE*pr); // don't know q, and don't need to, just do below

    // DEBUG:
    //    dualfprintf(fail_file,"S=%21.15g U[rho]=%21.15g pr[RHO]=%21.15g\n",Ugeomfree[ENTROPY],Ugeomfree[RHO],pr[RHO]);

    // \rho_0 U[ENTROPY]/U[RHO] is entropy density in expected form
    entropy=pr[RHO]*Ugeomfree[ENTROPY]/Ugeomfree[RHO];
    ufromentropy_calc(entropy, prother[DISSSIMPLEINVCO]);
    prother[DISSSIMPLEINVCO][UU]=prother[DISSSIMPLEINVCO][ENTROPY];
    // now prother[DISSSIMPLEINVCO][UU,ENTROPY] contains internal energy from simplified entropy inversion


    ////////////////////
    //
    // next obtain entropy from energy-evolution
    //
    ////////////////////

    // entropydiss is entropy from energy-evolution
    entropy_calc(pr, &entropydiss);




    ////////////////////////////////
    //
    // next obtain internal energy from full inversion (in general fluid elements being compared are actually different, so this can actually be less accurate)
    //
    ////////////////////////////////
    //    PALLLOOP(pl) prother[DISSFULLINVCO][pl]=pr[pl]; // guess
    PALLLOOP(pl) prother[DISSFULLINVCO][pl]=prbefore[pl]; // guess is better chosen from pre-energy evolution
    prother[DISSFULLINVCO][UU]=prother[DISSSIMPLEINVCO][UU]; // guess better with simple entropy version for internal energy
    prother[DISSFULLINVCO][ENTROPY]=prother[DISSSIMPLEINVCO][ENTROPY]; // guess much better with simple entropy version for internal energy
    // invert with entropy evolution
#if(AVOIDFULLINVERSION==0)
    // Noticed this uses too much computational power and may not even have solution causing Newton's method to reach large number of iterations
    Utoprimdiss(evolvetype, inputtype, U,  ptrgeom, prother[DISSFULLINVCO],&otherfail);
#else
    PALLLOOP(pl) prother[DISSFULLINVCO][pl] = prother[DISSSIMPLEINVCO][pl];
    otherfail=0;
#endif
    prother[DISSFULLINVCO][UU]=prother[DISSFULLINVCO][ENTROPY];
    // (prother[DISSFULLINVCO][UU,ENTROPY]) now both contain internal energy as derived from full entropy inversion
    // notice that all other quantities are could also be different (if doentropy==evolvefullentropy), hence the prother[DISSFULLINVCO] variable for temporary storage.




    ////////////////////////////////
    //
    // Get primitive state so have u_\mu for writing energy-momentum source (used for LAB1)
    //
    // Also obtain conserved quantities that assume internal energy is one from entropy evolution versions (used for LAB2)
    //
    ////////////////////////////////

    // assume u^\mu and b^\mu used from get_state() is same for all cases, so just get once
    if (get_stateforUdiss(pr, ptrgeom, &qenergy) >= 1) FAILSTATEMENT("utoprim.c:utoprim()", "get_state()", 1);
    qfull=qenergy; // copy entire q
    qsimple=qenergy; // copy entire q


    //    if (get_stateforUdiss(prother[DISSSIMPLEINVCO], ptrgeom, &qsimple) >= 1) FAILSTATEMENT("utoprim.c:utoprim()", "get_state()", 1);
    if (primtoU(UNOTHING, prother[DISSSIMPLEINVCO], &qenergy, ptrgeom, Uother[DISSSIMPLEINVLAB2]) >= 1) FAILSTATEMENT("utoprim.c:utoprim()", "primtoU()", 1);


    if(otherfail==UTOPRIMNOFAIL && AVOIDFULLINVERSION==0){
      // if didn't fail, then setup qfull and conserved quantity associated with fullinv version of primitive
      //      if (get_stateforUdiss(prother[DISSFULLINVCO], ptrgeom, &qfull) >= 1) FAILSTATEMENT("utoprim.c:utoprim()", "get_state()", 1);
      if (primtoU(UNOTHING, prother[DISSFULLINVCO], &qenergy, ptrgeom, Uother[DISSFULLINVLAB2]) >= 1) FAILSTATEMENT("utoprim.c:utoprim()", "primtoU()", 1);
    }
    else{
      PALLLOOP(pl) Uother[DISSFULLINVLAB2][pl]=Uother[DISSSIMPLEINVLAB2][pl];
    }


    

    // get state for final energy primitive
    if (primtoU(UNOTHING, pr, &qenergy, ptrgeom, Uenergy) >= 1) FAILSTATEMENT("utoprim.c:utoprim()", "primtoU()", 1);

    ///////////////////
    //
    //  choose to store entropy version of internal energy in primitive in case wanted
    //
    ///////////////////
    // just overwrite entropy primitive, leave rest same as from full energy equation
    //    pr[ENTROPY]=prother[DISSFULLINVCO][UU];
    pr[ENTROPY]=prother[DISSSIMPLEINVCO][ENTROPY]; // don't use this currently



    //////////////////////
    //
    // now compare pr[UU] and pr[ENTROPY] with some kind of diagnostic
    //
    //////////////////////
    if(evolvetype==EVOLVEUTOPRIM){
      // then during evolution and pr[UU]-pr[ENTROPY] is relevant to physical calculation
      // store difference

      // report failure to invert
      if(DODISS){
	if(otherfail!=UTOPRIMNOFAIL) dissfunpos[ptrgeom->i][ptrgeom->j][ptrgeom->k][DISSFAILUREINV]+=1.0;
      }


      //////////////
      //
      // obtain various versions of dissipation
      //
      //////////////
      for(loopinv=0;loopinv<NUMDISSVERSIONS;loopinv++){


	///////////////
	//
	// Choose which prother[] to use
	//
	// only use if inversion succeeded (otherwise assume entropy evolution wanted negative internal energy and so not a good solution)
	//
	//////////////
	if(
	   loopinv==DISSSIMPLEINVCO || loopinv==DISSSIMPLEINVCONOMAX 
	   || loopinv==DISSSIMPLEINVLAB1 || loopinv==DISSSIMPLEINVLAB1NOMAX
	   || loopinv==DISSSIMPLEINVLAB2 || loopinv==DISSSIMPLEINVLAB2NOMAX
	   || loopinv==DISSENTROPYCO || loopinv==DISSENTROPYCONOMAX
	   || loopinv==DISSENTROPYLAB1 || loopinv==DISSENTROPYLAB1NOMAX
	   || loopinv==DISSENTROPYLAB2 || loopinv==DISSENTROPYLAB2NOMAX
	   ){
	  choseninv=DISSSIMPLEINVCO;
	}
	else if(
		loopinv==DISSFULLINVCO || loopinv==DISSFULLINVCONOMAX
		|| loopinv==DISSFULLINVLAB1 || loopinv==DISSFULLINVLAB1NOMAX
		|| loopinv==DISSFULLINVLAB2 || loopinv==DISSFULLINVLAB2NOMAX
		){
	  if(otherfail==UTOPRIMNOFAIL) choseninv=DISSFULLINVCO;
	  else choseninv=DISSSIMPLEINVCO; // still want some result
	}
	else{
	  dualfprintf(fail_file,"In setting up choseninv: No such loopinv=%d\n",loopinv);
	  myexit(72698626);
	}


	if(! (choseninv==DISSSIMPLEINVCO||choseninv==DISSFULLINVCO)){
	  dualfprintf(fail_file,"In setting up choseninv: Chose bad choseninv=%d\n",choseninv);
	  myexit(1865728326);
	}


	////////////////
	//
	// Compute dissipation energy density over timestep
	// So final energy dissipated is found by multiplying by $\detg dV$ done after the below calculation
	//
	// GODMARK: The DISSVSR is only correct for diagonal coordinate Jacobian in which each i is same radius for all j                         
	//
	////////////////

	if(
	   loopinv==DISSSIMPLEINVCO || loopinv==DISSFULLINVCO || loopinv==DISSSIMPLEINVLAB1 || loopinv==DISSFULLINVLAB1
	   || loopinv==DISSSIMPLEINVCONOMAX || loopinv==DISSFULLINVCONOMAX || loopinv==DISSSIMPLEINVLAB1NOMAX || loopinv==DISSFULLINVLAB1NOMAX
	   ){
	  // dissipated energy is difference between entropy version of internal energy and energy version of internal energy
	  // so a shock has this as positive
	  // Notice Sod shock has negative dissipation sometimes at shock front
	  // So use MAX to NOT ALLOW negative dissipation (assume dissipation is really 0 and assume error in entropy calculation)
	  dissenergy[loopinv]=pr[UU]-prother[choseninv][UU];

	  // MAX versions:
	  if(loopinv==DISSSIMPLEINVCO || loopinv==DISSFULLINVCO || loopinv==DISSSIMPLEINVLAB1 || loopinv==DISSFULLINVLAB1){
	    dissenergy[loopinv]=MAX(dissenergy[loopinv],0.0);
	  }

 	  // LAB1 multiplications:
	  if(loopinv==DISSSIMPLEINVLAB1 || loopinv==DISSSIMPLEINVLAB1NOMAX){
	    // negative sign on u_t because dU/d\tau opposite sign compared to T^t_t
	    // here we use fluid velocity from simple evolution, which is same fluid velocity as from energy evolution
	    dissenergy[loopinv]*=(-qsimple.ucov[TT]);
	  }
	  else if(loopinv==DISSFULLINVLAB1 || loopinv==DISSFULLINVLAB1NOMAX ){
	    // here we use fluid velocity from entropy evolution as consistent with the meaning of "full" throughout
	    dissenergy[loopinv]*=(-qfull.ucov[TT]);
	  }

	}
	else if(loopinv==DISSENTROPYCO || loopinv==DISSENTROPYCONOMAX){
	  // assume rest are entropy generation quantities
	  // assume entropy must increase
	  dissenergy[loopinv]= entropydiss-entropy;

	  if(loopinv==DISSENTROPYCO){
	    dissenergy[loopinv]=MAX(dissenergy[loopinv],0.0);
	  }
	}
	else if(loopinv==DISSENTROPYLAB1 || loopinv==DISSENTROPYLAB1NOMAX){
	  // assume rest are entropy generation quantities
	  // assume entropy must increase
	  dissenergy[loopinv]= entropydiss-entropy; // same as comoving version as source term of entropy equation w.r.t. comparing fluid elements (not different times) That is, I set the d\rho_0/dt term to zero. GODMARK?

	  if(loopinv==DISSENTROPYLAB1){
	    dissenergy[loopinv]=MAX(dissenergy[loopinv],0.0);
	  }

	}
	else if(loopinv==DISSSIMPLEINVLAB2 || loopinv==DISSSIMPLEINVLAB2NOMAX){
	  // negative sign is because -T^t_t>0 for energy>0
	  // otherwise order of U's are so dissipation is positive
	  dissenergy[loopinv] = -(Uenergy[UU] - Uother[DISSSIMPLEINVLAB2][UU]); // could have used Ugeomfree[UU] too

	  if(loopinv==DISSSIMPLEINVLAB2){
	    dissenergy[loopinv]=MAX(dissenergy[loopinv],0.0);
	  }
	}
	else if(loopinv==DISSFULLINVLAB2 || loopinv==DISSFULLINVLAB2NOMAX){
	  // negative sign is because -T^t_t>0 for energy>0
	  // otherwise order of U's are so dissipation is positive
	  // leave-off MAX
	  // note stored U in Uother for SIMPLEINV for simplicity
	  dissenergy[loopinv] = -(Uenergy[UU] - Uother[DISSFULLINVLAB2][UU]); // could have used Ugeomfree[UU] too

	  if(loopinv==DISSFULLINVLAB2){
	    dissenergy[loopinv]=MAX(dissenergy[loopinv],0.0);
	  }
	}
	else if(loopinv==DISSENTROPYLAB2 || loopinv==DISSENTROPYLAB2NOMAX){
	  // computing S from U and S from Uother and taking difference
	  // order of U's are so dissipation is positive
	  // leave-off MAX
	  // using SIMPLEINV here like used SIMPLEINV to obtain entropy comoving version of dissipation
	  // Here Uenergy[ENTROPY] is conserved entropy as computed from energy-evolved pr's
	  dissenergy[loopinv] = (Uenergy[ENTROPY] - Uother[DISSSIMPLEINVLAB2][ENTROPY]);

	  if(loopinv==DISSENTROPYLAB2){
	    dissenergy[loopinv]=MAX(dissenergy[loopinv],0.0);
	  }

	}
	else{
	  dualfprintf(fail_file,"In computing dissenergy: No such loopinv=%d\n",loopinv);
	  myexit(72698627);	  
	}


	  
	// only for enerregion==0
	if(DODISSVSR) dissvsr[loopinv][startpos[1]+ptrgeom->i]+=dissenergy[loopinv]*ptrgeom->g * dVF;
	  
	if(DODISS){
	  // local integral
	  ENERREGIONLOOP(enerregion){
	    diss=dissreg[enerregion];
	    if( WITHINENERREGION(enerpos,ptrgeom->i,ptrgeom->j,ptrgeom->k) ){
	      diss[loopinv]+=dissenergy[loopinv]*ptrgeom->g * dVF; // actual energy
	    }
	  }
	    
	  // function over all space to be written as dump file
	  // energy density, which can be integrated in SM since grid is given
	  // multiply by gdet*dV since gdet may change in time!
	  dissfunpos[ptrgeom->i][ptrgeom->j][ptrgeom->k][loopinv]+=dissenergy[loopinv]*ptrgeom->g * dVF;
	}


	// DEBUG:
	//	dualfprintf(fail_file,"dissenergy[%d][%d][%d][%d]=%21.15g :: simple=%21.15g full=%21.15g energy=%21.15g\n",ptrgeom->i,ptrgeom->j,ptrgeom->k,loopinv,dissfunpos[ptrgeom->i][ptrgeom->j][ptrgeom->k][loopinv],prother[DISSSIMPLEINVCO][ENTROPY],prother[DISSFULLINVCO][ENTROPY],pr[UU]);
	  

  
      }// end loop over versions
    }// end if evolving
  }// end if doing comparison


  // now must redefine U[ENTROPY] so consistent with p(U[normalgrmhd])
  if(DOENOFLUX!=NOENOFLUX){
    //    primtoUcons=UENTROPY; // not UENTROPY since this doesn't have gdet factor!
    if (get_stateforUdiss(pr, ptrgeom, &q) >= 1) FAILSTATEMENT("utoprim.c:utoprim()", "get_state()", 1);
    if (primtoU(UEVOLVE, pr, &q, ptrgeom, Unew) >= 1) FAILSTATEMENT("utoprim.c:utoprim()", "primtoU()", 1);
    U[ENTROPY] = Unew[ENTROPY]; // now conserved entropy is consistent with real primitive state
  }
  // this is done automatically if doing NOENOFLUX since U is obtained again from new primitives.


  return(0);

}


