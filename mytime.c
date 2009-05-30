
#include "decs.h"

// global code timing function
int timecheck(int whichlocation, SFTYPE comptstart)
{
  static SFTYPE tlasttime;
  static long long int startnstep;
  long long int diffnstep;
  FTYPE ftemp;
  FILE* perfout;
  FILE *fileit;
  int i,ii;
  int badlog[NUMDTDS];
  char temps[MAXFILENAME];
  int numzones;
  int dostepout;
  int itemp;

  FTYPE var_ete,var_wt,var_tuperhr,var_eff;
  long long int var_zc;
  
  FTYPE var_lete,var_lwt,var_ltuperhr,var_leff;
  long long int var_lzc;
  
  ////////////////////////
  //
  // BEGIN TIMING STUFF
  //
#if(TIMEMETHOD==0)
  static time_t timestart,timestop;
  static time_t gtimestart,gtimestop,checktime;
#elif(TIMEMETHOD==1)
  static struct timeval timestart,timestop, gtimestart,gtimestop,checktime; struct timezone tz;
#elif((TIMEMETHOD==2)||(TIMEMETHOD==3))
  static clock_t timestart,timestop, gtimestart,gtimestop,checktime;
#endif
  static SFTYPE walltime=0,walltimelocal=0,walltot=0;
  // general time reports
  static clock_t usertmstimestart,usertmstimestop,systmstimestart,systmstimestop;
  static struct timeval wttimestart,wttimestop;
#if(TIMEMETHOD!=1)
  static struct timezone tz;
#endif
  // END TIMING STUFF
  //
  ////////////////////////



  if(myid==0){ // only report time/performance for myid==0

    
    if(whichlocation!=STARTTIME){
      diffnstep=nstep-startnstep;
    }


    if(whichlocation==STARTTIME){
      // SETUP TIME TRACKING
      GETTIME(&timestart);
      GETTIME(&gtimestart);
      microtime(&wttimestart);
      myustimes2(&usertmstimestart,&systmstimestart);
      startnstep=nstep; // since don't store the above starting times upon restart, but DO restart with correct nstep/realnstep, must offset
      diffnstep=0;
    }
    else if(whichlocation==STOPTIME){
      // get final time
      GETTIME(&timestop);
      microtime(&wttimestop);
      myustimes2(&usertmstimestop,&systmstimestop);
    }
    else if(whichlocation==CHECKTIME){


      // check up on how many timesteps per second so can calibrate period of step, perf, and gocheck
      if((diffnstep==1)||( (!(diffnstep%NTIMECHECK))) ){
	// no need to synch cpus since should be close and no MPI calls used so don't require exact synch
	//just use average time, instead of local time, to do a timestep
	GETTIME(&checktime);
	walltime=(SFTYPE) DELTATIME(checktime,timestart);
	if(walltime<1E-5){
	  dualfprintf(fail_file,"Warning: walltime=%21.15g < 1E-5\n",walltime);
	  walltime=1E-5;
	}
	// now calibrate everything
	NTIMECHECK=(int)((SFTYPE)(diffnstep)*DTtimecheck/walltime); // took walltime to go nsteps, so check every so steps corresponding to desired time.
	if(DTstep!=0.0) NDTCCHECK=(int)((SFTYPE)(diffnstep)*DTstep/walltime);
	if(DTstepdot!=0.0) NDTDOTCCHECK=(int)((SFTYPE)(diffnstep)*DTstepdot/walltime);
	if(DTperf!=0.0) NZCCHECK=(int)((SFTYPE)(diffnstep)*DTperf/walltime);
	if(DTgocheck!=0.0) NGOCHECK=(int)((SFTYPE)(diffnstep)*DTgocheck/walltime);
	
	if(NTIMECHECK<1) NTIMECHECK=1;
	if(NDTCCHECK<1) NDTCCHECK=1;
	if(NDTDOTCCHECK<1) NDTDOTCCHECK=1;
	if(NZCCHECK<1) NZCCHECK=1;
	if(NGOCHECK<1) NGOCHECK=1;
	
	// DEBUG:
	//fprintf(stderr,"WALLTIME: %21.15g :: %d %d %d %d %d\n",walltime,NTIMECHECK,NDTCCHECK,NDTDOTCCHECK,NZCCHECK,NGOCHECK);

      }


    }
    else if(whichlocation==SPEEDTIME){
    




      // speed check
      // setup so can plot in sm
      if(DOLOGPERF&&(!( (diffnstep)%NZCCHECK))){
	GETTIME(&gtimestop);
	// running average
	// running average zonecycle rate
	walltime=(SFTYPE) DELTATIME(gtimestop,timestart);
	if(walltime<1E-5) walltime=1E-5;
	// local zonecycle rate
	walltimelocal=(SFTYPE) DELTATIME(gtimestop,gtimestart);
	if(walltimelocal<1E-5) walltimelocal=1E-5;
	for(i=1;i<=1;i++){ // don't really want perf for i==0
	  if(i==0){
	    fileit=log_file;
	    dostepout=1;
	    numzones=N1*N2*N3; // GODMARK not really right, but not used
	  }
	  else if(i==1){
	    fileit=logperf_file;
	    dostepout=1;
	    numzones=realtotalzones;
	  }
	  if(dostepout){

	    var_ete=((tf-t+1.0E-6)/(t-comptstart+1.0E-6)*walltime*SEC2HOUR);
	    var_wt=walltime*SEC2HOUR;
	    var_zc=(int)((FTYPE)(numzones)*(FTYPE)(diffnstep)/walltime);
	    var_tuperhr=((t-comptstart)/(walltime*SEC2HOUR));
	    var_eff=var_zc/((FTYPE)(numprocs)*ZCPSESTIMATE); // estimate only as good as estimate of ZCPSESTIMATE
	  
	    var_lete=((tf-t+1.0E-6)/(t-tlasttime+1.0E-6)*walltimelocal*SEC2HOUR);
	    var_lwt=walltimelocal*SEC2HOUR;
	    var_lzc=(int)((FTYPE)(numzones)*(FTYPE)(NZCCHECK)/walltimelocal);
	    var_ltuperhr=((t-tlasttime)/(walltimelocal*SEC2HOUR));
	    var_leff=var_lzc/((FTYPE)(numprocs)*ZCPSESTIMATE); // estimate only as good as estimate of ZCPSESTIMATE

	    myfprintf(fileit,"#t, ete, n, wt, zc, tu/hr, Eff :  lete, ln, lwt, lzc, ltu/hr, lEff\n");
	    myfprintf(fileit,"%15.10g %15.10g %10ld %15.10g %10lld %10.5g %10.5g"
		    " %15.10g %5d %15.10g %10lld %10.5g %10.5g\n"
		    ,t
		    ,var_ete
		    ,nstep
		    ,var_wt
		    ,var_zc
		    ,var_tuperhr
		    ,var_eff
		    
		    ,var_lete
		    ,NZCCHECK
		    ,var_lwt
		    ,var_lzc
		    ,var_ltuperhr
		    ,var_leff
		    );
	    ftemp=walltimelocal/(t-tlasttime); // seconds per time unit

	    // how many seconds, below which if logging is taking place could impact performance and so is reported to user in the perf file
	    // comment shows what the time check is based upon
	    // the number inside is the # of seconds, the divisor is a scale factor(see main.c)
	  
	    // really want walltime(to do step)/dt << walltime(to do dump)/Dt
	    for(ii=0;ii<NUMDTDS;ii++){
	      if(ftemp<1/DTdumpgen[ii]) badlog[ii]=1; else badlog[ii]=0; 
	      // output to perf file if that log is impacting performance
	      if(badlog[ii]){
		myfprintf(fileit,"#LOGPERF!: %d\n",ii);
	      } 
	    }
	  }
	}
	GETTIME(&gtimestart);
	tlasttime=t;
      }// end if output speed
  

      if(PERFTEST){
	GETTIME(&checktime);
	walltime=(SFTYPE) DELTATIME(checktime,timestart);
	if(walltime<1E-5) walltime=1E-5;
	// setup so each turn is about the same WALLtime
	//itemp=(int)((SFTYPE)((SFTYPE)PERFWALLTIME/(walltime/(SFTYPE)(diffnstep))));
	// setup so each turn is same as estimated time and speed, but fixed timesteps for all runs: best for benchmark if you know ahead of time ZCPSESTIMATE, and use same value of this across all tests
	itemp=(int)((SFTYPE)(PERFWALLTIME*ZCPSESTIMATE) /( (SFTYPE)realtotalzones) );
	//    fprintf(stdout,"itemp: %d PWT: %d ZCPS: %d realtotalzones: %d\n",itemp,PERFWALLTIME,ZCPSESTIMATE,realtotalzones);
	if(itemp<1) itemp=1;
	// set final time so steps to desired # of steps
	tf=1.0*(SFTYPE)(itemp)/(SFTYPE)(diffnstep)*t;
	fprintf(stderr,"PERFTEST: nstep=%ld/%d t=%15.10g/%15.10g wt=%15.10g/%15.10g\n",nstep,itemp,t,tf,walltime,(SFTYPE)itemp*(SFTYPE)walltime/(SFTYPE)nstep); fflush(stderr);
      
	//    if(itemp>1000) itemp=1000;
	if(diffnstep==itemp) reallaststep=1;
	if(walltime>PERFWALLTIME) reallaststep=1;
	// GODMARK(commented) 
	//if(diffnstep==100) reallaststep=1;
      }
      // GODMARK(commented)
      //    if(diffnstep==1) reallaststep=1;




    }



    else if(whichlocation==REPORTTIME){



      if(PERFTEST){
	sprintf(temps,"%sfinalperf%s",DATADIR,".txt") ;
	if(!(perfout=fopen(temps,"at"))){
	  fprintf(fail_file,"Can't open %s\n",temps);
	  exit(1);
	}
	fprintf(stderr, "opened: %s\n", temps);
      }


      // running average zonecycle rate
      walltime=(SFTYPE) DELTATIME(timestop,timestart);
      if(walltime<1E-5) walltime=1E-5;
  

      fprintf(logfull_file,"#allproc: steps: %10ld wtime: %10.2g zcycles: %10d t: %10.2g\n",nstep,walltime*SEC2HOUR,(int)((FTYPE)(realtotalzones)*(FTYPE)diffnstep/walltime),(t-comptstart)) ;
#ifndef WIN32
      fprintf(logfull_file,"#(sec) walltime: %21.15g usertime: %21.15g systime: %21.15g\n",diffmicrotime(wttimestop,wttimestart),diffmyustimes(usertmstimestop,usertmstimestart),diffmyustimes(systmstimestop,systmstimestart));
#endif
      fprintf(stderr,"#(sec) walltime: %21.15g usertime: %21.15g systime: %21.15g\n",diffmicrotime(wttimestop,wttimestart),diffmyustimes(usertmstimestop,usertmstimestart),diffmyustimes(systmstimestop,systmstimestart));

      if(DOLOGPERF){
	myfprintf(logperf_file,"#done: steps: %10ld wtime: %10.2g zcycles: %10d t: %10.2g tu/hour: %10.5g\n",nstep,walltime*SEC2HOUR,(int)((FTYPE)(realtotalzones)*(FTYPE)diffnstep/walltime),(t-comptstart),(t-comptstart)/(walltime*SEC2HOUR)) ;
      }
      if(PERFTEST){
	myfprintf(perfout,"%10d\n",(int)((FTYPE)(realtotalzones)*(FTYPE)diffnstep/walltime)) ;
	myfprintf(stderr,"perf: N3: %d N2: %d N1: %d RTZ: %d ZCPS: %d steps: %ld walltime: %15.10g\n",N3,N2,N1,realtotalzones,(int)((FTYPE)(realtotalzones)*(FTYPE)diffnstep/walltime),nstep,walltime) ;
	fclose(perfout);
      }
      if(DOLOGSTEP){
	myfprintf(logstep_file,"#done: steps: %10ld wtime: %10.2g zcycles: %10d t: %10.2g\n",nstep,walltime*SEC2HOUR,(int)((FTYPE)(realtotalzones)*(FTYPE)diffnstep/walltime),(t-comptstart)) ;
      }

    }
  }// end if myid==0



  // global parameters of importance to all CPUs
#if(USEMPI)
  if(whichlocation==SPEEDTIME && PERFTEST){
    MPI_Bcast(&reallaststep,1,MPI_INT,0,MPI_COMM_GRMHD);
  }
#endif

  return(0);

}




void mycpuclock(clock_t *time)
{
  *time=clock();
}

#ifndef WIN32
void myustimes(clock_t *time) // returns number of microseconds
{
	struct tms mytimes;
	clock_t ret;
	long clockspersecond;

	clockspersecond=sysconf(_SC_CLK_TCK);
	ret=times(&mytimes);
	*time=(clock_t) (1000000.0*(SFTYPE)(mytimes.tms_utime+mytimes.tms_stime+mytimes.tms_cutime+mytimes.tms_cstime)/(SFTYPE)clockspersecond );
}
void myustimes2(clock_t *usertime,clock_t *systime) // returns number of microseconds
{
	struct tms mytimes;
	clock_t ret;
	long clockspersecond;

	clockspersecond=sysconf(_SC_CLK_TCK);
	ret=times(&mytimes);
	*usertime=(clock_t) (1000000.0*(SFTYPE)(mytimes.tms_utime+mytimes.tms_stime)/(SFTYPE)clockspersecond );
	*systime=(clock_t) (1000000.0*(SFTYPE)(mytimes.tms_cutime+mytimes.tms_cstime)/(SFTYPE)clockspersecond );
}
#else
void myustimes(clock_t *time) // returns number of microseconds
{
}
void myustimes2(clock_t *usertime,clock_t *systime) // returns number of microseconds
{
}
#endif



#ifdef WIN32

int gettimeofday(struct timeval *tv, struct timezone *tz)
{
  FILETIME ft;
  unsigned __int64 tmpres = 0;
  static int tzflag;
  if (NULL != tv)
    {
      GetSystemTimeAsFileTime(&ft);
      tmpres |= ft.dwHighDateTime;
      tmpres <<= 32;
      tmpres |= ft.dwLowDateTime;
      /*converting file time to unix epoch*/
      tmpres /= 10;  /*convert into microseconds*/
      tmpres -= DELTA_EPOCH_IN_MICROSECS;
      tv->tv_sec = (long)(tmpres / 1000000UL);
      tv->tv_usec = (long)(tmpres % 1000000UL);
    }
  if (NULL != tz)
    {
      if (!tzflag)
	{
	  _tzset();
	  tzflag++;
	}
      tz->tz_minuteswest = _timezone / 60;
      tz->tz_dsttime = _daylight;
    }
  return 0;
}

#endif







// check if time to output step/time/dt info
// setup so can plot in sm
// doensn't need starting time or nstep
int output_steptimedt_info(SFTYPE comptstart)
{
  int i;
  FILE *fileit;
  int dostepout;
  FTYPE strokecount;

  if(myid==0){ // only output step/dt info for myid==0
      
    if(DOLOGSTEP&&( (!(nstep%NDTDOTCCHECK)) )){
      myfprintf(logstep_file,".");
    }
    // GODMARK: comptstart was tstart in pnmhd code
    if(DOLOGSTEP&&( (!(nstep%NDTCCHECK))||(t>=tf-1.0E-7)||(t<=comptstart+1.0E-7) ) ){
      myfprintf(logstep_file,"\n");
      for(i=1;i<=1;i++){ // ==0 and ==2 not done since really not needed
	if(i==0){
	  fileit=log_file;
	  dostepout=1;
	  strokecount=1.0*nstroke/(2.0*N1*N2);
	}
	else if(i==1){
	  fileit=logstep_file;
	  dostepout=1;
	  strokecount=1.0*nstroke/(2.0*totalzones);
	}
	else if(i==2){
	  fileit=stderr;
	  dostepout=1;
	  strokecount=1.0*nstroke/(2.0*totalzones);
	}

	if(dostepout){
	  myfprintf(fileit,"#t dt cour nstep realnstep strokeperzone:\n"
		  "%21.15g %21.15g %21.15g %8ld %8ld %21.15g\n", t, dt, cour, nstep, realnstep,strokecount );
	
	  if(i==1) myfprintf(fileit,"#"); // for "." to be commented in SM
	}
      }// end over files
    }// end if outputting
  }

  return(0);
}

