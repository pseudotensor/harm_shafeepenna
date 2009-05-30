#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/////////////////
//
// chooseable:
//
/////////////////
#define HDF (1)
#define V5D (1)





/////////////////
//
// rest NOT chooseable unless you know what you are doing:
//
/////////////////
#ifndef _SWAP_ENDIAN
#define _SWAP_ENDIAN

// Macs and SGIs are Big-Endian; PCs are little endian
// returns TRUE if current machine is little endian
static int IsLittleEndian(void);

/******************************************************************************
  FUNCTION: SwapEndian
  PURPOSE: Swap the byte order of a structure
  EXAMPLE: float F=123.456;; SWAP_FLOAT(F);
******************************************************************************/

#define SWAP_SHORT(Var)  Var = *(short*)         SwapEndian((void*)&Var, sizeof(short))
#define SWAP_USHORT(Var) Var = *(unsigned short*)SwapEndian((void*)&Var, sizeof(short))
#define SWAP_LONG(Var)   Var = *(long*)          SwapEndian((void*)&Var, sizeof(long))
#define SWAP_ULONG(Var)  Var = *(unsigned long*) SwapEndian((void*)&Var, sizeof(long))
#define SWAP_RGB(Var)    Var = *(int*)           SwapEndian((void*)&Var, 3)
#define SWAP_FLOAT(Var)  Var = *(float*)         SwapEndian((void*)&Var, sizeof(float))
#define SWAP_DOUBLE(Var) Var = *(double*)        SwapEndian((void*)&Var, sizeof(double))

static void *SwapEndian(void* Addr, const int Nb);

#endif


#define LITTLE_ENDIAN_USER 0
#define BIG_ENDIAN_USER    1


// this program does not convert between data types, only between file formats
// for HDF we have to choose what datatype is read/written, and it's fixed for all output/input types

#if(HDF)
// note, for HDF, dim0,1,2,3,... means the array element array[dim0][dim1][dim2][dim3], which means dimN where N is the dimensionality of the vector, is most quickly iterated.
// so a loop like for(k) for(j) for(i) would be used on array[k][j][i]
// so i is actually dim2, j is dim1, and k is dim0
#include "hdf/hdf.h"
#include "hdf/mfhdf.h"

#define HDFTYPE DFNT_FLOAT32
//#define HDFTYPE DFNT_FLOAT64
// could be DFNT_BYTE DFNT_INT32 or DFNT_FLOAT64 for doubles
#define SDS_NAME      "Jonathan McKinney SDS"
#endif


#if(V5D)
#include <vis5d+/v5d.h>
#include <vis5d+/binio.h>
#endif

// note that for binary, should be concerned with big or little endianness
#define BYTE unsigned char
static void swap(BYTE *x, BYTE size);
#define swaporder16(A)  ( (((A) & 0xff00) >> 8) | (((A) & 0x00ff) << 8) )
#define swaporder32(A)  ( (((A) & 0xff000000) >> 24) | (((A) & 0x00ff0000) >> 8) | (((A) & 0x0000ff00) << 8) | (((A) & 0x000000ff) << 24))

#define DEBUG (0)


// debug compile                       gcc -O0 -g -Wall -o bin2txt bin2txt.c -lmfhdf -ldf -ljpeg -lz  -lv5d

// standard compile on ki-rh42:
// gcc -O2 -Wall -o bin2txt bin2txt.c -I /usr/include/hdf/ -L /usr/lib64/hdf/ -lmfhdf -ldf -ljpeg -lz -lv5d
 
// standard compile on ki-rh39:
// gcc -O2 -Wall -o bin2txt bin2txt.c -I /usr2/local/include/ -L /usr2/local/lib/ -lmfhdf -ldf -ljpeg -lz -lv5d

// JCM 10-4-01

// JCM 10-30-01 added hdf4 capability

#define BUFFERMAP ((k*N2*N1+j*N1+i)*numcolumns+nextbuf++)
#define BUFFERINIT nextbuf=0



#define MAXTERMS 100

int main(
         int argc,
         char *argv[],
         char *envp[]
         )
{
  int machineEndianness(void);
  int argstep=0;
  int NDIM;
  int i,j,k;
  int numrows,numcolumns,numterms,numheader,pipeheader,realnumheader;
  int bytesize,intsize,floatsize,doublesize;
  unsigned char dumb;
  unsigned short shortfordumb;
  int dumi;
  float dumf;
  double dumlf;
  char precision[MAXTERMS];
  int group[MAXTERMS];

  int N1,N2,N3,TS;
  int nextbuf;
  FILE * input;
  FILE * output;
  FILE * inputfilenamelist;
  FILE * outputfilenamelist;
  char INPUT_NAME[400];
  char OUTPUT_NAME[400];
  char INPUTLIST_NAME[400];
  char OUTPUTLIST_NAME[400];
  char V5DHEAD_NAME[400];
  int source,dest,ble;

  short *arrayb;
  int *arrayi;
  float *arrayf;
  double *arrayd;
  void * array;
  float *arrayvisf;
  void * arrayvis;
  float *arrayvisoutput;

#if(HDF)

  // hdf
  int32 sd_id, sds_id, sd_index,istat;
  intn  status;
  int32 rank;
  int32 dim_sizes[4], start[4], edges[4];
#endif
  int it;
#if(V5D)
  // vis5d
  int iv;
  int NumTimes;                      /* number of time steps */
  int NumVars;                       /* number of variables */
  int Nr, Nc, Nl[MAXVARS];           /* size of 3-D grids */
  char VarName[MAXVARS][MAXVARNAME];         /* names of variables */ 
  int TimeStamp[MAXTIMES];           /* real times for each time step */
  int DateStamp[MAXTIMES];           /* real dates for each time step */
  int CompressMode;                  /* number of bytes per grid */
  int Projection;                    /* a projection number */
  float ProjArgs[MAXPROJARGS];       /* the projection parameters */
  int Vertical;                      /* a vertical coord system number */
  float VertArgs[MAXVERTARGS];         /* the vertical coord sys parameters */
  FILE * vis5dheader;
  float pos[3+1][2+1];
  float minmax[2][MAXVARS];
  float a,b;
  int ijkvis5d,ijkjon;
#endif
  

  // r8 stuff
  char incommand[300],outcommand[300];
  int lname;
  int inputtype,outputtype;
  char ch;

  // begin

  bytesize=sizeof(unsigned char);
  intsize=sizeof(int);
  floatsize=sizeof(float);
  doublesize=sizeof(double);


  if(argc<11){
    fprintf(stderr,"non-V5D Usage: bin2txt SOURCE DEST BLE HEADERLINESTOSKIP NDIM N1 N2 N3 TS INPUTNAME OUTPUTNAME <format>\n");
    fprintf(stderr," SOURCE/DEST: 0=r8 (.gz) 1=binary 2=text 3=HDF4 4=HDF5 5=vis5d\n");
    fprintf(stderr," BLE: same endian (0), big to little or little to big (same process) (1), auto-conversion assuming ROMIO binary (-1)\n");
    fprintf(stderr," HEADERLINESTOSKIP: # of header text lines in r8 or text source or binary header (negative values mean pipe text header into new file)\n");
    fprintf(stderr," NDIM/N1/N2/N3/TS: # of dimensions and elements in each dimension (e.g. if 2D, N3's value doesn't matter but needs to be there) and # of timesteps to include TS must be >0\n");
    fprintf(stderr," INPUTNAME/OUTPUTNAME: input and output file names (not used if TS>1)\n");
    fprintf(stderr," <format> type is (for 1-4) b=byte i=integer, f=float, d=double, where you give <type1> #of1 <type2> #of2 ...\n");
    fprintf(stderr,"\ne.g. bin2txt 1 2 0 3 32 32 32 1 data32-3.txt data32-3.bin i 3 f 9\n\n");
    fprintf(stderr,"Vis5D with TS=1:\n");
    fprintf(stderr,"Usage (dest=5): bin2txt SOURCE DEST BLE HEADERLINESTOSKIP NDIM N1 N2 N3 TS HEADFILE INPUTNAME OUTPUTNAME <format>\n");
    fprintf(stderr,"HEADFILE contains return delimited list of names for all variables and their min and max values (min/max only used if source=0) (names MUST BE <=9 characters!!!)\n");
    fprintf(stderr,"Last line must have:\nx_start x_finish y_start y_finish z_start z_finish\n");
    fprintf(stderr,"e.g. for one variable:\n\ndensity 1E-4 1\n");
    fprintf(stderr,"0 1 0 1 0 1\n");
    fprintf(stderr,"V5D w/ TS>1 then also include:\n");
    fprintf(stderr,"Usage (dest=5): bin2txt SOURCE DEST BLE HEADERLINESTOSKIP NDIM N1 N2 N3 TS HEADFILE INPUTLIST OUTPUTNAME <format>\n");
    fprintf(stderr,"Usage (source=5): bin2txt SOURCE DEST BLE HEADERLINESTOSKIP NDIM N1 N2 N3 TS HEADFILE INPUTNAME OUTPUTLIST <format>\n");
    fprintf(stderr,"INTPUTLIST contains: sequence of return delimited input file names for each timeslice\n");
    fprintf(stderr,"OUTPUTLIST contains: sequence of return delimited output file names for each timeslice\n");
    fprintf(stderr,"vis5d e.g.: ./bin2txt 0 5 0 3 64 64 64 1 vis5d.image.head imx0-0-0-s1-0000.dat.r8.gz imx0-0-0-s1-0000.dat.v5d b 1\n");
    fprintf(stderr,"vis5d: bin2txtn 2 5 0 7 3 32 32 32 101 vis5d.dump.head dumplist.txt dump.v5d d 9\n");
    fprintf(stderr,"bin2txtn 2 5 0 0 3 32 32 32 101 vis5d.force.head forcelist.txt forcex.v5d d 8\n");

    exit(1);
  }

  // argv[0] is filename of program
  argstep=1;
  source=atoi(argv[argstep++]);
  dest=atoi(argv[argstep++]);

  ble=atoi(argv[argstep++]);
  if(ble==-1){
    if(machineEndianness()==LITTLE_ENDIAN_USER){
      ble=1; // need to convert
    }
    else{
      ble=0; // no need to convert
    }
  }
  // else use user-inputted ble
  if(ble!=0 && ble!=1){
    fprintf(stderr,"ble=%d undefined\n",ble);
    exit(1);
  }

  numheader=atoi(argv[argstep++]);
  
  if(numheader<0){
    pipeheader=1;
    realnumheader=-numheader;
  }
  else{
    pipeheader=0;
    realnumheader=numheader;
  }
  NDIM=atoi(argv[argstep++]);
  N1=atoi(argv[argstep++]);
  N2=atoi(argv[argstep++]);
  N3=atoi(argv[argstep++]);
  TS=atoi(argv[argstep++]);
  if(TS==0){
    fprintf(stderr,"Assumed by TS=0 you meant TS=1 really\n");
    TS=1;
  }
  fprintf(stderr,"source: %d dest: %d numheader: %d NDIM: %d N1: %d N2: %d N3: %d TS: %d\n",source,dest,numheader,NDIM,N1,N2,N3,TS);
  // modify N1,N2,N3 for various dimension
  if(NDIM==2){
    N3=1; // force
  }
  else if(NDIM==1){
    N2=N3=1; // force
  }

  numrows=N1*N2*N3;

  // get input-output file (or list)
  if((TS>1)&&(dest==5)){
    strcpy(V5DHEAD_NAME,argv[argstep++]); // for vis5d this is the header
    strcpy(INPUTLIST_NAME,argv[argstep++]); // for TS>1 this is list of files
    strcpy(OUTPUT_NAME,argv[argstep++]); // single output name
    fprintf(stderr,"V5DHEAD_NAME: %s INPUTLIST_NAME: %s  OUTPUT_NAME: %s\n",V5DHEAD_NAME,INPUTLIST_NAME,OUTPUT_NAME); fflush(stderr);
  }
  else if((TS>1)&&(source==5)){
    strcpy(V5DHEAD_NAME,argv[argstep++]); // for vis5d this is the header
    strcpy(INPUT_NAME,argv[argstep++]); // just single input name
    strcpy(OUTPUTLIST_NAME,argv[argstep++]); // for TS>1 this is list of files
    fprintf(stderr,"V5DHEAD_NAME: %s INPUT_NAME: %s OUTPUTLIST: %s \n",V5DHEAD_NAME,INPUT_NAME,OUTPUTLIST_NAME); fflush(stderr);
  }
  else if((dest==5)||(source==5)){
    strcpy(V5DHEAD_NAME,argv[argstep++]); // for vis5d this is the header
    strcpy(INPUT_NAME,argv[argstep++]); // just single input name
    strcpy(OUTPUT_NAME,argv[argstep++]); // single output name
    fprintf(stderr,"V5DHEAD_NAME: %s INPUT_NAME: %s OUTPUT_NAME: %s \n",V5DHEAD_NAME,INPUT_NAME,OUTPUT_NAME); fflush(stderr);
  }
  else{ // then don't need v5d header and don't need multi-file list
    strcpy(INPUT_NAME,argv[argstep++]); // just single input name
    strcpy(OUTPUT_NAME,argv[argstep++]); // single output name
    fprintf(stderr,"INPUT_NAME: %s OUTPUT_NAME: %s\n",INPUT_NAME,OUTPUT_NAME); fflush(stderr);
  }


  j=0;
  numcolumns=0;
  while(argstep<argc){  
    sscanf(argv[argstep++],"%c",&precision[j]);
    sscanf(argv[argstep++],"%d",&group[j]);
    numcolumns+=group[j];
    j++;
    if(j>MAXTERMS){
      fprintf(stderr,"too many terms!\n");
      exit(1);
    }
  }

  numterms=j;



  fprintf(stderr,"numrows: %d numcolumns: %d numterms: %d\n",numrows,numcolumns,numterms); fflush(stderr);



#if(V5D)
  /////////////////////////////
  //
  // vis5d stuff
  //
  //
  if(dest==5){
    NumTimes=TS;
    NumVars=numcolumns;
    Nr=N3; // (whine) Avery said use N3 instead of N1
    Nc=N2;  

    fprintf(stderr,"NumTimes=%d  NumVars=%d Nr=%d Nc=%d\n",NumTimes,NumVars,Nr,Nc); fflush(stderr);

    for(i=0;i<NumVars;i++){    
      Nl[i]=N1; // all variables have same # of N3 elements
      // (whine) Avery said use N1 instead of N3
      fprintf(stderr,"Nl[%d]=%d\n",i,Nl[i]);fflush(stderr);
    }

    // read special header file to setup vis5d
    if( (vis5dheader=fopen(V5DHEAD_NAME,"rt"))==NULL){
      fprintf(stderr,"can't open vis5d header file %s\n","vis5d.head");
      exit(1);
    }
    // header has in it:
    // list of names for all variables types in sequence (MUST BE <=9 characters!!!)
    // list of min max for each of the variables
    // x_start x_finish y_start y_finish z_start z_finish
    for(i=0;i<NumVars;i++){
      fscanf(vis5dheader,"%s %f %f",VarName[i],&minmax[0][i],&minmax[1][i]);
      fprintf(stderr,"VarName[%d]=%s min: %g max: %g\n",i,VarName[i],minmax[0][i],minmax[1][i]);  fflush(stderr);
    }    
    fscanf(vis5dheader,"%f",&pos[1][1]);
    fscanf(vis5dheader,"%f",&pos[1][2]);
    fscanf(vis5dheader,"%f",&pos[2][1]);
    fscanf(vis5dheader,"%f",&pos[2][2]);
    fscanf(vis5dheader,"%f",&pos[3][1]);
    fscanf(vis5dheader,"%f",&pos[3][2]);
    while(fgetc(vis5dheader)!='\n'); // skip rest of line
    for(i=0;i<NumTimes;i++){
      TimeStamp[i]=i%60+((i/60)%60)*100+(i/(60*60))*10000;
      //fprintf(stderr,"ts: %06d\n",TimeStamp[i]); fflush(stderr);
      DateStamp[i]=99036;
    }

  }

  if((dest==5)||(source==5)){

    CompressMode=1; // 1,2,4 bytes per grid point
    Projection=0; // 0=linear, rectangular, generic units 1=linear, rectangular,cylindrical-equidistant,2=Lambert Conformal, 3=Stereographic, 4=Rotated
    ProjArgs[0]=pos[1][2];
    ProjArgs[1]=pos[2][2];
    ProjArgs[2]=(pos[1][2]-pos[1][1])/(float)(Nr-1);
    ProjArgs[3]=(pos[2][2]-pos[2][1])/(float)(Nc-1);
    Vertical=0; // 0=equally spaced in generic units 1=equally spaced in km 2=unequally spaced in km 3=unequally spaced in mb
    VertArgs[0]=pos[3][1];
    VertArgs[1]=(pos[3][2]-pos[3][1])/(float)(Nl[0]-1);
    
    fprintf(stderr,"Args: %g %g : %g %g:  %g %g\n",ProjArgs[0],ProjArgs[1],ProjArgs[2],ProjArgs[3],VertArgs[0],VertArgs[1]); fflush(stderr);
    fclose(vis5dheader); // no longer needed
  }

  // see if TS>1 and get filename list
  if(TS>1){
    // then ignore argument file input/output names and get from list
    // order is ALL input names in 1 file, ALL output names in another file
    if(source!=5){ // otherwise 1 file
      fprintf(stderr,"Opening %s\n",INPUTLIST_NAME);
      if( (inputfilenamelist=fopen(INPUTLIST_NAME,"rt"))==NULL){
	fprintf(stderr,"can't open input filenamelist header file %s\n",INPUTLIST_NAME);
	exit(1);
      }
    }
    if(dest!=5){ // otherwise 1 file
      fprintf(stderr,"Opening %s\n",OUTPUTLIST_NAME);
      if( (outputfilenamelist=fopen(OUTPUTLIST_NAME,"rt"))==NULL){
	fprintf(stderr,"can't open output filenamelist header file %s\n",OUTPUTLIST_NAME);
	exit(1);
      }
    }
    // from now on input and output names will be formed from these lists
  }
#endif


#if(HDF)
  if((source==3)||(dest==3)||(source==4)||(dest==4)){ // use maximal holder array (doubles)
    if(HDFTYPE==DFNT_UINT8){
      arrayi=(int*)malloc(sizeof(int)*numcolumns*N1*N2*N3); // (DIM0->N1, DIM1->N2 where array[DIM0][DIM1])
      array=arrayi;
    }
    if(HDFTYPE==DFNT_INT32){
      arrayi=(int32*)malloc(sizeof(int32)*numcolumns*N1*N2*N3); // (DIM0->N1, DIM1->N2 where array[DIM0][DIM1])
      array=arrayi;
    }
    if(HDFTYPE==DFNT_FLOAT32){
      arrayf=(float32*)malloc(sizeof(float32)*numcolumns*N1*N2*N3); // (DIM0->N1, DIM1->N2 where array[DIM0][DIM1])
      array=arrayf;
    }
    if(HDFTYPE==DFNT_FLOAT64){
      arrayd=(float64*)malloc(sizeof(float64)*numcolumns*N1*N2*N3); // (DIM0->N1, DIM1->N2 where array[DIM0][DIM1])
      array=arrayd;
    }
    if(array==NULL){
      fprintf(stderr,"cannot allocate array hdf data\n");
      exit(1);
    }

    if(NDIM==3){
      dim_sizes[3] = numcolumns;
      dim_sizes[2] = N1;
      dim_sizes[1] = N2;
      dim_sizes[0] = N3;
      rank = 1+NDIM;
      
      edges[3] = numcolumns;
      edges[2] = N1;
      edges[1] = N2;
      edges[0] = N3;

      start[0]=start[1]=start[2]=start[3]=0;

    }
    else if(NDIM==2){
      // assum 2d
      dim_sizes[2] = numcolumns;
      dim_sizes[1] = N1;
      dim_sizes[0] = N2;
      rank = 1+NDIM;
      
      edges[2] = numcolumns;
      edges[1] = N1;
      edges[0] = N2;

      start[0]=start[1]=start[2]=0;
    }
    else if(NDIM==1){
      // assum 1d
      dim_sizes[1] = numcolumns;
      dim_sizes[0] = N1;
      rank = 1+NDIM;
      
      edges[1] = numcolumns;
      edges[0] = N1;

      start[0]=start[1]=0;
    }

  }
#endif

#if(V5D)
  if((source==5)||(dest==5)){
    fprintf(stderr,"Allocated memory for source=%d dest=%d\n",source,dest);
    arrayvisf=(float*)malloc(sizeof(float)*numcolumns*N1*N2*N3); // (DIM0->N1, DIM1->N2 where array[DIM0][DIM1])
    arrayvis=arrayvisf;
    if(arrayvis==NULL){
      fprintf(stderr,"cannot allocate array vis5d data\n");
      exit(1);
    }
    arrayvisoutput=(float*)malloc(sizeof(float)*N1*N2*N3); // (DIM0->N1, DIM1->N2 where array[DIM0][DIM1])
  }
#endif


  ///////////
  //
  //
  // BIG LOOP
  //
  ////////////

  for (it=0;it<TS;it++) { // loop over timeslices


    if(TS>1){ // otherwise already set
      if(source!=5) fscanf(inputfilenamelist,"%s",INPUT_NAME); // otherwise 1 file
      if(dest!=5) fscanf(outputfilenamelist,"%s",OUTPUT_NAME); // otherwise 1 file
      fprintf(stderr,"At TS=%d of %d using input file %s and output file %s\n",it,TS,INPUT_NAME,OUTPUT_NAME);
    }

    fprintf(stderr,"Open source=%d\n",source);
    
    // open source file
    if(source==1){
      if( (input=fopen(INPUT_NAME,"rb"))==NULL){
	fprintf(stderr,"cannot open %s\n",INPUT_NAME);
	exit(1);
      }
    }
    if(source==2){
      if( (input=fopen(INPUT_NAME,"rt"))==NULL){
	fprintf(stderr,"cannot open %s\n",INPUT_NAME);
	exit(1);
      }
    }
#if(HDF)
    if(source==3){
      // open HDF file
      sd_id=SDstart(INPUT_NAME,DFACC_READ);
      if (sd_id != FAIL)      printf ("Reading HDF file with READ access\n");
      
      sd_index = 0;
      sds_id = SDselect (sd_id, sd_index);
      
      istat = SDreaddata (sds_id, start, NULL, edges, (VOIDP) array);
    }
    if(source==4){
      // not yet
    }
#endif
#if(V5D)
    if((source==5)&&(it==0)){ // assume all input vis5d are multi-timed
      // not yet
    }
#endif

    if(source==0){
      
      // length of the entire unmodified input file name
      lname=strlen(INPUT_NAME);
      if( (INPUT_NAME[lname-1]=='z')&&(INPUT_NAME[lname-2]=='g')&&(INPUT_NAME[lname-3]=='.') ){
	inputtype=1;
	printf("input flagged as gzip\n");
      }
      else{
	inputtype=0;
      }
      
      if(inputtype==0){
	if( !(input=fopen(INPUT_NAME,"rb"))){
	  fprintf(stderr,"trouble opening input file: %s\n",INPUT_NAME);
	  exit(1);
	}
      }
      if(inputtype==1){
	sprintf(incommand,"gzip -d < %s",INPUT_NAME);
	if( !(input=popen(incommand,"r"))){
	  fprintf(stderr,"trouble opening input file: %s %s\n",INPUT_NAME,incommand);
	  exit(1);
	}
      }
      // assume header info provided and # lines of header provided (unlike in r8toras.c and block2tile.c)
      
    }

    //    fprintf(stderr,"Deal with header: pipeheader=%d\n",pipeheader);


    if( !(output=fopen(OUTPUT_NAME,"wt"))){
      fprintf(stderr,"trouble opening output file: %s\n",OUTPUT_NAME);
      exit(1);
      }
    fclose(output);


    /////////////////////////
    //  
    // deal with header
    //
    if(pipeheader){
      if( !(output=fopen(OUTPUT_NAME,"wt"))){
	fprintf(stderr,"trouble opening output file: %s\n",OUTPUT_NAME);
	exit(1);
      }
    }

    for(i=0;i<realnumheader;i++){
      //printf("headerline: %i\n",i);
      while((ch=fgetc(input))!='\n'){
	if(pipeheader) fputc(ch,output);
	//printf("%c",ch);
      }
    }

    // close header part if opened
    if(pipeheader){
      fprintf(output,"\n");
      fclose(output);
    }


    fprintf(stderr,"Open dest=%d\n",dest);
    // open destination file
    if(dest==0){
      lname=strlen(OUTPUT_NAME);
      if( (OUTPUT_NAME[lname-1]=='z')&&(OUTPUT_NAME[lname-2]=='g')&&(OUTPUT_NAME[lname-3]=='.') ){
	outputtype=1;
	printf("output flagged as gzip\n");
      }
      else{
	outputtype=0;
      }
      
      // open output file
      if(outputtype==0){
	if( !(output=fopen(OUTPUT_NAME,"ab"))){
	  fprintf(stderr,"trouble opening output file: %s\n",OUTPUT_NAME);
	  exit(1);
	}
      }
      if(outputtype==1){
	sprintf(outcommand,"gzip > %s",OUTPUT_NAME);
	if( !(output=popen(outcommand,"w"))){
	  fprintf(stderr,"trouble opening output file: %s %s\n",OUTPUT_NAME,outcommand);
	  exit(1);
	}
      }
    }
    if(dest==1){
      if( (output=fopen(OUTPUT_NAME,"at"))==NULL){
	fprintf(stderr,"cannot open %s\n",OUTPUT_NAME);
	exit(1);
      }
    }
    if(dest==2){
      if( (output=fopen(OUTPUT_NAME,"at"))==NULL){
	fprintf(stderr,"cannot open %s\n",OUTPUT_NAME);
	exit(1);
      }
    }
#if(HDF)
    if(dest==3){
      // open HDF file
      sd_id = SDstart (OUTPUT_NAME, DFACC_CREATE);
      
      // can change to other output types
      
      sds_id = SDcreate (sd_id, SDS_NAME, HDFTYPE, rank, dim_sizes);
      start[0]=start[1]=start[2]=start[3]=0;
    }
    if(dest==4){
      // not yet
    }
#endif
#if(V5D)
    if((dest==5)&&(it==0)){ // assume all output vis5d are multi-timed, so only open 1 file for all timeslices
      fprintf(stderr,"Opening vis5d file: %s %d %d %d %d %d %s\n",OUTPUT_NAME,NumTimes,NumVars,Nr,Nc,Nl[0],VarName[0]); fflush(stderr);

      /* use the v5dCreate call to create the v5d file and write the header */
      if (!v5dCreate( OUTPUT_NAME, NumTimes, NumVars, Nr, Nc, Nl,
		      (const char (*)[MAXVARNAME]) VarName,
		      TimeStamp, DateStamp, CompressMode,
		      Projection, ProjArgs, Vertical, VertArgs )) {
	printf("Error: couldn't create %s\n", OUTPUT_NAME );
	exit(1);
      }
    }
#endif

    
    
    
    
    


  fprintf(stderr,"reading file and putting into other format...\n"); fflush(stderr);

  BUFFERINIT;// use to fill/read array if HDF involved
  for(i=0;i<numrows;i++){
#if(DEBUG)
    fprintf(stderr,"rownumber: %d of %d\n",i,numrows-1); fflush(stderr);
#endif
    
    for(j=0;j<numterms;j++){ // over rows and each group
#if(DEBUG)
    fprintf(stderr,"termnumber: %d of %d\n",j,numterms-1); fflush(stderr);
#endif
      
      for(k=0;k<group[j];k++){ // over a group of same kind
#if(DEBUG)
    fprintf(stderr,"groupelements: %d of %d : precision: %c s: %d d: %d\n",k,group[j]-1,precision[j],source,dest); fflush(stderr);
#endif
	switch(source){	 
	case 0:
	  fread(&dumb,bytesize,1,input); precision[j]='b'; //forced
	  break;
	case 1:
	  if(precision[j]=='b')	            fread(&dumb,bytesize,1,input);
	  if(precision[j]=='i')	            fread(&dumi,intsize,1,input);
	  if(precision[j]=='f')             fread(&dumf,floatsize,1,input);
	  if(precision[j]=='d')             fread(&dumlf,doublesize,1,input);
	  break;
	case 2:
	  if(precision[j]=='b'){
	    fscanf(input,"%hu",&shortfordumb);
	    dumb=shortfordumb; // convert short to byte
	  }
	  
	  if(precision[j]=='i')	            fscanf(input,"%d",&dumi);
	  if(precision[j]=='f')             fscanf(input,"%f",&dumf);
	  if(precision[j]=='d')             fscanf(input,"%lf",&dumlf);
	  break;
#if(HDF)
	case 3:
	case 4:
	  if(precision[j]=='b') dumb=arrayb[nextbuf++];
	  if(precision[j]=='i') dumi=arrayi[nextbuf++];
	  if(precision[j]=='f') dumf=arrayf[nextbuf++];
	  if(precision[j]=='d') dumlf=arrayd[nextbuf++];
	  break;
#endif
#if(V5D)
	case 5:
	  // source is always float
	  dumf=arrayvisf[nextbuf++]; precision[j]='f'; // forced
	  break;
#endif
	default:
	  break;
	}
	switch(dest){
	case 0:
	  // always byte size output
	  if(precision[j]=='b'){dumb=dumb;  fwrite(&dumb,bytesize,1,output);}
	  if(precision[j]=='i'){dumb=dumi;  fwrite(&dumb,bytesize,1,output);}
	  if(precision[j]=='f'){dumb=dumf;  fwrite(&dumb,bytesize,1,output);}
	  if(precision[j]=='d'){dumb=dumlf; fwrite(&dumb,bytesize,1,output);}
	  break;
	case 1:
	  if(precision[j]=='b')		  fwrite(&dumb,bytesize,1,output);
	  if(precision[j]=='i'){if(ble){ SWAP_LONG(dumi);}     fwrite(&dumi,intsize,1,output);}
	  if(precision[j]=='f'){if(ble){ SWAP_FLOAT(dumf);}   fwrite(&dumf,floatsize,1,output);}
	  if(precision[j]=='d'){if(ble){ SWAP_DOUBLE(dumlf);} fwrite(&dumlf,doublesize,1,output);}
	  break;
	case 2:
	  if(precision[j]=='b'){
	    if(source==0 || source==1) fprintf(output,"%d ",dumb);
	    else fprintf(output,"%c ",dumb);
	  }
	  if(precision[j]=='i'){if(ble){SWAP_LONG(dumi);}     fprintf(output,"%d ",dumi);}
	  if(precision[j]=='f'){if(ble){SWAP_FLOAT(dumf);}   fprintf(output,"%17.10g ",dumf);}
	  if(precision[j]=='d'){if(ble){SWAP_DOUBLE(dumlf);} fprintf(output,"%26.20g ",dumlf);}
	  break;
#if(HDF)
	case 3:
	case 4:
	  if(precision[j]=='b') arrayb[nextbuf++]=dumb;
	  if(precision[j]=='i'){if(ble){SWAP_LONG(dumi);}     arrayi[nextbuf++]=dumi;}
	  if(precision[j]=='f'){if(ble){SWAP_FLOAT(dumf);}   arrayf[nextbuf++]=dumf;}
	  if(precision[j]=='d'){if(ble){SWAP_DOUBLE(dumlf);} arrayd[nextbuf++]=dumlf;}
	  break;
#endif
#if(V5D)
	case 5:
	  // same array, must be float
	  if(precision[j]=='b') arrayvisf[nextbuf++]=dumb;
	  if(precision[j]=='i'){if(ble){SWAP_LONG(dumi);}      arrayvisf[nextbuf++]=dumi;}
	  if(precision[j]=='f'){if(ble){SWAP_FLOAT(dumf);}    arrayvisf[nextbuf++]=dumf;}
	  if(precision[j]=='d'){if(ble){SWAP_DOUBLE(dumlf);}  arrayvisf[nextbuf++]=dumlf;}
	  break;
#endif
	default:
	  break;
	}
      }
    }
    // skip terms till carraige return (allows parsing out certain extra columns on end)
    if(source==2){ while(fgetc(input)!='\n'); }
    // output carraige return
    if(dest==2) fprintf(output,"\n");
  }

#if(HDF)
  if(dest==3){
    fprintf(stderr,"Writing HDF file ...\n"); fflush(stderr);
    status = SDwritedata (sds_id, start, NULL, edges, (VOIDP)array); 
    status = SDendaccess (sds_id);
    status = SDend (sd_id);
  }
#endif
  //#if(V5D&&0) // DEBUG GODMARK
#if(V5D)
  if(dest==5){
    fprintf(stderr,"Writing vis5d file ...\n"); fflush(stderr);

    for(iv=0;iv<NumVars;iv++){

      /**
       ** Read your 3-D grid data for timestep it and variable
       ** iv into the array g here.
       ** To help with 3-D array indexing we've defined a macro G.
       ** G(0,0,0) is the north-west-bottom corner, G(Nr-1,Nc-1,Nl-1) is
       ** the south-east-top corner.  If you want a value to be considered
       ** missing, assign it equal to the constant MISSING.  For example:
       ** G(ir,ic,il) = MISSING;
       **/
      //#define G(ROW, COLUMN, LEVEL)   g[ (ROW) + ((COLUMN) + (LEVEL) * Nc) * Nr ]
      
      
      if(source==0){ // then use min/max conversion for decent legend (at least for fixed scaled data)
	// also assumes linear legend!
	a=minmax[0][iv];
	b=minmax[1][iv];
      }
      else{ // to cancel change
	a=0.0;
	b=255.0;
      }
      // convert my format to vis5d format
      for(k=0;k<N3;k++) for(j=0;j<N2;j++) for(i=0;i<N1;i++){
	ijkjon=(i+(j+k*N2)*N1)*NumVars+iv;
	ijkvis5d=(N3-1-k) + ((j) + (i) * N2) * N3;
	arrayvisoutput[ijkvis5d]=(arrayvisf[ijkjon]/255.0)*(b-a)+a;

	// DEBUG:
	//	fprintf(stderr,"i=%d j=%d k=%d iv=%d ijkvis5d=%d a=%g b=%g value=%g\n",i,j,k,iv,ijkvis5d,a,b,arrayvisoutput[ijkvis5d]);
      }
      
      /* Write data to v5d file. */
      if (!v5dWrite( it+1, iv+1, arrayvisoutput )) {
	printf("Error while writing grid.  Disk full?\n");
	exit(1);
      }
    }
  }
#endif
  if((dest==1)||(dest==2)) fclose(output);
  else if(dest==0){
    if(outputtype==0) fclose(output);
    else if(outputtype==1) pclose(output);
  }

#if(HDF)
  if(source==3){
    /* Terminate access to the array. */
    istat = SDendaccess(sds_id);
    
    /* Terminate access to the SD interface and close the file. */
    istat = SDend(sd_id);
    if (istat != FAIL) printf("... file closed\n\n");
  }
  else if(source==4){
    // not yet
  }
#endif
  if((source==1)||(source==2)) fclose(input);
  else if(source==0){
    if(inputtype==0) fclose(input);
    else if(inputtype==1) pclose(input);
  }

  // END BIG LOOP
  } // over all timesteps


#if(V5D)
  // close v5d file which has entire time series in it
  if(dest==5){
    v5dClose();
  }
  if(source==5){ // close entire time series v5d
    // not yet
  }
#endif
  if(TS>1){
    if(source!=5){ // otherwise 1 file
      fclose(inputfilenamelist);
    }
    if(dest!=5){ // otherwise 1 file
      fclose(outputfilenamelist);
    }
  }

  

  fprintf(stderr,"done.\n");
  return(0);
  //  exit(0); // not reachable
}


// apparently doesn't work:
void swap(BYTE *x, BYTE size)
{
  unsigned char c;
  unsigned short s;
  unsigned long l;

  switch (size)
  {
      case 1: /* don't do anything */
	  break;
      case 2: /* swap two bytes */
	  c = *x;
	  *x = *(x+1);
	  *(x+1) = c;
	  break;
      case 4: /* swap two shorts (2-byte words) */
	  s = *(unsigned short *)x;
	  *(unsigned short *)x = *((unsigned short *)x + 1);
	  *((unsigned short *)x + 1) = s;
	  swap ((BYTE *)x, 2);
	  swap ((BYTE *)((unsigned short *)x+1), 2);
	  break;
      case 8: /* swap two longs (4-bytes words) */
	  l = *(unsigned long *)x;
	  *(unsigned long *)x = *((unsigned long *)x + 1);
	  *((unsigned long *)x + 1) = l;
	  swap ((BYTE *)x, 4);
	  swap ((BYTE *)((unsigned long *)x+1), 4);
	  break;
  }
}


int machineEndianness(void)
{
   int i = 1;
   char *p = (char *) &i;
   if (p[0] == 1){
     // Lowest address contains the least significant byte
     //     fprintf(stderr,"little\n");
     return LITTLE_ENDIAN_USER;
   }
   else{
     //     fprintf(stderr,"big\n");
     return BIG_ENDIAN_USER;
   }

   //   return(0); // not reachable
}



//#include "SwapEndian.h"



static long _TestEndian=1;

int IsLittleEndian(void) {
	return *(char*)&_TestEndian;
}

/******************************************************************************
  FUNCTION: SwapEndian
  PURPOSE: Swap the byte order of a structure
  EXAMPLE: float F=123.456;; SWAP_FLOAT(F);
******************************************************************************/

void *SwapEndian(void* Addr, const int Nb) {
	static char Swapped[16];
	switch (Nb) {
		case 2:	Swapped[0]=*((char*)Addr+1);
				Swapped[1]=*((char*)Addr  );
				break;
		case 3:	// As far as I know, 3 is used only with RGB images
				Swapped[0]=*((char*)Addr+2);
				Swapped[1]=*((char*)Addr+1);
				Swapped[2]=*((char*)Addr  );
				break;
		case 4:	Swapped[0]=*((char*)Addr+3);
				Swapped[1]=*((char*)Addr+2);
				Swapped[2]=*((char*)Addr+1);
				Swapped[3]=*((char*)Addr  );
				break;
		case 8:	Swapped[0]=*((char*)Addr+7);
				Swapped[1]=*((char*)Addr+6);
				Swapped[2]=*((char*)Addr+5);
				Swapped[3]=*((char*)Addr+4);
				Swapped[4]=*((char*)Addr+3);
				Swapped[5]=*((char*)Addr+2);
				Swapped[6]=*((char*)Addr+1);
				Swapped[7]=*((char*)Addr  );
				break;
		case 16:Swapped[0]=*((char*)Addr+15);
				Swapped[1]=*((char*)Addr+14);
				Swapped[2]=*((char*)Addr+13);
				Swapped[3]=*((char*)Addr+12);
				Swapped[4]=*((char*)Addr+11);
				Swapped[5]=*((char*)Addr+10);
				Swapped[6]=*((char*)Addr+9);
				Swapped[7]=*((char*)Addr+8);
				Swapped[8]=*((char*)Addr+7);
				Swapped[9]=*((char*)Addr+6);
				Swapped[10]=*((char*)Addr+5);
				Swapped[11]=*((char*)Addr+4);
				Swapped[12]=*((char*)Addr+3);
				Swapped[13]=*((char*)Addr+2);
				Swapped[14]=*((char*)Addr+1);
				Swapped[15]=*((char*)Addr  );
				break;
	}
	return (void*)Swapped;
}
