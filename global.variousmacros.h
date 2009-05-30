


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Various macros
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define MYDMIN(a,b) (mydminarg1=(a),mydminarg2=(b),(mydminarg1) < (mydminarg2) ?\
        (mydminarg1) : (mydminarg2))

#define delta(i,j) ((i == j) ? 1. : 0.)
#define dot(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3])

#define mink(I,J) (I != J ? (0.) : (I == 0 ? (-1.) : (1.)))

#define pfixupeach(pr,i,j,k,which,min) {if(pr[which]<min){ fladd[which]+=dV*MYGDET(i,j,k,CENT)*(min-pr[which]); pr[which]=min;}}

#define pfixup(pr,i,j,k) {pfixupeach(pr,i,j,k,RHO,RHOMIN); pfixupeach(pr,i,j,k,UU,UUMIN); }

// #define FAILSTATEMENT(file,function,number) {fprintf(fail_file,"%s
// %d-%s(): failure\n",file,number,function); fflush(fail_file);
// fprintf(fail_file,"rho[i][j][k]: %21.15g uu[i][j][k]: %21.15g rho2[i][j][k]:
// %21.15g uu2[i][j][k]: %21.15g i: %d j: %d pl:
// %d\n",p[i][j][k][RHO],p[i][j][k][UU],ph[i][j][k][RHO],ph[i][j][k][UU],i,j,pl);
// return(1);}

#define FAILSTATEMENT(file,function,number) {if(debugfail>=1){ dualfprintf(fail_file,"%s %d-%s(): failure\n",file,number,function); dualfprintf(fail_file,"i: %d j: %d k: %d p: %d\n",icurr,jcurr,kcurr,pcurr);} return(1);}

#define FAILSTATEMENTVOID(file,function,number) {if(debugfail>=1){ dualfprintf(fail_file,"%s %d-%s(): failure\n",file,number,function); dualfprintf(fail_file,"i: %d j: %d k: %d p: %d\n",icurr,jcurr,kcurr,pcurr);} }


#if(JONCHECKS2 && PRODUCTION==0)
#define MYFUN(fun,one,two,three) if(fun>=1){ FAILSTATEMENT(one,two,three);}
#else
// if PRODUCTION==1 then avoid if statement
#define MYFUN(fun,one,two,three) {fun;}
#endif
