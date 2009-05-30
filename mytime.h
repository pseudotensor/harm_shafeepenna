#ifndef _MYTIME_H
#define _MYTIME_H

#include <signal.h>
#include <time.h>

#ifndef WIN32
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <sys/times.h>
#else
#include <Winsock2.h>
#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif
struct timezone
{
  int  tz_minuteswest; /* minutes W of Greenwich */
  int  tz_dsttime;     /* type of dst correction */
};

extern int gettimeofday(struct timeval *tv, struct timezone *tz);

#endif

void mycpuclock(clock_t *time);
void myustimes(clock_t *time);
void myustimes2(clock_t *usertime,clock_t *systime);


#ifndef WIN32

#if(PERFTEST==1)
#define TIMEMETHOD 2 // SUPERMARK
#else
#define TIMEMETHOD 1 // used for measurements in wall time
#endif
// 0: old second accurate method for walltime
// 1: new microsecond accurate method for walltime
// 2: cpu(user+system) time, accurate to 1/CLOCKS_PER_SECOND seconds, not wall time
// 3: reports user/system and child user/system times (good for understanding if system is eating lot of your time (i.e. HD or net or whatever hardware device run by kernel))

#else
// can choose TIMEMETHOD == 0,1,2
#define TIMEMETHOD 1
#endif



// average time to completion is: T= ((zc/s)^{-1}) * (#zones) * (dt/tf) .  dt~dx/v -> T~dx^3 for 2D dx^4 for 3D



#define SEC2HOUR (2.77777777777777E-4)
#define microtime(time) gettimeofday(time,&tz)
#define diffmicrotime(timestop,timestart) ((SFTYPE)(timestop.tv_sec-timestart.tv_sec)+(SFTYPE)(timestop.tv_usec-timestart.tv_usec)*1E-6)
// can use clock() if time is < 71.5827882667 minutes for 32-bit systems.
// *time is used since sending pointer in general
#define TOTALTICKS (4294967296) // 32-bit system
//#define TOTALTICKS (18446744073709551616) // 64-bit system  (compiler complains even though shouldn't?!)
//#endif

#define cpuclock(time) mycpuclock(time)
#if(NCSA==0)
#define diffcpuclock(timestop,timestart) ((timestop>timestart) ? (SFTYPE)(timestop-timestart)/(SFTYPE)(CLOCKS_PER_SEC) : (SFTYPE)(timestop-(timestart-TOTALTICKS))/(SFTYPE)(CLOCKS_PER_SEC) )
#else
#define diffcpuclock(timestop,timestart)  ((SFTYPE)(timestop-timestart)/(SFTYPE)(CLOCKS_PER_SEC))  // can't trust above since not taking 64bit totaltick value
#endif
#define diffmyustimes(timestop,timestart) ((SFTYPE)(timestop-timestart)/1000000.0)

// can use nanosleep() or clock_nanosleep() to pause for 

#if(TIMEMETHOD==0)
#define GETTIME time
#define DELTATIME difftime
#elif(TIMEMETHOD==1)
#define GETTIME microtime
#define DELTATIME diffmicrotime
#elif(TIMEMETHOD==2)
#define GETTIME cpuclock
#define DELTATIME diffcpuclock
#elif(TIMEMETHOD==3)
#define GETTIME myustimes
#define DELTATIME diffmyustimes
#endif




#endif // end ifndef _MYTIME_H


