#ifndef _lu_output_darts_h_
#define _lu_output_darts_h_
#ifndef __DARTS_
#define __DARTS_
#endif
#include "npb-C.h"
#include "applu.h"
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <limits.h>
#include "utils.h"
#include <numa.h>
#include "darts.h"
#include "ompTP.h"
#include "tbb/concurrent_vector.h"
#include <mutex>
int main(int argc, char **argv) ;
class TP192;
class TP1;
typedef TP1 TP_blts;
class TP248;
/*Number of TPs to be used for the OMPFor in region TP248*/
#define NUMTPS248 NUMTPS
class TP2;
typedef TP2 TP_buts;
class TP1218;
/*Number of TPs to be used for the OMPFor in region TP1218*/
#define NUMTPS1218 NUMTPS
class TP2204;
class TP2223;
/*Number of TPs to be used for the OMPFor in region TP2223*/
#define NUMTPS2223 NUMTPS
class TP2274;
/*Number of TPs to be used for the OMPFor in region TP2274*/
#define NUMTPS2274 NUMTPS
class TP2409;
/*Number of TPs to be used for the OMPFor in region TP2409*/
#define NUMTPS2409 NUMTPS
class TP2558;
/*Number of TPs to be used for the OMPFor in region TP2558*/
#define NUMTPS2558 NUMTPS
class TP3196;
/*Number of TPs to be used for the OMPFor in region TP3196*/
#define NUMTPS3196 NUMTPS
class TP3345;
/*Number of TPs to be used for the OMPFor in region TP3345*/
#define NUMTPS3345 NUMTPS
class TP3980;
/*Number of TPs to be used for the OMPFor in region TP3980*/
#define NUMTPS3980 NUMTPS
class TP7;
typedef TP7 TP_jacld;
class TP4906;
/*Number of TPs to be used for the OMPFor in region TP4906*/
#define NUMTPS4906 NUMTPS
class TP8;
typedef TP8 TP_jacu;
class TP7428;
/*Number of TPs to be used for the OMPFor in region TP7428*/
#define NUMTPS7428 NUMTPS
class TP9879;
class TP9897;
/*Number of TPs to be used for the OMPFor in region TP9897*/
#define NUMTPS9897 NUMTPS
class TP10796;
class TP10811;
/*Number of TPs to be used for the OMPFor in region TP10811*/
#define NUMTPS10811 NUMTPS
class TP10871;
/*Number of TPs to be used for the OMPFor in region TP10871*/
#define NUMTPS10871 NUMTPS
class TP11020;
/*Number of TPs to be used for the OMPFor in region TP11020*/
#define NUMTPS11020 NUMTPS
class TP11660;
/*Number of TPs to be used for the OMPFor in region TP11660*/
#define NUMTPS11660 NUMTPS
class TP11809;
/*Number of TPs to be used for the OMPFor in region TP11809*/
#define NUMTPS11809 NUMTPS
class TP12446;
/*Number of TPs to be used for the OMPFor in region TP12446*/
#define NUMTPS12446 NUMTPS
class TP13197;
class TP13201;
/*Number of TPs to be used for the OMPFor in region TP13201*/
#define NUMTPS13201 NUMTPS
class TP13252;
/*Number of TPs to be used for the OMPFor in region TP13252*/
#define NUMTPS13252 NUMTPS
class TP13294;
/*Number of TPs to be used for the OMPFor in region TP13294*/
#define NUMTPS13294 NUMTPS
class TP13338;
/*Number of TPs to be used for the OMPFor in region TP13338*/
#define NUMTPS13338 NUMTPS
class TP13380;
/*Number of TPs to be used for the OMPFor in region TP13380*/
#define NUMTPS13380 NUMTPS
class TP13764;
class TP13770;
/*Number of TPs to be used for the OMPFor in region TP13770*/
#define NUMTPS13770 NUMTPS
class TP13896;
class TP13898;
/*Number of TPs to be used for the OMPFor in region TP13898*/
#define NUMTPS13898 NUMTPS
class TP13995;
class TP13997;
/*Number of TPs to be used for the OMPFor in region TP13997*/
#define NUMTPS13997 NUMTPS
class TP14054;
class TP14069;
class TP14080;
/*Number of TPs to be used for the OMPFor in region TP14080*/
#define NUMTPS14080 NUMTPS
extern int DARTS_CODELETS_MULT;
extern int NUMTPS;
extern size_t numOfCUs;
extern darts::Codelet* RuntimeFinalCodelet;
extern darts::ThreadAffinity *affin;
extern bool affinMaskRes;
extern darts::Runtime *myDARTSRuntime;
extern std::vector< std::vector < void* > > threadFunctionStack;
extern size_t ompNumThreads;
extern int ompSchedulePolicy;
extern int ompScheduleChunk;
extern void omp_set_num_threads(unsigned long numThreadsToSet);
extern int omp_get_num_threads();
extern int omp_get_max_threads();
extern int omp_get_num_procs();
extern double omp_get_wtime();
extern void omp_init_lock(omp_lock_t * lock);
extern void omp_destroy_lock(omp_lock_t * lock);
extern void omp_set_lock(omp_lock_t * lock);
extern void omp_unset_lock(omp_lock_t * lock);
/*TP192: OMPParallelDirective*/
class TP192 : public darts::ThreadedProcedure
{
 public:
class _barrierCodelets192 : public darts::Codelet
{
public:
TP192* inputsTPParent;
_barrierCodelets192():
darts::Codelet(){ }
_barrierCodelets192(uint32_t dep, uint32_t res, TP192* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets194:public darts::Codelet
{
public:
TP192* myTP;
TP192* inputsTPParent;
_checkInCodelets194():
darts::Codelet(){ }
_checkInCodelets194(uint32_t dep, uint32_t res, TP192* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
darts::Codelet* nextCodelet;
TP192 * TPParent;
TP192 * controlTPParent;
TP192* inputsTPParent;int *nthreads_darts192;/*OMP_SHARED - INPUT*/
int *nthreads_darts194;/*OMP_SHARED - INPUT*/
size_t TP194_alreadyLaunched;
_barrierCodelets192* barrierCodelets192;
_checkInCodelets194* checkInCodelets194;
TP192(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet, int *in_nthreads);
~TP192();
};
/*TP1: blts*/
class TP1:public ompTP
{
 public:
class _checkInCodelets248:public darts::Codelet
{
public:
TP1* myTP;
TP1* inputsTPParent;
_checkInCodelets248():
darts::Codelet(){ }
_checkInCodelets248(uint32_t dep, uint32_t res, TP1* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets248 : public darts::Codelet
{
public:
TP1* inputsTPParent;
_barrierCodelets248():
darts::Codelet(){ }
_barrierCodelets248(uint32_t dep, uint32_t res, TP1* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets352:public darts::Codelet
{
public:
TP1* myTP;
TP1* inputsTPParent;
_checkInCodelets352():
darts::Codelet(){ }
_checkInCodelets352(uint32_t dep, uint32_t res, TP1* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets354:public darts::Codelet
{
public:
TP1* myTP;
TP1* inputsTPParent;
_checkInCodelets354():
darts::Codelet(){ }
_checkInCodelets354(uint32_t dep, uint32_t res, TP1* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP1** ptrToThisFunctionTP;
TP1* inputsTPParent;
TP1* controlTPParent;
darts::Codelet** nextCodeletsblts;
darts::Codelet** nextSyncCodeletsblts;
int* nx_darts1/*VARIABLE*/;
int* ny_darts1/*VARIABLE*/;
int* nz_darts1/*VARIABLE*/;
int* k_darts1/*VARIABLE*/;
double* omega_darts1/*VARIABLE*/;
int* ist_darts1/*VARIABLE*/;
int* iend_darts1/*VARIABLE*/;
int* jst_darts1/*VARIABLE*/;
int* jend_darts1/*VARIABLE*/;
int* nx0_darts1/*VARIABLE*/;
int* ny0_darts1/*VARIABLE*/;
int* i_darts1/*VARIABLE*/;
int* j_darts1/*VARIABLE*/;
int* m_darts1/*VARIABLE*/;
double* tmp_darts1/*VARIABLE*/;
double* tmp1_darts1/*VARIABLE*/;
int i_darts354 ;
int* iend_darts354/*OMP_SHARED_PRIVATE - INPUT*/;
int* ist_darts354/*OMP_SHARED_PRIVATE - INPUT*/;
int j_darts354 ;
int* jend_darts354/*OMP_SHARED_PRIVATE - INPUT*/;
int* jst_darts354/*OMP_SHARED_PRIVATE - INPUT*/;
int* k_darts354/*OMP_SHARED_PRIVATE - INPUT*/;
int m_darts354 ;
double* omega_darts354/*OMP_SHARED_PRIVATE - INPUT*/;
double* tmp_darts354/*OMP_SHARED_PRIVATE - INPUT*/;
double* tmp1_darts354/*OMP_SHARED_PRIVATE - INPUT*/;
int id_darts354 /*VARIABLE*/;
TP248** TP248Ptr;
size_t *TP248_alreadyLaunched;
int numTPsSet248;
int numTPsReady248;
size_t TPsToUse248;
size_t codeletsPerTP248;
size_t totalCodelets248;
size_t TP352_alreadyLaunched;
size_t TP354_alreadyLaunched;
_checkInCodelets248* checkInCodelets248;
#if USE_SPIN_CODELETS == 0
_checkInCodelets248* firstCodelet;
#endif
_barrierCodelets248* barrierCodelets248;
_checkInCodelets352* checkInCodelets352;
_checkInCodelets354* checkInCodelets354;
TP1(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet, darts::Codelet* in_mainSyncCodelet, TP1** in_ptrToThisFunctionTP, int in_nx, int in_ny, int in_nz, int in_k, double in_omega, int in_ist, int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0);
~TP1();
void setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, int in_ist, int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID);};
/*TP248: OMPForDirective*/
class TP248 : public ompTP
{
 public:
class _barrierCodelets248 : public darts::Codelet
{
public:
TP248* inputsTPParent;
_barrierCodelets248():
darts::Codelet(){ }
_barrierCodelets248(uint32_t dep, uint32_t res, TP248* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations248(int* endRange, uint32_t codeletID);
class _checkInCodelets249:public darts::Codelet
{
public:
TP248* myTP;
TP248* inputsTPParent;
int endRange;
_checkInCodelets249():
darts::Codelet(){ }
_checkInCodelets249(uint32_t dep, uint32_t res, TP248* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP1 * TPParent;
TP248 * controlTPParent;
TP248* inputsTPParent;int *i_darts248/*OMP_PRIVATE - INPUT*/;
int **iend_darts248 /*OMP_SHARED_PRIVATE - INPUT*/;
int **ist_darts248 /*OMP_SHARED_PRIVATE - INPUT*/;
int *j_darts248/*OMP_PRIVATE - INPUT*/;
int **jend_darts248 /*OMP_SHARED_PRIVATE - INPUT*/;
int **jst_darts248 /*OMP_SHARED_PRIVATE - INPUT*/;
int **k_darts248 /*OMP_SHARED_PRIVATE - INPUT*/;
int *m_darts248/*OMP_PRIVATE - INPUT*/;
double **omega_darts248 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration248;
int lastIteration248;
int range248;
int rangePerCodelet248;
int minIteration248;
int remainderRange248;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets248* barrierCodelets248;
_checkInCodelets249* checkInCodelets249;
#if USE_SPIN_CODELETS == 0
_checkInCodelets249* firstCodelet;
#endif
TP248(int in_numThreads, int in_mainCodeletID, TP1* in_TPParent, int in_initIteration, int in_lastIteration, TP248** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP248();
};
/*TP2: buts*/
class TP2:public ompTP
{
 public:
class _checkInCodelets1218:public darts::Codelet
{
public:
TP2* myTP;
TP2* inputsTPParent;
_checkInCodelets1218():
darts::Codelet(){ }
_checkInCodelets1218(uint32_t dep, uint32_t res, TP2* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets1320:public darts::Codelet
{
public:
TP2* myTP;
TP2* inputsTPParent;
_checkInCodelets1320():
darts::Codelet(){ }
_checkInCodelets1320(uint32_t dep, uint32_t res, TP2* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP2** ptrToThisFunctionTP;
TP2* inputsTPParent;
TP2* controlTPParent;
darts::Codelet** nextCodeletsbuts;
darts::Codelet** nextSyncCodeletsbuts;
int* nx_darts2/*VARIABLE*/;
int* ny_darts2/*VARIABLE*/;
int* nz_darts2/*VARIABLE*/;
int* k_darts2/*VARIABLE*/;
double* omega_darts2/*VARIABLE*/;
int* ist_darts2/*VARIABLE*/;
int* iend_darts2/*VARIABLE*/;
int* jst_darts2/*VARIABLE*/;
int* jend_darts2/*VARIABLE*/;
int* nx0_darts2/*VARIABLE*/;
int* ny0_darts2/*VARIABLE*/;
int* i_darts2/*VARIABLE*/;
int* j_darts2/*VARIABLE*/;
int* m_darts2/*VARIABLE*/;
double* tmp_darts2/*VARIABLE*/;
double* tmp1_darts2/*VARIABLE*/;
int i_darts1320 ;
int* iend_darts1320/*OMP_SHARED_PRIVATE - INPUT*/;
int* ist_darts1320/*OMP_SHARED_PRIVATE - INPUT*/;
int j_darts1320 ;
int* jend_darts1320/*OMP_SHARED_PRIVATE - INPUT*/;
int* jst_darts1320/*OMP_SHARED_PRIVATE - INPUT*/;
int* k_darts1320/*OMP_SHARED_PRIVATE - INPUT*/;
int m_darts1320 ;
double* omega_darts1320/*OMP_SHARED_PRIVATE - INPUT*/;
double* tmp_darts1320/*OMP_SHARED_PRIVATE - INPUT*/;
double* tmp1_darts1320/*OMP_SHARED_PRIVATE - INPUT*/;
TP1218** TP1218Ptr;
size_t *TP1218_alreadyLaunched;
int numTPsSet1218;
int numTPsReady1218;
size_t TPsToUse1218;
size_t codeletsPerTP1218;
size_t totalCodelets1218;
size_t TP1320_alreadyLaunched;
_checkInCodelets1218* checkInCodelets1218;
#if USE_SPIN_CODELETS == 0
_checkInCodelets1218* firstCodelet;
#endif
_checkInCodelets1320* checkInCodelets1320;
TP2(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet, darts::Codelet* in_mainSyncCodelet, TP2** in_ptrToThisFunctionTP, int in_nx, int in_ny, int in_nz, int in_k, double in_omega, int in_ist, int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0);
~TP2();
void setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, int in_ist, int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID);};
/*TP1218: OMPForDirective*/
class TP1218 : public ompTP
{
 public:
bool requestNewRangeIterations1218(int* endRange, uint32_t codeletID);
class _checkInCodelets1219:public darts::Codelet
{
public:
TP1218* myTP;
TP1218* inputsTPParent;
int endRange;
_checkInCodelets1219():
darts::Codelet(){ }
_checkInCodelets1219(uint32_t dep, uint32_t res, TP1218* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP2 * TPParent;
TP1218 * controlTPParent;
TP1218* inputsTPParent;int *i_darts1218/*OMP_PRIVATE - INPUT*/;
int **iend_darts1218 /*OMP_SHARED_PRIVATE - INPUT*/;
int **ist_darts1218 /*OMP_SHARED_PRIVATE - INPUT*/;
int *j_darts1218/*OMP_PRIVATE - INPUT*/;
int **jend_darts1218 /*OMP_SHARED_PRIVATE - INPUT*/;
int **jst_darts1218 /*OMP_SHARED_PRIVATE - INPUT*/;
int **k_darts1218 /*OMP_SHARED_PRIVATE - INPUT*/;
int *m_darts1218/*OMP_PRIVATE - INPUT*/;
double **omega_darts1218 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration1218;
int lastIteration1218;
int range1218;
int rangePerCodelet1218;
int minIteration1218;
int remainderRange1218;
size_t readyCodelets;
int baseNumThreads;
int *signalNextReady;
_checkInCodelets1219* checkInCodelets1219;
#if USE_SPIN_CODELETS == 0
_checkInCodelets1219* firstCodelet;
#endif
TP1218(int in_numThreads, int in_mainCodeletID, TP2* in_TPParent, int in_initIteration, int in_lastIteration, TP1218** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP1218();
};
/*TP2204: OMPParallelDirective*/
class TP2204 : public darts::ThreadedProcedure
{
 public:
class _barrierCodelets2204 : public darts::Codelet
{
public:
TP2204* inputsTPParent;
_barrierCodelets2204():
darts::Codelet(){ }
_barrierCodelets2204(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets2206:public darts::Codelet
{
public:
TP2204* myTP;
TP2204* inputsTPParent;
_checkInCodelets2206():
darts::Codelet(){ }
_checkInCodelets2206(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets2223:public darts::Codelet
{
public:
TP2204* myTP;
TP2204* inputsTPParent;
_checkInCodelets2223():
darts::Codelet(){ }
_checkInCodelets2223(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets2223 : public darts::Codelet
{
public:
TP2204* inputsTPParent;
_barrierCodelets2223():
darts::Codelet(){ }
_barrierCodelets2223(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets2274:public darts::Codelet
{
public:
TP2204* myTP;
TP2204* inputsTPParent;
_checkInCodelets2274():
darts::Codelet(){ }
_checkInCodelets2274(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets2274 : public darts::Codelet
{
public:
TP2204* inputsTPParent;
_barrierCodelets2274():
darts::Codelet(){ }
_barrierCodelets2274(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets2406:public darts::Codelet
{
public:
TP2204* myTP;
TP2204* inputsTPParent;
_checkInCodelets2406():
darts::Codelet(){ }
_checkInCodelets2406(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets2409:public darts::Codelet
{
public:
TP2204* myTP;
TP2204* inputsTPParent;
_checkInCodelets2409():
darts::Codelet(){ }
_checkInCodelets2409(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets2409 : public darts::Codelet
{
public:
TP2204* inputsTPParent;
_barrierCodelets2409():
darts::Codelet(){ }
_barrierCodelets2409(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets2558:public darts::Codelet
{
public:
TP2204* myTP;
TP2204* inputsTPParent;
_checkInCodelets2558():
darts::Codelet(){ }
_checkInCodelets2558(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets2558 : public darts::Codelet
{
public:
TP2204* inputsTPParent;
_barrierCodelets2558():
darts::Codelet(){ }
_barrierCodelets2558(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets3193:public darts::Codelet
{
public:
TP2204* myTP;
TP2204* inputsTPParent;
_checkInCodelets3193():
darts::Codelet(){ }
_checkInCodelets3193(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets3196:public darts::Codelet
{
public:
TP2204* myTP;
TP2204* inputsTPParent;
_checkInCodelets3196():
darts::Codelet(){ }
_checkInCodelets3196(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets3196 : public darts::Codelet
{
public:
TP2204* inputsTPParent;
_barrierCodelets3196():
darts::Codelet(){ }
_barrierCodelets3196(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets3345:public darts::Codelet
{
public:
TP2204* myTP;
TP2204* inputsTPParent;
_checkInCodelets3345():
darts::Codelet(){ }
_checkInCodelets3345(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets3345 : public darts::Codelet
{
public:
TP2204* inputsTPParent;
_barrierCodelets3345():
darts::Codelet(){ }
_barrierCodelets3345(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets3980:public darts::Codelet
{
public:
TP2204* myTP;
TP2204* inputsTPParent;
_checkInCodelets3980():
darts::Codelet(){ }
_checkInCodelets3980(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets3980 : public darts::Codelet
{
public:
TP2204* inputsTPParent;
_barrierCodelets3980():
darts::Codelet(){ }
_barrierCodelets3980(uint32_t dep, uint32_t res, TP2204* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
darts::Codelet* nextCodelet;
TP2204 * TPParent;
TP2204 * controlTPParent;
TP2204* inputsTPParent;int *L1_darts2204/*VARIABLE*/;
int *L2_darts2204/*VARIABLE*/;
double *dsspm_darts2204/*VARIABLE*/;
double *eta_darts2204/*VARIABLE*/;
int *i_darts2204/*VARIABLE*/;
int *iend1_darts2204/*VARIABLE*/;
int *iglob_darts2204/*VARIABLE*/;
int *ist1_darts2204/*VARIABLE*/;
int *j_darts2204/*VARIABLE*/;
int *jend1_darts2204/*VARIABLE*/;
int *jglob_darts2204/*VARIABLE*/;
int *jst1_darts2204/*VARIABLE*/;
int *k_darts2204/*VARIABLE*/;
int *m_darts2204/*VARIABLE*/;
double *q_darts2204/*VARIABLE*/;
double *tmp_darts2204/*VARIABLE*/;
double *u21_darts2204/*VARIABLE*/;
double *u21i_darts2204/*VARIABLE*/;
double *u21im1_darts2204/*VARIABLE*/;
double *u21j_darts2204/*VARIABLE*/;
double *u21jm1_darts2204/*VARIABLE*/;
double *u21k_darts2204/*VARIABLE*/;
double *u21km1_darts2204/*VARIABLE*/;
double *u31_darts2204/*VARIABLE*/;
double *u31i_darts2204/*VARIABLE*/;
double *u31im1_darts2204/*VARIABLE*/;
double *u31j_darts2204/*VARIABLE*/;
double *u31jm1_darts2204/*VARIABLE*/;
double *u31k_darts2204/*VARIABLE*/;
double *u31km1_darts2204/*VARIABLE*/;
double *u41_darts2204/*VARIABLE*/;
double *u41i_darts2204/*VARIABLE*/;
double *u41im1_darts2204/*VARIABLE*/;
double *u41j_darts2204/*VARIABLE*/;
double *u41jm1_darts2204/*VARIABLE*/;
double *u41k_darts2204/*VARIABLE*/;
double *u41km1_darts2204/*VARIABLE*/;
double *u51i_darts2204/*VARIABLE*/;
double *u51im1_darts2204/*VARIABLE*/;
double *u51j_darts2204/*VARIABLE*/;
double *u51jm1_darts2204/*VARIABLE*/;
double *u51k_darts2204/*VARIABLE*/;
double *u51km1_darts2204/*VARIABLE*/;
double *xi_darts2204/*VARIABLE*/;
double *zeta_darts2204/*VARIABLE*/;
TP2223** TP2223Ptr;
size_t *TP2223_alreadyLaunched;
int numTPsSet2223;
int numTPsReady2223;
size_t TPsToUse2223;
size_t codeletsPerTP2223;
size_t totalCodelets2223;
TP2274** TP2274Ptr;
size_t *TP2274_alreadyLaunched;
int numTPsSet2274;
int numTPsReady2274;
size_t TPsToUse2274;
size_t codeletsPerTP2274;
size_t totalCodelets2274;
TP2409** TP2409Ptr;
size_t *TP2409_alreadyLaunched;
int numTPsSet2409;
int numTPsReady2409;
size_t TPsToUse2409;
size_t codeletsPerTP2409;
size_t totalCodelets2409;
TP2558** TP2558Ptr;
size_t *TP2558_alreadyLaunched;
int numTPsSet2558;
int numTPsReady2558;
size_t TPsToUse2558;
size_t codeletsPerTP2558;
size_t totalCodelets2558;
TP3196** TP3196Ptr;
size_t *TP3196_alreadyLaunched;
int numTPsSet3196;
int numTPsReady3196;
size_t TPsToUse3196;
size_t codeletsPerTP3196;
size_t totalCodelets3196;
TP3345** TP3345Ptr;
size_t *TP3345_alreadyLaunched;
int numTPsSet3345;
int numTPsReady3345;
size_t TPsToUse3345;
size_t codeletsPerTP3345;
size_t totalCodelets3345;
TP3980** TP3980Ptr;
size_t *TP3980_alreadyLaunched;
int numTPsSet3980;
int numTPsReady3980;
size_t TPsToUse3980;
size_t codeletsPerTP3980;
size_t totalCodelets3980;
_barrierCodelets2204* barrierCodelets2204;
_checkInCodelets2206* checkInCodelets2206;
_checkInCodelets2223* checkInCodelets2223;
_barrierCodelets2223* barrierCodelets2223;
_checkInCodelets2274* checkInCodelets2274;
_barrierCodelets2274* barrierCodelets2274;
_checkInCodelets2406* checkInCodelets2406;
_checkInCodelets2409* checkInCodelets2409;
_barrierCodelets2409* barrierCodelets2409;
_checkInCodelets2558* checkInCodelets2558;
_barrierCodelets2558* barrierCodelets2558;
_checkInCodelets3193* checkInCodelets3193;
_checkInCodelets3196* checkInCodelets3196;
_barrierCodelets3196* barrierCodelets3196;
_checkInCodelets3345* checkInCodelets3345;
_barrierCodelets3345* barrierCodelets3345;
_checkInCodelets3980* checkInCodelets3980;
_barrierCodelets3980* barrierCodelets3980;
TP2204(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
~TP2204();
};
/*TP2223: OMPForDirective*/
class TP2223 : public ompTP
{
 public:
class _barrierCodelets2223 : public darts::Codelet
{
public:
TP2223* inputsTPParent;
_barrierCodelets2223():
darts::Codelet(){ }
_barrierCodelets2223(uint32_t dep, uint32_t res, TP2223* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations2223(int* endRange, uint32_t codeletID);
class _checkInCodelets2224:public darts::Codelet
{
public:
TP2223* myTP;
TP2223* inputsTPParent;
int endRange;
_checkInCodelets2224():
darts::Codelet(){ }
_checkInCodelets2224(uint32_t dep, uint32_t res, TP2223* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP2204 * TPParent;
TP2223 * controlTPParent;
TP2223* inputsTPParent;int *i_darts2223/*OMP_PRIVATE - INPUT*/;
int *j_darts2223/*OMP_PRIVATE - INPUT*/;
int *k_darts2223/*OMP_PRIVATE - INPUT*/;
int *m_darts2223/*OMP_PRIVATE - INPUT*/;
int initIteration2223;
int lastIteration2223;
int range2223;
int rangePerCodelet2223;
int minIteration2223;
int remainderRange2223;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets2223* barrierCodelets2223;
_checkInCodelets2224* checkInCodelets2224;
#if USE_SPIN_CODELETS == 0
_checkInCodelets2224* firstCodelet;
#endif
TP2223(int in_numThreads, int in_mainCodeletID, TP2204* in_TPParent, int in_initIteration, int in_lastIteration, TP2223** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP2223();
};
/*TP2274: OMPForDirective*/
class TP2274 : public ompTP
{
 public:
class _barrierCodelets2274 : public darts::Codelet
{
public:
TP2274* inputsTPParent;
_barrierCodelets2274():
darts::Codelet(){ }
_barrierCodelets2274(uint32_t dep, uint32_t res, TP2274* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations2274(int* endRange, uint32_t codeletID);
class _checkInCodelets2275:public darts::Codelet
{
public:
TP2274* myTP;
TP2274* inputsTPParent;
int endRange;
_checkInCodelets2275():
darts::Codelet(){ }
_checkInCodelets2275(uint32_t dep, uint32_t res, TP2274* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP2204 * TPParent;
TP2274 * controlTPParent;
TP2274* inputsTPParent;double **eta_darts2274 /*OMP_SHARED_PRIVATE - INPUT*/;
int *i_darts2274/*OMP_PRIVATE - INPUT*/;
int **iglob_darts2274 /*OMP_SHARED_PRIVATE - INPUT*/;
int *j_darts2274/*OMP_PRIVATE - INPUT*/;
int **jglob_darts2274 /*OMP_SHARED_PRIVATE - INPUT*/;
int *k_darts2274/*OMP_PRIVATE - INPUT*/;
int *m_darts2274/*OMP_PRIVATE - INPUT*/;
double **xi_darts2274 /*OMP_SHARED_PRIVATE - INPUT*/;
double **zeta_darts2274 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration2274;
int lastIteration2274;
int range2274;
int rangePerCodelet2274;
int minIteration2274;
int remainderRange2274;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets2274* barrierCodelets2274;
_checkInCodelets2275* checkInCodelets2275;
#if USE_SPIN_CODELETS == 0
_checkInCodelets2275* firstCodelet;
#endif
TP2274(int in_numThreads, int in_mainCodeletID, TP2204* in_TPParent, int in_initIteration, int in_lastIteration, TP2274** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP2274();
};
/*TP2409: OMPForDirective*/
class TP2409 : public ompTP
{
 public:
class _barrierCodelets2409 : public darts::Codelet
{
public:
TP2409* inputsTPParent;
_barrierCodelets2409():
darts::Codelet(){ }
_barrierCodelets2409(uint32_t dep, uint32_t res, TP2409* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations2409(int* endRange, uint32_t codeletID);
class _checkInCodelets2410:public darts::Codelet
{
public:
TP2409* myTP;
TP2409* inputsTPParent;
int endRange;
_checkInCodelets2410():
darts::Codelet(){ }
_checkInCodelets2410(uint32_t dep, uint32_t res, TP2409* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP2204 * TPParent;
TP2409 * controlTPParent;
TP2409* inputsTPParent;int **L1_darts2409 /*OMP_SHARED_PRIVATE - INPUT*/;
int **L2_darts2409 /*OMP_SHARED_PRIVATE - INPUT*/;
int *i_darts2409/*OMP_PRIVATE - INPUT*/;
int *j_darts2409/*OMP_PRIVATE - INPUT*/;
int *k_darts2409/*OMP_PRIVATE - INPUT*/;
double **q_darts2409 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u21_darts2409 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration2409;
int lastIteration2409;
int range2409;
int rangePerCodelet2409;
int minIteration2409;
int remainderRange2409;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets2409* barrierCodelets2409;
_checkInCodelets2410* checkInCodelets2410;
#if USE_SPIN_CODELETS == 0
_checkInCodelets2410* firstCodelet;
#endif
TP2409(int in_numThreads, int in_mainCodeletID, TP2204* in_TPParent, int in_initIteration, int in_lastIteration, TP2409** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP2409();
};
/*TP2558: OMPForDirective*/
class TP2558 : public ompTP
{
 public:
class _barrierCodelets2558 : public darts::Codelet
{
public:
TP2558* inputsTPParent;
_barrierCodelets2558():
darts::Codelet(){ }
_barrierCodelets2558(uint32_t dep, uint32_t res, TP2558* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations2558(int* endRange, uint32_t codeletID);
class _checkInCodelets2559:public darts::Codelet
{
public:
TP2558* myTP;
TP2558* inputsTPParent;
int endRange;
_checkInCodelets2559():
darts::Codelet(){ }
_checkInCodelets2559(uint32_t dep, uint32_t res, TP2558* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP2204 * TPParent;
TP2558 * controlTPParent;
TP2558* inputsTPParent;int **L2_darts2558 /*OMP_SHARED_PRIVATE - INPUT*/;
double **dsspm_darts2558 /*OMP_SHARED_PRIVATE - INPUT*/;
int *i_darts2558/*OMP_PRIVATE - INPUT*/;
int **iend1_darts2558 /*OMP_SHARED_PRIVATE - INPUT*/;
int **ist1_darts2558 /*OMP_SHARED_PRIVATE - INPUT*/;
int *j_darts2558/*OMP_PRIVATE - INPUT*/;
int *k_darts2558/*OMP_PRIVATE - INPUT*/;
int *m_darts2558/*OMP_PRIVATE - INPUT*/;
double **tmp_darts2558 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u21i_darts2558 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u21im1_darts2558 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u31i_darts2558 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u31im1_darts2558 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u41i_darts2558 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u41im1_darts2558 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u51i_darts2558 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u51im1_darts2558 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration2558;
int lastIteration2558;
int range2558;
int rangePerCodelet2558;
int minIteration2558;
int remainderRange2558;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets2558* barrierCodelets2558;
_checkInCodelets2559* checkInCodelets2559;
#if USE_SPIN_CODELETS == 0
_checkInCodelets2559* firstCodelet;
#endif
TP2558(int in_numThreads, int in_mainCodeletID, TP2204* in_TPParent, int in_initIteration, int in_lastIteration, TP2558** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP2558();
};
/*TP3196: OMPForDirective*/
class TP3196 : public ompTP
{
 public:
class _barrierCodelets3196 : public darts::Codelet
{
public:
TP3196* inputsTPParent;
_barrierCodelets3196():
darts::Codelet(){ }
_barrierCodelets3196(uint32_t dep, uint32_t res, TP3196* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations3196(int* endRange, uint32_t codeletID);
class _checkInCodelets3197:public darts::Codelet
{
public:
TP3196* myTP;
TP3196* inputsTPParent;
int endRange;
_checkInCodelets3197():
darts::Codelet(){ }
_checkInCodelets3197(uint32_t dep, uint32_t res, TP3196* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP2204 * TPParent;
TP3196 * controlTPParent;
TP3196* inputsTPParent;int **L1_darts3196 /*OMP_SHARED_PRIVATE - INPUT*/;
int **L2_darts3196 /*OMP_SHARED_PRIVATE - INPUT*/;
int *i_darts3196/*OMP_PRIVATE - INPUT*/;
int *j_darts3196/*OMP_PRIVATE - INPUT*/;
int *k_darts3196/*OMP_PRIVATE - INPUT*/;
double **q_darts3196 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u31_darts3196 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration3196;
int lastIteration3196;
int range3196;
int rangePerCodelet3196;
int minIteration3196;
int remainderRange3196;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets3196* barrierCodelets3196;
_checkInCodelets3197* checkInCodelets3197;
#if USE_SPIN_CODELETS == 0
_checkInCodelets3197* firstCodelet;
#endif
TP3196(int in_numThreads, int in_mainCodeletID, TP2204* in_TPParent, int in_initIteration, int in_lastIteration, TP3196** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP3196();
};
/*TP3345: OMPForDirective*/
class TP3345 : public ompTP
{
 public:
class _barrierCodelets3345 : public darts::Codelet
{
public:
TP3345* inputsTPParent;
_barrierCodelets3345():
darts::Codelet(){ }
_barrierCodelets3345(uint32_t dep, uint32_t res, TP3345* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations3345(int* endRange, uint32_t codeletID);
class _checkInCodelets3346:public darts::Codelet
{
public:
TP3345* myTP;
TP3345* inputsTPParent;
int endRange;
_checkInCodelets3346():
darts::Codelet(){ }
_checkInCodelets3346(uint32_t dep, uint32_t res, TP3345* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP2204 * TPParent;
TP3345 * controlTPParent;
TP3345* inputsTPParent;int **L2_darts3345 /*OMP_SHARED_PRIVATE - INPUT*/;
double **dsspm_darts3345 /*OMP_SHARED_PRIVATE - INPUT*/;
int *i_darts3345/*OMP_PRIVATE - INPUT*/;
int *j_darts3345/*OMP_PRIVATE - INPUT*/;
int **jend1_darts3345 /*OMP_SHARED_PRIVATE - INPUT*/;
int **jst1_darts3345 /*OMP_SHARED_PRIVATE - INPUT*/;
int *k_darts3345/*OMP_PRIVATE - INPUT*/;
int *m_darts3345/*OMP_PRIVATE - INPUT*/;
double **tmp_darts3345 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u21j_darts3345 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u21jm1_darts3345 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u31j_darts3345 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u31jm1_darts3345 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u41j_darts3345 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u41jm1_darts3345 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u51j_darts3345 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u51jm1_darts3345 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration3345;
int lastIteration3345;
int range3345;
int rangePerCodelet3345;
int minIteration3345;
int remainderRange3345;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets3345* barrierCodelets3345;
_checkInCodelets3346* checkInCodelets3346;
#if USE_SPIN_CODELETS == 0
_checkInCodelets3346* firstCodelet;
#endif
TP3345(int in_numThreads, int in_mainCodeletID, TP2204* in_TPParent, int in_initIteration, int in_lastIteration, TP3345** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP3345();
};
/*TP3980: OMPForDirective*/
class TP3980 : public ompTP
{
 public:
class _barrierCodelets3980 : public darts::Codelet
{
public:
TP3980* inputsTPParent;
_barrierCodelets3980():
darts::Codelet(){ }
_barrierCodelets3980(uint32_t dep, uint32_t res, TP3980* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations3980(int* endRange, uint32_t codeletID);
class _checkInCodelets3981:public darts::Codelet
{
public:
TP3980* myTP;
TP3980* inputsTPParent;
int endRange;
_checkInCodelets3981():
darts::Codelet(){ }
_checkInCodelets3981(uint32_t dep, uint32_t res, TP3980* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP2204 * TPParent;
TP3980 * controlTPParent;
TP3980* inputsTPParent;double **dsspm_darts3980 /*OMP_SHARED_PRIVATE - INPUT*/;
int *i_darts3980/*OMP_PRIVATE - INPUT*/;
int *j_darts3980/*OMP_PRIVATE - INPUT*/;
int *k_darts3980/*OMP_PRIVATE - INPUT*/;
int *m_darts3980/*OMP_PRIVATE - INPUT*/;
double **q_darts3980 /*OMP_SHARED_PRIVATE - INPUT*/;
double **tmp_darts3980 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u21k_darts3980 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u21km1_darts3980 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u31k_darts3980 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u31km1_darts3980 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u41_darts3980 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u41k_darts3980 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u41km1_darts3980 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u51k_darts3980 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u51km1_darts3980 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration3980;
int lastIteration3980;
int range3980;
int rangePerCodelet3980;
int minIteration3980;
int remainderRange3980;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets3980* barrierCodelets3980;
_checkInCodelets3981* checkInCodelets3981;
#if USE_SPIN_CODELETS == 0
_checkInCodelets3981* firstCodelet;
#endif
TP3980(int in_numThreads, int in_mainCodeletID, TP2204* in_TPParent, int in_initIteration, int in_lastIteration, TP3980** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP3980();
};
/*TP7: jacld*/
class TP7:public ompTP
{
 public:
class _checkInCodelets4885:public darts::Codelet
{
public:
TP7* myTP;
TP7* inputsTPParent;
_checkInCodelets4885():
darts::Codelet(){ }
_checkInCodelets4885(uint32_t dep, uint32_t res, TP7* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets4906:public darts::Codelet
{
public:
TP7* myTP;
TP7* inputsTPParent;
_checkInCodelets4906():
darts::Codelet(){ }
_checkInCodelets4906(uint32_t dep, uint32_t res, TP7* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets4906 : public darts::Codelet
{
public:
TP7* inputsTPParent;
_barrierCodelets4906():
darts::Codelet(){ }
_barrierCodelets4906(uint32_t dep, uint32_t res, TP7* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP7** ptrToThisFunctionTP;
TP7* inputsTPParent;
TP7* controlTPParent;
darts::Codelet** nextCodeletsjacld;
darts::Codelet** nextSyncCodeletsjacld;
int* k_darts7/*VARIABLE*/;
double* c1345_darts7/*VARIABLE*/;
double* c34_darts7/*VARIABLE*/;
int* i_darts7/*VARIABLE*/;
int* j_darts7/*VARIABLE*/;
double* r43_darts7/*VARIABLE*/;
double* tmp1_darts7/*VARIABLE*/;
double* tmp2_darts7/*VARIABLE*/;
double* tmp3_darts7/*VARIABLE*/;
TP4906** TP4906Ptr;
size_t *TP4906_alreadyLaunched;
int numTPsSet4906;
int numTPsReady4906;
size_t TPsToUse4906;
size_t codeletsPerTP4906;
size_t totalCodelets4906;
_checkInCodelets4885* checkInCodelets4885;
#if USE_SPIN_CODELETS == 0
_checkInCodelets4885* firstCodelet;
#endif
_checkInCodelets4906* checkInCodelets4906;
_barrierCodelets4906* barrierCodelets4906;
TP7(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet, darts::Codelet* in_mainSyncCodelet, TP7** in_ptrToThisFunctionTP, int in_k);
~TP7();
void setNewInputs(int in_k, size_t codeletID);};
/*TP4906: OMPForDirective*/
class TP4906 : public ompTP
{
 public:
class _barrierCodelets4906 : public darts::Codelet
{
public:
TP4906* inputsTPParent;
_barrierCodelets4906():
darts::Codelet(){ }
_barrierCodelets4906(uint32_t dep, uint32_t res, TP4906* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations4906(int* endRange, uint32_t codeletID);
class _checkInCodelets4907:public darts::Codelet
{
public:
TP4906* myTP;
TP4906* inputsTPParent;
int endRange;
_checkInCodelets4907():
darts::Codelet(){ }
_checkInCodelets4907(uint32_t dep, uint32_t res, TP4906* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP7 * TPParent;
TP4906 * controlTPParent;
TP4906* inputsTPParent;double **c1345_darts4906 /*OMP_SHARED_PRIVATE - INPUT*/;
double **c34_darts4906 /*OMP_SHARED_PRIVATE - INPUT*/;
int *i_darts4906/*OMP_PRIVATE - INPUT*/;
int *j_darts4906/*OMP_PRIVATE - INPUT*/;
int **k_darts4906 /*OMP_SHARED_PRIVATE - INPUT*/;
double **r43_darts4906 /*OMP_SHARED_PRIVATE - INPUT*/;
double **tmp1_darts4906 /*OMP_SHARED_PRIVATE - INPUT*/;
double **tmp2_darts4906 /*OMP_SHARED_PRIVATE - INPUT*/;
double **tmp3_darts4906 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration4906;
int lastIteration4906;
int range4906;
int rangePerCodelet4906;
int minIteration4906;
int remainderRange4906;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets4906* barrierCodelets4906;
_checkInCodelets4907* checkInCodelets4907;
#if USE_SPIN_CODELETS == 0
_checkInCodelets4907* firstCodelet;
#endif
TP4906(int in_numThreads, int in_mainCodeletID, TP7* in_TPParent, int in_initIteration, int in_lastIteration, TP4906** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP4906();
};
/*TP8: jacu*/
class TP8:public ompTP
{
 public:
class _checkInCodelets7407:public darts::Codelet
{
public:
TP8* myTP;
TP8* inputsTPParent;
_checkInCodelets7407():
darts::Codelet(){ }
_checkInCodelets7407(uint32_t dep, uint32_t res, TP8* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets7428:public darts::Codelet
{
public:
TP8* myTP;
TP8* inputsTPParent;
_checkInCodelets7428():
darts::Codelet(){ }
_checkInCodelets7428(uint32_t dep, uint32_t res, TP8* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets7428 : public darts::Codelet
{
public:
TP8* inputsTPParent;
_barrierCodelets7428():
darts::Codelet(){ }
_barrierCodelets7428(uint32_t dep, uint32_t res, TP8* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP8** ptrToThisFunctionTP;
TP8* inputsTPParent;
TP8* controlTPParent;
darts::Codelet** nextCodeletsjacu;
darts::Codelet** nextSyncCodeletsjacu;
int* k_darts8/*VARIABLE*/;
double* c1345_darts8/*VARIABLE*/;
double* c34_darts8/*VARIABLE*/;
int* i_darts8/*VARIABLE*/;
int* j_darts8/*VARIABLE*/;
double* r43_darts8/*VARIABLE*/;
double* tmp1_darts8/*VARIABLE*/;
double* tmp2_darts8/*VARIABLE*/;
double* tmp3_darts8/*VARIABLE*/;
TP7428** TP7428Ptr;
size_t *TP7428_alreadyLaunched;
int numTPsSet7428;
int numTPsReady7428;
size_t TPsToUse7428;
size_t codeletsPerTP7428;
size_t totalCodelets7428;
_checkInCodelets7407* checkInCodelets7407;
#if USE_SPIN_CODELETS == 0
_checkInCodelets7407* firstCodelet;
#endif
_checkInCodelets7428* checkInCodelets7428;
_barrierCodelets7428* barrierCodelets7428;
TP8(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet, darts::Codelet* in_mainSyncCodelet, TP8** in_ptrToThisFunctionTP, int in_k);
~TP8();
void setNewInputs(int in_k, size_t codeletID);};
/*TP7428: OMPForDirective*/
class TP7428 : public ompTP
{
 public:
class _barrierCodelets7428 : public darts::Codelet
{
public:
TP7428* inputsTPParent;
_barrierCodelets7428():
darts::Codelet(){ }
_barrierCodelets7428(uint32_t dep, uint32_t res, TP7428* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations7428(int* endRange, uint32_t codeletID);
class _checkInCodelets7429:public darts::Codelet
{
public:
TP7428* myTP;
TP7428* inputsTPParent;
int endRange;
_checkInCodelets7429():
darts::Codelet(){ }
_checkInCodelets7429(uint32_t dep, uint32_t res, TP7428* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP8 * TPParent;
TP7428 * controlTPParent;
TP7428* inputsTPParent;double **c1345_darts7428 /*OMP_SHARED_PRIVATE - INPUT*/;
double **c34_darts7428 /*OMP_SHARED_PRIVATE - INPUT*/;
int *i_darts7428/*OMP_PRIVATE - INPUT*/;
int *j_darts7428/*OMP_PRIVATE - INPUT*/;
int **k_darts7428 /*OMP_SHARED_PRIVATE - INPUT*/;
double **r43_darts7428 /*OMP_SHARED_PRIVATE - INPUT*/;
double **tmp1_darts7428 /*OMP_SHARED_PRIVATE - INPUT*/;
double **tmp2_darts7428 /*OMP_SHARED_PRIVATE - INPUT*/;
double **tmp3_darts7428 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration7428;
int lastIteration7428;
int range7428;
int rangePerCodelet7428;
int minIteration7428;
int remainderRange7428;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets7428* barrierCodelets7428;
_checkInCodelets7429* checkInCodelets7429;
#if USE_SPIN_CODELETS == 0
_checkInCodelets7429* firstCodelet;
#endif
TP7428(int in_numThreads, int in_mainCodeletID, TP8* in_TPParent, int in_initIteration, int in_lastIteration, TP7428** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP7428();
};
/*TP9879: OMPParallelDirective*/
class TP9879 : public darts::ThreadedProcedure
{
 public:
class _barrierCodelets9879 : public darts::Codelet
{
public:
TP9879* inputsTPParent;
_barrierCodelets9879():
darts::Codelet(){ }
_barrierCodelets9879(uint32_t dep, uint32_t res, TP9879* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets9881:public darts::Codelet
{
public:
TP9879* myTP;
TP9879* inputsTPParent;
_checkInCodelets9881():
darts::Codelet(){ }
_checkInCodelets9881(uint32_t dep, uint32_t res, TP9879* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets9888:public darts::Codelet
{
public:
TP9879* myTP;
TP9879* inputsTPParent;
_checkInCodelets9888():
darts::Codelet(){ }
_checkInCodelets9888(uint32_t dep, uint32_t res, TP9879* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets9888 : public darts::Codelet
{
public:
TP9879* inputsTPParent;
_barrierCodelets9888():
darts::Codelet(){ }
_barrierCodelets9888(uint32_t dep, uint32_t res, TP9879* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets9897:public darts::Codelet
{
public:
TP9879* myTP;
TP9879* inputsTPParent;
_checkInCodelets9897():
darts::Codelet(){ }
_checkInCodelets9897(uint32_t dep, uint32_t res, TP9879* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets9994:public darts::Codelet
{
public:
TP9879* myTP;
TP9879* inputsTPParent;
_checkInCodelets9994():
darts::Codelet(){ }
_checkInCodelets9994(uint32_t dep, uint32_t res, TP9879* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets10004 : public darts::Codelet
{
public:
TP9879* inputsTPParent;
_barrierCodelets10004():
darts::Codelet(){ }
_barrierCodelets10004(uint32_t dep, uint32_t res, TP9879* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets10005:public darts::Codelet
{
public:
TP9879* myTP;
TP9879* inputsTPParent;
_checkInCodelets10005():
darts::Codelet(){ }
_checkInCodelets10005(uint32_t dep, uint32_t res, TP9879* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets10005 : public darts::Codelet
{
public:
TP9879* inputsTPParent;
_barrierCodelets10005():
darts::Codelet(){ }
_barrierCodelets10005(uint32_t dep, uint32_t res, TP9879* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
darts::Codelet* nextCodelet;
TP9879 * TPParent;
TP9879 * controlTPParent;
TP9879* inputsTPParent;int *iend_darts9879;/*OMP_SHARED - INPUT*/
int *ist_darts9879;/*OMP_SHARED - INPUT*/
int *jend_darts9879;/*OMP_SHARED - INPUT*/
int *jst_darts9879;/*OMP_SHARED - INPUT*/
int *nx0_darts9879;/*OMP_SHARED - INPUT*/
int *ny0_darts9879;/*OMP_SHARED - INPUT*/
int *nz0_darts9879;/*OMP_SHARED - INPUT*/
double * *sum_darts9879;/*OMP_SHARED - INPUT*/
int *i_darts9879/*VARIABLE*/;
int *j_darts9879/*VARIABLE*/;
int *k_darts9879/*VARIABLE*/;
int *m_darts9879/*VARIABLE*/;
double *sum0_darts9879/*VARIABLE*/;
double *sum1_darts9879/*VARIABLE*/;
double *sum2_darts9879/*VARIABLE*/;
double *sum3_darts9879/*VARIABLE*/;
double *sum4_darts9879/*VARIABLE*/;
int m_darts10005 ;
int *nx0_darts10005;/*OMP_SHARED - INPUT*/
int *ny0_darts10005;/*OMP_SHARED - INPUT*/
int *nz0_darts10005;/*OMP_SHARED - INPUT*/
double * *sum_darts10005;/*OMP_SHARED - INPUT*/
int m_darts9888 ;
double * *sum_darts9888;/*OMP_SHARED - INPUT*/
size_t TP9888_alreadyLaunched;
TP9897** TP9897Ptr;
size_t *TP9897_alreadyLaunched;
int numTPsSet9897;
int numTPsReady9897;
size_t TPsToUse9897;
size_t codeletsPerTP9897;
size_t totalCodelets9897;
size_t TP10005_alreadyLaunched;
_barrierCodelets9879* barrierCodelets9879;
_checkInCodelets9881* checkInCodelets9881;
_checkInCodelets9888* checkInCodelets9888;
_barrierCodelets9888* barrierCodelets9888;
_checkInCodelets9897* checkInCodelets9897;
_checkInCodelets9994* checkInCodelets9994;
_barrierCodelets10004* barrierCodelets10004;
_checkInCodelets10005* checkInCodelets10005;
_barrierCodelets10005* barrierCodelets10005;
TP9879(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet, int *in_iend, int *in_ist, int *in_jend, int *in_jst, int *in_nx0, int *in_ny0, int *in_nz0, double * *in_sum);
~TP9879();
};
/*TP9897: OMPForDirective*/
class TP9897 : public ompTP
{
 public:
bool requestNewRangeIterations9897(int* endRange, uint32_t codeletID);
class _checkInCodelets9898:public darts::Codelet
{
public:
TP9897* myTP;
TP9897* inputsTPParent;
int endRange;
_checkInCodelets9898():
darts::Codelet(){ }
_checkInCodelets9898(uint32_t dep, uint32_t res, TP9897* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP9879 * TPParent;
TP9897 * controlTPParent;
TP9897* inputsTPParent;int *i_darts9897/*OMP_PRIVATE - INPUT*/;
int *iend_darts9897;/*OMP_SHARED - INPUT*/
int *ist_darts9897;/*OMP_SHARED - INPUT*/
int *j_darts9897/*OMP_PRIVATE - INPUT*/;
int *jend_darts9897;/*OMP_SHARED - INPUT*/
int *jst_darts9897;/*OMP_SHARED - INPUT*/
int *k_darts9897/*OMP_PRIVATE - INPUT*/;
int *nz0_darts9897;/*OMP_SHARED - INPUT*/
double **sum0_darts9897 /*OMP_SHARED_PRIVATE - INPUT*/;
double **sum1_darts9897 /*OMP_SHARED_PRIVATE - INPUT*/;
double **sum2_darts9897 /*OMP_SHARED_PRIVATE - INPUT*/;
double **sum3_darts9897 /*OMP_SHARED_PRIVATE - INPUT*/;
double **sum4_darts9897 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration9897;
int lastIteration9897;
int range9897;
int rangePerCodelet9897;
int minIteration9897;
int remainderRange9897;
size_t readyCodelets;
int baseNumThreads;
int *signalNextReady;
_checkInCodelets9898* checkInCodelets9898;
#if USE_SPIN_CODELETS == 0
_checkInCodelets9898* firstCodelet;
#endif
TP9897(int in_numThreads, int in_mainCodeletID, TP9879* in_TPParent, int in_initIteration, int in_lastIteration, int *in_iend, int *in_ist, int *in_jend, int *in_jst, int *in_nz0, TP9897** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP9897();
};
/*TP10796: OMPParallelDirective*/
class TP10796 : public darts::ThreadedProcedure
{
 public:
class _barrierCodelets10796 : public darts::Codelet
{
public:
TP10796* inputsTPParent;
_barrierCodelets10796():
darts::Codelet(){ }
_barrierCodelets10796(uint32_t dep, uint32_t res, TP10796* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets10811:public darts::Codelet
{
public:
TP10796* myTP;
TP10796* inputsTPParent;
_checkInCodelets10811():
darts::Codelet(){ }
_checkInCodelets10811(uint32_t dep, uint32_t res, TP10796* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets10811 : public darts::Codelet
{
public:
TP10796* inputsTPParent;
_barrierCodelets10811():
darts::Codelet(){ }
_barrierCodelets10811(uint32_t dep, uint32_t res, TP10796* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets10868:public darts::Codelet
{
public:
TP10796* myTP;
TP10796* inputsTPParent;
_checkInCodelets10868():
darts::Codelet(){ }
_checkInCodelets10868(uint32_t dep, uint32_t res, TP10796* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets10871:public darts::Codelet
{
public:
TP10796* myTP;
TP10796* inputsTPParent;
_checkInCodelets10871():
darts::Codelet(){ }
_checkInCodelets10871(uint32_t dep, uint32_t res, TP10796* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets10871 : public darts::Codelet
{
public:
TP10796* inputsTPParent;
_barrierCodelets10871():
darts::Codelet(){ }
_barrierCodelets10871(uint32_t dep, uint32_t res, TP10796* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets11020:public darts::Codelet
{
public:
TP10796* myTP;
TP10796* inputsTPParent;
_checkInCodelets11020():
darts::Codelet(){ }
_checkInCodelets11020(uint32_t dep, uint32_t res, TP10796* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets11020 : public darts::Codelet
{
public:
TP10796* inputsTPParent;
_barrierCodelets11020():
darts::Codelet(){ }
_barrierCodelets11020(uint32_t dep, uint32_t res, TP10796* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets11657:public darts::Codelet
{
public:
TP10796* myTP;
TP10796* inputsTPParent;
_checkInCodelets11657():
darts::Codelet(){ }
_checkInCodelets11657(uint32_t dep, uint32_t res, TP10796* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets11660:public darts::Codelet
{
public:
TP10796* myTP;
TP10796* inputsTPParent;
_checkInCodelets11660():
darts::Codelet(){ }
_checkInCodelets11660(uint32_t dep, uint32_t res, TP10796* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets11660 : public darts::Codelet
{
public:
TP10796* inputsTPParent;
_barrierCodelets11660():
darts::Codelet(){ }
_barrierCodelets11660(uint32_t dep, uint32_t res, TP10796* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets11809:public darts::Codelet
{
public:
TP10796* myTP;
TP10796* inputsTPParent;
_checkInCodelets11809():
darts::Codelet(){ }
_checkInCodelets11809(uint32_t dep, uint32_t res, TP10796* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets11809 : public darts::Codelet
{
public:
TP10796* inputsTPParent;
_barrierCodelets11809():
darts::Codelet(){ }
_barrierCodelets11809(uint32_t dep, uint32_t res, TP10796* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets12446:public darts::Codelet
{
public:
TP10796* myTP;
TP10796* inputsTPParent;
_checkInCodelets12446():
darts::Codelet(){ }
_checkInCodelets12446(uint32_t dep, uint32_t res, TP10796* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets12446 : public darts::Codelet
{
public:
TP10796* inputsTPParent;
_barrierCodelets12446():
darts::Codelet(){ }
_barrierCodelets12446(uint32_t dep, uint32_t res, TP10796* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
darts::Codelet* nextCodelet;
TP10796 * TPParent;
TP10796 * controlTPParent;
TP10796* inputsTPParent;int *L1_darts10796/*VARIABLE*/;
int *L2_darts10796/*VARIABLE*/;
int *i_darts10796/*VARIABLE*/;
int *iend1_darts10796/*VARIABLE*/;
int *ist1_darts10796/*VARIABLE*/;
int *j_darts10796/*VARIABLE*/;
int *jend1_darts10796/*VARIABLE*/;
int *jst1_darts10796/*VARIABLE*/;
int *k_darts10796/*VARIABLE*/;
int *m_darts10796/*VARIABLE*/;
double *q_darts10796/*VARIABLE*/;
double *tmp_darts10796/*VARIABLE*/;
double *u21_darts10796/*VARIABLE*/;
double *u21i_darts10796/*VARIABLE*/;
double *u21im1_darts10796/*VARIABLE*/;
double *u21j_darts10796/*VARIABLE*/;
double *u21jm1_darts10796/*VARIABLE*/;
double *u21k_darts10796/*VARIABLE*/;
double *u21km1_darts10796/*VARIABLE*/;
double *u31_darts10796/*VARIABLE*/;
double *u31i_darts10796/*VARIABLE*/;
double *u31im1_darts10796/*VARIABLE*/;
double *u31j_darts10796/*VARIABLE*/;
double *u31jm1_darts10796/*VARIABLE*/;
double *u31k_darts10796/*VARIABLE*/;
double *u31km1_darts10796/*VARIABLE*/;
double *u41_darts10796/*VARIABLE*/;
double *u41i_darts10796/*VARIABLE*/;
double *u41im1_darts10796/*VARIABLE*/;
double *u41j_darts10796/*VARIABLE*/;
double *u41jm1_darts10796/*VARIABLE*/;
double *u41k_darts10796/*VARIABLE*/;
double *u41km1_darts10796/*VARIABLE*/;
double *u51i_darts10796/*VARIABLE*/;
double *u51im1_darts10796/*VARIABLE*/;
double *u51j_darts10796/*VARIABLE*/;
double *u51jm1_darts10796/*VARIABLE*/;
double *u51k_darts10796/*VARIABLE*/;
double *u51km1_darts10796/*VARIABLE*/;
TP10811** TP10811Ptr;
size_t *TP10811_alreadyLaunched;
int numTPsSet10811;
int numTPsReady10811;
size_t TPsToUse10811;
size_t codeletsPerTP10811;
size_t totalCodelets10811;
TP10871** TP10871Ptr;
size_t *TP10871_alreadyLaunched;
int numTPsSet10871;
int numTPsReady10871;
size_t TPsToUse10871;
size_t codeletsPerTP10871;
size_t totalCodelets10871;
TP11020** TP11020Ptr;
size_t *TP11020_alreadyLaunched;
int numTPsSet11020;
int numTPsReady11020;
size_t TPsToUse11020;
size_t codeletsPerTP11020;
size_t totalCodelets11020;
TP11660** TP11660Ptr;
size_t *TP11660_alreadyLaunched;
int numTPsSet11660;
int numTPsReady11660;
size_t TPsToUse11660;
size_t codeletsPerTP11660;
size_t totalCodelets11660;
TP11809** TP11809Ptr;
size_t *TP11809_alreadyLaunched;
int numTPsSet11809;
int numTPsReady11809;
size_t TPsToUse11809;
size_t codeletsPerTP11809;
size_t totalCodelets11809;
TP12446** TP12446Ptr;
size_t *TP12446_alreadyLaunched;
int numTPsSet12446;
int numTPsReady12446;
size_t TPsToUse12446;
size_t codeletsPerTP12446;
size_t totalCodelets12446;
_barrierCodelets10796* barrierCodelets10796;
_checkInCodelets10811* checkInCodelets10811;
_barrierCodelets10811* barrierCodelets10811;
_checkInCodelets10868* checkInCodelets10868;
_checkInCodelets10871* checkInCodelets10871;
_barrierCodelets10871* barrierCodelets10871;
_checkInCodelets11020* checkInCodelets11020;
_barrierCodelets11020* barrierCodelets11020;
_checkInCodelets11657* checkInCodelets11657;
_checkInCodelets11660* checkInCodelets11660;
_barrierCodelets11660* barrierCodelets11660;
_checkInCodelets11809* checkInCodelets11809;
_barrierCodelets11809* barrierCodelets11809;
_checkInCodelets12446* checkInCodelets12446;
_barrierCodelets12446* barrierCodelets12446;
TP10796(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
~TP10796();
};
/*TP10811: OMPForDirective*/
class TP10811 : public ompTP
{
 public:
class _barrierCodelets10811 : public darts::Codelet
{
public:
TP10811* inputsTPParent;
_barrierCodelets10811():
darts::Codelet(){ }
_barrierCodelets10811(uint32_t dep, uint32_t res, TP10811* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations10811(int* endRange, uint32_t codeletID);
class _checkInCodelets10812:public darts::Codelet
{
public:
TP10811* myTP;
TP10811* inputsTPParent;
int endRange;
_checkInCodelets10812():
darts::Codelet(){ }
_checkInCodelets10812(uint32_t dep, uint32_t res, TP10811* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP10796 * TPParent;
TP10811 * controlTPParent;
TP10811* inputsTPParent;int *i_darts10811/*OMP_PRIVATE - INPUT*/;
int *j_darts10811/*OMP_PRIVATE - INPUT*/;
int *k_darts10811/*OMP_PRIVATE - INPUT*/;
int *m_darts10811/*OMP_PRIVATE - INPUT*/;
int initIteration10811;
int lastIteration10811;
int range10811;
int rangePerCodelet10811;
int minIteration10811;
int remainderRange10811;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets10811* barrierCodelets10811;
_checkInCodelets10812* checkInCodelets10812;
#if USE_SPIN_CODELETS == 0
_checkInCodelets10812* firstCodelet;
#endif
TP10811(int in_numThreads, int in_mainCodeletID, TP10796* in_TPParent, int in_initIteration, int in_lastIteration, TP10811** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP10811();
};
/*TP10871: OMPForDirective*/
class TP10871 : public ompTP
{
 public:
class _barrierCodelets10871 : public darts::Codelet
{
public:
TP10871* inputsTPParent;
_barrierCodelets10871():
darts::Codelet(){ }
_barrierCodelets10871(uint32_t dep, uint32_t res, TP10871* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations10871(int* endRange, uint32_t codeletID);
class _checkInCodelets10872:public darts::Codelet
{
public:
TP10871* myTP;
TP10871* inputsTPParent;
int endRange;
_checkInCodelets10872():
darts::Codelet(){ }
_checkInCodelets10872(uint32_t dep, uint32_t res, TP10871* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP10796 * TPParent;
TP10871 * controlTPParent;
TP10871* inputsTPParent;int **L1_darts10871 /*OMP_SHARED_PRIVATE - INPUT*/;
int **L2_darts10871 /*OMP_SHARED_PRIVATE - INPUT*/;
int *i_darts10871/*OMP_PRIVATE - INPUT*/;
int *j_darts10871/*OMP_PRIVATE - INPUT*/;
int *k_darts10871/*OMP_PRIVATE - INPUT*/;
double **q_darts10871 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u21_darts10871 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration10871;
int lastIteration10871;
int range10871;
int rangePerCodelet10871;
int minIteration10871;
int remainderRange10871;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets10871* barrierCodelets10871;
_checkInCodelets10872* checkInCodelets10872;
#if USE_SPIN_CODELETS == 0
_checkInCodelets10872* firstCodelet;
#endif
TP10871(int in_numThreads, int in_mainCodeletID, TP10796* in_TPParent, int in_initIteration, int in_lastIteration, TP10871** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP10871();
};
/*TP11020: OMPForDirective*/
class TP11020 : public ompTP
{
 public:
class _barrierCodelets11020 : public darts::Codelet
{
public:
TP11020* inputsTPParent;
_barrierCodelets11020():
darts::Codelet(){ }
_barrierCodelets11020(uint32_t dep, uint32_t res, TP11020* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations11020(int* endRange, uint32_t codeletID);
class _checkInCodelets11021:public darts::Codelet
{
public:
TP11020* myTP;
TP11020* inputsTPParent;
int endRange;
_checkInCodelets11021():
darts::Codelet(){ }
_checkInCodelets11021(uint32_t dep, uint32_t res, TP11020* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP10796 * TPParent;
TP11020 * controlTPParent;
TP11020* inputsTPParent;int **L2_darts11020 /*OMP_SHARED_PRIVATE - INPUT*/;
int *i_darts11020/*OMP_PRIVATE - INPUT*/;
int **iend1_darts11020 /*OMP_SHARED_PRIVATE - INPUT*/;
int **ist1_darts11020 /*OMP_SHARED_PRIVATE - INPUT*/;
int *j_darts11020/*OMP_PRIVATE - INPUT*/;
int *k_darts11020/*OMP_PRIVATE - INPUT*/;
int *m_darts11020/*OMP_PRIVATE - INPUT*/;
double **tmp_darts11020 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u21i_darts11020 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u21im1_darts11020 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u31i_darts11020 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u31im1_darts11020 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u41i_darts11020 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u41im1_darts11020 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u51i_darts11020 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u51im1_darts11020 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration11020;
int lastIteration11020;
int range11020;
int rangePerCodelet11020;
int minIteration11020;
int remainderRange11020;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets11020* barrierCodelets11020;
_checkInCodelets11021* checkInCodelets11021;
#if USE_SPIN_CODELETS == 0
_checkInCodelets11021* firstCodelet;
#endif
TP11020(int in_numThreads, int in_mainCodeletID, TP10796* in_TPParent, int in_initIteration, int in_lastIteration, TP11020** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP11020();
};
/*TP11660: OMPForDirective*/
class TP11660 : public ompTP
{
 public:
class _barrierCodelets11660 : public darts::Codelet
{
public:
TP11660* inputsTPParent;
_barrierCodelets11660():
darts::Codelet(){ }
_barrierCodelets11660(uint32_t dep, uint32_t res, TP11660* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations11660(int* endRange, uint32_t codeletID);
class _checkInCodelets11661:public darts::Codelet
{
public:
TP11660* myTP;
TP11660* inputsTPParent;
int endRange;
_checkInCodelets11661():
darts::Codelet(){ }
_checkInCodelets11661(uint32_t dep, uint32_t res, TP11660* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP10796 * TPParent;
TP11660 * controlTPParent;
TP11660* inputsTPParent;int **L1_darts11660 /*OMP_SHARED_PRIVATE - INPUT*/;
int **L2_darts11660 /*OMP_SHARED_PRIVATE - INPUT*/;
int *i_darts11660/*OMP_PRIVATE - INPUT*/;
int *j_darts11660/*OMP_PRIVATE - INPUT*/;
int *k_darts11660/*OMP_PRIVATE - INPUT*/;
double **q_darts11660 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u31_darts11660 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration11660;
int lastIteration11660;
int range11660;
int rangePerCodelet11660;
int minIteration11660;
int remainderRange11660;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets11660* barrierCodelets11660;
_checkInCodelets11661* checkInCodelets11661;
#if USE_SPIN_CODELETS == 0
_checkInCodelets11661* firstCodelet;
#endif
TP11660(int in_numThreads, int in_mainCodeletID, TP10796* in_TPParent, int in_initIteration, int in_lastIteration, TP11660** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP11660();
};
/*TP11809: OMPForDirective*/
class TP11809 : public ompTP
{
 public:
class _barrierCodelets11809 : public darts::Codelet
{
public:
TP11809* inputsTPParent;
_barrierCodelets11809():
darts::Codelet(){ }
_barrierCodelets11809(uint32_t dep, uint32_t res, TP11809* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations11809(int* endRange, uint32_t codeletID);
class _checkInCodelets11810:public darts::Codelet
{
public:
TP11809* myTP;
TP11809* inputsTPParent;
int endRange;
_checkInCodelets11810():
darts::Codelet(){ }
_checkInCodelets11810(uint32_t dep, uint32_t res, TP11809* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP10796 * TPParent;
TP11809 * controlTPParent;
TP11809* inputsTPParent;int **L2_darts11809 /*OMP_SHARED_PRIVATE - INPUT*/;
int *i_darts11809/*OMP_PRIVATE - INPUT*/;
int *j_darts11809/*OMP_PRIVATE - INPUT*/;
int **jend1_darts11809 /*OMP_SHARED_PRIVATE - INPUT*/;
int **jst1_darts11809 /*OMP_SHARED_PRIVATE - INPUT*/;
int *k_darts11809/*OMP_PRIVATE - INPUT*/;
int *m_darts11809/*OMP_PRIVATE - INPUT*/;
double **tmp_darts11809 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u21j_darts11809 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u21jm1_darts11809 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u31j_darts11809 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u31jm1_darts11809 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u41j_darts11809 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u41jm1_darts11809 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u51j_darts11809 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u51jm1_darts11809 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration11809;
int lastIteration11809;
int range11809;
int rangePerCodelet11809;
int minIteration11809;
int remainderRange11809;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets11809* barrierCodelets11809;
_checkInCodelets11810* checkInCodelets11810;
#if USE_SPIN_CODELETS == 0
_checkInCodelets11810* firstCodelet;
#endif
TP11809(int in_numThreads, int in_mainCodeletID, TP10796* in_TPParent, int in_initIteration, int in_lastIteration, TP11809** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP11809();
};
/*TP12446: OMPForDirective*/
class TP12446 : public ompTP
{
 public:
class _barrierCodelets12446 : public darts::Codelet
{
public:
TP12446* inputsTPParent;
_barrierCodelets12446():
darts::Codelet(){ }
_barrierCodelets12446(uint32_t dep, uint32_t res, TP12446* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations12446(int* endRange, uint32_t codeletID);
class _checkInCodelets12447:public darts::Codelet
{
public:
TP12446* myTP;
TP12446* inputsTPParent;
int endRange;
_checkInCodelets12447():
darts::Codelet(){ }
_checkInCodelets12447(uint32_t dep, uint32_t res, TP12446* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP10796 * TPParent;
TP12446 * controlTPParent;
TP12446* inputsTPParent;int *i_darts12446/*OMP_PRIVATE - INPUT*/;
int *j_darts12446/*OMP_PRIVATE - INPUT*/;
int *k_darts12446/*OMP_PRIVATE - INPUT*/;
int *m_darts12446/*OMP_PRIVATE - INPUT*/;
double **q_darts12446 /*OMP_SHARED_PRIVATE - INPUT*/;
double **tmp_darts12446 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u21k_darts12446 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u21km1_darts12446 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u31k_darts12446 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u31km1_darts12446 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u41_darts12446 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u41k_darts12446 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u41km1_darts12446 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u51k_darts12446 /*OMP_SHARED_PRIVATE - INPUT*/;
double **u51km1_darts12446 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration12446;
int lastIteration12446;
int range12446;
int rangePerCodelet12446;
int minIteration12446;
int remainderRange12446;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets12446* barrierCodelets12446;
_checkInCodelets12447* checkInCodelets12447;
#if USE_SPIN_CODELETS == 0
_checkInCodelets12447* firstCodelet;
#endif
TP12446(int in_numThreads, int in_mainCodeletID, TP10796* in_TPParent, int in_initIteration, int in_lastIteration, TP12446** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP12446();
};
/*TP13197: OMPParallelDirective*/
class TP13197 : public darts::ThreadedProcedure
{
 public:
class _barrierCodelets13197 : public darts::Codelet
{
public:
TP13197* inputsTPParent;
_barrierCodelets13197():
darts::Codelet(){ }
_barrierCodelets13197(uint32_t dep, uint32_t res, TP13197* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets13201:public darts::Codelet
{
public:
TP13197* myTP;
TP13197* inputsTPParent;
_checkInCodelets13201():
darts::Codelet(){ }
_checkInCodelets13201(uint32_t dep, uint32_t res, TP13197* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets13201 : public darts::Codelet
{
public:
TP13197* inputsTPParent;
_barrierCodelets13201():
darts::Codelet(){ }
_barrierCodelets13201(uint32_t dep, uint32_t res, TP13197* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets13252:public darts::Codelet
{
public:
TP13197* myTP;
TP13197* inputsTPParent;
_checkInCodelets13252():
darts::Codelet(){ }
_checkInCodelets13252(uint32_t dep, uint32_t res, TP13197* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets13252 : public darts::Codelet
{
public:
TP13197* inputsTPParent;
_barrierCodelets13252():
darts::Codelet(){ }
_barrierCodelets13252(uint32_t dep, uint32_t res, TP13197* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets13294:public darts::Codelet
{
public:
TP13197* myTP;
TP13197* inputsTPParent;
_checkInCodelets13294():
darts::Codelet(){ }
_checkInCodelets13294(uint32_t dep, uint32_t res, TP13197* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets13294 : public darts::Codelet
{
public:
TP13197* inputsTPParent;
_barrierCodelets13294():
darts::Codelet(){ }
_barrierCodelets13294(uint32_t dep, uint32_t res, TP13197* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets13338:public darts::Codelet
{
public:
TP13197* myTP;
TP13197* inputsTPParent;
_checkInCodelets13338():
darts::Codelet(){ }
_checkInCodelets13338(uint32_t dep, uint32_t res, TP13197* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets13338 : public darts::Codelet
{
public:
TP13197* inputsTPParent;
_barrierCodelets13338():
darts::Codelet(){ }
_barrierCodelets13338(uint32_t dep, uint32_t res, TP13197* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets13380:public darts::Codelet
{
public:
TP13197* myTP;
TP13197* inputsTPParent;
_checkInCodelets13380():
darts::Codelet(){ }
_checkInCodelets13380(uint32_t dep, uint32_t res, TP13197* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets13380 : public darts::Codelet
{
public:
TP13197* inputsTPParent;
_barrierCodelets13380():
darts::Codelet(){ }
_barrierCodelets13380(uint32_t dep, uint32_t res, TP13197* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
darts::Codelet* nextCodelet;
TP13197 * TPParent;
TP13197 * controlTPParent;
TP13197* inputsTPParent;int *i_darts13197/*VARIABLE*/;
int *iglob_darts13197/*VARIABLE*/;
int *j_darts13197/*VARIABLE*/;
int *jglob_darts13197/*VARIABLE*/;
int *k_darts13197/*VARIABLE*/;
TP13201** TP13201Ptr;
size_t *TP13201_alreadyLaunched;
int numTPsSet13201;
int numTPsReady13201;
size_t TPsToUse13201;
size_t codeletsPerTP13201;
size_t totalCodelets13201;
TP13252** TP13252Ptr;
size_t *TP13252_alreadyLaunched;
int numTPsSet13252;
int numTPsReady13252;
size_t TPsToUse13252;
size_t codeletsPerTP13252;
size_t totalCodelets13252;
TP13294** TP13294Ptr;
size_t *TP13294_alreadyLaunched;
int numTPsSet13294;
int numTPsReady13294;
size_t TPsToUse13294;
size_t codeletsPerTP13294;
size_t totalCodelets13294;
TP13338** TP13338Ptr;
size_t *TP13338_alreadyLaunched;
int numTPsSet13338;
int numTPsReady13338;
size_t TPsToUse13338;
size_t codeletsPerTP13338;
size_t totalCodelets13338;
TP13380** TP13380Ptr;
size_t *TP13380_alreadyLaunched;
int numTPsSet13380;
int numTPsReady13380;
size_t TPsToUse13380;
size_t codeletsPerTP13380;
size_t totalCodelets13380;
_barrierCodelets13197* barrierCodelets13197;
_checkInCodelets13201* checkInCodelets13201;
_barrierCodelets13201* barrierCodelets13201;
_checkInCodelets13252* checkInCodelets13252;
_barrierCodelets13252* barrierCodelets13252;
_checkInCodelets13294* checkInCodelets13294;
_barrierCodelets13294* barrierCodelets13294;
_checkInCodelets13338* checkInCodelets13338;
_barrierCodelets13338* barrierCodelets13338;
_checkInCodelets13380* checkInCodelets13380;
_barrierCodelets13380* barrierCodelets13380;
TP13197(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
~TP13197();
};
/*TP13201: OMPForDirective*/
class TP13201 : public ompTP
{
 public:
class _barrierCodelets13201 : public darts::Codelet
{
public:
TP13201* inputsTPParent;
_barrierCodelets13201():
darts::Codelet(){ }
_barrierCodelets13201(uint32_t dep, uint32_t res, TP13201* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations13201(int* endRange, uint32_t codeletID);
class _checkInCodelets13202:public darts::Codelet
{
public:
TP13201* myTP;
TP13201* inputsTPParent;
int endRange;
_checkInCodelets13202():
darts::Codelet(){ }
_checkInCodelets13202(uint32_t dep, uint32_t res, TP13201* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP13197 * TPParent;
TP13201 * controlTPParent;
TP13201* inputsTPParent;int *i_darts13201/*OMP_PRIVATE - INPUT*/;
int **iglob_darts13201 /*OMP_SHARED_PRIVATE - INPUT*/;
int *j_darts13201/*OMP_PRIVATE - INPUT*/;
int **jglob_darts13201 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration13201;
int lastIteration13201;
int range13201;
int rangePerCodelet13201;
int minIteration13201;
int remainderRange13201;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets13201* barrierCodelets13201;
_checkInCodelets13202* checkInCodelets13202;
#if USE_SPIN_CODELETS == 0
_checkInCodelets13202* firstCodelet;
#endif
TP13201(int in_numThreads, int in_mainCodeletID, TP13197* in_TPParent, int in_initIteration, int in_lastIteration, TP13201** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP13201();
};
/*TP13252: OMPForDirective*/
class TP13252 : public ompTP
{
 public:
class _barrierCodelets13252 : public darts::Codelet
{
public:
TP13252* inputsTPParent;
_barrierCodelets13252():
darts::Codelet(){ }
_barrierCodelets13252(uint32_t dep, uint32_t res, TP13252* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations13252(int* endRange, uint32_t codeletID);
class _checkInCodelets13253:public darts::Codelet
{
public:
TP13252* myTP;
TP13252* inputsTPParent;
int endRange;
_checkInCodelets13253():
darts::Codelet(){ }
_checkInCodelets13253(uint32_t dep, uint32_t res, TP13252* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP13197 * TPParent;
TP13252 * controlTPParent;
TP13252* inputsTPParent;int *i_darts13252/*OMP_PRIVATE - INPUT*/;
int **iglob_darts13252 /*OMP_SHARED_PRIVATE - INPUT*/;
int *k_darts13252/*OMP_PRIVATE - INPUT*/;
int initIteration13252;
int lastIteration13252;
int range13252;
int rangePerCodelet13252;
int minIteration13252;
int remainderRange13252;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets13252* barrierCodelets13252;
_checkInCodelets13253* checkInCodelets13253;
#if USE_SPIN_CODELETS == 0
_checkInCodelets13253* firstCodelet;
#endif
TP13252(int in_numThreads, int in_mainCodeletID, TP13197* in_TPParent, int in_initIteration, int in_lastIteration, TP13252** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP13252();
};
/*TP13294: OMPForDirective*/
class TP13294 : public ompTP
{
 public:
class _barrierCodelets13294 : public darts::Codelet
{
public:
TP13294* inputsTPParent;
_barrierCodelets13294():
darts::Codelet(){ }
_barrierCodelets13294(uint32_t dep, uint32_t res, TP13294* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations13294(int* endRange, uint32_t codeletID);
class _checkInCodelets13295:public darts::Codelet
{
public:
TP13294* myTP;
TP13294* inputsTPParent;
int endRange;
_checkInCodelets13295():
darts::Codelet(){ }
_checkInCodelets13295(uint32_t dep, uint32_t res, TP13294* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP13197 * TPParent;
TP13294 * controlTPParent;
TP13294* inputsTPParent;int *i_darts13294/*OMP_PRIVATE - INPUT*/;
int **iglob_darts13294 /*OMP_SHARED_PRIVATE - INPUT*/;
int *k_darts13294/*OMP_PRIVATE - INPUT*/;
int initIteration13294;
int lastIteration13294;
int range13294;
int rangePerCodelet13294;
int minIteration13294;
int remainderRange13294;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets13294* barrierCodelets13294;
_checkInCodelets13295* checkInCodelets13295;
#if USE_SPIN_CODELETS == 0
_checkInCodelets13295* firstCodelet;
#endif
TP13294(int in_numThreads, int in_mainCodeletID, TP13197* in_TPParent, int in_initIteration, int in_lastIteration, TP13294** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP13294();
};
/*TP13338: OMPForDirective*/
class TP13338 : public ompTP
{
 public:
class _barrierCodelets13338 : public darts::Codelet
{
public:
TP13338* inputsTPParent;
_barrierCodelets13338():
darts::Codelet(){ }
_barrierCodelets13338(uint32_t dep, uint32_t res, TP13338* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations13338(int* endRange, uint32_t codeletID);
class _checkInCodelets13339:public darts::Codelet
{
public:
TP13338* myTP;
TP13338* inputsTPParent;
int endRange;
_checkInCodelets13339():
darts::Codelet(){ }
_checkInCodelets13339(uint32_t dep, uint32_t res, TP13338* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP13197 * TPParent;
TP13338 * controlTPParent;
TP13338* inputsTPParent;int *j_darts13338/*OMP_PRIVATE - INPUT*/;
int **jglob_darts13338 /*OMP_SHARED_PRIVATE - INPUT*/;
int *k_darts13338/*OMP_PRIVATE - INPUT*/;
int initIteration13338;
int lastIteration13338;
int range13338;
int rangePerCodelet13338;
int minIteration13338;
int remainderRange13338;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets13338* barrierCodelets13338;
_checkInCodelets13339* checkInCodelets13339;
#if USE_SPIN_CODELETS == 0
_checkInCodelets13339* firstCodelet;
#endif
TP13338(int in_numThreads, int in_mainCodeletID, TP13197* in_TPParent, int in_initIteration, int in_lastIteration, TP13338** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP13338();
};
/*TP13380: OMPForDirective*/
class TP13380 : public ompTP
{
 public:
class _barrierCodelets13380 : public darts::Codelet
{
public:
TP13380* inputsTPParent;
_barrierCodelets13380():
darts::Codelet(){ }
_barrierCodelets13380(uint32_t dep, uint32_t res, TP13380* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations13380(int* endRange, uint32_t codeletID);
class _checkInCodelets13381:public darts::Codelet
{
public:
TP13380* myTP;
TP13380* inputsTPParent;
int endRange;
_checkInCodelets13381():
darts::Codelet(){ }
_checkInCodelets13381(uint32_t dep, uint32_t res, TP13380* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP13197 * TPParent;
TP13380 * controlTPParent;
TP13380* inputsTPParent;int *j_darts13380/*OMP_PRIVATE - INPUT*/;
int **jglob_darts13380 /*OMP_SHARED_PRIVATE - INPUT*/;
int *k_darts13380/*OMP_PRIVATE - INPUT*/;
int initIteration13380;
int lastIteration13380;
int range13380;
int rangePerCodelet13380;
int minIteration13380;
int remainderRange13380;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets13380* barrierCodelets13380;
_checkInCodelets13381* checkInCodelets13381;
#if USE_SPIN_CODELETS == 0
_checkInCodelets13381* firstCodelet;
#endif
TP13380(int in_numThreads, int in_mainCodeletID, TP13197* in_TPParent, int in_initIteration, int in_lastIteration, TP13380** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP13380();
};
/*TP13764: OMPParallelDirective*/
class TP13764 : public darts::ThreadedProcedure
{
 public:
class _barrierCodelets13764 : public darts::Codelet
{
public:
TP13764* inputsTPParent;
_barrierCodelets13764():
darts::Codelet(){ }
_barrierCodelets13764(uint32_t dep, uint32_t res, TP13764* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets13770:public darts::Codelet
{
public:
TP13764* myTP;
TP13764* inputsTPParent;
_checkInCodelets13770():
darts::Codelet(){ }
_checkInCodelets13770(uint32_t dep, uint32_t res, TP13764* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets13770 : public darts::Codelet
{
public:
TP13764* inputsTPParent;
_barrierCodelets13770():
darts::Codelet(){ }
_barrierCodelets13770(uint32_t dep, uint32_t res, TP13764* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
darts::Codelet* nextCodelet;
TP13764 * TPParent;
TP13764 * controlTPParent;
TP13764* inputsTPParent;double *eta_darts13764/*VARIABLE*/;
int *i_darts13764/*VARIABLE*/;
int *iglob_darts13764/*VARIABLE*/;
int *j_darts13764/*VARIABLE*/;
int *jglob_darts13764/*VARIABLE*/;
int *k_darts13764/*VARIABLE*/;
int *m_darts13764/*VARIABLE*/;
double *peta_darts13764/*VARIABLE*/;
double *pxi_darts13764/*VARIABLE*/;
double *pzeta_darts13764/*VARIABLE*/;
double *xi_darts13764/*VARIABLE*/;
double *zeta_darts13764/*VARIABLE*/;
TP13770** TP13770Ptr;
size_t *TP13770_alreadyLaunched;
int numTPsSet13770;
int numTPsReady13770;
size_t TPsToUse13770;
size_t codeletsPerTP13770;
size_t totalCodelets13770;
_barrierCodelets13764* barrierCodelets13764;
_checkInCodelets13770* checkInCodelets13770;
_barrierCodelets13770* barrierCodelets13770;
TP13764(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
~TP13764();
};
/*TP13770: OMPForDirective*/
class TP13770 : public ompTP
{
 public:
class _barrierCodelets13770 : public darts::Codelet
{
public:
TP13770* inputsTPParent;
_barrierCodelets13770():
darts::Codelet(){ }
_barrierCodelets13770(uint32_t dep, uint32_t res, TP13770* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations13770(int* endRange, uint32_t codeletID);
class _checkInCodelets13771:public darts::Codelet
{
public:
TP13770* myTP;
TP13770* inputsTPParent;
int endRange;
_checkInCodelets13771():
darts::Codelet(){ }
_checkInCodelets13771(uint32_t dep, uint32_t res, TP13770* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP13764 * TPParent;
TP13770 * controlTPParent;
TP13770* inputsTPParent;double **eta_darts13770 /*OMP_SHARED_PRIVATE - INPUT*/;
int *i_darts13770/*OMP_PRIVATE - INPUT*/;
int **iglob_darts13770 /*OMP_SHARED_PRIVATE - INPUT*/;
int *j_darts13770/*OMP_PRIVATE - INPUT*/;
int **jglob_darts13770 /*OMP_SHARED_PRIVATE - INPUT*/;
int *k_darts13770/*OMP_PRIVATE - INPUT*/;
int *m_darts13770/*OMP_PRIVATE - INPUT*/;
double **peta_darts13770 /*OMP_SHARED_PRIVATE - INPUT*/;
double **pxi_darts13770 /*OMP_SHARED_PRIVATE - INPUT*/;
double **pzeta_darts13770 /*OMP_SHARED_PRIVATE - INPUT*/;
double **xi_darts13770 /*OMP_SHARED_PRIVATE - INPUT*/;
double **zeta_darts13770 /*OMP_SHARED_PRIVATE - INPUT*/;
int initIteration13770;
int lastIteration13770;
int range13770;
int rangePerCodelet13770;
int minIteration13770;
int remainderRange13770;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets13770* barrierCodelets13770;
_checkInCodelets13771* checkInCodelets13771;
#if USE_SPIN_CODELETS == 0
_checkInCodelets13771* firstCodelet;
#endif
TP13770(int in_numThreads, int in_mainCodeletID, TP13764* in_TPParent, int in_initIteration, int in_lastIteration, TP13770** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP13770();
};
/*TP13896: OMPParallelDirective*/
class TP13896 : public darts::ThreadedProcedure
{
 public:
class _barrierCodelets13896 : public darts::Codelet
{
public:
TP13896* inputsTPParent;
_barrierCodelets13896():
darts::Codelet(){ }
_barrierCodelets13896(uint32_t dep, uint32_t res, TP13896* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets13898:public darts::Codelet
{
public:
TP13896* myTP;
TP13896* inputsTPParent;
_checkInCodelets13898():
darts::Codelet(){ }
_checkInCodelets13898(uint32_t dep, uint32_t res, TP13896* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets13898 : public darts::Codelet
{
public:
TP13896* inputsTPParent;
_barrierCodelets13898():
darts::Codelet(){ }
_barrierCodelets13898(uint32_t dep, uint32_t res, TP13896* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
darts::Codelet* nextCodelet;
TP13896 * TPParent;
TP13896 * controlTPParent;
TP13896* inputsTPParent;int *i_darts13896/*OMP_PRIVATE - INPUT*/;
int *j_darts13896/*OMP_PRIVATE - INPUT*/;
int *k_darts13896/*OMP_PRIVATE - INPUT*/;
int *m_darts13896/*OMP_PRIVATE - INPUT*/;
TP13898** TP13898Ptr;
size_t *TP13898_alreadyLaunched;
int numTPsSet13898;
int numTPsReady13898;
size_t TPsToUse13898;
size_t codeletsPerTP13898;
size_t totalCodelets13898;
_barrierCodelets13896* barrierCodelets13896;
_checkInCodelets13898* checkInCodelets13898;
_barrierCodelets13898* barrierCodelets13898;
TP13896(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
~TP13896();
};
/*TP13898: OMPForDirective*/
class TP13898 : public ompTP
{
 public:
class _barrierCodelets13898 : public darts::Codelet
{
public:
TP13898* inputsTPParent;
_barrierCodelets13898():
darts::Codelet(){ }
_barrierCodelets13898(uint32_t dep, uint32_t res, TP13898* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations13898(int* endRange, uint32_t codeletID);
class _checkInCodelets13899:public darts::Codelet
{
public:
TP13898* myTP;
TP13898* inputsTPParent;
int endRange;
_checkInCodelets13899():
darts::Codelet(){ }
_checkInCodelets13899(uint32_t dep, uint32_t res, TP13898* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP13896 * TPParent;
TP13898 * controlTPParent;
TP13898* inputsTPParent;int *i_darts13898/*OMP_PRIVATE - INPUT*/;
int *j_darts13898/*OMP_PRIVATE - INPUT*/;
int *k_darts13898/*OMP_PRIVATE - INPUT*/;
int *m_darts13898/*OMP_PRIVATE - INPUT*/;
int initIteration13898;
int lastIteration13898;
int range13898;
int rangePerCodelet13898;
int minIteration13898;
int remainderRange13898;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets13898* barrierCodelets13898;
_checkInCodelets13899* checkInCodelets13899;
#if USE_SPIN_CODELETS == 0
_checkInCodelets13899* firstCodelet;
#endif
TP13898(int in_numThreads, int in_mainCodeletID, TP13896* in_TPParent, int in_initIteration, int in_lastIteration, TP13898** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP13898();
};
/*TP13995: OMPParallelDirective*/
class TP13995 : public darts::ThreadedProcedure
{
 public:
class _barrierCodelets13995 : public darts::Codelet
{
public:
TP13995* inputsTPParent;
_barrierCodelets13995():
darts::Codelet(){ }
_barrierCodelets13995(uint32_t dep, uint32_t res, TP13995* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets13997:public darts::Codelet
{
public:
TP13995* myTP;
TP13995* inputsTPParent;
_checkInCodelets13997():
darts::Codelet(){ }
_checkInCodelets13997(uint32_t dep, uint32_t res, TP13995* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets13997 : public darts::Codelet
{
public:
TP13995* inputsTPParent;
_barrierCodelets13997():
darts::Codelet(){ }
_barrierCodelets13997(uint32_t dep, uint32_t res, TP13995* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets14052:public darts::Codelet
{
public:
TP13995* myTP;
TP13995* inputsTPParent;
_checkInCodelets14052():
darts::Codelet(){ }
_checkInCodelets14052(uint32_t dep, uint32_t res, TP13995* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets14055:public darts::Codelet
{
public:
TP13995* myTP;
TP13995* inputsTPParent;
_checkInCodelets14055():
darts::Codelet(){ }
_checkInCodelets14055(uint32_t dep, uint32_t res, TP13995* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets14054:public darts::Codelet
{
public:
TP13995* myTP;
TP13995* inputsTPParent;
_checkInCodelets14054():
darts::Codelet(){ }
_checkInCodelets14054(uint32_t dep, uint32_t res, TP13995* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets14058:public darts::Codelet
{
public:
TP13995* myTP;
TP13995* inputsTPParent;
_checkInCodelets14058():
darts::Codelet(){ }
_checkInCodelets14058(uint32_t dep, uint32_t res, TP13995* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets14066 : public darts::Codelet
{
public:
TP13995* inputsTPParent;
_barrierCodelets14066():
darts::Codelet(){ }
_barrierCodelets14066(uint32_t dep, uint32_t res, TP13995* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets14067:public darts::Codelet
{
public:
TP13995* myTP;
TP13995* inputsTPParent;
_checkInCodelets14067():
darts::Codelet(){ }
_checkInCodelets14067(uint32_t dep, uint32_t res, TP13995* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets14070:public darts::Codelet
{
public:
TP13995* myTP;
TP13995* inputsTPParent;
_checkInCodelets14070():
darts::Codelet(){ }
_checkInCodelets14070(uint32_t dep, uint32_t res, TP13995* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets14069:public darts::Codelet
{
public:
TP13995* myTP;
TP13995* inputsTPParent;
_checkInCodelets14069():
darts::Codelet(){ }
_checkInCodelets14069(uint32_t dep, uint32_t res, TP13995* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets14073:public darts::Codelet
{
public:
TP13995* myTP;
TP13995* inputsTPParent;
_checkInCodelets14073():
darts::Codelet(){ }
_checkInCodelets14073(uint32_t dep, uint32_t res, TP13995* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets14077 : public darts::Codelet
{
public:
TP13995* inputsTPParent;
_barrierCodelets14077():
darts::Codelet(){ }
_barrierCodelets14077(uint32_t dep, uint32_t res, TP13995* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets14078:public darts::Codelet
{
public:
TP13995* myTP;
TP13995* inputsTPParent;
_checkInCodelets14078():
darts::Codelet(){ }
_checkInCodelets14078(uint32_t dep, uint32_t res, TP13995* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets14080:public darts::Codelet
{
public:
TP13995* myTP;
TP13995* inputsTPParent;
_checkInCodelets14080():
darts::Codelet(){ }
_checkInCodelets14080(uint32_t dep, uint32_t res, TP13995* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _barrierCodelets14080 : public darts::Codelet
{
public:
TP13995* inputsTPParent;
_barrierCodelets14080():
darts::Codelet(){ }
_barrierCodelets14080(uint32_t dep, uint32_t res, TP13995* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
darts::Codelet* nextCodelet;
TP13995 * TPParent;
TP13995 * controlTPParent;
TP13995* inputsTPParent;int *i_darts13995/*OMP_PRIVATE - INPUT*/;
int *istep_darts13995/*OMP_PRIVATE - INPUT*/;
int *j_darts13995/*OMP_PRIVATE - INPUT*/;
int *k_darts13995/*OMP_PRIVATE - INPUT*/;
int *m_darts13995/*OMP_PRIVATE - INPUT*/;
double *tmp_darts13995;/*OMP_SHARED - INPUT*/
int* k_darts14064/*OMP_SHARED_PRIVATE - INPUT*/;
int* k_darts14061/*OMP_SHARED_PRIVATE - INPUT*/;
TP13997** TP13997Ptr;
size_t *TP13997_alreadyLaunched;
int numTPsSet13997;
int numTPsReady13997;
size_t TPsToUse13997;
size_t codeletsPerTP13997;
size_t totalCodelets13997;
size_t TP14052_alreadyLaunched;
unsigned int TP14054_LoopCounter;
unsigned int *TP14054_LoopCounterPerThread;
tbb::concurrent_vector<TP14054*> TP14054PtrVec;
size_t TP14067_alreadyLaunched;
unsigned int TP14069_LoopCounter;
unsigned int *TP14069_LoopCounterPerThread;
tbb::concurrent_vector<TP14069*> TP14069PtrVec;
size_t TP14078_alreadyLaunched;
TP14080** TP14080Ptr;
size_t *TP14080_alreadyLaunched;
int numTPsSet14080;
int numTPsReady14080;
size_t TPsToUse14080;
size_t codeletsPerTP14080;
size_t totalCodelets14080;
_barrierCodelets13995* barrierCodelets13995;
_checkInCodelets13997* checkInCodelets13997;
_barrierCodelets13997* barrierCodelets13997;
_checkInCodelets14052* checkInCodelets14052;
_checkInCodelets14055* checkInCodelets14055;
_checkInCodelets14054* checkInCodelets14054;
_checkInCodelets14058* checkInCodelets14058;
_barrierCodelets14066* barrierCodelets14066;
_checkInCodelets14067* checkInCodelets14067;
_checkInCodelets14070* checkInCodelets14070;
_checkInCodelets14069* checkInCodelets14069;
_checkInCodelets14073* checkInCodelets14073;
_barrierCodelets14077* barrierCodelets14077;
_checkInCodelets14078* checkInCodelets14078;
_checkInCodelets14080* checkInCodelets14080;
_barrierCodelets14080* barrierCodelets14080;
TP13995(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet, double *in_tmp);
~TP13995();
};
/*TP13997: OMPForDirective*/
class TP13997 : public ompTP
{
 public:
class _barrierCodelets13997 : public darts::Codelet
{
public:
TP13997* inputsTPParent;
_barrierCodelets13997():
darts::Codelet(){ }
_barrierCodelets13997(uint32_t dep, uint32_t res, TP13997* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations13997(int* endRange, uint32_t codeletID);
class _checkInCodelets13998:public darts::Codelet
{
public:
TP13997* myTP;
TP13997* inputsTPParent;
int endRange;
_checkInCodelets13998():
darts::Codelet(){ }
_checkInCodelets13998(uint32_t dep, uint32_t res, TP13997* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP13995 * TPParent;
TP13997 * controlTPParent;
TP13997* inputsTPParent;int *i_darts13997/*OMP_PRIVATE - INPUT*/;
int *j_darts13997/*OMP_PRIVATE - INPUT*/;
int *k_darts13997/*OMP_PRIVATE - INPUT*/;
int *m_darts13997/*OMP_PRIVATE - INPUT*/;
int initIteration13997;
int lastIteration13997;
int range13997;
int rangePerCodelet13997;
int minIteration13997;
int remainderRange13997;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets13997* barrierCodelets13997;
_checkInCodelets13998* checkInCodelets13998;
#if USE_SPIN_CODELETS == 0
_checkInCodelets13998* firstCodelet;
#endif
TP13997(int in_numThreads, int in_mainCodeletID, TP13995* in_TPParent, int in_initIteration, int in_lastIteration, TP13997** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP13997();
};
/*TP14054: ForStmt*/
class TP14054 : public ompTP
{
 public:
class _checkInCodelets14060:public darts::Codelet
{
public:
TP14054* myTP;
TP13995* inputsTPParent;
_checkInCodelets14060():
darts::Codelet(){ }
_checkInCodelets14060(uint32_t dep, uint32_t res, TP14054* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets14061:public darts::Codelet
{
public:
TP14054* myTP;
TP13995* inputsTPParent;
_checkInCodelets14061():
darts::Codelet(){ }
_checkInCodelets14061(uint32_t dep, uint32_t res, TP14054* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets14063:public darts::Codelet
{
public:
TP14054* myTP;
TP13995* inputsTPParent;
_checkInCodelets14063():
darts::Codelet(){ }
_checkInCodelets14063(uint32_t dep, uint32_t res, TP14054* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets14064:public darts::Codelet
{
public:
TP14054* myTP;
TP13995* inputsTPParent;
_checkInCodelets14064():
darts::Codelet(){ }
_checkInCodelets14064(uint32_t dep, uint32_t res, TP14054* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP13995 * TPParent;
TP14054 * controlTPParent;
TP13995* inputsTPParent;TP14054** ptrToThisTP;
TP_jacld* TP14060Ptr;
int TP14060_alreadyLaunched;
size_t TP14061_alreadyLaunched;
TP_blts* TP14063Ptr;
int TP14063_alreadyLaunched;
size_t TP14064_alreadyLaunched;
_checkInCodelets14060* checkInCodelets14060;
#if USE_SPIN_CODELETS == 0
_checkInCodelets14060* firstCodelet;
#endif
_checkInCodelets14061* checkInCodelets14061;
_checkInCodelets14063* checkInCodelets14063;
_checkInCodelets14064* checkInCodelets14064;
TP14054(int in_numThreads, int in_mainCodeletID, TP13995* in_TPParent, TP13995* in_inputsTPParent, TP14054** in_ptrToThisTP);
~TP14054();
};
/*TP14069: ForStmt*/
class TP14069 : public ompTP
{
 public:
class _checkInCodelets14075:public darts::Codelet
{
public:
TP14069* myTP;
TP13995* inputsTPParent;
_checkInCodelets14075():
darts::Codelet(){ }
_checkInCodelets14075(uint32_t dep, uint32_t res, TP14069* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
class _checkInCodelets14076:public darts::Codelet
{
public:
TP14069* myTP;
TP13995* inputsTPParent;
_checkInCodelets14076():
darts::Codelet(){ }
_checkInCodelets14076(uint32_t dep, uint32_t res, TP14069* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP13995 * TPParent;
TP14069 * controlTPParent;
TP13995* inputsTPParent;TP14069** ptrToThisTP;
TP_jacu* TP14075Ptr;
int TP14075_alreadyLaunched;
TP_buts* TP14076Ptr;
int TP14076_alreadyLaunched;
_checkInCodelets14075* checkInCodelets14075;
#if USE_SPIN_CODELETS == 0
_checkInCodelets14075* firstCodelet;
#endif
_checkInCodelets14076* checkInCodelets14076;
TP14069(int in_numThreads, int in_mainCodeletID, TP13995* in_TPParent, TP13995* in_inputsTPParent, TP14069** in_ptrToThisTP);
~TP14069();
};
/*TP14080: OMPForDirective*/
class TP14080 : public ompTP
{
 public:
class _barrierCodelets14080 : public darts::Codelet
{
public:
TP14080* inputsTPParent;
_barrierCodelets14080():
darts::Codelet(){ }
_barrierCodelets14080(uint32_t dep, uint32_t res, TP14080* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
bool requestNewRangeIterations14080(int* endRange, uint32_t codeletID);
class _checkInCodelets14081:public darts::Codelet
{
public:
TP14080* myTP;
TP14080* inputsTPParent;
int endRange;
_checkInCodelets14081():
darts::Codelet(){ }
_checkInCodelets14081(uint32_t dep, uint32_t res, TP14080* myTP, uint32_t id):
darts::Codelet(dep,res,myTP,LONGWAIT, id), myTP(myTP), inputsTPParent(myTP->inputsTPParent){ }
void fire(void);
};
TP13995 * TPParent;
TP14080 * controlTPParent;
TP14080* inputsTPParent;int *i_darts14080/*OMP_PRIVATE - INPUT*/;
int *j_darts14080/*OMP_PRIVATE - INPUT*/;
int *k_darts14080/*OMP_PRIVATE - INPUT*/;
int *m_darts14080/*OMP_PRIVATE - INPUT*/;
double *tmp_darts14080;/*OMP_SHARED - INPUT*/
int initIteration14080;
int lastIteration14080;
int range14080;
int rangePerCodelet14080;
int minIteration14080;
int remainderRange14080;
size_t readyCodelets;
int baseNumThreads;
_barrierCodelets14080* barrierCodelets14080;
_checkInCodelets14081* checkInCodelets14081;
#if USE_SPIN_CODELETS == 0
_checkInCodelets14081* firstCodelet;
#endif
TP14080(int in_numThreads, int in_mainCodeletID, TP13995* in_TPParent, int in_initIteration, int in_lastIteration, double *in_tmp, TP14080** in_ptrToThisTP);
void inline dispatchCodelet(size_t codeletID);
~TP14080();
};
#endif
