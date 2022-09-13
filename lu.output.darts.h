#ifndef _lu_output_darts_h_
#define _lu_output_darts_h_
#ifndef __DARTS_
#define __DARTS_
#endif
#include "applu.h"
#include "darts.h"
#include "npb-C.h"
#include "ompTP.h"
#include "tbb/concurrent_vector.h"
#include "utils.h"
#include <limits.h>
#include <math.h>
#include <mutex>
#include <numa.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
int main(int argc, char** argv);
class TP192;
class TP1;
typedef TP1 TP_blts;
class TP248;
/*Number of TPs to be used for the OMPFor in region TP248*/
#define NUMTPS248 NUMTPS
class TP2;
typedef TP2 TP_buts;
class TP1213;
/*Number of TPs to be used for the OMPFor in region TP1213*/
#define NUMTPS1213 NUMTPS
class TP2199;
class TP2218;
/*Number of TPs to be used for the OMPFor in region TP2218*/
#define NUMTPS2218 NUMTPS
class TP2269;
/*Number of TPs to be used for the OMPFor in region TP2269*/
#define NUMTPS2269 NUMTPS
class TP2404;
/*Number of TPs to be used for the OMPFor in region TP2404*/
#define NUMTPS2404 NUMTPS
class TP2553;
/*Number of TPs to be used for the OMPFor in region TP2553*/
#define NUMTPS2553 NUMTPS
class TP3191;
/*Number of TPs to be used for the OMPFor in region TP3191*/
#define NUMTPS3191 NUMTPS
class TP3340;
/*Number of TPs to be used for the OMPFor in region TP3340*/
#define NUMTPS3340 NUMTPS
class TP3975;
/*Number of TPs to be used for the OMPFor in region TP3975*/
#define NUMTPS3975 NUMTPS
class TP7;
typedef TP7 TP_jacld;
class TP4901;
/*Number of TPs to be used for the OMPFor in region TP4901*/
#define NUMTPS4901 NUMTPS
class TP8;
typedef TP8 TP_jacu;
class TP7420;
/*Number of TPs to be used for the OMPFor in region TP7420*/
#define NUMTPS7420 NUMTPS
class TP9871;
class TP9889;
/*Number of TPs to be used for the OMPFor in region TP9889*/
#define NUMTPS9889 NUMTPS
class TP10788;
class TP10803;
/*Number of TPs to be used for the OMPFor in region TP10803*/
#define NUMTPS10803 NUMTPS
class TP10863;
/*Number of TPs to be used for the OMPFor in region TP10863*/
#define NUMTPS10863 NUMTPS
class TP11012;
/*Number of TPs to be used for the OMPFor in region TP11012*/
#define NUMTPS11012 NUMTPS
class TP11652;
/*Number of TPs to be used for the OMPFor in region TP11652*/
#define NUMTPS11652 NUMTPS
class TP11801;
/*Number of TPs to be used for the OMPFor in region TP11801*/
#define NUMTPS11801 NUMTPS
class TP12438;
/*Number of TPs to be used for the OMPFor in region TP12438*/
#define NUMTPS12438 NUMTPS
class TP13189;
class TP13193;
/*Number of TPs to be used for the OMPFor in region TP13193*/
#define NUMTPS13193 NUMTPS
class TP13244;
/*Number of TPs to be used for the OMPFor in region TP13244*/
#define NUMTPS13244 NUMTPS
class TP13286;
/*Number of TPs to be used for the OMPFor in region TP13286*/
#define NUMTPS13286 NUMTPS
class TP13330;
/*Number of TPs to be used for the OMPFor in region TP13330*/
#define NUMTPS13330 NUMTPS
class TP13372;
/*Number of TPs to be used for the OMPFor in region TP13372*/
#define NUMTPS13372 NUMTPS
class TP13756;
class TP13762;
/*Number of TPs to be used for the OMPFor in region TP13762*/
#define NUMTPS13762 NUMTPS
class TP13888;
class TP13890;
/*Number of TPs to be used for the OMPFor in region TP13890*/
#define NUMTPS13890 NUMTPS
class TP13987;
class TP13989;
/*Number of TPs to be used for the OMPFor in region TP13989*/
#define NUMTPS13989 NUMTPS
class TP14044;
class TP14053;
class TP14062;
/*Number of TPs to be used for the OMPFor in region TP14062*/
#define NUMTPS14062 NUMTPS
extern int DARTS_CODELETS_MULT;
extern int NUMTPS;
extern size_t numOfCUs;
extern darts::Codelet* RuntimeFinalCodelet;
extern darts::ThreadAffinity* affin;
extern bool affinMaskRes;
extern darts::Runtime* myDARTSRuntime;
extern std::vector<std::vector<void*>> threadFunctionStack;
extern size_t ompNumThreads;
extern int ompSchedulePolicy;
extern int ompScheduleChunk;
extern void omp_set_num_threads(unsigned long numThreadsToSet);
extern int omp_get_num_threads();
extern int omp_get_max_threads();
extern int omp_get_num_procs();
extern double omp_get_wtime();
extern void omp_init_lock(omp_lock_t* lock);
extern void omp_destroy_lock(omp_lock_t* lock);
extern void omp_set_lock(omp_lock_t* lock);
extern void omp_unset_lock(omp_lock_t* lock);
/*TP192: OMPParallelDirective*/
class TP192 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets192 : public darts::Codelet {
    public:
        TP192* inputsTPParent;
        _barrierCodelets192()
            : darts::Codelet()
        {
        }
        _barrierCodelets192(uint32_t dep, uint32_t res, TP192* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets194 : public darts::Codelet {
    public:
        TP192* myTP;
        TP192* inputsTPParent;
        _checkInCodelets194()
            : darts::Codelet()
        {
        }
        _checkInCodelets194(uint32_t dep, uint32_t res, TP192* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP192* TPParent;
    TP192* controlTPParent;
    TP192* inputsTPParent;
    int* nthreads_darts192; /*OMP_SHARED - INPUT*/
    int* nthreads_darts194; /*OMP_SHARED - INPUT*/
    size_t TP194_alreadyLaunched;
    _barrierCodelets192* barrierCodelets192;
    _checkInCodelets194* checkInCodelets194;
    TP192(
        int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet, int* in_nthreads);
    ~TP192();
};
/*TP1: blts*/
class TP1 : public ompTP {
public:
    class _checkInCodelets248 : public darts::Codelet {
    public:
        TP1* myTP;
        TP1* inputsTPParent;
        _checkInCodelets248()
            : darts::Codelet()
        {
        }
        _checkInCodelets248(uint32_t dep, uint32_t res, TP1* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets248 : public darts::Codelet {
    public:
        TP1* inputsTPParent;
        _barrierCodelets248()
            : darts::Codelet()
        {
        }
        _barrierCodelets248(uint32_t dep, uint32_t res, TP1* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets352 : public darts::Codelet {
    public:
        TP1* myTP;
        TP1* inputsTPParent;
        _checkInCodelets352()
            : darts::Codelet()
        {
        }
        _checkInCodelets352(uint32_t dep, uint32_t res, TP1* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP1** ptrToThisFunctionTP;
    TP1* inputsTPParent;
    TP1* controlTPParent;
    darts::Codelet** nextCodeletsblts;
    darts::Codelet** nextSyncCodeletsblts;
    int* nx_darts1 /*VARIABLE*/;
    int* ny_darts1 /*VARIABLE*/;
    int* nz_darts1 /*VARIABLE*/;
    int* k_darts1 /*VARIABLE*/;
    double* omega_darts1 /*VARIABLE*/;
    int* ist_darts1 /*VARIABLE*/;
    int* iend_darts1 /*VARIABLE*/;
    int* jst_darts1 /*VARIABLE*/;
    int* jend_darts1 /*VARIABLE*/;
    int* nx0_darts1 /*VARIABLE*/;
    int* ny0_darts1 /*VARIABLE*/;
    int* i_darts1 /*VARIABLE*/;
    int* j_darts1 /*VARIABLE*/;
    int* m_darts1 /*VARIABLE*/;
    double* tmp_darts1 /*VARIABLE*/;
    double* tmp1_darts1 /*VARIABLE*/;
    int i_darts352;
    int* iend_darts352 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* ist_darts352 /*OMP_SHARED_PRIVATE - INPUT*/;
    int j_darts352;
    int* jend_darts352 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* jst_darts352 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts352 /*OMP_SHARED_PRIVATE - INPUT*/;
    int m_darts352;
    double* omega_darts352 /*OMP_SHARED_PRIVATE - INPUT*/;
    double* tmp_darts352 /*OMP_SHARED_PRIVATE - INPUT*/;
    double* tmp1_darts352 /*OMP_SHARED_PRIVATE - INPUT*/;
    TP248** TP248Ptr;
    size_t* TP248_alreadyLaunched;
    int numTPsSet248;
    int numTPsReady248;
    size_t TPsToUse248;
    size_t codeletsPerTP248;
    size_t totalCodelets248;
    size_t TP352_alreadyLaunched;
    _checkInCodelets248* checkInCodelets248;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets248* firstCodelet;
#endif
    _barrierCodelets248* barrierCodelets248;
    _checkInCodelets352* checkInCodelets352;
    TP1(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
        darts::Codelet* in_mainSyncCodelet, TP1** in_ptrToThisFunctionTP, int in_nx, int in_ny,
        int in_nz, int in_k, double in_omega, int in_ist, int in_iend, int in_jst, int in_jend,
        int in_nx0, int in_ny0);
    ~TP1();
    void setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, int in_ist,
        int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID);
};
/*TP248: OMPForDirective*/
class TP248 : public ompTP {
public:
    class _barrierCodelets248 : public darts::Codelet {
    public:
        TP248* inputsTPParent;
        _barrierCodelets248()
            : darts::Codelet()
        {
        }
        _barrierCodelets248(uint32_t dep, uint32_t res, TP248* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations248(int* endRange, uint32_t codeletID);
    class _checkInCodelets249 : public darts::Codelet {
    public:
        TP248* myTP;
        TP248* inputsTPParent;
        int endRange;
        _checkInCodelets249()
            : darts::Codelet()
        {
        }
        _checkInCodelets249(uint32_t dep, uint32_t res, TP248* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP1* TPParent;
    TP248* controlTPParent;
    TP248* inputsTPParent;
    int* i_darts248 /*OMP_PRIVATE - INPUT*/;
    int** iend_darts248 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** ist_darts248 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts248 /*OMP_PRIVATE - INPUT*/;
    int** jend_darts248 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** jst_darts248 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** k_darts248 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* m_darts248 /*OMP_PRIVATE - INPUT*/;
    double** omega_darts248 /*OMP_SHARED_PRIVATE - INPUT*/;
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
    TP248(int in_numThreads, int in_mainCodeletID, TP1* in_TPParent, int in_initIteration,
        int in_lastIteration, TP248** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP248();
};
/*TP2: buts*/
class TP2 : public ompTP {
public:
    class _checkInCodelets1213 : public darts::Codelet {
    public:
        TP2* myTP;
        TP2* inputsTPParent;
        _checkInCodelets1213()
            : darts::Codelet()
        {
        }
        _checkInCodelets1213(uint32_t dep, uint32_t res, TP2* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets1213 : public darts::Codelet {
    public:
        TP2* inputsTPParent;
        _barrierCodelets1213()
            : darts::Codelet()
        {
        }
        _barrierCodelets1213(uint32_t dep, uint32_t res, TP2* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets1315 : public darts::Codelet {
    public:
        TP2* myTP;
        TP2* inputsTPParent;
        _checkInCodelets1315()
            : darts::Codelet()
        {
        }
        _checkInCodelets1315(uint32_t dep, uint32_t res, TP2* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2** ptrToThisFunctionTP;
    TP2* inputsTPParent;
    TP2* controlTPParent;
    darts::Codelet** nextCodeletsbuts;
    darts::Codelet** nextSyncCodeletsbuts;
    int* nx_darts2 /*VARIABLE*/;
    int* ny_darts2 /*VARIABLE*/;
    int* nz_darts2 /*VARIABLE*/;
    int* k_darts2 /*VARIABLE*/;
    double* omega_darts2 /*VARIABLE*/;
    int* ist_darts2 /*VARIABLE*/;
    int* iend_darts2 /*VARIABLE*/;
    int* jst_darts2 /*VARIABLE*/;
    int* jend_darts2 /*VARIABLE*/;
    int* nx0_darts2 /*VARIABLE*/;
    int* ny0_darts2 /*VARIABLE*/;
    int* i_darts2 /*VARIABLE*/;
    int* j_darts2 /*VARIABLE*/;
    int* m_darts2 /*VARIABLE*/;
    double* tmp_darts2 /*VARIABLE*/;
    double* tmp1_darts2 /*VARIABLE*/;
    int i_darts1315;
    int* iend_darts1315 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* ist_darts1315 /*OMP_SHARED_PRIVATE - INPUT*/;
    int j_darts1315;
    int* jend_darts1315 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* jst_darts1315 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts1315 /*OMP_SHARED_PRIVATE - INPUT*/;
    int m_darts1315;
    double* omega_darts1315 /*OMP_SHARED_PRIVATE - INPUT*/;
    double* tmp_darts1315 /*OMP_SHARED_PRIVATE - INPUT*/;
    double* tmp1_darts1315 /*OMP_SHARED_PRIVATE - INPUT*/;
    TP1213** TP1213Ptr;
    size_t* TP1213_alreadyLaunched;
    int numTPsSet1213;
    int numTPsReady1213;
    size_t TPsToUse1213;
    size_t codeletsPerTP1213;
    size_t totalCodelets1213;
    size_t TP1315_alreadyLaunched;
    _checkInCodelets1213* checkInCodelets1213;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets1213* firstCodelet;
#endif
    _barrierCodelets1213* barrierCodelets1213;
    _checkInCodelets1315* checkInCodelets1315;
    TP2(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
        darts::Codelet* in_mainSyncCodelet, TP2** in_ptrToThisFunctionTP, int in_nx, int in_ny,
        int in_nz, int in_k, double in_omega, int in_ist, int in_iend, int in_jst, int in_jend,
        int in_nx0, int in_ny0);
    ~TP2();
    void setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, int in_ist,
        int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID);
};
/*TP1213: OMPForDirective*/
class TP1213 : public ompTP {
public:
    class _barrierCodelets1213 : public darts::Codelet {
    public:
        TP1213* inputsTPParent;
        _barrierCodelets1213()
            : darts::Codelet()
        {
        }
        _barrierCodelets1213(uint32_t dep, uint32_t res, TP1213* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations1213(int* endRange, uint32_t codeletID);
    class _checkInCodelets1214 : public darts::Codelet {
    public:
        TP1213* myTP;
        TP1213* inputsTPParent;
        int endRange;
        _checkInCodelets1214()
            : darts::Codelet()
        {
        }
        _checkInCodelets1214(uint32_t dep, uint32_t res, TP1213* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2* TPParent;
    TP1213* controlTPParent;
    TP1213* inputsTPParent;
    int* i_darts1213 /*OMP_PRIVATE - INPUT*/;
    int** iend_darts1213 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** ist_darts1213 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts1213 /*OMP_PRIVATE - INPUT*/;
    int** jend_darts1213 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** jst_darts1213 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** k_darts1213 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* m_darts1213 /*OMP_PRIVATE - INPUT*/;
    double** omega_darts1213 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration1213;
    int lastIteration1213;
    int range1213;
    int rangePerCodelet1213;
    int minIteration1213;
    int remainderRange1213;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets1213* barrierCodelets1213;
    _checkInCodelets1214* checkInCodelets1214;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets1214* firstCodelet;
#endif
    TP1213(int in_numThreads, int in_mainCodeletID, TP2* in_TPParent, int in_initIteration,
        int in_lastIteration, TP1213** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP1213();
};
/*TP2199: OMPParallelDirective*/
class TP2199 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets2199 : public darts::Codelet {
    public:
        TP2199* inputsTPParent;
        _barrierCodelets2199()
            : darts::Codelet()
        {
        }
        _barrierCodelets2199(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2201 : public darts::Codelet {
    public:
        TP2199* myTP;
        TP2199* inputsTPParent;
        _checkInCodelets2201()
            : darts::Codelet()
        {
        }
        _checkInCodelets2201(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2218 : public darts::Codelet {
    public:
        TP2199* myTP;
        TP2199* inputsTPParent;
        _checkInCodelets2218()
            : darts::Codelet()
        {
        }
        _checkInCodelets2218(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets2218 : public darts::Codelet {
    public:
        TP2199* inputsTPParent;
        _barrierCodelets2218()
            : darts::Codelet()
        {
        }
        _barrierCodelets2218(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2269 : public darts::Codelet {
    public:
        TP2199* myTP;
        TP2199* inputsTPParent;
        _checkInCodelets2269()
            : darts::Codelet()
        {
        }
        _checkInCodelets2269(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets2269 : public darts::Codelet {
    public:
        TP2199* inputsTPParent;
        _barrierCodelets2269()
            : darts::Codelet()
        {
        }
        _barrierCodelets2269(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2401 : public darts::Codelet {
    public:
        TP2199* myTP;
        TP2199* inputsTPParent;
        _checkInCodelets2401()
            : darts::Codelet()
        {
        }
        _checkInCodelets2401(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2404 : public darts::Codelet {
    public:
        TP2199* myTP;
        TP2199* inputsTPParent;
        _checkInCodelets2404()
            : darts::Codelet()
        {
        }
        _checkInCodelets2404(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets2404 : public darts::Codelet {
    public:
        TP2199* inputsTPParent;
        _barrierCodelets2404()
            : darts::Codelet()
        {
        }
        _barrierCodelets2404(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2553 : public darts::Codelet {
    public:
        TP2199* myTP;
        TP2199* inputsTPParent;
        _checkInCodelets2553()
            : darts::Codelet()
        {
        }
        _checkInCodelets2553(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets2553 : public darts::Codelet {
    public:
        TP2199* inputsTPParent;
        _barrierCodelets2553()
            : darts::Codelet()
        {
        }
        _barrierCodelets2553(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets3188 : public darts::Codelet {
    public:
        TP2199* myTP;
        TP2199* inputsTPParent;
        _checkInCodelets3188()
            : darts::Codelet()
        {
        }
        _checkInCodelets3188(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets3191 : public darts::Codelet {
    public:
        TP2199* myTP;
        TP2199* inputsTPParent;
        _checkInCodelets3191()
            : darts::Codelet()
        {
        }
        _checkInCodelets3191(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets3191 : public darts::Codelet {
    public:
        TP2199* inputsTPParent;
        _barrierCodelets3191()
            : darts::Codelet()
        {
        }
        _barrierCodelets3191(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets3340 : public darts::Codelet {
    public:
        TP2199* myTP;
        TP2199* inputsTPParent;
        _checkInCodelets3340()
            : darts::Codelet()
        {
        }
        _checkInCodelets3340(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets3340 : public darts::Codelet {
    public:
        TP2199* inputsTPParent;
        _barrierCodelets3340()
            : darts::Codelet()
        {
        }
        _barrierCodelets3340(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets3975 : public darts::Codelet {
    public:
        TP2199* myTP;
        TP2199* inputsTPParent;
        _checkInCodelets3975()
            : darts::Codelet()
        {
        }
        _checkInCodelets3975(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets3975 : public darts::Codelet {
    public:
        TP2199* inputsTPParent;
        _barrierCodelets3975()
            : darts::Codelet()
        {
        }
        _barrierCodelets3975(uint32_t dep, uint32_t res, TP2199* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP2199* TPParent;
    TP2199* controlTPParent;
    TP2199* inputsTPParent;
    int* L1_darts2199 /*VARIABLE*/;
    int* L2_darts2199 /*VARIABLE*/;
    double* dsspm_darts2199 /*VARIABLE*/;
    double* eta_darts2199 /*VARIABLE*/;
    int* i_darts2199 /*VARIABLE*/;
    int* iend1_darts2199 /*VARIABLE*/;
    int* iglob_darts2199 /*VARIABLE*/;
    int* ist1_darts2199 /*VARIABLE*/;
    int* j_darts2199 /*VARIABLE*/;
    int* jend1_darts2199 /*VARIABLE*/;
    int* jglob_darts2199 /*VARIABLE*/;
    int* jst1_darts2199 /*VARIABLE*/;
    int* k_darts2199 /*VARIABLE*/;
    int* m_darts2199 /*VARIABLE*/;
    double* q_darts2199 /*VARIABLE*/;
    double* tmp_darts2199 /*VARIABLE*/;
    double* u21_darts2199 /*VARIABLE*/;
    double* u21i_darts2199 /*VARIABLE*/;
    double* u21im1_darts2199 /*VARIABLE*/;
    double* u21j_darts2199 /*VARIABLE*/;
    double* u21jm1_darts2199 /*VARIABLE*/;
    double* u21k_darts2199 /*VARIABLE*/;
    double* u21km1_darts2199 /*VARIABLE*/;
    double* u31_darts2199 /*VARIABLE*/;
    double* u31i_darts2199 /*VARIABLE*/;
    double* u31im1_darts2199 /*VARIABLE*/;
    double* u31j_darts2199 /*VARIABLE*/;
    double* u31jm1_darts2199 /*VARIABLE*/;
    double* u31k_darts2199 /*VARIABLE*/;
    double* u31km1_darts2199 /*VARIABLE*/;
    double* u41_darts2199 /*VARIABLE*/;
    double* u41i_darts2199 /*VARIABLE*/;
    double* u41im1_darts2199 /*VARIABLE*/;
    double* u41j_darts2199 /*VARIABLE*/;
    double* u41jm1_darts2199 /*VARIABLE*/;
    double* u41k_darts2199 /*VARIABLE*/;
    double* u41km1_darts2199 /*VARIABLE*/;
    double* u51i_darts2199 /*VARIABLE*/;
    double* u51im1_darts2199 /*VARIABLE*/;
    double* u51j_darts2199 /*VARIABLE*/;
    double* u51jm1_darts2199 /*VARIABLE*/;
    double* u51k_darts2199 /*VARIABLE*/;
    double* u51km1_darts2199 /*VARIABLE*/;
    double* xi_darts2199 /*VARIABLE*/;
    double* zeta_darts2199 /*VARIABLE*/;
    TP2218** TP2218Ptr;
    size_t* TP2218_alreadyLaunched;
    int numTPsSet2218;
    int numTPsReady2218;
    size_t TPsToUse2218;
    size_t codeletsPerTP2218;
    size_t totalCodelets2218;
    TP2269** TP2269Ptr;
    size_t* TP2269_alreadyLaunched;
    int numTPsSet2269;
    int numTPsReady2269;
    size_t TPsToUse2269;
    size_t codeletsPerTP2269;
    size_t totalCodelets2269;
    TP2404** TP2404Ptr;
    size_t* TP2404_alreadyLaunched;
    int numTPsSet2404;
    int numTPsReady2404;
    size_t TPsToUse2404;
    size_t codeletsPerTP2404;
    size_t totalCodelets2404;
    TP2553** TP2553Ptr;
    size_t* TP2553_alreadyLaunched;
    int numTPsSet2553;
    int numTPsReady2553;
    size_t TPsToUse2553;
    size_t codeletsPerTP2553;
    size_t totalCodelets2553;
    TP3191** TP3191Ptr;
    size_t* TP3191_alreadyLaunched;
    int numTPsSet3191;
    int numTPsReady3191;
    size_t TPsToUse3191;
    size_t codeletsPerTP3191;
    size_t totalCodelets3191;
    TP3340** TP3340Ptr;
    size_t* TP3340_alreadyLaunched;
    int numTPsSet3340;
    int numTPsReady3340;
    size_t TPsToUse3340;
    size_t codeletsPerTP3340;
    size_t totalCodelets3340;
    TP3975** TP3975Ptr;
    size_t* TP3975_alreadyLaunched;
    int numTPsSet3975;
    int numTPsReady3975;
    size_t TPsToUse3975;
    size_t codeletsPerTP3975;
    size_t totalCodelets3975;
    _barrierCodelets2199* barrierCodelets2199;
    _checkInCodelets2201* checkInCodelets2201;
    _checkInCodelets2218* checkInCodelets2218;
    _barrierCodelets2218* barrierCodelets2218;
    _checkInCodelets2269* checkInCodelets2269;
    _barrierCodelets2269* barrierCodelets2269;
    _checkInCodelets2401* checkInCodelets2401;
    _checkInCodelets2404* checkInCodelets2404;
    _barrierCodelets2404* barrierCodelets2404;
    _checkInCodelets2553* checkInCodelets2553;
    _barrierCodelets2553* barrierCodelets2553;
    _checkInCodelets3188* checkInCodelets3188;
    _checkInCodelets3191* checkInCodelets3191;
    _barrierCodelets3191* barrierCodelets3191;
    _checkInCodelets3340* checkInCodelets3340;
    _barrierCodelets3340* barrierCodelets3340;
    _checkInCodelets3975* checkInCodelets3975;
    _barrierCodelets3975* barrierCodelets3975;
    TP2199(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
    ~TP2199();
};
/*TP2218: OMPForDirective*/
class TP2218 : public ompTP {
public:
    class _barrierCodelets2218 : public darts::Codelet {
    public:
        TP2218* inputsTPParent;
        _barrierCodelets2218()
            : darts::Codelet()
        {
        }
        _barrierCodelets2218(uint32_t dep, uint32_t res, TP2218* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations2218(int* endRange, uint32_t codeletID);
    class _checkInCodelets2219 : public darts::Codelet {
    public:
        TP2218* myTP;
        TP2218* inputsTPParent;
        int endRange;
        _checkInCodelets2219()
            : darts::Codelet()
        {
        }
        _checkInCodelets2219(uint32_t dep, uint32_t res, TP2218* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2199* TPParent;
    TP2218* controlTPParent;
    TP2218* inputsTPParent;
    int* i_darts2218 /*OMP_PRIVATE - INPUT*/;
    int* j_darts2218 /*OMP_PRIVATE - INPUT*/;
    int* k_darts2218 /*OMP_PRIVATE - INPUT*/;
    int* m_darts2218 /*OMP_PRIVATE - INPUT*/;
    int initIteration2218;
    int lastIteration2218;
    int range2218;
    int rangePerCodelet2218;
    int minIteration2218;
    int remainderRange2218;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets2218* barrierCodelets2218;
    _checkInCodelets2219* checkInCodelets2219;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets2219* firstCodelet;
#endif
    TP2218(int in_numThreads, int in_mainCodeletID, TP2199* in_TPParent, int in_initIteration,
        int in_lastIteration, TP2218** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP2218();
};
/*TP2269: OMPForDirective*/
class TP2269 : public ompTP {
public:
    class _barrierCodelets2269 : public darts::Codelet {
    public:
        TP2269* inputsTPParent;
        _barrierCodelets2269()
            : darts::Codelet()
        {
        }
        _barrierCodelets2269(uint32_t dep, uint32_t res, TP2269* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations2269(int* endRange, uint32_t codeletID);
    class _checkInCodelets2270 : public darts::Codelet {
    public:
        TP2269* myTP;
        TP2269* inputsTPParent;
        int endRange;
        _checkInCodelets2270()
            : darts::Codelet()
        {
        }
        _checkInCodelets2270(uint32_t dep, uint32_t res, TP2269* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2199* TPParent;
    TP2269* controlTPParent;
    TP2269* inputsTPParent;
    double** eta_darts2269 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts2269 /*OMP_PRIVATE - INPUT*/;
    int** iglob_darts2269 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts2269 /*OMP_PRIVATE - INPUT*/;
    int** jglob_darts2269 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts2269 /*OMP_PRIVATE - INPUT*/;
    int* m_darts2269 /*OMP_PRIVATE - INPUT*/;
    double** xi_darts2269 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** zeta_darts2269 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration2269;
    int lastIteration2269;
    int range2269;
    int rangePerCodelet2269;
    int minIteration2269;
    int remainderRange2269;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets2269* barrierCodelets2269;
    _checkInCodelets2270* checkInCodelets2270;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets2270* firstCodelet;
#endif
    TP2269(int in_numThreads, int in_mainCodeletID, TP2199* in_TPParent, int in_initIteration,
        int in_lastIteration, TP2269** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP2269();
};
/*TP2404: OMPForDirective*/
class TP2404 : public ompTP {
public:
    class _barrierCodelets2404 : public darts::Codelet {
    public:
        TP2404* inputsTPParent;
        _barrierCodelets2404()
            : darts::Codelet()
        {
        }
        _barrierCodelets2404(uint32_t dep, uint32_t res, TP2404* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations2404(int* endRange, uint32_t codeletID);
    class _checkInCodelets2405 : public darts::Codelet {
    public:
        TP2404* myTP;
        TP2404* inputsTPParent;
        int endRange;
        _checkInCodelets2405()
            : darts::Codelet()
        {
        }
        _checkInCodelets2405(uint32_t dep, uint32_t res, TP2404* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2199* TPParent;
    TP2404* controlTPParent;
    TP2404* inputsTPParent;
    int** L1_darts2404 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** L2_darts2404 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts2404 /*OMP_PRIVATE - INPUT*/;
    int* j_darts2404 /*OMP_PRIVATE - INPUT*/;
    int* k_darts2404 /*OMP_PRIVATE - INPUT*/;
    double** q_darts2404 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21_darts2404 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration2404;
    int lastIteration2404;
    int range2404;
    int rangePerCodelet2404;
    int minIteration2404;
    int remainderRange2404;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets2404* barrierCodelets2404;
    _checkInCodelets2405* checkInCodelets2405;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets2405* firstCodelet;
#endif
    TP2404(int in_numThreads, int in_mainCodeletID, TP2199* in_TPParent, int in_initIteration,
        int in_lastIteration, TP2404** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP2404();
};
/*TP2553: OMPForDirective*/
class TP2553 : public ompTP {
public:
    class _barrierCodelets2553 : public darts::Codelet {
    public:
        TP2553* inputsTPParent;
        _barrierCodelets2553()
            : darts::Codelet()
        {
        }
        _barrierCodelets2553(uint32_t dep, uint32_t res, TP2553* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations2553(int* endRange, uint32_t codeletID);
    class _checkInCodelets2554 : public darts::Codelet {
    public:
        TP2553* myTP;
        TP2553* inputsTPParent;
        int endRange;
        _checkInCodelets2554()
            : darts::Codelet()
        {
        }
        _checkInCodelets2554(uint32_t dep, uint32_t res, TP2553* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2199* TPParent;
    TP2553* controlTPParent;
    TP2553* inputsTPParent;
    int** L2_darts2553 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** dsspm_darts2553 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts2553 /*OMP_PRIVATE - INPUT*/;
    int** iend1_darts2553 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** ist1_darts2553 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts2553 /*OMP_PRIVATE - INPUT*/;
    int* k_darts2553 /*OMP_PRIVATE - INPUT*/;
    int* m_darts2553 /*OMP_PRIVATE - INPUT*/;
    double** tmp_darts2553 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21i_darts2553 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21im1_darts2553 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31i_darts2553 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31im1_darts2553 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41i_darts2553 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41im1_darts2553 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51i_darts2553 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51im1_darts2553 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration2553;
    int lastIteration2553;
    int range2553;
    int rangePerCodelet2553;
    int minIteration2553;
    int remainderRange2553;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets2553* barrierCodelets2553;
    _checkInCodelets2554* checkInCodelets2554;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets2554* firstCodelet;
#endif
    TP2553(int in_numThreads, int in_mainCodeletID, TP2199* in_TPParent, int in_initIteration,
        int in_lastIteration, TP2553** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP2553();
};
/*TP3191: OMPForDirective*/
class TP3191 : public ompTP {
public:
    class _barrierCodelets3191 : public darts::Codelet {
    public:
        TP3191* inputsTPParent;
        _barrierCodelets3191()
            : darts::Codelet()
        {
        }
        _barrierCodelets3191(uint32_t dep, uint32_t res, TP3191* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations3191(int* endRange, uint32_t codeletID);
    class _checkInCodelets3192 : public darts::Codelet {
    public:
        TP3191* myTP;
        TP3191* inputsTPParent;
        int endRange;
        _checkInCodelets3192()
            : darts::Codelet()
        {
        }
        _checkInCodelets3192(uint32_t dep, uint32_t res, TP3191* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2199* TPParent;
    TP3191* controlTPParent;
    TP3191* inputsTPParent;
    int** L1_darts3191 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** L2_darts3191 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts3191 /*OMP_PRIVATE - INPUT*/;
    int* j_darts3191 /*OMP_PRIVATE - INPUT*/;
    int* k_darts3191 /*OMP_PRIVATE - INPUT*/;
    double** q_darts3191 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31_darts3191 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration3191;
    int lastIteration3191;
    int range3191;
    int rangePerCodelet3191;
    int minIteration3191;
    int remainderRange3191;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets3191* barrierCodelets3191;
    _checkInCodelets3192* checkInCodelets3192;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets3192* firstCodelet;
#endif
    TP3191(int in_numThreads, int in_mainCodeletID, TP2199* in_TPParent, int in_initIteration,
        int in_lastIteration, TP3191** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP3191();
};
/*TP3340: OMPForDirective*/
class TP3340 : public ompTP {
public:
    class _barrierCodelets3340 : public darts::Codelet {
    public:
        TP3340* inputsTPParent;
        _barrierCodelets3340()
            : darts::Codelet()
        {
        }
        _barrierCodelets3340(uint32_t dep, uint32_t res, TP3340* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations3340(int* endRange, uint32_t codeletID);
    class _checkInCodelets3341 : public darts::Codelet {
    public:
        TP3340* myTP;
        TP3340* inputsTPParent;
        int endRange;
        _checkInCodelets3341()
            : darts::Codelet()
        {
        }
        _checkInCodelets3341(uint32_t dep, uint32_t res, TP3340* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2199* TPParent;
    TP3340* controlTPParent;
    TP3340* inputsTPParent;
    int** L2_darts3340 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** dsspm_darts3340 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts3340 /*OMP_PRIVATE - INPUT*/;
    int* j_darts3340 /*OMP_PRIVATE - INPUT*/;
    int** jend1_darts3340 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** jst1_darts3340 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts3340 /*OMP_PRIVATE - INPUT*/;
    int* m_darts3340 /*OMP_PRIVATE - INPUT*/;
    double** tmp_darts3340 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21j_darts3340 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21jm1_darts3340 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31j_darts3340 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31jm1_darts3340 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41j_darts3340 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41jm1_darts3340 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51j_darts3340 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51jm1_darts3340 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration3340;
    int lastIteration3340;
    int range3340;
    int rangePerCodelet3340;
    int minIteration3340;
    int remainderRange3340;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets3340* barrierCodelets3340;
    _checkInCodelets3341* checkInCodelets3341;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets3341* firstCodelet;
#endif
    TP3340(int in_numThreads, int in_mainCodeletID, TP2199* in_TPParent, int in_initIteration,
        int in_lastIteration, TP3340** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP3340();
};
/*TP3975: OMPForDirective*/
class TP3975 : public ompTP {
public:
    class _barrierCodelets3975 : public darts::Codelet {
    public:
        TP3975* inputsTPParent;
        _barrierCodelets3975()
            : darts::Codelet()
        {
        }
        _barrierCodelets3975(uint32_t dep, uint32_t res, TP3975* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations3975(int* endRange, uint32_t codeletID);
    class _checkInCodelets3976 : public darts::Codelet {
    public:
        TP3975* myTP;
        TP3975* inputsTPParent;
        int endRange;
        _checkInCodelets3976()
            : darts::Codelet()
        {
        }
        _checkInCodelets3976(uint32_t dep, uint32_t res, TP3975* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2199* TPParent;
    TP3975* controlTPParent;
    TP3975* inputsTPParent;
    double** dsspm_darts3975 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts3975 /*OMP_PRIVATE - INPUT*/;
    int* j_darts3975 /*OMP_PRIVATE - INPUT*/;
    int* k_darts3975 /*OMP_PRIVATE - INPUT*/;
    int* m_darts3975 /*OMP_PRIVATE - INPUT*/;
    double** q_darts3975 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp_darts3975 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21k_darts3975 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21km1_darts3975 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31k_darts3975 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31km1_darts3975 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41_darts3975 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41k_darts3975 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41km1_darts3975 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51k_darts3975 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51km1_darts3975 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration3975;
    int lastIteration3975;
    int range3975;
    int rangePerCodelet3975;
    int minIteration3975;
    int remainderRange3975;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets3975* barrierCodelets3975;
    _checkInCodelets3976* checkInCodelets3976;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets3976* firstCodelet;
#endif
    TP3975(int in_numThreads, int in_mainCodeletID, TP2199* in_TPParent, int in_initIteration,
        int in_lastIteration, TP3975** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP3975();
};
/*TP7: jacld*/
class TP7 : public ompTP {
public:
    class _checkInCodelets4880 : public darts::Codelet {
    public:
        TP7* myTP;
        TP7* inputsTPParent;
        _checkInCodelets4880()
            : darts::Codelet()
        {
        }
        _checkInCodelets4880(uint32_t dep, uint32_t res, TP7* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets4901 : public darts::Codelet {
    public:
        TP7* myTP;
        TP7* inputsTPParent;
        _checkInCodelets4901()
            : darts::Codelet()
        {
        }
        _checkInCodelets4901(uint32_t dep, uint32_t res, TP7* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets4901 : public darts::Codelet {
    public:
        TP7* inputsTPParent;
        _barrierCodelets4901()
            : darts::Codelet()
        {
        }
        _barrierCodelets4901(uint32_t dep, uint32_t res, TP7* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP7** ptrToThisFunctionTP;
    TP7* inputsTPParent;
    TP7* controlTPParent;
    darts::Codelet** nextCodeletsjacld;
    darts::Codelet** nextSyncCodeletsjacld;
    int* k_darts7 /*VARIABLE*/;
    double* c1345_darts7 /*VARIABLE*/;
    double* c34_darts7 /*VARIABLE*/;
    int* i_darts7 /*VARIABLE*/;
    int* j_darts7 /*VARIABLE*/;
    double* r43_darts7 /*VARIABLE*/;
    double* tmp1_darts7 /*VARIABLE*/;
    double* tmp2_darts7 /*VARIABLE*/;
    double* tmp3_darts7 /*VARIABLE*/;
    TP4901** TP4901Ptr;
    size_t* TP4901_alreadyLaunched;
    int numTPsSet4901;
    int numTPsReady4901;
    size_t TPsToUse4901;
    size_t codeletsPerTP4901;
    size_t totalCodelets4901;
    _checkInCodelets4880* checkInCodelets4880;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets4880* firstCodelet;
#endif
    _checkInCodelets4901* checkInCodelets4901;
    _barrierCodelets4901* barrierCodelets4901;
    TP7(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
        darts::Codelet* in_mainSyncCodelet, TP7** in_ptrToThisFunctionTP, int in_k);
    ~TP7();
    void setNewInputs(int in_k, size_t codeletID);
};
/*TP4901: OMPForDirective*/
class TP4901 : public ompTP {
public:
    class _barrierCodelets4901 : public darts::Codelet {
    public:
        TP4901* inputsTPParent;
        _barrierCodelets4901()
            : darts::Codelet()
        {
        }
        _barrierCodelets4901(uint32_t dep, uint32_t res, TP4901* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations4901(int* endRange, uint32_t codeletID);
    class _checkInCodelets4902 : public darts::Codelet {
    public:
        TP4901* myTP;
        TP4901* inputsTPParent;
        int endRange;
        _checkInCodelets4902()
            : darts::Codelet()
        {
        }
        _checkInCodelets4902(uint32_t dep, uint32_t res, TP4901* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP7* TPParent;
    TP4901* controlTPParent;
    TP4901* inputsTPParent;
    double** c1345_darts4901 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** c34_darts4901 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts4901 /*OMP_PRIVATE - INPUT*/;
    int* j_darts4901 /*OMP_PRIVATE - INPUT*/;
    int** k_darts4901 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** r43_darts4901 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp1_darts4901 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp2_darts4901 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp3_darts4901 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration4901;
    int lastIteration4901;
    int range4901;
    int rangePerCodelet4901;
    int minIteration4901;
    int remainderRange4901;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets4901* barrierCodelets4901;
    _checkInCodelets4902* checkInCodelets4902;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets4902* firstCodelet;
#endif
    TP4901(int in_numThreads, int in_mainCodeletID, TP7* in_TPParent, int in_initIteration,
        int in_lastIteration, TP4901** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP4901();
};
/*TP8: jacu*/
class TP8 : public ompTP {
public:
    class _checkInCodelets7399 : public darts::Codelet {
    public:
        TP8* myTP;
        TP8* inputsTPParent;
        _checkInCodelets7399()
            : darts::Codelet()
        {
        }
        _checkInCodelets7399(uint32_t dep, uint32_t res, TP8* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets7420 : public darts::Codelet {
    public:
        TP8* myTP;
        TP8* inputsTPParent;
        _checkInCodelets7420()
            : darts::Codelet()
        {
        }
        _checkInCodelets7420(uint32_t dep, uint32_t res, TP8* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets7420 : public darts::Codelet {
    public:
        TP8* inputsTPParent;
        _barrierCodelets7420()
            : darts::Codelet()
        {
        }
        _barrierCodelets7420(uint32_t dep, uint32_t res, TP8* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP8** ptrToThisFunctionTP;
    TP8* inputsTPParent;
    TP8* controlTPParent;
    darts::Codelet** nextCodeletsjacu;
    darts::Codelet** nextSyncCodeletsjacu;
    int* k_darts8 /*VARIABLE*/;
    double* c1345_darts8 /*VARIABLE*/;
    double* c34_darts8 /*VARIABLE*/;
    int* i_darts8 /*VARIABLE*/;
    int* j_darts8 /*VARIABLE*/;
    double* r43_darts8 /*VARIABLE*/;
    double* tmp1_darts8 /*VARIABLE*/;
    double* tmp2_darts8 /*VARIABLE*/;
    double* tmp3_darts8 /*VARIABLE*/;
    TP7420** TP7420Ptr;
    size_t* TP7420_alreadyLaunched;
    int numTPsSet7420;
    int numTPsReady7420;
    size_t TPsToUse7420;
    size_t codeletsPerTP7420;
    size_t totalCodelets7420;
    _checkInCodelets7399* checkInCodelets7399;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets7399* firstCodelet;
#endif
    _checkInCodelets7420* checkInCodelets7420;
    _barrierCodelets7420* barrierCodelets7420;
    TP8(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
        darts::Codelet* in_mainSyncCodelet, TP8** in_ptrToThisFunctionTP, int in_k);
    ~TP8();
    void setNewInputs(int in_k, size_t codeletID);
};
/*TP7420: OMPForDirective*/
class TP7420 : public ompTP {
public:
    class _barrierCodelets7420 : public darts::Codelet {
    public:
        TP7420* inputsTPParent;
        _barrierCodelets7420()
            : darts::Codelet()
        {
        }
        _barrierCodelets7420(uint32_t dep, uint32_t res, TP7420* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations7420(int* endRange, uint32_t codeletID);
    class _checkInCodelets7421 : public darts::Codelet {
    public:
        TP7420* myTP;
        TP7420* inputsTPParent;
        int endRange;
        _checkInCodelets7421()
            : darts::Codelet()
        {
        }
        _checkInCodelets7421(uint32_t dep, uint32_t res, TP7420* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP8* TPParent;
    TP7420* controlTPParent;
    TP7420* inputsTPParent;
    double** c1345_darts7420 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** c34_darts7420 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts7420 /*OMP_PRIVATE - INPUT*/;
    int* j_darts7420 /*OMP_PRIVATE - INPUT*/;
    int** k_darts7420 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** r43_darts7420 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp1_darts7420 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp2_darts7420 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp3_darts7420 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration7420;
    int lastIteration7420;
    int range7420;
    int rangePerCodelet7420;
    int minIteration7420;
    int remainderRange7420;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets7420* barrierCodelets7420;
    _checkInCodelets7421* checkInCodelets7421;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets7421* firstCodelet;
#endif
    TP7420(int in_numThreads, int in_mainCodeletID, TP8* in_TPParent, int in_initIteration,
        int in_lastIteration, TP7420** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP7420();
};
/*TP9871: OMPParallelDirective*/
class TP9871 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets9871 : public darts::Codelet {
    public:
        TP9871* inputsTPParent;
        _barrierCodelets9871()
            : darts::Codelet()
        {
        }
        _barrierCodelets9871(uint32_t dep, uint32_t res, TP9871* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets9873 : public darts::Codelet {
    public:
        TP9871* myTP;
        TP9871* inputsTPParent;
        _checkInCodelets9873()
            : darts::Codelet()
        {
        }
        _checkInCodelets9873(uint32_t dep, uint32_t res, TP9871* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets9880 : public darts::Codelet {
    public:
        TP9871* myTP;
        TP9871* inputsTPParent;
        _checkInCodelets9880()
            : darts::Codelet()
        {
        }
        _checkInCodelets9880(uint32_t dep, uint32_t res, TP9871* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets9880 : public darts::Codelet {
    public:
        TP9871* inputsTPParent;
        _barrierCodelets9880()
            : darts::Codelet()
        {
        }
        _barrierCodelets9880(uint32_t dep, uint32_t res, TP9871* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets9889 : public darts::Codelet {
    public:
        TP9871* myTP;
        TP9871* inputsTPParent;
        _checkInCodelets9889()
            : darts::Codelet()
        {
        }
        _checkInCodelets9889(uint32_t dep, uint32_t res, TP9871* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets9986 : public darts::Codelet {
    public:
        TP9871* myTP;
        TP9871* inputsTPParent;
        _checkInCodelets9986()
            : darts::Codelet()
        {
        }
        _checkInCodelets9986(uint32_t dep, uint32_t res, TP9871* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets9996 : public darts::Codelet {
    public:
        TP9871* inputsTPParent;
        _barrierCodelets9996()
            : darts::Codelet()
        {
        }
        _barrierCodelets9996(uint32_t dep, uint32_t res, TP9871* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets9997 : public darts::Codelet {
    public:
        TP9871* myTP;
        TP9871* inputsTPParent;
        _checkInCodelets9997()
            : darts::Codelet()
        {
        }
        _checkInCodelets9997(uint32_t dep, uint32_t res, TP9871* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets9997 : public darts::Codelet {
    public:
        TP9871* inputsTPParent;
        _barrierCodelets9997()
            : darts::Codelet()
        {
        }
        _barrierCodelets9997(uint32_t dep, uint32_t res, TP9871* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP9871* TPParent;
    TP9871* controlTPParent;
    TP9871* inputsTPParent;
    int* iend_darts9871; /*OMP_SHARED - INPUT*/
    int* ist_darts9871; /*OMP_SHARED - INPUT*/
    int* jend_darts9871; /*OMP_SHARED - INPUT*/
    int* jst_darts9871; /*OMP_SHARED - INPUT*/
    int* nx0_darts9871; /*OMP_SHARED - INPUT*/
    int* ny0_darts9871; /*OMP_SHARED - INPUT*/
    int* nz0_darts9871; /*OMP_SHARED - INPUT*/
    double** sum_darts9871; /*OMP_SHARED - INPUT*/
    int* i_darts9871 /*VARIABLE*/;
    int* j_darts9871 /*VARIABLE*/;
    int* k_darts9871 /*VARIABLE*/;
    int* m_darts9871 /*VARIABLE*/;
    double* sum0_darts9871 /*VARIABLE*/;
    double* sum1_darts9871 /*VARIABLE*/;
    double* sum2_darts9871 /*VARIABLE*/;
    double* sum3_darts9871 /*VARIABLE*/;
    double* sum4_darts9871 /*VARIABLE*/;
    int m_darts9997;
    int* nx0_darts9997; /*OMP_SHARED - INPUT*/
    int* ny0_darts9997; /*OMP_SHARED - INPUT*/
    int* nz0_darts9997; /*OMP_SHARED - INPUT*/
    double** sum_darts9997; /*OMP_SHARED - INPUT*/
    int m_darts9880;
    double** sum_darts9880; /*OMP_SHARED - INPUT*/
    size_t TP9880_alreadyLaunched;
    TP9889** TP9889Ptr;
    size_t* TP9889_alreadyLaunched;
    int numTPsSet9889;
    int numTPsReady9889;
    size_t TPsToUse9889;
    size_t codeletsPerTP9889;
    size_t totalCodelets9889;
    size_t TP9997_alreadyLaunched;
    _barrierCodelets9871* barrierCodelets9871;
    _checkInCodelets9873* checkInCodelets9873;
    _checkInCodelets9880* checkInCodelets9880;
    _barrierCodelets9880* barrierCodelets9880;
    _checkInCodelets9889* checkInCodelets9889;
    _checkInCodelets9986* checkInCodelets9986;
    _barrierCodelets9996* barrierCodelets9996;
    _checkInCodelets9997* checkInCodelets9997;
    _barrierCodelets9997* barrierCodelets9997;
    TP9871(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet, int* in_iend,
        int* in_ist, int* in_jend, int* in_jst, int* in_nx0, int* in_ny0, int* in_nz0,
        double** in_sum);
    ~TP9871();
};
/*TP9889: OMPForDirective*/
class TP9889 : public ompTP {
public:
    bool requestNewRangeIterations9889(int* endRange, uint32_t codeletID);
    class _checkInCodelets9890 : public darts::Codelet {
    public:
        TP9889* myTP;
        TP9889* inputsTPParent;
        int endRange;
        _checkInCodelets9890()
            : darts::Codelet()
        {
        }
        _checkInCodelets9890(uint32_t dep, uint32_t res, TP9889* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP9871* TPParent;
    TP9889* controlTPParent;
    TP9889* inputsTPParent;
    int* i_darts9889 /*OMP_PRIVATE - INPUT*/;
    int* iend_darts9889; /*OMP_SHARED - INPUT*/
    int* ist_darts9889; /*OMP_SHARED - INPUT*/
    int* j_darts9889 /*OMP_PRIVATE - INPUT*/;
    int* jend_darts9889; /*OMP_SHARED - INPUT*/
    int* jst_darts9889; /*OMP_SHARED - INPUT*/
    int* k_darts9889 /*OMP_PRIVATE - INPUT*/;
    int* nz0_darts9889; /*OMP_SHARED - INPUT*/
    double** sum0_darts9889 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** sum1_darts9889 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** sum2_darts9889 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** sum3_darts9889 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** sum4_darts9889 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration9889;
    int lastIteration9889;
    int range9889;
    int rangePerCodelet9889;
    int minIteration9889;
    int remainderRange9889;
    size_t readyCodelets;
    int baseNumThreads;
    int* signalNextReady;
    _checkInCodelets9890* checkInCodelets9890;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets9890* firstCodelet;
#endif
    TP9889(int in_numThreads, int in_mainCodeletID, TP9871* in_TPParent, int in_initIteration,
        int in_lastIteration, int* in_iend, int* in_ist, int* in_jend, int* in_jst, int* in_nz0,
        TP9889** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP9889();
};
/*TP10788: OMPParallelDirective*/
class TP10788 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets10788 : public darts::Codelet {
    public:
        TP10788* inputsTPParent;
        _barrierCodelets10788()
            : darts::Codelet()
        {
        }
        _barrierCodelets10788(uint32_t dep, uint32_t res, TP10788* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets10803 : public darts::Codelet {
    public:
        TP10788* myTP;
        TP10788* inputsTPParent;
        _checkInCodelets10803()
            : darts::Codelet()
        {
        }
        _checkInCodelets10803(uint32_t dep, uint32_t res, TP10788* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets10803 : public darts::Codelet {
    public:
        TP10788* inputsTPParent;
        _barrierCodelets10803()
            : darts::Codelet()
        {
        }
        _barrierCodelets10803(uint32_t dep, uint32_t res, TP10788* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets10860 : public darts::Codelet {
    public:
        TP10788* myTP;
        TP10788* inputsTPParent;
        _checkInCodelets10860()
            : darts::Codelet()
        {
        }
        _checkInCodelets10860(uint32_t dep, uint32_t res, TP10788* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets10863 : public darts::Codelet {
    public:
        TP10788* myTP;
        TP10788* inputsTPParent;
        _checkInCodelets10863()
            : darts::Codelet()
        {
        }
        _checkInCodelets10863(uint32_t dep, uint32_t res, TP10788* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets10863 : public darts::Codelet {
    public:
        TP10788* inputsTPParent;
        _barrierCodelets10863()
            : darts::Codelet()
        {
        }
        _barrierCodelets10863(uint32_t dep, uint32_t res, TP10788* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets11012 : public darts::Codelet {
    public:
        TP10788* myTP;
        TP10788* inputsTPParent;
        _checkInCodelets11012()
            : darts::Codelet()
        {
        }
        _checkInCodelets11012(uint32_t dep, uint32_t res, TP10788* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets11012 : public darts::Codelet {
    public:
        TP10788* inputsTPParent;
        _barrierCodelets11012()
            : darts::Codelet()
        {
        }
        _barrierCodelets11012(uint32_t dep, uint32_t res, TP10788* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets11649 : public darts::Codelet {
    public:
        TP10788* myTP;
        TP10788* inputsTPParent;
        _checkInCodelets11649()
            : darts::Codelet()
        {
        }
        _checkInCodelets11649(uint32_t dep, uint32_t res, TP10788* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets11652 : public darts::Codelet {
    public:
        TP10788* myTP;
        TP10788* inputsTPParent;
        _checkInCodelets11652()
            : darts::Codelet()
        {
        }
        _checkInCodelets11652(uint32_t dep, uint32_t res, TP10788* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets11652 : public darts::Codelet {
    public:
        TP10788* inputsTPParent;
        _barrierCodelets11652()
            : darts::Codelet()
        {
        }
        _barrierCodelets11652(uint32_t dep, uint32_t res, TP10788* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets11801 : public darts::Codelet {
    public:
        TP10788* myTP;
        TP10788* inputsTPParent;
        _checkInCodelets11801()
            : darts::Codelet()
        {
        }
        _checkInCodelets11801(uint32_t dep, uint32_t res, TP10788* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets11801 : public darts::Codelet {
    public:
        TP10788* inputsTPParent;
        _barrierCodelets11801()
            : darts::Codelet()
        {
        }
        _barrierCodelets11801(uint32_t dep, uint32_t res, TP10788* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets12438 : public darts::Codelet {
    public:
        TP10788* myTP;
        TP10788* inputsTPParent;
        _checkInCodelets12438()
            : darts::Codelet()
        {
        }
        _checkInCodelets12438(uint32_t dep, uint32_t res, TP10788* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets12438 : public darts::Codelet {
    public:
        TP10788* inputsTPParent;
        _barrierCodelets12438()
            : darts::Codelet()
        {
        }
        _barrierCodelets12438(uint32_t dep, uint32_t res, TP10788* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP10788* TPParent;
    TP10788* controlTPParent;
    TP10788* inputsTPParent;
    int* L1_darts10788 /*VARIABLE*/;
    int* L2_darts10788 /*VARIABLE*/;
    int* i_darts10788 /*VARIABLE*/;
    int* iend1_darts10788 /*VARIABLE*/;
    int* ist1_darts10788 /*VARIABLE*/;
    int* j_darts10788 /*VARIABLE*/;
    int* jend1_darts10788 /*VARIABLE*/;
    int* jst1_darts10788 /*VARIABLE*/;
    int* k_darts10788 /*VARIABLE*/;
    int* m_darts10788 /*VARIABLE*/;
    double* q_darts10788 /*VARIABLE*/;
    double* tmp_darts10788 /*VARIABLE*/;
    double* u21_darts10788 /*VARIABLE*/;
    double* u21i_darts10788 /*VARIABLE*/;
    double* u21im1_darts10788 /*VARIABLE*/;
    double* u21j_darts10788 /*VARIABLE*/;
    double* u21jm1_darts10788 /*VARIABLE*/;
    double* u21k_darts10788 /*VARIABLE*/;
    double* u21km1_darts10788 /*VARIABLE*/;
    double* u31_darts10788 /*VARIABLE*/;
    double* u31i_darts10788 /*VARIABLE*/;
    double* u31im1_darts10788 /*VARIABLE*/;
    double* u31j_darts10788 /*VARIABLE*/;
    double* u31jm1_darts10788 /*VARIABLE*/;
    double* u31k_darts10788 /*VARIABLE*/;
    double* u31km1_darts10788 /*VARIABLE*/;
    double* u41_darts10788 /*VARIABLE*/;
    double* u41i_darts10788 /*VARIABLE*/;
    double* u41im1_darts10788 /*VARIABLE*/;
    double* u41j_darts10788 /*VARIABLE*/;
    double* u41jm1_darts10788 /*VARIABLE*/;
    double* u41k_darts10788 /*VARIABLE*/;
    double* u41km1_darts10788 /*VARIABLE*/;
    double* u51i_darts10788 /*VARIABLE*/;
    double* u51im1_darts10788 /*VARIABLE*/;
    double* u51j_darts10788 /*VARIABLE*/;
    double* u51jm1_darts10788 /*VARIABLE*/;
    double* u51k_darts10788 /*VARIABLE*/;
    double* u51km1_darts10788 /*VARIABLE*/;
    TP10803** TP10803Ptr;
    size_t* TP10803_alreadyLaunched;
    int numTPsSet10803;
    int numTPsReady10803;
    size_t TPsToUse10803;
    size_t codeletsPerTP10803;
    size_t totalCodelets10803;
    TP10863** TP10863Ptr;
    size_t* TP10863_alreadyLaunched;
    int numTPsSet10863;
    int numTPsReady10863;
    size_t TPsToUse10863;
    size_t codeletsPerTP10863;
    size_t totalCodelets10863;
    TP11012** TP11012Ptr;
    size_t* TP11012_alreadyLaunched;
    int numTPsSet11012;
    int numTPsReady11012;
    size_t TPsToUse11012;
    size_t codeletsPerTP11012;
    size_t totalCodelets11012;
    TP11652** TP11652Ptr;
    size_t* TP11652_alreadyLaunched;
    int numTPsSet11652;
    int numTPsReady11652;
    size_t TPsToUse11652;
    size_t codeletsPerTP11652;
    size_t totalCodelets11652;
    TP11801** TP11801Ptr;
    size_t* TP11801_alreadyLaunched;
    int numTPsSet11801;
    int numTPsReady11801;
    size_t TPsToUse11801;
    size_t codeletsPerTP11801;
    size_t totalCodelets11801;
    TP12438** TP12438Ptr;
    size_t* TP12438_alreadyLaunched;
    int numTPsSet12438;
    int numTPsReady12438;
    size_t TPsToUse12438;
    size_t codeletsPerTP12438;
    size_t totalCodelets12438;
    _barrierCodelets10788* barrierCodelets10788;
    _checkInCodelets10803* checkInCodelets10803;
    _barrierCodelets10803* barrierCodelets10803;
    _checkInCodelets10860* checkInCodelets10860;
    _checkInCodelets10863* checkInCodelets10863;
    _barrierCodelets10863* barrierCodelets10863;
    _checkInCodelets11012* checkInCodelets11012;
    _barrierCodelets11012* barrierCodelets11012;
    _checkInCodelets11649* checkInCodelets11649;
    _checkInCodelets11652* checkInCodelets11652;
    _barrierCodelets11652* barrierCodelets11652;
    _checkInCodelets11801* checkInCodelets11801;
    _barrierCodelets11801* barrierCodelets11801;
    _checkInCodelets12438* checkInCodelets12438;
    _barrierCodelets12438* barrierCodelets12438;
    TP10788(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
    ~TP10788();
};
/*TP10803: OMPForDirective*/
class TP10803 : public ompTP {
public:
    class _barrierCodelets10803 : public darts::Codelet {
    public:
        TP10803* inputsTPParent;
        _barrierCodelets10803()
            : darts::Codelet()
        {
        }
        _barrierCodelets10803(uint32_t dep, uint32_t res, TP10803* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations10803(int* endRange, uint32_t codeletID);
    class _checkInCodelets10804 : public darts::Codelet {
    public:
        TP10803* myTP;
        TP10803* inputsTPParent;
        int endRange;
        _checkInCodelets10804()
            : darts::Codelet()
        {
        }
        _checkInCodelets10804(uint32_t dep, uint32_t res, TP10803* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10788* TPParent;
    TP10803* controlTPParent;
    TP10803* inputsTPParent;
    int* i_darts10803 /*OMP_PRIVATE - INPUT*/;
    int* j_darts10803 /*OMP_PRIVATE - INPUT*/;
    int* k_darts10803 /*OMP_PRIVATE - INPUT*/;
    int* m_darts10803 /*OMP_PRIVATE - INPUT*/;
    int initIteration10803;
    int lastIteration10803;
    int range10803;
    int rangePerCodelet10803;
    int minIteration10803;
    int remainderRange10803;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets10803* barrierCodelets10803;
    _checkInCodelets10804* checkInCodelets10804;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets10804* firstCodelet;
#endif
    TP10803(int in_numThreads, int in_mainCodeletID, TP10788* in_TPParent, int in_initIteration,
        int in_lastIteration, TP10803** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP10803();
};
/*TP10863: OMPForDirective*/
class TP10863 : public ompTP {
public:
    class _barrierCodelets10863 : public darts::Codelet {
    public:
        TP10863* inputsTPParent;
        _barrierCodelets10863()
            : darts::Codelet()
        {
        }
        _barrierCodelets10863(uint32_t dep, uint32_t res, TP10863* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations10863(int* endRange, uint32_t codeletID);
    class _checkInCodelets10864 : public darts::Codelet {
    public:
        TP10863* myTP;
        TP10863* inputsTPParent;
        int endRange;
        _checkInCodelets10864()
            : darts::Codelet()
        {
        }
        _checkInCodelets10864(uint32_t dep, uint32_t res, TP10863* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10788* TPParent;
    TP10863* controlTPParent;
    TP10863* inputsTPParent;
    int** L1_darts10863 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** L2_darts10863 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts10863 /*OMP_PRIVATE - INPUT*/;
    int* j_darts10863 /*OMP_PRIVATE - INPUT*/;
    int* k_darts10863 /*OMP_PRIVATE - INPUT*/;
    double** q_darts10863 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21_darts10863 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration10863;
    int lastIteration10863;
    int range10863;
    int rangePerCodelet10863;
    int minIteration10863;
    int remainderRange10863;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets10863* barrierCodelets10863;
    _checkInCodelets10864* checkInCodelets10864;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets10864* firstCodelet;
#endif
    TP10863(int in_numThreads, int in_mainCodeletID, TP10788* in_TPParent, int in_initIteration,
        int in_lastIteration, TP10863** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP10863();
};
/*TP11012: OMPForDirective*/
class TP11012 : public ompTP {
public:
    class _barrierCodelets11012 : public darts::Codelet {
    public:
        TP11012* inputsTPParent;
        _barrierCodelets11012()
            : darts::Codelet()
        {
        }
        _barrierCodelets11012(uint32_t dep, uint32_t res, TP11012* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations11012(int* endRange, uint32_t codeletID);
    class _checkInCodelets11013 : public darts::Codelet {
    public:
        TP11012* myTP;
        TP11012* inputsTPParent;
        int endRange;
        _checkInCodelets11013()
            : darts::Codelet()
        {
        }
        _checkInCodelets11013(uint32_t dep, uint32_t res, TP11012* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10788* TPParent;
    TP11012* controlTPParent;
    TP11012* inputsTPParent;
    int** L2_darts11012 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts11012 /*OMP_PRIVATE - INPUT*/;
    int** iend1_darts11012 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** ist1_darts11012 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts11012 /*OMP_PRIVATE - INPUT*/;
    int* k_darts11012 /*OMP_PRIVATE - INPUT*/;
    int* m_darts11012 /*OMP_PRIVATE - INPUT*/;
    double** tmp_darts11012 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21i_darts11012 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21im1_darts11012 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31i_darts11012 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31im1_darts11012 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41i_darts11012 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41im1_darts11012 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51i_darts11012 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51im1_darts11012 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration11012;
    int lastIteration11012;
    int range11012;
    int rangePerCodelet11012;
    int minIteration11012;
    int remainderRange11012;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets11012* barrierCodelets11012;
    _checkInCodelets11013* checkInCodelets11013;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets11013* firstCodelet;
#endif
    TP11012(int in_numThreads, int in_mainCodeletID, TP10788* in_TPParent, int in_initIteration,
        int in_lastIteration, TP11012** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP11012();
};
/*TP11652: OMPForDirective*/
class TP11652 : public ompTP {
public:
    class _barrierCodelets11652 : public darts::Codelet {
    public:
        TP11652* inputsTPParent;
        _barrierCodelets11652()
            : darts::Codelet()
        {
        }
        _barrierCodelets11652(uint32_t dep, uint32_t res, TP11652* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations11652(int* endRange, uint32_t codeletID);
    class _checkInCodelets11653 : public darts::Codelet {
    public:
        TP11652* myTP;
        TP11652* inputsTPParent;
        int endRange;
        _checkInCodelets11653()
            : darts::Codelet()
        {
        }
        _checkInCodelets11653(uint32_t dep, uint32_t res, TP11652* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10788* TPParent;
    TP11652* controlTPParent;
    TP11652* inputsTPParent;
    int** L1_darts11652 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** L2_darts11652 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts11652 /*OMP_PRIVATE - INPUT*/;
    int* j_darts11652 /*OMP_PRIVATE - INPUT*/;
    int* k_darts11652 /*OMP_PRIVATE - INPUT*/;
    double** q_darts11652 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31_darts11652 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration11652;
    int lastIteration11652;
    int range11652;
    int rangePerCodelet11652;
    int minIteration11652;
    int remainderRange11652;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets11652* barrierCodelets11652;
    _checkInCodelets11653* checkInCodelets11653;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets11653* firstCodelet;
#endif
    TP11652(int in_numThreads, int in_mainCodeletID, TP10788* in_TPParent, int in_initIteration,
        int in_lastIteration, TP11652** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP11652();
};
/*TP11801: OMPForDirective*/
class TP11801 : public ompTP {
public:
    class _barrierCodelets11801 : public darts::Codelet {
    public:
        TP11801* inputsTPParent;
        _barrierCodelets11801()
            : darts::Codelet()
        {
        }
        _barrierCodelets11801(uint32_t dep, uint32_t res, TP11801* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations11801(int* endRange, uint32_t codeletID);
    class _checkInCodelets11802 : public darts::Codelet {
    public:
        TP11801* myTP;
        TP11801* inputsTPParent;
        int endRange;
        _checkInCodelets11802()
            : darts::Codelet()
        {
        }
        _checkInCodelets11802(uint32_t dep, uint32_t res, TP11801* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10788* TPParent;
    TP11801* controlTPParent;
    TP11801* inputsTPParent;
    int** L2_darts11801 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts11801 /*OMP_PRIVATE - INPUT*/;
    int* j_darts11801 /*OMP_PRIVATE - INPUT*/;
    int** jend1_darts11801 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** jst1_darts11801 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts11801 /*OMP_PRIVATE - INPUT*/;
    int* m_darts11801 /*OMP_PRIVATE - INPUT*/;
    double** tmp_darts11801 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21j_darts11801 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21jm1_darts11801 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31j_darts11801 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31jm1_darts11801 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41j_darts11801 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41jm1_darts11801 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51j_darts11801 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51jm1_darts11801 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration11801;
    int lastIteration11801;
    int range11801;
    int rangePerCodelet11801;
    int minIteration11801;
    int remainderRange11801;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets11801* barrierCodelets11801;
    _checkInCodelets11802* checkInCodelets11802;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets11802* firstCodelet;
#endif
    TP11801(int in_numThreads, int in_mainCodeletID, TP10788* in_TPParent, int in_initIteration,
        int in_lastIteration, TP11801** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP11801();
};
/*TP12438: OMPForDirective*/
class TP12438 : public ompTP {
public:
    class _barrierCodelets12438 : public darts::Codelet {
    public:
        TP12438* inputsTPParent;
        _barrierCodelets12438()
            : darts::Codelet()
        {
        }
        _barrierCodelets12438(uint32_t dep, uint32_t res, TP12438* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations12438(int* endRange, uint32_t codeletID);
    class _checkInCodelets12439 : public darts::Codelet {
    public:
        TP12438* myTP;
        TP12438* inputsTPParent;
        int endRange;
        _checkInCodelets12439()
            : darts::Codelet()
        {
        }
        _checkInCodelets12439(uint32_t dep, uint32_t res, TP12438* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10788* TPParent;
    TP12438* controlTPParent;
    TP12438* inputsTPParent;
    int* i_darts12438 /*OMP_PRIVATE - INPUT*/;
    int* j_darts12438 /*OMP_PRIVATE - INPUT*/;
    int* k_darts12438 /*OMP_PRIVATE - INPUT*/;
    int* m_darts12438 /*OMP_PRIVATE - INPUT*/;
    double** q_darts12438 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp_darts12438 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21k_darts12438 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21km1_darts12438 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31k_darts12438 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31km1_darts12438 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41_darts12438 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41k_darts12438 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41km1_darts12438 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51k_darts12438 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51km1_darts12438 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration12438;
    int lastIteration12438;
    int range12438;
    int rangePerCodelet12438;
    int minIteration12438;
    int remainderRange12438;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets12438* barrierCodelets12438;
    _checkInCodelets12439* checkInCodelets12439;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets12439* firstCodelet;
#endif
    TP12438(int in_numThreads, int in_mainCodeletID, TP10788* in_TPParent, int in_initIteration,
        int in_lastIteration, TP12438** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP12438();
};
/*TP13189: OMPParallelDirective*/
class TP13189 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets13189 : public darts::Codelet {
    public:
        TP13189* inputsTPParent;
        _barrierCodelets13189()
            : darts::Codelet()
        {
        }
        _barrierCodelets13189(uint32_t dep, uint32_t res, TP13189* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13193 : public darts::Codelet {
    public:
        TP13189* myTP;
        TP13189* inputsTPParent;
        _checkInCodelets13193()
            : darts::Codelet()
        {
        }
        _checkInCodelets13193(uint32_t dep, uint32_t res, TP13189* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13193 : public darts::Codelet {
    public:
        TP13189* inputsTPParent;
        _barrierCodelets13193()
            : darts::Codelet()
        {
        }
        _barrierCodelets13193(uint32_t dep, uint32_t res, TP13189* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13244 : public darts::Codelet {
    public:
        TP13189* myTP;
        TP13189* inputsTPParent;
        _checkInCodelets13244()
            : darts::Codelet()
        {
        }
        _checkInCodelets13244(uint32_t dep, uint32_t res, TP13189* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13244 : public darts::Codelet {
    public:
        TP13189* inputsTPParent;
        _barrierCodelets13244()
            : darts::Codelet()
        {
        }
        _barrierCodelets13244(uint32_t dep, uint32_t res, TP13189* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13286 : public darts::Codelet {
    public:
        TP13189* myTP;
        TP13189* inputsTPParent;
        _checkInCodelets13286()
            : darts::Codelet()
        {
        }
        _checkInCodelets13286(uint32_t dep, uint32_t res, TP13189* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13286 : public darts::Codelet {
    public:
        TP13189* inputsTPParent;
        _barrierCodelets13286()
            : darts::Codelet()
        {
        }
        _barrierCodelets13286(uint32_t dep, uint32_t res, TP13189* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13330 : public darts::Codelet {
    public:
        TP13189* myTP;
        TP13189* inputsTPParent;
        _checkInCodelets13330()
            : darts::Codelet()
        {
        }
        _checkInCodelets13330(uint32_t dep, uint32_t res, TP13189* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13330 : public darts::Codelet {
    public:
        TP13189* inputsTPParent;
        _barrierCodelets13330()
            : darts::Codelet()
        {
        }
        _barrierCodelets13330(uint32_t dep, uint32_t res, TP13189* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13372 : public darts::Codelet {
    public:
        TP13189* myTP;
        TP13189* inputsTPParent;
        _checkInCodelets13372()
            : darts::Codelet()
        {
        }
        _checkInCodelets13372(uint32_t dep, uint32_t res, TP13189* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13372 : public darts::Codelet {
    public:
        TP13189* inputsTPParent;
        _barrierCodelets13372()
            : darts::Codelet()
        {
        }
        _barrierCodelets13372(uint32_t dep, uint32_t res, TP13189* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP13189* TPParent;
    TP13189* controlTPParent;
    TP13189* inputsTPParent;
    int* i_darts13189 /*VARIABLE*/;
    int* iglob_darts13189 /*VARIABLE*/;
    int* j_darts13189 /*VARIABLE*/;
    int* jglob_darts13189 /*VARIABLE*/;
    int* k_darts13189 /*VARIABLE*/;
    TP13193** TP13193Ptr;
    size_t* TP13193_alreadyLaunched;
    int numTPsSet13193;
    int numTPsReady13193;
    size_t TPsToUse13193;
    size_t codeletsPerTP13193;
    size_t totalCodelets13193;
    TP13244** TP13244Ptr;
    size_t* TP13244_alreadyLaunched;
    int numTPsSet13244;
    int numTPsReady13244;
    size_t TPsToUse13244;
    size_t codeletsPerTP13244;
    size_t totalCodelets13244;
    TP13286** TP13286Ptr;
    size_t* TP13286_alreadyLaunched;
    int numTPsSet13286;
    int numTPsReady13286;
    size_t TPsToUse13286;
    size_t codeletsPerTP13286;
    size_t totalCodelets13286;
    TP13330** TP13330Ptr;
    size_t* TP13330_alreadyLaunched;
    int numTPsSet13330;
    int numTPsReady13330;
    size_t TPsToUse13330;
    size_t codeletsPerTP13330;
    size_t totalCodelets13330;
    TP13372** TP13372Ptr;
    size_t* TP13372_alreadyLaunched;
    int numTPsSet13372;
    int numTPsReady13372;
    size_t TPsToUse13372;
    size_t codeletsPerTP13372;
    size_t totalCodelets13372;
    _barrierCodelets13189* barrierCodelets13189;
    _checkInCodelets13193* checkInCodelets13193;
    _barrierCodelets13193* barrierCodelets13193;
    _checkInCodelets13244* checkInCodelets13244;
    _barrierCodelets13244* barrierCodelets13244;
    _checkInCodelets13286* checkInCodelets13286;
    _barrierCodelets13286* barrierCodelets13286;
    _checkInCodelets13330* checkInCodelets13330;
    _barrierCodelets13330* barrierCodelets13330;
    _checkInCodelets13372* checkInCodelets13372;
    _barrierCodelets13372* barrierCodelets13372;
    TP13189(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
    ~TP13189();
};
/*TP13193: OMPForDirective*/
class TP13193 : public ompTP {
public:
    class _barrierCodelets13193 : public darts::Codelet {
    public:
        TP13193* inputsTPParent;
        _barrierCodelets13193()
            : darts::Codelet()
        {
        }
        _barrierCodelets13193(uint32_t dep, uint32_t res, TP13193* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13193(int* endRange, uint32_t codeletID);
    class _checkInCodelets13194 : public darts::Codelet {
    public:
        TP13193* myTP;
        TP13193* inputsTPParent;
        int endRange;
        _checkInCodelets13194()
            : darts::Codelet()
        {
        }
        _checkInCodelets13194(uint32_t dep, uint32_t res, TP13193* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13189* TPParent;
    TP13193* controlTPParent;
    TP13193* inputsTPParent;
    int* i_darts13193 /*OMP_PRIVATE - INPUT*/;
    int** iglob_darts13193 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts13193 /*OMP_PRIVATE - INPUT*/;
    int** jglob_darts13193 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration13193;
    int lastIteration13193;
    int range13193;
    int rangePerCodelet13193;
    int minIteration13193;
    int remainderRange13193;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13193* barrierCodelets13193;
    _checkInCodelets13194* checkInCodelets13194;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13194* firstCodelet;
#endif
    TP13193(int in_numThreads, int in_mainCodeletID, TP13189* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13193** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13193();
};
/*TP13244: OMPForDirective*/
class TP13244 : public ompTP {
public:
    class _barrierCodelets13244 : public darts::Codelet {
    public:
        TP13244* inputsTPParent;
        _barrierCodelets13244()
            : darts::Codelet()
        {
        }
        _barrierCodelets13244(uint32_t dep, uint32_t res, TP13244* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13244(int* endRange, uint32_t codeletID);
    class _checkInCodelets13245 : public darts::Codelet {
    public:
        TP13244* myTP;
        TP13244* inputsTPParent;
        int endRange;
        _checkInCodelets13245()
            : darts::Codelet()
        {
        }
        _checkInCodelets13245(uint32_t dep, uint32_t res, TP13244* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13189* TPParent;
    TP13244* controlTPParent;
    TP13244* inputsTPParent;
    int* i_darts13244 /*OMP_PRIVATE - INPUT*/;
    int** iglob_darts13244 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts13244 /*OMP_PRIVATE - INPUT*/;
    int initIteration13244;
    int lastIteration13244;
    int range13244;
    int rangePerCodelet13244;
    int minIteration13244;
    int remainderRange13244;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13244* barrierCodelets13244;
    _checkInCodelets13245* checkInCodelets13245;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13245* firstCodelet;
#endif
    TP13244(int in_numThreads, int in_mainCodeletID, TP13189* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13244** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13244();
};
/*TP13286: OMPForDirective*/
class TP13286 : public ompTP {
public:
    class _barrierCodelets13286 : public darts::Codelet {
    public:
        TP13286* inputsTPParent;
        _barrierCodelets13286()
            : darts::Codelet()
        {
        }
        _barrierCodelets13286(uint32_t dep, uint32_t res, TP13286* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13286(int* endRange, uint32_t codeletID);
    class _checkInCodelets13287 : public darts::Codelet {
    public:
        TP13286* myTP;
        TP13286* inputsTPParent;
        int endRange;
        _checkInCodelets13287()
            : darts::Codelet()
        {
        }
        _checkInCodelets13287(uint32_t dep, uint32_t res, TP13286* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13189* TPParent;
    TP13286* controlTPParent;
    TP13286* inputsTPParent;
    int* i_darts13286 /*OMP_PRIVATE - INPUT*/;
    int** iglob_darts13286 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts13286 /*OMP_PRIVATE - INPUT*/;
    int initIteration13286;
    int lastIteration13286;
    int range13286;
    int rangePerCodelet13286;
    int minIteration13286;
    int remainderRange13286;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13286* barrierCodelets13286;
    _checkInCodelets13287* checkInCodelets13287;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13287* firstCodelet;
#endif
    TP13286(int in_numThreads, int in_mainCodeletID, TP13189* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13286** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13286();
};
/*TP13330: OMPForDirective*/
class TP13330 : public ompTP {
public:
    class _barrierCodelets13330 : public darts::Codelet {
    public:
        TP13330* inputsTPParent;
        _barrierCodelets13330()
            : darts::Codelet()
        {
        }
        _barrierCodelets13330(uint32_t dep, uint32_t res, TP13330* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13330(int* endRange, uint32_t codeletID);
    class _checkInCodelets13331 : public darts::Codelet {
    public:
        TP13330* myTP;
        TP13330* inputsTPParent;
        int endRange;
        _checkInCodelets13331()
            : darts::Codelet()
        {
        }
        _checkInCodelets13331(uint32_t dep, uint32_t res, TP13330* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13189* TPParent;
    TP13330* controlTPParent;
    TP13330* inputsTPParent;
    int* j_darts13330 /*OMP_PRIVATE - INPUT*/;
    int** jglob_darts13330 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts13330 /*OMP_PRIVATE - INPUT*/;
    int initIteration13330;
    int lastIteration13330;
    int range13330;
    int rangePerCodelet13330;
    int minIteration13330;
    int remainderRange13330;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13330* barrierCodelets13330;
    _checkInCodelets13331* checkInCodelets13331;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13331* firstCodelet;
#endif
    TP13330(int in_numThreads, int in_mainCodeletID, TP13189* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13330** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13330();
};
/*TP13372: OMPForDirective*/
class TP13372 : public ompTP {
public:
    class _barrierCodelets13372 : public darts::Codelet {
    public:
        TP13372* inputsTPParent;
        _barrierCodelets13372()
            : darts::Codelet()
        {
        }
        _barrierCodelets13372(uint32_t dep, uint32_t res, TP13372* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13372(int* endRange, uint32_t codeletID);
    class _checkInCodelets13373 : public darts::Codelet {
    public:
        TP13372* myTP;
        TP13372* inputsTPParent;
        int endRange;
        _checkInCodelets13373()
            : darts::Codelet()
        {
        }
        _checkInCodelets13373(uint32_t dep, uint32_t res, TP13372* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13189* TPParent;
    TP13372* controlTPParent;
    TP13372* inputsTPParent;
    int* j_darts13372 /*OMP_PRIVATE - INPUT*/;
    int** jglob_darts13372 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts13372 /*OMP_PRIVATE - INPUT*/;
    int initIteration13372;
    int lastIteration13372;
    int range13372;
    int rangePerCodelet13372;
    int minIteration13372;
    int remainderRange13372;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13372* barrierCodelets13372;
    _checkInCodelets13373* checkInCodelets13373;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13373* firstCodelet;
#endif
    TP13372(int in_numThreads, int in_mainCodeletID, TP13189* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13372** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13372();
};
/*TP13756: OMPParallelDirective*/
class TP13756 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets13756 : public darts::Codelet {
    public:
        TP13756* inputsTPParent;
        _barrierCodelets13756()
            : darts::Codelet()
        {
        }
        _barrierCodelets13756(uint32_t dep, uint32_t res, TP13756* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13762 : public darts::Codelet {
    public:
        TP13756* myTP;
        TP13756* inputsTPParent;
        _checkInCodelets13762()
            : darts::Codelet()
        {
        }
        _checkInCodelets13762(uint32_t dep, uint32_t res, TP13756* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13762 : public darts::Codelet {
    public:
        TP13756* inputsTPParent;
        _barrierCodelets13762()
            : darts::Codelet()
        {
        }
        _barrierCodelets13762(uint32_t dep, uint32_t res, TP13756* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP13756* TPParent;
    TP13756* controlTPParent;
    TP13756* inputsTPParent;
    double* eta_darts13756 /*VARIABLE*/;
    int* i_darts13756 /*VARIABLE*/;
    int* iglob_darts13756 /*VARIABLE*/;
    int* j_darts13756 /*VARIABLE*/;
    int* jglob_darts13756 /*VARIABLE*/;
    int* k_darts13756 /*VARIABLE*/;
    int* m_darts13756 /*VARIABLE*/;
    double* peta_darts13756 /*VARIABLE*/;
    double* pxi_darts13756 /*VARIABLE*/;
    double* pzeta_darts13756 /*VARIABLE*/;
    double* xi_darts13756 /*VARIABLE*/;
    double* zeta_darts13756 /*VARIABLE*/;
    TP13762** TP13762Ptr;
    size_t* TP13762_alreadyLaunched;
    int numTPsSet13762;
    int numTPsReady13762;
    size_t TPsToUse13762;
    size_t codeletsPerTP13762;
    size_t totalCodelets13762;
    _barrierCodelets13756* barrierCodelets13756;
    _checkInCodelets13762* checkInCodelets13762;
    _barrierCodelets13762* barrierCodelets13762;
    TP13756(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
    ~TP13756();
};
/*TP13762: OMPForDirective*/
class TP13762 : public ompTP {
public:
    class _barrierCodelets13762 : public darts::Codelet {
    public:
        TP13762* inputsTPParent;
        _barrierCodelets13762()
            : darts::Codelet()
        {
        }
        _barrierCodelets13762(uint32_t dep, uint32_t res, TP13762* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13762(int* endRange, uint32_t codeletID);
    class _checkInCodelets13763 : public darts::Codelet {
    public:
        TP13762* myTP;
        TP13762* inputsTPParent;
        int endRange;
        _checkInCodelets13763()
            : darts::Codelet()
        {
        }
        _checkInCodelets13763(uint32_t dep, uint32_t res, TP13762* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13756* TPParent;
    TP13762* controlTPParent;
    TP13762* inputsTPParent;
    double** eta_darts13762 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts13762 /*OMP_PRIVATE - INPUT*/;
    int** iglob_darts13762 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts13762 /*OMP_PRIVATE - INPUT*/;
    int** jglob_darts13762 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts13762 /*OMP_PRIVATE - INPUT*/;
    int* m_darts13762 /*OMP_PRIVATE - INPUT*/;
    double** peta_darts13762 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** pxi_darts13762 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** pzeta_darts13762 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** xi_darts13762 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** zeta_darts13762 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration13762;
    int lastIteration13762;
    int range13762;
    int rangePerCodelet13762;
    int minIteration13762;
    int remainderRange13762;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13762* barrierCodelets13762;
    _checkInCodelets13763* checkInCodelets13763;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13763* firstCodelet;
#endif
    TP13762(int in_numThreads, int in_mainCodeletID, TP13756* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13762** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13762();
};
/*TP13888: OMPParallelDirective*/
class TP13888 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets13888 : public darts::Codelet {
    public:
        TP13888* inputsTPParent;
        _barrierCodelets13888()
            : darts::Codelet()
        {
        }
        _barrierCodelets13888(uint32_t dep, uint32_t res, TP13888* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13890 : public darts::Codelet {
    public:
        TP13888* myTP;
        TP13888* inputsTPParent;
        _checkInCodelets13890()
            : darts::Codelet()
        {
        }
        _checkInCodelets13890(uint32_t dep, uint32_t res, TP13888* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13890 : public darts::Codelet {
    public:
        TP13888* inputsTPParent;
        _barrierCodelets13890()
            : darts::Codelet()
        {
        }
        _barrierCodelets13890(uint32_t dep, uint32_t res, TP13888* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP13888* TPParent;
    TP13888* controlTPParent;
    TP13888* inputsTPParent;
    int* i_darts13888 /*OMP_PRIVATE - INPUT*/;
    int* j_darts13888 /*OMP_PRIVATE - INPUT*/;
    int* k_darts13888 /*OMP_PRIVATE - INPUT*/;
    int* m_darts13888 /*OMP_PRIVATE - INPUT*/;
    TP13890** TP13890Ptr;
    size_t* TP13890_alreadyLaunched;
    int numTPsSet13890;
    int numTPsReady13890;
    size_t TPsToUse13890;
    size_t codeletsPerTP13890;
    size_t totalCodelets13890;
    _barrierCodelets13888* barrierCodelets13888;
    _checkInCodelets13890* checkInCodelets13890;
    _barrierCodelets13890* barrierCodelets13890;
    TP13888(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
    ~TP13888();
};
/*TP13890: OMPForDirective*/
class TP13890 : public ompTP {
public:
    class _barrierCodelets13890 : public darts::Codelet {
    public:
        TP13890* inputsTPParent;
        _barrierCodelets13890()
            : darts::Codelet()
        {
        }
        _barrierCodelets13890(uint32_t dep, uint32_t res, TP13890* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13890(int* endRange, uint32_t codeletID);
    class _checkInCodelets13891 : public darts::Codelet {
    public:
        TP13890* myTP;
        TP13890* inputsTPParent;
        int endRange;
        _checkInCodelets13891()
            : darts::Codelet()
        {
        }
        _checkInCodelets13891(uint32_t dep, uint32_t res, TP13890* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13888* TPParent;
    TP13890* controlTPParent;
    TP13890* inputsTPParent;
    int* i_darts13890 /*OMP_PRIVATE - INPUT*/;
    int* j_darts13890 /*OMP_PRIVATE - INPUT*/;
    int* k_darts13890 /*OMP_PRIVATE - INPUT*/;
    int* m_darts13890 /*OMP_PRIVATE - INPUT*/;
    int initIteration13890;
    int lastIteration13890;
    int range13890;
    int rangePerCodelet13890;
    int minIteration13890;
    int remainderRange13890;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13890* barrierCodelets13890;
    _checkInCodelets13891* checkInCodelets13891;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13891* firstCodelet;
#endif
    TP13890(int in_numThreads, int in_mainCodeletID, TP13888* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13890** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13890();
};
/*TP13987: OMPParallelDirective*/
class TP13987 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets13987 : public darts::Codelet {
    public:
        TP13987* inputsTPParent;
        _barrierCodelets13987()
            : darts::Codelet()
        {
        }
        _barrierCodelets13987(uint32_t dep, uint32_t res, TP13987* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13989 : public darts::Codelet {
    public:
        TP13987* myTP;
        TP13987* inputsTPParent;
        _checkInCodelets13989()
            : darts::Codelet()
        {
        }
        _checkInCodelets13989(uint32_t dep, uint32_t res, TP13987* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13989 : public darts::Codelet {
    public:
        TP13987* inputsTPParent;
        _barrierCodelets13989()
            : darts::Codelet()
        {
        }
        _barrierCodelets13989(uint32_t dep, uint32_t res, TP13987* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14045 : public darts::Codelet {
    public:
        TP13987* myTP;
        TP13987* inputsTPParent;
        _checkInCodelets14045()
            : darts::Codelet()
        {
        }
        _checkInCodelets14045(uint32_t dep, uint32_t res, TP13987* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14044 : public darts::Codelet {
    public:
        TP13987* myTP;
        TP13987* inputsTPParent;
        _checkInCodelets14044()
            : darts::Codelet()
        {
        }
        _checkInCodelets14044(uint32_t dep, uint32_t res, TP13987* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14048 : public darts::Codelet {
    public:
        TP13987* myTP;
        TP13987* inputsTPParent;
        _checkInCodelets14048()
            : darts::Codelet()
        {
        }
        _checkInCodelets14048(uint32_t dep, uint32_t res, TP13987* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets14052 : public darts::Codelet {
    public:
        TP13987* inputsTPParent;
        _barrierCodelets14052()
            : darts::Codelet()
        {
        }
        _barrierCodelets14052(uint32_t dep, uint32_t res, TP13987* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14054 : public darts::Codelet {
    public:
        TP13987* myTP;
        TP13987* inputsTPParent;
        _checkInCodelets14054()
            : darts::Codelet()
        {
        }
        _checkInCodelets14054(uint32_t dep, uint32_t res, TP13987* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14053 : public darts::Codelet {
    public:
        TP13987* myTP;
        TP13987* inputsTPParent;
        _checkInCodelets14053()
            : darts::Codelet()
        {
        }
        _checkInCodelets14053(uint32_t dep, uint32_t res, TP13987* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14057 : public darts::Codelet {
    public:
        TP13987* myTP;
        TP13987* inputsTPParent;
        _checkInCodelets14057()
            : darts::Codelet()
        {
        }
        _checkInCodelets14057(uint32_t dep, uint32_t res, TP13987* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets14061 : public darts::Codelet {
    public:
        TP13987* inputsTPParent;
        _barrierCodelets14061()
            : darts::Codelet()
        {
        }
        _barrierCodelets14061(uint32_t dep, uint32_t res, TP13987* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14062 : public darts::Codelet {
    public:
        TP13987* myTP;
        TP13987* inputsTPParent;
        _checkInCodelets14062()
            : darts::Codelet()
        {
        }
        _checkInCodelets14062(uint32_t dep, uint32_t res, TP13987* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets14062 : public darts::Codelet {
    public:
        TP13987* inputsTPParent;
        _barrierCodelets14062()
            : darts::Codelet()
        {
        }
        _barrierCodelets14062(uint32_t dep, uint32_t res, TP13987* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP13987* TPParent;
    TP13987* controlTPParent;
    TP13987* inputsTPParent;
    int* i_darts13987 /*OMP_PRIVATE - INPUT*/;
    int* istep_darts13987 /*OMP_PRIVATE - INPUT*/;
    int* j_darts13987 /*OMP_PRIVATE - INPUT*/;
    int* k_darts13987 /*OMP_PRIVATE - INPUT*/;
    int* m_darts13987 /*OMP_PRIVATE - INPUT*/;
    double* tmp_darts13987; /*OMP_SHARED - INPUT*/
    TP13989** TP13989Ptr;
    size_t* TP13989_alreadyLaunched;
    int numTPsSet13989;
    int numTPsReady13989;
    size_t TPsToUse13989;
    size_t codeletsPerTP13989;
    size_t totalCodelets13989;
    unsigned int TP14044_LoopCounter;
    unsigned int* TP14044_LoopCounterPerThread;
    tbb::concurrent_vector<TP14044*> TP14044PtrVec;
    unsigned int TP14053_LoopCounter;
    unsigned int* TP14053_LoopCounterPerThread;
    tbb::concurrent_vector<TP14053*> TP14053PtrVec;
    TP14062** TP14062Ptr;
    size_t* TP14062_alreadyLaunched;
    int numTPsSet14062;
    int numTPsReady14062;
    size_t TPsToUse14062;
    size_t codeletsPerTP14062;
    size_t totalCodelets14062;
    _barrierCodelets13987* barrierCodelets13987;
    _checkInCodelets13989* checkInCodelets13989;
    _barrierCodelets13989* barrierCodelets13989;
    _checkInCodelets14045* checkInCodelets14045;
    _checkInCodelets14044* checkInCodelets14044;
    _checkInCodelets14048* checkInCodelets14048;
    _barrierCodelets14052* barrierCodelets14052;
    _checkInCodelets14054* checkInCodelets14054;
    _checkInCodelets14053* checkInCodelets14053;
    _checkInCodelets14057* checkInCodelets14057;
    _barrierCodelets14061* barrierCodelets14061;
    _checkInCodelets14062* checkInCodelets14062;
    _barrierCodelets14062* barrierCodelets14062;
    TP13987(
        int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet, double* in_tmp);
    ~TP13987();
};
/*TP13989: OMPForDirective*/
class TP13989 : public ompTP {
public:
    class _barrierCodelets13989 : public darts::Codelet {
    public:
        TP13989* inputsTPParent;
        _barrierCodelets13989()
            : darts::Codelet()
        {
        }
        _barrierCodelets13989(uint32_t dep, uint32_t res, TP13989* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13989(int* endRange, uint32_t codeletID);
    class _checkInCodelets13990 : public darts::Codelet {
    public:
        TP13989* myTP;
        TP13989* inputsTPParent;
        int endRange;
        _checkInCodelets13990()
            : darts::Codelet()
        {
        }
        _checkInCodelets13990(uint32_t dep, uint32_t res, TP13989* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13987* TPParent;
    TP13989* controlTPParent;
    TP13989* inputsTPParent;
    int* i_darts13989 /*OMP_PRIVATE - INPUT*/;
    int* j_darts13989 /*OMP_PRIVATE - INPUT*/;
    int* k_darts13989 /*OMP_PRIVATE - INPUT*/;
    int* m_darts13989 /*OMP_PRIVATE - INPUT*/;
    int initIteration13989;
    int lastIteration13989;
    int range13989;
    int rangePerCodelet13989;
    int minIteration13989;
    int remainderRange13989;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13989* barrierCodelets13989;
    _checkInCodelets13990* checkInCodelets13990;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13990* firstCodelet;
#endif
    TP13989(int in_numThreads, int in_mainCodeletID, TP13987* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13989** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13989();
};
/*TP14044: ForStmt*/
class TP14044 : public ompTP {
public:
    class _checkInCodelets14050 : public darts::Codelet {
    public:
        TP14044* myTP;
        TP13987* inputsTPParent;
        _checkInCodelets14050()
            : darts::Codelet()
        {
        }
        _checkInCodelets14050(uint32_t dep, uint32_t res, TP14044* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14051 : public darts::Codelet {
    public:
        TP14044* myTP;
        TP13987* inputsTPParent;
        _checkInCodelets14051()
            : darts::Codelet()
        {
        }
        _checkInCodelets14051(uint32_t dep, uint32_t res, TP14044* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13987* TPParent;
    TP14044* controlTPParent;
    TP13987* inputsTPParent;
    TP14044** ptrToThisTP;
    TP_jacld* TP14050Ptr;
    int TP14050_alreadyLaunched;
    TP_blts* TP14051Ptr;
    int TP14051_alreadyLaunched;
    _checkInCodelets14050* checkInCodelets14050;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets14050* firstCodelet;
#endif
    _checkInCodelets14051* checkInCodelets14051;
    TP14044(int in_numThreads, int in_mainCodeletID, TP13987* in_TPParent,
        TP13987* in_inputsTPParent, TP14044** in_ptrToThisTP);
    ~TP14044();
};
/*TP14053: ForStmt*/
class TP14053 : public ompTP {
public:
    class _checkInCodelets14059 : public darts::Codelet {
    public:
        TP14053* myTP;
        TP13987* inputsTPParent;
        _checkInCodelets14059()
            : darts::Codelet()
        {
        }
        _checkInCodelets14059(uint32_t dep, uint32_t res, TP14053* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14060 : public darts::Codelet {
    public:
        TP14053* myTP;
        TP13987* inputsTPParent;
        _checkInCodelets14060()
            : darts::Codelet()
        {
        }
        _checkInCodelets14060(uint32_t dep, uint32_t res, TP14053* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13987* TPParent;
    TP14053* controlTPParent;
    TP13987* inputsTPParent;
    TP14053** ptrToThisTP;
    TP_jacu* TP14059Ptr;
    int TP14059_alreadyLaunched;
    TP_buts* TP14060Ptr;
    int TP14060_alreadyLaunched;
    _checkInCodelets14059* checkInCodelets14059;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets14059* firstCodelet;
#endif
    _checkInCodelets14060* checkInCodelets14060;
    TP14053(int in_numThreads, int in_mainCodeletID, TP13987* in_TPParent,
        TP13987* in_inputsTPParent, TP14053** in_ptrToThisTP);
    ~TP14053();
};
/*TP14062: OMPForDirective*/
class TP14062 : public ompTP {
public:
    class _barrierCodelets14062 : public darts::Codelet {
    public:
        TP14062* inputsTPParent;
        _barrierCodelets14062()
            : darts::Codelet()
        {
        }
        _barrierCodelets14062(uint32_t dep, uint32_t res, TP14062* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations14062(int* endRange, uint32_t codeletID);
    class _checkInCodelets14063 : public darts::Codelet {
    public:
        TP14062* myTP;
        TP14062* inputsTPParent;
        int endRange;
        _checkInCodelets14063()
            : darts::Codelet()
        {
        }
        _checkInCodelets14063(uint32_t dep, uint32_t res, TP14062* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13987* TPParent;
    TP14062* controlTPParent;
    TP14062* inputsTPParent;
    int* i_darts14062 /*OMP_PRIVATE - INPUT*/;
    int* j_darts14062 /*OMP_PRIVATE - INPUT*/;
    int* k_darts14062 /*OMP_PRIVATE - INPUT*/;
    int* m_darts14062 /*OMP_PRIVATE - INPUT*/;
    double* tmp_darts14062; /*OMP_SHARED - INPUT*/
    int initIteration14062;
    int lastIteration14062;
    int range14062;
    int rangePerCodelet14062;
    int minIteration14062;
    int remainderRange14062;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets14062* barrierCodelets14062;
    _checkInCodelets14063* checkInCodelets14063;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets14063* firstCodelet;
#endif
    TP14062(int in_numThreads, int in_mainCodeletID, TP13987* in_TPParent, int in_initIteration,
        int in_lastIteration, double* in_tmp, TP14062** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP14062();
};
#endif
