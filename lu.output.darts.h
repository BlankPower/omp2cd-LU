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
class TP186;
class TP1;
typedef TP1 TP_blts;
class TP238;
/*Number of TPs to be used for the OMPFor in region TP238*/
#define NUMTPS238 NUMTPS
class TP344;
/*Number of TPs to be used for the OMPFor in region TP344*/
#define NUMTPS344 NUMTPS
class TP2;
typedef TP2 TP_buts;
class TP1228;
/*Number of TPs to be used for the OMPFor in region TP1228*/
#define NUMTPS1228 NUMTPS
class TP1330;
/*Number of TPs to be used for the OMPFor in region TP1330*/
#define NUMTPS1330 NUMTPS
class TP2241;
class TP2260;
/*Number of TPs to be used for the OMPFor in region TP2260*/
#define NUMTPS2260 NUMTPS
class TP2311;
/*Number of TPs to be used for the OMPFor in region TP2311*/
#define NUMTPS2311 NUMTPS
class TP2446;
/*Number of TPs to be used for the OMPFor in region TP2446*/
#define NUMTPS2446 NUMTPS
class TP2595;
/*Number of TPs to be used for the OMPFor in region TP2595*/
#define NUMTPS2595 NUMTPS
class TP3233;
/*Number of TPs to be used for the OMPFor in region TP3233*/
#define NUMTPS3233 NUMTPS
class TP3382;
/*Number of TPs to be used for the OMPFor in region TP3382*/
#define NUMTPS3382 NUMTPS
class TP4017;
/*Number of TPs to be used for the OMPFor in region TP4017*/
#define NUMTPS4017 NUMTPS
class TP7;
typedef TP7 TP_jacld;
class TP4943;
/*Number of TPs to be used for the OMPFor in region TP4943*/
#define NUMTPS4943 NUMTPS
class TP8;
typedef TP8 TP_jacu;
class TP7462;
/*Number of TPs to be used for the OMPFor in region TP7462*/
#define NUMTPS7462 NUMTPS
class TP9913;
class TP9931;
/*Number of TPs to be used for the OMPFor in region TP9931*/
#define NUMTPS9931 NUMTPS
class TP10830;
class TP10845;
/*Number of TPs to be used for the OMPFor in region TP10845*/
#define NUMTPS10845 NUMTPS
class TP10905;
/*Number of TPs to be used for the OMPFor in region TP10905*/
#define NUMTPS10905 NUMTPS
class TP11054;
/*Number of TPs to be used for the OMPFor in region TP11054*/
#define NUMTPS11054 NUMTPS
class TP11694;
/*Number of TPs to be used for the OMPFor in region TP11694*/
#define NUMTPS11694 NUMTPS
class TP11843;
/*Number of TPs to be used for the OMPFor in region TP11843*/
#define NUMTPS11843 NUMTPS
class TP12480;
/*Number of TPs to be used for the OMPFor in region TP12480*/
#define NUMTPS12480 NUMTPS
class TP13231;
class TP13235;
/*Number of TPs to be used for the OMPFor in region TP13235*/
#define NUMTPS13235 NUMTPS
class TP13286;
/*Number of TPs to be used for the OMPFor in region TP13286*/
#define NUMTPS13286 NUMTPS
class TP13328;
/*Number of TPs to be used for the OMPFor in region TP13328*/
#define NUMTPS13328 NUMTPS
class TP13372;
/*Number of TPs to be used for the OMPFor in region TP13372*/
#define NUMTPS13372 NUMTPS
class TP13414;
/*Number of TPs to be used for the OMPFor in region TP13414*/
#define NUMTPS13414 NUMTPS
class TP13903;
class TP13905;
/*Number of TPs to be used for the OMPFor in region TP13905*/
#define NUMTPS13905 NUMTPS
class TP14002;
class TP14004;
/*Number of TPs to be used for the OMPFor in region TP14004*/
#define NUMTPS14004 NUMTPS
class TP14059;
class TP14072;
class TP14081;
/*Number of TPs to be used for the OMPFor in region TP14081*/
#define NUMTPS14081 NUMTPS
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
/*TP186: OMPParallelDirective*/
class TP186 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets186 : public darts::Codelet {
    public:
        TP186* inputsTPParent;
        _barrierCodelets186()
            : darts::Codelet()
        {
        }
        _barrierCodelets186(uint32_t dep, uint32_t res, TP186* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets188 : public darts::Codelet {
    public:
        TP186* myTP;
        TP186* inputsTPParent;
        _checkInCodelets188()
            : darts::Codelet()
        {
        }
        _checkInCodelets188(uint32_t dep, uint32_t res, TP186* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP186* TPParent;
    TP186* controlTPParent;
    TP186* inputsTPParent;
    int* nthreads_darts186; /*OMP_SHARED - INPUT*/
    int* nthreads_darts188; /*OMP_SHARED - INPUT*/
    size_t TP188_alreadyLaunched;
    _barrierCodelets186* barrierCodelets186;
    _checkInCodelets188* checkInCodelets188;
    TP186(
        int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet, int* in_nthreads);
    ~TP186();
};
/*TP1: blts*/
class TP1 : public ompTP {
public:
    class _checkInCodelets238 : public darts::Codelet {
    public:
        TP1* myTP;
        TP1* inputsTPParent;
        _checkInCodelets238()
            : darts::Codelet()
        {
        }
        _checkInCodelets238(uint32_t dep, uint32_t res, TP1* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets238 : public darts::Codelet {
    public:
        TP1* inputsTPParent;
        _barrierCodelets238()
            : darts::Codelet()
        {
        }
        _barrierCodelets238(uint32_t dep, uint32_t res, TP1* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets342 : public darts::Codelet {
    public:
        TP1* myTP;
        TP1* inputsTPParent;
        _checkInCodelets342()
            : darts::Codelet()
        {
        }
        _checkInCodelets342(uint32_t dep, uint32_t res, TP1* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets344 : public darts::Codelet {
    public:
        TP1* myTP;
        TP1* inputsTPParent;
        _checkInCodelets344()
            : darts::Codelet()
        {
        }
        _checkInCodelets344(uint32_t dep, uint32_t res, TP1* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets344 : public darts::Codelet {
    public:
        TP1* inputsTPParent;
        _barrierCodelets344()
            : darts::Codelet()
        {
        }
        _barrierCodelets344(uint32_t dep, uint32_t res, TP1* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
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
    TP238** TP238Ptr;
    size_t* TP238_alreadyLaunched;
    int numTPsSet238;
    int numTPsReady238;
    size_t TPsToUse238;
    size_t codeletsPerTP238;
    size_t totalCodelets238;
    size_t TP342_alreadyLaunched;
    TP344** TP344Ptr;
    size_t* TP344_alreadyLaunched;
    int numTPsSet344;
    int numTPsReady344;
    size_t TPsToUse344;
    size_t codeletsPerTP344;
    size_t totalCodelets344;
    _checkInCodelets238* checkInCodelets238;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets238* firstCodelet;
#endif
    _barrierCodelets238* barrierCodelets238;
    _checkInCodelets342* checkInCodelets342;
    _checkInCodelets344* checkInCodelets344;
    _barrierCodelets344* barrierCodelets344;
    TP1(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
        darts::Codelet* in_mainSyncCodelet, TP1** in_ptrToThisFunctionTP, int in_nx, int in_ny,
        int in_nz, int in_k, double in_omega, int in_ist, int in_iend, int in_jst, int in_jend,
        int in_nx0, int in_ny0);
    ~TP1();
    void setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, int in_ist,
        int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID);
};
/*TP238: OMPForDirective*/
class TP238 : public ompTP {
public:
    class _barrierCodelets238 : public darts::Codelet {
    public:
        TP238* inputsTPParent;
        _barrierCodelets238()
            : darts::Codelet()
        {
        }
        _barrierCodelets238(uint32_t dep, uint32_t res, TP238* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations238(int* endRange, uint32_t codeletID);
    class _checkInCodelets239 : public darts::Codelet {
    public:
        TP238* myTP;
        TP238* inputsTPParent;
        int endRange;
        _checkInCodelets239()
            : darts::Codelet()
        {
        }
        _checkInCodelets239(uint32_t dep, uint32_t res, TP238* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP1* TPParent;
    TP238* controlTPParent;
    TP238* inputsTPParent;
    int* i_darts238 /*OMP_PRIVATE - INPUT*/;
    int** iend_darts238 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** ist_darts238 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts238 /*OMP_PRIVATE - INPUT*/;
    int** jend_darts238 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** jst_darts238 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** k_darts238 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* m_darts238 /*OMP_PRIVATE - INPUT*/;
    double** omega_darts238 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration238;
    int lastIteration238;
    int range238;
    int rangePerCodelet238;
    int minIteration238;
    int remainderRange238;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets238* barrierCodelets238;
    _checkInCodelets239* checkInCodelets239;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets239* firstCodelet;
#endif
    TP238(int in_numThreads, int in_mainCodeletID, TP1* in_TPParent, int in_initIteration,
        int in_lastIteration, TP238** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP238();
};
/*TP344: OMPForDirective*/
class TP344 : public ompTP {
public:
    class _barrierCodelets344 : public darts::Codelet {
    public:
        TP344* inputsTPParent;
        _barrierCodelets344()
            : darts::Codelet()
        {
        }
        _barrierCodelets344(uint32_t dep, uint32_t res, TP344* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations344(int* endRange, uint32_t codeletID);
    class _checkInCodelets345 : public darts::Codelet {
    public:
        TP344* myTP;
        TP344* inputsTPParent;
        int endRange;
        _checkInCodelets345()
            : darts::Codelet()
        {
        }
        _checkInCodelets345(uint32_t dep, uint32_t res, TP344* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP1* TPParent;
    TP344* controlTPParent;
    TP344* inputsTPParent;
    int* i_darts344 /*OMP_PRIVATE - INPUT*/;
    int** iend_darts344 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** ist_darts344 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts344 /*OMP_PRIVATE - INPUT*/;
    int** jend_darts344 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** jst_darts344 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** k_darts344 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* m_darts344 /*OMP_PRIVATE - INPUT*/;
    double** omega_darts344 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp_darts344 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp1_darts344 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration344;
    int lastIteration344;
    int range344;
    int rangePerCodelet344;
    int minIteration344;
    int remainderRange344;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets344* barrierCodelets344;
    _checkInCodelets345* checkInCodelets345;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets345* firstCodelet;
#endif
    TP344(int in_numThreads, int in_mainCodeletID, TP1* in_TPParent, int in_initIteration,
        int in_lastIteration, TP344** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP344();
};
/*TP2: buts*/
class TP2 : public ompTP {
public:
    class _checkInCodelets1228 : public darts::Codelet {
    public:
        TP2* myTP;
        TP2* inputsTPParent;
        _checkInCodelets1228()
            : darts::Codelet()
        {
        }
        _checkInCodelets1228(uint32_t dep, uint32_t res, TP2* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets1330 : public darts::Codelet {
    public:
        TP2* myTP;
        TP2* inputsTPParent;
        _checkInCodelets1330()
            : darts::Codelet()
        {
        }
        _checkInCodelets1330(uint32_t dep, uint32_t res, TP2* myTP, uint32_t id)
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
    TP1228** TP1228Ptr;
    size_t* TP1228_alreadyLaunched;
    int numTPsSet1228;
    int numTPsReady1228;
    size_t TPsToUse1228;
    size_t codeletsPerTP1228;
    size_t totalCodelets1228;
    TP1330** TP1330Ptr;
    size_t* TP1330_alreadyLaunched;
    int numTPsSet1330;
    int numTPsReady1330;
    size_t TPsToUse1330;
    size_t codeletsPerTP1330;
    size_t totalCodelets1330;
    _checkInCodelets1228* checkInCodelets1228;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets1228* firstCodelet;
#endif
    _checkInCodelets1330* checkInCodelets1330;
    TP2(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
        darts::Codelet* in_mainSyncCodelet, TP2** in_ptrToThisFunctionTP, int in_nx, int in_ny,
        int in_nz, int in_k, double in_omega, int in_ist, int in_iend, int in_jst, int in_jend,
        int in_nx0, int in_ny0);
    ~TP2();
    void setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, int in_ist,
        int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID);
};
/*TP1228: OMPForDirective*/
class TP1228 : public ompTP {
public:
    bool requestNewRangeIterations1228(int* endRange, uint32_t codeletID);
    class _checkInCodelets1229 : public darts::Codelet {
    public:
        TP1228* myTP;
        TP1228* inputsTPParent;
        int endRange;
        _checkInCodelets1229()
            : darts::Codelet()
        {
        }
        _checkInCodelets1229(uint32_t dep, uint32_t res, TP1228* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2* TPParent;
    TP1228* controlTPParent;
    TP1228* inputsTPParent;
    int* i_darts1228 /*OMP_PRIVATE - INPUT*/;
    int** iend_darts1228 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** ist_darts1228 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts1228 /*OMP_PRIVATE - INPUT*/;
    int** jend_darts1228 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** jst_darts1228 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** k_darts1228 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* m_darts1228 /*OMP_PRIVATE - INPUT*/;
    double** omega_darts1228 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration1228;
    int lastIteration1228;
    int range1228;
    int rangePerCodelet1228;
    int minIteration1228;
    int remainderRange1228;
    size_t readyCodelets;
    int baseNumThreads;
    int* signalNextReady;
    _checkInCodelets1229* checkInCodelets1229;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets1229* firstCodelet;
#endif
    TP1228(int in_numThreads, int in_mainCodeletID, TP2* in_TPParent, int in_initIteration,
        int in_lastIteration, TP1228** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP1228();
};
/*TP1330: OMPForDirective*/
class TP1330 : public ompTP {
public:
    bool requestNewRangeIterations1330(int* endRange, uint32_t codeletID);
    class _checkInCodelets1331 : public darts::Codelet {
    public:
        TP1330* myTP;
        TP1330* inputsTPParent;
        int endRange;
        _checkInCodelets1331()
            : darts::Codelet()
        {
        }
        _checkInCodelets1331(uint32_t dep, uint32_t res, TP1330* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2* TPParent;
    TP1330* controlTPParent;
    TP1330* inputsTPParent;
    int* i_darts1330 /*OMP_PRIVATE - INPUT*/;
    int** iend_darts1330 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** ist_darts1330 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts1330 /*OMP_PRIVATE - INPUT*/;
    int** jend_darts1330 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** jst_darts1330 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** k_darts1330 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* m_darts1330 /*OMP_PRIVATE - INPUT*/;
    double** omega_darts1330 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp_darts1330 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp1_darts1330 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration1330;
    int lastIteration1330;
    int range1330;
    int rangePerCodelet1330;
    int minIteration1330;
    int remainderRange1330;
    size_t readyCodelets;
    int baseNumThreads;
    int* signalNextReady;
    _checkInCodelets1331* checkInCodelets1331;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets1331* firstCodelet;
#endif
    TP1330(int in_numThreads, int in_mainCodeletID, TP2* in_TPParent, int in_initIteration,
        int in_lastIteration, TP1330** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP1330();
};
/*TP2241: OMPParallelDirective*/
class TP2241 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets2241 : public darts::Codelet {
    public:
        TP2241* inputsTPParent;
        _barrierCodelets2241()
            : darts::Codelet()
        {
        }
        _barrierCodelets2241(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2243 : public darts::Codelet {
    public:
        TP2241* myTP;
        TP2241* inputsTPParent;
        _checkInCodelets2243()
            : darts::Codelet()
        {
        }
        _checkInCodelets2243(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2260 : public darts::Codelet {
    public:
        TP2241* myTP;
        TP2241* inputsTPParent;
        _checkInCodelets2260()
            : darts::Codelet()
        {
        }
        _checkInCodelets2260(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets2260 : public darts::Codelet {
    public:
        TP2241* inputsTPParent;
        _barrierCodelets2260()
            : darts::Codelet()
        {
        }
        _barrierCodelets2260(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2311 : public darts::Codelet {
    public:
        TP2241* myTP;
        TP2241* inputsTPParent;
        _checkInCodelets2311()
            : darts::Codelet()
        {
        }
        _checkInCodelets2311(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets2311 : public darts::Codelet {
    public:
        TP2241* inputsTPParent;
        _barrierCodelets2311()
            : darts::Codelet()
        {
        }
        _barrierCodelets2311(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2443 : public darts::Codelet {
    public:
        TP2241* myTP;
        TP2241* inputsTPParent;
        _checkInCodelets2443()
            : darts::Codelet()
        {
        }
        _checkInCodelets2443(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2446 : public darts::Codelet {
    public:
        TP2241* myTP;
        TP2241* inputsTPParent;
        _checkInCodelets2446()
            : darts::Codelet()
        {
        }
        _checkInCodelets2446(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets2446 : public darts::Codelet {
    public:
        TP2241* inputsTPParent;
        _barrierCodelets2446()
            : darts::Codelet()
        {
        }
        _barrierCodelets2446(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2595 : public darts::Codelet {
    public:
        TP2241* myTP;
        TP2241* inputsTPParent;
        _checkInCodelets2595()
            : darts::Codelet()
        {
        }
        _checkInCodelets2595(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets2595 : public darts::Codelet {
    public:
        TP2241* inputsTPParent;
        _barrierCodelets2595()
            : darts::Codelet()
        {
        }
        _barrierCodelets2595(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets3230 : public darts::Codelet {
    public:
        TP2241* myTP;
        TP2241* inputsTPParent;
        _checkInCodelets3230()
            : darts::Codelet()
        {
        }
        _checkInCodelets3230(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets3233 : public darts::Codelet {
    public:
        TP2241* myTP;
        TP2241* inputsTPParent;
        _checkInCodelets3233()
            : darts::Codelet()
        {
        }
        _checkInCodelets3233(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets3233 : public darts::Codelet {
    public:
        TP2241* inputsTPParent;
        _barrierCodelets3233()
            : darts::Codelet()
        {
        }
        _barrierCodelets3233(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets3382 : public darts::Codelet {
    public:
        TP2241* myTP;
        TP2241* inputsTPParent;
        _checkInCodelets3382()
            : darts::Codelet()
        {
        }
        _checkInCodelets3382(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets3382 : public darts::Codelet {
    public:
        TP2241* inputsTPParent;
        _barrierCodelets3382()
            : darts::Codelet()
        {
        }
        _barrierCodelets3382(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets4017 : public darts::Codelet {
    public:
        TP2241* myTP;
        TP2241* inputsTPParent;
        _checkInCodelets4017()
            : darts::Codelet()
        {
        }
        _checkInCodelets4017(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets4017 : public darts::Codelet {
    public:
        TP2241* inputsTPParent;
        _barrierCodelets4017()
            : darts::Codelet()
        {
        }
        _barrierCodelets4017(uint32_t dep, uint32_t res, TP2241* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP2241* TPParent;
    TP2241* controlTPParent;
    TP2241* inputsTPParent;
    int* L1_darts2241 /*VARIABLE*/;
    int* L2_darts2241 /*VARIABLE*/;
    double* dsspm_darts2241 /*VARIABLE*/;
    double* eta_darts2241 /*VARIABLE*/;
    int* i_darts2241 /*VARIABLE*/;
    int* iend1_darts2241 /*VARIABLE*/;
    int* iglob_darts2241 /*VARIABLE*/;
    int* ist1_darts2241 /*VARIABLE*/;
    int* j_darts2241 /*VARIABLE*/;
    int* jend1_darts2241 /*VARIABLE*/;
    int* jglob_darts2241 /*VARIABLE*/;
    int* jst1_darts2241 /*VARIABLE*/;
    int* k_darts2241 /*VARIABLE*/;
    int* m_darts2241 /*VARIABLE*/;
    double* q_darts2241 /*VARIABLE*/;
    double* tmp_darts2241 /*VARIABLE*/;
    double* u21_darts2241 /*VARIABLE*/;
    double* u21i_darts2241 /*VARIABLE*/;
    double* u21im1_darts2241 /*VARIABLE*/;
    double* u21j_darts2241 /*VARIABLE*/;
    double* u21jm1_darts2241 /*VARIABLE*/;
    double* u21k_darts2241 /*VARIABLE*/;
    double* u21km1_darts2241 /*VARIABLE*/;
    double* u31_darts2241 /*VARIABLE*/;
    double* u31i_darts2241 /*VARIABLE*/;
    double* u31im1_darts2241 /*VARIABLE*/;
    double* u31j_darts2241 /*VARIABLE*/;
    double* u31jm1_darts2241 /*VARIABLE*/;
    double* u31k_darts2241 /*VARIABLE*/;
    double* u31km1_darts2241 /*VARIABLE*/;
    double* u41_darts2241 /*VARIABLE*/;
    double* u41i_darts2241 /*VARIABLE*/;
    double* u41im1_darts2241 /*VARIABLE*/;
    double* u41j_darts2241 /*VARIABLE*/;
    double* u41jm1_darts2241 /*VARIABLE*/;
    double* u41k_darts2241 /*VARIABLE*/;
    double* u41km1_darts2241 /*VARIABLE*/;
    double* u51i_darts2241 /*VARIABLE*/;
    double* u51im1_darts2241 /*VARIABLE*/;
    double* u51j_darts2241 /*VARIABLE*/;
    double* u51jm1_darts2241 /*VARIABLE*/;
    double* u51k_darts2241 /*VARIABLE*/;
    double* u51km1_darts2241 /*VARIABLE*/;
    double* xi_darts2241 /*VARIABLE*/;
    double* zeta_darts2241 /*VARIABLE*/;
    TP2260** TP2260Ptr;
    size_t* TP2260_alreadyLaunched;
    int numTPsSet2260;
    int numTPsReady2260;
    size_t TPsToUse2260;
    size_t codeletsPerTP2260;
    size_t totalCodelets2260;
    TP2311** TP2311Ptr;
    size_t* TP2311_alreadyLaunched;
    int numTPsSet2311;
    int numTPsReady2311;
    size_t TPsToUse2311;
    size_t codeletsPerTP2311;
    size_t totalCodelets2311;
    TP2446** TP2446Ptr;
    size_t* TP2446_alreadyLaunched;
    int numTPsSet2446;
    int numTPsReady2446;
    size_t TPsToUse2446;
    size_t codeletsPerTP2446;
    size_t totalCodelets2446;
    TP2595** TP2595Ptr;
    size_t* TP2595_alreadyLaunched;
    int numTPsSet2595;
    int numTPsReady2595;
    size_t TPsToUse2595;
    size_t codeletsPerTP2595;
    size_t totalCodelets2595;
    TP3233** TP3233Ptr;
    size_t* TP3233_alreadyLaunched;
    int numTPsSet3233;
    int numTPsReady3233;
    size_t TPsToUse3233;
    size_t codeletsPerTP3233;
    size_t totalCodelets3233;
    TP3382** TP3382Ptr;
    size_t* TP3382_alreadyLaunched;
    int numTPsSet3382;
    int numTPsReady3382;
    size_t TPsToUse3382;
    size_t codeletsPerTP3382;
    size_t totalCodelets3382;
    TP4017** TP4017Ptr;
    size_t* TP4017_alreadyLaunched;
    int numTPsSet4017;
    int numTPsReady4017;
    size_t TPsToUse4017;
    size_t codeletsPerTP4017;
    size_t totalCodelets4017;
    _barrierCodelets2241* barrierCodelets2241;
    _checkInCodelets2243* checkInCodelets2243;
    _checkInCodelets2260* checkInCodelets2260;
    _barrierCodelets2260* barrierCodelets2260;
    _checkInCodelets2311* checkInCodelets2311;
    _barrierCodelets2311* barrierCodelets2311;
    _checkInCodelets2443* checkInCodelets2443;
    _checkInCodelets2446* checkInCodelets2446;
    _barrierCodelets2446* barrierCodelets2446;
    _checkInCodelets2595* checkInCodelets2595;
    _barrierCodelets2595* barrierCodelets2595;
    _checkInCodelets3230* checkInCodelets3230;
    _checkInCodelets3233* checkInCodelets3233;
    _barrierCodelets3233* barrierCodelets3233;
    _checkInCodelets3382* checkInCodelets3382;
    _barrierCodelets3382* barrierCodelets3382;
    _checkInCodelets4017* checkInCodelets4017;
    _barrierCodelets4017* barrierCodelets4017;
    TP2241(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
    ~TP2241();
};
/*TP2260: OMPForDirective*/
class TP2260 : public ompTP {
public:
    class _barrierCodelets2260 : public darts::Codelet {
    public:
        TP2260* inputsTPParent;
        _barrierCodelets2260()
            : darts::Codelet()
        {
        }
        _barrierCodelets2260(uint32_t dep, uint32_t res, TP2260* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations2260(int* endRange, uint32_t codeletID);
    class _checkInCodelets2261 : public darts::Codelet {
    public:
        TP2260* myTP;
        TP2260* inputsTPParent;
        int endRange;
        _checkInCodelets2261()
            : darts::Codelet()
        {
        }
        _checkInCodelets2261(uint32_t dep, uint32_t res, TP2260* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2241* TPParent;
    TP2260* controlTPParent;
    TP2260* inputsTPParent;
    int* i_darts2260 /*OMP_PRIVATE - INPUT*/;
    int* j_darts2260 /*OMP_PRIVATE - INPUT*/;
    int* k_darts2260 /*OMP_PRIVATE - INPUT*/;
    int* m_darts2260 /*OMP_PRIVATE - INPUT*/;
    int initIteration2260;
    int lastIteration2260;
    int range2260;
    int rangePerCodelet2260;
    int minIteration2260;
    int remainderRange2260;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets2260* barrierCodelets2260;
    _checkInCodelets2261* checkInCodelets2261;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets2261* firstCodelet;
#endif
    TP2260(int in_numThreads, int in_mainCodeletID, TP2241* in_TPParent, int in_initIteration,
        int in_lastIteration, TP2260** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP2260();
};
/*TP2311: OMPForDirective*/
class TP2311 : public ompTP {
public:
    class _barrierCodelets2311 : public darts::Codelet {
    public:
        TP2311* inputsTPParent;
        _barrierCodelets2311()
            : darts::Codelet()
        {
        }
        _barrierCodelets2311(uint32_t dep, uint32_t res, TP2311* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations2311(int* endRange, uint32_t codeletID);
    class _checkInCodelets2312 : public darts::Codelet {
    public:
        TP2311* myTP;
        TP2311* inputsTPParent;
        int endRange;
        _checkInCodelets2312()
            : darts::Codelet()
        {
        }
        _checkInCodelets2312(uint32_t dep, uint32_t res, TP2311* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2241* TPParent;
    TP2311* controlTPParent;
    TP2311* inputsTPParent;
    double** eta_darts2311 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts2311 /*OMP_PRIVATE - INPUT*/;
    int** iglob_darts2311 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts2311 /*OMP_PRIVATE - INPUT*/;
    int** jglob_darts2311 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts2311 /*OMP_PRIVATE - INPUT*/;
    int* m_darts2311 /*OMP_PRIVATE - INPUT*/;
    double** xi_darts2311 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** zeta_darts2311 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration2311;
    int lastIteration2311;
    int range2311;
    int rangePerCodelet2311;
    int minIteration2311;
    int remainderRange2311;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets2311* barrierCodelets2311;
    _checkInCodelets2312* checkInCodelets2312;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets2312* firstCodelet;
#endif
    TP2311(int in_numThreads, int in_mainCodeletID, TP2241* in_TPParent, int in_initIteration,
        int in_lastIteration, TP2311** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP2311();
};
/*TP2446: OMPForDirective*/
class TP2446 : public ompTP {
public:
    class _barrierCodelets2446 : public darts::Codelet {
    public:
        TP2446* inputsTPParent;
        _barrierCodelets2446()
            : darts::Codelet()
        {
        }
        _barrierCodelets2446(uint32_t dep, uint32_t res, TP2446* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations2446(int* endRange, uint32_t codeletID);
    class _checkInCodelets2447 : public darts::Codelet {
    public:
        TP2446* myTP;
        TP2446* inputsTPParent;
        int endRange;
        _checkInCodelets2447()
            : darts::Codelet()
        {
        }
        _checkInCodelets2447(uint32_t dep, uint32_t res, TP2446* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2241* TPParent;
    TP2446* controlTPParent;
    TP2446* inputsTPParent;
    int** L1_darts2446 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** L2_darts2446 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts2446 /*OMP_PRIVATE - INPUT*/;
    int* j_darts2446 /*OMP_PRIVATE - INPUT*/;
    int* k_darts2446 /*OMP_PRIVATE - INPUT*/;
    double** q_darts2446 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21_darts2446 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration2446;
    int lastIteration2446;
    int range2446;
    int rangePerCodelet2446;
    int minIteration2446;
    int remainderRange2446;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets2446* barrierCodelets2446;
    _checkInCodelets2447* checkInCodelets2447;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets2447* firstCodelet;
#endif
    TP2446(int in_numThreads, int in_mainCodeletID, TP2241* in_TPParent, int in_initIteration,
        int in_lastIteration, TP2446** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP2446();
};
/*TP2595: OMPForDirective*/
class TP2595 : public ompTP {
public:
    class _barrierCodelets2595 : public darts::Codelet {
    public:
        TP2595* inputsTPParent;
        _barrierCodelets2595()
            : darts::Codelet()
        {
        }
        _barrierCodelets2595(uint32_t dep, uint32_t res, TP2595* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations2595(int* endRange, uint32_t codeletID);
    class _checkInCodelets2596 : public darts::Codelet {
    public:
        TP2595* myTP;
        TP2595* inputsTPParent;
        int endRange;
        _checkInCodelets2596()
            : darts::Codelet()
        {
        }
        _checkInCodelets2596(uint32_t dep, uint32_t res, TP2595* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2241* TPParent;
    TP2595* controlTPParent;
    TP2595* inputsTPParent;
    int** L2_darts2595 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** dsspm_darts2595 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts2595 /*OMP_PRIVATE - INPUT*/;
    int** iend1_darts2595 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** ist1_darts2595 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts2595 /*OMP_PRIVATE - INPUT*/;
    int* k_darts2595 /*OMP_PRIVATE - INPUT*/;
    int* m_darts2595 /*OMP_PRIVATE - INPUT*/;
    double** tmp_darts2595 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21i_darts2595 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21im1_darts2595 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31i_darts2595 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31im1_darts2595 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41i_darts2595 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41im1_darts2595 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51i_darts2595 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51im1_darts2595 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration2595;
    int lastIteration2595;
    int range2595;
    int rangePerCodelet2595;
    int minIteration2595;
    int remainderRange2595;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets2595* barrierCodelets2595;
    _checkInCodelets2596* checkInCodelets2596;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets2596* firstCodelet;
#endif
    TP2595(int in_numThreads, int in_mainCodeletID, TP2241* in_TPParent, int in_initIteration,
        int in_lastIteration, TP2595** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP2595();
};
/*TP3233: OMPForDirective*/
class TP3233 : public ompTP {
public:
    class _barrierCodelets3233 : public darts::Codelet {
    public:
        TP3233* inputsTPParent;
        _barrierCodelets3233()
            : darts::Codelet()
        {
        }
        _barrierCodelets3233(uint32_t dep, uint32_t res, TP3233* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations3233(int* endRange, uint32_t codeletID);
    class _checkInCodelets3234 : public darts::Codelet {
    public:
        TP3233* myTP;
        TP3233* inputsTPParent;
        int endRange;
        _checkInCodelets3234()
            : darts::Codelet()
        {
        }
        _checkInCodelets3234(uint32_t dep, uint32_t res, TP3233* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2241* TPParent;
    TP3233* controlTPParent;
    TP3233* inputsTPParent;
    int** L1_darts3233 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** L2_darts3233 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts3233 /*OMP_PRIVATE - INPUT*/;
    int* j_darts3233 /*OMP_PRIVATE - INPUT*/;
    int* k_darts3233 /*OMP_PRIVATE - INPUT*/;
    double** q_darts3233 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31_darts3233 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration3233;
    int lastIteration3233;
    int range3233;
    int rangePerCodelet3233;
    int minIteration3233;
    int remainderRange3233;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets3233* barrierCodelets3233;
    _checkInCodelets3234* checkInCodelets3234;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets3234* firstCodelet;
#endif
    TP3233(int in_numThreads, int in_mainCodeletID, TP2241* in_TPParent, int in_initIteration,
        int in_lastIteration, TP3233** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP3233();
};
/*TP3382: OMPForDirective*/
class TP3382 : public ompTP {
public:
    class _barrierCodelets3382 : public darts::Codelet {
    public:
        TP3382* inputsTPParent;
        _barrierCodelets3382()
            : darts::Codelet()
        {
        }
        _barrierCodelets3382(uint32_t dep, uint32_t res, TP3382* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations3382(int* endRange, uint32_t codeletID);
    class _checkInCodelets3383 : public darts::Codelet {
    public:
        TP3382* myTP;
        TP3382* inputsTPParent;
        int endRange;
        _checkInCodelets3383()
            : darts::Codelet()
        {
        }
        _checkInCodelets3383(uint32_t dep, uint32_t res, TP3382* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2241* TPParent;
    TP3382* controlTPParent;
    TP3382* inputsTPParent;
    int** L2_darts3382 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** dsspm_darts3382 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts3382 /*OMP_PRIVATE - INPUT*/;
    int* j_darts3382 /*OMP_PRIVATE - INPUT*/;
    int** jend1_darts3382 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** jst1_darts3382 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts3382 /*OMP_PRIVATE - INPUT*/;
    int* m_darts3382 /*OMP_PRIVATE - INPUT*/;
    double** tmp_darts3382 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21j_darts3382 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21jm1_darts3382 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31j_darts3382 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31jm1_darts3382 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41j_darts3382 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41jm1_darts3382 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51j_darts3382 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51jm1_darts3382 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration3382;
    int lastIteration3382;
    int range3382;
    int rangePerCodelet3382;
    int minIteration3382;
    int remainderRange3382;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets3382* barrierCodelets3382;
    _checkInCodelets3383* checkInCodelets3383;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets3383* firstCodelet;
#endif
    TP3382(int in_numThreads, int in_mainCodeletID, TP2241* in_TPParent, int in_initIteration,
        int in_lastIteration, TP3382** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP3382();
};
/*TP4017: OMPForDirective*/
class TP4017 : public ompTP {
public:
    class _barrierCodelets4017 : public darts::Codelet {
    public:
        TP4017* inputsTPParent;
        _barrierCodelets4017()
            : darts::Codelet()
        {
        }
        _barrierCodelets4017(uint32_t dep, uint32_t res, TP4017* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations4017(int* endRange, uint32_t codeletID);
    class _checkInCodelets4018 : public darts::Codelet {
    public:
        TP4017* myTP;
        TP4017* inputsTPParent;
        int endRange;
        _checkInCodelets4018()
            : darts::Codelet()
        {
        }
        _checkInCodelets4018(uint32_t dep, uint32_t res, TP4017* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2241* TPParent;
    TP4017* controlTPParent;
    TP4017* inputsTPParent;
    double** dsspm_darts4017 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts4017 /*OMP_PRIVATE - INPUT*/;
    int* j_darts4017 /*OMP_PRIVATE - INPUT*/;
    int* k_darts4017 /*OMP_PRIVATE - INPUT*/;
    int* m_darts4017 /*OMP_PRIVATE - INPUT*/;
    double** q_darts4017 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp_darts4017 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21k_darts4017 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21km1_darts4017 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31k_darts4017 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31km1_darts4017 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41_darts4017 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41k_darts4017 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41km1_darts4017 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51k_darts4017 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51km1_darts4017 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration4017;
    int lastIteration4017;
    int range4017;
    int rangePerCodelet4017;
    int minIteration4017;
    int remainderRange4017;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets4017* barrierCodelets4017;
    _checkInCodelets4018* checkInCodelets4018;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets4018* firstCodelet;
#endif
    TP4017(int in_numThreads, int in_mainCodeletID, TP2241* in_TPParent, int in_initIteration,
        int in_lastIteration, TP4017** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP4017();
};
/*TP7: jacld*/
class TP7 : public ompTP {
public:
    class _checkInCodelets4922 : public darts::Codelet {
    public:
        TP7* myTP;
        TP7* inputsTPParent;
        _checkInCodelets4922()
            : darts::Codelet()
        {
        }
        _checkInCodelets4922(uint32_t dep, uint32_t res, TP7* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets4943 : public darts::Codelet {
    public:
        TP7* myTP;
        TP7* inputsTPParent;
        _checkInCodelets4943()
            : darts::Codelet()
        {
        }
        _checkInCodelets4943(uint32_t dep, uint32_t res, TP7* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets4943 : public darts::Codelet {
    public:
        TP7* inputsTPParent;
        _barrierCodelets4943()
            : darts::Codelet()
        {
        }
        _barrierCodelets4943(uint32_t dep, uint32_t res, TP7* myTP, uint32_t id)
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
    TP4943** TP4943Ptr;
    size_t* TP4943_alreadyLaunched;
    int numTPsSet4943;
    int numTPsReady4943;
    size_t TPsToUse4943;
    size_t codeletsPerTP4943;
    size_t totalCodelets4943;
    _checkInCodelets4922* checkInCodelets4922;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets4922* firstCodelet;
#endif
    _checkInCodelets4943* checkInCodelets4943;
    _barrierCodelets4943* barrierCodelets4943;
    TP7(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
        darts::Codelet* in_mainSyncCodelet, TP7** in_ptrToThisFunctionTP, int in_k);
    ~TP7();
    void setNewInputs(int in_k, size_t codeletID);
};
/*TP4943: OMPForDirective*/
class TP4943 : public ompTP {
public:
    class _barrierCodelets4943 : public darts::Codelet {
    public:
        TP4943* inputsTPParent;
        _barrierCodelets4943()
            : darts::Codelet()
        {
        }
        _barrierCodelets4943(uint32_t dep, uint32_t res, TP4943* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations4943(int* endRange, uint32_t codeletID);
    class _checkInCodelets4944 : public darts::Codelet {
    public:
        TP4943* myTP;
        TP4943* inputsTPParent;
        int endRange;
        _checkInCodelets4944()
            : darts::Codelet()
        {
        }
        _checkInCodelets4944(uint32_t dep, uint32_t res, TP4943* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP7* TPParent;
    TP4943* controlTPParent;
    TP4943* inputsTPParent;
    double** c1345_darts4943 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** c34_darts4943 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts4943 /*OMP_PRIVATE - INPUT*/;
    int* j_darts4943 /*OMP_PRIVATE - INPUT*/;
    int** k_darts4943 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** r43_darts4943 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp1_darts4943 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp2_darts4943 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp3_darts4943 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration4943;
    int lastIteration4943;
    int range4943;
    int rangePerCodelet4943;
    int minIteration4943;
    int remainderRange4943;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets4943* barrierCodelets4943;
    _checkInCodelets4944* checkInCodelets4944;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets4944* firstCodelet;
#endif
    TP4943(int in_numThreads, int in_mainCodeletID, TP7* in_TPParent, int in_initIteration,
        int in_lastIteration, TP4943** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP4943();
};
/*TP8: jacu*/
class TP8 : public ompTP {
public:
    class _checkInCodelets7441 : public darts::Codelet {
    public:
        TP8* myTP;
        TP8* inputsTPParent;
        _checkInCodelets7441()
            : darts::Codelet()
        {
        }
        _checkInCodelets7441(uint32_t dep, uint32_t res, TP8* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets7462 : public darts::Codelet {
    public:
        TP8* myTP;
        TP8* inputsTPParent;
        _checkInCodelets7462()
            : darts::Codelet()
        {
        }
        _checkInCodelets7462(uint32_t dep, uint32_t res, TP8* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
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
    TP7462** TP7462Ptr;
    size_t* TP7462_alreadyLaunched;
    int numTPsSet7462;
    int numTPsReady7462;
    size_t TPsToUse7462;
    size_t codeletsPerTP7462;
    size_t totalCodelets7462;
    _checkInCodelets7441* checkInCodelets7441;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets7441* firstCodelet;
#endif
    _checkInCodelets7462* checkInCodelets7462;
    TP8(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
        darts::Codelet* in_mainSyncCodelet, TP8** in_ptrToThisFunctionTP, int in_k);
    ~TP8();
    void setNewInputs(int in_k, size_t codeletID);
};
/*TP7462: OMPForDirective*/
class TP7462 : public ompTP {
public:
    bool requestNewRangeIterations7462(int* endRange, uint32_t codeletID);
    class _checkInCodelets7463 : public darts::Codelet {
    public:
        TP7462* myTP;
        TP7462* inputsTPParent;
        int endRange;
        _checkInCodelets7463()
            : darts::Codelet()
        {
        }
        _checkInCodelets7463(uint32_t dep, uint32_t res, TP7462* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP8* TPParent;
    TP7462* controlTPParent;
    TP7462* inputsTPParent;
    double** c1345_darts7462 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** c34_darts7462 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts7462 /*OMP_PRIVATE - INPUT*/;
    int* j_darts7462 /*OMP_PRIVATE - INPUT*/;
    int** k_darts7462 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** r43_darts7462 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp1_darts7462 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp2_darts7462 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp3_darts7462 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration7462;
    int lastIteration7462;
    int range7462;
    int rangePerCodelet7462;
    int minIteration7462;
    int remainderRange7462;
    size_t readyCodelets;
    int baseNumThreads;
    int* signalNextReady;
    _checkInCodelets7463* checkInCodelets7463;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets7463* firstCodelet;
#endif
    TP7462(int in_numThreads, int in_mainCodeletID, TP8* in_TPParent, int in_initIteration,
        int in_lastIteration, TP7462** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP7462();
};
/*TP9913: OMPParallelDirective*/
class TP9913 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets9913 : public darts::Codelet {
    public:
        TP9913* inputsTPParent;
        _barrierCodelets9913()
            : darts::Codelet()
        {
        }
        _barrierCodelets9913(uint32_t dep, uint32_t res, TP9913* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets9915 : public darts::Codelet {
    public:
        TP9913* myTP;
        TP9913* inputsTPParent;
        _checkInCodelets9915()
            : darts::Codelet()
        {
        }
        _checkInCodelets9915(uint32_t dep, uint32_t res, TP9913* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets9922 : public darts::Codelet {
    public:
        TP9913* myTP;
        TP9913* inputsTPParent;
        _checkInCodelets9922()
            : darts::Codelet()
        {
        }
        _checkInCodelets9922(uint32_t dep, uint32_t res, TP9913* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets9922 : public darts::Codelet {
    public:
        TP9913* inputsTPParent;
        _barrierCodelets9922()
            : darts::Codelet()
        {
        }
        _barrierCodelets9922(uint32_t dep, uint32_t res, TP9913* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets9931 : public darts::Codelet {
    public:
        TP9913* myTP;
        TP9913* inputsTPParent;
        _checkInCodelets9931()
            : darts::Codelet()
        {
        }
        _checkInCodelets9931(uint32_t dep, uint32_t res, TP9913* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets10028 : public darts::Codelet {
    public:
        TP9913* myTP;
        TP9913* inputsTPParent;
        _checkInCodelets10028()
            : darts::Codelet()
        {
        }
        _checkInCodelets10028(uint32_t dep, uint32_t res, TP9913* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets10038 : public darts::Codelet {
    public:
        TP9913* inputsTPParent;
        _barrierCodelets10038()
            : darts::Codelet()
        {
        }
        _barrierCodelets10038(uint32_t dep, uint32_t res, TP9913* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets10039 : public darts::Codelet {
    public:
        TP9913* myTP;
        TP9913* inputsTPParent;
        _checkInCodelets10039()
            : darts::Codelet()
        {
        }
        _checkInCodelets10039(uint32_t dep, uint32_t res, TP9913* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets10039 : public darts::Codelet {
    public:
        TP9913* inputsTPParent;
        _barrierCodelets10039()
            : darts::Codelet()
        {
        }
        _barrierCodelets10039(uint32_t dep, uint32_t res, TP9913* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP9913* TPParent;
    TP9913* controlTPParent;
    TP9913* inputsTPParent;
    int* iend_darts9913; /*OMP_SHARED - INPUT*/
    int* ist_darts9913; /*OMP_SHARED - INPUT*/
    int* jend_darts9913; /*OMP_SHARED - INPUT*/
    int* jst_darts9913; /*OMP_SHARED - INPUT*/
    int* nx0_darts9913; /*OMP_SHARED - INPUT*/
    int* ny0_darts9913; /*OMP_SHARED - INPUT*/
    int* nz0_darts9913; /*OMP_SHARED - INPUT*/
    double** sum_darts9913; /*OMP_SHARED - INPUT*/
    int* i_darts9913 /*VARIABLE*/;
    int* j_darts9913 /*VARIABLE*/;
    int* k_darts9913 /*VARIABLE*/;
    int* m_darts9913 /*VARIABLE*/;
    double* sum0_darts9913 /*VARIABLE*/;
    double* sum1_darts9913 /*VARIABLE*/;
    double* sum2_darts9913 /*VARIABLE*/;
    double* sum3_darts9913 /*VARIABLE*/;
    double* sum4_darts9913 /*VARIABLE*/;
    int m_darts10039;
    int* nx0_darts10039; /*OMP_SHARED - INPUT*/
    int* ny0_darts10039; /*OMP_SHARED - INPUT*/
    int* nz0_darts10039; /*OMP_SHARED - INPUT*/
    double** sum_darts10039; /*OMP_SHARED - INPUT*/
    int m_darts9922;
    double** sum_darts9922; /*OMP_SHARED - INPUT*/
    size_t TP9922_alreadyLaunched;
    TP9931** TP9931Ptr;
    size_t* TP9931_alreadyLaunched;
    int numTPsSet9931;
    int numTPsReady9931;
    size_t TPsToUse9931;
    size_t codeletsPerTP9931;
    size_t totalCodelets9931;
    size_t TP10039_alreadyLaunched;
    _barrierCodelets9913* barrierCodelets9913;
    _checkInCodelets9915* checkInCodelets9915;
    _checkInCodelets9922* checkInCodelets9922;
    _barrierCodelets9922* barrierCodelets9922;
    _checkInCodelets9931* checkInCodelets9931;
    _checkInCodelets10028* checkInCodelets10028;
    _barrierCodelets10038* barrierCodelets10038;
    _checkInCodelets10039* checkInCodelets10039;
    _barrierCodelets10039* barrierCodelets10039;
    TP9913(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet, int* in_iend,
        int* in_ist, int* in_jend, int* in_jst, int* in_nx0, int* in_ny0, int* in_nz0,
        double** in_sum);
    ~TP9913();
};
/*TP9931: OMPForDirective*/
class TP9931 : public ompTP {
public:
    bool requestNewRangeIterations9931(int* endRange, uint32_t codeletID);
    class _checkInCodelets9932 : public darts::Codelet {
    public:
        TP9931* myTP;
        TP9931* inputsTPParent;
        int endRange;
        _checkInCodelets9932()
            : darts::Codelet()
        {
        }
        _checkInCodelets9932(uint32_t dep, uint32_t res, TP9931* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP9913* TPParent;
    TP9931* controlTPParent;
    TP9931* inputsTPParent;
    int* i_darts9931 /*OMP_PRIVATE - INPUT*/;
    int* iend_darts9931; /*OMP_SHARED - INPUT*/
    int* ist_darts9931; /*OMP_SHARED - INPUT*/
    int* j_darts9931 /*OMP_PRIVATE - INPUT*/;
    int* jend_darts9931; /*OMP_SHARED - INPUT*/
    int* jst_darts9931; /*OMP_SHARED - INPUT*/
    int* k_darts9931 /*OMP_PRIVATE - INPUT*/;
    int* nz0_darts9931; /*OMP_SHARED - INPUT*/
    double** sum0_darts9931 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** sum1_darts9931 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** sum2_darts9931 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** sum3_darts9931 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** sum4_darts9931 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration9931;
    int lastIteration9931;
    int range9931;
    int rangePerCodelet9931;
    int minIteration9931;
    int remainderRange9931;
    size_t readyCodelets;
    int baseNumThreads;
    int* signalNextReady;
    _checkInCodelets9932* checkInCodelets9932;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets9932* firstCodelet;
#endif
    TP9931(int in_numThreads, int in_mainCodeletID, TP9913* in_TPParent, int in_initIteration,
        int in_lastIteration, int* in_iend, int* in_ist, int* in_jend, int* in_jst, int* in_nz0,
        TP9931** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP9931();
};
/*TP10830: OMPParallelDirective*/
class TP10830 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets10830 : public darts::Codelet {
    public:
        TP10830* inputsTPParent;
        _barrierCodelets10830()
            : darts::Codelet()
        {
        }
        _barrierCodelets10830(uint32_t dep, uint32_t res, TP10830* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets10845 : public darts::Codelet {
    public:
        TP10830* myTP;
        TP10830* inputsTPParent;
        _checkInCodelets10845()
            : darts::Codelet()
        {
        }
        _checkInCodelets10845(uint32_t dep, uint32_t res, TP10830* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets10845 : public darts::Codelet {
    public:
        TP10830* inputsTPParent;
        _barrierCodelets10845()
            : darts::Codelet()
        {
        }
        _barrierCodelets10845(uint32_t dep, uint32_t res, TP10830* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets10902 : public darts::Codelet {
    public:
        TP10830* myTP;
        TP10830* inputsTPParent;
        _checkInCodelets10902()
            : darts::Codelet()
        {
        }
        _checkInCodelets10902(uint32_t dep, uint32_t res, TP10830* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets10905 : public darts::Codelet {
    public:
        TP10830* myTP;
        TP10830* inputsTPParent;
        _checkInCodelets10905()
            : darts::Codelet()
        {
        }
        _checkInCodelets10905(uint32_t dep, uint32_t res, TP10830* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets10905 : public darts::Codelet {
    public:
        TP10830* inputsTPParent;
        _barrierCodelets10905()
            : darts::Codelet()
        {
        }
        _barrierCodelets10905(uint32_t dep, uint32_t res, TP10830* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets11054 : public darts::Codelet {
    public:
        TP10830* myTP;
        TP10830* inputsTPParent;
        _checkInCodelets11054()
            : darts::Codelet()
        {
        }
        _checkInCodelets11054(uint32_t dep, uint32_t res, TP10830* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets11054 : public darts::Codelet {
    public:
        TP10830* inputsTPParent;
        _barrierCodelets11054()
            : darts::Codelet()
        {
        }
        _barrierCodelets11054(uint32_t dep, uint32_t res, TP10830* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets11691 : public darts::Codelet {
    public:
        TP10830* myTP;
        TP10830* inputsTPParent;
        _checkInCodelets11691()
            : darts::Codelet()
        {
        }
        _checkInCodelets11691(uint32_t dep, uint32_t res, TP10830* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets11694 : public darts::Codelet {
    public:
        TP10830* myTP;
        TP10830* inputsTPParent;
        _checkInCodelets11694()
            : darts::Codelet()
        {
        }
        _checkInCodelets11694(uint32_t dep, uint32_t res, TP10830* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets11694 : public darts::Codelet {
    public:
        TP10830* inputsTPParent;
        _barrierCodelets11694()
            : darts::Codelet()
        {
        }
        _barrierCodelets11694(uint32_t dep, uint32_t res, TP10830* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets11843 : public darts::Codelet {
    public:
        TP10830* myTP;
        TP10830* inputsTPParent;
        _checkInCodelets11843()
            : darts::Codelet()
        {
        }
        _checkInCodelets11843(uint32_t dep, uint32_t res, TP10830* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets11843 : public darts::Codelet {
    public:
        TP10830* inputsTPParent;
        _barrierCodelets11843()
            : darts::Codelet()
        {
        }
        _barrierCodelets11843(uint32_t dep, uint32_t res, TP10830* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets12480 : public darts::Codelet {
    public:
        TP10830* myTP;
        TP10830* inputsTPParent;
        _checkInCodelets12480()
            : darts::Codelet()
        {
        }
        _checkInCodelets12480(uint32_t dep, uint32_t res, TP10830* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets12480 : public darts::Codelet {
    public:
        TP10830* inputsTPParent;
        _barrierCodelets12480()
            : darts::Codelet()
        {
        }
        _barrierCodelets12480(uint32_t dep, uint32_t res, TP10830* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP10830* TPParent;
    TP10830* controlTPParent;
    TP10830* inputsTPParent;
    int* L1_darts10830 /*VARIABLE*/;
    int* L2_darts10830 /*VARIABLE*/;
    int* i_darts10830 /*VARIABLE*/;
    int* iend1_darts10830 /*VARIABLE*/;
    int* ist1_darts10830 /*VARIABLE*/;
    int* j_darts10830 /*VARIABLE*/;
    int* jend1_darts10830 /*VARIABLE*/;
    int* jst1_darts10830 /*VARIABLE*/;
    int* k_darts10830 /*VARIABLE*/;
    int* m_darts10830 /*VARIABLE*/;
    double* q_darts10830 /*VARIABLE*/;
    double* tmp_darts10830 /*VARIABLE*/;
    double* u21_darts10830 /*VARIABLE*/;
    double* u21i_darts10830 /*VARIABLE*/;
    double* u21im1_darts10830 /*VARIABLE*/;
    double* u21j_darts10830 /*VARIABLE*/;
    double* u21jm1_darts10830 /*VARIABLE*/;
    double* u21k_darts10830 /*VARIABLE*/;
    double* u21km1_darts10830 /*VARIABLE*/;
    double* u31_darts10830 /*VARIABLE*/;
    double* u31i_darts10830 /*VARIABLE*/;
    double* u31im1_darts10830 /*VARIABLE*/;
    double* u31j_darts10830 /*VARIABLE*/;
    double* u31jm1_darts10830 /*VARIABLE*/;
    double* u31k_darts10830 /*VARIABLE*/;
    double* u31km1_darts10830 /*VARIABLE*/;
    double* u41_darts10830 /*VARIABLE*/;
    double* u41i_darts10830 /*VARIABLE*/;
    double* u41im1_darts10830 /*VARIABLE*/;
    double* u41j_darts10830 /*VARIABLE*/;
    double* u41jm1_darts10830 /*VARIABLE*/;
    double* u41k_darts10830 /*VARIABLE*/;
    double* u41km1_darts10830 /*VARIABLE*/;
    double* u51i_darts10830 /*VARIABLE*/;
    double* u51im1_darts10830 /*VARIABLE*/;
    double* u51j_darts10830 /*VARIABLE*/;
    double* u51jm1_darts10830 /*VARIABLE*/;
    double* u51k_darts10830 /*VARIABLE*/;
    double* u51km1_darts10830 /*VARIABLE*/;
    TP10845** TP10845Ptr;
    size_t* TP10845_alreadyLaunched;
    int numTPsSet10845;
    int numTPsReady10845;
    size_t TPsToUse10845;
    size_t codeletsPerTP10845;
    size_t totalCodelets10845;
    TP10905** TP10905Ptr;
    size_t* TP10905_alreadyLaunched;
    int numTPsSet10905;
    int numTPsReady10905;
    size_t TPsToUse10905;
    size_t codeletsPerTP10905;
    size_t totalCodelets10905;
    TP11054** TP11054Ptr;
    size_t* TP11054_alreadyLaunched;
    int numTPsSet11054;
    int numTPsReady11054;
    size_t TPsToUse11054;
    size_t codeletsPerTP11054;
    size_t totalCodelets11054;
    TP11694** TP11694Ptr;
    size_t* TP11694_alreadyLaunched;
    int numTPsSet11694;
    int numTPsReady11694;
    size_t TPsToUse11694;
    size_t codeletsPerTP11694;
    size_t totalCodelets11694;
    TP11843** TP11843Ptr;
    size_t* TP11843_alreadyLaunched;
    int numTPsSet11843;
    int numTPsReady11843;
    size_t TPsToUse11843;
    size_t codeletsPerTP11843;
    size_t totalCodelets11843;
    TP12480** TP12480Ptr;
    size_t* TP12480_alreadyLaunched;
    int numTPsSet12480;
    int numTPsReady12480;
    size_t TPsToUse12480;
    size_t codeletsPerTP12480;
    size_t totalCodelets12480;
    _barrierCodelets10830* barrierCodelets10830;
    _checkInCodelets10845* checkInCodelets10845;
    _barrierCodelets10845* barrierCodelets10845;
    _checkInCodelets10902* checkInCodelets10902;
    _checkInCodelets10905* checkInCodelets10905;
    _barrierCodelets10905* barrierCodelets10905;
    _checkInCodelets11054* checkInCodelets11054;
    _barrierCodelets11054* barrierCodelets11054;
    _checkInCodelets11691* checkInCodelets11691;
    _checkInCodelets11694* checkInCodelets11694;
    _barrierCodelets11694* barrierCodelets11694;
    _checkInCodelets11843* checkInCodelets11843;
    _barrierCodelets11843* barrierCodelets11843;
    _checkInCodelets12480* checkInCodelets12480;
    _barrierCodelets12480* barrierCodelets12480;
    TP10830(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
    ~TP10830();
};
/*TP10845: OMPForDirective*/
class TP10845 : public ompTP {
public:
    class _barrierCodelets10845 : public darts::Codelet {
    public:
        TP10845* inputsTPParent;
        _barrierCodelets10845()
            : darts::Codelet()
        {
        }
        _barrierCodelets10845(uint32_t dep, uint32_t res, TP10845* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations10845(int* endRange, uint32_t codeletID);
    class _checkInCodelets10846 : public darts::Codelet {
    public:
        TP10845* myTP;
        TP10845* inputsTPParent;
        int endRange;
        _checkInCodelets10846()
            : darts::Codelet()
        {
        }
        _checkInCodelets10846(uint32_t dep, uint32_t res, TP10845* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10830* TPParent;
    TP10845* controlTPParent;
    TP10845* inputsTPParent;
    int* i_darts10845 /*OMP_PRIVATE - INPUT*/;
    int* j_darts10845 /*OMP_PRIVATE - INPUT*/;
    int* k_darts10845 /*OMP_PRIVATE - INPUT*/;
    int* m_darts10845 /*OMP_PRIVATE - INPUT*/;
    int initIteration10845;
    int lastIteration10845;
    int range10845;
    int rangePerCodelet10845;
    int minIteration10845;
    int remainderRange10845;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets10845* barrierCodelets10845;
    _checkInCodelets10846* checkInCodelets10846;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets10846* firstCodelet;
#endif
    TP10845(int in_numThreads, int in_mainCodeletID, TP10830* in_TPParent, int in_initIteration,
        int in_lastIteration, TP10845** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP10845();
};
/*TP10905: OMPForDirective*/
class TP10905 : public ompTP {
public:
    class _barrierCodelets10905 : public darts::Codelet {
    public:
        TP10905* inputsTPParent;
        _barrierCodelets10905()
            : darts::Codelet()
        {
        }
        _barrierCodelets10905(uint32_t dep, uint32_t res, TP10905* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations10905(int* endRange, uint32_t codeletID);
    class _checkInCodelets10906 : public darts::Codelet {
    public:
        TP10905* myTP;
        TP10905* inputsTPParent;
        int endRange;
        _checkInCodelets10906()
            : darts::Codelet()
        {
        }
        _checkInCodelets10906(uint32_t dep, uint32_t res, TP10905* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10830* TPParent;
    TP10905* controlTPParent;
    TP10905* inputsTPParent;
    int** L1_darts10905 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** L2_darts10905 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts10905 /*OMP_PRIVATE - INPUT*/;
    int* j_darts10905 /*OMP_PRIVATE - INPUT*/;
    int* k_darts10905 /*OMP_PRIVATE - INPUT*/;
    double** q_darts10905 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21_darts10905 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration10905;
    int lastIteration10905;
    int range10905;
    int rangePerCodelet10905;
    int minIteration10905;
    int remainderRange10905;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets10905* barrierCodelets10905;
    _checkInCodelets10906* checkInCodelets10906;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets10906* firstCodelet;
#endif
    TP10905(int in_numThreads, int in_mainCodeletID, TP10830* in_TPParent, int in_initIteration,
        int in_lastIteration, TP10905** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP10905();
};
/*TP11054: OMPForDirective*/
class TP11054 : public ompTP {
public:
    class _barrierCodelets11054 : public darts::Codelet {
    public:
        TP11054* inputsTPParent;
        _barrierCodelets11054()
            : darts::Codelet()
        {
        }
        _barrierCodelets11054(uint32_t dep, uint32_t res, TP11054* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations11054(int* endRange, uint32_t codeletID);
    class _checkInCodelets11055 : public darts::Codelet {
    public:
        TP11054* myTP;
        TP11054* inputsTPParent;
        int endRange;
        _checkInCodelets11055()
            : darts::Codelet()
        {
        }
        _checkInCodelets11055(uint32_t dep, uint32_t res, TP11054* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10830* TPParent;
    TP11054* controlTPParent;
    TP11054* inputsTPParent;
    int** L2_darts11054 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts11054 /*OMP_PRIVATE - INPUT*/;
    int** iend1_darts11054 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** ist1_darts11054 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts11054 /*OMP_PRIVATE - INPUT*/;
    int* k_darts11054 /*OMP_PRIVATE - INPUT*/;
    int* m_darts11054 /*OMP_PRIVATE - INPUT*/;
    double** tmp_darts11054 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21i_darts11054 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21im1_darts11054 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31i_darts11054 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31im1_darts11054 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41i_darts11054 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41im1_darts11054 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51i_darts11054 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51im1_darts11054 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration11054;
    int lastIteration11054;
    int range11054;
    int rangePerCodelet11054;
    int minIteration11054;
    int remainderRange11054;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets11054* barrierCodelets11054;
    _checkInCodelets11055* checkInCodelets11055;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets11055* firstCodelet;
#endif
    TP11054(int in_numThreads, int in_mainCodeletID, TP10830* in_TPParent, int in_initIteration,
        int in_lastIteration, TP11054** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP11054();
};
/*TP11694: OMPForDirective*/
class TP11694 : public ompTP {
public:
    class _barrierCodelets11694 : public darts::Codelet {
    public:
        TP11694* inputsTPParent;
        _barrierCodelets11694()
            : darts::Codelet()
        {
        }
        _barrierCodelets11694(uint32_t dep, uint32_t res, TP11694* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations11694(int* endRange, uint32_t codeletID);
    class _checkInCodelets11695 : public darts::Codelet {
    public:
        TP11694* myTP;
        TP11694* inputsTPParent;
        int endRange;
        _checkInCodelets11695()
            : darts::Codelet()
        {
        }
        _checkInCodelets11695(uint32_t dep, uint32_t res, TP11694* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10830* TPParent;
    TP11694* controlTPParent;
    TP11694* inputsTPParent;
    int** L1_darts11694 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** L2_darts11694 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts11694 /*OMP_PRIVATE - INPUT*/;
    int* j_darts11694 /*OMP_PRIVATE - INPUT*/;
    int* k_darts11694 /*OMP_PRIVATE - INPUT*/;
    double** q_darts11694 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31_darts11694 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration11694;
    int lastIteration11694;
    int range11694;
    int rangePerCodelet11694;
    int minIteration11694;
    int remainderRange11694;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets11694* barrierCodelets11694;
    _checkInCodelets11695* checkInCodelets11695;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets11695* firstCodelet;
#endif
    TP11694(int in_numThreads, int in_mainCodeletID, TP10830* in_TPParent, int in_initIteration,
        int in_lastIteration, TP11694** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP11694();
};
/*TP11843: OMPForDirective*/
class TP11843 : public ompTP {
public:
    class _barrierCodelets11843 : public darts::Codelet {
    public:
        TP11843* inputsTPParent;
        _barrierCodelets11843()
            : darts::Codelet()
        {
        }
        _barrierCodelets11843(uint32_t dep, uint32_t res, TP11843* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations11843(int* endRange, uint32_t codeletID);
    class _checkInCodelets11844 : public darts::Codelet {
    public:
        TP11843* myTP;
        TP11843* inputsTPParent;
        int endRange;
        _checkInCodelets11844()
            : darts::Codelet()
        {
        }
        _checkInCodelets11844(uint32_t dep, uint32_t res, TP11843* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10830* TPParent;
    TP11843* controlTPParent;
    TP11843* inputsTPParent;
    int** L2_darts11843 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts11843 /*OMP_PRIVATE - INPUT*/;
    int* j_darts11843 /*OMP_PRIVATE - INPUT*/;
    int** jend1_darts11843 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** jst1_darts11843 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts11843 /*OMP_PRIVATE - INPUT*/;
    int* m_darts11843 /*OMP_PRIVATE - INPUT*/;
    double** tmp_darts11843 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21j_darts11843 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21jm1_darts11843 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31j_darts11843 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31jm1_darts11843 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41j_darts11843 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41jm1_darts11843 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51j_darts11843 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51jm1_darts11843 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration11843;
    int lastIteration11843;
    int range11843;
    int rangePerCodelet11843;
    int minIteration11843;
    int remainderRange11843;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets11843* barrierCodelets11843;
    _checkInCodelets11844* checkInCodelets11844;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets11844* firstCodelet;
#endif
    TP11843(int in_numThreads, int in_mainCodeletID, TP10830* in_TPParent, int in_initIteration,
        int in_lastIteration, TP11843** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP11843();
};
/*TP12480: OMPForDirective*/
class TP12480 : public ompTP {
public:
    class _barrierCodelets12480 : public darts::Codelet {
    public:
        TP12480* inputsTPParent;
        _barrierCodelets12480()
            : darts::Codelet()
        {
        }
        _barrierCodelets12480(uint32_t dep, uint32_t res, TP12480* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations12480(int* endRange, uint32_t codeletID);
    class _checkInCodelets12481 : public darts::Codelet {
    public:
        TP12480* myTP;
        TP12480* inputsTPParent;
        int endRange;
        _checkInCodelets12481()
            : darts::Codelet()
        {
        }
        _checkInCodelets12481(uint32_t dep, uint32_t res, TP12480* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10830* TPParent;
    TP12480* controlTPParent;
    TP12480* inputsTPParent;
    int* i_darts12480 /*OMP_PRIVATE - INPUT*/;
    int* j_darts12480 /*OMP_PRIVATE - INPUT*/;
    int* k_darts12480 /*OMP_PRIVATE - INPUT*/;
    int* m_darts12480 /*OMP_PRIVATE - INPUT*/;
    double** q_darts12480 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp_darts12480 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21k_darts12480 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21km1_darts12480 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31k_darts12480 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31km1_darts12480 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41_darts12480 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41k_darts12480 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41km1_darts12480 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51k_darts12480 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51km1_darts12480 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration12480;
    int lastIteration12480;
    int range12480;
    int rangePerCodelet12480;
    int minIteration12480;
    int remainderRange12480;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets12480* barrierCodelets12480;
    _checkInCodelets12481* checkInCodelets12481;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets12481* firstCodelet;
#endif
    TP12480(int in_numThreads, int in_mainCodeletID, TP10830* in_TPParent, int in_initIteration,
        int in_lastIteration, TP12480** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP12480();
};
/*TP13231: OMPParallelDirective*/
class TP13231 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets13231 : public darts::Codelet {
    public:
        TP13231* inputsTPParent;
        _barrierCodelets13231()
            : darts::Codelet()
        {
        }
        _barrierCodelets13231(uint32_t dep, uint32_t res, TP13231* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13235 : public darts::Codelet {
    public:
        TP13231* myTP;
        TP13231* inputsTPParent;
        _checkInCodelets13235()
            : darts::Codelet()
        {
        }
        _checkInCodelets13235(uint32_t dep, uint32_t res, TP13231* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13235 : public darts::Codelet {
    public:
        TP13231* inputsTPParent;
        _barrierCodelets13235()
            : darts::Codelet()
        {
        }
        _barrierCodelets13235(uint32_t dep, uint32_t res, TP13231* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13286 : public darts::Codelet {
    public:
        TP13231* myTP;
        TP13231* inputsTPParent;
        _checkInCodelets13286()
            : darts::Codelet()
        {
        }
        _checkInCodelets13286(uint32_t dep, uint32_t res, TP13231* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13286 : public darts::Codelet {
    public:
        TP13231* inputsTPParent;
        _barrierCodelets13286()
            : darts::Codelet()
        {
        }
        _barrierCodelets13286(uint32_t dep, uint32_t res, TP13231* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13328 : public darts::Codelet {
    public:
        TP13231* myTP;
        TP13231* inputsTPParent;
        _checkInCodelets13328()
            : darts::Codelet()
        {
        }
        _checkInCodelets13328(uint32_t dep, uint32_t res, TP13231* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13328 : public darts::Codelet {
    public:
        TP13231* inputsTPParent;
        _barrierCodelets13328()
            : darts::Codelet()
        {
        }
        _barrierCodelets13328(uint32_t dep, uint32_t res, TP13231* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13372 : public darts::Codelet {
    public:
        TP13231* myTP;
        TP13231* inputsTPParent;
        _checkInCodelets13372()
            : darts::Codelet()
        {
        }
        _checkInCodelets13372(uint32_t dep, uint32_t res, TP13231* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13372 : public darts::Codelet {
    public:
        TP13231* inputsTPParent;
        _barrierCodelets13372()
            : darts::Codelet()
        {
        }
        _barrierCodelets13372(uint32_t dep, uint32_t res, TP13231* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13414 : public darts::Codelet {
    public:
        TP13231* myTP;
        TP13231* inputsTPParent;
        _checkInCodelets13414()
            : darts::Codelet()
        {
        }
        _checkInCodelets13414(uint32_t dep, uint32_t res, TP13231* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13414 : public darts::Codelet {
    public:
        TP13231* inputsTPParent;
        _barrierCodelets13414()
            : darts::Codelet()
        {
        }
        _barrierCodelets13414(uint32_t dep, uint32_t res, TP13231* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP13231* TPParent;
    TP13231* controlTPParent;
    TP13231* inputsTPParent;
    int* i_darts13231 /*VARIABLE*/;
    int* iglob_darts13231 /*VARIABLE*/;
    int* j_darts13231 /*VARIABLE*/;
    int* jglob_darts13231 /*VARIABLE*/;
    int* k_darts13231 /*VARIABLE*/;
    TP13235** TP13235Ptr;
    size_t* TP13235_alreadyLaunched;
    int numTPsSet13235;
    int numTPsReady13235;
    size_t TPsToUse13235;
    size_t codeletsPerTP13235;
    size_t totalCodelets13235;
    TP13286** TP13286Ptr;
    size_t* TP13286_alreadyLaunched;
    int numTPsSet13286;
    int numTPsReady13286;
    size_t TPsToUse13286;
    size_t codeletsPerTP13286;
    size_t totalCodelets13286;
    TP13328** TP13328Ptr;
    size_t* TP13328_alreadyLaunched;
    int numTPsSet13328;
    int numTPsReady13328;
    size_t TPsToUse13328;
    size_t codeletsPerTP13328;
    size_t totalCodelets13328;
    TP13372** TP13372Ptr;
    size_t* TP13372_alreadyLaunched;
    int numTPsSet13372;
    int numTPsReady13372;
    size_t TPsToUse13372;
    size_t codeletsPerTP13372;
    size_t totalCodelets13372;
    TP13414** TP13414Ptr;
    size_t* TP13414_alreadyLaunched;
    int numTPsSet13414;
    int numTPsReady13414;
    size_t TPsToUse13414;
    size_t codeletsPerTP13414;
    size_t totalCodelets13414;
    _barrierCodelets13231* barrierCodelets13231;
    _checkInCodelets13235* checkInCodelets13235;
    _barrierCodelets13235* barrierCodelets13235;
    _checkInCodelets13286* checkInCodelets13286;
    _barrierCodelets13286* barrierCodelets13286;
    _checkInCodelets13328* checkInCodelets13328;
    _barrierCodelets13328* barrierCodelets13328;
    _checkInCodelets13372* checkInCodelets13372;
    _barrierCodelets13372* barrierCodelets13372;
    _checkInCodelets13414* checkInCodelets13414;
    _barrierCodelets13414* barrierCodelets13414;
    TP13231(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
    ~TP13231();
};
/*TP13235: OMPForDirective*/
class TP13235 : public ompTP {
public:
    class _barrierCodelets13235 : public darts::Codelet {
    public:
        TP13235* inputsTPParent;
        _barrierCodelets13235()
            : darts::Codelet()
        {
        }
        _barrierCodelets13235(uint32_t dep, uint32_t res, TP13235* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13235(int* endRange, uint32_t codeletID);
    class _checkInCodelets13236 : public darts::Codelet {
    public:
        TP13235* myTP;
        TP13235* inputsTPParent;
        int endRange;
        _checkInCodelets13236()
            : darts::Codelet()
        {
        }
        _checkInCodelets13236(uint32_t dep, uint32_t res, TP13235* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13231* TPParent;
    TP13235* controlTPParent;
    TP13235* inputsTPParent;
    int* i_darts13235 /*OMP_PRIVATE - INPUT*/;
    int** iglob_darts13235 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts13235 /*OMP_PRIVATE - INPUT*/;
    int** jglob_darts13235 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration13235;
    int lastIteration13235;
    int range13235;
    int rangePerCodelet13235;
    int minIteration13235;
    int remainderRange13235;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13235* barrierCodelets13235;
    _checkInCodelets13236* checkInCodelets13236;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13236* firstCodelet;
#endif
    TP13235(int in_numThreads, int in_mainCodeletID, TP13231* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13235** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13235();
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
    TP13231* TPParent;
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
    TP13286(int in_numThreads, int in_mainCodeletID, TP13231* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13286** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13286();
};
/*TP13328: OMPForDirective*/
class TP13328 : public ompTP {
public:
    class _barrierCodelets13328 : public darts::Codelet {
    public:
        TP13328* inputsTPParent;
        _barrierCodelets13328()
            : darts::Codelet()
        {
        }
        _barrierCodelets13328(uint32_t dep, uint32_t res, TP13328* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13328(int* endRange, uint32_t codeletID);
    class _checkInCodelets13329 : public darts::Codelet {
    public:
        TP13328* myTP;
        TP13328* inputsTPParent;
        int endRange;
        _checkInCodelets13329()
            : darts::Codelet()
        {
        }
        _checkInCodelets13329(uint32_t dep, uint32_t res, TP13328* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13231* TPParent;
    TP13328* controlTPParent;
    TP13328* inputsTPParent;
    int* i_darts13328 /*OMP_PRIVATE - INPUT*/;
    int** iglob_darts13328 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts13328 /*OMP_PRIVATE - INPUT*/;
    int initIteration13328;
    int lastIteration13328;
    int range13328;
    int rangePerCodelet13328;
    int minIteration13328;
    int remainderRange13328;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13328* barrierCodelets13328;
    _checkInCodelets13329* checkInCodelets13329;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13329* firstCodelet;
#endif
    TP13328(int in_numThreads, int in_mainCodeletID, TP13231* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13328** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13328();
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
    TP13231* TPParent;
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
    TP13372(int in_numThreads, int in_mainCodeletID, TP13231* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13372** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13372();
};
/*TP13414: OMPForDirective*/
class TP13414 : public ompTP {
public:
    class _barrierCodelets13414 : public darts::Codelet {
    public:
        TP13414* inputsTPParent;
        _barrierCodelets13414()
            : darts::Codelet()
        {
        }
        _barrierCodelets13414(uint32_t dep, uint32_t res, TP13414* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13414(int* endRange, uint32_t codeletID);
    class _checkInCodelets13415 : public darts::Codelet {
    public:
        TP13414* myTP;
        TP13414* inputsTPParent;
        int endRange;
        _checkInCodelets13415()
            : darts::Codelet()
        {
        }
        _checkInCodelets13415(uint32_t dep, uint32_t res, TP13414* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13231* TPParent;
    TP13414* controlTPParent;
    TP13414* inputsTPParent;
    int* j_darts13414 /*OMP_PRIVATE - INPUT*/;
    int** jglob_darts13414 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts13414 /*OMP_PRIVATE - INPUT*/;
    int initIteration13414;
    int lastIteration13414;
    int range13414;
    int rangePerCodelet13414;
    int minIteration13414;
    int remainderRange13414;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13414* barrierCodelets13414;
    _checkInCodelets13415* checkInCodelets13415;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13415* firstCodelet;
#endif
    TP13414(int in_numThreads, int in_mainCodeletID, TP13231* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13414** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13414();
};
/*TP13903: OMPParallelDirective*/
class TP13903 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets13903 : public darts::Codelet {
    public:
        TP13903* inputsTPParent;
        _barrierCodelets13903()
            : darts::Codelet()
        {
        }
        _barrierCodelets13903(uint32_t dep, uint32_t res, TP13903* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13905 : public darts::Codelet {
    public:
        TP13903* myTP;
        TP13903* inputsTPParent;
        _checkInCodelets13905()
            : darts::Codelet()
        {
        }
        _checkInCodelets13905(uint32_t dep, uint32_t res, TP13903* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13905 : public darts::Codelet {
    public:
        TP13903* inputsTPParent;
        _barrierCodelets13905()
            : darts::Codelet()
        {
        }
        _barrierCodelets13905(uint32_t dep, uint32_t res, TP13903* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP13903* TPParent;
    TP13903* controlTPParent;
    TP13903* inputsTPParent;
    int* i_darts13903 /*OMP_PRIVATE - INPUT*/;
    int* j_darts13903 /*OMP_PRIVATE - INPUT*/;
    int* k_darts13903 /*OMP_PRIVATE - INPUT*/;
    int* m_darts13903 /*OMP_PRIVATE - INPUT*/;
    TP13905** TP13905Ptr;
    size_t* TP13905_alreadyLaunched;
    int numTPsSet13905;
    int numTPsReady13905;
    size_t TPsToUse13905;
    size_t codeletsPerTP13905;
    size_t totalCodelets13905;
    _barrierCodelets13903* barrierCodelets13903;
    _checkInCodelets13905* checkInCodelets13905;
    _barrierCodelets13905* barrierCodelets13905;
    TP13903(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
    ~TP13903();
};
/*TP13905: OMPForDirective*/
class TP13905 : public ompTP {
public:
    class _barrierCodelets13905 : public darts::Codelet {
    public:
        TP13905* inputsTPParent;
        _barrierCodelets13905()
            : darts::Codelet()
        {
        }
        _barrierCodelets13905(uint32_t dep, uint32_t res, TP13905* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13905(int* endRange, uint32_t codeletID);
    class _checkInCodelets13906 : public darts::Codelet {
    public:
        TP13905* myTP;
        TP13905* inputsTPParent;
        int endRange;
        _checkInCodelets13906()
            : darts::Codelet()
        {
        }
        _checkInCodelets13906(uint32_t dep, uint32_t res, TP13905* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13903* TPParent;
    TP13905* controlTPParent;
    TP13905* inputsTPParent;
    int* i_darts13905 /*OMP_PRIVATE - INPUT*/;
    int* j_darts13905 /*OMP_PRIVATE - INPUT*/;
    int* k_darts13905 /*OMP_PRIVATE - INPUT*/;
    int* m_darts13905 /*OMP_PRIVATE - INPUT*/;
    int initIteration13905;
    int lastIteration13905;
    int range13905;
    int rangePerCodelet13905;
    int minIteration13905;
    int remainderRange13905;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13905* barrierCodelets13905;
    _checkInCodelets13906* checkInCodelets13906;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13906* firstCodelet;
#endif
    TP13905(int in_numThreads, int in_mainCodeletID, TP13903* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13905** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13905();
};
/*TP14002: OMPParallelDirective*/
class TP14002 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets14002 : public darts::Codelet {
    public:
        TP14002* inputsTPParent;
        _barrierCodelets14002()
            : darts::Codelet()
        {
        }
        _barrierCodelets14002(uint32_t dep, uint32_t res, TP14002* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14004 : public darts::Codelet {
    public:
        TP14002* myTP;
        TP14002* inputsTPParent;
        _checkInCodelets14004()
            : darts::Codelet()
        {
        }
        _checkInCodelets14004(uint32_t dep, uint32_t res, TP14002* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets14004 : public darts::Codelet {
    public:
        TP14002* inputsTPParent;
        _barrierCodelets14004()
            : darts::Codelet()
        {
        }
        _barrierCodelets14004(uint32_t dep, uint32_t res, TP14002* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14060 : public darts::Codelet {
    public:
        TP14002* myTP;
        TP14002* inputsTPParent;
        _checkInCodelets14060()
            : darts::Codelet()
        {
        }
        _checkInCodelets14060(uint32_t dep, uint32_t res, TP14002* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14059 : public darts::Codelet {
    public:
        TP14002* myTP;
        TP14002* inputsTPParent;
        _checkInCodelets14059()
            : darts::Codelet()
        {
        }
        _checkInCodelets14059(uint32_t dep, uint32_t res, TP14002* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14063 : public darts::Codelet {
    public:
        TP14002* myTP;
        TP14002* inputsTPParent;
        _checkInCodelets14063()
            : darts::Codelet()
        {
        }
        _checkInCodelets14063(uint32_t dep, uint32_t res, TP14002* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets14071 : public darts::Codelet {
    public:
        TP14002* inputsTPParent;
        _barrierCodelets14071()
            : darts::Codelet()
        {
        }
        _barrierCodelets14071(uint32_t dep, uint32_t res, TP14002* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14073 : public darts::Codelet {
    public:
        TP14002* myTP;
        TP14002* inputsTPParent;
        _checkInCodelets14073()
            : darts::Codelet()
        {
        }
        _checkInCodelets14073(uint32_t dep, uint32_t res, TP14002* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14072 : public darts::Codelet {
    public:
        TP14002* myTP;
        TP14002* inputsTPParent;
        _checkInCodelets14072()
            : darts::Codelet()
        {
        }
        _checkInCodelets14072(uint32_t dep, uint32_t res, TP14002* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14076 : public darts::Codelet {
    public:
        TP14002* myTP;
        TP14002* inputsTPParent;
        _checkInCodelets14076()
            : darts::Codelet()
        {
        }
        _checkInCodelets14076(uint32_t dep, uint32_t res, TP14002* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets14080 : public darts::Codelet {
    public:
        TP14002* inputsTPParent;
        _barrierCodelets14080()
            : darts::Codelet()
        {
        }
        _barrierCodelets14080(uint32_t dep, uint32_t res, TP14002* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14081 : public darts::Codelet {
    public:
        TP14002* myTP;
        TP14002* inputsTPParent;
        _checkInCodelets14081()
            : darts::Codelet()
        {
        }
        _checkInCodelets14081(uint32_t dep, uint32_t res, TP14002* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets14081 : public darts::Codelet {
    public:
        TP14002* inputsTPParent;
        _barrierCodelets14081()
            : darts::Codelet()
        {
        }
        _barrierCodelets14081(uint32_t dep, uint32_t res, TP14002* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP14002* TPParent;
    TP14002* controlTPParent;
    TP14002* inputsTPParent;
    int* i_darts14002 /*OMP_PRIVATE - INPUT*/;
    int* istep_darts14002 /*OMP_PRIVATE - INPUT*/;
    int* j_darts14002 /*OMP_PRIVATE - INPUT*/;
    int* k_darts14002 /*OMP_PRIVATE - INPUT*/;
    int* m_darts14002 /*OMP_PRIVATE - INPUT*/;
    double* tmp_darts14002; /*OMP_SHARED - INPUT*/
    int* k_darts14069 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts14066 /*OMP_SHARED_PRIVATE - INPUT*/;
    TP14004** TP14004Ptr;
    size_t* TP14004_alreadyLaunched;
    int numTPsSet14004;
    int numTPsReady14004;
    size_t TPsToUse14004;
    size_t codeletsPerTP14004;
    size_t totalCodelets14004;
    unsigned int TP14059_LoopCounter;
    unsigned int* TP14059_LoopCounterPerThread;
    tbb::concurrent_vector<TP14059*> TP14059PtrVec;
    unsigned int TP14072_LoopCounter;
    unsigned int* TP14072_LoopCounterPerThread;
    tbb::concurrent_vector<TP14072*> TP14072PtrVec;
    TP14081** TP14081Ptr;
    size_t* TP14081_alreadyLaunched;
    int numTPsSet14081;
    int numTPsReady14081;
    size_t TPsToUse14081;
    size_t codeletsPerTP14081;
    size_t totalCodelets14081;
    _barrierCodelets14002* barrierCodelets14002;
    _checkInCodelets14004* checkInCodelets14004;
    _barrierCodelets14004* barrierCodelets14004;
    _checkInCodelets14060* checkInCodelets14060;
    _checkInCodelets14059* checkInCodelets14059;
    _checkInCodelets14063* checkInCodelets14063;
    _barrierCodelets14071* barrierCodelets14071;
    _checkInCodelets14073* checkInCodelets14073;
    _checkInCodelets14072* checkInCodelets14072;
    _checkInCodelets14076* checkInCodelets14076;
    _barrierCodelets14080* barrierCodelets14080;
    _checkInCodelets14081* checkInCodelets14081;
    _barrierCodelets14081* barrierCodelets14081;
    TP14002(
        int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet, double* in_tmp);
    ~TP14002();
};
/*TP14004: OMPForDirective*/
class TP14004 : public ompTP {
public:
    class _barrierCodelets14004 : public darts::Codelet {
    public:
        TP14004* inputsTPParent;
        _barrierCodelets14004()
            : darts::Codelet()
        {
        }
        _barrierCodelets14004(uint32_t dep, uint32_t res, TP14004* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations14004(int* endRange, uint32_t codeletID);
    class _checkInCodelets14005 : public darts::Codelet {
    public:
        TP14004* myTP;
        TP14004* inputsTPParent;
        int endRange;
        _checkInCodelets14005()
            : darts::Codelet()
        {
        }
        _checkInCodelets14005(uint32_t dep, uint32_t res, TP14004* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP14002* TPParent;
    TP14004* controlTPParent;
    TP14004* inputsTPParent;
    int* i_darts14004 /*OMP_PRIVATE - INPUT*/;
    int* j_darts14004 /*OMP_PRIVATE - INPUT*/;
    int* k_darts14004 /*OMP_PRIVATE - INPUT*/;
    int* m_darts14004 /*OMP_PRIVATE - INPUT*/;
    int initIteration14004;
    int lastIteration14004;
    int range14004;
    int rangePerCodelet14004;
    int minIteration14004;
    int remainderRange14004;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets14004* barrierCodelets14004;
    _checkInCodelets14005* checkInCodelets14005;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets14005* firstCodelet;
#endif
    TP14004(int in_numThreads, int in_mainCodeletID, TP14002* in_TPParent, int in_initIteration,
        int in_lastIteration, TP14004** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP14004();
};
/*TP14059: ForStmt*/
class TP14059 : public ompTP {
public:
    class _checkInCodelets14065 : public darts::Codelet {
    public:
        TP14059* myTP;
        TP14002* inputsTPParent;
        _checkInCodelets14065()
            : darts::Codelet()
        {
        }
        _checkInCodelets14065(uint32_t dep, uint32_t res, TP14059* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14066 : public darts::Codelet {
    public:
        TP14059* myTP;
        TP14002* inputsTPParent;
        _checkInCodelets14066()
            : darts::Codelet()
        {
        }
        _checkInCodelets14066(uint32_t dep, uint32_t res, TP14059* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14068 : public darts::Codelet {
    public:
        TP14059* myTP;
        TP14002* inputsTPParent;
        _checkInCodelets14068()
            : darts::Codelet()
        {
        }
        _checkInCodelets14068(uint32_t dep, uint32_t res, TP14059* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14069 : public darts::Codelet {
    public:
        TP14059* myTP;
        TP14002* inputsTPParent;
        _checkInCodelets14069()
            : darts::Codelet()
        {
        }
        _checkInCodelets14069(uint32_t dep, uint32_t res, TP14059* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP14002* TPParent;
    TP14059* controlTPParent;
    TP14002* inputsTPParent;
    TP14059** ptrToThisTP;
    TP_jacld* TP14065Ptr;
    int TP14065_alreadyLaunched;
    size_t TP14066_alreadyLaunched;
    TP_blts* TP14068Ptr;
    int TP14068_alreadyLaunched;
    size_t TP14069_alreadyLaunched;
    _checkInCodelets14065* checkInCodelets14065;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets14065* firstCodelet;
#endif
    _checkInCodelets14066* checkInCodelets14066;
    _checkInCodelets14068* checkInCodelets14068;
    _checkInCodelets14069* checkInCodelets14069;
    TP14059(int in_numThreads, int in_mainCodeletID, TP14002* in_TPParent,
        TP14002* in_inputsTPParent, TP14059** in_ptrToThisTP);
    ~TP14059();
};
/*TP14072: ForStmt*/
class TP14072 : public ompTP {
public:
    class _checkInCodelets14078 : public darts::Codelet {
    public:
        TP14072* myTP;
        TP14002* inputsTPParent;
        _checkInCodelets14078()
            : darts::Codelet()
        {
        }
        _checkInCodelets14078(uint32_t dep, uint32_t res, TP14072* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14079 : public darts::Codelet {
    public:
        TP14072* myTP;
        TP14002* inputsTPParent;
        _checkInCodelets14079()
            : darts::Codelet()
        {
        }
        _checkInCodelets14079(uint32_t dep, uint32_t res, TP14072* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP14002* TPParent;
    TP14072* controlTPParent;
    TP14002* inputsTPParent;
    TP14072** ptrToThisTP;
    TP_jacu* TP14078Ptr;
    int TP14078_alreadyLaunched;
    TP_buts* TP14079Ptr;
    int TP14079_alreadyLaunched;
    _checkInCodelets14078* checkInCodelets14078;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets14078* firstCodelet;
#endif
    _checkInCodelets14079* checkInCodelets14079;
    TP14072(int in_numThreads, int in_mainCodeletID, TP14002* in_TPParent,
        TP14002* in_inputsTPParent, TP14072** in_ptrToThisTP);
    ~TP14072();
};
/*TP14081: OMPForDirective*/
class TP14081 : public ompTP {
public:
    class _barrierCodelets14081 : public darts::Codelet {
    public:
        TP14081* inputsTPParent;
        _barrierCodelets14081()
            : darts::Codelet()
        {
        }
        _barrierCodelets14081(uint32_t dep, uint32_t res, TP14081* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations14081(int* endRange, uint32_t codeletID);
    class _checkInCodelets14082 : public darts::Codelet {
    public:
        TP14081* myTP;
        TP14081* inputsTPParent;
        int endRange;
        _checkInCodelets14082()
            : darts::Codelet()
        {
        }
        _checkInCodelets14082(uint32_t dep, uint32_t res, TP14081* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP14002* TPParent;
    TP14081* controlTPParent;
    TP14081* inputsTPParent;
    int* i_darts14081 /*OMP_PRIVATE - INPUT*/;
    int* j_darts14081 /*OMP_PRIVATE - INPUT*/;
    int* k_darts14081 /*OMP_PRIVATE - INPUT*/;
    int* m_darts14081 /*OMP_PRIVATE - INPUT*/;
    double* tmp_darts14081; /*OMP_SHARED - INPUT*/
    int initIteration14081;
    int lastIteration14081;
    int range14081;
    int rangePerCodelet14081;
    int minIteration14081;
    int remainderRange14081;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets14081* barrierCodelets14081;
    _checkInCodelets14082* checkInCodelets14082;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets14082* firstCodelet;
#endif
    TP14081(int in_numThreads, int in_mainCodeletID, TP14002* in_TPParent, int in_initIteration,
        int in_lastIteration, double* in_tmp, TP14081** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP14081();
};
#endif
