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
class TP1;
typedef TP1 TP_blts;
class TP2;
typedef TP2 TP_buts;
class TP234;
class TP1;
typedef TP1 TP_blts;
class TP2;
typedef TP2 TP_buts;
class TP2301;
class TP2320;
/*Number of TPs to be used for the OMPFor in region TP2320*/
#define NUMTPS2320 NUMTPS
class TP2371;
/*Number of TPs to be used for the OMPFor in region TP2371*/
#define NUMTPS2371 NUMTPS
class TP2506;
/*Number of TPs to be used for the OMPFor in region TP2506*/
#define NUMTPS2506 NUMTPS
class TP2655;
/*Number of TPs to be used for the OMPFor in region TP2655*/
#define NUMTPS2655 NUMTPS
class TP3293;
/*Number of TPs to be used for the OMPFor in region TP3293*/
#define NUMTPS3293 NUMTPS
class TP3442;
/*Number of TPs to be used for the OMPFor in region TP3442*/
#define NUMTPS3442 NUMTPS
class TP4077;
/*Number of TPs to be used for the OMPFor in region TP4077*/
#define NUMTPS4077 NUMTPS
class TP7;
typedef TP7 TP_jacld;
class TP5003;
/*Number of TPs to be used for the OMPFor in region TP5003*/
#define NUMTPS5003 NUMTPS
class TP8;
typedef TP8 TP_jacu;
class TP7522;
/*Number of TPs to be used for the OMPFor in region TP7522*/
#define NUMTPS7522 NUMTPS
class TP10896;
class TP10911;
/*Number of TPs to be used for the OMPFor in region TP10911*/
#define NUMTPS10911 NUMTPS
class TP10971;
/*Number of TPs to be used for the OMPFor in region TP10971*/
#define NUMTPS10971 NUMTPS
class TP11120;
/*Number of TPs to be used for the OMPFor in region TP11120*/
#define NUMTPS11120 NUMTPS
class TP11760;
/*Number of TPs to be used for the OMPFor in region TP11760*/
#define NUMTPS11760 NUMTPS
class TP11909;
/*Number of TPs to be used for the OMPFor in region TP11909*/
#define NUMTPS11909 NUMTPS
class TP12546;
/*Number of TPs to be used for the OMPFor in region TP12546*/
#define NUMTPS12546 NUMTPS
class TP13297;
class TP13301;
/*Number of TPs to be used for the OMPFor in region TP13301*/
#define NUMTPS13301 NUMTPS
class TP13352;
/*Number of TPs to be used for the OMPFor in region TP13352*/
#define NUMTPS13352 NUMTPS
class TP13394;
/*Number of TPs to be used for the OMPFor in region TP13394*/
#define NUMTPS13394 NUMTPS
class TP13438;
/*Number of TPs to be used for the OMPFor in region TP13438*/
#define NUMTPS13438 NUMTPS
class TP13480;
/*Number of TPs to be used for the OMPFor in region TP13480*/
#define NUMTPS13480 NUMTPS
class TP13864;
class TP13871;
/*Number of TPs to be used for the OMPFor in region TP13871*/
#define NUMTPS13871 NUMTPS
class TP13997;
class TP13999;
/*Number of TPs to be used for the OMPFor in region TP13999*/
#define NUMTPS13999 NUMTPS
class TP14096;
class TP14098;
/*Number of TPs to be used for the OMPFor in region TP14098*/
#define NUMTPS14098 NUMTPS
class TP14153;
class TP14162;
class TP14171;
/*Number of TPs to be used for the OMPFor in region TP14171*/
#define NUMTPS14171 NUMTPS
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
/*TP1: blts*/
class TP1 : public ompTP {
public:
    class _checkInCodelets205 : public darts::Codelet {
    public:
        TP1* myTP;
        TP1* inputsTPParent;
        _checkInCodelets205()
            : darts::Codelet()
        {
        }
        _checkInCodelets205(uint32_t dep, uint32_t res, TP1* myTP, uint32_t id)
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
    double(*) * *v_darts1 /*VARIABLE*/;
    int v_outer1_size;
    double(*) * *ldz_darts1 /*VARIABLE*/;
    int ldz_outer1_size;
    double(*) * *ldy_darts1 /*VARIABLE*/;
    int ldy_outer1_size;
    double(*) * *ldx_darts1 /*VARIABLE*/;
    int ldx_outer1_size;
    double(*) * *d_darts1 /*VARIABLE*/;
    int d_outer1_size;
    int* ist_darts1 /*VARIABLE*/;
    int* iend_darts1 /*VARIABLE*/;
    int* jst_darts1 /*VARIABLE*/;
    int* jend_darts1 /*VARIABLE*/;
    int* nx0_darts1 /*VARIABLE*/;
    int* ny0_darts1 /*VARIABLE*/;
    _checkInCodelets205* checkInCodelets205;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets205* firstCodelet;
#endif
    TP1(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
        darts::Codelet* in_mainSyncCodelet, TP1** in_ptrToThisFunctionTP, int in_nx, int in_ny,
        int in_nz, int in_k, double in_omega, double(*) in_v, double(*) in_ldz, double(*) in_ldy,
        double(*) in_ldx, double(*) in_d, int in_ist, int in_iend, int in_jst, int in_jend,
        int in_nx0, int in_ny0);
    ~TP1();
    void setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, double(*) in_v,
        double(*) in_ldz, double(*) in_ldy, double(*) in_ldx, double(*) in_d, int in_ist,
        int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID);
};
/*TP2: buts*/
class TP2 : public ompTP {
public:
    class _checkInCodelets211 : public darts::Codelet {
    public:
        TP2* myTP;
        TP2* inputsTPParent;
        _checkInCodelets211()
            : darts::Codelet()
        {
        }
        _checkInCodelets211(uint32_t dep, uint32_t res, TP2* myTP, uint32_t id)
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
    double(*) * *v_darts2 /*VARIABLE*/;
    int v_outer2_size;
    double(*) * *tv_darts2 /*VARIABLE*/;
    int tv_outer2_size;
    double(*) * *d_darts2 /*VARIABLE*/;
    int d_outer2_size;
    double(*) * *udx_darts2 /*VARIABLE*/;
    int udx_outer2_size;
    double(*) * *udy_darts2 /*VARIABLE*/;
    int udy_outer2_size;
    double(*) * *udz_darts2 /*VARIABLE*/;
    int udz_outer2_size;
    int* ist_darts2 /*VARIABLE*/;
    int* iend_darts2 /*VARIABLE*/;
    int* jst_darts2 /*VARIABLE*/;
    int* jend_darts2 /*VARIABLE*/;
    int* nx0_darts2 /*VARIABLE*/;
    int* ny0_darts2 /*VARIABLE*/;
    _checkInCodelets211* checkInCodelets211;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets211* firstCodelet;
#endif
    TP2(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
        darts::Codelet* in_mainSyncCodelet, TP2** in_ptrToThisFunctionTP, int in_nx, int in_ny,
        int in_nz, int in_k, double in_omega, double(*) in_v, double(*) in_tv, double(*) in_d,
        double(*) in_udx, double(*) in_udy, double(*) in_udz, int in_ist, int in_iend, int in_jst,
        int in_jend, int in_nx0, int in_ny0);
    ~TP2();
    void setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, double(*) in_v,
        double(*) in_tv, double(*) in_d, double(*) in_udx, double(*) in_udy, double(*) in_udz,
        int in_ist, int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID);
};
/*TP234: OMPParallelDirective*/
class TP234 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets234 : public darts::Codelet {
    public:
        TP234* inputsTPParent;
        _barrierCodelets234()
            : darts::Codelet()
        {
        }
        _barrierCodelets234(uint32_t dep, uint32_t res, TP234* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets236 : public darts::Codelet {
    public:
        TP234* myTP;
        TP234* inputsTPParent;
        _checkInCodelets236()
            : darts::Codelet()
        {
        }
        _checkInCodelets236(uint32_t dep, uint32_t res, TP234* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP234* TPParent;
    TP234* controlTPParent;
    TP234* inputsTPParent;
    int* nthreads_darts234; /*OMP_SHARED - INPUT*/
    int* nthreads_darts236; /*OMP_SHARED - INPUT*/
    size_t TP236_alreadyLaunched;
    _barrierCodelets234* barrierCodelets234;
    _checkInCodelets236* checkInCodelets236;
    TP234(
        int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet, int* in_nthreads);
    ~TP234();
};
/*TP1: blts*/
class TP1 : public ompTP {
public:
    class _checkInCodelets205 : public darts::Codelet {
    public:
        TP1* myTP;
        TP1* inputsTPParent;
        _checkInCodelets205()
            : darts::Codelet()
        {
        }
        _checkInCodelets205(uint32_t dep, uint32_t res, TP1* myTP, uint32_t id)
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
    double(*) * *v_darts1 /*VARIABLE*/;
    int v_outer1_size;
    double(*) * *ldz_darts1 /*VARIABLE*/;
    int ldz_outer1_size;
    double(*) * *ldy_darts1 /*VARIABLE*/;
    int ldy_outer1_size;
    double(*) * *ldx_darts1 /*VARIABLE*/;
    int ldx_outer1_size;
    double(*) * *d_darts1 /*VARIABLE*/;
    int d_outer1_size;
    int* ist_darts1 /*VARIABLE*/;
    int* iend_darts1 /*VARIABLE*/;
    int* jst_darts1 /*VARIABLE*/;
    int* jend_darts1 /*VARIABLE*/;
    int* nx0_darts1 /*VARIABLE*/;
    int* ny0_darts1 /*VARIABLE*/;
    _checkInCodelets205* checkInCodelets205;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets205* firstCodelet;
#endif
    TP1(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
        darts::Codelet* in_mainSyncCodelet, TP1** in_ptrToThisFunctionTP, int in_nx, int in_ny,
        int in_nz, int in_k, double in_omega, double(*) in_v, double(*) in_ldz, double(*) in_ldy,
        double(*) in_ldx, double(*) in_d, int in_ist, int in_iend, int in_jst, int in_jend,
        int in_nx0, int in_ny0);
    ~TP1();
    void setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, double(*) in_v,
        double(*) in_ldz, double(*) in_ldy, double(*) in_ldx, double(*) in_d, int in_ist,
        int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID);
};
/*TP2: buts*/
class TP2 : public ompTP {
public:
    class _checkInCodelets211 : public darts::Codelet {
    public:
        TP2* myTP;
        TP2* inputsTPParent;
        _checkInCodelets211()
            : darts::Codelet()
        {
        }
        _checkInCodelets211(uint32_t dep, uint32_t res, TP2* myTP, uint32_t id)
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
    double(*) * *v_darts2 /*VARIABLE*/;
    int v_outer2_size;
    double(*) * *tv_darts2 /*VARIABLE*/;
    int tv_outer2_size;
    double(*) * *d_darts2 /*VARIABLE*/;
    int d_outer2_size;
    double(*) * *udx_darts2 /*VARIABLE*/;
    int udx_outer2_size;
    double(*) * *udy_darts2 /*VARIABLE*/;
    int udy_outer2_size;
    double(*) * *udz_darts2 /*VARIABLE*/;
    int udz_outer2_size;
    int* ist_darts2 /*VARIABLE*/;
    int* iend_darts2 /*VARIABLE*/;
    int* jst_darts2 /*VARIABLE*/;
    int* jend_darts2 /*VARIABLE*/;
    int* nx0_darts2 /*VARIABLE*/;
    int* ny0_darts2 /*VARIABLE*/;
    _checkInCodelets211* checkInCodelets211;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets211* firstCodelet;
#endif
    TP2(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
        darts::Codelet* in_mainSyncCodelet, TP2** in_ptrToThisFunctionTP, int in_nx, int in_ny,
        int in_nz, int in_k, double in_omega, double(*) in_v, double(*) in_tv, double(*) in_d,
        double(*) in_udx, double(*) in_udy, double(*) in_udz, int in_ist, int in_iend, int in_jst,
        int in_jend, int in_nx0, int in_ny0);
    ~TP2();
    void setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, double(*) in_v,
        double(*) in_tv, double(*) in_d, double(*) in_udx, double(*) in_udy, double(*) in_udz,
        int in_ist, int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID);
};
/*TP2301: OMPParallelDirective*/
class TP2301 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets2301 : public darts::Codelet {
    public:
        TP2301* inputsTPParent;
        _barrierCodelets2301()
            : darts::Codelet()
        {
        }
        _barrierCodelets2301(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2303 : public darts::Codelet {
    public:
        TP2301* myTP;
        TP2301* inputsTPParent;
        _checkInCodelets2303()
            : darts::Codelet()
        {
        }
        _checkInCodelets2303(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2320 : public darts::Codelet {
    public:
        TP2301* myTP;
        TP2301* inputsTPParent;
        _checkInCodelets2320()
            : darts::Codelet()
        {
        }
        _checkInCodelets2320(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets2320 : public darts::Codelet {
    public:
        TP2301* inputsTPParent;
        _barrierCodelets2320()
            : darts::Codelet()
        {
        }
        _barrierCodelets2320(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2371 : public darts::Codelet {
    public:
        TP2301* myTP;
        TP2301* inputsTPParent;
        _checkInCodelets2371()
            : darts::Codelet()
        {
        }
        _checkInCodelets2371(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets2371 : public darts::Codelet {
    public:
        TP2301* inputsTPParent;
        _barrierCodelets2371()
            : darts::Codelet()
        {
        }
        _barrierCodelets2371(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2503 : public darts::Codelet {
    public:
        TP2301* myTP;
        TP2301* inputsTPParent;
        _checkInCodelets2503()
            : darts::Codelet()
        {
        }
        _checkInCodelets2503(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2506 : public darts::Codelet {
    public:
        TP2301* myTP;
        TP2301* inputsTPParent;
        _checkInCodelets2506()
            : darts::Codelet()
        {
        }
        _checkInCodelets2506(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets2506 : public darts::Codelet {
    public:
        TP2301* inputsTPParent;
        _barrierCodelets2506()
            : darts::Codelet()
        {
        }
        _barrierCodelets2506(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets2655 : public darts::Codelet {
    public:
        TP2301* myTP;
        TP2301* inputsTPParent;
        _checkInCodelets2655()
            : darts::Codelet()
        {
        }
        _checkInCodelets2655(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets2655 : public darts::Codelet {
    public:
        TP2301* inputsTPParent;
        _barrierCodelets2655()
            : darts::Codelet()
        {
        }
        _barrierCodelets2655(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets3290 : public darts::Codelet {
    public:
        TP2301* myTP;
        TP2301* inputsTPParent;
        _checkInCodelets3290()
            : darts::Codelet()
        {
        }
        _checkInCodelets3290(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets3293 : public darts::Codelet {
    public:
        TP2301* myTP;
        TP2301* inputsTPParent;
        _checkInCodelets3293()
            : darts::Codelet()
        {
        }
        _checkInCodelets3293(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets3293 : public darts::Codelet {
    public:
        TP2301* inputsTPParent;
        _barrierCodelets3293()
            : darts::Codelet()
        {
        }
        _barrierCodelets3293(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets3442 : public darts::Codelet {
    public:
        TP2301* myTP;
        TP2301* inputsTPParent;
        _checkInCodelets3442()
            : darts::Codelet()
        {
        }
        _checkInCodelets3442(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets3442 : public darts::Codelet {
    public:
        TP2301* inputsTPParent;
        _barrierCodelets3442()
            : darts::Codelet()
        {
        }
        _barrierCodelets3442(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets4077 : public darts::Codelet {
    public:
        TP2301* myTP;
        TP2301* inputsTPParent;
        _checkInCodelets4077()
            : darts::Codelet()
        {
        }
        _checkInCodelets4077(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets4077 : public darts::Codelet {
    public:
        TP2301* inputsTPParent;
        _barrierCodelets4077()
            : darts::Codelet()
        {
        }
        _barrierCodelets4077(uint32_t dep, uint32_t res, TP2301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP2301* TPParent;
    TP2301* controlTPParent;
    TP2301* inputsTPParent;
    int* L1_darts2301 /*VARIABLE*/;
    int* L2_darts2301 /*VARIABLE*/;
    double* dsspm_darts2301 /*VARIABLE*/;
    double* eta_darts2301 /*VARIABLE*/;
    int* i_darts2301 /*VARIABLE*/;
    int* iend1_darts2301 /*VARIABLE*/;
    int* iglob_darts2301 /*VARIABLE*/;
    int* ist1_darts2301 /*VARIABLE*/;
    int* j_darts2301 /*VARIABLE*/;
    int* jend1_darts2301 /*VARIABLE*/;
    int* jglob_darts2301 /*VARIABLE*/;
    int* jst1_darts2301 /*VARIABLE*/;
    int* k_darts2301 /*VARIABLE*/;
    int* m_darts2301 /*VARIABLE*/;
    double* q_darts2301 /*VARIABLE*/;
    double* tmp_darts2301 /*VARIABLE*/;
    double* u21_darts2301 /*VARIABLE*/;
    double* u21i_darts2301 /*VARIABLE*/;
    double* u21im1_darts2301 /*VARIABLE*/;
    double* u21j_darts2301 /*VARIABLE*/;
    double* u21jm1_darts2301 /*VARIABLE*/;
    double* u21k_darts2301 /*VARIABLE*/;
    double* u21km1_darts2301 /*VARIABLE*/;
    double* u31_darts2301 /*VARIABLE*/;
    double* u31i_darts2301 /*VARIABLE*/;
    double* u31im1_darts2301 /*VARIABLE*/;
    double* u31j_darts2301 /*VARIABLE*/;
    double* u31jm1_darts2301 /*VARIABLE*/;
    double* u31k_darts2301 /*VARIABLE*/;
    double* u31km1_darts2301 /*VARIABLE*/;
    double* u41_darts2301 /*VARIABLE*/;
    double* u41i_darts2301 /*VARIABLE*/;
    double* u41im1_darts2301 /*VARIABLE*/;
    double* u41j_darts2301 /*VARIABLE*/;
    double* u41jm1_darts2301 /*VARIABLE*/;
    double* u41k_darts2301 /*VARIABLE*/;
    double* u41km1_darts2301 /*VARIABLE*/;
    double* u51i_darts2301 /*VARIABLE*/;
    double* u51im1_darts2301 /*VARIABLE*/;
    double* u51j_darts2301 /*VARIABLE*/;
    double* u51jm1_darts2301 /*VARIABLE*/;
    double* u51k_darts2301 /*VARIABLE*/;
    double* u51km1_darts2301 /*VARIABLE*/;
    double* xi_darts2301 /*VARIABLE*/;
    double* zeta_darts2301 /*VARIABLE*/;
    TP2320** TP2320Ptr;
    size_t* TP2320_alreadyLaunched;
    int numTPsSet2320;
    int numTPsReady2320;
    size_t TPsToUse2320;
    size_t codeletsPerTP2320;
    size_t totalCodelets2320;
    TP2371** TP2371Ptr;
    size_t* TP2371_alreadyLaunched;
    int numTPsSet2371;
    int numTPsReady2371;
    size_t TPsToUse2371;
    size_t codeletsPerTP2371;
    size_t totalCodelets2371;
    TP2506** TP2506Ptr;
    size_t* TP2506_alreadyLaunched;
    int numTPsSet2506;
    int numTPsReady2506;
    size_t TPsToUse2506;
    size_t codeletsPerTP2506;
    size_t totalCodelets2506;
    TP2655** TP2655Ptr;
    size_t* TP2655_alreadyLaunched;
    int numTPsSet2655;
    int numTPsReady2655;
    size_t TPsToUse2655;
    size_t codeletsPerTP2655;
    size_t totalCodelets2655;
    TP3293** TP3293Ptr;
    size_t* TP3293_alreadyLaunched;
    int numTPsSet3293;
    int numTPsReady3293;
    size_t TPsToUse3293;
    size_t codeletsPerTP3293;
    size_t totalCodelets3293;
    TP3442** TP3442Ptr;
    size_t* TP3442_alreadyLaunched;
    int numTPsSet3442;
    int numTPsReady3442;
    size_t TPsToUse3442;
    size_t codeletsPerTP3442;
    size_t totalCodelets3442;
    TP4077** TP4077Ptr;
    size_t* TP4077_alreadyLaunched;
    int numTPsSet4077;
    int numTPsReady4077;
    size_t TPsToUse4077;
    size_t codeletsPerTP4077;
    size_t totalCodelets4077;
    _barrierCodelets2301* barrierCodelets2301;
    _checkInCodelets2303* checkInCodelets2303;
    _checkInCodelets2320* checkInCodelets2320;
    _barrierCodelets2320* barrierCodelets2320;
    _checkInCodelets2371* checkInCodelets2371;
    _barrierCodelets2371* barrierCodelets2371;
    _checkInCodelets2503* checkInCodelets2503;
    _checkInCodelets2506* checkInCodelets2506;
    _barrierCodelets2506* barrierCodelets2506;
    _checkInCodelets2655* checkInCodelets2655;
    _barrierCodelets2655* barrierCodelets2655;
    _checkInCodelets3290* checkInCodelets3290;
    _checkInCodelets3293* checkInCodelets3293;
    _barrierCodelets3293* barrierCodelets3293;
    _checkInCodelets3442* checkInCodelets3442;
    _barrierCodelets3442* barrierCodelets3442;
    _checkInCodelets4077* checkInCodelets4077;
    _barrierCodelets4077* barrierCodelets4077;
    TP2301(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
    ~TP2301();
};
/*TP2320: OMPForDirective*/
class TP2320 : public ompTP {
public:
    class _barrierCodelets2320 : public darts::Codelet {
    public:
        TP2320* inputsTPParent;
        _barrierCodelets2320()
            : darts::Codelet()
        {
        }
        _barrierCodelets2320(uint32_t dep, uint32_t res, TP2320* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations2320(int* endRange, uint32_t codeletID);
    class _checkInCodelets2321 : public darts::Codelet {
    public:
        TP2320* myTP;
        TP2320* inputsTPParent;
        int endRange;
        _checkInCodelets2321()
            : darts::Codelet()
        {
        }
        _checkInCodelets2321(uint32_t dep, uint32_t res, TP2320* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2301* TPParent;
    TP2320* controlTPParent;
    TP2320* inputsTPParent;
    int* i_darts2320 /*OMP_PRIVATE - INPUT*/;
    int* j_darts2320 /*OMP_PRIVATE - INPUT*/;
    int* k_darts2320 /*OMP_PRIVATE - INPUT*/;
    int* m_darts2320 /*OMP_PRIVATE - INPUT*/;
    int initIteration2320;
    int lastIteration2320;
    int range2320;
    int rangePerCodelet2320;
    int minIteration2320;
    int remainderRange2320;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets2320* barrierCodelets2320;
    _checkInCodelets2321* checkInCodelets2321;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets2321* firstCodelet;
#endif
    TP2320(int in_numThreads, int in_mainCodeletID, TP2301* in_TPParent, int in_initIteration,
        int in_lastIteration, TP2320** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP2320();
};
/*TP2371: OMPForDirective*/
class TP2371 : public ompTP {
public:
    class _barrierCodelets2371 : public darts::Codelet {
    public:
        TP2371* inputsTPParent;
        _barrierCodelets2371()
            : darts::Codelet()
        {
        }
        _barrierCodelets2371(uint32_t dep, uint32_t res, TP2371* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations2371(int* endRange, uint32_t codeletID);
    class _checkInCodelets2372 : public darts::Codelet {
    public:
        TP2371* myTP;
        TP2371* inputsTPParent;
        int endRange;
        _checkInCodelets2372()
            : darts::Codelet()
        {
        }
        _checkInCodelets2372(uint32_t dep, uint32_t res, TP2371* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2301* TPParent;
    TP2371* controlTPParent;
    TP2371* inputsTPParent;
    double** eta_darts2371 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts2371 /*OMP_PRIVATE - INPUT*/;
    int** iglob_darts2371 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts2371 /*OMP_PRIVATE - INPUT*/;
    int** jglob_darts2371 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts2371 /*OMP_PRIVATE - INPUT*/;
    int* m_darts2371 /*OMP_PRIVATE - INPUT*/;
    double** xi_darts2371 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** zeta_darts2371 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration2371;
    int lastIteration2371;
    int range2371;
    int rangePerCodelet2371;
    int minIteration2371;
    int remainderRange2371;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets2371* barrierCodelets2371;
    _checkInCodelets2372* checkInCodelets2372;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets2372* firstCodelet;
#endif
    TP2371(int in_numThreads, int in_mainCodeletID, TP2301* in_TPParent, int in_initIteration,
        int in_lastIteration, TP2371** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP2371();
};
/*TP2506: OMPForDirective*/
class TP2506 : public ompTP {
public:
    class _barrierCodelets2506 : public darts::Codelet {
    public:
        TP2506* inputsTPParent;
        _barrierCodelets2506()
            : darts::Codelet()
        {
        }
        _barrierCodelets2506(uint32_t dep, uint32_t res, TP2506* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations2506(int* endRange, uint32_t codeletID);
    class _checkInCodelets2507 : public darts::Codelet {
    public:
        TP2506* myTP;
        TP2506* inputsTPParent;
        int endRange;
        _checkInCodelets2507()
            : darts::Codelet()
        {
        }
        _checkInCodelets2507(uint32_t dep, uint32_t res, TP2506* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2301* TPParent;
    TP2506* controlTPParent;
    TP2506* inputsTPParent;
    int** L1_darts2506 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** L2_darts2506 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts2506 /*OMP_PRIVATE - INPUT*/;
    int* j_darts2506 /*OMP_PRIVATE - INPUT*/;
    int* k_darts2506 /*OMP_PRIVATE - INPUT*/;
    double** q_darts2506 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21_darts2506 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration2506;
    int lastIteration2506;
    int range2506;
    int rangePerCodelet2506;
    int minIteration2506;
    int remainderRange2506;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets2506* barrierCodelets2506;
    _checkInCodelets2507* checkInCodelets2507;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets2507* firstCodelet;
#endif
    TP2506(int in_numThreads, int in_mainCodeletID, TP2301* in_TPParent, int in_initIteration,
        int in_lastIteration, TP2506** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP2506();
};
/*TP2655: OMPForDirective*/
class TP2655 : public ompTP {
public:
    class _barrierCodelets2655 : public darts::Codelet {
    public:
        TP2655* inputsTPParent;
        _barrierCodelets2655()
            : darts::Codelet()
        {
        }
        _barrierCodelets2655(uint32_t dep, uint32_t res, TP2655* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations2655(int* endRange, uint32_t codeletID);
    class _checkInCodelets2656 : public darts::Codelet {
    public:
        TP2655* myTP;
        TP2655* inputsTPParent;
        int endRange;
        _checkInCodelets2656()
            : darts::Codelet()
        {
        }
        _checkInCodelets2656(uint32_t dep, uint32_t res, TP2655* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2301* TPParent;
    TP2655* controlTPParent;
    TP2655* inputsTPParent;
    int** L2_darts2655 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** dsspm_darts2655 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts2655 /*OMP_PRIVATE - INPUT*/;
    int** iend1_darts2655 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** ist1_darts2655 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts2655 /*OMP_PRIVATE - INPUT*/;
    int* k_darts2655 /*OMP_PRIVATE - INPUT*/;
    int* m_darts2655 /*OMP_PRIVATE - INPUT*/;
    double** tmp_darts2655 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21i_darts2655 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21im1_darts2655 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31i_darts2655 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31im1_darts2655 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41i_darts2655 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41im1_darts2655 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51i_darts2655 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51im1_darts2655 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration2655;
    int lastIteration2655;
    int range2655;
    int rangePerCodelet2655;
    int minIteration2655;
    int remainderRange2655;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets2655* barrierCodelets2655;
    _checkInCodelets2656* checkInCodelets2656;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets2656* firstCodelet;
#endif
    TP2655(int in_numThreads, int in_mainCodeletID, TP2301* in_TPParent, int in_initIteration,
        int in_lastIteration, TP2655** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP2655();
};
/*TP3293: OMPForDirective*/
class TP3293 : public ompTP {
public:
    class _barrierCodelets3293 : public darts::Codelet {
    public:
        TP3293* inputsTPParent;
        _barrierCodelets3293()
            : darts::Codelet()
        {
        }
        _barrierCodelets3293(uint32_t dep, uint32_t res, TP3293* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations3293(int* endRange, uint32_t codeletID);
    class _checkInCodelets3294 : public darts::Codelet {
    public:
        TP3293* myTP;
        TP3293* inputsTPParent;
        int endRange;
        _checkInCodelets3294()
            : darts::Codelet()
        {
        }
        _checkInCodelets3294(uint32_t dep, uint32_t res, TP3293* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2301* TPParent;
    TP3293* controlTPParent;
    TP3293* inputsTPParent;
    int** L1_darts3293 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** L2_darts3293 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts3293 /*OMP_PRIVATE - INPUT*/;
    int* j_darts3293 /*OMP_PRIVATE - INPUT*/;
    int* k_darts3293 /*OMP_PRIVATE - INPUT*/;
    double** q_darts3293 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31_darts3293 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration3293;
    int lastIteration3293;
    int range3293;
    int rangePerCodelet3293;
    int minIteration3293;
    int remainderRange3293;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets3293* barrierCodelets3293;
    _checkInCodelets3294* checkInCodelets3294;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets3294* firstCodelet;
#endif
    TP3293(int in_numThreads, int in_mainCodeletID, TP2301* in_TPParent, int in_initIteration,
        int in_lastIteration, TP3293** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP3293();
};
/*TP3442: OMPForDirective*/
class TP3442 : public ompTP {
public:
    class _barrierCodelets3442 : public darts::Codelet {
    public:
        TP3442* inputsTPParent;
        _barrierCodelets3442()
            : darts::Codelet()
        {
        }
        _barrierCodelets3442(uint32_t dep, uint32_t res, TP3442* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations3442(int* endRange, uint32_t codeletID);
    class _checkInCodelets3443 : public darts::Codelet {
    public:
        TP3442* myTP;
        TP3442* inputsTPParent;
        int endRange;
        _checkInCodelets3443()
            : darts::Codelet()
        {
        }
        _checkInCodelets3443(uint32_t dep, uint32_t res, TP3442* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2301* TPParent;
    TP3442* controlTPParent;
    TP3442* inputsTPParent;
    int** L2_darts3442 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** dsspm_darts3442 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts3442 /*OMP_PRIVATE - INPUT*/;
    int* j_darts3442 /*OMP_PRIVATE - INPUT*/;
    int** jend1_darts3442 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** jst1_darts3442 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts3442 /*OMP_PRIVATE - INPUT*/;
    int* m_darts3442 /*OMP_PRIVATE - INPUT*/;
    double** tmp_darts3442 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21j_darts3442 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21jm1_darts3442 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31j_darts3442 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31jm1_darts3442 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41j_darts3442 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41jm1_darts3442 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51j_darts3442 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51jm1_darts3442 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration3442;
    int lastIteration3442;
    int range3442;
    int rangePerCodelet3442;
    int minIteration3442;
    int remainderRange3442;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets3442* barrierCodelets3442;
    _checkInCodelets3443* checkInCodelets3443;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets3443* firstCodelet;
#endif
    TP3442(int in_numThreads, int in_mainCodeletID, TP2301* in_TPParent, int in_initIteration,
        int in_lastIteration, TP3442** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP3442();
};
/*TP4077: OMPForDirective*/
class TP4077 : public ompTP {
public:
    class _barrierCodelets4077 : public darts::Codelet {
    public:
        TP4077* inputsTPParent;
        _barrierCodelets4077()
            : darts::Codelet()
        {
        }
        _barrierCodelets4077(uint32_t dep, uint32_t res, TP4077* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations4077(int* endRange, uint32_t codeletID);
    class _checkInCodelets4078 : public darts::Codelet {
    public:
        TP4077* myTP;
        TP4077* inputsTPParent;
        int endRange;
        _checkInCodelets4078()
            : darts::Codelet()
        {
        }
        _checkInCodelets4078(uint32_t dep, uint32_t res, TP4077* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP2301* TPParent;
    TP4077* controlTPParent;
    TP4077* inputsTPParent;
    double** dsspm_darts4077 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts4077 /*OMP_PRIVATE - INPUT*/;
    int* j_darts4077 /*OMP_PRIVATE - INPUT*/;
    int* k_darts4077 /*OMP_PRIVATE - INPUT*/;
    int* m_darts4077 /*OMP_PRIVATE - INPUT*/;
    double** q_darts4077 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp_darts4077 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21k_darts4077 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21km1_darts4077 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31k_darts4077 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31km1_darts4077 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41_darts4077 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41k_darts4077 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41km1_darts4077 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51k_darts4077 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51km1_darts4077 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration4077;
    int lastIteration4077;
    int range4077;
    int rangePerCodelet4077;
    int minIteration4077;
    int remainderRange4077;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets4077* barrierCodelets4077;
    _checkInCodelets4078* checkInCodelets4078;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets4078* firstCodelet;
#endif
    TP4077(int in_numThreads, int in_mainCodeletID, TP2301* in_TPParent, int in_initIteration,
        int in_lastIteration, TP4077** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP4077();
};
/*TP7: jacld*/
class TP7 : public ompTP {
public:
    class _checkInCodelets4982 : public darts::Codelet {
    public:
        TP7* myTP;
        TP7* inputsTPParent;
        _checkInCodelets4982()
            : darts::Codelet()
        {
        }
        _checkInCodelets4982(uint32_t dep, uint32_t res, TP7* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets5003 : public darts::Codelet {
    public:
        TP7* myTP;
        TP7* inputsTPParent;
        _checkInCodelets5003()
            : darts::Codelet()
        {
        }
        _checkInCodelets5003(uint32_t dep, uint32_t res, TP7* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
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
    TP5003** TP5003Ptr;
    size_t* TP5003_alreadyLaunched;
    int numTPsSet5003;
    int numTPsReady5003;
    size_t TPsToUse5003;
    size_t codeletsPerTP5003;
    size_t totalCodelets5003;
    _checkInCodelets4982* checkInCodelets4982;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets4982* firstCodelet;
#endif
    _checkInCodelets5003* checkInCodelets5003;
    TP7(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
        darts::Codelet* in_mainSyncCodelet, TP7** in_ptrToThisFunctionTP, int in_k);
    ~TP7();
    void setNewInputs(int in_k, size_t codeletID);
};
/*TP5003: OMPForDirective*/
class TP5003 : public ompTP {
public:
    bool requestNewRangeIterations5003(int* endRange, uint32_t codeletID);
    class _checkInCodelets5004 : public darts::Codelet {
    public:
        TP5003* myTP;
        TP5003* inputsTPParent;
        int endRange;
        _checkInCodelets5004()
            : darts::Codelet()
        {
        }
        _checkInCodelets5004(uint32_t dep, uint32_t res, TP5003* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP7* TPParent;
    TP5003* controlTPParent;
    TP5003* inputsTPParent;
    double** c1345_darts5003 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** c34_darts5003 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts5003 /*OMP_PRIVATE - INPUT*/;
    int* j_darts5003 /*OMP_PRIVATE - INPUT*/;
    int** k_darts5003 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** r43_darts5003 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp1_darts5003 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp2_darts5003 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp3_darts5003 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration5003;
    int lastIteration5003;
    int range5003;
    int rangePerCodelet5003;
    int minIteration5003;
    int remainderRange5003;
    size_t readyCodelets;
    int baseNumThreads;
    int* signalNextReady;
    _checkInCodelets5004* checkInCodelets5004;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets5004* firstCodelet;
#endif
    TP5003(int in_numThreads, int in_mainCodeletID, TP7* in_TPParent, int in_initIteration,
        int in_lastIteration, TP5003** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP5003();
};
/*TP8: jacu*/
class TP8 : public ompTP {
public:
    class _checkInCodelets7501 : public darts::Codelet {
    public:
        TP8* myTP;
        TP8* inputsTPParent;
        _checkInCodelets7501()
            : darts::Codelet()
        {
        }
        _checkInCodelets7501(uint32_t dep, uint32_t res, TP8* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets7522 : public darts::Codelet {
    public:
        TP8* myTP;
        TP8* inputsTPParent;
        _checkInCodelets7522()
            : darts::Codelet()
        {
        }
        _checkInCodelets7522(uint32_t dep, uint32_t res, TP8* myTP, uint32_t id)
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
    TP7522** TP7522Ptr;
    size_t* TP7522_alreadyLaunched;
    int numTPsSet7522;
    int numTPsReady7522;
    size_t TPsToUse7522;
    size_t codeletsPerTP7522;
    size_t totalCodelets7522;
    _checkInCodelets7501* checkInCodelets7501;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets7501* firstCodelet;
#endif
    _checkInCodelets7522* checkInCodelets7522;
    TP8(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
        darts::Codelet* in_mainSyncCodelet, TP8** in_ptrToThisFunctionTP, int in_k);
    ~TP8();
    void setNewInputs(int in_k, size_t codeletID);
};
/*TP7522: OMPForDirective*/
class TP7522 : public ompTP {
public:
    bool requestNewRangeIterations7522(int* endRange, uint32_t codeletID);
    class _checkInCodelets7523 : public darts::Codelet {
    public:
        TP7522* myTP;
        TP7522* inputsTPParent;
        int endRange;
        _checkInCodelets7523()
            : darts::Codelet()
        {
        }
        _checkInCodelets7523(uint32_t dep, uint32_t res, TP7522* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP8* TPParent;
    TP7522* controlTPParent;
    TP7522* inputsTPParent;
    double** c1345_darts7522 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** c34_darts7522 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts7522 /*OMP_PRIVATE - INPUT*/;
    int* j_darts7522 /*OMP_PRIVATE - INPUT*/;
    int** k_darts7522 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** r43_darts7522 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp1_darts7522 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp2_darts7522 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp3_darts7522 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration7522;
    int lastIteration7522;
    int range7522;
    int rangePerCodelet7522;
    int minIteration7522;
    int remainderRange7522;
    size_t readyCodelets;
    int baseNumThreads;
    int* signalNextReady;
    _checkInCodelets7523* checkInCodelets7523;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets7523* firstCodelet;
#endif
    TP7522(int in_numThreads, int in_mainCodeletID, TP8* in_TPParent, int in_initIteration,
        int in_lastIteration, TP7522** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP7522();
};
/*TP10896: OMPParallelDirective*/
class TP10896 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets10896 : public darts::Codelet {
    public:
        TP10896* inputsTPParent;
        _barrierCodelets10896()
            : darts::Codelet()
        {
        }
        _barrierCodelets10896(uint32_t dep, uint32_t res, TP10896* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets10911 : public darts::Codelet {
    public:
        TP10896* myTP;
        TP10896* inputsTPParent;
        _checkInCodelets10911()
            : darts::Codelet()
        {
        }
        _checkInCodelets10911(uint32_t dep, uint32_t res, TP10896* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets10911 : public darts::Codelet {
    public:
        TP10896* inputsTPParent;
        _barrierCodelets10911()
            : darts::Codelet()
        {
        }
        _barrierCodelets10911(uint32_t dep, uint32_t res, TP10896* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets10968 : public darts::Codelet {
    public:
        TP10896* myTP;
        TP10896* inputsTPParent;
        _checkInCodelets10968()
            : darts::Codelet()
        {
        }
        _checkInCodelets10968(uint32_t dep, uint32_t res, TP10896* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets10971 : public darts::Codelet {
    public:
        TP10896* myTP;
        TP10896* inputsTPParent;
        _checkInCodelets10971()
            : darts::Codelet()
        {
        }
        _checkInCodelets10971(uint32_t dep, uint32_t res, TP10896* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets10971 : public darts::Codelet {
    public:
        TP10896* inputsTPParent;
        _barrierCodelets10971()
            : darts::Codelet()
        {
        }
        _barrierCodelets10971(uint32_t dep, uint32_t res, TP10896* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets11120 : public darts::Codelet {
    public:
        TP10896* myTP;
        TP10896* inputsTPParent;
        _checkInCodelets11120()
            : darts::Codelet()
        {
        }
        _checkInCodelets11120(uint32_t dep, uint32_t res, TP10896* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets11120 : public darts::Codelet {
    public:
        TP10896* inputsTPParent;
        _barrierCodelets11120()
            : darts::Codelet()
        {
        }
        _barrierCodelets11120(uint32_t dep, uint32_t res, TP10896* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets11757 : public darts::Codelet {
    public:
        TP10896* myTP;
        TP10896* inputsTPParent;
        _checkInCodelets11757()
            : darts::Codelet()
        {
        }
        _checkInCodelets11757(uint32_t dep, uint32_t res, TP10896* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets11760 : public darts::Codelet {
    public:
        TP10896* myTP;
        TP10896* inputsTPParent;
        _checkInCodelets11760()
            : darts::Codelet()
        {
        }
        _checkInCodelets11760(uint32_t dep, uint32_t res, TP10896* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets11760 : public darts::Codelet {
    public:
        TP10896* inputsTPParent;
        _barrierCodelets11760()
            : darts::Codelet()
        {
        }
        _barrierCodelets11760(uint32_t dep, uint32_t res, TP10896* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets11909 : public darts::Codelet {
    public:
        TP10896* myTP;
        TP10896* inputsTPParent;
        _checkInCodelets11909()
            : darts::Codelet()
        {
        }
        _checkInCodelets11909(uint32_t dep, uint32_t res, TP10896* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets11909 : public darts::Codelet {
    public:
        TP10896* inputsTPParent;
        _barrierCodelets11909()
            : darts::Codelet()
        {
        }
        _barrierCodelets11909(uint32_t dep, uint32_t res, TP10896* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets12546 : public darts::Codelet {
    public:
        TP10896* myTP;
        TP10896* inputsTPParent;
        _checkInCodelets12546()
            : darts::Codelet()
        {
        }
        _checkInCodelets12546(uint32_t dep, uint32_t res, TP10896* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets12546 : public darts::Codelet {
    public:
        TP10896* inputsTPParent;
        _barrierCodelets12546()
            : darts::Codelet()
        {
        }
        _barrierCodelets12546(uint32_t dep, uint32_t res, TP10896* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP10896* TPParent;
    TP10896* controlTPParent;
    TP10896* inputsTPParent;
    int* L1_darts10896 /*VARIABLE*/;
    int* L2_darts10896 /*VARIABLE*/;
    int* i_darts10896 /*VARIABLE*/;
    int* iend1_darts10896 /*VARIABLE*/;
    int* ist1_darts10896 /*VARIABLE*/;
    int* j_darts10896 /*VARIABLE*/;
    int* jend1_darts10896 /*VARIABLE*/;
    int* jst1_darts10896 /*VARIABLE*/;
    int* k_darts10896 /*VARIABLE*/;
    int* m_darts10896 /*VARIABLE*/;
    double* q_darts10896 /*VARIABLE*/;
    double* tmp_darts10896 /*VARIABLE*/;
    double* u21_darts10896 /*VARIABLE*/;
    double* u21i_darts10896 /*VARIABLE*/;
    double* u21im1_darts10896 /*VARIABLE*/;
    double* u21j_darts10896 /*VARIABLE*/;
    double* u21jm1_darts10896 /*VARIABLE*/;
    double* u21k_darts10896 /*VARIABLE*/;
    double* u21km1_darts10896 /*VARIABLE*/;
    double* u31_darts10896 /*VARIABLE*/;
    double* u31i_darts10896 /*VARIABLE*/;
    double* u31im1_darts10896 /*VARIABLE*/;
    double* u31j_darts10896 /*VARIABLE*/;
    double* u31jm1_darts10896 /*VARIABLE*/;
    double* u31k_darts10896 /*VARIABLE*/;
    double* u31km1_darts10896 /*VARIABLE*/;
    double* u41_darts10896 /*VARIABLE*/;
    double* u41i_darts10896 /*VARIABLE*/;
    double* u41im1_darts10896 /*VARIABLE*/;
    double* u41j_darts10896 /*VARIABLE*/;
    double* u41jm1_darts10896 /*VARIABLE*/;
    double* u41k_darts10896 /*VARIABLE*/;
    double* u41km1_darts10896 /*VARIABLE*/;
    double* u51i_darts10896 /*VARIABLE*/;
    double* u51im1_darts10896 /*VARIABLE*/;
    double* u51j_darts10896 /*VARIABLE*/;
    double* u51jm1_darts10896 /*VARIABLE*/;
    double* u51k_darts10896 /*VARIABLE*/;
    double* u51km1_darts10896 /*VARIABLE*/;
    TP10911** TP10911Ptr;
    size_t* TP10911_alreadyLaunched;
    int numTPsSet10911;
    int numTPsReady10911;
    size_t TPsToUse10911;
    size_t codeletsPerTP10911;
    size_t totalCodelets10911;
    TP10971** TP10971Ptr;
    size_t* TP10971_alreadyLaunched;
    int numTPsSet10971;
    int numTPsReady10971;
    size_t TPsToUse10971;
    size_t codeletsPerTP10971;
    size_t totalCodelets10971;
    TP11120** TP11120Ptr;
    size_t* TP11120_alreadyLaunched;
    int numTPsSet11120;
    int numTPsReady11120;
    size_t TPsToUse11120;
    size_t codeletsPerTP11120;
    size_t totalCodelets11120;
    TP11760** TP11760Ptr;
    size_t* TP11760_alreadyLaunched;
    int numTPsSet11760;
    int numTPsReady11760;
    size_t TPsToUse11760;
    size_t codeletsPerTP11760;
    size_t totalCodelets11760;
    TP11909** TP11909Ptr;
    size_t* TP11909_alreadyLaunched;
    int numTPsSet11909;
    int numTPsReady11909;
    size_t TPsToUse11909;
    size_t codeletsPerTP11909;
    size_t totalCodelets11909;
    TP12546** TP12546Ptr;
    size_t* TP12546_alreadyLaunched;
    int numTPsSet12546;
    int numTPsReady12546;
    size_t TPsToUse12546;
    size_t codeletsPerTP12546;
    size_t totalCodelets12546;
    _barrierCodelets10896* barrierCodelets10896;
    _checkInCodelets10911* checkInCodelets10911;
    _barrierCodelets10911* barrierCodelets10911;
    _checkInCodelets10968* checkInCodelets10968;
    _checkInCodelets10971* checkInCodelets10971;
    _barrierCodelets10971* barrierCodelets10971;
    _checkInCodelets11120* checkInCodelets11120;
    _barrierCodelets11120* barrierCodelets11120;
    _checkInCodelets11757* checkInCodelets11757;
    _checkInCodelets11760* checkInCodelets11760;
    _barrierCodelets11760* barrierCodelets11760;
    _checkInCodelets11909* checkInCodelets11909;
    _barrierCodelets11909* barrierCodelets11909;
    _checkInCodelets12546* checkInCodelets12546;
    _barrierCodelets12546* barrierCodelets12546;
    TP10896(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
    ~TP10896();
};
/*TP10911: OMPForDirective*/
class TP10911 : public ompTP {
public:
    class _barrierCodelets10911 : public darts::Codelet {
    public:
        TP10911* inputsTPParent;
        _barrierCodelets10911()
            : darts::Codelet()
        {
        }
        _barrierCodelets10911(uint32_t dep, uint32_t res, TP10911* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations10911(int* endRange, uint32_t codeletID);
    class _checkInCodelets10912 : public darts::Codelet {
    public:
        TP10911* myTP;
        TP10911* inputsTPParent;
        int endRange;
        _checkInCodelets10912()
            : darts::Codelet()
        {
        }
        _checkInCodelets10912(uint32_t dep, uint32_t res, TP10911* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10896* TPParent;
    TP10911* controlTPParent;
    TP10911* inputsTPParent;
    int* i_darts10911 /*OMP_PRIVATE - INPUT*/;
    int* j_darts10911 /*OMP_PRIVATE - INPUT*/;
    int* k_darts10911 /*OMP_PRIVATE - INPUT*/;
    int* m_darts10911 /*OMP_PRIVATE - INPUT*/;
    int initIteration10911;
    int lastIteration10911;
    int range10911;
    int rangePerCodelet10911;
    int minIteration10911;
    int remainderRange10911;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets10911* barrierCodelets10911;
    _checkInCodelets10912* checkInCodelets10912;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets10912* firstCodelet;
#endif
    TP10911(int in_numThreads, int in_mainCodeletID, TP10896* in_TPParent, int in_initIteration,
        int in_lastIteration, TP10911** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP10911();
};
/*TP10971: OMPForDirective*/
class TP10971 : public ompTP {
public:
    class _barrierCodelets10971 : public darts::Codelet {
    public:
        TP10971* inputsTPParent;
        _barrierCodelets10971()
            : darts::Codelet()
        {
        }
        _barrierCodelets10971(uint32_t dep, uint32_t res, TP10971* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations10971(int* endRange, uint32_t codeletID);
    class _checkInCodelets10972 : public darts::Codelet {
    public:
        TP10971* myTP;
        TP10971* inputsTPParent;
        int endRange;
        _checkInCodelets10972()
            : darts::Codelet()
        {
        }
        _checkInCodelets10972(uint32_t dep, uint32_t res, TP10971* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10896* TPParent;
    TP10971* controlTPParent;
    TP10971* inputsTPParent;
    int** L1_darts10971 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** L2_darts10971 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts10971 /*OMP_PRIVATE - INPUT*/;
    int* j_darts10971 /*OMP_PRIVATE - INPUT*/;
    int* k_darts10971 /*OMP_PRIVATE - INPUT*/;
    double** q_darts10971 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21_darts10971 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration10971;
    int lastIteration10971;
    int range10971;
    int rangePerCodelet10971;
    int minIteration10971;
    int remainderRange10971;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets10971* barrierCodelets10971;
    _checkInCodelets10972* checkInCodelets10972;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets10972* firstCodelet;
#endif
    TP10971(int in_numThreads, int in_mainCodeletID, TP10896* in_TPParent, int in_initIteration,
        int in_lastIteration, TP10971** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP10971();
};
/*TP11120: OMPForDirective*/
class TP11120 : public ompTP {
public:
    class _barrierCodelets11120 : public darts::Codelet {
    public:
        TP11120* inputsTPParent;
        _barrierCodelets11120()
            : darts::Codelet()
        {
        }
        _barrierCodelets11120(uint32_t dep, uint32_t res, TP11120* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations11120(int* endRange, uint32_t codeletID);
    class _checkInCodelets11121 : public darts::Codelet {
    public:
        TP11120* myTP;
        TP11120* inputsTPParent;
        int endRange;
        _checkInCodelets11121()
            : darts::Codelet()
        {
        }
        _checkInCodelets11121(uint32_t dep, uint32_t res, TP11120* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10896* TPParent;
    TP11120* controlTPParent;
    TP11120* inputsTPParent;
    int** L2_darts11120 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts11120 /*OMP_PRIVATE - INPUT*/;
    int** iend1_darts11120 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** ist1_darts11120 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts11120 /*OMP_PRIVATE - INPUT*/;
    int* k_darts11120 /*OMP_PRIVATE - INPUT*/;
    int* m_darts11120 /*OMP_PRIVATE - INPUT*/;
    double** tmp_darts11120 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21i_darts11120 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21im1_darts11120 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31i_darts11120 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31im1_darts11120 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41i_darts11120 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41im1_darts11120 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51i_darts11120 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51im1_darts11120 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration11120;
    int lastIteration11120;
    int range11120;
    int rangePerCodelet11120;
    int minIteration11120;
    int remainderRange11120;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets11120* barrierCodelets11120;
    _checkInCodelets11121* checkInCodelets11121;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets11121* firstCodelet;
#endif
    TP11120(int in_numThreads, int in_mainCodeletID, TP10896* in_TPParent, int in_initIteration,
        int in_lastIteration, TP11120** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP11120();
};
/*TP11760: OMPForDirective*/
class TP11760 : public ompTP {
public:
    class _barrierCodelets11760 : public darts::Codelet {
    public:
        TP11760* inputsTPParent;
        _barrierCodelets11760()
            : darts::Codelet()
        {
        }
        _barrierCodelets11760(uint32_t dep, uint32_t res, TP11760* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations11760(int* endRange, uint32_t codeletID);
    class _checkInCodelets11761 : public darts::Codelet {
    public:
        TP11760* myTP;
        TP11760* inputsTPParent;
        int endRange;
        _checkInCodelets11761()
            : darts::Codelet()
        {
        }
        _checkInCodelets11761(uint32_t dep, uint32_t res, TP11760* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10896* TPParent;
    TP11760* controlTPParent;
    TP11760* inputsTPParent;
    int** L1_darts11760 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** L2_darts11760 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts11760 /*OMP_PRIVATE - INPUT*/;
    int* j_darts11760 /*OMP_PRIVATE - INPUT*/;
    int* k_darts11760 /*OMP_PRIVATE - INPUT*/;
    double** q_darts11760 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31_darts11760 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration11760;
    int lastIteration11760;
    int range11760;
    int rangePerCodelet11760;
    int minIteration11760;
    int remainderRange11760;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets11760* barrierCodelets11760;
    _checkInCodelets11761* checkInCodelets11761;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets11761* firstCodelet;
#endif
    TP11760(int in_numThreads, int in_mainCodeletID, TP10896* in_TPParent, int in_initIteration,
        int in_lastIteration, TP11760** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP11760();
};
/*TP11909: OMPForDirective*/
class TP11909 : public ompTP {
public:
    class _barrierCodelets11909 : public darts::Codelet {
    public:
        TP11909* inputsTPParent;
        _barrierCodelets11909()
            : darts::Codelet()
        {
        }
        _barrierCodelets11909(uint32_t dep, uint32_t res, TP11909* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations11909(int* endRange, uint32_t codeletID);
    class _checkInCodelets11910 : public darts::Codelet {
    public:
        TP11909* myTP;
        TP11909* inputsTPParent;
        int endRange;
        _checkInCodelets11910()
            : darts::Codelet()
        {
        }
        _checkInCodelets11910(uint32_t dep, uint32_t res, TP11909* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10896* TPParent;
    TP11909* controlTPParent;
    TP11909* inputsTPParent;
    int** L2_darts11909 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts11909 /*OMP_PRIVATE - INPUT*/;
    int* j_darts11909 /*OMP_PRIVATE - INPUT*/;
    int** jend1_darts11909 /*OMP_SHARED_PRIVATE - INPUT*/;
    int** jst1_darts11909 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts11909 /*OMP_PRIVATE - INPUT*/;
    int* m_darts11909 /*OMP_PRIVATE - INPUT*/;
    double** tmp_darts11909 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21j_darts11909 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21jm1_darts11909 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31j_darts11909 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31jm1_darts11909 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41j_darts11909 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41jm1_darts11909 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51j_darts11909 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51jm1_darts11909 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration11909;
    int lastIteration11909;
    int range11909;
    int rangePerCodelet11909;
    int minIteration11909;
    int remainderRange11909;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets11909* barrierCodelets11909;
    _checkInCodelets11910* checkInCodelets11910;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets11910* firstCodelet;
#endif
    TP11909(int in_numThreads, int in_mainCodeletID, TP10896* in_TPParent, int in_initIteration,
        int in_lastIteration, TP11909** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP11909();
};
/*TP12546: OMPForDirective*/
class TP12546 : public ompTP {
public:
    class _barrierCodelets12546 : public darts::Codelet {
    public:
        TP12546* inputsTPParent;
        _barrierCodelets12546()
            : darts::Codelet()
        {
        }
        _barrierCodelets12546(uint32_t dep, uint32_t res, TP12546* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations12546(int* endRange, uint32_t codeletID);
    class _checkInCodelets12547 : public darts::Codelet {
    public:
        TP12546* myTP;
        TP12546* inputsTPParent;
        int endRange;
        _checkInCodelets12547()
            : darts::Codelet()
        {
        }
        _checkInCodelets12547(uint32_t dep, uint32_t res, TP12546* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP10896* TPParent;
    TP12546* controlTPParent;
    TP12546* inputsTPParent;
    int* i_darts12546 /*OMP_PRIVATE - INPUT*/;
    int* j_darts12546 /*OMP_PRIVATE - INPUT*/;
    int* k_darts12546 /*OMP_PRIVATE - INPUT*/;
    int* m_darts12546 /*OMP_PRIVATE - INPUT*/;
    double** q_darts12546 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** tmp_darts12546 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21k_darts12546 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u21km1_darts12546 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31k_darts12546 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u31km1_darts12546 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41_darts12546 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41k_darts12546 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u41km1_darts12546 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51k_darts12546 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** u51km1_darts12546 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration12546;
    int lastIteration12546;
    int range12546;
    int rangePerCodelet12546;
    int minIteration12546;
    int remainderRange12546;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets12546* barrierCodelets12546;
    _checkInCodelets12547* checkInCodelets12547;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets12547* firstCodelet;
#endif
    TP12546(int in_numThreads, int in_mainCodeletID, TP10896* in_TPParent, int in_initIteration,
        int in_lastIteration, TP12546** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP12546();
};
/*TP13297: OMPParallelDirective*/
class TP13297 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets13297 : public darts::Codelet {
    public:
        TP13297* inputsTPParent;
        _barrierCodelets13297()
            : darts::Codelet()
        {
        }
        _barrierCodelets13297(uint32_t dep, uint32_t res, TP13297* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13301 : public darts::Codelet {
    public:
        TP13297* myTP;
        TP13297* inputsTPParent;
        _checkInCodelets13301()
            : darts::Codelet()
        {
        }
        _checkInCodelets13301(uint32_t dep, uint32_t res, TP13297* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13301 : public darts::Codelet {
    public:
        TP13297* inputsTPParent;
        _barrierCodelets13301()
            : darts::Codelet()
        {
        }
        _barrierCodelets13301(uint32_t dep, uint32_t res, TP13297* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13352 : public darts::Codelet {
    public:
        TP13297* myTP;
        TP13297* inputsTPParent;
        _checkInCodelets13352()
            : darts::Codelet()
        {
        }
        _checkInCodelets13352(uint32_t dep, uint32_t res, TP13297* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13352 : public darts::Codelet {
    public:
        TP13297* inputsTPParent;
        _barrierCodelets13352()
            : darts::Codelet()
        {
        }
        _barrierCodelets13352(uint32_t dep, uint32_t res, TP13297* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13394 : public darts::Codelet {
    public:
        TP13297* myTP;
        TP13297* inputsTPParent;
        _checkInCodelets13394()
            : darts::Codelet()
        {
        }
        _checkInCodelets13394(uint32_t dep, uint32_t res, TP13297* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13394 : public darts::Codelet {
    public:
        TP13297* inputsTPParent;
        _barrierCodelets13394()
            : darts::Codelet()
        {
        }
        _barrierCodelets13394(uint32_t dep, uint32_t res, TP13297* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13438 : public darts::Codelet {
    public:
        TP13297* myTP;
        TP13297* inputsTPParent;
        _checkInCodelets13438()
            : darts::Codelet()
        {
        }
        _checkInCodelets13438(uint32_t dep, uint32_t res, TP13297* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13438 : public darts::Codelet {
    public:
        TP13297* inputsTPParent;
        _barrierCodelets13438()
            : darts::Codelet()
        {
        }
        _barrierCodelets13438(uint32_t dep, uint32_t res, TP13297* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13480 : public darts::Codelet {
    public:
        TP13297* myTP;
        TP13297* inputsTPParent;
        _checkInCodelets13480()
            : darts::Codelet()
        {
        }
        _checkInCodelets13480(uint32_t dep, uint32_t res, TP13297* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13480 : public darts::Codelet {
    public:
        TP13297* inputsTPParent;
        _barrierCodelets13480()
            : darts::Codelet()
        {
        }
        _barrierCodelets13480(uint32_t dep, uint32_t res, TP13297* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP13297* TPParent;
    TP13297* controlTPParent;
    TP13297* inputsTPParent;
    int* i_darts13297 /*VARIABLE*/;
    int* iglob_darts13297 /*VARIABLE*/;
    int* j_darts13297 /*VARIABLE*/;
    int* jglob_darts13297 /*VARIABLE*/;
    int* k_darts13297 /*VARIABLE*/;
    TP13301** TP13301Ptr;
    size_t* TP13301_alreadyLaunched;
    int numTPsSet13301;
    int numTPsReady13301;
    size_t TPsToUse13301;
    size_t codeletsPerTP13301;
    size_t totalCodelets13301;
    TP13352** TP13352Ptr;
    size_t* TP13352_alreadyLaunched;
    int numTPsSet13352;
    int numTPsReady13352;
    size_t TPsToUse13352;
    size_t codeletsPerTP13352;
    size_t totalCodelets13352;
    TP13394** TP13394Ptr;
    size_t* TP13394_alreadyLaunched;
    int numTPsSet13394;
    int numTPsReady13394;
    size_t TPsToUse13394;
    size_t codeletsPerTP13394;
    size_t totalCodelets13394;
    TP13438** TP13438Ptr;
    size_t* TP13438_alreadyLaunched;
    int numTPsSet13438;
    int numTPsReady13438;
    size_t TPsToUse13438;
    size_t codeletsPerTP13438;
    size_t totalCodelets13438;
    TP13480** TP13480Ptr;
    size_t* TP13480_alreadyLaunched;
    int numTPsSet13480;
    int numTPsReady13480;
    size_t TPsToUse13480;
    size_t codeletsPerTP13480;
    size_t totalCodelets13480;
    _barrierCodelets13297* barrierCodelets13297;
    _checkInCodelets13301* checkInCodelets13301;
    _barrierCodelets13301* barrierCodelets13301;
    _checkInCodelets13352* checkInCodelets13352;
    _barrierCodelets13352* barrierCodelets13352;
    _checkInCodelets13394* checkInCodelets13394;
    _barrierCodelets13394* barrierCodelets13394;
    _checkInCodelets13438* checkInCodelets13438;
    _barrierCodelets13438* barrierCodelets13438;
    _checkInCodelets13480* checkInCodelets13480;
    _barrierCodelets13480* barrierCodelets13480;
    TP13297(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
    ~TP13297();
};
/*TP13301: OMPForDirective*/
class TP13301 : public ompTP {
public:
    class _barrierCodelets13301 : public darts::Codelet {
    public:
        TP13301* inputsTPParent;
        _barrierCodelets13301()
            : darts::Codelet()
        {
        }
        _barrierCodelets13301(uint32_t dep, uint32_t res, TP13301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13301(int* endRange, uint32_t codeletID);
    class _checkInCodelets13302 : public darts::Codelet {
    public:
        TP13301* myTP;
        TP13301* inputsTPParent;
        int endRange;
        _checkInCodelets13302()
            : darts::Codelet()
        {
        }
        _checkInCodelets13302(uint32_t dep, uint32_t res, TP13301* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13297* TPParent;
    TP13301* controlTPParent;
    TP13301* inputsTPParent;
    int* i_darts13301 /*OMP_PRIVATE - INPUT*/;
    int** iglob_darts13301 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts13301 /*OMP_PRIVATE - INPUT*/;
    int** jglob_darts13301 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration13301;
    int lastIteration13301;
    int range13301;
    int rangePerCodelet13301;
    int minIteration13301;
    int remainderRange13301;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13301* barrierCodelets13301;
    _checkInCodelets13302* checkInCodelets13302;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13302* firstCodelet;
#endif
    TP13301(int in_numThreads, int in_mainCodeletID, TP13297* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13301** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13301();
};
/*TP13352: OMPForDirective*/
class TP13352 : public ompTP {
public:
    class _barrierCodelets13352 : public darts::Codelet {
    public:
        TP13352* inputsTPParent;
        _barrierCodelets13352()
            : darts::Codelet()
        {
        }
        _barrierCodelets13352(uint32_t dep, uint32_t res, TP13352* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13352(int* endRange, uint32_t codeletID);
    class _checkInCodelets13353 : public darts::Codelet {
    public:
        TP13352* myTP;
        TP13352* inputsTPParent;
        int endRange;
        _checkInCodelets13353()
            : darts::Codelet()
        {
        }
        _checkInCodelets13353(uint32_t dep, uint32_t res, TP13352* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13297* TPParent;
    TP13352* controlTPParent;
    TP13352* inputsTPParent;
    int* i_darts13352 /*OMP_PRIVATE - INPUT*/;
    int** iglob_darts13352 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts13352 /*OMP_PRIVATE - INPUT*/;
    int initIteration13352;
    int lastIteration13352;
    int range13352;
    int rangePerCodelet13352;
    int minIteration13352;
    int remainderRange13352;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13352* barrierCodelets13352;
    _checkInCodelets13353* checkInCodelets13353;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13353* firstCodelet;
#endif
    TP13352(int in_numThreads, int in_mainCodeletID, TP13297* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13352** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13352();
};
/*TP13394: OMPForDirective*/
class TP13394 : public ompTP {
public:
    class _barrierCodelets13394 : public darts::Codelet {
    public:
        TP13394* inputsTPParent;
        _barrierCodelets13394()
            : darts::Codelet()
        {
        }
        _barrierCodelets13394(uint32_t dep, uint32_t res, TP13394* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13394(int* endRange, uint32_t codeletID);
    class _checkInCodelets13395 : public darts::Codelet {
    public:
        TP13394* myTP;
        TP13394* inputsTPParent;
        int endRange;
        _checkInCodelets13395()
            : darts::Codelet()
        {
        }
        _checkInCodelets13395(uint32_t dep, uint32_t res, TP13394* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13297* TPParent;
    TP13394* controlTPParent;
    TP13394* inputsTPParent;
    int* i_darts13394 /*OMP_PRIVATE - INPUT*/;
    int** iglob_darts13394 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts13394 /*OMP_PRIVATE - INPUT*/;
    int initIteration13394;
    int lastIteration13394;
    int range13394;
    int rangePerCodelet13394;
    int minIteration13394;
    int remainderRange13394;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13394* barrierCodelets13394;
    _checkInCodelets13395* checkInCodelets13395;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13395* firstCodelet;
#endif
    TP13394(int in_numThreads, int in_mainCodeletID, TP13297* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13394** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13394();
};
/*TP13438: OMPForDirective*/
class TP13438 : public ompTP {
public:
    class _barrierCodelets13438 : public darts::Codelet {
    public:
        TP13438* inputsTPParent;
        _barrierCodelets13438()
            : darts::Codelet()
        {
        }
        _barrierCodelets13438(uint32_t dep, uint32_t res, TP13438* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13438(int* endRange, uint32_t codeletID);
    class _checkInCodelets13439 : public darts::Codelet {
    public:
        TP13438* myTP;
        TP13438* inputsTPParent;
        int endRange;
        _checkInCodelets13439()
            : darts::Codelet()
        {
        }
        _checkInCodelets13439(uint32_t dep, uint32_t res, TP13438* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13297* TPParent;
    TP13438* controlTPParent;
    TP13438* inputsTPParent;
    int* j_darts13438 /*OMP_PRIVATE - INPUT*/;
    int** jglob_darts13438 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts13438 /*OMP_PRIVATE - INPUT*/;
    int initIteration13438;
    int lastIteration13438;
    int range13438;
    int rangePerCodelet13438;
    int minIteration13438;
    int remainderRange13438;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13438* barrierCodelets13438;
    _checkInCodelets13439* checkInCodelets13439;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13439* firstCodelet;
#endif
    TP13438(int in_numThreads, int in_mainCodeletID, TP13297* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13438** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13438();
};
/*TP13480: OMPForDirective*/
class TP13480 : public ompTP {
public:
    class _barrierCodelets13480 : public darts::Codelet {
    public:
        TP13480* inputsTPParent;
        _barrierCodelets13480()
            : darts::Codelet()
        {
        }
        _barrierCodelets13480(uint32_t dep, uint32_t res, TP13480* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13480(int* endRange, uint32_t codeletID);
    class _checkInCodelets13481 : public darts::Codelet {
    public:
        TP13480* myTP;
        TP13480* inputsTPParent;
        int endRange;
        _checkInCodelets13481()
            : darts::Codelet()
        {
        }
        _checkInCodelets13481(uint32_t dep, uint32_t res, TP13480* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13297* TPParent;
    TP13480* controlTPParent;
    TP13480* inputsTPParent;
    int* j_darts13480 /*OMP_PRIVATE - INPUT*/;
    int** jglob_darts13480 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts13480 /*OMP_PRIVATE - INPUT*/;
    int initIteration13480;
    int lastIteration13480;
    int range13480;
    int rangePerCodelet13480;
    int minIteration13480;
    int remainderRange13480;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13480* barrierCodelets13480;
    _checkInCodelets13481* checkInCodelets13481;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13481* firstCodelet;
#endif
    TP13480(int in_numThreads, int in_mainCodeletID, TP13297* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13480** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13480();
};
/*TP13864: OMPParallelDirective*/
class TP13864 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets13864 : public darts::Codelet {
    public:
        TP13864* inputsTPParent;
        _barrierCodelets13864()
            : darts::Codelet()
        {
        }
        _barrierCodelets13864(uint32_t dep, uint32_t res, TP13864* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13871 : public darts::Codelet {
    public:
        TP13864* myTP;
        TP13864* inputsTPParent;
        _checkInCodelets13871()
            : darts::Codelet()
        {
        }
        _checkInCodelets13871(uint32_t dep, uint32_t res, TP13864* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13871 : public darts::Codelet {
    public:
        TP13864* inputsTPParent;
        _barrierCodelets13871()
            : darts::Codelet()
        {
        }
        _barrierCodelets13871(uint32_t dep, uint32_t res, TP13864* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP13864* TPParent;
    TP13864* controlTPParent;
    TP13864* inputsTPParent;
    double* eta_darts13864 /*VARIABLE*/;
    int* i_darts13864 /*VARIABLE*/;
    int* iglob_darts13864 /*VARIABLE*/;
    int* j_darts13864 /*VARIABLE*/;
    int* jglob_darts13864 /*VARIABLE*/;
    int* k_darts13864 /*VARIABLE*/;
    int* m_darts13864 /*VARIABLE*/;
    double* peta_darts13864 /*VARIABLE*/;
    double* pxi_darts13864 /*VARIABLE*/;
    double* pzeta_darts13864 /*VARIABLE*/;
    double** ue_1jk_darts13864 /*VARIABLE*/;
    uint64_t ue_1jk_outer13864_size;
    double** ue_i1k_darts13864 /*VARIABLE*/;
    uint64_t ue_i1k_outer13864_size;
    double** ue_ij1_darts13864 /*VARIABLE*/;
    uint64_t ue_ij1_outer13864_size;
    double** ue_ijnz_darts13864 /*VARIABLE*/;
    uint64_t ue_ijnz_outer13864_size;
    double** ue_iny0k_darts13864 /*VARIABLE*/;
    uint64_t ue_iny0k_outer13864_size;
    double** ue_nx0jk_darts13864 /*VARIABLE*/;
    uint64_t ue_nx0jk_outer13864_size;
    double* xi_darts13864 /*VARIABLE*/;
    double* zeta_darts13864 /*VARIABLE*/;
    TP13871** TP13871Ptr;
    size_t* TP13871_alreadyLaunched;
    int numTPsSet13871;
    int numTPsReady13871;
    size_t TPsToUse13871;
    size_t codeletsPerTP13871;
    size_t totalCodelets13871;
    _barrierCodelets13864* barrierCodelets13864;
    _checkInCodelets13871* checkInCodelets13871;
    _barrierCodelets13871* barrierCodelets13871;
    TP13864(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
    ~TP13864();
};
/*TP13871: OMPForDirective*/
class TP13871 : public ompTP {
public:
    class _barrierCodelets13871 : public darts::Codelet {
    public:
        TP13871* inputsTPParent;
        _barrierCodelets13871()
            : darts::Codelet()
        {
        }
        _barrierCodelets13871(uint32_t dep, uint32_t res, TP13871* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13871(int* endRange, uint32_t codeletID);
    class _checkInCodelets13872 : public darts::Codelet {
    public:
        TP13871* myTP;
        TP13871* inputsTPParent;
        int endRange;
        _checkInCodelets13872()
            : darts::Codelet()
        {
        }
        _checkInCodelets13872(uint32_t dep, uint32_t res, TP13871* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13864* TPParent;
    TP13871* controlTPParent;
    TP13871* inputsTPParent;
    double** eta_darts13871 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* i_darts13871 /*OMP_PRIVATE - INPUT*/;
    int** iglob_darts13871 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* j_darts13871 /*OMP_PRIVATE - INPUT*/;
    int** jglob_darts13871 /*OMP_SHARED_PRIVATE - INPUT*/;
    int* k_darts13871 /*OMP_PRIVATE - INPUT*/;
    int* m_darts13871 /*OMP_PRIVATE - INPUT*/;
    double** peta_darts13871 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** pxi_darts13871 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** pzeta_darts13871 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** ue_1jk_darts13871 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** ue_i1k_darts13871 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** ue_ij1_darts13871 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** ue_ijnz_darts13871 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** ue_iny0k_darts13871 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** ue_nx0jk_darts13871 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** xi_darts13871 /*OMP_SHARED_PRIVATE - INPUT*/;
    double** zeta_darts13871 /*OMP_SHARED_PRIVATE - INPUT*/;
    int initIteration13871;
    int lastIteration13871;
    int range13871;
    int rangePerCodelet13871;
    int minIteration13871;
    int remainderRange13871;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13871* barrierCodelets13871;
    _checkInCodelets13872* checkInCodelets13872;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets13872* firstCodelet;
#endif
    TP13871(int in_numThreads, int in_mainCodeletID, TP13864* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13871** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13871();
};
/*TP13997: OMPParallelDirective*/
class TP13997 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets13997 : public darts::Codelet {
    public:
        TP13997* inputsTPParent;
        _barrierCodelets13997()
            : darts::Codelet()
        {
        }
        _barrierCodelets13997(uint32_t dep, uint32_t res, TP13997* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets13999 : public darts::Codelet {
    public:
        TP13997* myTP;
        TP13997* inputsTPParent;
        _checkInCodelets13999()
            : darts::Codelet()
        {
        }
        _checkInCodelets13999(uint32_t dep, uint32_t res, TP13997* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets13999 : public darts::Codelet {
    public:
        TP13997* inputsTPParent;
        _barrierCodelets13999()
            : darts::Codelet()
        {
        }
        _barrierCodelets13999(uint32_t dep, uint32_t res, TP13997* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP13997* TPParent;
    TP13997* controlTPParent;
    TP13997* inputsTPParent;
    int* i_darts13997 /*OMP_PRIVATE - INPUT*/;
    int* j_darts13997 /*OMP_PRIVATE - INPUT*/;
    int* k_darts13997 /*OMP_PRIVATE - INPUT*/;
    int* m_darts13997 /*OMP_PRIVATE - INPUT*/;
    TP13999** TP13999Ptr;
    size_t* TP13999_alreadyLaunched;
    int numTPsSet13999;
    int numTPsReady13999;
    size_t TPsToUse13999;
    size_t codeletsPerTP13999;
    size_t totalCodelets13999;
    _barrierCodelets13997* barrierCodelets13997;
    _checkInCodelets13999* checkInCodelets13999;
    _barrierCodelets13999* barrierCodelets13999;
    TP13997(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet);
    ~TP13997();
};
/*TP13999: OMPForDirective*/
class TP13999 : public ompTP {
public:
    class _barrierCodelets13999 : public darts::Codelet {
    public:
        TP13999* inputsTPParent;
        _barrierCodelets13999()
            : darts::Codelet()
        {
        }
        _barrierCodelets13999(uint32_t dep, uint32_t res, TP13999* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations13999(int* endRange, uint32_t codeletID);
    class _checkInCodelets14000 : public darts::Codelet {
    public:
        TP13999* myTP;
        TP13999* inputsTPParent;
        int endRange;
        _checkInCodelets14000()
            : darts::Codelet()
        {
        }
        _checkInCodelets14000(uint32_t dep, uint32_t res, TP13999* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP13997* TPParent;
    TP13999* controlTPParent;
    TP13999* inputsTPParent;
    int* i_darts13999 /*OMP_PRIVATE - INPUT*/;
    int* j_darts13999 /*OMP_PRIVATE - INPUT*/;
    int* k_darts13999 /*OMP_PRIVATE - INPUT*/;
    int* m_darts13999 /*OMP_PRIVATE - INPUT*/;
    int initIteration13999;
    int lastIteration13999;
    int range13999;
    int rangePerCodelet13999;
    int minIteration13999;
    int remainderRange13999;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets13999* barrierCodelets13999;
    _checkInCodelets14000* checkInCodelets14000;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets14000* firstCodelet;
#endif
    TP13999(int in_numThreads, int in_mainCodeletID, TP13997* in_TPParent, int in_initIteration,
        int in_lastIteration, TP13999** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP13999();
};
/*TP14096: OMPParallelDirective*/
class TP14096 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets14096 : public darts::Codelet {
    public:
        TP14096* inputsTPParent;
        _barrierCodelets14096()
            : darts::Codelet()
        {
        }
        _barrierCodelets14096(uint32_t dep, uint32_t res, TP14096* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14098 : public darts::Codelet {
    public:
        TP14096* myTP;
        TP14096* inputsTPParent;
        _checkInCodelets14098()
            : darts::Codelet()
        {
        }
        _checkInCodelets14098(uint32_t dep, uint32_t res, TP14096* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets14098 : public darts::Codelet {
    public:
        TP14096* inputsTPParent;
        _barrierCodelets14098()
            : darts::Codelet()
        {
        }
        _barrierCodelets14098(uint32_t dep, uint32_t res, TP14096* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14154 : public darts::Codelet {
    public:
        TP14096* myTP;
        TP14096* inputsTPParent;
        _checkInCodelets14154()
            : darts::Codelet()
        {
        }
        _checkInCodelets14154(uint32_t dep, uint32_t res, TP14096* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14153 : public darts::Codelet {
    public:
        TP14096* myTP;
        TP14096* inputsTPParent;
        _checkInCodelets14153()
            : darts::Codelet()
        {
        }
        _checkInCodelets14153(uint32_t dep, uint32_t res, TP14096* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14157 : public darts::Codelet {
    public:
        TP14096* myTP;
        TP14096* inputsTPParent;
        _checkInCodelets14157()
            : darts::Codelet()
        {
        }
        _checkInCodelets14157(uint32_t dep, uint32_t res, TP14096* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets14161 : public darts::Codelet {
    public:
        TP14096* inputsTPParent;
        _barrierCodelets14161()
            : darts::Codelet()
        {
        }
        _barrierCodelets14161(uint32_t dep, uint32_t res, TP14096* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14163 : public darts::Codelet {
    public:
        TP14096* myTP;
        TP14096* inputsTPParent;
        _checkInCodelets14163()
            : darts::Codelet()
        {
        }
        _checkInCodelets14163(uint32_t dep, uint32_t res, TP14096* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14162 : public darts::Codelet {
    public:
        TP14096* myTP;
        TP14096* inputsTPParent;
        _checkInCodelets14162()
            : darts::Codelet()
        {
        }
        _checkInCodelets14162(uint32_t dep, uint32_t res, TP14096* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14166 : public darts::Codelet {
    public:
        TP14096* myTP;
        TP14096* inputsTPParent;
        _checkInCodelets14166()
            : darts::Codelet()
        {
        }
        _checkInCodelets14166(uint32_t dep, uint32_t res, TP14096* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets14170 : public darts::Codelet {
    public:
        TP14096* inputsTPParent;
        _barrierCodelets14170()
            : darts::Codelet()
        {
        }
        _barrierCodelets14170(uint32_t dep, uint32_t res, TP14096* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14171 : public darts::Codelet {
    public:
        TP14096* myTP;
        TP14096* inputsTPParent;
        _checkInCodelets14171()
            : darts::Codelet()
        {
        }
        _checkInCodelets14171(uint32_t dep, uint32_t res, TP14096* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets14171 : public darts::Codelet {
    public:
        TP14096* inputsTPParent;
        _barrierCodelets14171()
            : darts::Codelet()
        {
        }
        _barrierCodelets14171(uint32_t dep, uint32_t res, TP14096* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP14096* TPParent;
    TP14096* controlTPParent;
    TP14096* inputsTPParent;
    int* i_darts14096 /*OMP_PRIVATE - INPUT*/;
    int* istep_darts14096 /*OMP_PRIVATE - INPUT*/;
    int* j_darts14096 /*OMP_PRIVATE - INPUT*/;
    int* k_darts14096 /*OMP_PRIVATE - INPUT*/;
    int* m_darts14096 /*OMP_PRIVATE - INPUT*/;
    double* tmp_darts14096; /*OMP_SHARED - INPUT*/
    double* tv_darts14096; /*OMP_SHARED - INPUT*/
    uint64_t tv_outer14096_size;
    TP14098** TP14098Ptr;
    size_t* TP14098_alreadyLaunched;
    int numTPsSet14098;
    int numTPsReady14098;
    size_t TPsToUse14098;
    size_t codeletsPerTP14098;
    size_t totalCodelets14098;
    unsigned int TP14153_LoopCounter;
    unsigned int* TP14153_LoopCounterPerThread;
    tbb::concurrent_vector<TP14153*> TP14153PtrVec;
    unsigned int TP14162_LoopCounter;
    unsigned int* TP14162_LoopCounterPerThread;
    tbb::concurrent_vector<TP14162*> TP14162PtrVec;
    TP14171** TP14171Ptr;
    size_t* TP14171_alreadyLaunched;
    int numTPsSet14171;
    int numTPsReady14171;
    size_t TPsToUse14171;
    size_t codeletsPerTP14171;
    size_t totalCodelets14171;
    _barrierCodelets14096* barrierCodelets14096;
    _checkInCodelets14098* checkInCodelets14098;
    _barrierCodelets14098* barrierCodelets14098;
    _checkInCodelets14154* checkInCodelets14154;
    _checkInCodelets14153* checkInCodelets14153;
    _checkInCodelets14157* checkInCodelets14157;
    _barrierCodelets14161* barrierCodelets14161;
    _checkInCodelets14163* checkInCodelets14163;
    _checkInCodelets14162* checkInCodelets14162;
    _checkInCodelets14166* checkInCodelets14166;
    _barrierCodelets14170* barrierCodelets14170;
    _checkInCodelets14171* checkInCodelets14171;
    _barrierCodelets14171* barrierCodelets14171;
    TP14096(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet, double* in_tmp,
        double* in_tv, int in_tv_outer14096_size);
    ~TP14096();
};
/*TP14098: OMPForDirective*/
class TP14098 : public ompTP {
public:
    class _barrierCodelets14098 : public darts::Codelet {
    public:
        TP14098* inputsTPParent;
        _barrierCodelets14098()
            : darts::Codelet()
        {
        }
        _barrierCodelets14098(uint32_t dep, uint32_t res, TP14098* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations14098(int* endRange, uint32_t codeletID);
    class _checkInCodelets14099 : public darts::Codelet {
    public:
        TP14098* myTP;
        TP14098* inputsTPParent;
        int endRange;
        _checkInCodelets14099()
            : darts::Codelet()
        {
        }
        _checkInCodelets14099(uint32_t dep, uint32_t res, TP14098* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP14096* TPParent;
    TP14098* controlTPParent;
    TP14098* inputsTPParent;
    int* i_darts14098 /*OMP_PRIVATE - INPUT*/;
    int* j_darts14098 /*OMP_PRIVATE - INPUT*/;
    int* k_darts14098 /*OMP_PRIVATE - INPUT*/;
    int* m_darts14098 /*OMP_PRIVATE - INPUT*/;
    int initIteration14098;
    int lastIteration14098;
    int range14098;
    int rangePerCodelet14098;
    int minIteration14098;
    int remainderRange14098;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets14098* barrierCodelets14098;
    _checkInCodelets14099* checkInCodelets14099;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets14099* firstCodelet;
#endif
    TP14098(int in_numThreads, int in_mainCodeletID, TP14096* in_TPParent, int in_initIteration,
        int in_lastIteration, TP14098** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP14098();
};
/*TP14153: ForStmt*/
class TP14153 : public ompTP {
public:
    class _checkInCodelets14159 : public darts::Codelet {
    public:
        TP14153* myTP;
        TP14096* inputsTPParent;
        _checkInCodelets14159()
            : darts::Codelet()
        {
        }
        _checkInCodelets14159(uint32_t dep, uint32_t res, TP14153* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14160 : public darts::Codelet {
    public:
        TP14153* myTP;
        TP14096* inputsTPParent;
        _checkInCodelets14160()
            : darts::Codelet()
        {
        }
        _checkInCodelets14160(uint32_t dep, uint32_t res, TP14153* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP14096* TPParent;
    TP14153* controlTPParent;
    TP14096* inputsTPParent;
    TP14153** ptrToThisTP;
    TP_jacld* TP14159Ptr;
    int TP14159_alreadyLaunched;
    TP_blts* TP14160Ptr;
    int TP14160_alreadyLaunched;
    _checkInCodelets14159* checkInCodelets14159;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets14159* firstCodelet;
#endif
    _checkInCodelets14160* checkInCodelets14160;
    TP14153(int in_numThreads, int in_mainCodeletID, TP14096* in_TPParent,
        TP14096* in_inputsTPParent, TP14153** in_ptrToThisTP);
    ~TP14153();
};
/*TP14162: ForStmt*/
class TP14162 : public ompTP {
public:
    class _checkInCodelets14168 : public darts::Codelet {
    public:
        TP14162* myTP;
        TP14096* inputsTPParent;
        _checkInCodelets14168()
            : darts::Codelet()
        {
        }
        _checkInCodelets14168(uint32_t dep, uint32_t res, TP14162* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets14169 : public darts::Codelet {
    public:
        TP14162* myTP;
        TP14096* inputsTPParent;
        _checkInCodelets14169()
            : darts::Codelet()
        {
        }
        _checkInCodelets14169(uint32_t dep, uint32_t res, TP14162* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP14096* TPParent;
    TP14162* controlTPParent;
    TP14096* inputsTPParent;
    TP14162** ptrToThisTP;
    TP_jacu* TP14168Ptr;
    int TP14168_alreadyLaunched;
    TP_buts* TP14169Ptr;
    int TP14169_alreadyLaunched;
    _checkInCodelets14168* checkInCodelets14168;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets14168* firstCodelet;
#endif
    _checkInCodelets14169* checkInCodelets14169;
    TP14162(int in_numThreads, int in_mainCodeletID, TP14096* in_TPParent,
        TP14096* in_inputsTPParent, TP14162** in_ptrToThisTP);
    ~TP14162();
};
/*TP14171: OMPForDirective*/
class TP14171 : public ompTP {
public:
    class _barrierCodelets14171 : public darts::Codelet {
    public:
        TP14171* inputsTPParent;
        _barrierCodelets14171()
            : darts::Codelet()
        {
        }
        _barrierCodelets14171(uint32_t dep, uint32_t res, TP14171* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations14171(int* endRange, uint32_t codeletID);
    class _checkInCodelets14172 : public darts::Codelet {
    public:
        TP14171* myTP;
        TP14171* inputsTPParent;
        int endRange;
        _checkInCodelets14172()
            : darts::Codelet()
        {
        }
        _checkInCodelets14172(uint32_t dep, uint32_t res, TP14171* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP14096* TPParent;
    TP14171* controlTPParent;
    TP14171* inputsTPParent;
    int* i_darts14171 /*OMP_PRIVATE - INPUT*/;
    int* j_darts14171 /*OMP_PRIVATE - INPUT*/;
    int* k_darts14171 /*OMP_PRIVATE - INPUT*/;
    int* m_darts14171 /*OMP_PRIVATE - INPUT*/;
    double* tmp_darts14171; /*OMP_SHARED - INPUT*/
    int initIteration14171;
    int lastIteration14171;
    int range14171;
    int rangePerCodelet14171;
    int minIteration14171;
    int remainderRange14171;
    size_t readyCodelets;
    int baseNumThreads;
    _barrierCodelets14171* barrierCodelets14171;
    _checkInCodelets14172* checkInCodelets14172;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets14172* firstCodelet;
#endif
    TP14171(int in_numThreads, int in_mainCodeletID, TP14096* in_TPParent, int in_initIteration,
        int in_lastIteration, double* in_tmp, TP14171** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP14171();
};
#endif
