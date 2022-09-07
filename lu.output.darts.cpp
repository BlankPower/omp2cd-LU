#include "lu.output.darts.h"
using namespace darts;
using namespace std;
static boolean flag[163];
static void verify(double xcr[5], double xce[5], double xci, char* class, boolean* verified);
static void ssor();
static void setcoeff();
static void setbv();
static void setiv();
static void rhs();
static void l2norm(int nx0, int ny0, int nz0, int ist, int iend, int jst, int jend,
    double v[162][163][163][5], double sum[5]);
static void error();
static void pintgr();
static void domain();
static void exact(int i, int j, int k, double u000ijk[5]);
static void erhs();
static void read_input();
/*Function: l2norm, ID: 9*/
static void l2norm(int nx0, int ny0, int nz0, int ist, int iend, int jst, int jend,
    double v[162][163][163][5], double sum[5])
{
    /*l2norm:9*/
    162 / 2 * 2 + 1;
    162 / 2 * 2 + 1;
}
/*Function: main, ID: 18*/
int main(int argc, char** argv)
{
    getOMPNumThreads();
    getOMPSchedulePolicy();
    getTPLoopThresholds();
    getNumTPs();
    affin = new ThreadAffinity(
        ompNumThreads / NUMTPS - 1, NUMTPS, COMPACT, getDARTSTPPolicy(), getDARTSMCPolicy());
    affinMaskRes = affin->generateMask();
    myDARTSRuntime = new Runtime(affin);
    RuntimeFinalCodelet = &(myDARTSRuntime->finalSignal);
    /*main:18*/
    /*CompoundStmt:223*/
    char class;
    boolean verified;
    double mflops;
    int nthreads = 1;
    read_input();
    domain();
    setcoeff();
    setbv();
    setiv();
    erhs();
    if (affinMaskRes) {
        myDARTSRuntime->run(launch<TP234>(
            ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet, (int*)&((nthreads))));
    }
    ssor();
    error();
    pintgr();
    verify(rsdnm, errnm, frc, &class, &verified);
    mflops = (double)itmax
        * (1984.77 * (double)nx0 * (double)ny0 * (double)nz0
            - 10923.299999999999
                * (((double)(nx0 + ny0 + nz0) / 3.) * ((double)(nx0 + ny0 + nz0) / 3.))
            + 27770.900000000001 * (double)(nx0 + ny0 + nz0) / 3. - 144010.)
        / (maxtime * 1.0E+6);
    c_print_results("LU", class, nx0, ny0, nz0, itmax, nthreads, maxtime, mflops,
        "          floating point", verified, "3.0 structured", "02 Sep 2022", "gcc", "gcc",
        "-lm -fopenmp", "-I../common", "-O3 -fopenmp", "(none)", "(none)");
}
/*Function: domain, ID: 3*/
static void domain()
{
    /*domain:3*/
    /*CompoundStmt:2280*/
    nx = nx0;
    ny = ny0;
    nz = nz0;
    if (nx < 4 || ny < 4 || nz < 4) {
        printf("     SUBDOMAIN SIZE IS TOO SMALL - \n     ADJUST PROBLEM SIZE OR NUMBER OF "
               "PROCESSORS\n     SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL\n     TO 4 THEY "
               "ARE CURRENTLY%3d%3d%3d\n",
            nx, ny, nz);
        exit(1);
    }
    if (nx > 162 || ny > 162 || nz > 162) {
        printf("     SUBDOMAIN SIZE IS TOO LARGE - \n     ADJUST PROBLEM SIZE OR NUMBER OF "
               "PROCESSORS\n     SO THAT NX, NY AND NZ ARE LESS THAN OR EQUAL TO \n     ISIZ1, "
               "ISIZ2 AND ISIZ3 RESPECTIVELY.  THEY ARE\n     CURRENTLY%4d%4d%4d\n",
            nx, ny, nz);
        exit(1);
    }
    ist = 1;
    iend = nx - 2;
    jst = 1;
    jend = ny - 2;
}
/*Function: erhs, ID: 4*/
static void erhs()
{
    /*erhs:4*/
    /*CompoundStmt:2300*/
    if (affinMaskRes) {
        myDARTSRuntime->run(
            launch<TP2301>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
    }
}
/*Function: error, ID: 5*/
static void error()
{
    /*error:5*/
    /*CompoundStmt:4827*/
    int i, j, k, m;
    int iglob, jglob;
    double tmp;
    double u000ijk[5];
    for (m = 0; m < 5; m++) {
        errnm[m] = 0.;
    }
    for (i = ist; i <= iend; i++) {
        iglob = i;
        for (j = jst; j <= jend; j++) {
            jglob = j;
            for (k = 1; k <= nz - 2; k++) {
                exact(iglob, jglob, k, u000ijk);
                for (m = 0; m < 5; m++) {
                    tmp = (u000ijk[m] - u[i][j][k][m]);
                    errnm[m] = errnm[m] + tmp * tmp;
                }
            }
        }
    }
    for (m = 0; m < 5; m++) {
        errnm[m] = sqrt(errnm[m] / ((nx0 - 2) * (ny0 - 2) * (nz0 - 2)));
    }
}
/*Function: exact, ID: 6*/
static void exact(int i, int j, int k, double u000ijk[5])
{
    /*exact:6*/
    /*CompoundStmt:4891*/
    int m;
    double xi, eta, zeta;
    xi = ((double)i) / (nx0 - 1);
    eta = ((double)j) / (ny0 - 1);
    zeta = ((double)k) / (nz - 1);
    for (m = 0; m < 5; m++) {
        u000ijk[m] = ce[m][0] + ce[m][1] * xi + ce[m][2] * eta + ce[m][3] * zeta
            + ce[m][4] * xi * xi + ce[m][5] * eta * eta + ce[m][6] * zeta * zeta
            + ce[m][7] * xi * xi * xi + ce[m][8] * eta * eta * eta + ce[m][9] * zeta * zeta * zeta
            + ce[m][10] * xi * xi * xi * xi + ce[m][11] * eta * eta * eta * eta
            + ce[m][12] * zeta * zeta * zeta * zeta;
    }
}
/*Function: l2norm, ID: 9*/
static void l2norm(int nx0, int ny0, int nz0, int ist, int iend, int jst, int jend,
    double v[162][163][163][5], double sum[5])
{
    /*l2norm:9*/
    162 / 2 * 2 + 1;
    162 / 2 * 2 + 1;
}
/*Function: pintgr, ID: 10*/
static void pintgr()
{
    /*pintgr:10*/
    /*CompoundStmt:10121*/
    int i, j, k;
    int ibeg, ifin, ifin1;
    int jbeg, jfin, jfin1;
    int iglob, iglob1, iglob2;
    int jglob, jglob1, jglob2;
    double phi1[164][164];
    double phi2[164][164];
    double frc1, frc2, frc3;
    ibeg = nx;
    ifin = 0;
    iglob1 = -1;
    iglob2 = nx - 1;
    if (iglob1 >= ii1 && iglob2 < ii2 + nx)
        ibeg = 0;
    if (iglob1 >= ii1 - nx && iglob2 <= ii2)
        ifin = nx;
    if (ii1 >= iglob1 && ii1 <= iglob2)
        ibeg = ii1;
    if (ii2 >= iglob1 && ii2 <= iglob2)
        ifin = ii2;
    jbeg = ny;
    jfin = -1;
    jglob1 = 0;
    jglob2 = ny - 1;
    if (jglob1 >= ji1 && jglob2 < ji2 + ny)
        jbeg = 0;
    if (jglob1 > ji1 - ny && jglob2 <= ji2)
        jfin = ny;
    if (ji1 >= jglob1 && ji1 <= jglob2)
        jbeg = ji1;
    if (ji2 >= jglob1 && ji2 <= jglob2)
        jfin = ji2;
    ifin1 = ifin;
    jfin1 = jfin;
    if (ifin1 == ii2)
        ifin1 = ifin - 1;
    if (jfin1 == ji2)
        jfin1 = jfin - 1;
    for (i = 0; i <= 162 + 1; i++) {
        for (k = 0; k <= 162 + 1; k++) {
            phi1[i][k] = 0.;
            phi2[i][k] = 0.;
        }
    }
    for (i = ibeg; i <= ifin; i++) {
        iglob = i;
        for (j = jbeg; j <= jfin; j++) {
            jglob = j;
            k = ki1;
            phi1[i][j] = 0.40000000000000002
                * (u[i][j][k][4]
                    - 0.5
                        * (((u[i][j][k][1]) * (u[i][j][k][1])) + ((u[i][j][k][2]) * (u[i][j][k][2]))
                            + ((u[i][j][k][3]) * (u[i][j][k][3])))
                        / u[i][j][k][0]);
            k = ki2;
            phi2[i][j] = 0.40000000000000002
                * (u[i][j][k][4]
                    - 0.5
                        * (((u[i][j][k][1]) * (u[i][j][k][1])) + ((u[i][j][k][2]) * (u[i][j][k][2]))
                            + ((u[i][j][k][3]) * (u[i][j][k][3])))
                        / u[i][j][k][0]);
        }
    }
    frc1 = 0.;
    for (i = ibeg; i <= ifin1; i++) {
        for (j = jbeg; j <= jfin1; j++) {
            frc1 = frc1
                + (phi1[i][j] + phi1[i + 1][j] + phi1[i][j + 1] + phi1[i + 1][j + 1] + phi2[i][j]
                    + phi2[i + 1][j] + phi2[i][j + 1] + phi2[i + 1][j + 1]);
        }
    }
    frc1 = dxi * deta * frc1;
    for (i = 0; i <= 162 + 1; i++) {
        for (k = 0; k <= 162 + 1; k++) {
            phi1[i][k] = 0.;
            phi2[i][k] = 0.;
        }
    }
    jglob = jbeg;
    if (jglob == ji1) {
        for (i = ibeg; i <= ifin; i++) {
            iglob = i;
            for (k = ki1; k <= ki2; k++) {
                phi1[i][k] = 0.40000000000000002
                    * (u[i][jbeg][k][4]
                        - 0.5
                            * (((u[i][jbeg][k][1]) * (u[i][jbeg][k][1]))
                                + ((u[i][jbeg][k][2]) * (u[i][jbeg][k][2]))
                                + ((u[i][jbeg][k][3]) * (u[i][jbeg][k][3])))
                            / u[i][jbeg][k][0]);
            }
        }
    }
    jglob = jfin;
    if (jglob == ji2) {
        for (i = ibeg; i <= ifin; i++) {
            iglob = i;
            for (k = ki1; k <= ki2; k++) {
                phi2[i][k] = 0.40000000000000002
                    * (u[i][jfin][k][4]
                        - 0.5
                            * (((u[i][jfin][k][1]) * (u[i][jfin][k][1]))
                                + ((u[i][jfin][k][2]) * (u[i][jfin][k][2]))
                                + ((u[i][jfin][k][3]) * (u[i][jfin][k][3])))
                            / u[i][jfin][k][0]);
            }
        }
    }
    frc2 = 0.;
    for (i = ibeg; i <= ifin1; i++) {
        for (k = ki1; k <= ki2 - 1; k++) {
            frc2 = frc2
                + (phi1[i][k] + phi1[i + 1][k] + phi1[i][k + 1] + phi1[i + 1][k + 1] + phi2[i][k]
                    + phi2[i + 1][k] + phi2[i][k + 1] + phi2[i + 1][k + 1]);
        }
    }
    frc2 = dxi * dzeta * frc2;
    for (i = 0; i <= 162 + 1; i++) {
        for (k = 0; k <= 162 + 1; k++) {
            phi1[i][k] = 0.;
            phi2[i][k] = 0.;
        }
    }
    iglob = ibeg;
    if (iglob == ii1) {
        for (j = jbeg; j <= jfin; j++) {
            jglob = j;
            for (k = ki1; k <= ki2; k++) {
                phi1[j][k] = 0.40000000000000002
                    * (u[ibeg][j][k][4]
                        - 0.5
                            * (((u[ibeg][j][k][1]) * (u[ibeg][j][k][1]))
                                + ((u[ibeg][j][k][2]) * (u[ibeg][j][k][2]))
                                + ((u[ibeg][j][k][3]) * (u[ibeg][j][k][3])))
                            / u[ibeg][j][k][0]);
            }
        }
    }
    iglob = ifin;
    if (iglob == ii2) {
        for (j = jbeg; j <= jfin; j++) {
            jglob = j;
            for (k = ki1; k <= ki2; k++) {
                phi2[j][k] = 0.40000000000000002
                    * (u[ifin][j][k][4]
                        - 0.5
                            * (((u[ifin][j][k][1]) * (u[ifin][j][k][1]))
                                + ((u[ifin][j][k][2]) * (u[ifin][j][k][2]))
                                + ((u[ifin][j][k][3]) * (u[ifin][j][k][3])))
                            / u[ifin][j][k][0]);
            }
        }
    }
    frc3 = 0.;
    for (j = jbeg; j <= jfin1; j++) {
        for (k = ki1; k <= ki2 - 1; k++) {
            frc3 = frc3
                + (phi1[j][k] + phi1[j + 1][k] + phi1[j][k + 1] + phi1[j + 1][k + 1] + phi2[j][k]
                    + phi2[j + 1][k] + phi2[j][k + 1] + phi2[j + 1][k + 1]);
        }
    }
    frc3 = deta * dzeta * frc3;
    frc = 0.25 * (frc1 + frc2 + frc3);
}
/*Function: read_input, ID: 11*/
static void read_input()
{
    /*read_input:11*/
    /*CompoundStmt:10751*/
    FILE* fp;
    printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version - LU Benchmark\n\n");
    fp = fopen("inputlu.data", "r");
    if (fp != ((void*)0)) {
        printf(" Reading from input file inputlu.data\n");
        while (fgetc(fp) != '\n')
            ;
        while (fgetc(fp) != '\n')
            ;
        fscanf(fp, "%d%d", &ipr, &inorm);
        while (fgetc(fp) != '\n')
            ;
        while (fgetc(fp) != '\n')
            ;
        while (fgetc(fp) != '\n')
            ;
        fscanf(fp, "%d", &itmax);
        while (fgetc(fp) != '\n')
            ;
        while (fgetc(fp) != '\n')
            ;
        while (fgetc(fp) != '\n')
            ;
        fscanf(fp, "%lf", &dt);
        while (fgetc(fp) != '\n')
            ;
        while (fgetc(fp) != '\n')
            ;
        while (fgetc(fp) != '\n')
            ;
        fscanf(fp, "%lf", &omega);
        while (fgetc(fp) != '\n')
            ;
        while (fgetc(fp) != '\n')
            ;
        while (fgetc(fp) != '\n')
            ;
        fscanf(fp, "%lf%lf%lf%lf%lf", &tolrsd[0], &tolrsd[1], &tolrsd[2], &tolrsd[3], &tolrsd[4]);
        while (fgetc(fp) != '\n')
            ;
        while (fgetc(fp) != '\n')
            ;
        while (fgetc(fp) != '\n')
            ;
        fscanf(fp, "%d%d%d", &nx0, &ny0, &nz0);
        while (fgetc(fp) != '\n')
            ;
        fclose(fp);
    } else {
        ipr = 1;
        inorm = 250;
        itmax = 250;
        dt = 2.;
        omega = 1.2;
        tolrsd[0] = 1.0E-8;
        tolrsd[1] = 1.0E-8;
        tolrsd[2] = 1.0E-8;
        tolrsd[3] = 1.0E-8;
        tolrsd[4] = 1.0E-8;
        nx0 = 162;
        ny0 = 162;
        nz0 = 162;
    }
    if (nx0 < 4 || ny0 < 4 || nz0 < 4) {
        printf("     PROBLEM SIZE IS TOO SMALL - \n     SET EACH OF NX, NY AND NZ AT LEAST EQUAL "
               "TO 5\n");
        exit(1);
    }
    if (nx0 > 162 || ny0 > 162 || nz0 > 162) {
        printf("     PROBLEM SIZE IS TOO LARGE - \n     NX, NY AND NZ SHOULD BE EQUAL TO \n     "
               "ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY\n");
        exit(1);
    }
    printf(" Size: %3dx%3dx%3d\n", nx0, ny0, nz0);
    printf(" Iterations: %3d\n", itmax);
}
/*Function: rhs, ID: 12*/
static void rhs()
{
    /*rhs:12*/
    /*CompoundStmt:10895*/
    if (affinMaskRes) {
        myDARTSRuntime->run(
            launch<TP10896>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
    }
}
/*Function: setbv, ID: 13*/
static void setbv()
{
    /*setbv:13*/
    /*CompoundStmt:13296*/
    if (affinMaskRes) {
        myDARTSRuntime->run(
            launch<TP13297>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
    }
}
/*Function: setcoeff, ID: 14*/
static void setcoeff()
{
    /*setcoeff:14*/
    /*CompoundStmt:13524*/
    dxi = 1. / (nx0 - 1);
    deta = 1. / (ny0 - 1);
    dzeta = 1. / (nz0 - 1);
    tx1 = 1. / (dxi * dxi);
    tx2 = 1. / (2. * dxi);
    tx3 = 1. / dxi;
    ty1 = 1. / (deta * deta);
    ty2 = 1. / (2. * deta);
    ty3 = 1. / deta;
    tz1 = 1. / (dzeta * dzeta);
    tz2 = 1. / (2. * dzeta);
    tz3 = 1. / dzeta;
    ii1 = 1;
    ii2 = nx0 - 2;
    ji1 = 1;
    ji2 = ny0 - 3;
    ki1 = 2;
    ki2 = nz0 - 2;
    dx1 = 0.75;
    dx2 = dx1;
    dx3 = dx1;
    dx4 = dx1;
    dx5 = dx1;
    dy1 = 0.75;
    dy2 = dy1;
    dy3 = dy1;
    dy4 = dy1;
    dy5 = dy1;
    dz1 = 1.;
    dz2 = dz1;
    dz3 = dz1;
    dz4 = dz1;
    dz5 = dz1;
    dssp = dz1 / 4.;
    ce[0][0] = 2.;
    ce[0][1] = 0.;
    ce[0][2] = 0.;
    ce[0][3] = 4.;
    ce[0][4] = 5.;
    ce[0][5] = 3.;
    ce[0][6] = 0.5;
    ce[0][7] = 0.02;
    ce[0][8] = 0.01;
    ce[0][9] = 0.029999999999999999;
    ce[0][10] = 0.5;
    ce[0][11] = 0.40000000000000002;
    ce[0][12] = 0.29999999999999999;
    ce[1][0] = 1.;
    ce[1][1] = 0.;
    ce[1][2] = 0.;
    ce[1][3] = 0.;
    ce[1][4] = 1.;
    ce[1][5] = 2.;
    ce[1][6] = 3.;
    ce[1][7] = 0.01;
    ce[1][8] = 0.029999999999999999;
    ce[1][9] = 0.02;
    ce[1][10] = 0.40000000000000002;
    ce[1][11] = 0.29999999999999999;
    ce[1][12] = 0.5;
    ce[2][0] = 2.;
    ce[2][1] = 2.;
    ce[2][2] = 0.;
    ce[2][3] = 0.;
    ce[2][4] = 0.;
    ce[2][5] = 2.;
    ce[2][6] = 3.;
    ce[2][7] = 0.040000000000000001;
    ce[2][8] = 0.029999999999999999;
    ce[2][9] = 0.050000000000000003;
    ce[2][10] = 0.29999999999999999;
    ce[2][11] = 0.5;
    ce[2][12] = 0.40000000000000002;
    ce[3][0] = 2.;
    ce[3][1] = 2.;
    ce[3][2] = 0.;
    ce[3][3] = 0.;
    ce[3][4] = 0.;
    ce[3][5] = 2.;
    ce[3][6] = 3.;
    ce[3][7] = 0.029999999999999999;
    ce[3][8] = 0.050000000000000003;
    ce[3][9] = 0.040000000000000001;
    ce[3][10] = 0.20000000000000001;
    ce[3][11] = 0.10000000000000001;
    ce[3][12] = 0.29999999999999999;
    ce[4][0] = 5.;
    ce[4][1] = 4.;
    ce[4][2] = 3.;
    ce[4][3] = 2.;
    ce[4][4] = 0.10000000000000001;
    ce[4][5] = 0.40000000000000002;
    ce[4][6] = 0.29999999999999999;
    ce[4][7] = 0.050000000000000003;
    ce[4][8] = 0.040000000000000001;
    ce[4][9] = 0.029999999999999999;
    ce[4][10] = 0.10000000000000001;
    ce[4][11] = 0.29999999999999999;
    ce[4][12] = 0.20000000000000001;
}
/*Function: setiv, ID: 15*/
static void setiv()
{
    /*setiv:15*/
    /*CompoundStmt:13863*/
    if (affinMaskRes) {
        myDARTSRuntime->run(
            launch<TP13864>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
    }
}
/*Function: ssor, ID: 16*/
static void ssor()
{
    /*ssor:16*/
    /*CompoundStmt:13986*/
    int i, j, k, m;
    int istep;
    double tmp;
    double delunm[5], tv[162][162][5];
    tmp = 1. / (omega * (2. - omega));
    if (affinMaskRes) {
        myDARTSRuntime->run(
            launch<TP13997>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
    }
    rhs();
    l2norm(nx0, ny0, nz0, ist, iend, jst, jend, rsd, rsdnm);
    timer_clear(1);
    timer_start(1);
    for (istep = 1; istep <= itmax; istep++) {
        /*CompoundStmt:14090*/
        if (istep % 20 == 0 || istep == itmax || istep == 1) {
            //#pragma omp master
            printf(" Time step %4d\n", istep);
        }
        if (affinMaskRes) {
myDARTSRuntime->run(launch<TP14096>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet, (double*)&((tmp)), (double *)((tv)), (162][162][5)));
        }
        if (istep % inorm == 0) {
            l2norm(nx0, ny0, nz0, ist, iend, jst, jend, rsd, delunm);
        }
        rhs();
        if ((istep % inorm == 0) || (istep == itmax)) {
            l2norm(nx0, ny0, nz0, ist, iend, jst, jend, rsd, rsdnm);
        }
        if ((rsdnm[0] < tolrsd[0]) && (rsdnm[1] < tolrsd[1]) && (rsdnm[2] < tolrsd[2])
            && (rsdnm[3] < tolrsd[3]) && (rsdnm[4] < tolrsd[4])) {
            exit(1);
        }
    }
    timer_stop(1);
    maxtime = timer_read(1);
}
/*Function: verify, ID: 17*/
static void verify(double xcr[5], double xce[5], double xci, char* class, boolean* verified)
{
    /*verify:17*/
    /*CompoundStmt:14247*/
    double xcrref[5], xceref[5], xciref, xcrdif[5], xcedif[5], xcidif, epsilon, dtref;
    int m;
    epsilon = 1.0E-8;
    *class = 'U';
    *verified = 1;
    for (m = 0; m < 5; m++) {
        xcrref[m] = 1.;
        xceref[m] = 1.;
    }
    xciref = 1.;
    if (nx0 == 12 && ny0 == 12 && nz0 == 12 && itmax == 50) {
        *class = 'S';
        dtref = 0.5;
        xcrref[0] = 0.016196343210976703;
        xcrref[1] = 0.002197674516482132;
        xcrref[2] = 0.0015179927653399185;
        xcrref[3] = 0.0015029584435994323;
        xcrref[4] = 0.034264073155896461;
        xceref[0] = 6.4223319957960922E-4;
        xceref[1] = 8.4144342047347924E-5;
        xceref[2] = 5.8588269616485187E-5;
        xceref[3] = 5.847422259515735E-5;
        xceref[4] = 0.0013103347914111294;
        xciref = 7.8418928865937083;
    } else if (nx0 == 33 && ny0 == 33 && nz0 == 33 && itmax == 300) {
        *class = 'W';
        dtref = 0.0015;
        xcrref[0] = 12.36511638192;
        xcrref[1] = 1.317228477799;
        xcrref[2] = 2.5501207130950001;
        xcrref[3] = 2.3261877502520001;
        xcrref[4] = 28.26799444189;
        xceref[0] = 0.48678771442160002;
        xceref[1] = 0.050646528809819999;
        xceref[2] = 0.092818181019600002;
        xceref[3] = 0.085701265427330003;
        xceref[4] = 1.084277417792;
        xciref = 11.61399311023;
    } else if (nx0 == 64 && ny0 == 64 && nz0 == 64 && itmax == 250) {
        *class = 'A';
        dtref = 2.;
        xcrref[0] = 779.02107606689367;
        xcrref[1] = 63.40276525969287;
        xcrref[2] = 194.99249727292479;
        xcrref[3] = 178.45301160418538;
        xcrref[4] = 1838.4760349464248;
        xceref[0] = 29.964085685471943;
        xceref[1] = 2.819457636500335;
        xceref[2] = 7.3473412698774743;
        xceref[3] = 6.7139225687777051;
        xceref[4] = 70.715315688392579;
        xciref = 26.030925604886278;
    } else if (nx0 == 102 && ny0 == 102 && nz0 == 102 && itmax == 250) {
        *class = 'B';
        dtref = 2.;
        xcrref[0] = 3553.2672969982737;
        xcrref[1] = 262.14750795310692;
        xcrref[2] = 883.3372185095219;
        xcrref[3] = 778.12774739425265;
        xcrref[4] = 7308.7969592545314;
        xceref[0] = 114.01176380212709;
        xceref[1] = 8.1098963655421574;
        xceref[2] = 28.480597317698308;
        xceref[3] = 25.905394567832939;
        xceref[4] = 260.54907504857414;
        xciref = 47.887162703308228;
    } else if (nx0 == 162 && ny0 == 162 && nz0 == 162 && itmax == 250) {
        *class = 'C';
        dtref = 2.;
        xcrref[0] = 10376.698032353785;
        xcrref[1] = 892.21245880100855;
        xcrref[2] = 2562.3881458266087;
        xcrref[3] = 2191.9434385783143;
        xcrref[4] = 17807.805726106119;
        xceref[0] = 215.98639971694928;
        xceref[1] = 15.57895592398636;
        xceref[2] = 54.131886307720777;
        xceref[3] = 48.226264315404542;
        xceref[4] = 455.90291004325036;
        xciref = 66.64045535721813;
    } else {
        *verified = 0;
    }
    for (m = 0; m < 5; m++) {
        xcrdif[m] = fabs((xcr[m] - xcrref[m]) / xcrref[m]);
        xcedif[m] = fabs((xce[m] - xceref[m]) / xceref[m]);
    }
    xcidif = fabs((xci - xciref) / xciref);
    if (*class != 'U') {
        printf("\n Verification being performed for class %1c\n", *class);
        printf(" Accuracy setting for epsilon = %20.13e\n", epsilon);
        if (fabs(dt - dtref) > epsilon) {
            *verified = 0;
            *class = 'U';
            printf(" DT does not match the reference value of %15.8e\n", dtref);
        }
    } else {
        printf(" Unknown class\n");
    }
    if (*class != 'U') {
        printf(" Comparison of RMS-norms of residual\n");
    } else {
        printf(" RMS-norms of residual\n");
    }
    for (m = 0; m < 5; m++) {
        if (*class == 'U') {
            printf("          %2d  %20.13e\n", m, xcr[m]);
        } else if (xcrdif[m] > epsilon) {
            *verified = 0;
            printf(" FAILURE: %2d  %20.13e%20.13e%20.13e\n", m, xcr[m], xcrref[m], xcrdif[m]);
        } else {
            printf("          %2d  %20.13e%20.13e%20.13e\n", m, xcr[m], xcrref[m], xcrdif[m]);
        }
    }
    if (*class != 'U') {
        printf(" Comparison of RMS-norms of solution error\n");
    } else {
        printf(" RMS-norms of solution error\n");
    }
    for (m = 0; m < 5; m++) {
        if (*class == 'U') {
            printf("          %2d  %20.13e\n", m, xce[m]);
        } else if (xcedif[m] > epsilon) {
            *verified = 0;
            printf(" FAILURE: %2d  %20.13e%20.13e%20.13e\n", m, xce[m], xceref[m], xcedif[m]);
        } else {
            printf("          %2d  %20.13e%20.13e%20.13e\n", m, xce[m], xceref[m], xcedif[m]);
        }
    }
    if (*class != 'U') {
        printf(" Comparison of surface integral\n");
    } else {
        printf(" Surface integral\n");
    }
    if (*class == 'U') {
        printf("              %20.13e\n", xci);
    } else if (xcidif > epsilon) {
        *verified = 0;
        printf(" FAILURE:     %20.13e%20.13e%20.13e\n", xci, xciref, xcidif);
    } else {
        printf("              %20.13e%20.13e%20.13e\n", xci, xciref, xcidif);
    }
    if (*class == 'U') {
        printf(" No reference values provided\n");
        printf(" No verification performed\n");
    } else if (*verified) {
        printf(" Verification Successful\n");
    } else {
        printf(" Verification failed\n");
    }
}
/*TP1: TP_blts*/
void TP1::_checkInCodelets205::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Allocate each array var in the codelet */
this->inputsTPParent->d_outer1_size = 162][162][5][5;
this->inputsTPParent->d_darts1[this->getID()] = (double *)malloc(sizeof(double ) * this->inputsTPParent->d_outer1_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->ldx_outer1_size = 162][5][5;
this->inputsTPParent->ldx_darts1[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->ldx_outer1_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->ldy_outer1_size = 162][5][5;
this->inputsTPParent->ldy_darts1[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->ldy_outer1_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->ldz_outer1_size = 162][5][5;
this->inputsTPParent->ldz_darts1[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->ldz_outer1_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->v_outer1_size = 163][163][5;
this->inputsTPParent->v_darts1[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->v_outer1_size);

/*printing node 205: BinaryOperator*/
162 / 2 * 2 + 1;

/*printing node 208: BinaryOperator*/
162 / 2 * 2 + 1;
/*Signaling next codelet from last stmt in the codelet*/
/*Find and signal the next codelet*/

myTP->controlTPParent->nextCodeletsblts[this->getID()]->decDep();
}
TP1::TP1(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
    darts::Codelet* in_mainSyncCodelet, TP1** in_ptrToThisFunctionTP, int in_nx, int in_ny,
    int in_nz, int in_k, double in_omega, double(*) in_v, double(*) in_ldz, double(*) in_ldy,
    double(*) in_ldx, double(*) in_d, int in_ist, int in_iend, int in_jst, int in_jend, int in_nx0,
    int in_ny0)
    : ompTP(in_numThreads, in_mainCodeletID)
    , ptrToThisFunctionTP(in_ptrToThisFunctionTP)
    , inputsTPParent(this)
    , controlTPParent(this)
    , nextCodeletsblts(new Codelet*[in_numThreads])
    , nextSyncCodeletsblts(new Codelet*[in_numThreads])
    , nx_darts1(new int[this->numThreads])
    , ny_darts1(new int[this->numThreads])
    , nz_darts1(new int[this->numThreads])
    , k_darts1(new int[this->numThreads])
    , omega_darts1(new double[this->numThreads])
    , v_darts1(new double(*) *[this->numThreads])
    , ldz_darts1(new double(*) *[this->numThreads])
    , ldy_darts1(new double(*) *[this->numThreads])
    , ldx_darts1(new double(*) *[this->numThreads])
    , d_darts1(new double(*) *[this->numThreads])
    , ist_darts1(new int[this->numThreads])
    , iend_darts1(new int[this->numThreads])
    , jst_darts1(new int[this->numThreads])
    , jend_darts1(new int[this->numThreads])
    , nx0_darts1(new int[this->numThreads])
    , ny0_darts1(new int[this->numThreads])
    , checkInCodelets205(new _checkInCodelets205[this->numThreads])
{
    _checkInCodelets205* checkInCodelets205Ptr = (this->checkInCodelets205);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets205);
#endif
    for (size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets205Ptr) = _checkInCodelets205(2, 1, this, codeletCounter);
#else
        (*checkInCodelets205Ptr) = _checkInCodelets205(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets205Ptr).decDep();
        checkInCodelets205Ptr++;
    }
    if (this->numThreads == 1) {
        this->nextCodeletsblts[0] = in_mainNextCodelet;
        this->nextSyncCodeletsblts[0] = in_mainSyncCodelet;
        this->nx_darts1[0] = in_nx;
        this->ny_darts1[0] = in_ny;
        this->nz_darts1[0] = in_nz;
        this->k_darts1[0] = in_k;
        this->omega_darts1[0] = in_omega;
        this->v_darts1[0] = in_v;
        this->ldz_darts1[0] = in_ldz;
        this->ldy_darts1[0] = in_ldy;
        this->ldx_darts1[0] = in_ldx;
        this->d_darts1[0] = in_d;
        this->ist_darts1[0] = in_ist;
        this->iend_darts1[0] = in_iend;
        this->jst_darts1[0] = in_jst;
        this->jend_darts1[0] = in_jend;
        this->nx0_darts1[0] = in_nx0;
        this->ny0_darts1[0] = in_ny0;
        this->availableCodelets[0] = 1;
    } else {
        this->nx_darts1[this->mainCodeletID] = in_nx;
        this->ny_darts1[this->mainCodeletID] = in_ny;
        this->nz_darts1[this->mainCodeletID] = in_nz;
        this->k_darts1[this->mainCodeletID] = in_k;
        this->omega_darts1[this->mainCodeletID] = in_omega;
        this->v_darts1[this->mainCodeletID] = in_v;
        this->ldz_darts1[this->mainCodeletID] = in_ldz;
        this->ldy_darts1[this->mainCodeletID] = in_ldy;
        this->ldx_darts1[this->mainCodeletID] = in_ldx;
        this->d_darts1[this->mainCodeletID] = in_d;
        this->ist_darts1[this->mainCodeletID] = in_ist;
        this->iend_darts1[this->mainCodeletID] = in_iend;
        this->jst_darts1[this->mainCodeletID] = in_jst;
        this->jend_darts1[this->mainCodeletID] = in_jend;
        this->nx0_darts1[this->mainCodeletID] = in_nx0;
        this->ny0_darts1[this->mainCodeletID] = in_ny0;
        this->nextCodeletsblts[in_mainCodeletID] = in_mainNextCodelet;
        this->nextSyncCodeletsblts[in_mainCodeletID] = in_mainSyncCodelet;
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[this->mainCodeletID].decDep();
#else
        this->availableCodelets[this->mainCodeletID] = 1;
#endif
        *(this->ptrToThisFunctionTP) = this;
    }
}
TP1::~TP1()
{
    delete[] checkInCodelets205;
    delete[] nextCodeletsblts;
    delete[] nextSyncCodeletsblts;
    delete[] nx_darts1;
    delete[] ny_darts1;
    delete[] nz_darts1;
    delete[] k_darts1;
    delete[] omega_darts1;
    delete[] v_darts1;
    delete[] ldz_darts1;
    delete[] ldy_darts1;
    delete[] ldx_darts1;
    delete[] d_darts1;
    delete[] ist_darts1;
    delete[] iend_darts1;
    delete[] jst_darts1;
    delete[] jend_darts1;
    delete[] nx0_darts1;
    delete[] ny0_darts1;
}
void TP1::setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, double(*) in_v,
    double(*) in_ldz, double(*) in_ldy, double(*) in_ldx, double(*) in_d, int in_ist, int in_iend,
    int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID)
{
    this->nx_darts1[codeletID] = in_nx;
    this->ny_darts1[codeletID] = in_ny;
    this->nz_darts1[codeletID] = in_nz;
    this->k_darts1[codeletID] = in_k;
    this->omega_darts1[codeletID] = in_omega;
    this->v_darts1[codeletID] = in_v;
    this->ldz_darts1[codeletID] = in_ldz;
    this->ldy_darts1[codeletID] = in_ldy;
    this->ldx_darts1[codeletID] = in_ldx;
    this->d_darts1[codeletID] = in_d;
    this->ist_darts1[codeletID] = in_ist;
    this->iend_darts1[codeletID] = in_iend;
    this->jst_darts1[codeletID] = in_jst;
    this->jend_darts1[codeletID] = in_jend;
    this->nx0_darts1[codeletID] = in_nx0;
    this->ny0_darts1[codeletID] = in_ny0;
}
/*TP2: TP_buts*/
void TP2::_checkInCodelets211::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Allocate each array var in the codelet */
this->inputsTPParent->d_outer2_size = 162][162][5][5;
this->inputsTPParent->d_darts2[this->getID()] = (double *)malloc(sizeof(double ) * this->inputsTPParent->d_outer2_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->tv_outer2_size = 162][5;
this->inputsTPParent->tv_darts2[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->tv_outer2_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->udx_outer2_size = 162][5][5;
this->inputsTPParent->udx_darts2[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->udx_outer2_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->udy_outer2_size = 162][5][5;
this->inputsTPParent->udy_darts2[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->udy_outer2_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->udz_outer2_size = 162][5][5;
this->inputsTPParent->udz_darts2[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->udz_outer2_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->v_outer2_size = 163][163][5;
this->inputsTPParent->v_darts2[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->v_outer2_size);

/*printing node 211: BinaryOperator*/
162 / 2 * 2 + 1;

/*printing node 214: BinaryOperator*/
162 / 2 * 2 + 1;
/*Signaling next codelet from last stmt in the codelet*/
/*Find and signal the next codelet*/

myTP->controlTPParent->nextCodeletsbuts[this->getID()]->decDep();
}
TP2::TP2(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
    darts::Codelet* in_mainSyncCodelet, TP2** in_ptrToThisFunctionTP, int in_nx, int in_ny,
    int in_nz, int in_k, double in_omega, double(*) in_v, double(*) in_tv, double(*) in_d,
    double(*) in_udx, double(*) in_udy, double(*) in_udz, int in_ist, int in_iend, int in_jst,
    int in_jend, int in_nx0, int in_ny0)
    : ompTP(in_numThreads, in_mainCodeletID)
    , ptrToThisFunctionTP(in_ptrToThisFunctionTP)
    , inputsTPParent(this)
    , controlTPParent(this)
    , nextCodeletsbuts(new Codelet*[in_numThreads])
    , nextSyncCodeletsbuts(new Codelet*[in_numThreads])
    , nx_darts2(new int[this->numThreads])
    , ny_darts2(new int[this->numThreads])
    , nz_darts2(new int[this->numThreads])
    , k_darts2(new int[this->numThreads])
    , omega_darts2(new double[this->numThreads])
    , v_darts2(new double(*) *[this->numThreads])
    , tv_darts2(new double(*) *[this->numThreads])
    , d_darts2(new double(*) *[this->numThreads])
    , udx_darts2(new double(*) *[this->numThreads])
    , udy_darts2(new double(*) *[this->numThreads])
    , udz_darts2(new double(*) *[this->numThreads])
    , ist_darts2(new int[this->numThreads])
    , iend_darts2(new int[this->numThreads])
    , jst_darts2(new int[this->numThreads])
    , jend_darts2(new int[this->numThreads])
    , nx0_darts2(new int[this->numThreads])
    , ny0_darts2(new int[this->numThreads])
    , checkInCodelets211(new _checkInCodelets211[this->numThreads])
{
    _checkInCodelets211* checkInCodelets211Ptr = (this->checkInCodelets211);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets211);
#endif
    for (size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets211Ptr) = _checkInCodelets211(2, 1, this, codeletCounter);
#else
        (*checkInCodelets211Ptr) = _checkInCodelets211(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets211Ptr).decDep();
        checkInCodelets211Ptr++;
    }
    if (this->numThreads == 1) {
        this->nextCodeletsbuts[0] = in_mainNextCodelet;
        this->nextSyncCodeletsbuts[0] = in_mainSyncCodelet;
        this->nx_darts2[0] = in_nx;
        this->ny_darts2[0] = in_ny;
        this->nz_darts2[0] = in_nz;
        this->k_darts2[0] = in_k;
        this->omega_darts2[0] = in_omega;
        this->v_darts2[0] = in_v;
        this->tv_darts2[0] = in_tv;
        this->d_darts2[0] = in_d;
        this->udx_darts2[0] = in_udx;
        this->udy_darts2[0] = in_udy;
        this->udz_darts2[0] = in_udz;
        this->ist_darts2[0] = in_ist;
        this->iend_darts2[0] = in_iend;
        this->jst_darts2[0] = in_jst;
        this->jend_darts2[0] = in_jend;
        this->nx0_darts2[0] = in_nx0;
        this->ny0_darts2[0] = in_ny0;
        this->availableCodelets[0] = 1;
    } else {
        this->nx_darts2[this->mainCodeletID] = in_nx;
        this->ny_darts2[this->mainCodeletID] = in_ny;
        this->nz_darts2[this->mainCodeletID] = in_nz;
        this->k_darts2[this->mainCodeletID] = in_k;
        this->omega_darts2[this->mainCodeletID] = in_omega;
        this->v_darts2[this->mainCodeletID] = in_v;
        this->tv_darts2[this->mainCodeletID] = in_tv;
        this->d_darts2[this->mainCodeletID] = in_d;
        this->udx_darts2[this->mainCodeletID] = in_udx;
        this->udy_darts2[this->mainCodeletID] = in_udy;
        this->udz_darts2[this->mainCodeletID] = in_udz;
        this->ist_darts2[this->mainCodeletID] = in_ist;
        this->iend_darts2[this->mainCodeletID] = in_iend;
        this->jst_darts2[this->mainCodeletID] = in_jst;
        this->jend_darts2[this->mainCodeletID] = in_jend;
        this->nx0_darts2[this->mainCodeletID] = in_nx0;
        this->ny0_darts2[this->mainCodeletID] = in_ny0;
        this->nextCodeletsbuts[in_mainCodeletID] = in_mainNextCodelet;
        this->nextSyncCodeletsbuts[in_mainCodeletID] = in_mainSyncCodelet;
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[this->mainCodeletID].decDep();
#else
        this->availableCodelets[this->mainCodeletID] = 1;
#endif
        *(this->ptrToThisFunctionTP) = this;
    }
}
TP2::~TP2()
{
    delete[] checkInCodelets211;
    delete[] nextCodeletsbuts;
    delete[] nextSyncCodeletsbuts;
    delete[] nx_darts2;
    delete[] ny_darts2;
    delete[] nz_darts2;
    delete[] k_darts2;
    delete[] omega_darts2;
    delete[] v_darts2;
    delete[] tv_darts2;
    delete[] d_darts2;
    delete[] udx_darts2;
    delete[] udy_darts2;
    delete[] udz_darts2;
    delete[] ist_darts2;
    delete[] iend_darts2;
    delete[] jst_darts2;
    delete[] jend_darts2;
    delete[] nx0_darts2;
    delete[] ny0_darts2;
}
void TP2::setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, double(*) in_v,
    double(*) in_tv, double(*) in_d, double(*) in_udx, double(*) in_udy, double(*) in_udz,
    int in_ist, int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID)
{
    this->nx_darts2[codeletID] = in_nx;
    this->ny_darts2[codeletID] = in_ny;
    this->nz_darts2[codeletID] = in_nz;
    this->k_darts2[codeletID] = in_k;
    this->omega_darts2[codeletID] = in_omega;
    this->v_darts2[codeletID] = in_v;
    this->tv_darts2[codeletID] = in_tv;
    this->d_darts2[codeletID] = in_d;
    this->udx_darts2[codeletID] = in_udx;
    this->udy_darts2[codeletID] = in_udy;
    this->udz_darts2[codeletID] = in_udz;
    this->ist_darts2[codeletID] = in_ist;
    this->iend_darts2[codeletID] = in_iend;
    this->jst_darts2[codeletID] = in_jst;
    this->jend_darts2[codeletID] = in_jend;
    this->nx0_darts2[codeletID] = in_nx0;
    this->ny0_darts2[codeletID] = in_ny0;
}
/*TP234: OMPParallelDirective*/
void TP234::_barrierCodelets234::fire(void)
{
    TP234* myTP = static_cast<TP234*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP234::_checkInCodelets236::fire(void)
{
    if (this->getID() == 0) {
        /*Init the vars for this region*/
        /*Initialize the vars of the inlined region*/
        this->inputsTPParent->nthreads_darts236
            = (this->inputsTPParent->nthreads_darts234) /*OMP_SHARED - VAR INLINED*/;

        /*printing node 237: BinaryOperator*/
        (*(this->inputsTPParent->nthreads_darts236)) = omp_get_num_threads();
        /*Signaling next codelet from last stmt in the codelet*/
        /*Find and signal the next codelet*/
        myTP->controlTPParent->TPParent->barrierCodelets234[0].decDep();
    } else {
        /*Find and signal the next codelet*/
        myTP->TPParent->barrierCodelets234[0].decDep();
    }
}
TP234::TP234(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet, int* in_nthreads)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , nthreads_darts234(in_nthreads) /*OMP_SHARED - INPUT*/
    , TP236_alreadyLaunched(0)
    , barrierCodelets234(new _barrierCodelets234[1])
    , checkInCodelets236(new _checkInCodelets236[this->numThreads])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets234[0] = _barrierCodelets234(ompNumThreads, ompNumThreads, this, 0);
    _checkInCodelets236* checkInCodelets236Ptr = (this->checkInCodelets236);
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets236Ptr) = _checkInCodelets236(1, 1, this, codeletCounter);
        (*checkInCodelets236Ptr).decDep();
        checkInCodelets236Ptr++;
    }
}
TP234::~TP234()
{
    delete[] barrierCodelets234;
    delete[] checkInCodelets236;
}
/*TP1: TP_blts*/
void TP1::_checkInCodelets205::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Allocate each array var in the codelet */
this->inputsTPParent->d_outer1_size = 162][162][5][5;
this->inputsTPParent->d_darts1[this->getID()] = (double *)malloc(sizeof(double ) * this->inputsTPParent->d_outer1_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->ldx_outer1_size = 162][5][5;
this->inputsTPParent->ldx_darts1[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->ldx_outer1_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->ldy_outer1_size = 162][5][5;
this->inputsTPParent->ldy_darts1[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->ldy_outer1_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->ldz_outer1_size = 162][5][5;
this->inputsTPParent->ldz_darts1[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->ldz_outer1_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->v_outer1_size = 163][163][5;
this->inputsTPParent->v_darts1[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->v_outer1_size);

/*printing node 205: BinaryOperator*/
162 / 2 * 2 + 1;
162 / 2 * 2 + 1;

/*printing node 208: BinaryOperator*/
162 / 2 * 2 + 1;
162 / 2 * 2 + 1;
/*Signaling next codelet from last stmt in the codelet*/
/*Find and signal the next codelet*/

myTP->controlTPParent->nextCodeletsblts[this->getID()]->decDep();
}
TP1::TP1(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
    darts::Codelet* in_mainSyncCodelet, TP1** in_ptrToThisFunctionTP, int in_nx, int in_ny,
    int in_nz, int in_k, double in_omega, double(*) in_v, double(*) in_ldz, double(*) in_ldy,
    double(*) in_ldx, double(*) in_d, int in_ist, int in_iend, int in_jst, int in_jend, int in_nx0,
    int in_ny0)
    : ompTP(in_numThreads, in_mainCodeletID)
    , ptrToThisFunctionTP(in_ptrToThisFunctionTP)
    , inputsTPParent(this)
    , controlTPParent(this)
    , nextCodeletsblts(new Codelet*[in_numThreads])
    , nextSyncCodeletsblts(new Codelet*[in_numThreads])
    , nx_darts1(new int[this->numThreads])
    , ny_darts1(new int[this->numThreads])
    , nz_darts1(new int[this->numThreads])
    , k_darts1(new int[this->numThreads])
    , omega_darts1(new double[this->numThreads])
    , v_darts1(new double(*) *[this->numThreads])
    , ldz_darts1(new double(*) *[this->numThreads])
    , ldy_darts1(new double(*) *[this->numThreads])
    , ldx_darts1(new double(*) *[this->numThreads])
    , d_darts1(new double(*) *[this->numThreads])
    , ist_darts1(new int[this->numThreads])
    , iend_darts1(new int[this->numThreads])
    , jst_darts1(new int[this->numThreads])
    , jend_darts1(new int[this->numThreads])
    , nx0_darts1(new int[this->numThreads])
    , ny0_darts1(new int[this->numThreads])
    , checkInCodelets205(new _checkInCodelets205[this->numThreads])
{
    _checkInCodelets205* checkInCodelets205Ptr = (this->checkInCodelets205);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets205);
#endif
    for (size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets205Ptr) = _checkInCodelets205(2, 1, this, codeletCounter);
#else
        (*checkInCodelets205Ptr) = _checkInCodelets205(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets205Ptr).decDep();
        checkInCodelets205Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets205Ptr) = _checkInCodelets205(2, 1, this, codeletCounter);
#else
        (*checkInCodelets205Ptr) = _checkInCodelets205(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets205Ptr).decDep();
        checkInCodelets205Ptr++;
    }
    if (this->numThreads == 1) {
        this->nextCodeletsblts[0] = in_mainNextCodelet;
        this->nextSyncCodeletsblts[0] = in_mainSyncCodelet;
        this->nx_darts1[0] = in_nx;
        this->ny_darts1[0] = in_ny;
        this->nz_darts1[0] = in_nz;
        this->k_darts1[0] = in_k;
        this->omega_darts1[0] = in_omega;
        this->v_darts1[0] = in_v;
        this->ldz_darts1[0] = in_ldz;
        this->ldy_darts1[0] = in_ldy;
        this->ldx_darts1[0] = in_ldx;
        this->d_darts1[0] = in_d;
        this->ist_darts1[0] = in_ist;
        this->iend_darts1[0] = in_iend;
        this->jst_darts1[0] = in_jst;
        this->jend_darts1[0] = in_jend;
        this->nx0_darts1[0] = in_nx0;
        this->ny0_darts1[0] = in_ny0;
        this->availableCodelets[0] = 1;
    } else {
        this->nx_darts1[this->mainCodeletID] = in_nx;
        this->ny_darts1[this->mainCodeletID] = in_ny;
        this->nz_darts1[this->mainCodeletID] = in_nz;
        this->k_darts1[this->mainCodeletID] = in_k;
        this->omega_darts1[this->mainCodeletID] = in_omega;
        this->v_darts1[this->mainCodeletID] = in_v;
        this->ldz_darts1[this->mainCodeletID] = in_ldz;
        this->ldy_darts1[this->mainCodeletID] = in_ldy;
        this->ldx_darts1[this->mainCodeletID] = in_ldx;
        this->d_darts1[this->mainCodeletID] = in_d;
        this->ist_darts1[this->mainCodeletID] = in_ist;
        this->iend_darts1[this->mainCodeletID] = in_iend;
        this->jst_darts1[this->mainCodeletID] = in_jst;
        this->jend_darts1[this->mainCodeletID] = in_jend;
        this->nx0_darts1[this->mainCodeletID] = in_nx0;
        this->ny0_darts1[this->mainCodeletID] = in_ny0;
        this->nextCodeletsblts[in_mainCodeletID] = in_mainNextCodelet;
        this->nextSyncCodeletsblts[in_mainCodeletID] = in_mainSyncCodelet;
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[this->mainCodeletID].decDep();
#else
        this->availableCodelets[this->mainCodeletID] = 1;
#endif
        *(this->ptrToThisFunctionTP) = this;
    }
}
TP1::~TP1()
{
    delete[] checkInCodelets205;
    delete[] nextCodeletsblts;
    delete[] nextSyncCodeletsblts;
    delete[] nx_darts1;
    delete[] ny_darts1;
    delete[] nz_darts1;
    delete[] k_darts1;
    delete[] omega_darts1;
    delete[] v_darts1;
    delete[] ldz_darts1;
    delete[] ldy_darts1;
    delete[] ldx_darts1;
    delete[] d_darts1;
    delete[] ist_darts1;
    delete[] iend_darts1;
    delete[] jst_darts1;
    delete[] jend_darts1;
    delete[] nx0_darts1;
    delete[] ny0_darts1;
}
void TP1::setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, double(*) in_v,
    double(*) in_ldz, double(*) in_ldy, double(*) in_ldx, double(*) in_d, int in_ist, int in_iend,
    int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID)
{
    this->nx_darts1[codeletID] = in_nx;
    this->ny_darts1[codeletID] = in_ny;
    this->nz_darts1[codeletID] = in_nz;
    this->k_darts1[codeletID] = in_k;
    this->omega_darts1[codeletID] = in_omega;
    this->v_darts1[codeletID] = in_v;
    this->ldz_darts1[codeletID] = in_ldz;
    this->ldy_darts1[codeletID] = in_ldy;
    this->ldx_darts1[codeletID] = in_ldx;
    this->d_darts1[codeletID] = in_d;
    this->ist_darts1[codeletID] = in_ist;
    this->iend_darts1[codeletID] = in_iend;
    this->jst_darts1[codeletID] = in_jst;
    this->jend_darts1[codeletID] = in_jend;
    this->nx0_darts1[codeletID] = in_nx0;
    this->ny0_darts1[codeletID] = in_ny0;
}
/*TP2: TP_buts*/
void TP2::_checkInCodelets211::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Allocate each array var in the codelet */
this->inputsTPParent->d_outer2_size = 162][162][5][5;
this->inputsTPParent->d_darts2[this->getID()] = (double *)malloc(sizeof(double ) * this->inputsTPParent->d_outer2_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->tv_outer2_size = 162][5;
this->inputsTPParent->tv_darts2[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->tv_outer2_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->udx_outer2_size = 162][5][5;
this->inputsTPParent->udx_darts2[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->udx_outer2_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->udy_outer2_size = 162][5][5;
this->inputsTPParent->udy_darts2[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->udy_outer2_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->udz_outer2_size = 162][5][5;
this->inputsTPParent->udz_darts2[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->udz_outer2_size);
/*Allocate each array var in the codelet */
this->inputsTPParent->v_outer2_size = 163][163][5;
this->inputsTPParent->v_darts2[this->getID()] = (double (*)*)malloc(sizeof(double (*)) * this->inputsTPParent->v_outer2_size);

/*printing node 211: BinaryOperator*/
162 / 2 * 2 + 1;
162 / 2 * 2 + 1;

/*printing node 214: BinaryOperator*/
162 / 2 * 2 + 1;
162 / 2 * 2 + 1;
/*Signaling next codelet from last stmt in the codelet*/
/*Find and signal the next codelet*/

myTP->controlTPParent->nextCodeletsbuts[this->getID()]->decDep();
}
TP2::TP2(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
    darts::Codelet* in_mainSyncCodelet, TP2** in_ptrToThisFunctionTP, int in_nx, int in_ny,
    int in_nz, int in_k, double in_omega, double(*) in_v, double(*) in_tv, double(*) in_d,
    double(*) in_udx, double(*) in_udy, double(*) in_udz, int in_ist, int in_iend, int in_jst,
    int in_jend, int in_nx0, int in_ny0)
    : ompTP(in_numThreads, in_mainCodeletID)
    , ptrToThisFunctionTP(in_ptrToThisFunctionTP)
    , inputsTPParent(this)
    , controlTPParent(this)
    , nextCodeletsbuts(new Codelet*[in_numThreads])
    , nextSyncCodeletsbuts(new Codelet*[in_numThreads])
    , nx_darts2(new int[this->numThreads])
    , ny_darts2(new int[this->numThreads])
    , nz_darts2(new int[this->numThreads])
    , k_darts2(new int[this->numThreads])
    , omega_darts2(new double[this->numThreads])
    , v_darts2(new double(*) *[this->numThreads])
    , tv_darts2(new double(*) *[this->numThreads])
    , d_darts2(new double(*) *[this->numThreads])
    , udx_darts2(new double(*) *[this->numThreads])
    , udy_darts2(new double(*) *[this->numThreads])
    , udz_darts2(new double(*) *[this->numThreads])
    , ist_darts2(new int[this->numThreads])
    , iend_darts2(new int[this->numThreads])
    , jst_darts2(new int[this->numThreads])
    , jend_darts2(new int[this->numThreads])
    , nx0_darts2(new int[this->numThreads])
    , ny0_darts2(new int[this->numThreads])
    , checkInCodelets211(new _checkInCodelets211[this->numThreads])
{
    _checkInCodelets211* checkInCodelets211Ptr = (this->checkInCodelets211);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets211);
#endif
    for (size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets211Ptr) = _checkInCodelets211(2, 1, this, codeletCounter);
#else
        (*checkInCodelets211Ptr) = _checkInCodelets211(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets211Ptr).decDep();
        checkInCodelets211Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets211Ptr) = _checkInCodelets211(2, 1, this, codeletCounter);
#else
        (*checkInCodelets211Ptr) = _checkInCodelets211(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets211Ptr).decDep();
        checkInCodelets211Ptr++;
    }
    if (this->numThreads == 1) {
        this->nextCodeletsbuts[0] = in_mainNextCodelet;
        this->nextSyncCodeletsbuts[0] = in_mainSyncCodelet;
        this->nx_darts2[0] = in_nx;
        this->ny_darts2[0] = in_ny;
        this->nz_darts2[0] = in_nz;
        this->k_darts2[0] = in_k;
        this->omega_darts2[0] = in_omega;
        this->v_darts2[0] = in_v;
        this->tv_darts2[0] = in_tv;
        this->d_darts2[0] = in_d;
        this->udx_darts2[0] = in_udx;
        this->udy_darts2[0] = in_udy;
        this->udz_darts2[0] = in_udz;
        this->ist_darts2[0] = in_ist;
        this->iend_darts2[0] = in_iend;
        this->jst_darts2[0] = in_jst;
        this->jend_darts2[0] = in_jend;
        this->nx0_darts2[0] = in_nx0;
        this->ny0_darts2[0] = in_ny0;
        this->availableCodelets[0] = 1;
    } else {
        this->nx_darts2[this->mainCodeletID] = in_nx;
        this->ny_darts2[this->mainCodeletID] = in_ny;
        this->nz_darts2[this->mainCodeletID] = in_nz;
        this->k_darts2[this->mainCodeletID] = in_k;
        this->omega_darts2[this->mainCodeletID] = in_omega;
        this->v_darts2[this->mainCodeletID] = in_v;
        this->tv_darts2[this->mainCodeletID] = in_tv;
        this->d_darts2[this->mainCodeletID] = in_d;
        this->udx_darts2[this->mainCodeletID] = in_udx;
        this->udy_darts2[this->mainCodeletID] = in_udy;
        this->udz_darts2[this->mainCodeletID] = in_udz;
        this->ist_darts2[this->mainCodeletID] = in_ist;
        this->iend_darts2[this->mainCodeletID] = in_iend;
        this->jst_darts2[this->mainCodeletID] = in_jst;
        this->jend_darts2[this->mainCodeletID] = in_jend;
        this->nx0_darts2[this->mainCodeletID] = in_nx0;
        this->ny0_darts2[this->mainCodeletID] = in_ny0;
        this->nextCodeletsbuts[in_mainCodeletID] = in_mainNextCodelet;
        this->nextSyncCodeletsbuts[in_mainCodeletID] = in_mainSyncCodelet;
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[this->mainCodeletID].decDep();
#else
        this->availableCodelets[this->mainCodeletID] = 1;
#endif
        *(this->ptrToThisFunctionTP) = this;
    }
}
TP2::~TP2()
{
    delete[] checkInCodelets211;
    delete[] nextCodeletsbuts;
    delete[] nextSyncCodeletsbuts;
    delete[] nx_darts2;
    delete[] ny_darts2;
    delete[] nz_darts2;
    delete[] k_darts2;
    delete[] omega_darts2;
    delete[] v_darts2;
    delete[] tv_darts2;
    delete[] d_darts2;
    delete[] udx_darts2;
    delete[] udy_darts2;
    delete[] udz_darts2;
    delete[] ist_darts2;
    delete[] iend_darts2;
    delete[] jst_darts2;
    delete[] jend_darts2;
    delete[] nx0_darts2;
    delete[] ny0_darts2;
}
void TP2::setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, double(*) in_v,
    double(*) in_tv, double(*) in_d, double(*) in_udx, double(*) in_udy, double(*) in_udz,
    int in_ist, int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID)
{
    this->nx_darts2[codeletID] = in_nx;
    this->ny_darts2[codeletID] = in_ny;
    this->nz_darts2[codeletID] = in_nz;
    this->k_darts2[codeletID] = in_k;
    this->omega_darts2[codeletID] = in_omega;
    this->v_darts2[codeletID] = in_v;
    this->tv_darts2[codeletID] = in_tv;
    this->d_darts2[codeletID] = in_d;
    this->udx_darts2[codeletID] = in_udx;
    this->udy_darts2[codeletID] = in_udy;
    this->udz_darts2[codeletID] = in_udz;
    this->ist_darts2[codeletID] = in_ist;
    this->iend_darts2[codeletID] = in_iend;
    this->jst_darts2[codeletID] = in_jst;
    this->jend_darts2[codeletID] = in_jend;
    this->nx0_darts2[codeletID] = in_nx0;
    this->ny0_darts2[codeletID] = in_ny0;
}
/*TP2301: OMPParallelDirective*/
void TP2301::_barrierCodelets2301::fire(void)
{
    TP2301* myTP = static_cast<TP2301*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP2301::_checkInCodelets2303::fire(void)
{
    /*Init the vars for this region*/

    /*printing node 2303: DeclStmt*/

    /*printing node 2304: DeclStmt*/

    /*printing node 2305: DeclStmt*/

    /*printing node 2306: DeclStmt*/

    /*printing node 2307: DeclStmt*/

    /*printing node 2308: DeclStmt*/

    /*printing node 2309: DeclStmt*/

    /*printing node 2310: DeclStmt*/

    /*printing node 2311: DeclStmt*/

    /*printing node 2312: DeclStmt*/

    /*printing node 2313: DeclStmt*/

    /*printing node 2314: DeclStmt*/

    /*printing node 2315: DeclStmt*/

    /*printing node 2316: DeclStmt*/

    /*printing node 2317: DeclStmt*/

    /*printing node 2318: DeclStmt*/

    /*printing node 2319: BinaryOperator*/
    (this->inputsTPParent->dsspm_darts2301[this->getID()]) = dssp;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 2303 nextRegion: 2320 */
    myTP->controlTPParent->checkInCodelets2320[this->getID()].decDep();
}
void TP2301::_checkInCodelets2320::fire(void)
{
    /*region 2320 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP2320;
    if (idx < myTP->TPsToUse2320) {
        if (!__sync_val_compare_and_swap(&(myTP->TP2320_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(nx - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse2320;
            int minIteration = min<int>(nx, 0);
            int remainderRange = range % myTP->TPsToUse2320;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (0 < nx) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse2320 - 1) {
                lastIteration = nx;
            }
#if USEINVOKE == 1
            invoke<TP2320>(myTP, myTP->codeletsPerTP2320 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP2320Ptr[idx]));
#else
            place<TP2320>(idx, myTP, myTP->codeletsPerTP2320 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP2320Ptr[idx]));
#endif
        } else {
            if (myTP->TP2320Ptr[idx] != nullptr) {
                myTP->TP2320Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2301::_barrierCodelets2320::fire(void)
{
    TP2301* myTP = static_cast<TP2301*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets2371[codeletsCounter].decDep();
        }
    }
}
void TP2301::_checkInCodelets2371::fire(void)
{
    /*region 2371 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP2371;
    if (idx < myTP->TPsToUse2371) {
        if (!__sync_val_compare_and_swap(&(myTP->TP2371_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(nx - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse2371;
            int minIteration = min<int>(nx, 0);
            int remainderRange = range % myTP->TPsToUse2371;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (0 < nx) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse2371 - 1) {
                lastIteration = nx;
            }
#if USEINVOKE == 1
            invoke<TP2371>(myTP, myTP->codeletsPerTP2371 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP2371Ptr[idx]));
#else
            place<TP2371>(idx, myTP, myTP->codeletsPerTP2371 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP2371Ptr[idx]));
#endif
        } else {
            if (myTP->TP2371Ptr[idx] != nullptr) {
                myTP->TP2371Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2301::_barrierCodelets2371::fire(void)
{
    TP2301* myTP = static_cast<TP2301*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets2503[codeletsCounter].decDep();
        }
    }
}
void TP2301::_checkInCodelets2503::fire(void)
{

    /*printing node 2503: BinaryOperator*/
    (this->inputsTPParent->L1_darts2301[this->getID()]) = 0;

    /*printing node 2504: BinaryOperator*/
    (this->inputsTPParent->L2_darts2301[this->getID()]) = nx - 1;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 2503 nextRegion: 2506 */
    myTP->controlTPParent->checkInCodelets2506[this->getID()].decDep();
}
void TP2301::_checkInCodelets2506::fire(void)
{
    /*region 2506 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP2506;
    if (idx < myTP->TPsToUse2506) {
        if (!__sync_val_compare_and_swap(&(myTP->TP2506_alreadyLaunched[idx]), 0, 1)) {
            int range = abs((this->inputsTPParent->L2_darts2301[this->getID()])
                            - (this->inputsTPParent->L1_darts2301[this->getID()]))
                / 1;
            int rangePerCodelet = range / myTP->TPsToUse2506;
            int minIteration = min<int>((this->inputsTPParent->L2_darts2301[this->getID()]),
                (this->inputsTPParent->L1_darts2301[this->getID()]));
            int remainderRange = range % myTP->TPsToUse2506;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if ((this->inputsTPParent->L1_darts2301[this->getID()])
                < (this->inputsTPParent->L2_darts2301[this->getID()])) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse2506 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse2506 - 1) {
                lastIteration = (this->inputsTPParent->L2_darts2301[this->getID()]);
            }
#if USEINVOKE == 1
            invoke<TP2506>(myTP, myTP->codeletsPerTP2506 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP2506Ptr[idx]));
#else
            place<TP2506>(idx, myTP, myTP->codeletsPerTP2506 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP2506Ptr[idx]));
#endif
        } else {
            if (myTP->TP2506Ptr[idx] != nullptr) {
                myTP->TP2506Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2301::_barrierCodelets2506::fire(void)
{
    TP2301* myTP = static_cast<TP2301*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets2655[codeletsCounter].decDep();
        }
    }
}
void TP2301::_checkInCodelets2655::fire(void)
{
    /*region 2655 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP2655;
    if (idx < myTP->TPsToUse2655) {
        if (!__sync_val_compare_and_swap(&(myTP->TP2655_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(jend - jst) / 1;
            int rangePerCodelet = range / myTP->TPsToUse2655;
            int minIteration = min<int>(jend, jst);
            int remainderRange = range % myTP->TPsToUse2655;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (jst < jend) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse2655 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse2655 - 1) {
                lastIteration = jend;
            }
#if USEINVOKE == 1
            invoke<TP2655>(myTP, myTP->codeletsPerTP2655 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP2655Ptr[idx]));
#else
            place<TP2655>(idx, myTP, myTP->codeletsPerTP2655 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP2655Ptr[idx]));
#endif
        } else {
            if (myTP->TP2655Ptr[idx] != nullptr) {
                myTP->TP2655Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2301::_barrierCodelets2655::fire(void)
{
    TP2301* myTP = static_cast<TP2301*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets3290[codeletsCounter].decDep();
        }
    }
}
void TP2301::_checkInCodelets3290::fire(void)
{

    /*printing node 3290: BinaryOperator*/
    (this->inputsTPParent->L1_darts2301[this->getID()]) = 0;

    /*printing node 3291: BinaryOperator*/
    (this->inputsTPParent->L2_darts2301[this->getID()]) = ny - 1;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 3290 nextRegion: 3293 */
    myTP->controlTPParent->checkInCodelets3293[this->getID()].decDep();
}
void TP2301::_checkInCodelets3293::fire(void)
{
    /*region 3293 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP3293;
    if (idx < myTP->TPsToUse3293) {
        if (!__sync_val_compare_and_swap(&(myTP->TP3293_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse3293;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse3293;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (ist < iend) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse3293 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse3293 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP3293>(myTP, myTP->codeletsPerTP3293 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP3293Ptr[idx]));
#else
            place<TP3293>(idx, myTP, myTP->codeletsPerTP3293 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP3293Ptr[idx]));
#endif
        } else {
            if (myTP->TP3293Ptr[idx] != nullptr) {
                myTP->TP3293Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2301::_barrierCodelets3293::fire(void)
{
    TP2301* myTP = static_cast<TP2301*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets3442[codeletsCounter].decDep();
        }
    }
}
void TP2301::_checkInCodelets3442::fire(void)
{
    /*region 3442 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP3442;
    if (idx < myTP->TPsToUse3442) {
        if (!__sync_val_compare_and_swap(&(myTP->TP3442_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse3442;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse3442;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (ist < iend) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse3442 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse3442 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP3442>(myTP, myTP->codeletsPerTP3442 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP3442Ptr[idx]));
#else
            place<TP3442>(idx, myTP, myTP->codeletsPerTP3442 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP3442Ptr[idx]));
#endif
        } else {
            if (myTP->TP3442Ptr[idx] != nullptr) {
                myTP->TP3442Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2301::_barrierCodelets3442::fire(void)
{
    TP2301* myTP = static_cast<TP2301*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets4077[codeletsCounter].decDep();
        }
    }
}
void TP2301::_checkInCodelets4077::fire(void)
{
    /*region 4077 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP4077;
    if (idx < myTP->TPsToUse4077) {
        if (!__sync_val_compare_and_swap(&(myTP->TP4077_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse4077;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse4077;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (ist < iend) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse4077 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse4077 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP4077>(myTP, myTP->codeletsPerTP4077 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP4077Ptr[idx]));
#else
            place<TP4077>(idx, myTP, myTP->codeletsPerTP4077 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP4077Ptr[idx]));
#endif
        } else {
            if (myTP->TP4077Ptr[idx] != nullptr) {
                myTP->TP4077Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2301::_barrierCodelets4077::fire(void)
{
    TP2301* myTP = static_cast<TP2301*>(myTP_);
    myTP->TPParent->barrierCodelets2301[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets2301[0]));
}
TP2301::TP2301(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts2301(new int[this->numThreads]) /*VARIABLE*/
    , L2_darts2301(new int[this->numThreads]) /*VARIABLE*/
    , dsspm_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , eta_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , i_darts2301(new int[this->numThreads]) /*VARIABLE*/
    , iend1_darts2301(new int[this->numThreads]) /*VARIABLE*/
    , iglob_darts2301(new int[this->numThreads]) /*VARIABLE*/
    , ist1_darts2301(new int[this->numThreads]) /*VARIABLE*/
    , j_darts2301(new int[this->numThreads]) /*VARIABLE*/
    , jend1_darts2301(new int[this->numThreads]) /*VARIABLE*/
    , jglob_darts2301(new int[this->numThreads]) /*VARIABLE*/
    , jst1_darts2301(new int[this->numThreads]) /*VARIABLE*/
    , k_darts2301(new int[this->numThreads]) /*VARIABLE*/
    , m_darts2301(new int[this->numThreads]) /*VARIABLE*/
    , q_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , tmp_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u21_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u21i_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u21im1_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u21j_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u21jm1_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u21k_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u21km1_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u31_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u31i_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u31im1_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u31j_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u31jm1_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u31k_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u31km1_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u41_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u41i_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u41im1_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u41j_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u41jm1_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u41k_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u41km1_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u51i_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u51im1_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u51j_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u51jm1_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u51k_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , u51km1_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , xi_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , zeta_darts2301(new double[this->numThreads]) /*VARIABLE*/
    , TP2320Ptr(new TP2320*[NUMTPS2320])
    , TP2320_alreadyLaunched(new size_t[NUMTPS2320])
    , numTPsSet2320(0)
    , numTPsReady2320(0)
    , TPsToUse2320(NUMTPS2320)
    , codeletsPerTP2320(this->numThreads / NUMTPS2320)
    , totalCodelets2320(this->TPsToUse2320 * this->codeletsPerTP2320)
    , TP2371Ptr(new TP2371*[NUMTPS2371])
    , TP2371_alreadyLaunched(new size_t[NUMTPS2371])
    , numTPsSet2371(0)
    , numTPsReady2371(0)
    , TPsToUse2371(NUMTPS2371)
    , codeletsPerTP2371(this->numThreads / NUMTPS2371)
    , totalCodelets2371(this->TPsToUse2371 * this->codeletsPerTP2371)
    , TP2506Ptr(new TP2506*[NUMTPS2506])
    , TP2506_alreadyLaunched(new size_t[NUMTPS2506])
    , numTPsSet2506(0)
    , numTPsReady2506(0)
    , TPsToUse2506(NUMTPS2506)
    , codeletsPerTP2506(this->numThreads / NUMTPS2506)
    , totalCodelets2506(this->TPsToUse2506 * this->codeletsPerTP2506)
    , TP2655Ptr(new TP2655*[NUMTPS2655])
    , TP2655_alreadyLaunched(new size_t[NUMTPS2655])
    , numTPsSet2655(0)
    , numTPsReady2655(0)
    , TPsToUse2655(NUMTPS2655)
    , codeletsPerTP2655(this->numThreads / NUMTPS2655)
    , totalCodelets2655(this->TPsToUse2655 * this->codeletsPerTP2655)
    , TP3293Ptr(new TP3293*[NUMTPS3293])
    , TP3293_alreadyLaunched(new size_t[NUMTPS3293])
    , numTPsSet3293(0)
    , numTPsReady3293(0)
    , TPsToUse3293(NUMTPS3293)
    , codeletsPerTP3293(this->numThreads / NUMTPS3293)
    , totalCodelets3293(this->TPsToUse3293 * this->codeletsPerTP3293)
    , TP3442Ptr(new TP3442*[NUMTPS3442])
    , TP3442_alreadyLaunched(new size_t[NUMTPS3442])
    , numTPsSet3442(0)
    , numTPsReady3442(0)
    , TPsToUse3442(NUMTPS3442)
    , codeletsPerTP3442(this->numThreads / NUMTPS3442)
    , totalCodelets3442(this->TPsToUse3442 * this->codeletsPerTP3442)
    , TP4077Ptr(new TP4077*[NUMTPS4077])
    , TP4077_alreadyLaunched(new size_t[NUMTPS4077])
    , numTPsSet4077(0)
    , numTPsReady4077(0)
    , TPsToUse4077(NUMTPS4077)
    , codeletsPerTP4077(this->numThreads / NUMTPS4077)
    , totalCodelets4077(this->TPsToUse4077 * this->codeletsPerTP4077)
    , barrierCodelets2301(new _barrierCodelets2301[1])
    , checkInCodelets2303(new _checkInCodelets2303[this->numThreads])
    , checkInCodelets2320(new _checkInCodelets2320[this->numThreads])
    , barrierCodelets2320(new _barrierCodelets2320[1])
    , checkInCodelets2371(new _checkInCodelets2371[this->numThreads])
    , barrierCodelets2371(new _barrierCodelets2371[1])
    , checkInCodelets2503(new _checkInCodelets2503[this->numThreads])
    , checkInCodelets2506(new _checkInCodelets2506[this->numThreads])
    , barrierCodelets2506(new _barrierCodelets2506[1])
    , checkInCodelets2655(new _checkInCodelets2655[this->numThreads])
    , barrierCodelets2655(new _barrierCodelets2655[1])
    , checkInCodelets3290(new _checkInCodelets3290[this->numThreads])
    , checkInCodelets3293(new _checkInCodelets3293[this->numThreads])
    , barrierCodelets3293(new _barrierCodelets3293[1])
    , checkInCodelets3442(new _checkInCodelets3442[this->numThreads])
    , barrierCodelets3442(new _barrierCodelets3442[1])
    , checkInCodelets4077(new _checkInCodelets4077[this->numThreads])
    , barrierCodelets4077(new _barrierCodelets4077[1])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets2301[0] = _barrierCodelets2301(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets4077[0] = _barrierCodelets4077(NUMTPS4077, NUMTPS4077, this, 0);
    barrierCodelets3442[0] = _barrierCodelets3442(NUMTPS3442, NUMTPS3442, this, 0);
    barrierCodelets3293[0] = _barrierCodelets3293(NUMTPS3293, NUMTPS3293, this, 0);
    barrierCodelets2655[0] = _barrierCodelets2655(NUMTPS2655, NUMTPS2655, this, 0);
    barrierCodelets2506[0] = _barrierCodelets2506(NUMTPS2506, NUMTPS2506, this, 0);
    barrierCodelets2371[0] = _barrierCodelets2371(NUMTPS2371, NUMTPS2371, this, 0);
    barrierCodelets2320[0] = _barrierCodelets2320(NUMTPS2320, NUMTPS2320, this, 0);
    _checkInCodelets4077* checkInCodelets4077Ptr = (this->checkInCodelets4077);
    for (int i = 0; i < NUMTPS4077; i++) {
        TP4077Ptr[i] = nullptr;
        TP4077_alreadyLaunched[i] = 0;
    }
    _checkInCodelets3442* checkInCodelets3442Ptr = (this->checkInCodelets3442);
    for (int i = 0; i < NUMTPS3442; i++) {
        TP3442Ptr[i] = nullptr;
        TP3442_alreadyLaunched[i] = 0;
    }
    _checkInCodelets3293* checkInCodelets3293Ptr = (this->checkInCodelets3293);
    for (int i = 0; i < NUMTPS3293; i++) {
        TP3293Ptr[i] = nullptr;
        TP3293_alreadyLaunched[i] = 0;
    }
    _checkInCodelets3290* checkInCodelets3290Ptr = (this->checkInCodelets3290);
    _checkInCodelets2655* checkInCodelets2655Ptr = (this->checkInCodelets2655);
    for (int i = 0; i < NUMTPS2655; i++) {
        TP2655Ptr[i] = nullptr;
        TP2655_alreadyLaunched[i] = 0;
    }
    _checkInCodelets2506* checkInCodelets2506Ptr = (this->checkInCodelets2506);
    for (int i = 0; i < NUMTPS2506; i++) {
        TP2506Ptr[i] = nullptr;
        TP2506_alreadyLaunched[i] = 0;
    }
    _checkInCodelets2503* checkInCodelets2503Ptr = (this->checkInCodelets2503);
    _checkInCodelets2371* checkInCodelets2371Ptr = (this->checkInCodelets2371);
    for (int i = 0; i < NUMTPS2371; i++) {
        TP2371Ptr[i] = nullptr;
        TP2371_alreadyLaunched[i] = 0;
    }
    _checkInCodelets2320* checkInCodelets2320Ptr = (this->checkInCodelets2320);
    for (int i = 0; i < NUMTPS2320; i++) {
        TP2320Ptr[i] = nullptr;
        TP2320_alreadyLaunched[i] = 0;
    }
    _checkInCodelets2303* checkInCodelets2303Ptr = (this->checkInCodelets2303);
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets4077Ptr) = _checkInCodelets4077(1, 1, this, codeletCounter);
        checkInCodelets4077Ptr++;
        (*checkInCodelets3442Ptr) = _checkInCodelets3442(1, 1, this, codeletCounter);
        checkInCodelets3442Ptr++;
        (*checkInCodelets3293Ptr) = _checkInCodelets3293(1, 1, this, codeletCounter);
        checkInCodelets3293Ptr++;
        (*checkInCodelets3290Ptr) = _checkInCodelets3290(1, 1, this, codeletCounter);
        checkInCodelets3290Ptr++;
        (*checkInCodelets2655Ptr) = _checkInCodelets2655(1, 1, this, codeletCounter);
        checkInCodelets2655Ptr++;
        (*checkInCodelets2506Ptr) = _checkInCodelets2506(1, 1, this, codeletCounter);
        checkInCodelets2506Ptr++;
        (*checkInCodelets2503Ptr) = _checkInCodelets2503(1, 1, this, codeletCounter);
        checkInCodelets2503Ptr++;
        (*checkInCodelets2371Ptr) = _checkInCodelets2371(1, 1, this, codeletCounter);
        checkInCodelets2371Ptr++;
        (*checkInCodelets2320Ptr) = _checkInCodelets2320(1, 1, this, codeletCounter);
        checkInCodelets2320Ptr++;
        (*checkInCodelets2303Ptr) = _checkInCodelets2303(1, 1, this, codeletCounter);
        (*checkInCodelets2303Ptr).decDep();
        checkInCodelets2303Ptr++;
    }
}
TP2301::~TP2301()
{
    delete[] L1_darts2301;
    delete[] L2_darts2301;
    delete[] dsspm_darts2301;
    delete[] eta_darts2301;
    delete[] i_darts2301;
    delete[] iend1_darts2301;
    delete[] iglob_darts2301;
    delete[] ist1_darts2301;
    delete[] j_darts2301;
    delete[] jend1_darts2301;
    delete[] jglob_darts2301;
    delete[] jst1_darts2301;
    delete[] k_darts2301;
    delete[] m_darts2301;
    delete[] q_darts2301;
    delete[] tmp_darts2301;
    delete[] u21_darts2301;
    delete[] u21i_darts2301;
    delete[] u21im1_darts2301;
    delete[] u21j_darts2301;
    delete[] u21jm1_darts2301;
    delete[] u21k_darts2301;
    delete[] u21km1_darts2301;
    delete[] u31_darts2301;
    delete[] u31i_darts2301;
    delete[] u31im1_darts2301;
    delete[] u31j_darts2301;
    delete[] u31jm1_darts2301;
    delete[] u31k_darts2301;
    delete[] u31km1_darts2301;
    delete[] u41_darts2301;
    delete[] u41i_darts2301;
    delete[] u41im1_darts2301;
    delete[] u41j_darts2301;
    delete[] u41jm1_darts2301;
    delete[] u41k_darts2301;
    delete[] u41km1_darts2301;
    delete[] u51i_darts2301;
    delete[] u51im1_darts2301;
    delete[] u51j_darts2301;
    delete[] u51jm1_darts2301;
    delete[] u51k_darts2301;
    delete[] u51km1_darts2301;
    delete[] xi_darts2301;
    delete[] zeta_darts2301;
    delete[] barrierCodelets2301;
    delete[] barrierCodelets4077;
    delete[] checkInCodelets4077;
    delete[] barrierCodelets3442;
    delete[] checkInCodelets3442;
    delete[] barrierCodelets3293;
    delete[] checkInCodelets3293;
    delete[] checkInCodelets3290;
    delete[] barrierCodelets2655;
    delete[] checkInCodelets2655;
    delete[] barrierCodelets2506;
    delete[] checkInCodelets2506;
    delete[] checkInCodelets2503;
    delete[] barrierCodelets2371;
    delete[] checkInCodelets2371;
    delete[] barrierCodelets2320;
    delete[] checkInCodelets2320;
    delete[] checkInCodelets2303;
}
/*TP2320: OMPForDirective*/
void TP2320::_barrierCodelets2320::fire(void)
{
    TP2320* myTP = static_cast<TP2320*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets2320[0].decDep();
}
bool TP2320::requestNewRangeIterations2320(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet2320 * codeletID;
        int tempEndRange = rangePerCodelet2320 * (codeletID + 1);
        if (remainderRange2320 != 0) {
            if (codeletID < (uint32_t)remainderRange2320) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange2320;
                tempEndRange += remainderRange2320;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration2320;
        tempEndRange = tempEndRange * 1 + minIteration2320;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration2320 < lastIteration2320) {
            (this->inputsTPParent->i_darts2320[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts2320[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration2320;
        }
    }
    return isThereNewIteration;
}
void TP2320::_checkInCodelets2321::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 2321: ForStmt*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    int* i = &(this->inputsTPParent->i_darts2320[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts2320[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts2320[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts2320[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations2320(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets2320[0].decDep();
        return;
    }
    for (int i_darts_counter_temp2320 = (*i); i_darts_counter_temp2320 < endRange
         && i_darts_counter_temp2320 < this->inputsTPParent->lastIteration2320;
         i_darts_counter_temp2320++) {
        {
            {
                /*Loop's init*/
                (*j) = 0;
                int j_darts_counter_temp2320 = (*j);
                for (; j_darts_counter_temp2320 < ny; j_darts_counter_temp2320++) {
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp2320 = (*k);
                        for (; k_darts_counter_temp2320 < nz; k_darts_counter_temp2320++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp2320 = (*m);
                                for (; m_darts_counter_temp2320 < 5; m_darts_counter_temp2320++) {
                                    frct[(i_darts_counter_temp2320)][j_darts_counter_temp2320]
                                        [k_darts_counter_temp2320][m_darts_counter_temp2320]
                                        = 0.;
                                }
                                (*m) = m_darts_counter_temp2320;
                            }
                        }
                        (*k) = k_darts_counter_temp2320;
                    }
                }
                (*j) = j_darts_counter_temp2320;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets2320[0].decDep();
}
TP2320::TP2320(int in_numThreads, int in_mainCodeletID, TP2301* in_TPParent, int in_initIteration,
    int in_lastIteration, TP2320** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts2320(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts2320(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts2320(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts2320(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration2320(in_initIteration)
    , lastIteration2320(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets2320(new _barrierCodelets2320[1])
    , checkInCodelets2321(new _checkInCodelets2321[this->numThreads])
{
    /*Initialize the loop parameters*/
    range2320 = abs(lastIteration2320 - initIteration2320) / 1;
    rangePerCodelet2320 = range2320 / numThreads;
    minIteration2320 = min<int>(lastIteration2320, initIteration2320);
    remainderRange2320 = range2320 % numThreads;
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets2320[0] = _barrierCodelets2320(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets2321* checkInCodelets2321Ptr = (this->checkInCodelets2321);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets2321);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets2321Ptr) = _checkInCodelets2321(2, 1, this, codeletCounter);
#else
        (*checkInCodelets2321Ptr) = _checkInCodelets2321(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets2321Ptr).decDep();
        checkInCodelets2321Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP2320::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets2321[localID].setID(codeletID);
    this->checkInCodelets2321[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets2321[localID + this->baseNumThreads * i]
            = _checkInCodelets2321(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets2321[localID + this->baseNumThreads * i]
            = _checkInCodelets2321(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets2321[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets2321[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP2320::~TP2320()
{
    delete[] barrierCodelets2320;
    delete[] checkInCodelets2321;
}
/*TP2371: OMPForDirective*/
void TP2371::_barrierCodelets2371::fire(void)
{
    TP2371* myTP = static_cast<TP2371*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets2371[0].decDep();
}
bool TP2371::requestNewRangeIterations2371(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet2371 * codeletID;
        int tempEndRange = rangePerCodelet2371 * (codeletID + 1);
        if (remainderRange2371 != 0) {
            if (codeletID < (uint32_t)remainderRange2371) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange2371;
                tempEndRange += remainderRange2371;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration2371;
        tempEndRange = tempEndRange * 1 + minIteration2371;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration2371 < lastIteration2371) {
            (this->inputsTPParent->i_darts2371[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts2371[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration2371;
        }
    }
    return isThereNewIteration;
}
void TP2371::_checkInCodelets2372::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->eta_darts2371[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->eta_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->iglob_darts2371[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jglob_darts2371[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->xi_darts2371[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->xi_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->zeta_darts2371[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->zeta_darts2301[this->getID()]);

    /*printing node 2372: ForStmt*/
    /*var: eta*/
    /*var: i*/
    /*var: iglob*/
    /*var: j*/
    /*var: jglob*/
    /*var: k*/
    /*var: m*/
    /*var: xi*/
    /*var: zeta*/
    double** eta = &(this->inputsTPParent->eta_darts2371[this->getLocalID()]);
    (void)eta /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts2371[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iglob = &(this->inputsTPParent->iglob_darts2371[this->getLocalID()]);
    (void)iglob /*OMP_SHARED_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts2371[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jglob = &(this->inputsTPParent->jglob_darts2371[this->getLocalID()]);
    (void)jglob /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts2371[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts2371[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** xi = &(this->inputsTPParent->xi_darts2371[this->getLocalID()]);
    (void)xi /*OMP_SHARED_PRIVATE*/;
    double** zeta = &(this->inputsTPParent->zeta_darts2371[this->getLocalID()]);
    (void)zeta /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations2371(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets2371[0].decDep();
        return;
    }
    for (int i_darts_counter_temp2371 = (*i); i_darts_counter_temp2371 < endRange
         && i_darts_counter_temp2371 < this->inputsTPParent->lastIteration2371;
         i_darts_counter_temp2371++) {
        {
            (*(*iglob)) = (i_darts_counter_temp2371);
            (*(*xi)) = ((double)((*(*iglob)))) / (nx0 - 1);
            {
                /*Loop's init*/
                (*j) = 0;
                int j_darts_counter_temp2371 = (*j);
                for (; j_darts_counter_temp2371 < ny; j_darts_counter_temp2371++) {
                    (*(*jglob)) = j_darts_counter_temp2371;
                    (*(*eta)) = ((double)((*(*jglob)))) / (ny0 - 1);
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp2371 = (*k);
                        for (; k_darts_counter_temp2371 < nz; k_darts_counter_temp2371++) {
                            (*(*zeta)) = ((double)(k_darts_counter_temp2371)) / (nz - 1);
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp2371 = (*m);
                                for (; m_darts_counter_temp2371 < 5; m_darts_counter_temp2371++) {
                                    rsd[(i_darts_counter_temp2371)][j_darts_counter_temp2371]
                                       [k_darts_counter_temp2371][m_darts_counter_temp2371]
                                        = ce[m_darts_counter_temp2371][0]
                                        + ce[m_darts_counter_temp2371][1] * (*(*xi))
                                        + ce[m_darts_counter_temp2371][2] * (*(*eta))
                                        + ce[m_darts_counter_temp2371][3] * (*(*zeta))
                                        + ce[m_darts_counter_temp2371][4] * (*(*xi)) * (*(*xi))
                                        + ce[m_darts_counter_temp2371][5] * (*(*eta)) * (*(*eta))
                                        + ce[m_darts_counter_temp2371][6] * (*(*zeta)) * (*(*zeta))
                                        + ce[m_darts_counter_temp2371][7] * (*(*xi)) * (*(*xi))
                                            * (*(*xi))
                                        + ce[m_darts_counter_temp2371][8] * (*(*eta)) * (*(*eta))
                                            * (*(*eta))
                                        + ce[m_darts_counter_temp2371][9] * (*(*zeta)) * (*(*zeta))
                                            * (*(*zeta))
                                        + ce[m_darts_counter_temp2371][10] * (*(*xi)) * (*(*xi))
                                            * (*(*xi)) * (*(*xi))
                                        + ce[m_darts_counter_temp2371][11] * (*(*eta)) * (*(*eta))
                                            * (*(*eta)) * (*(*eta))
                                        + ce[m_darts_counter_temp2371][12] * (*(*zeta)) * (*(*zeta))
                                            * (*(*zeta)) * (*(*zeta));
                                }
                                (*m) = m_darts_counter_temp2371;
                            }
                        }
                        (*k) = k_darts_counter_temp2371;
                    }
                }
                (*j) = j_darts_counter_temp2371;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets2371[0].decDep();
}
TP2371::TP2371(int in_numThreads, int in_mainCodeletID, TP2301* in_TPParent, int in_initIteration,
    int in_lastIteration, TP2371** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , eta_darts2371(new double*[this->numThreads])
    , i_darts2371(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iglob_darts2371(new int*[this->numThreads])
    , j_darts2371(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jglob_darts2371(new int*[this->numThreads])
    , k_darts2371(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts2371(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , xi_darts2371(new double*[this->numThreads])
    , zeta_darts2371(new double*[this->numThreads])
    , initIteration2371(in_initIteration)
    , lastIteration2371(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets2371(new _barrierCodelets2371[1])
    , checkInCodelets2372(new _checkInCodelets2372[this->numThreads])
{
    /*Initialize the loop parameters*/
    range2371 = abs(lastIteration2371 - initIteration2371) / 1;
    rangePerCodelet2371 = range2371 / numThreads;
    minIteration2371 = min<int>(lastIteration2371, initIteration2371);
    remainderRange2371 = range2371 % numThreads;
    /*Initialize inputs and vars.*/
    this->eta_darts2371
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->iglob_darts2371 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jglob_darts2371 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->xi_darts2371
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->zeta_darts2371
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets2371[0] = _barrierCodelets2371(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets2372* checkInCodelets2372Ptr = (this->checkInCodelets2372);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets2372);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets2372Ptr) = _checkInCodelets2372(2, 1, this, codeletCounter);
#else
        (*checkInCodelets2372Ptr) = _checkInCodelets2372(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets2372Ptr).decDep();
        checkInCodelets2372Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP2371::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets2372[localID].setID(codeletID);
    this->checkInCodelets2372[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets2372[localID + this->baseNumThreads * i]
            = _checkInCodelets2372(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets2372[localID + this->baseNumThreads * i]
            = _checkInCodelets2372(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets2372[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets2372[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP2371::~TP2371()
{
    delete[] eta_darts2371;
    delete[] iglob_darts2371;
    delete[] jglob_darts2371;
    delete[] xi_darts2371;
    delete[] zeta_darts2371;
    delete[] barrierCodelets2371;
    delete[] checkInCodelets2372;
}
/*TP2506: OMPForDirective*/
void TP2506::_barrierCodelets2506::fire(void)
{
    TP2506* myTP = static_cast<TP2506*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets2506[0].decDep();
}
bool TP2506::requestNewRangeIterations2506(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet2506 * codeletID;
        int tempEndRange = rangePerCodelet2506 * (codeletID + 1);
        if (remainderRange2506 != 0) {
            if (codeletID < (uint32_t)remainderRange2506) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange2506;
                tempEndRange += remainderRange2506;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration2506;
        tempEndRange = tempEndRange * 1 + minIteration2506;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration2506 < lastIteration2506) {
            (this->inputsTPParent->i_darts2506[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts2506[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration2506;
        }
    }
    return isThereNewIteration;
}
void TP2506::_checkInCodelets2507::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L1_darts2506[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L1_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L2_darts2506[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->q_darts2506[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21_darts2506[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21_darts2301[this->getID()]);

    /*printing node 2507: ForStmt*/
    /*var: L1*/
    /*var: L2*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: q*/
    /*var: u21*/
    int* i = &(this->inputsTPParent->i_darts2506[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts2506[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts2506[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts2506[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** u21 = &(this->inputsTPParent->u21_darts2506[this->getLocalID()]);
    (void)u21 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations2506(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets2506[0].decDep();
        return;
    }
    for (int i_darts_counter_temp2506 = (*i); i_darts_counter_temp2506 <= endRange
         && i_darts_counter_temp2506 <= this->inputsTPParent->lastIteration2506;
         i_darts_counter_temp2506++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp2506 = (*j);
                for (; j_darts_counter_temp2506 <= jend; j_darts_counter_temp2506++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp2506 = (*k);
                        for (; k_darts_counter_temp2506 < nz - 1; k_darts_counter_temp2506++) {
                            flux[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                [k_darts_counter_temp2506][0]
                                = rsd[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                     [k_darts_counter_temp2506][1];
                            (*(*u21)) = rsd[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                           [k_darts_counter_temp2506][1]
                                / rsd[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                     [k_darts_counter_temp2506][0];
                            (*(*q)) = 0.5
                                * (rsd[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                      [k_darts_counter_temp2506][1]
                                        * rsd[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                             [k_darts_counter_temp2506][1]
                                    + rsd[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                         [k_darts_counter_temp2506][2]
                                        * rsd[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                             [k_darts_counter_temp2506][2]
                                    + rsd[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                         [k_darts_counter_temp2506][3]
                                        * rsd[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                             [k_darts_counter_temp2506][3])
                                / rsd[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                     [k_darts_counter_temp2506][0];
                            flux[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                [k_darts_counter_temp2506][1]
                                = rsd[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                     [k_darts_counter_temp2506][1]
                                    * (*(*u21))
                                + 0.40000000000000002
                                    * (rsd[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                          [k_darts_counter_temp2506][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                [k_darts_counter_temp2506][2]
                                = rsd[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                     [k_darts_counter_temp2506][2]
                                * (*(*u21));
                            flux[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                [k_darts_counter_temp2506][3]
                                = rsd[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                     [k_darts_counter_temp2506][3]
                                * (*(*u21));
                            flux[(i_darts_counter_temp2506)][j_darts_counter_temp2506]
                                [k_darts_counter_temp2506][4]
                                = (1.3999999999999999
                                          * rsd[(i_darts_counter_temp2506)]
                                               [j_darts_counter_temp2506][k_darts_counter_temp2506]
                                               [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u21));
                        }
                        (*k) = k_darts_counter_temp2506;
                    }
                }
                (*j) = j_darts_counter_temp2506;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets2506[0].decDep();
}
TP2506::TP2506(int in_numThreads, int in_mainCodeletID, TP2301* in_TPParent, int in_initIteration,
    int in_lastIteration, TP2506** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts2506(new int*[this->numThreads])
    , L2_darts2506(new int*[this->numThreads])
    , i_darts2506(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts2506(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts2506(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts2506(new double*[this->numThreads])
    , u21_darts2506(new double*[this->numThreads])
    , initIteration2506(in_initIteration)
    , lastIteration2506(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets2506(new _barrierCodelets2506[1])
    , checkInCodelets2507(new _checkInCodelets2507[this->numThreads])
{
    /*Initialize the loop parameters*/
    range2506 = abs(lastIteration2506 - initIteration2506) / 1;
    rangePerCodelet2506 = range2506 / numThreads;
    minIteration2506 = min<int>(lastIteration2506, initIteration2506);
    remainderRange2506 = range2506 % numThreads;
    /*Initialize inputs and vars.*/
    this->L1_darts2506 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->L2_darts2506 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->q_darts2506 = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21_darts2506
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets2506[0] = _barrierCodelets2506(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets2507* checkInCodelets2507Ptr = (this->checkInCodelets2507);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets2507);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets2507Ptr) = _checkInCodelets2507(2, 1, this, codeletCounter);
#else
        (*checkInCodelets2507Ptr) = _checkInCodelets2507(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets2507Ptr).decDep();
        checkInCodelets2507Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP2506::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets2507[localID].setID(codeletID);
    this->checkInCodelets2507[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets2507[localID + this->baseNumThreads * i]
            = _checkInCodelets2507(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets2507[localID + this->baseNumThreads * i]
            = _checkInCodelets2507(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets2507[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets2507[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP2506::~TP2506()
{
    delete[] L1_darts2506;
    delete[] L2_darts2506;
    delete[] q_darts2506;
    delete[] u21_darts2506;
    delete[] barrierCodelets2506;
    delete[] checkInCodelets2507;
}
/*TP2655: OMPForDirective*/
void TP2655::_barrierCodelets2655::fire(void)
{
    TP2655* myTP = static_cast<TP2655*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets2655[0].decDep();
}
bool TP2655::requestNewRangeIterations2655(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet2655 * codeletID;
        int tempEndRange = rangePerCodelet2655 * (codeletID + 1);
        if (remainderRange2655 != 0) {
            if (codeletID < (uint32_t)remainderRange2655) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange2655;
                tempEndRange += remainderRange2655;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration2655;
        tempEndRange = tempEndRange * 1 + minIteration2655;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration2655 < lastIteration2655) {
            (this->inputsTPParent->j_darts2655[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->j_darts2655[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration2655;
        }
    }
    return isThereNewIteration;
}
void TP2655::_checkInCodelets2656::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L2_darts2655[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->dsspm_darts2655[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->dsspm_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->iend1_darts2655[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iend1_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ist1_darts2655[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->ist1_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts2655[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21i_darts2655[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21i_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21im1_darts2655[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21im1_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31i_darts2655[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31i_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31im1_darts2655[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31im1_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41i_darts2655[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41i_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41im1_darts2655[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41im1_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51i_darts2655[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51i_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51im1_darts2655[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51im1_darts2301[this->getID()]);

    /*printing node 2656: ForStmt*/
    /*var: L2*/
    /*var: dsspm*/
    /*var: i*/
    /*var: iend1*/
    /*var: ist1*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    /*var: tmp*/
    /*var: u21i*/
    /*var: u21im1*/
    /*var: u31i*/
    /*var: u31im1*/
    /*var: u41i*/
    /*var: u41im1*/
    /*var: u51i*/
    /*var: u51im1*/
    int** L2 = &(this->inputsTPParent->L2_darts2655[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    double** dsspm = &(this->inputsTPParent->dsspm_darts2655[this->getLocalID()]);
    (void)dsspm /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts2655[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iend1 = &(this->inputsTPParent->iend1_darts2655[this->getLocalID()]);
    (void)iend1 /*OMP_SHARED_PRIVATE*/;
    int** ist1 = &(this->inputsTPParent->ist1_darts2655[this->getLocalID()]);
    (void)ist1 /*OMP_SHARED_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts2655[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts2655[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts2655[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts2655[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21i = &(this->inputsTPParent->u21i_darts2655[this->getLocalID()]);
    (void)u21i /*OMP_SHARED_PRIVATE*/;
    double** u21im1 = &(this->inputsTPParent->u21im1_darts2655[this->getLocalID()]);
    (void)u21im1 /*OMP_SHARED_PRIVATE*/;
    double** u31i = &(this->inputsTPParent->u31i_darts2655[this->getLocalID()]);
    (void)u31i /*OMP_SHARED_PRIVATE*/;
    double** u31im1 = &(this->inputsTPParent->u31im1_darts2655[this->getLocalID()]);
    (void)u31im1 /*OMP_SHARED_PRIVATE*/;
    double** u41i = &(this->inputsTPParent->u41i_darts2655[this->getLocalID()]);
    (void)u41i /*OMP_SHARED_PRIVATE*/;
    double** u41im1 = &(this->inputsTPParent->u41im1_darts2655[this->getLocalID()]);
    (void)u41im1 /*OMP_SHARED_PRIVATE*/;
    double** u51i = &(this->inputsTPParent->u51i_darts2655[this->getLocalID()]);
    (void)u51i /*OMP_SHARED_PRIVATE*/;
    double** u51im1 = &(this->inputsTPParent->u51im1_darts2655[this->getLocalID()]);
    (void)u51im1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations2655(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets2655[0].decDep();
        return;
    }
    for (int j_darts_counter_temp2655 = (*j); j_darts_counter_temp2655 <= endRange
         && j_darts_counter_temp2655 <= this->inputsTPParent->lastIteration2655;
         j_darts_counter_temp2655++) {
        {
            {
                /*Loop's init*/
                (*k) = 1;
                int k_darts_counter_temp2655 = (*k);
                for (; k_darts_counter_temp2655 <= nz - 2; k_darts_counter_temp2655++) {
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp2655 = (*i);
                        for (; i_darts_counter_temp2655 <= iend; i_darts_counter_temp2655++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp2655 = (*m);
                                for (; m_darts_counter_temp2655 < 5; m_darts_counter_temp2655++) {
                                    frct[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                        [k_darts_counter_temp2655][m_darts_counter_temp2655]
                                        = frct[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                              [k_darts_counter_temp2655][m_darts_counter_temp2655]
                                        - tx2
                                            * (flux[i_darts_counter_temp2655 + 1]
                                                   [(j_darts_counter_temp2655)]
                                                   [k_darts_counter_temp2655]
                                                   [m_darts_counter_temp2655]
                                                - flux[i_darts_counter_temp2655 - 1]
                                                      [(j_darts_counter_temp2655)]
                                                      [k_darts_counter_temp2655]
                                                      [m_darts_counter_temp2655]);
                                }
                                (*m) = m_darts_counter_temp2655;
                            }
                        }
                        (*i) = i_darts_counter_temp2655;
                    }
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp2655 = (*i);
                        for (; i_darts_counter_temp2655 <= (*(*L2)); i_darts_counter_temp2655++) {
                            (*(*tmp)) = 1.
                                / rsd[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                     [k_darts_counter_temp2655][0];
                            (*(*u21i)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                     [k_darts_counter_temp2655][1];
                            (*(*u31i)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                     [k_darts_counter_temp2655][2];
                            (*(*u41i)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                     [k_darts_counter_temp2655][3];
                            (*(*u51i)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                     [k_darts_counter_temp2655][4];
                            (*(*tmp)) = 1.
                                / rsd[i_darts_counter_temp2655 - 1][(j_darts_counter_temp2655)]
                                     [k_darts_counter_temp2655][0];
                            (*(*u21im1)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2655 - 1][(j_darts_counter_temp2655)]
                                     [k_darts_counter_temp2655][1];
                            (*(*u31im1)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2655 - 1][(j_darts_counter_temp2655)]
                                     [k_darts_counter_temp2655][2];
                            (*(*u41im1)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2655 - 1][(j_darts_counter_temp2655)]
                                     [k_darts_counter_temp2655][3];
                            (*(*u51im1)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2655 - 1][(j_darts_counter_temp2655)]
                                     [k_darts_counter_temp2655][4];
                            flux[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                [k_darts_counter_temp2655][1]
                                = (4. / 3.) * tx3 * ((*(*u21i)) - (*(*u21im1)));
                            flux[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                [k_darts_counter_temp2655][2]
                                = tx3 * ((*(*u31i)) - (*(*u31im1)));
                            flux[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                [k_darts_counter_temp2655][3]
                                = tx3 * ((*(*u41i)) - (*(*u41im1)));
                            flux[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                [k_darts_counter_temp2655][4]
                                = 0.5 * (1. - 1.3999999999999999 * 1.3999999999999999) * tx3
                                    * (((*(*u21i)) * (*(*u21i)) + (*(*u31i)) * (*(*u31i))
                                           + (*(*u41i)) * (*(*u41i)))
                                        - ((*(*u21im1)) * (*(*u21im1)) + (*(*u31im1)) * (*(*u31im1))
                                            + (*(*u41im1)) * (*(*u41im1))))
                                + (1. / 6.) * tx3
                                    * ((*(*u21i)) * (*(*u21i)) - (*(*u21im1)) * (*(*u21im1)))
                                + 1.3999999999999999 * 1.3999999999999999 * tx3
                                    * ((*(*u51i)) - (*(*u51im1)));
                        }
                        (*i) = i_darts_counter_temp2655;
                    }
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp2655 = (*i);
                        for (; i_darts_counter_temp2655 <= iend; i_darts_counter_temp2655++) {
                            frct[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                [k_darts_counter_temp2655][0]
                                = frct[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                      [k_darts_counter_temp2655][0]
                                + dx1 * tx1
                                    * (rsd[i_darts_counter_temp2655 - 1][(j_darts_counter_temp2655)]
                                          [k_darts_counter_temp2655][0]
                                        - 2.
                                            * rsd[i_darts_counter_temp2655]
                                                 [(j_darts_counter_temp2655)]
                                                 [k_darts_counter_temp2655][0]
                                        + rsd[i_darts_counter_temp2655 + 1]
                                             [(j_darts_counter_temp2655)][k_darts_counter_temp2655]
                                             [0]);
                            frct[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                [k_darts_counter_temp2655][1]
                                = frct[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                      [k_darts_counter_temp2655][1]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp2655 + 1]
                                           [(j_darts_counter_temp2655)][k_darts_counter_temp2655][1]
                                        - flux[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                              [k_darts_counter_temp2655][1])
                                + dx2 * tx1
                                    * (rsd[i_darts_counter_temp2655 - 1][(j_darts_counter_temp2655)]
                                          [k_darts_counter_temp2655][1]
                                        - 2.
                                            * rsd[i_darts_counter_temp2655]
                                                 [(j_darts_counter_temp2655)]
                                                 [k_darts_counter_temp2655][1]
                                        + rsd[i_darts_counter_temp2655 + 1]
                                             [(j_darts_counter_temp2655)][k_darts_counter_temp2655]
                                             [1]);
                            frct[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                [k_darts_counter_temp2655][2]
                                = frct[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                      [k_darts_counter_temp2655][2]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp2655 + 1]
                                           [(j_darts_counter_temp2655)][k_darts_counter_temp2655][2]
                                        - flux[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                              [k_darts_counter_temp2655][2])
                                + dx3 * tx1
                                    * (rsd[i_darts_counter_temp2655 - 1][(j_darts_counter_temp2655)]
                                          [k_darts_counter_temp2655][2]
                                        - 2.
                                            * rsd[i_darts_counter_temp2655]
                                                 [(j_darts_counter_temp2655)]
                                                 [k_darts_counter_temp2655][2]
                                        + rsd[i_darts_counter_temp2655 + 1]
                                             [(j_darts_counter_temp2655)][k_darts_counter_temp2655]
                                             [2]);
                            frct[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                [k_darts_counter_temp2655][3]
                                = frct[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                      [k_darts_counter_temp2655][3]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp2655 + 1]
                                           [(j_darts_counter_temp2655)][k_darts_counter_temp2655][3]
                                        - flux[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                              [k_darts_counter_temp2655][3])
                                + dx4 * tx1
                                    * (rsd[i_darts_counter_temp2655 - 1][(j_darts_counter_temp2655)]
                                          [k_darts_counter_temp2655][3]
                                        - 2.
                                            * rsd[i_darts_counter_temp2655]
                                                 [(j_darts_counter_temp2655)]
                                                 [k_darts_counter_temp2655][3]
                                        + rsd[i_darts_counter_temp2655 + 1]
                                             [(j_darts_counter_temp2655)][k_darts_counter_temp2655]
                                             [3]);
                            frct[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                [k_darts_counter_temp2655][4]
                                = frct[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                      [k_darts_counter_temp2655][4]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp2655 + 1]
                                           [(j_darts_counter_temp2655)][k_darts_counter_temp2655][4]
                                        - flux[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                              [k_darts_counter_temp2655][4])
                                + dx5 * tx1
                                    * (rsd[i_darts_counter_temp2655 - 1][(j_darts_counter_temp2655)]
                                          [k_darts_counter_temp2655][4]
                                        - 2.
                                            * rsd[i_darts_counter_temp2655]
                                                 [(j_darts_counter_temp2655)]
                                                 [k_darts_counter_temp2655][4]
                                        + rsd[i_darts_counter_temp2655 + 1]
                                             [(j_darts_counter_temp2655)][k_darts_counter_temp2655]
                                             [4]);
                        }
                        (*i) = i_darts_counter_temp2655;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp2655 = (*m);
                        for (; m_darts_counter_temp2655 < 5; m_darts_counter_temp2655++) {
                            frct[1][(j_darts_counter_temp2655)][k_darts_counter_temp2655]
                                [m_darts_counter_temp2655]
                                = frct[1][(j_darts_counter_temp2655)][k_darts_counter_temp2655]
                                      [m_darts_counter_temp2655]
                                - (*(*dsspm))
                                    * (+5.
                                            * rsd[1][(j_darts_counter_temp2655)]
                                                 [k_darts_counter_temp2655]
                                                 [m_darts_counter_temp2655]
                                        - 4.
                                            * rsd[2][(j_darts_counter_temp2655)]
                                                 [k_darts_counter_temp2655]
                                                 [m_darts_counter_temp2655]
                                        + rsd[3][(j_darts_counter_temp2655)]
                                             [k_darts_counter_temp2655][m_darts_counter_temp2655]);
                            frct[2][(j_darts_counter_temp2655)][k_darts_counter_temp2655]
                                [m_darts_counter_temp2655]
                                = frct[2][(j_darts_counter_temp2655)][k_darts_counter_temp2655]
                                      [m_darts_counter_temp2655]
                                - (*(*dsspm))
                                    * (-4.
                                            * rsd[1][(j_darts_counter_temp2655)]
                                                 [k_darts_counter_temp2655]
                                                 [m_darts_counter_temp2655]
                                        + 6.
                                            * rsd[2][(j_darts_counter_temp2655)]
                                                 [k_darts_counter_temp2655]
                                                 [m_darts_counter_temp2655]
                                        - 4.
                                            * rsd[3][(j_darts_counter_temp2655)]
                                                 [k_darts_counter_temp2655]
                                                 [m_darts_counter_temp2655]
                                        + rsd[4][(j_darts_counter_temp2655)]
                                             [k_darts_counter_temp2655][m_darts_counter_temp2655]);
                        }
                        (*m) = m_darts_counter_temp2655;
                    }
                    (*(*ist1)) = 3;
                    (*(*iend1)) = nx - 4;
                    {
                        /*Loop's init*/
                        (*i) = (*(*ist1));
                        int i_darts_counter_temp2655 = (*i);
                        for (; i_darts_counter_temp2655 <= (*(*iend1));
                             i_darts_counter_temp2655++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp2655 = (*m);
                                for (; m_darts_counter_temp2655 < 5; m_darts_counter_temp2655++) {
                                    frct[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                        [k_darts_counter_temp2655][m_darts_counter_temp2655]
                                        = frct[i_darts_counter_temp2655][(j_darts_counter_temp2655)]
                                              [k_darts_counter_temp2655][m_darts_counter_temp2655]
                                        - (*(*dsspm))
                                            * (rsd[i_darts_counter_temp2655 - 2]
                                                  [(j_darts_counter_temp2655)]
                                                  [k_darts_counter_temp2655]
                                                  [m_darts_counter_temp2655]
                                                - 4.
                                                    * rsd[i_darts_counter_temp2655 - 1]
                                                         [(j_darts_counter_temp2655)]
                                                         [k_darts_counter_temp2655]
                                                         [m_darts_counter_temp2655]
                                                + 6.
                                                    * rsd[i_darts_counter_temp2655]
                                                         [(j_darts_counter_temp2655)]
                                                         [k_darts_counter_temp2655]
                                                         [m_darts_counter_temp2655]
                                                - 4.
                                                    * rsd[i_darts_counter_temp2655 + 1]
                                                         [(j_darts_counter_temp2655)]
                                                         [k_darts_counter_temp2655]
                                                         [m_darts_counter_temp2655]
                                                + rsd[i_darts_counter_temp2655 + 2]
                                                     [(j_darts_counter_temp2655)]
                                                     [k_darts_counter_temp2655]
                                                     [m_darts_counter_temp2655]);
                                }
                                (*m) = m_darts_counter_temp2655;
                            }
                        }
                        (*i) = i_darts_counter_temp2655;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp2655 = (*m);
                        for (; m_darts_counter_temp2655 < 5; m_darts_counter_temp2655++) {
                            frct[nx - 3][(j_darts_counter_temp2655)][k_darts_counter_temp2655]
                                [m_darts_counter_temp2655]
                                = frct[nx - 3][(j_darts_counter_temp2655)][k_darts_counter_temp2655]
                                      [m_darts_counter_temp2655]
                                - (*(*dsspm))
                                    * (rsd[nx - 5][(j_darts_counter_temp2655)]
                                          [k_darts_counter_temp2655][m_darts_counter_temp2655]
                                        - 4.
                                            * rsd[nx - 4][(j_darts_counter_temp2655)]
                                                 [k_darts_counter_temp2655]
                                                 [m_darts_counter_temp2655]
                                        + 6.
                                            * rsd[nx - 3][(j_darts_counter_temp2655)]
                                                 [k_darts_counter_temp2655]
                                                 [m_darts_counter_temp2655]
                                        - 4.
                                            * rsd[nx - 2][(j_darts_counter_temp2655)]
                                                 [k_darts_counter_temp2655]
                                                 [m_darts_counter_temp2655]);
                            frct[nx - 2][(j_darts_counter_temp2655)][k_darts_counter_temp2655]
                                [m_darts_counter_temp2655]
                                = frct[nx - 2][(j_darts_counter_temp2655)][k_darts_counter_temp2655]
                                      [m_darts_counter_temp2655]
                                - (*(*dsspm))
                                    * (rsd[nx - 4][(j_darts_counter_temp2655)]
                                          [k_darts_counter_temp2655][m_darts_counter_temp2655]
                                        - 4.
                                            * rsd[nx - 3][(j_darts_counter_temp2655)]
                                                 [k_darts_counter_temp2655]
                                                 [m_darts_counter_temp2655]
                                        + 5.
                                            * rsd[nx - 2][(j_darts_counter_temp2655)]
                                                 [k_darts_counter_temp2655]
                                                 [m_darts_counter_temp2655]);
                        }
                        (*m) = m_darts_counter_temp2655;
                    }
                }
                (*k) = k_darts_counter_temp2655;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets2655[0].decDep();
}
TP2655::TP2655(int in_numThreads, int in_mainCodeletID, TP2301* in_TPParent, int in_initIteration,
    int in_lastIteration, TP2655** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L2_darts2655(new int*[this->numThreads])
    , dsspm_darts2655(new double*[this->numThreads])
    , i_darts2655(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iend1_darts2655(new int*[this->numThreads])
    , ist1_darts2655(new int*[this->numThreads])
    , j_darts2655(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts2655(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts2655(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts2655(new double*[this->numThreads])
    , u21i_darts2655(new double*[this->numThreads])
    , u21im1_darts2655(new double*[this->numThreads])
    , u31i_darts2655(new double*[this->numThreads])
    , u31im1_darts2655(new double*[this->numThreads])
    , u41i_darts2655(new double*[this->numThreads])
    , u41im1_darts2655(new double*[this->numThreads])
    , u51i_darts2655(new double*[this->numThreads])
    , u51im1_darts2655(new double*[this->numThreads])
    , initIteration2655(in_initIteration)
    , lastIteration2655(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets2655(new _barrierCodelets2655[1])
    , checkInCodelets2656(new _checkInCodelets2656[this->numThreads])
{
    /*Initialize the loop parameters*/
    range2655 = abs(lastIteration2655 - initIteration2655) / 1;
    rangePerCodelet2655 = range2655 / numThreads;
    minIteration2655 = min<int>(lastIteration2655, initIteration2655);
    remainderRange2655 = range2655 % numThreads;
    /*Initialize inputs and vars.*/
    this->L2_darts2655 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->dsspm_darts2655
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->iend1_darts2655 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ist1_darts2655 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts2655
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21i_darts2655
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21im1_darts2655
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31i_darts2655
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31im1_darts2655
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41i_darts2655
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41im1_darts2655
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51i_darts2655
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51im1_darts2655
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets2655[0] = _barrierCodelets2655(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets2656* checkInCodelets2656Ptr = (this->checkInCodelets2656);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets2656);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets2656Ptr) = _checkInCodelets2656(2, 1, this, codeletCounter);
#else
        (*checkInCodelets2656Ptr) = _checkInCodelets2656(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets2656Ptr).decDep();
        checkInCodelets2656Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP2655::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets2656[localID].setID(codeletID);
    this->checkInCodelets2656[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets2656[localID + this->baseNumThreads * i]
            = _checkInCodelets2656(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets2656[localID + this->baseNumThreads * i]
            = _checkInCodelets2656(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets2656[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets2656[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP2655::~TP2655()
{
    delete[] L2_darts2655;
    delete[] dsspm_darts2655;
    delete[] iend1_darts2655;
    delete[] ist1_darts2655;
    delete[] tmp_darts2655;
    delete[] u21i_darts2655;
    delete[] u21im1_darts2655;
    delete[] u31i_darts2655;
    delete[] u31im1_darts2655;
    delete[] u41i_darts2655;
    delete[] u41im1_darts2655;
    delete[] u51i_darts2655;
    delete[] u51im1_darts2655;
    delete[] barrierCodelets2655;
    delete[] checkInCodelets2656;
}
/*TP3293: OMPForDirective*/
void TP3293::_barrierCodelets3293::fire(void)
{
    TP3293* myTP = static_cast<TP3293*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets3293[0].decDep();
}
bool TP3293::requestNewRangeIterations3293(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet3293 * codeletID;
        int tempEndRange = rangePerCodelet3293 * (codeletID + 1);
        if (remainderRange3293 != 0) {
            if (codeletID < (uint32_t)remainderRange3293) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange3293;
                tempEndRange += remainderRange3293;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration3293;
        tempEndRange = tempEndRange * 1 + minIteration3293;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration3293 < lastIteration3293) {
            (this->inputsTPParent->i_darts3293[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts3293[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration3293;
        }
    }
    return isThereNewIteration;
}
void TP3293::_checkInCodelets3294::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L1_darts3293[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L1_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L2_darts3293[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->q_darts3293[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31_darts3293[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31_darts2301[this->getID()]);

    /*printing node 3294: ForStmt*/
    /*var: L1*/
    /*var: L2*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: q*/
    /*var: u31*/
    int** L1 = &(this->inputsTPParent->L1_darts3293[this->getLocalID()]);
    (void)L1 /*OMP_SHARED_PRIVATE*/;
    int** L2 = &(this->inputsTPParent->L2_darts3293[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts3293[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts3293[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts3293[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts3293[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** u31 = &(this->inputsTPParent->u31_darts3293[this->getLocalID()]);
    (void)u31 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations3293(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets3293[0].decDep();
        return;
    }
    for (int i_darts_counter_temp3293 = (*i); i_darts_counter_temp3293 <= endRange
         && i_darts_counter_temp3293 <= this->inputsTPParent->lastIteration3293;
         i_darts_counter_temp3293++) {
        {
            {
                /*Loop's init*/
                (*j) = (*(*L1));
                int j_darts_counter_temp3293 = (*j);
                for (; j_darts_counter_temp3293 <= (*(*L2)); j_darts_counter_temp3293++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp3293 = (*k);
                        for (; k_darts_counter_temp3293 <= nz - 2; k_darts_counter_temp3293++) {
                            flux[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                [k_darts_counter_temp3293][0]
                                = rsd[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                     [k_darts_counter_temp3293][2];
                            (*(*u31)) = rsd[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                           [k_darts_counter_temp3293][2]
                                / rsd[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                     [k_darts_counter_temp3293][0];
                            (*(*q)) = 0.5
                                * (rsd[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                      [k_darts_counter_temp3293][1]
                                        * rsd[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                             [k_darts_counter_temp3293][1]
                                    + rsd[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                         [k_darts_counter_temp3293][2]
                                        * rsd[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                             [k_darts_counter_temp3293][2]
                                    + rsd[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                         [k_darts_counter_temp3293][3]
                                        * rsd[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                             [k_darts_counter_temp3293][3])
                                / rsd[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                     [k_darts_counter_temp3293][0];
                            flux[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                [k_darts_counter_temp3293][1]
                                = rsd[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                     [k_darts_counter_temp3293][1]
                                * (*(*u31));
                            flux[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                [k_darts_counter_temp3293][2]
                                = rsd[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                     [k_darts_counter_temp3293][2]
                                    * (*(*u31))
                                + 0.40000000000000002
                                    * (rsd[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                          [k_darts_counter_temp3293][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                [k_darts_counter_temp3293][3]
                                = rsd[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                     [k_darts_counter_temp3293][3]
                                * (*(*u31));
                            flux[(i_darts_counter_temp3293)][j_darts_counter_temp3293]
                                [k_darts_counter_temp3293][4]
                                = (1.3999999999999999
                                          * rsd[(i_darts_counter_temp3293)]
                                               [j_darts_counter_temp3293][k_darts_counter_temp3293]
                                               [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u31));
                        }
                        (*k) = k_darts_counter_temp3293;
                    }
                }
                (*j) = j_darts_counter_temp3293;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets3293[0].decDep();
}
TP3293::TP3293(int in_numThreads, int in_mainCodeletID, TP2301* in_TPParent, int in_initIteration,
    int in_lastIteration, TP3293** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts3293(new int*[this->numThreads])
    , L2_darts3293(new int*[this->numThreads])
    , i_darts3293(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts3293(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts3293(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts3293(new double*[this->numThreads])
    , u31_darts3293(new double*[this->numThreads])
    , initIteration3293(in_initIteration)
    , lastIteration3293(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets3293(new _barrierCodelets3293[1])
    , checkInCodelets3294(new _checkInCodelets3294[this->numThreads])
{
    /*Initialize the loop parameters*/
    range3293 = abs(lastIteration3293 - initIteration3293) / 1;
    rangePerCodelet3293 = range3293 / numThreads;
    minIteration3293 = min<int>(lastIteration3293, initIteration3293);
    remainderRange3293 = range3293 % numThreads;
    /*Initialize inputs and vars.*/
    this->L1_darts3293 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->L2_darts3293 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->q_darts3293 = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31_darts3293
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets3293[0] = _barrierCodelets3293(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets3294* checkInCodelets3294Ptr = (this->checkInCodelets3294);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets3294);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets3294Ptr) = _checkInCodelets3294(2, 1, this, codeletCounter);
#else
        (*checkInCodelets3294Ptr) = _checkInCodelets3294(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets3294Ptr).decDep();
        checkInCodelets3294Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP3293::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets3294[localID].setID(codeletID);
    this->checkInCodelets3294[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets3294[localID + this->baseNumThreads * i]
            = _checkInCodelets3294(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets3294[localID + this->baseNumThreads * i]
            = _checkInCodelets3294(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets3294[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets3294[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP3293::~TP3293()
{
    delete[] L1_darts3293;
    delete[] L2_darts3293;
    delete[] q_darts3293;
    delete[] u31_darts3293;
    delete[] barrierCodelets3293;
    delete[] checkInCodelets3294;
}
/*TP3442: OMPForDirective*/
void TP3442::_barrierCodelets3442::fire(void)
{
    TP3442* myTP = static_cast<TP3442*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets3442[0].decDep();
}
bool TP3442::requestNewRangeIterations3442(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet3442 * codeletID;
        int tempEndRange = rangePerCodelet3442 * (codeletID + 1);
        if (remainderRange3442 != 0) {
            if (codeletID < (uint32_t)remainderRange3442) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange3442;
                tempEndRange += remainderRange3442;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration3442;
        tempEndRange = tempEndRange * 1 + minIteration3442;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration3442 < lastIteration3442) {
            (this->inputsTPParent->i_darts3442[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts3442[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration3442;
        }
    }
    return isThereNewIteration;
}
void TP3442::_checkInCodelets3443::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L2_darts3442[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->dsspm_darts3442[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->dsspm_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jend1_darts3442[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jend1_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jst1_darts3442[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jst1_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts3442[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21j_darts3442[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21j_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21jm1_darts3442[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21jm1_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31j_darts3442[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31j_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31jm1_darts3442[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31jm1_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41j_darts3442[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41j_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41jm1_darts3442[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41jm1_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51j_darts3442[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51j_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51jm1_darts3442[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51jm1_darts2301[this->getID()]);

    /*printing node 3443: ForStmt*/
    /*var: L2*/
    /*var: dsspm*/
    /*var: i*/
    /*var: j*/
    /*var: jend1*/
    /*var: jst1*/
    /*var: k*/
    /*var: m*/
    /*var: tmp*/
    /*var: u21j*/
    /*var: u21jm1*/
    /*var: u31j*/
    /*var: u31jm1*/
    /*var: u41j*/
    /*var: u41jm1*/
    /*var: u51j*/
    /*var: u51jm1*/
    int** L2 = &(this->inputsTPParent->L2_darts3442[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    double** dsspm = &(this->inputsTPParent->dsspm_darts3442[this->getLocalID()]);
    (void)dsspm /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts3442[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts3442[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jend1 = &(this->inputsTPParent->jend1_darts3442[this->getLocalID()]);
    (void)jend1 /*OMP_SHARED_PRIVATE*/;
    int** jst1 = &(this->inputsTPParent->jst1_darts3442[this->getLocalID()]);
    (void)jst1 /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts3442[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts3442[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts3442[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21j = &(this->inputsTPParent->u21j_darts3442[this->getLocalID()]);
    (void)u21j /*OMP_SHARED_PRIVATE*/;
    double** u21jm1 = &(this->inputsTPParent->u21jm1_darts3442[this->getLocalID()]);
    (void)u21jm1 /*OMP_SHARED_PRIVATE*/;
    double** u31j = &(this->inputsTPParent->u31j_darts3442[this->getLocalID()]);
    (void)u31j /*OMP_SHARED_PRIVATE*/;
    double** u31jm1 = &(this->inputsTPParent->u31jm1_darts3442[this->getLocalID()]);
    (void)u31jm1 /*OMP_SHARED_PRIVATE*/;
    double** u41j = &(this->inputsTPParent->u41j_darts3442[this->getLocalID()]);
    (void)u41j /*OMP_SHARED_PRIVATE*/;
    double** u41jm1 = &(this->inputsTPParent->u41jm1_darts3442[this->getLocalID()]);
    (void)u41jm1 /*OMP_SHARED_PRIVATE*/;
    double** u51j = &(this->inputsTPParent->u51j_darts3442[this->getLocalID()]);
    (void)u51j /*OMP_SHARED_PRIVATE*/;
    double** u51jm1 = &(this->inputsTPParent->u51jm1_darts3442[this->getLocalID()]);
    (void)u51jm1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations3442(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets3442[0].decDep();
        return;
    }
    for (int i_darts_counter_temp3442 = (*i); i_darts_counter_temp3442 <= endRange
         && i_darts_counter_temp3442 <= this->inputsTPParent->lastIteration3442;
         i_darts_counter_temp3442++) {
        {
            {
                /*Loop's init*/
                (*k) = 1;
                int k_darts_counter_temp3442 = (*k);
                for (; k_darts_counter_temp3442 <= nz - 2; k_darts_counter_temp3442++) {
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp3442 = (*j);
                        for (; j_darts_counter_temp3442 <= jend; j_darts_counter_temp3442++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp3442 = (*m);
                                for (; m_darts_counter_temp3442 < 5; m_darts_counter_temp3442++) {
                                    frct[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                        [k_darts_counter_temp3442][m_darts_counter_temp3442]
                                        = frct[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                              [k_darts_counter_temp3442][m_darts_counter_temp3442]
                                        - ty2
                                            * (flux[(i_darts_counter_temp3442)]
                                                   [j_darts_counter_temp3442 + 1]
                                                   [k_darts_counter_temp3442]
                                                   [m_darts_counter_temp3442]
                                                - flux[(i_darts_counter_temp3442)]
                                                      [j_darts_counter_temp3442 - 1]
                                                      [k_darts_counter_temp3442]
                                                      [m_darts_counter_temp3442]);
                                }
                                (*m) = m_darts_counter_temp3442;
                            }
                        }
                        (*j) = j_darts_counter_temp3442;
                    }
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp3442 = (*j);
                        for (; j_darts_counter_temp3442 <= (*(*L2)); j_darts_counter_temp3442++) {
                            (*(*tmp)) = 1.
                                / rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                     [k_darts_counter_temp3442][0];
                            (*(*u21j)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                     [k_darts_counter_temp3442][1];
                            (*(*u31j)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                     [k_darts_counter_temp3442][2];
                            (*(*u41j)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                     [k_darts_counter_temp3442][3];
                            (*(*u51j)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                     [k_darts_counter_temp3442][4];
                            (*(*tmp)) = 1.
                                / rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442 - 1]
                                     [k_darts_counter_temp3442][0];
                            (*(*u21jm1)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442 - 1]
                                     [k_darts_counter_temp3442][1];
                            (*(*u31jm1)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442 - 1]
                                     [k_darts_counter_temp3442][2];
                            (*(*u41jm1)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442 - 1]
                                     [k_darts_counter_temp3442][3];
                            (*(*u51jm1)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442 - 1]
                                     [k_darts_counter_temp3442][4];
                            flux[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                [k_darts_counter_temp3442][1]
                                = ty3 * ((*(*u21j)) - (*(*u21jm1)));
                            flux[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                [k_darts_counter_temp3442][2]
                                = (4. / 3.) * ty3 * ((*(*u31j)) - (*(*u31jm1)));
                            flux[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                [k_darts_counter_temp3442][3]
                                = ty3 * ((*(*u41j)) - (*(*u41jm1)));
                            flux[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                [k_darts_counter_temp3442][4]
                                = 0.5 * (1. - 1.3999999999999999 * 1.3999999999999999) * ty3
                                    * (((*(*u21j)) * (*(*u21j)) + (*(*u31j)) * (*(*u31j))
                                           + (*(*u41j)) * (*(*u41j)))
                                        - ((*(*u21jm1)) * (*(*u21jm1)) + (*(*u31jm1)) * (*(*u31jm1))
                                            + (*(*u41jm1)) * (*(*u41jm1))))
                                + (1. / 6.) * ty3
                                    * ((*(*u31j)) * (*(*u31j)) - (*(*u31jm1)) * (*(*u31jm1)))
                                + 1.3999999999999999 * 1.3999999999999999 * ty3
                                    * ((*(*u51j)) - (*(*u51jm1)));
                        }
                        (*j) = j_darts_counter_temp3442;
                    }
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp3442 = (*j);
                        for (; j_darts_counter_temp3442 <= jend; j_darts_counter_temp3442++) {
                            frct[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                [k_darts_counter_temp3442][0]
                                = frct[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                      [k_darts_counter_temp3442][0]
                                + dy1 * ty1
                                    * (rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442 - 1]
                                          [k_darts_counter_temp3442][0]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3442)]
                                                 [j_darts_counter_temp3442]
                                                 [k_darts_counter_temp3442][0]
                                        + rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442
                                            + 1][k_darts_counter_temp3442][0]);
                            frct[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                [k_darts_counter_temp3442][1]
                                = frct[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                      [k_darts_counter_temp3442][1]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp3442)][j_darts_counter_temp3442
                                           + 1][k_darts_counter_temp3442][1]
                                        - flux[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                              [k_darts_counter_temp3442][1])
                                + dy2 * ty1
                                    * (rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442 - 1]
                                          [k_darts_counter_temp3442][1]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3442)]
                                                 [j_darts_counter_temp3442]
                                                 [k_darts_counter_temp3442][1]
                                        + rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442
                                            + 1][k_darts_counter_temp3442][1]);
                            frct[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                [k_darts_counter_temp3442][2]
                                = frct[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                      [k_darts_counter_temp3442][2]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp3442)][j_darts_counter_temp3442
                                           + 1][k_darts_counter_temp3442][2]
                                        - flux[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                              [k_darts_counter_temp3442][2])
                                + dy3 * ty1
                                    * (rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442 - 1]
                                          [k_darts_counter_temp3442][2]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3442)]
                                                 [j_darts_counter_temp3442]
                                                 [k_darts_counter_temp3442][2]
                                        + rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442
                                            + 1][k_darts_counter_temp3442][2]);
                            frct[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                [k_darts_counter_temp3442][3]
                                = frct[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                      [k_darts_counter_temp3442][3]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp3442)][j_darts_counter_temp3442
                                           + 1][k_darts_counter_temp3442][3]
                                        - flux[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                              [k_darts_counter_temp3442][3])
                                + dy4 * ty1
                                    * (rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442 - 1]
                                          [k_darts_counter_temp3442][3]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3442)]
                                                 [j_darts_counter_temp3442]
                                                 [k_darts_counter_temp3442][3]
                                        + rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442
                                            + 1][k_darts_counter_temp3442][3]);
                            frct[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                [k_darts_counter_temp3442][4]
                                = frct[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                      [k_darts_counter_temp3442][4]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp3442)][j_darts_counter_temp3442
                                           + 1][k_darts_counter_temp3442][4]
                                        - flux[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                              [k_darts_counter_temp3442][4])
                                + dy5 * ty1
                                    * (rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442 - 1]
                                          [k_darts_counter_temp3442][4]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3442)]
                                                 [j_darts_counter_temp3442]
                                                 [k_darts_counter_temp3442][4]
                                        + rsd[(i_darts_counter_temp3442)][j_darts_counter_temp3442
                                            + 1][k_darts_counter_temp3442][4]);
                        }
                        (*j) = j_darts_counter_temp3442;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp3442 = (*m);
                        for (; m_darts_counter_temp3442 < 5; m_darts_counter_temp3442++) {
                            frct[(i_darts_counter_temp3442)][1][k_darts_counter_temp3442]
                                [m_darts_counter_temp3442]
                                = frct[(i_darts_counter_temp3442)][1][k_darts_counter_temp3442]
                                      [m_darts_counter_temp3442]
                                - (*(*dsspm))
                                    * (+5.
                                            * rsd[(i_darts_counter_temp3442)][1]
                                                 [k_darts_counter_temp3442]
                                                 [m_darts_counter_temp3442]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3442)][2]
                                                 [k_darts_counter_temp3442]
                                                 [m_darts_counter_temp3442]
                                        + rsd[(i_darts_counter_temp3442)][3]
                                             [k_darts_counter_temp3442][m_darts_counter_temp3442]);
                            frct[(i_darts_counter_temp3442)][2][k_darts_counter_temp3442]
                                [m_darts_counter_temp3442]
                                = frct[(i_darts_counter_temp3442)][2][k_darts_counter_temp3442]
                                      [m_darts_counter_temp3442]
                                - (*(*dsspm))
                                    * (-4.
                                            * rsd[(i_darts_counter_temp3442)][1]
                                                 [k_darts_counter_temp3442]
                                                 [m_darts_counter_temp3442]
                                        + 6.
                                            * rsd[(i_darts_counter_temp3442)][2]
                                                 [k_darts_counter_temp3442]
                                                 [m_darts_counter_temp3442]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3442)][3]
                                                 [k_darts_counter_temp3442]
                                                 [m_darts_counter_temp3442]
                                        + rsd[(i_darts_counter_temp3442)][4]
                                             [k_darts_counter_temp3442][m_darts_counter_temp3442]);
                        }
                        (*m) = m_darts_counter_temp3442;
                    }
                    (*(*jst1)) = 3;
                    (*(*jend1)) = ny - 4;
                    {
                        /*Loop's init*/
                        (*j) = (*(*jst1));
                        int j_darts_counter_temp3442 = (*j);
                        for (; j_darts_counter_temp3442 <= (*(*jend1));
                             j_darts_counter_temp3442++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp3442 = (*m);
                                for (; m_darts_counter_temp3442 < 5; m_darts_counter_temp3442++) {
                                    frct[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                        [k_darts_counter_temp3442][m_darts_counter_temp3442]
                                        = frct[(i_darts_counter_temp3442)][j_darts_counter_temp3442]
                                              [k_darts_counter_temp3442][m_darts_counter_temp3442]
                                        - (*(*dsspm))
                                            * (rsd[(i_darts_counter_temp3442)]
                                                  [j_darts_counter_temp3442 - 2]
                                                  [k_darts_counter_temp3442]
                                                  [m_darts_counter_temp3442]
                                                - 4.
                                                    * rsd[(i_darts_counter_temp3442)]
                                                         [j_darts_counter_temp3442 - 1]
                                                         [k_darts_counter_temp3442]
                                                         [m_darts_counter_temp3442]
                                                + 6.
                                                    * rsd[(i_darts_counter_temp3442)]
                                                         [j_darts_counter_temp3442]
                                                         [k_darts_counter_temp3442]
                                                         [m_darts_counter_temp3442]
                                                - 4.
                                                    * rsd[(i_darts_counter_temp3442)]
                                                         [j_darts_counter_temp3442 + 1]
                                                         [k_darts_counter_temp3442]
                                                         [m_darts_counter_temp3442]
                                                + rsd[(i_darts_counter_temp3442)]
                                                     [j_darts_counter_temp3442 + 2]
                                                     [k_darts_counter_temp3442]
                                                     [m_darts_counter_temp3442]);
                                }
                                (*m) = m_darts_counter_temp3442;
                            }
                        }
                        (*j) = j_darts_counter_temp3442;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp3442 = (*m);
                        for (; m_darts_counter_temp3442 < 5; m_darts_counter_temp3442++) {
                            frct[(i_darts_counter_temp3442)][ny - 3][k_darts_counter_temp3442]
                                [m_darts_counter_temp3442]
                                = frct[(i_darts_counter_temp3442)][ny - 3][k_darts_counter_temp3442]
                                      [m_darts_counter_temp3442]
                                - (*(*dsspm))
                                    * (rsd[(i_darts_counter_temp3442)][ny - 5]
                                          [k_darts_counter_temp3442][m_darts_counter_temp3442]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3442)][ny - 4]
                                                 [k_darts_counter_temp3442]
                                                 [m_darts_counter_temp3442]
                                        + 6.
                                            * rsd[(i_darts_counter_temp3442)][ny - 3]
                                                 [k_darts_counter_temp3442]
                                                 [m_darts_counter_temp3442]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3442)][ny - 2]
                                                 [k_darts_counter_temp3442]
                                                 [m_darts_counter_temp3442]);
                            frct[(i_darts_counter_temp3442)][ny - 2][k_darts_counter_temp3442]
                                [m_darts_counter_temp3442]
                                = frct[(i_darts_counter_temp3442)][ny - 2][k_darts_counter_temp3442]
                                      [m_darts_counter_temp3442]
                                - (*(*dsspm))
                                    * (rsd[(i_darts_counter_temp3442)][ny - 4]
                                          [k_darts_counter_temp3442][m_darts_counter_temp3442]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3442)][ny - 3]
                                                 [k_darts_counter_temp3442]
                                                 [m_darts_counter_temp3442]
                                        + 5.
                                            * rsd[(i_darts_counter_temp3442)][ny - 2]
                                                 [k_darts_counter_temp3442]
                                                 [m_darts_counter_temp3442]);
                        }
                        (*m) = m_darts_counter_temp3442;
                    }
                }
                (*k) = k_darts_counter_temp3442;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets3442[0].decDep();
}
TP3442::TP3442(int in_numThreads, int in_mainCodeletID, TP2301* in_TPParent, int in_initIteration,
    int in_lastIteration, TP3442** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L2_darts3442(new int*[this->numThreads])
    , dsspm_darts3442(new double*[this->numThreads])
    , i_darts3442(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts3442(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jend1_darts3442(new int*[this->numThreads])
    , jst1_darts3442(new int*[this->numThreads])
    , k_darts3442(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts3442(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts3442(new double*[this->numThreads])
    , u21j_darts3442(new double*[this->numThreads])
    , u21jm1_darts3442(new double*[this->numThreads])
    , u31j_darts3442(new double*[this->numThreads])
    , u31jm1_darts3442(new double*[this->numThreads])
    , u41j_darts3442(new double*[this->numThreads])
    , u41jm1_darts3442(new double*[this->numThreads])
    , u51j_darts3442(new double*[this->numThreads])
    , u51jm1_darts3442(new double*[this->numThreads])
    , initIteration3442(in_initIteration)
    , lastIteration3442(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets3442(new _barrierCodelets3442[1])
    , checkInCodelets3443(new _checkInCodelets3443[this->numThreads])
{
    /*Initialize the loop parameters*/
    range3442 = abs(lastIteration3442 - initIteration3442) / 1;
    rangePerCodelet3442 = range3442 / numThreads;
    minIteration3442 = min<int>(lastIteration3442, initIteration3442);
    remainderRange3442 = range3442 % numThreads;
    /*Initialize inputs and vars.*/
    this->L2_darts3442 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->dsspm_darts3442
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jend1_darts3442 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jst1_darts3442 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts3442
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21j_darts3442
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21jm1_darts3442
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31j_darts3442
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31jm1_darts3442
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41j_darts3442
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41jm1_darts3442
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51j_darts3442
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51jm1_darts3442
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets3442[0] = _barrierCodelets3442(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets3443* checkInCodelets3443Ptr = (this->checkInCodelets3443);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets3443);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets3443Ptr) = _checkInCodelets3443(2, 1, this, codeletCounter);
#else
        (*checkInCodelets3443Ptr) = _checkInCodelets3443(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets3443Ptr).decDep();
        checkInCodelets3443Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP3442::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets3443[localID].setID(codeletID);
    this->checkInCodelets3443[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets3443[localID + this->baseNumThreads * i]
            = _checkInCodelets3443(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets3443[localID + this->baseNumThreads * i]
            = _checkInCodelets3443(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets3443[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets3443[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP3442::~TP3442()
{
    delete[] L2_darts3442;
    delete[] dsspm_darts3442;
    delete[] jend1_darts3442;
    delete[] jst1_darts3442;
    delete[] tmp_darts3442;
    delete[] u21j_darts3442;
    delete[] u21jm1_darts3442;
    delete[] u31j_darts3442;
    delete[] u31jm1_darts3442;
    delete[] u41j_darts3442;
    delete[] u41jm1_darts3442;
    delete[] u51j_darts3442;
    delete[] u51jm1_darts3442;
    delete[] barrierCodelets3442;
    delete[] checkInCodelets3443;
}
/*TP4077: OMPForDirective*/
void TP4077::_barrierCodelets4077::fire(void)
{
    TP4077* myTP = static_cast<TP4077*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets4077[0].decDep();
}
bool TP4077::requestNewRangeIterations4077(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet4077 * codeletID;
        int tempEndRange = rangePerCodelet4077 * (codeletID + 1);
        if (remainderRange4077 != 0) {
            if (codeletID < (uint32_t)remainderRange4077) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange4077;
                tempEndRange += remainderRange4077;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration4077;
        tempEndRange = tempEndRange * 1 + minIteration4077;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration4077 < lastIteration4077) {
            (this->inputsTPParent->i_darts4077[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts4077[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration4077;
        }
    }
    return isThereNewIteration;
}
void TP4077::_checkInCodelets4078::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->dsspm_darts4077[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->dsspm_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->q_darts4077[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts4077[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21k_darts4077[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21k_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21km1_darts4077[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21km1_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31k_darts4077[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31k_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31km1_darts4077[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31km1_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41_darts4077[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41k_darts4077[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41k_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41km1_darts4077[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41km1_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51k_darts4077[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51k_darts2301[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51km1_darts4077[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51km1_darts2301[this->getID()]);

    /*printing node 4078: ForStmt*/
    /*var: dsspm*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    /*var: q*/
    /*var: tmp*/
    /*var: u21k*/
    /*var: u21km1*/
    /*var: u31k*/
    /*var: u31km1*/
    /*var: u41*/
    /*var: u41k*/
    /*var: u41km1*/
    /*var: u51k*/
    /*var: u51km1*/
    double** dsspm = &(this->inputsTPParent->dsspm_darts4077[this->getLocalID()]);
    (void)dsspm /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts4077[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts4077[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts4077[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts4077[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts4077[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts4077[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21k = &(this->inputsTPParent->u21k_darts4077[this->getLocalID()]);
    (void)u21k /*OMP_SHARED_PRIVATE*/;
    double** u21km1 = &(this->inputsTPParent->u21km1_darts4077[this->getLocalID()]);
    (void)u21km1 /*OMP_SHARED_PRIVATE*/;
    double** u31k = &(this->inputsTPParent->u31k_darts4077[this->getLocalID()]);
    (void)u31k /*OMP_SHARED_PRIVATE*/;
    double** u31km1 = &(this->inputsTPParent->u31km1_darts4077[this->getLocalID()]);
    (void)u31km1 /*OMP_SHARED_PRIVATE*/;
    double** u41 = &(this->inputsTPParent->u41_darts4077[this->getLocalID()]);
    (void)u41 /*OMP_SHARED_PRIVATE*/;
    double** u41k = &(this->inputsTPParent->u41k_darts4077[this->getLocalID()]);
    (void)u41k /*OMP_SHARED_PRIVATE*/;
    double** u41km1 = &(this->inputsTPParent->u41km1_darts4077[this->getLocalID()]);
    (void)u41km1 /*OMP_SHARED_PRIVATE*/;
    double** u51k = &(this->inputsTPParent->u51k_darts4077[this->getLocalID()]);
    (void)u51k /*OMP_SHARED_PRIVATE*/;
    double** u51km1 = &(this->inputsTPParent->u51km1_darts4077[this->getLocalID()]);
    (void)u51km1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations4077(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets4077[0].decDep();
        return;
    }
    for (int i_darts_counter_temp4077 = (*i); i_darts_counter_temp4077 <= endRange
         && i_darts_counter_temp4077 <= this->inputsTPParent->lastIteration4077;
         i_darts_counter_temp4077++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp4077 = (*j);
                for (; j_darts_counter_temp4077 <= jend; j_darts_counter_temp4077++) {
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp4077 = (*k);
                        for (; k_darts_counter_temp4077 <= nz - 1; k_darts_counter_temp4077++) {
                            flux[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                [k_darts_counter_temp4077][0]
                                = rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                     [k_darts_counter_temp4077][3];
                            (*(*u41)) = rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                           [k_darts_counter_temp4077][3]
                                / rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                     [k_darts_counter_temp4077][0];
                            (*(*q)) = 0.5
                                * (rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                      [k_darts_counter_temp4077][1]
                                        * rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                             [k_darts_counter_temp4077][1]
                                    + rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                         [k_darts_counter_temp4077][2]
                                        * rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                             [k_darts_counter_temp4077][2]
                                    + rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                         [k_darts_counter_temp4077][3]
                                        * rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                             [k_darts_counter_temp4077][3])
                                / rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                     [k_darts_counter_temp4077][0];
                            flux[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                [k_darts_counter_temp4077][1]
                                = rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                     [k_darts_counter_temp4077][1]
                                * (*(*u41));
                            flux[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                [k_darts_counter_temp4077][2]
                                = rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                     [k_darts_counter_temp4077][2]
                                * (*(*u41));
                            flux[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                [k_darts_counter_temp4077][3]
                                = rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                     [k_darts_counter_temp4077][3]
                                    * (*(*u41))
                                + 0.40000000000000002
                                    * (rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                          [k_darts_counter_temp4077][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                [k_darts_counter_temp4077][4]
                                = (1.3999999999999999
                                          * rsd[(i_darts_counter_temp4077)]
                                               [j_darts_counter_temp4077][k_darts_counter_temp4077]
                                               [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u41));
                        }
                        (*k) = k_darts_counter_temp4077;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp4077 = (*k);
                        for (; k_darts_counter_temp4077 <= nz - 2; k_darts_counter_temp4077++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp4077 = (*m);
                                for (; m_darts_counter_temp4077 < 5; m_darts_counter_temp4077++) {
                                    frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                        [k_darts_counter_temp4077][m_darts_counter_temp4077]
                                        = frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                              [k_darts_counter_temp4077][m_darts_counter_temp4077]
                                        - tz2
                                            * (flux[(i_darts_counter_temp4077)]
                                                   [j_darts_counter_temp4077]
                                                   [k_darts_counter_temp4077 + 1]
                                                   [m_darts_counter_temp4077]
                                                - flux[k_darts_counter_temp4077 - 1]
                                                      [(i_darts_counter_temp4077)]
                                                      [j_darts_counter_temp4077]
                                                      [m_darts_counter_temp4077]);
                                }
                                (*m) = m_darts_counter_temp4077;
                            }
                        }
                        (*k) = k_darts_counter_temp4077;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp4077 = (*k);
                        for (; k_darts_counter_temp4077 <= nz - 1; k_darts_counter_temp4077++) {
                            (*(*tmp)) = 1.
                                / rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                     [k_darts_counter_temp4077][0];
                            (*(*u21k)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                     [k_darts_counter_temp4077][1];
                            (*(*u31k)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                     [k_darts_counter_temp4077][2];
                            (*(*u41k)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                     [k_darts_counter_temp4077][3];
                            (*(*u51k)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                     [k_darts_counter_temp4077][4];
                            (*(*tmp)) = 1.
                                / rsd[k_darts_counter_temp4077 - 1][(i_darts_counter_temp4077)]
                                     [j_darts_counter_temp4077][0];
                            (*(*u21km1)) = (*(*tmp))
                                * rsd[k_darts_counter_temp4077 - 1][(i_darts_counter_temp4077)]
                                     [j_darts_counter_temp4077][1];
                            (*(*u31km1)) = (*(*tmp))
                                * rsd[k_darts_counter_temp4077 - 1][(i_darts_counter_temp4077)]
                                     [j_darts_counter_temp4077][2];
                            (*(*u41km1)) = (*(*tmp))
                                * rsd[k_darts_counter_temp4077 - 1][(i_darts_counter_temp4077)]
                                     [j_darts_counter_temp4077][3];
                            (*(*u51km1)) = (*(*tmp))
                                * rsd[k_darts_counter_temp4077 - 1][(i_darts_counter_temp4077)]
                                     [j_darts_counter_temp4077][4];
                            flux[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                [k_darts_counter_temp4077][1]
                                = tz3 * ((*(*u21k)) - (*(*u21km1)));
                            flux[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                [k_darts_counter_temp4077][2]
                                = tz3 * ((*(*u31k)) - (*(*u31km1)));
                            flux[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                [k_darts_counter_temp4077][3]
                                = (4. / 3.) * tz3 * ((*(*u41k)) - (*(*u41km1)));
                            flux[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                [k_darts_counter_temp4077][4]
                                = 0.5 * (1. - 1.3999999999999999 * 1.3999999999999999) * tz3
                                    * (((*(*u21k)) * (*(*u21k)) + (*(*u31k)) * (*(*u31k))
                                           + (*(*u41k)) * (*(*u41k)))
                                        - ((*(*u21km1)) * (*(*u21km1)) + (*(*u31km1)) * (*(*u31km1))
                                            + (*(*u41km1)) * (*(*u41km1))))
                                + (1. / 6.) * tz3
                                    * ((*(*u41k)) * (*(*u41k)) - (*(*u41km1)) * (*(*u41km1)))
                                + 1.3999999999999999 * 1.3999999999999999 * tz3
                                    * ((*(*u51k)) - (*(*u51km1)));
                        }
                        (*k) = k_darts_counter_temp4077;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp4077 = (*k);
                        for (; k_darts_counter_temp4077 <= nz - 2; k_darts_counter_temp4077++) {
                            frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                [k_darts_counter_temp4077][0]
                                = frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                      [k_darts_counter_temp4077][0]
                                + dz1 * tz1
                                    * (rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                          [k_darts_counter_temp4077 + 1][0]
                                        - 2.
                                            * rsd[(i_darts_counter_temp4077)]
                                                 [j_darts_counter_temp4077]
                                                 [k_darts_counter_temp4077][0]
                                        + rsd[k_darts_counter_temp4077 - 1]
                                             [(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                             [0]);
                            frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                [k_darts_counter_temp4077][1]
                                = frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                      [k_darts_counter_temp4077][1]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                           [k_darts_counter_temp4077 + 1][1]
                                        - flux[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                              [k_darts_counter_temp4077][1])
                                + dz2 * tz1
                                    * (rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                          [k_darts_counter_temp4077 + 1][1]
                                        - 2.
                                            * rsd[(i_darts_counter_temp4077)]
                                                 [j_darts_counter_temp4077]
                                                 [k_darts_counter_temp4077][1]
                                        + rsd[k_darts_counter_temp4077 - 1]
                                             [(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                             [1]);
                            frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                [k_darts_counter_temp4077][2]
                                = frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                      [k_darts_counter_temp4077][2]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                           [k_darts_counter_temp4077 + 1][2]
                                        - flux[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                              [k_darts_counter_temp4077][2])
                                + dz3 * tz1
                                    * (rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                          [k_darts_counter_temp4077 + 1][2]
                                        - 2.
                                            * rsd[(i_darts_counter_temp4077)]
                                                 [j_darts_counter_temp4077]
                                                 [k_darts_counter_temp4077][2]
                                        + rsd[k_darts_counter_temp4077 - 1]
                                             [(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                             [2]);
                            frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                [k_darts_counter_temp4077][3]
                                = frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                      [k_darts_counter_temp4077][3]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                           [k_darts_counter_temp4077 + 1][3]
                                        - flux[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                              [k_darts_counter_temp4077][3])
                                + dz4 * tz1
                                    * (rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                          [k_darts_counter_temp4077 + 1][3]
                                        - 2.
                                            * rsd[(i_darts_counter_temp4077)]
                                                 [j_darts_counter_temp4077]
                                                 [k_darts_counter_temp4077][3]
                                        + rsd[k_darts_counter_temp4077 - 1]
                                             [(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                             [3]);
                            frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                [k_darts_counter_temp4077][4]
                                = frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                      [k_darts_counter_temp4077][4]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                           [k_darts_counter_temp4077 + 1][4]
                                        - flux[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                              [k_darts_counter_temp4077][4])
                                + dz5 * tz1
                                    * (rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                          [k_darts_counter_temp4077 + 1][4]
                                        - 2.
                                            * rsd[(i_darts_counter_temp4077)]
                                                 [j_darts_counter_temp4077]
                                                 [k_darts_counter_temp4077][4]
                                        + rsd[k_darts_counter_temp4077 - 1]
                                             [(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                             [4]);
                        }
                        (*k) = k_darts_counter_temp4077;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp4077 = (*m);
                        for (; m_darts_counter_temp4077 < 5; m_darts_counter_temp4077++) {
                            frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077][1]
                                [m_darts_counter_temp4077]
                                = frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077][1]
                                      [m_darts_counter_temp4077]
                                - (*(*dsspm))
                                    * (+5.
                                            * rsd[(i_darts_counter_temp4077)]
                                                 [j_darts_counter_temp4077][1]
                                                 [m_darts_counter_temp4077]
                                        - 4.
                                            * rsd[(i_darts_counter_temp4077)]
                                                 [j_darts_counter_temp4077][2]
                                                 [m_darts_counter_temp4077]
                                        + rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                             [3][m_darts_counter_temp4077]);
                            frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077][2]
                                [m_darts_counter_temp4077]
                                = frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077][2]
                                      [m_darts_counter_temp4077]
                                - (*(*dsspm))
                                    * (-4.
                                            * rsd[(i_darts_counter_temp4077)]
                                                 [j_darts_counter_temp4077][1]
                                                 [m_darts_counter_temp4077]
                                        + 6.
                                            * rsd[(i_darts_counter_temp4077)]
                                                 [j_darts_counter_temp4077][2]
                                                 [m_darts_counter_temp4077]
                                        - 4.
                                            * rsd[(i_darts_counter_temp4077)]
                                                 [j_darts_counter_temp4077][3]
                                                 [m_darts_counter_temp4077]
                                        + rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                             [4][m_darts_counter_temp4077]);
                        }
                        (*m) = m_darts_counter_temp4077;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 3;
                        int k_darts_counter_temp4077 = (*k);
                        for (; k_darts_counter_temp4077 <= nz - 4; k_darts_counter_temp4077++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp4077 = (*m);
                                for (; m_darts_counter_temp4077 < 5; m_darts_counter_temp4077++) {
                                    frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                        [k_darts_counter_temp4077][m_darts_counter_temp4077]
                                        = frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                              [k_darts_counter_temp4077][m_darts_counter_temp4077]
                                        - (*(*dsspm))
                                            * (rsd[(i_darts_counter_temp4077)]
                                                  [j_darts_counter_temp4077]
                                                  [k_darts_counter_temp4077 - 2]
                                                  [m_darts_counter_temp4077]
                                                - 4.
                                                    * rsd[k_darts_counter_temp4077 - 1]
                                                         [(i_darts_counter_temp4077)]
                                                         [j_darts_counter_temp4077]
                                                         [m_darts_counter_temp4077]
                                                + 6.
                                                    * rsd[(i_darts_counter_temp4077)]
                                                         [j_darts_counter_temp4077]
                                                         [k_darts_counter_temp4077]
                                                         [m_darts_counter_temp4077]
                                                - 4.
                                                    * rsd[(i_darts_counter_temp4077)]
                                                         [j_darts_counter_temp4077]
                                                         [k_darts_counter_temp4077 + 1]
                                                         [m_darts_counter_temp4077]
                                                + rsd[(i_darts_counter_temp4077)]
                                                     [j_darts_counter_temp4077]
                                                     [k_darts_counter_temp4077 + 2]
                                                     [m_darts_counter_temp4077]);
                                }
                                (*m) = m_darts_counter_temp4077;
                            }
                        }
                        (*k) = k_darts_counter_temp4077;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp4077 = (*m);
                        for (; m_darts_counter_temp4077 < 5; m_darts_counter_temp4077++) {
                            frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077][nz - 3]
                                [m_darts_counter_temp4077]
                                = frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077][nz - 3]
                                      [m_darts_counter_temp4077]
                                - (*(*dsspm))
                                    * (rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                          [nz - 5][m_darts_counter_temp4077]
                                        - 4.
                                            * rsd[(i_darts_counter_temp4077)]
                                                 [j_darts_counter_temp4077][nz - 4]
                                                 [m_darts_counter_temp4077]
                                        + 6.
                                            * rsd[(i_darts_counter_temp4077)]
                                                 [j_darts_counter_temp4077][nz - 3]
                                                 [m_darts_counter_temp4077]
                                        - 4.
                                            * rsd[(i_darts_counter_temp4077)]
                                                 [j_darts_counter_temp4077][nz - 2]
                                                 [m_darts_counter_temp4077]);
                            frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077][nz - 2]
                                [m_darts_counter_temp4077]
                                = frct[(i_darts_counter_temp4077)][j_darts_counter_temp4077][nz - 2]
                                      [m_darts_counter_temp4077]
                                - (*(*dsspm))
                                    * (rsd[(i_darts_counter_temp4077)][j_darts_counter_temp4077]
                                          [nz - 4][m_darts_counter_temp4077]
                                        - 4.
                                            * rsd[(i_darts_counter_temp4077)]
                                                 [j_darts_counter_temp4077][nz - 3]
                                                 [m_darts_counter_temp4077]
                                        + 5.
                                            * rsd[(i_darts_counter_temp4077)]
                                                 [j_darts_counter_temp4077][nz - 2]
                                                 [m_darts_counter_temp4077]);
                        }
                        (*m) = m_darts_counter_temp4077;
                    }
                }
                (*j) = j_darts_counter_temp4077;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets4077[0].decDep();
}
TP4077::TP4077(int in_numThreads, int in_mainCodeletID, TP2301* in_TPParent, int in_initIteration,
    int in_lastIteration, TP4077** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , dsspm_darts4077(new double*[this->numThreads])
    , i_darts4077(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts4077(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts4077(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts4077(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts4077(new double*[this->numThreads])
    , tmp_darts4077(new double*[this->numThreads])
    , u21k_darts4077(new double*[this->numThreads])
    , u21km1_darts4077(new double*[this->numThreads])
    , u31k_darts4077(new double*[this->numThreads])
    , u31km1_darts4077(new double*[this->numThreads])
    , u41_darts4077(new double*[this->numThreads])
    , u41k_darts4077(new double*[this->numThreads])
    , u41km1_darts4077(new double*[this->numThreads])
    , u51k_darts4077(new double*[this->numThreads])
    , u51km1_darts4077(new double*[this->numThreads])
    , initIteration4077(in_initIteration)
    , lastIteration4077(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets4077(new _barrierCodelets4077[1])
    , checkInCodelets4078(new _checkInCodelets4078[this->numThreads])
{
    /*Initialize the loop parameters*/
    range4077 = abs(lastIteration4077 - initIteration4077) / 1;
    rangePerCodelet4077 = range4077 / numThreads;
    minIteration4077 = min<int>(lastIteration4077, initIteration4077);
    remainderRange4077 = range4077 % numThreads;
    /*Initialize inputs and vars.*/
    this->dsspm_darts4077
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->q_darts4077 = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts4077
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21k_darts4077
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21km1_darts4077
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31k_darts4077
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31km1_darts4077
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41_darts4077
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41k_darts4077
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41km1_darts4077
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51k_darts4077
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51km1_darts4077
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets4077[0] = _barrierCodelets4077(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets4078* checkInCodelets4078Ptr = (this->checkInCodelets4078);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets4078);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets4078Ptr) = _checkInCodelets4078(2, 1, this, codeletCounter);
#else
        (*checkInCodelets4078Ptr) = _checkInCodelets4078(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets4078Ptr).decDep();
        checkInCodelets4078Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP4077::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets4078[localID].setID(codeletID);
    this->checkInCodelets4078[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets4078[localID + this->baseNumThreads * i]
            = _checkInCodelets4078(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets4078[localID + this->baseNumThreads * i]
            = _checkInCodelets4078(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets4078[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets4078[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP4077::~TP4077()
{
    delete[] dsspm_darts4077;
    delete[] q_darts4077;
    delete[] tmp_darts4077;
    delete[] u21k_darts4077;
    delete[] u21km1_darts4077;
    delete[] u31k_darts4077;
    delete[] u31km1_darts4077;
    delete[] u41_darts4077;
    delete[] u41k_darts4077;
    delete[] u41km1_darts4077;
    delete[] u51k_darts4077;
    delete[] u51km1_darts4077;
    delete[] barrierCodelets4077;
    delete[] checkInCodelets4078;
}
/*TP7: TP_jacld*/
void TP7::_checkInCodelets4982::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 4982: DeclStmt*/

    /*printing node 4983: DeclStmt*/

    /*printing node 4984: DeclStmt*/

    /*printing node 4985: DeclStmt*/

    /*printing node 4986: DeclStmt*/

    /*printing node 4987: BinaryOperator*/
    (this->inputsTPParent->r43_darts7[this->getID()]) = (4. / 3.);

    /*printing node 4991: BinaryOperator*/
    (this->inputsTPParent->c1345_darts7[this->getID()])
        = 1.3999999999999999 * 0.10000000000000001 * 1. * 1.3999999999999999;

    /*printing node 4999: BinaryOperator*/
    (this->inputsTPParent->c34_darts7[this->getID()]) = 0.10000000000000001 * 1.;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 4982 nextRegion: 5003 */
    myTP->controlTPParent->checkInCodelets5003[this->getID()].decDep();
}
void TP7::_checkInCodelets5003::fire(void)
{
    /*region 5003 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP5003;
    if (idx < myTP->TPsToUse5003) {
        if (!__sync_val_compare_and_swap(&(myTP->TP5003_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse5003;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse5003;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (ist < iend) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse5003 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse5003 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP5003>(myTP, myTP->codeletsPerTP5003 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP5003Ptr[idx]));
#else
            place<TP5003>(idx, myTP, myTP->codeletsPerTP5003 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP5003Ptr[idx]));
#endif
        } else {
            if (myTP->TP5003Ptr[idx] != nullptr) {
                myTP->TP5003Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    } else {
        /*Find and signal the next codelet*/

        myTP->controlTPParent->nextCodeletsjacld[this->getID()]->decDep();
    }
}
TP7::TP7(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
    darts::Codelet* in_mainSyncCodelet, TP7** in_ptrToThisFunctionTP, int in_k)
    : ompTP(in_numThreads, in_mainCodeletID)
    , ptrToThisFunctionTP(in_ptrToThisFunctionTP)
    , inputsTPParent(this)
    , controlTPParent(this)
    , nextCodeletsjacld(new Codelet*[in_numThreads])
    , nextSyncCodeletsjacld(new Codelet*[in_numThreads])
    , k_darts7(new int[this->numThreads])
    , c1345_darts7(new double[this->numThreads])
    , c34_darts7(new double[this->numThreads])
    , i_darts7(new int[this->numThreads])
    , j_darts7(new int[this->numThreads])
    , r43_darts7(new double[this->numThreads])
    , tmp1_darts7(new double[this->numThreads])
    , tmp2_darts7(new double[this->numThreads])
    , tmp3_darts7(new double[this->numThreads])
    , TP5003Ptr(new TP5003*[NUMTPS5003])
    , TP5003_alreadyLaunched(new size_t[NUMTPS5003])
    , numTPsSet5003(0)
    , numTPsReady5003(0)
    , TPsToUse5003(NUMTPS5003)
    , codeletsPerTP5003(this->numThreads / NUMTPS5003)
    , totalCodelets5003(this->TPsToUse5003 * this->codeletsPerTP5003)
    , checkInCodelets4982(new _checkInCodelets4982[this->numThreads])
    , checkInCodelets5003(new _checkInCodelets5003[this->numThreads])
{
    _checkInCodelets5003* checkInCodelets5003Ptr = (this->checkInCodelets5003);
    for (int i = 0; i < NUMTPS5003; i++) {
        TP5003Ptr[i] = nullptr;
        TP5003_alreadyLaunched[i] = 0;
    }
    _checkInCodelets4982* checkInCodelets4982Ptr = (this->checkInCodelets4982);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets4982);
#endif
    for (size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++) {
        (*checkInCodelets5003Ptr) = _checkInCodelets5003(1, 1, this, codeletCounter);
        checkInCodelets5003Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets4982Ptr) = _checkInCodelets4982(2, 1, this, codeletCounter);
#else
        (*checkInCodelets4982Ptr) = _checkInCodelets4982(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets4982Ptr).decDep();
        checkInCodelets4982Ptr++;
    }
    if (this->numThreads == 1) {
        this->nextCodeletsjacld[0] = in_mainNextCodelet;
        this->nextSyncCodeletsjacld[0] = in_mainSyncCodelet;
        this->k_darts7[0] = in_k;
        this->availableCodelets[0] = 1;
    } else {
        this->k_darts7[this->mainCodeletID] = in_k;
        this->nextCodeletsjacld[in_mainCodeletID] = in_mainNextCodelet;
        this->nextSyncCodeletsjacld[in_mainCodeletID] = in_mainSyncCodelet;
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[this->mainCodeletID].decDep();
#else
        this->availableCodelets[this->mainCodeletID] = 1;
#endif
        *(this->ptrToThisFunctionTP) = this;
    }
}
TP7::~TP7()
{
    delete[] checkInCodelets5003;
    delete[] checkInCodelets4982;
    delete[] nextCodeletsjacld;
    delete[] nextSyncCodeletsjacld;
    delete[] k_darts7;
    delete[] c1345_darts7;
    delete[] c34_darts7;
    delete[] i_darts7;
    delete[] j_darts7;
    delete[] r43_darts7;
    delete[] tmp1_darts7;
    delete[] tmp2_darts7;
    delete[] tmp3_darts7;
}
void TP7::setNewInputs(int in_k, size_t codeletID) { this->k_darts7[codeletID] = in_k; }
/*TP5003: OMPForDirective*/
bool TP5003::requestNewRangeIterations5003(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 1*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet5003 * codeletID;
        int tempEndRange = rangePerCodelet5003 * (codeletID + 1);
        if (remainderRange5003 != 0) {
            if (codeletID < (uint32_t)remainderRange5003) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange5003;
                tempEndRange += remainderRange5003;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration5003;
        tempEndRange = tempEndRange * 1 + minIteration5003;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration5003 < lastIteration5003) {
            (this->inputsTPParent->i_darts5003[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts5003[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration5003;
        }
    }
    return isThereNewIteration;
}
void TP5003::_checkInCodelets5004::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->c1345_darts5003[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->c1345_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->c34_darts5003[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->c34_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->k_darts5003[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->k_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->r43_darts5003[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->r43_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp1_darts5003[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp1_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp2_darts5003[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp2_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp3_darts5003[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp3_darts7[this->getID()]);

    /*printing node 5004: ForStmt*/
    /*var: c1345*/
    /*var: c34*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: r43*/
    /*var: tmp1*/
    /*var: tmp2*/
    /*var: tmp3*/
    double** c1345 = &(this->inputsTPParent->c1345_darts5003[this->getLocalID()]);
    (void)c1345 /*OMP_SHARED_PRIVATE*/;
    double** c34 = &(this->inputsTPParent->c34_darts5003[this->getLocalID()]);
    (void)c34 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts5003[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts5003[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** k = &(this->inputsTPParent->k_darts5003[this->getLocalID()]);
    (void)k /*OMP_SHARED_PRIVATE*/;
    double** r43 = &(this->inputsTPParent->r43_darts5003[this->getLocalID()]);
    (void)r43 /*OMP_SHARED_PRIVATE*/;
    double** tmp1 = &(this->inputsTPParent->tmp1_darts5003[this->getLocalID()]);
    (void)tmp1 /*OMP_SHARED_PRIVATE*/;
    double** tmp2 = &(this->inputsTPParent->tmp2_darts5003[this->getLocalID()]);
    (void)tmp2 /*OMP_SHARED_PRIVATE*/;
    double** tmp3 = &(this->inputsTPParent->tmp3_darts5003[this->getLocalID()]);
    (void)tmp3 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations5003(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Find and signal the next codelet*/
        return;
    }
    for (int i_darts_counter_temp5003 = (*i); i_darts_counter_temp5003 <= endRange
         && i_darts_counter_temp5003 <= this->inputsTPParent->lastIteration5003;
         i_darts_counter_temp5003++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp5003 = (*j);
                for (; j_darts_counter_temp5003 <= jend; j_darts_counter_temp5003++) {
                    (*(*tmp1))
                        = 1. / u[(i_darts_counter_temp5003)][j_darts_counter_temp5003][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][0]
                        = 1. + dt * 2. * (tx1 * dx1 + ty1 * dy1 + tz1 * dz1);
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][1] = 0.;
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][2] = 0.;
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][3] = 0.;
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][4] = 0.;
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][0] = dt * 2.
                        * (tx1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                       [(*(*k))][1])
                            + ty1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                       [(*(*k))][1])
                            + tz1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                       [(*(*k))][1]));
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][1] = 1.
                        + dt * 2.
                            * (tx1 * (*(*r43)) * (*(*c34)) * (*(*tmp1))
                                + ty1 * (*(*c34)) * (*(*tmp1)) + tz1 * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx2 + ty1 * dy2 + tz1 * dz2);
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][2] = 0.;
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][3] = 0.;
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][4] = 0.;
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][0] = dt * 2.
                        * (tx1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                       [(*(*k))][2])
                            + ty1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                       [(*(*k))][2])
                            + tz1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                       [(*(*k))][2]));
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][1] = 0.;
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][2] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c34)) * (*(*tmp1))
                                + ty1 * (*(*r43)) * (*(*c34)) * (*(*tmp1))
                                + tz1 * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx3 + ty1 * dy3 + tz1 * dz3);
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][3] = 0.;
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][4] = 0.;
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][0] = dt * 2.
                        * (tx1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                       [(*(*k))][3])
                            + ty1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                       [(*(*k))][3])
                            + tz1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                       [(*(*k))][3]));
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][1] = 0.;
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][2] = 0.;
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][3] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c34)) * (*(*tmp1)) + ty1 * (*(*c34)) * (*(*tmp1))
                                + tz1 * (*(*r43)) * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx4 + ty1 * dy4 + tz1 * dz4);
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][4] = 0.;
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][0] = dt * 2.
                        * (tx1
                                * (-((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp5003)]
                                                [j_darts_counter_temp5003][(*(*k))][1])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp5003)]
                                                [j_darts_counter_temp5003][(*(*k))][2])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp5003)]
                                                [j_darts_counter_temp5003][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                           [(*(*k))][4])
                            + ty1
                                * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp5003)]
                                                [j_darts_counter_temp5003][(*(*k))][1])))
                                    - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp5003)]
                                                [j_darts_counter_temp5003][(*(*k))][2])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp5003)]
                                                [j_darts_counter_temp5003][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                           [(*(*k))][4])
                            + tz1
                                * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp5003)]
                                                [j_darts_counter_temp5003][(*(*k))][1])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp5003)]
                                                [j_darts_counter_temp5003][(*(*k))][2])))
                                    - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp5003)]
                                                [j_darts_counter_temp5003][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                           [(*(*k))][4]));
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][1] = dt * 2.
                        * (tx1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003][(*(*k))]
                                   [1]
                            + ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003][(*(*k))]
                                   [1]
                            + tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003][(*(*k))]
                                   [1]);
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][2] = dt * 2.
                        * (tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003][(*(*k))]
                                   [2]
                            + ty1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003][(*(*k))]
                                   [2]
                            + tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003][(*(*k))]
                                   [2]);
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][3] = dt * 2.
                        * (tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003][(*(*k))]
                                   [3]
                            + ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003][(*(*k))]
                                   [3]
                            + tz1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003][(*(*k))]
                                   [3]);
                    d[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][4] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c1345)) * (*(*tmp1)) + ty1 * (*(*c1345)) * (*(*tmp1))
                                + tz1 * (*(*c1345)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx5 + ty1 * dy5 + tz1 * dz5);
                    (*(*tmp1)) = 1.
                        / u[(*(*k)) - 1][(i_darts_counter_temp5003)][j_darts_counter_temp5003][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][0] = -dt * tz1 * dz1;
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][1] = 0.;
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][2] = 0.;
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][3] = -dt * tz2;
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][4] = 0.;
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][0] = -dt * tz2
                            * (-(u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                  [j_darts_counter_temp5003][1]
                                   * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                      [j_darts_counter_temp5003][3])
                                * (*(*tmp2)))
                        - dt * tz1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                   [j_darts_counter_temp5003][1]);
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][1] = -dt * tz2
                            * (u[(*(*k)) - 1][(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * (*(*c34)) * (*(*tmp1)) - dt * tz1 * dz2;
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][2] = 0.;
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][3] = -dt * tz2
                        * (u[(*(*k)) - 1][(i_darts_counter_temp5003)][j_darts_counter_temp5003][1]
                            * (*(*tmp1)));
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][4] = 0.;
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][0] = -dt * tz2
                            * (-(u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                  [j_darts_counter_temp5003][2]
                                   * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                      [j_darts_counter_temp5003][3])
                                * (*(*tmp2)))
                        - dt * tz1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                   [j_darts_counter_temp5003][2]);
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][1] = 0.;
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][2] = -dt * tz2
                            * (u[(*(*k)) - 1][(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * ((*(*c34)) * (*(*tmp1))) - dt * tz1 * dz3;
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][3] = -dt * tz2
                        * (u[(*(*k)) - 1][(i_darts_counter_temp5003)][j_darts_counter_temp5003][2]
                            * (*(*tmp1)));
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][4] = 0.;
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][0] = -dt * tz2
                            * (-(u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                  [j_darts_counter_temp5003][3]
                                   * (*(*tmp1)))
                                    * (u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                        [j_darts_counter_temp5003][3]
                                        * (*(*tmp1)))
                                + 0.5 * 0.40000000000000002
                                    * ((u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                         [j_darts_counter_temp5003][1]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003][1]
                                           + u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                              [j_darts_counter_temp5003][2]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003][2]
                                           + u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                              [j_darts_counter_temp5003][3]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003][3])
                                        * (*(*tmp2))))
                        - dt * tz1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                   [j_darts_counter_temp5003][3]);
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][1] = -dt * tz2
                        * (-0.40000000000000002
                            * (u[(*(*k)) - 1][(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                [1]
                                * (*(*tmp1))));
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][2] = -dt * tz2
                        * (-0.40000000000000002
                            * (u[(*(*k)) - 1][(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                [2]
                                * (*(*tmp1))));
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][3] = -dt * tz2
                            * (2. - 0.40000000000000002)
                            * (u[(*(*k)) - 1][(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * tz1 * dz4;
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][4]
                        = -dt * tz2 * 0.40000000000000002;
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][0] = -dt * tz2
                            * ((0.40000000000000002
                                       * (u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                           [j_darts_counter_temp5003][1]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003][1]
                                           + u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                              [j_darts_counter_temp5003][2]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003][2]
                                           + u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                              [j_darts_counter_temp5003][3]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                           [j_darts_counter_temp5003][4]
                                           * (*(*tmp1))))
                                * (u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                    [j_darts_counter_temp5003][3]
                                    * (*(*tmp1))))
                        - dt * tz1
                            * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                        [j_darts_counter_temp5003][1]
                                        * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                           [j_darts_counter_temp5003][1])
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                        [j_darts_counter_temp5003][2]
                                        * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                           [j_darts_counter_temp5003][2])
                                - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                        [j_darts_counter_temp5003][3]
                                        * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                           [j_darts_counter_temp5003][3])
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                       [j_darts_counter_temp5003][4]);
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][1] = -dt * tz2
                            * (-0.40000000000000002
                                * (u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                    [j_darts_counter_temp5003][1]
                                    * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                       [j_darts_counter_temp5003][3])
                                * (*(*tmp2)))
                        - dt * tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(*(*k)) - 1][(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                               [1];
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][2] = -dt * tz2
                            * (-0.40000000000000002
                                * (u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                    [j_darts_counter_temp5003][2]
                                    * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                       [j_darts_counter_temp5003][3])
                                * (*(*tmp2)))
                        - dt * tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(*(*k)) - 1][(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                               [2];
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][3] = -dt * tz2
                            * (1.3999999999999999
                                    * (u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                        [j_darts_counter_temp5003][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                         [j_darts_counter_temp5003][1]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003][1]
                                           + u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                              [j_darts_counter_temp5003][2]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003][2]
                                           + 3.
                                               * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003][3]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003][3])
                                        * (*(*tmp2))))
                        - dt * tz1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(*(*k)) - 1][(i_darts_counter_temp5003)][j_darts_counter_temp5003]
                               [3];
                    a[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][4] = -dt * tz2
                            * (1.3999999999999999
                                * (u[(*(*k)) - 1][(i_darts_counter_temp5003)]
                                    [j_darts_counter_temp5003][3]
                                    * (*(*tmp1))))
                        - dt * tz1 * (*(*c1345)) * (*(*tmp1)) - dt * tz1 * dz5;
                    (*(*tmp1)) = 1.
                        / u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][0] = -dt * ty1 * dy1;
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][1] = 0.;
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][2] = -dt * ty2;
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][3] = 0.;
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][4] = 0.;
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][0] = -dt * ty2
                            * (-(u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                  [(*(*k))][1]
                                   * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                      [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * ty1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                   [(*(*k))][1]);
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][1] = -dt * ty2
                            * (u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1][(*(*k))]
                                [2]
                                * (*(*tmp1)))
                        - dt * ty1 * ((*(*c34)) * (*(*tmp1))) - dt * ty1 * dy2;
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][2] = -dt * ty2
                        * (u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1][(*(*k))][1]
                            * (*(*tmp1)));
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][3] = 0.;
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][4] = 0.;
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][0] = -dt * ty2
                            * (-(u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                  [(*(*k))][2]
                                   * (*(*tmp1)))
                                    * (u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                        [(*(*k))][2]
                                        * (*(*tmp1)))
                                + 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                         [(*(*k))][1]
                                               * u[(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003 - 1][(*(*k))][1]
                                           + u[(i_darts_counter_temp5003)]
                                              [j_darts_counter_temp5003 - 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003 - 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp5003)]
                                              [j_darts_counter_temp5003 - 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003 - 1][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * ty1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                   [(*(*k))][2]);
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][1] = -dt * ty2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1][(*(*k))]
                                [1]
                                * (*(*tmp1))));
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][2] = -dt * ty2
                            * ((2. - 0.40000000000000002)
                                * (u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * ty1 * dy3;
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][3] = -dt * ty2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1][(*(*k))]
                                [3]
                                * (*(*tmp1))));
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][4]
                        = -dt * ty2 * 0.40000000000000002;
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][0] = -dt * ty2
                            * (-(u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                  [(*(*k))][2]
                                   * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                      [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * ty1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                   [(*(*k))][3]);
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][1] = 0.;
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][2] = -dt * ty2
                        * (u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1][(*(*k))][3]
                            * (*(*tmp1)));
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][3] = -dt * ty2
                            * (u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1][(*(*k))]
                                [2]
                                * (*(*tmp1)))
                        - dt * ty1 * ((*(*c34)) * (*(*tmp1))) - dt * ty1 * dy4;
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][4] = 0.;
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][0] = -dt * ty2
                            * ((0.40000000000000002
                                       * (u[(i_darts_counter_temp5003)]
                                           [j_darts_counter_temp5003 - 1][(*(*k))][1]
                                               * u[(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003 - 1][(*(*k))][1]
                                           + u[(i_darts_counter_temp5003)]
                                              [j_darts_counter_temp5003 - 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003 - 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp5003)]
                                              [j_darts_counter_temp5003 - 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003 - 1][(*(*k))][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(i_darts_counter_temp5003)]
                                           [j_darts_counter_temp5003 - 1][(*(*k))][4]
                                           * (*(*tmp1))))
                                * (u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1
                            * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                          [(*(*k))][1])
                                        * (u[(i_darts_counter_temp5003)]
                                            [j_darts_counter_temp5003 - 1][(*(*k))][1])))
                                - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                          [(*(*k))][2])
                                        * (u[(i_darts_counter_temp5003)]
                                            [j_darts_counter_temp5003 - 1][(*(*k))][2])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                          [(*(*k))][3])
                                        * (u[(i_darts_counter_temp5003)]
                                            [j_darts_counter_temp5003 - 1][(*(*k))][3])))
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                       [(*(*k))][4]);
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][1] = -dt * ty2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                    [(*(*k))][1]
                                    * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                       [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1][(*(*k))]
                               [1];
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][2] = -dt * ty2
                            * (1.3999999999999999
                                    * (u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                        [(*(*k))][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                         [(*(*k))][1]
                                               * u[(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003 - 1][(*(*k))][1]
                                           + 3.
                                               * u[(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003 - 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003 - 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp5003)]
                                              [j_darts_counter_temp5003 - 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp5003)]
                                                  [j_darts_counter_temp5003 - 1][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * ty1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1][(*(*k))]
                               [2];
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][3] = -dt * ty2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                    [(*(*k))][2]
                                    * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                       [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1][(*(*k))]
                               [3];
                    b[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][4] = -dt * ty2
                            * (1.3999999999999999
                                * (u[(i_darts_counter_temp5003)][j_darts_counter_temp5003 - 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1 * (*(*c1345)) * (*(*tmp1)) - dt * ty1 * dy5;
                    (*(*tmp1)) = 1.
                        / u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][0] = -dt * tx1 * dx1;
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][1] = -dt * tx2;
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][2] = 0.;
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][3] = 0.;
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][0][4] = 0.;
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][0] = -dt * tx2
                            * (-(u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003][(*(*k))]
                                  [1]
                                   * (*(*tmp1)))
                                    * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                        [(*(*k))][1]
                                        * (*(*tmp1)))
                                + 0.40000000000000002 * 0.5
                                    * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                        [(*(*k))][1]
                                            * u[(i_darts_counter_temp5003)-1]
                                               [j_darts_counter_temp5003][(*(*k))][1]
                                        + u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                           [(*(*k))][2]
                                            * u[(i_darts_counter_temp5003)-1]
                                               [j_darts_counter_temp5003][(*(*k))][2]
                                        + u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                           [(*(*k))][3]
                                            * u[(i_darts_counter_temp5003)-1]
                                               [j_darts_counter_temp5003][(*(*k))][3])
                                    * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003][(*(*k))]
                                   [1]);
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][1] = -dt * tx2
                            * ((2. - 0.40000000000000002)
                                * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * tx1 * dx2;
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][2] = -dt * tx2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003][(*(*k))][2]
                                * (*(*tmp1))));
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][3] = -dt * tx2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003][(*(*k))][3]
                                * (*(*tmp1))));
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][1][4]
                        = -dt * tx2 * 0.40000000000000002;
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][0] = -dt * tx2
                            * (-(u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003][(*(*k))]
                                  [1]
                                   * u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                      [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003][(*(*k))]
                                   [2]);
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][1] = -dt * tx2
                        * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003][(*(*k))][2]
                            * (*(*tmp1)));
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][2] = -dt * tx2
                            * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003][(*(*k))][1]
                                * (*(*tmp1)))
                        - dt * tx1 * ((*(*c34)) * (*(*tmp1))) - dt * tx1 * dx3;
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][3] = 0.;
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][2][4] = 0.;
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][0] = -dt * tx2
                            * (-(u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003][(*(*k))]
                                  [1]
                                   * u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                      [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003][(*(*k))]
                                   [3]);
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][1] = -dt * tx2
                        * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003][(*(*k))][3]
                            * (*(*tmp1)));
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][2] = 0.;
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][3] = -dt * tx2
                            * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003][(*(*k))][1]
                                * (*(*tmp1)))
                        - dt * tx1 * ((*(*c34)) * (*(*tmp1))) - dt * tx1 * dx4;
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][3][4] = 0.;
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][0] = -dt * tx2
                            * ((0.40000000000000002
                                       * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                           [(*(*k))][1]
                                               * u[(i_darts_counter_temp5003)-1]
                                                  [j_darts_counter_temp5003][(*(*k))][1]
                                           + u[(i_darts_counter_temp5003)-1]
                                              [j_darts_counter_temp5003][(*(*k))][2]
                                               * u[(i_darts_counter_temp5003)-1]
                                                  [j_darts_counter_temp5003][(*(*k))][2]
                                           + u[(i_darts_counter_temp5003)-1]
                                              [j_darts_counter_temp5003][(*(*k))][3]
                                               * u[(i_darts_counter_temp5003)-1]
                                                  [j_darts_counter_temp5003][(*(*k))][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                           [(*(*k))][4]
                                           * (*(*tmp1))))
                                * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1
                            * (-((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                          [(*(*k))][1])
                                        * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                            [(*(*k))][1])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                          [(*(*k))][2])
                                        * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                            [(*(*k))][2])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                          [(*(*k))][3])
                                        * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                            [(*(*k))][3])))
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                       [(*(*k))][4]);
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][1] = -dt * tx2
                            * (1.3999999999999999
                                    * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                        [(*(*k))][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((3.
                                               * u[(i_darts_counter_temp5003)-1]
                                                  [j_darts_counter_temp5003][(*(*k))][1]
                                               * u[(i_darts_counter_temp5003)-1]
                                                  [j_darts_counter_temp5003][(*(*k))][1]
                                           + u[(i_darts_counter_temp5003)-1]
                                              [j_darts_counter_temp5003][(*(*k))][2]
                                               * u[(i_darts_counter_temp5003)-1]
                                                  [j_darts_counter_temp5003][(*(*k))][2]
                                           + u[(i_darts_counter_temp5003)-1]
                                              [j_darts_counter_temp5003][(*(*k))][3]
                                               * u[(i_darts_counter_temp5003)-1]
                                                  [j_darts_counter_temp5003][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * tx1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003][(*(*k))][1];
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][2] = -dt * tx2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                    [(*(*k))][2]
                                    * u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                       [(*(*k))][1])
                                * (*(*tmp2)))
                        - dt * tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003][(*(*k))][2];
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][3] = -dt * tx2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                    [(*(*k))][3]
                                    * u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                       [(*(*k))][1])
                                * (*(*tmp2)))
                        - dt * tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003][(*(*k))][3];
                    c[(i_darts_counter_temp5003)][j_darts_counter_temp5003][4][4] = -dt * tx2
                            * (1.3999999999999999
                                * (u[(i_darts_counter_temp5003)-1][j_darts_counter_temp5003]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1 * (*(*c1345)) * (*(*tmp1)) - dt * tx1 * dx5;
                }
                (*j) = j_darts_counter_temp5003;
            }
        }
    }
    /*If this omp for has no barrier,
    check if all the codelets
    replicated from the same
    global ID has finished and
    signal the next codelet.
    Otherwise, return.*/
    uint32_t completedMultCodelet = __sync_fetch_and_add(
        &(myTP->signalNextReady[this->getLocalID() % myTP->baseNumThreads]), 1);
    if (completedMultCodelet < (uint32_t)(DARTS_CODELETS_MULT - 1))
        return;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Find and signal the next codelet*/
}
TP5003::TP5003(int in_numThreads, int in_mainCodeletID, TP7* in_TPParent, int in_initIteration,
    int in_lastIteration, TP5003** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , c1345_darts5003(new double*[this->numThreads])
    , c34_darts5003(new double*[this->numThreads])
    , i_darts5003(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts5003(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts5003(new int*[this->numThreads])
    , r43_darts5003(new double*[this->numThreads])
    , tmp1_darts5003(new double*[this->numThreads])
    , tmp2_darts5003(new double*[this->numThreads])
    , tmp3_darts5003(new double*[this->numThreads])
    , initIteration5003(in_initIteration)
    , lastIteration5003(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , signalNextReady(new int[baseNumThreads])
    , checkInCodelets5004(new _checkInCodelets5004[this->numThreads])
{
    /*Initialize the loop parameters*/
    range5003 = abs(lastIteration5003 - initIteration5003) / 1;
    rangePerCodelet5003 = range5003 / numThreads;
    minIteration5003 = min<int>(lastIteration5003, initIteration5003);
    remainderRange5003 = range5003 % numThreads;
    /*Initialize inputs and vars.*/
    this->c1345_darts5003
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->c34_darts5003
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->k_darts5003 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->r43_darts5003
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp1_darts5003
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp2_darts5003
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp3_darts5003
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    _checkInCodelets5004* checkInCodelets5004Ptr = (this->checkInCodelets5004);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets5004);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
        this->signalNextReady[codeletCounter] = 0;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets5004Ptr) = _checkInCodelets5004(2, 1, this, codeletCounter);
#else
        (*checkInCodelets5004Ptr) = _checkInCodelets5004(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets5004Ptr).decDep();
        checkInCodelets5004Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP5003::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets5004[localID].setID(codeletID);
    this->checkInCodelets5004[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets5004[localID + this->baseNumThreads * i]
            = _checkInCodelets5004(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets5004[localID + this->baseNumThreads * i]
            = _checkInCodelets5004(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets5004[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets5004[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP5003::~TP5003()
{
    delete[] c1345_darts5003;
    delete[] c34_darts5003;
    delete[] k_darts5003;
    delete[] r43_darts5003;
    delete[] tmp1_darts5003;
    delete[] tmp2_darts5003;
    delete[] tmp3_darts5003;
    delete[] checkInCodelets5004;
}
/*TP8: TP_jacu*/
void TP8::_checkInCodelets7501::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 7501: DeclStmt*/

    /*printing node 7502: DeclStmt*/

    /*printing node 7503: DeclStmt*/

    /*printing node 7504: DeclStmt*/

    /*printing node 7505: DeclStmt*/

    /*printing node 7506: BinaryOperator*/
    (this->inputsTPParent->r43_darts8[this->getID()]) = (4. / 3.);

    /*printing node 7510: BinaryOperator*/
    (this->inputsTPParent->c1345_darts8[this->getID()])
        = 1.3999999999999999 * 0.10000000000000001 * 1. * 1.3999999999999999;

    /*printing node 7518: BinaryOperator*/
    (this->inputsTPParent->c34_darts8[this->getID()]) = 0.10000000000000001 * 1.;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 7501 nextRegion: 7522 */
    myTP->controlTPParent->checkInCodelets7522[this->getID()].decDep();
}
void TP8::_checkInCodelets7522::fire(void)
{
    /*region 7522 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP7522;
    if (idx < myTP->TPsToUse7522) {
        if (!__sync_val_compare_and_swap(&(myTP->TP7522_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(ist - iend) / 1;
            int rangePerCodelet = range / myTP->TPsToUse7522;
            int minIteration = min<int>(ist, iend);
            int remainderRange = range % myTP->TPsToUse7522;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (iend < ist) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == 0) {
                lastIteration = lastIteration - 1;
            }
            if (idx == myTP->TPsToUse7522 - 1) {
                lastIteration = ist;
            }
#if USEINVOKE == 1
            invoke<TP7522>(myTP, myTP->codeletsPerTP7522 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP7522Ptr[idx]));
#else
            place<TP7522>(idx, myTP, myTP->codeletsPerTP7522 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP7522Ptr[idx]));
#endif
        } else {
            if (myTP->TP7522Ptr[idx] != nullptr) {
                myTP->TP7522Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    } else {
        /*Find and signal the next codelet*/

        myTP->controlTPParent->nextCodeletsjacu[this->getID()]->decDep();
    }
}
TP8::TP8(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
    darts::Codelet* in_mainSyncCodelet, TP8** in_ptrToThisFunctionTP, int in_k)
    : ompTP(in_numThreads, in_mainCodeletID)
    , ptrToThisFunctionTP(in_ptrToThisFunctionTP)
    , inputsTPParent(this)
    , controlTPParent(this)
    , nextCodeletsjacu(new Codelet*[in_numThreads])
    , nextSyncCodeletsjacu(new Codelet*[in_numThreads])
    , k_darts8(new int[this->numThreads])
    , c1345_darts8(new double[this->numThreads])
    , c34_darts8(new double[this->numThreads])
    , i_darts8(new int[this->numThreads])
    , j_darts8(new int[this->numThreads])
    , r43_darts8(new double[this->numThreads])
    , tmp1_darts8(new double[this->numThreads])
    , tmp2_darts8(new double[this->numThreads])
    , tmp3_darts8(new double[this->numThreads])
    , TP7522Ptr(new TP7522*[NUMTPS7522])
    , TP7522_alreadyLaunched(new size_t[NUMTPS7522])
    , numTPsSet7522(0)
    , numTPsReady7522(0)
    , TPsToUse7522(NUMTPS7522)
    , codeletsPerTP7522(this->numThreads / NUMTPS7522)
    , totalCodelets7522(this->TPsToUse7522 * this->codeletsPerTP7522)
    , checkInCodelets7501(new _checkInCodelets7501[this->numThreads])
    , checkInCodelets7522(new _checkInCodelets7522[this->numThreads])
{
    _checkInCodelets7522* checkInCodelets7522Ptr = (this->checkInCodelets7522);
    for (int i = 0; i < NUMTPS7522; i++) {
        TP7522Ptr[i] = nullptr;
        TP7522_alreadyLaunched[i] = 0;
    }
    _checkInCodelets7501* checkInCodelets7501Ptr = (this->checkInCodelets7501);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets7501);
#endif
    for (size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++) {
        (*checkInCodelets7522Ptr) = _checkInCodelets7522(1, 1, this, codeletCounter);
        checkInCodelets7522Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets7501Ptr) = _checkInCodelets7501(2, 1, this, codeletCounter);
#else
        (*checkInCodelets7501Ptr) = _checkInCodelets7501(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets7501Ptr).decDep();
        checkInCodelets7501Ptr++;
    }
    if (this->numThreads == 1) {
        this->nextCodeletsjacu[0] = in_mainNextCodelet;
        this->nextSyncCodeletsjacu[0] = in_mainSyncCodelet;
        this->k_darts8[0] = in_k;
        this->availableCodelets[0] = 1;
    } else {
        this->k_darts8[this->mainCodeletID] = in_k;
        this->nextCodeletsjacu[in_mainCodeletID] = in_mainNextCodelet;
        this->nextSyncCodeletsjacu[in_mainCodeletID] = in_mainSyncCodelet;
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[this->mainCodeletID].decDep();
#else
        this->availableCodelets[this->mainCodeletID] = 1;
#endif
        *(this->ptrToThisFunctionTP) = this;
    }
}
TP8::~TP8()
{
    delete[] checkInCodelets7522;
    delete[] checkInCodelets7501;
    delete[] nextCodeletsjacu;
    delete[] nextSyncCodeletsjacu;
    delete[] k_darts8;
    delete[] c1345_darts8;
    delete[] c34_darts8;
    delete[] i_darts8;
    delete[] j_darts8;
    delete[] r43_darts8;
    delete[] tmp1_darts8;
    delete[] tmp2_darts8;
    delete[] tmp3_darts8;
}
void TP8::setNewInputs(int in_k, size_t codeletID) { this->k_darts8[codeletID] = in_k; }
/*TP7522: OMPForDirective*/
bool TP7522::requestNewRangeIterations7522(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 1*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet7522 * codeletID;
        int tempEndRange = rangePerCodelet7522 * (codeletID + 1);
        if (remainderRange7522 != 0) {
            if (codeletID < (uint32_t)remainderRange7522) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange7522;
                tempEndRange += remainderRange7522;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration7522;
        tempEndRange = tempEndRange * 1 + minIteration7522;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration7522 < lastIteration7522) {
            (this->inputsTPParent->i_darts7522[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts7522[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == 0) {
            *endRange = *endRange - 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration7522;
        }
    }
    return isThereNewIteration;
}
void TP7522::_checkInCodelets7523::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->c1345_darts7522[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->c1345_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->c34_darts7522[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->c34_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->k_darts7522[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->k_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->r43_darts7522[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->r43_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp1_darts7522[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp1_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp2_darts7522[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp2_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp3_darts7522[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp3_darts8[this->getID()]);

    /*printing node 7523: ForStmt*/
    /*var: c1345*/
    /*var: c34*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: r43*/
    /*var: tmp1*/
    /*var: tmp2*/
    /*var: tmp3*/
    double** c1345 = &(this->inputsTPParent->c1345_darts7522[this->getLocalID()]);
    (void)c1345 /*OMP_SHARED_PRIVATE*/;
    double** c34 = &(this->inputsTPParent->c34_darts7522[this->getLocalID()]);
    (void)c34 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts7522[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts7522[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** k = &(this->inputsTPParent->k_darts7522[this->getLocalID()]);
    (void)k /*OMP_SHARED_PRIVATE*/;
    double** r43 = &(this->inputsTPParent->r43_darts7522[this->getLocalID()]);
    (void)r43 /*OMP_SHARED_PRIVATE*/;
    double** tmp1 = &(this->inputsTPParent->tmp1_darts7522[this->getLocalID()]);
    (void)tmp1 /*OMP_SHARED_PRIVATE*/;
    double** tmp2 = &(this->inputsTPParent->tmp2_darts7522[this->getLocalID()]);
    (void)tmp2 /*OMP_SHARED_PRIVATE*/;
    double** tmp3 = &(this->inputsTPParent->tmp3_darts7522[this->getLocalID()]);
    (void)tmp3 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations7522(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Find and signal the next codelet*/
        return;
    }
    for (int i_darts_counter_temp7522 = (*i); i_darts_counter_temp7522 >= endRange
         && i_darts_counter_temp7522 >= this->inputsTPParent->lastIteration7522;
         i_darts_counter_temp7522--) {
        {
            {
                /*Loop's init*/
                (*j) = jend;
                int j_darts_counter_temp7522 = (*j);
                for (; j_darts_counter_temp7522 >= jst; j_darts_counter_temp7522--) {
                    (*(*tmp1))
                        = 1. / u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][0]
                        = 1. + dt * 2. * (tx1 * dx1 + ty1 * dy1 + tz1 * dz1);
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][1] = 0.;
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][2] = 0.;
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][3] = 0.;
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][4] = 0.;
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][0] = dt * 2.
                        * (tx1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                       [(*(*k))][1])
                            + ty1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                       [(*(*k))][1])
                            + tz1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                       [(*(*k))][1]));
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][1] = 1.
                        + dt * 2.
                            * (tx1 * (*(*r43)) * (*(*c34)) * (*(*tmp1))
                                + ty1 * (*(*c34)) * (*(*tmp1)) + tz1 * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx2 + ty1 * dy2 + tz1 * dz2);
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][2] = 0.;
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][3] = 0.;
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][4] = 0.;
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][0] = dt * 2.
                        * (tx1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                       [(*(*k))][2])
                            + ty1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                       [(*(*k))][2])
                            + tz1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                       [(*(*k))][2]));
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][1] = 0.;
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][2] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c34)) * (*(*tmp1))
                                + ty1 * (*(*r43)) * (*(*c34)) * (*(*tmp1))
                                + tz1 * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx3 + ty1 * dy3 + tz1 * dz3);
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][3] = 0.;
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][4] = 0.;
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][0] = dt * 2.
                        * (tx1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                       [(*(*k))][3])
                            + ty1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                       [(*(*k))][3])
                            + tz1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                       [(*(*k))][3]));
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][1] = 0.;
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][2] = 0.;
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][3] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c34)) * (*(*tmp1)) + ty1 * (*(*c34)) * (*(*tmp1))
                                + tz1 * (*(*r43)) * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx4 + ty1 * dy4 + tz1 * dz4);
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][4] = 0.;
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][0] = dt * 2.
                        * (tx1
                                * (-((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp7522)]
                                                [j_darts_counter_temp7522][(*(*k))][1])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp7522)]
                                                [j_darts_counter_temp7522][(*(*k))][2])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp7522)]
                                                [j_darts_counter_temp7522][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                           [(*(*k))][4])
                            + ty1
                                * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp7522)]
                                                [j_darts_counter_temp7522][(*(*k))][1])))
                                    - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp7522)]
                                                [j_darts_counter_temp7522][(*(*k))][2])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp7522)]
                                                [j_darts_counter_temp7522][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                           [(*(*k))][4])
                            + tz1
                                * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp7522)]
                                                [j_darts_counter_temp7522][(*(*k))][1])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp7522)]
                                                [j_darts_counter_temp7522][(*(*k))][2])))
                                    - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp7522)]
                                                [j_darts_counter_temp7522][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                           [(*(*k))][4]));
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][1] = dt * 2.
                        * (tx1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k))]
                                   [1]
                            + ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k))]
                                   [1]
                            + tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k))]
                                   [1]);
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][2] = dt * 2.
                        * (tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k))]
                                   [2]
                            + ty1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k))]
                                   [2]
                            + tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k))]
                                   [2]);
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][3] = dt * 2.
                        * (tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k))]
                                   [3]
                            + ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k))]
                                   [3]
                            + tz1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k))]
                                   [3]);
                    d[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][4] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c1345)) * (*(*tmp1)) + ty1 * (*(*c1345)) * (*(*tmp1))
                                + tz1 * (*(*c1345)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx5 + ty1 * dy5 + tz1 * dz5);
                    (*(*tmp1)) = 1.
                        / u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][0] = -dt * tx1 * dx1;
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][1] = dt * tx2;
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][2] = 0.;
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][3] = 0.;
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][4] = 0.;
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][0] = dt * tx2
                            * (-(u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                  [(*(*k))][1]
                                   * (*(*tmp1)))
                                    * (u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                        [(*(*k))][1]
                                        * (*(*tmp1)))
                                + 0.40000000000000002 * 0.5
                                    * (u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                        [(*(*k))][1]
                                            * u[(i_darts_counter_temp7522) + 1]
                                               [j_darts_counter_temp7522][(*(*k))][1]
                                        + u[(i_darts_counter_temp7522) + 1]
                                           [j_darts_counter_temp7522][(*(*k))][2]
                                            * u[(i_darts_counter_temp7522) + 1]
                                               [j_darts_counter_temp7522][(*(*k))][2]
                                        + u[(i_darts_counter_temp7522) + 1]
                                           [j_darts_counter_temp7522][(*(*k))][3]
                                            * u[(i_darts_counter_temp7522) + 1]
                                               [j_darts_counter_temp7522][(*(*k))][3])
                                    * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                   [(*(*k))][1]);
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][1] = dt * tx2
                            * ((2. - 0.40000000000000002)
                                * (u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * tx1 * dx2;
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][2] = dt * tx2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522][(*(*k))]
                                [2]
                                * (*(*tmp1))));
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][3] = dt * tx2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522][(*(*k))]
                                [3]
                                * (*(*tmp1))));
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][4]
                        = dt * tx2 * 0.40000000000000002;
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][0] = dt * tx2
                            * (-(u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                  [(*(*k))][1]
                                   * u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                      [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                   [(*(*k))][2]);
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][1] = dt * tx2
                        * (u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522][(*(*k))][2]
                            * (*(*tmp1)));
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][2] = dt * tx2
                            * (u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522][(*(*k))]
                                [1]
                                * (*(*tmp1)))
                        - dt * tx1 * ((*(*c34)) * (*(*tmp1))) - dt * tx1 * dx3;
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][3] = 0.;
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][4] = 0.;
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][0] = dt * tx2
                            * (-(u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                  [(*(*k))][1]
                                   * u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                      [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                   [(*(*k))][3]);
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][1] = dt * tx2
                        * (u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522][(*(*k))][3]
                            * (*(*tmp1)));
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][2] = 0.;
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][3] = dt * tx2
                            * (u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522][(*(*k))]
                                [1]
                                * (*(*tmp1)))
                        - dt * tx1 * ((*(*c34)) * (*(*tmp1))) - dt * tx1 * dx4;
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][4] = 0.;
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][0] = dt * tx2
                            * ((0.40000000000000002
                                       * (u[(i_darts_counter_temp7522) + 1]
                                           [j_darts_counter_temp7522][(*(*k))][1]
                                               * u[(i_darts_counter_temp7522) + 1]
                                                  [j_darts_counter_temp7522][(*(*k))][1]
                                           + u[(i_darts_counter_temp7522) + 1]
                                              [j_darts_counter_temp7522][(*(*k))][2]
                                               * u[(i_darts_counter_temp7522) + 1]
                                                  [j_darts_counter_temp7522][(*(*k))][2]
                                           + u[(i_darts_counter_temp7522) + 1]
                                              [j_darts_counter_temp7522][(*(*k))][3]
                                               * u[(i_darts_counter_temp7522) + 1]
                                                  [j_darts_counter_temp7522][(*(*k))][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(i_darts_counter_temp7522) + 1]
                                           [j_darts_counter_temp7522][(*(*k))][4]
                                           * (*(*tmp1))))
                                * (u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1
                            * (-((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                          [(*(*k))][1])
                                        * (u[(i_darts_counter_temp7522) + 1]
                                            [j_darts_counter_temp7522][(*(*k))][1])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                          [(*(*k))][2])
                                        * (u[(i_darts_counter_temp7522) + 1]
                                            [j_darts_counter_temp7522][(*(*k))][2])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                          [(*(*k))][3])
                                        * (u[(i_darts_counter_temp7522) + 1]
                                            [j_darts_counter_temp7522][(*(*k))][3])))
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                       [(*(*k))][4]);
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][1] = dt * tx2
                            * (1.3999999999999999
                                    * (u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                        [(*(*k))][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((3.
                                               * u[(i_darts_counter_temp7522) + 1]
                                                  [j_darts_counter_temp7522][(*(*k))][1]
                                               * u[(i_darts_counter_temp7522) + 1]
                                                  [j_darts_counter_temp7522][(*(*k))][1]
                                           + u[(i_darts_counter_temp7522) + 1]
                                              [j_darts_counter_temp7522][(*(*k))][2]
                                               * u[(i_darts_counter_temp7522) + 1]
                                                  [j_darts_counter_temp7522][(*(*k))][2]
                                           + u[(i_darts_counter_temp7522) + 1]
                                              [j_darts_counter_temp7522][(*(*k))][3]
                                               * u[(i_darts_counter_temp7522) + 1]
                                                  [j_darts_counter_temp7522][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * tx1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522][(*(*k))]
                               [1];
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][2] = dt * tx2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                    [(*(*k))][2]
                                    * u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                       [(*(*k))][1])
                                * (*(*tmp2)))
                        - dt * tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522][(*(*k))]
                               [2];
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][3] = dt * tx2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                    [(*(*k))][3]
                                    * u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                       [(*(*k))][1])
                                * (*(*tmp2)))
                        - dt * tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522][(*(*k))]
                               [3];
                    a[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][4] = dt * tx2
                            * (1.3999999999999999
                                * (u[(i_darts_counter_temp7522) + 1][j_darts_counter_temp7522]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1 * (*(*c1345)) * (*(*tmp1)) - dt * tx1 * dx5;
                    (*(*tmp1)) = 1.
                        / u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][0] = -dt * ty1 * dy1;
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][1] = 0.;
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][2] = dt * ty2;
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][3] = 0.;
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][4] = 0.;
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][0] = dt * ty2
                            * (-(u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                  [(*(*k))][1]
                                   * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                      [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * ty1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                   [(*(*k))][1]);
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][1] = dt * ty2
                            * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1][(*(*k))]
                                [2]
                                * (*(*tmp1)))
                        - dt * ty1 * ((*(*c34)) * (*(*tmp1))) - dt * ty1 * dy2;
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][2] = dt * ty2
                        * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1][(*(*k))][1]
                            * (*(*tmp1)));
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][3] = 0.;
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][4] = 0.;
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][0] = dt * ty2
                            * (-(u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                  [(*(*k))][2]
                                   * (*(*tmp1)))
                                    * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                        [(*(*k))][2]
                                        * (*(*tmp1)))
                                + 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                         [(*(*k))][1]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522 + 1][(*(*k))][1]
                                           + u[(i_darts_counter_temp7522)]
                                              [j_darts_counter_temp7522 + 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522 + 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp7522)]
                                              [j_darts_counter_temp7522 + 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522 + 1][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * ty1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                   [(*(*k))][2]);
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][1] = dt * ty2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1][(*(*k))]
                                [1]
                                * (*(*tmp1))));
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][2] = dt * ty2
                            * ((2. - 0.40000000000000002)
                                * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * ty1 * dy3;
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][3] = dt * ty2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1][(*(*k))]
                                [3]
                                * (*(*tmp1))));
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][4]
                        = dt * ty2 * 0.40000000000000002;
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][0] = dt * ty2
                            * (-(u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                  [(*(*k))][2]
                                   * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                      [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * ty1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                   [(*(*k))][3]);
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][1] = 0.;
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][2] = dt * ty2
                        * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1][(*(*k))][3]
                            * (*(*tmp1)));
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][3] = dt * ty2
                            * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1][(*(*k))]
                                [2]
                                * (*(*tmp1)))
                        - dt * ty1 * ((*(*c34)) * (*(*tmp1))) - dt * ty1 * dy4;
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][4] = 0.;
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][0] = dt * ty2
                            * ((0.40000000000000002
                                       * (u[(i_darts_counter_temp7522)]
                                           [j_darts_counter_temp7522 + 1][(*(*k))][1]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522 + 1][(*(*k))][1]
                                           + u[(i_darts_counter_temp7522)]
                                              [j_darts_counter_temp7522 + 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522 + 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp7522)]
                                              [j_darts_counter_temp7522 + 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522 + 1][(*(*k))][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(i_darts_counter_temp7522)]
                                           [j_darts_counter_temp7522 + 1][(*(*k))][4]
                                           * (*(*tmp1))))
                                * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1
                            * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                          [(*(*k))][1])
                                        * (u[(i_darts_counter_temp7522)]
                                            [j_darts_counter_temp7522 + 1][(*(*k))][1])))
                                - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                          [(*(*k))][2])
                                        * (u[(i_darts_counter_temp7522)]
                                            [j_darts_counter_temp7522 + 1][(*(*k))][2])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                          [(*(*k))][3])
                                        * (u[(i_darts_counter_temp7522)]
                                            [j_darts_counter_temp7522 + 1][(*(*k))][3])))
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                       [(*(*k))][4]);
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][1] = dt * ty2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                    [(*(*k))][1]
                                    * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                       [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1][(*(*k))]
                               [1];
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][2] = dt * ty2
                            * (1.3999999999999999
                                    * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                        [(*(*k))][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                         [(*(*k))][1]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522 + 1][(*(*k))][1]
                                           + 3.
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522 + 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522 + 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp7522)]
                                              [j_darts_counter_temp7522 + 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522 + 1][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * ty1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1][(*(*k))]
                               [2];
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][3] = dt * ty2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                    [(*(*k))][2]
                                    * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                       [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1][(*(*k))]
                               [3];
                    b[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][4] = dt * ty2
                            * (1.3999999999999999
                                * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522 + 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1 * (*(*c1345)) * (*(*tmp1)) - dt * ty1 * dy5;
                    (*(*tmp1)) = 1.
                        / u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k)) + 1][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][0] = -dt * tz1 * dz1;
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][1] = 0.;
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][2] = 0.;
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][3] = dt * tz2;
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][0][4] = 0.;
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][0] = dt * tz2
                            * (-(u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                  [(*(*k)) + 1][1]
                                   * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                      [(*(*k)) + 1][3])
                                * (*(*tmp2)))
                        - dt * tz1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                   [(*(*k)) + 1][1]);
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][1] = dt * tz2
                            * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k)) + 1]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * (*(*c34)) * (*(*tmp1)) - dt * tz1 * dz2;
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][2] = 0.;
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][3] = dt * tz2
                        * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k)) + 1][1]
                            * (*(*tmp1)));
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][1][4] = 0.;
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][0] = dt * tz2
                            * (-(u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                  [(*(*k)) + 1][2]
                                   * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                      [(*(*k)) + 1][3])
                                * (*(*tmp2)))
                        - dt * tz1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                   [(*(*k)) + 1][2]);
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][1] = 0.;
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][2] = dt * tz2
                            * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k)) + 1]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * ((*(*c34)) * (*(*tmp1))) - dt * tz1 * dz3;
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][3] = dt * tz2
                        * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k)) + 1][2]
                            * (*(*tmp1)));
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][2][4] = 0.;
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][0] = dt * tz2
                            * (-(u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                  [(*(*k)) + 1][3]
                                   * (*(*tmp1)))
                                    * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                        [(*(*k)) + 1][3]
                                        * (*(*tmp1)))
                                + 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                         [(*(*k)) + 1][1]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522][(*(*k)) + 1][1]
                                           + u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                              [(*(*k)) + 1][2]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522][(*(*k)) + 1][2]
                                           + u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                              [(*(*k)) + 1][3]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522][(*(*k)) + 1][3])
                                        * (*(*tmp2))))
                        - dt * tz1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                   [(*(*k)) + 1][3]);
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][1] = dt * tz2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k)) + 1]
                                [1]
                                * (*(*tmp1))));
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][2] = dt * tz2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k)) + 1]
                                [2]
                                * (*(*tmp1))));
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][3] = dt * tz2
                            * (2. - 0.40000000000000002)
                            * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k)) + 1]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * tz1 * dz4;
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][3][4]
                        = dt * tz2 * 0.40000000000000002;
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][0] = dt * tz2
                            * ((0.40000000000000002
                                       * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                           [(*(*k)) + 1][1]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522][(*(*k)) + 1][1]
                                           + u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                              [(*(*k)) + 1][2]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522][(*(*k)) + 1][2]
                                           + u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                              [(*(*k)) + 1][3]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522][(*(*k)) + 1][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                           [(*(*k)) + 1][4]
                                           * (*(*tmp1))))
                                * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                    [(*(*k)) + 1][3]
                                    * (*(*tmp1))))
                        - dt * tz1
                            * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                          [(*(*k)) + 1][1])
                                        * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                            [(*(*k)) + 1][1])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                          [(*(*k)) + 1][2])
                                        * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                            [(*(*k)) + 1][2])))
                                - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                          [(*(*k)) + 1][3])
                                        * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                            [(*(*k)) + 1][3])))
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                       [(*(*k)) + 1][4]);
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][1] = dt * tz2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                    [(*(*k)) + 1][1]
                                    * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                       [(*(*k)) + 1][3])
                                * (*(*tmp2)))
                        - dt * tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k)) + 1]
                               [1];
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][2] = dt * tz2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                    [(*(*k)) + 1][2]
                                    * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                       [(*(*k)) + 1][3])
                                * (*(*tmp2)))
                        - dt * tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k)) + 1]
                               [2];
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][3] = dt * tz2
                            * (1.3999999999999999
                                    * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                        [(*(*k)) + 1][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                         [(*(*k)) + 1][1]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522][(*(*k)) + 1][1]
                                           + u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                              [(*(*k)) + 1][2]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522][(*(*k)) + 1][2]
                                           + 3.
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522][(*(*k)) + 1][3]
                                               * u[(i_darts_counter_temp7522)]
                                                  [j_darts_counter_temp7522][(*(*k)) + 1][3])
                                        * (*(*tmp2))))
                        - dt * tz1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7522)][j_darts_counter_temp7522][(*(*k)) + 1]
                               [3];
                    c[(i_darts_counter_temp7522)][j_darts_counter_temp7522][4][4] = dt * tz2
                            * (1.3999999999999999
                                * (u[(i_darts_counter_temp7522)][j_darts_counter_temp7522]
                                    [(*(*k)) + 1][3]
                                    * (*(*tmp1))))
                        - dt * tz1 * (*(*c1345)) * (*(*tmp1)) - dt * tz1 * dz5;
                }
                (*j) = j_darts_counter_temp7522;
            }
        }
    }
    /*If this omp for has no barrier,
    check if all the codelets
    replicated from the same
    global ID has finished and
    signal the next codelet.
    Otherwise, return.*/
    uint32_t completedMultCodelet = __sync_fetch_and_add(
        &(myTP->signalNextReady[this->getLocalID() % myTP->baseNumThreads]), 1);
    if (completedMultCodelet < (uint32_t)(DARTS_CODELETS_MULT - 1))
        return;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Find and signal the next codelet*/
}
TP7522::TP7522(int in_numThreads, int in_mainCodeletID, TP8* in_TPParent, int in_initIteration,
    int in_lastIteration, TP7522** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , c1345_darts7522(new double*[this->numThreads])
    , c34_darts7522(new double*[this->numThreads])
    , i_darts7522(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts7522(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts7522(new int*[this->numThreads])
    , r43_darts7522(new double*[this->numThreads])
    , tmp1_darts7522(new double*[this->numThreads])
    , tmp2_darts7522(new double*[this->numThreads])
    , tmp3_darts7522(new double*[this->numThreads])
    , initIteration7522(in_initIteration)
    , lastIteration7522(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , signalNextReady(new int[baseNumThreads])
    , checkInCodelets7523(new _checkInCodelets7523[this->numThreads])
{
    /*Initialize the loop parameters*/
    range7522 = abs(lastIteration7522 - initIteration7522) / 1;
    rangePerCodelet7522 = range7522 / numThreads;
    minIteration7522 = min<int>(lastIteration7522, initIteration7522);
    remainderRange7522 = range7522 % numThreads;
    /*Initialize inputs and vars.*/
    this->c1345_darts7522
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->c34_darts7522
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->k_darts7522 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->r43_darts7522
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp1_darts7522
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp2_darts7522
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp3_darts7522
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    _checkInCodelets7523* checkInCodelets7523Ptr = (this->checkInCodelets7523);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets7523);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
        this->signalNextReady[codeletCounter] = 0;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets7523Ptr) = _checkInCodelets7523(2, 1, this, codeletCounter);
#else
        (*checkInCodelets7523Ptr) = _checkInCodelets7523(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets7523Ptr).decDep();
        checkInCodelets7523Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP7522::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets7523[localID].setID(codeletID);
    this->checkInCodelets7523[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets7523[localID + this->baseNumThreads * i]
            = _checkInCodelets7523(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets7523[localID + this->baseNumThreads * i]
            = _checkInCodelets7523(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets7523[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets7523[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP7522::~TP7522()
{
    delete[] c1345_darts7522;
    delete[] c34_darts7522;
    delete[] k_darts7522;
    delete[] r43_darts7522;
    delete[] tmp1_darts7522;
    delete[] tmp2_darts7522;
    delete[] tmp3_darts7522;
    delete[] checkInCodelets7523;
}
/*TP10896: OMPParallelDirective*/
void TP10896::_barrierCodelets10896::fire(void)
{
    TP10896* myTP = static_cast<TP10896*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP10896::_checkInCodelets10911::fire(void)
{
    /*region 10911 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP10911;
    if (idx < myTP->TPsToUse10911) {
        if (!__sync_val_compare_and_swap(&(myTP->TP10911_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(nx - 1 - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse10911;
            int minIteration = min<int>(nx - 1, 0);
            int remainderRange = range % myTP->TPsToUse10911;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (0 < nx - 1) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse10911 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse10911 - 1) {
                lastIteration = nx - 1;
            }
#if USEINVOKE == 1
            invoke<TP10911>(myTP, myTP->codeletsPerTP10911 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP10911Ptr[idx]));
#else
            place<TP10911>(idx, myTP, myTP->codeletsPerTP10911 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP10911Ptr[idx]));
#endif
        } else {
            if (myTP->TP10911Ptr[idx] != nullptr) {
                myTP->TP10911Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10896::_barrierCodelets10911::fire(void)
{
    TP10896* myTP = static_cast<TP10896*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets10968[codeletsCounter].decDep();
        }
    }
}
void TP10896::_checkInCodelets10968::fire(void)
{

    /*printing node 10968: BinaryOperator*/
    (this->inputsTPParent->L1_darts10896[this->getID()]) = 0;

    /*printing node 10969: BinaryOperator*/
    (this->inputsTPParent->L2_darts10896[this->getID()]) = nx - 1;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 10968 nextRegion: 10971 */
    myTP->controlTPParent->checkInCodelets10971[this->getID()].decDep();
}
void TP10896::_checkInCodelets10971::fire(void)
{
    /*region 10971 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP10971;
    if (idx < myTP->TPsToUse10971) {
        if (!__sync_val_compare_and_swap(&(myTP->TP10971_alreadyLaunched[idx]), 0, 1)) {
            int range = abs((this->inputsTPParent->L2_darts10896[this->getID()])
                            - (this->inputsTPParent->L1_darts10896[this->getID()]))
                / 1;
            int rangePerCodelet = range / myTP->TPsToUse10971;
            int minIteration = min<int>((this->inputsTPParent->L2_darts10896[this->getID()]),
                (this->inputsTPParent->L1_darts10896[this->getID()]));
            int remainderRange = range % myTP->TPsToUse10971;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if ((this->inputsTPParent->L1_darts10896[this->getID()])
                < (this->inputsTPParent->L2_darts10896[this->getID()])) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse10971 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse10971 - 1) {
                lastIteration = (this->inputsTPParent->L2_darts10896[this->getID()]);
            }
#if USEINVOKE == 1
            invoke<TP10971>(myTP, myTP->codeletsPerTP10971 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP10971Ptr[idx]));
#else
            place<TP10971>(idx, myTP, myTP->codeletsPerTP10971 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP10971Ptr[idx]));
#endif
        } else {
            if (myTP->TP10971Ptr[idx] != nullptr) {
                myTP->TP10971Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10896::_barrierCodelets10971::fire(void)
{
    TP10896* myTP = static_cast<TP10896*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets11120[codeletsCounter].decDep();
        }
    }
}
void TP10896::_checkInCodelets11120::fire(void)
{
    /*region 11120 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP11120;
    if (idx < myTP->TPsToUse11120) {
        if (!__sync_val_compare_and_swap(&(myTP->TP11120_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(jend - jst) / 1;
            int rangePerCodelet = range / myTP->TPsToUse11120;
            int minIteration = min<int>(jend, jst);
            int remainderRange = range % myTP->TPsToUse11120;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (jst < jend) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse11120 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse11120 - 1) {
                lastIteration = jend;
            }
#if USEINVOKE == 1
            invoke<TP11120>(myTP, myTP->codeletsPerTP11120 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11120Ptr[idx]));
#else
            place<TP11120>(idx, myTP, myTP->codeletsPerTP11120 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11120Ptr[idx]));
#endif
        } else {
            if (myTP->TP11120Ptr[idx] != nullptr) {
                myTP->TP11120Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10896::_barrierCodelets11120::fire(void)
{
    TP10896* myTP = static_cast<TP10896*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets11757[codeletsCounter].decDep();
        }
    }
}
void TP10896::_checkInCodelets11757::fire(void)
{

    /*printing node 11757: BinaryOperator*/
    (this->inputsTPParent->L1_darts10896[this->getID()]) = 0;

    /*printing node 11758: BinaryOperator*/
    (this->inputsTPParent->L2_darts10896[this->getID()]) = ny - 1;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 11757 nextRegion: 11760 */
    myTP->controlTPParent->checkInCodelets11760[this->getID()].decDep();
}
void TP10896::_checkInCodelets11760::fire(void)
{
    /*region 11760 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP11760;
    if (idx < myTP->TPsToUse11760) {
        if (!__sync_val_compare_and_swap(&(myTP->TP11760_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse11760;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse11760;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (ist < iend) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse11760 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse11760 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP11760>(myTP, myTP->codeletsPerTP11760 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11760Ptr[idx]));
#else
            place<TP11760>(idx, myTP, myTP->codeletsPerTP11760 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11760Ptr[idx]));
#endif
        } else {
            if (myTP->TP11760Ptr[idx] != nullptr) {
                myTP->TP11760Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10896::_barrierCodelets11760::fire(void)
{
    TP10896* myTP = static_cast<TP10896*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets11909[codeletsCounter].decDep();
        }
    }
}
void TP10896::_checkInCodelets11909::fire(void)
{
    /*region 11909 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP11909;
    if (idx < myTP->TPsToUse11909) {
        if (!__sync_val_compare_and_swap(&(myTP->TP11909_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse11909;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse11909;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (ist < iend) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse11909 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse11909 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP11909>(myTP, myTP->codeletsPerTP11909 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11909Ptr[idx]));
#else
            place<TP11909>(idx, myTP, myTP->codeletsPerTP11909 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11909Ptr[idx]));
#endif
        } else {
            if (myTP->TP11909Ptr[idx] != nullptr) {
                myTP->TP11909Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10896::_barrierCodelets11909::fire(void)
{
    TP10896* myTP = static_cast<TP10896*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets12546[codeletsCounter].decDep();
        }
    }
}
void TP10896::_checkInCodelets12546::fire(void)
{
    /*region 12546 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP12546;
    if (idx < myTP->TPsToUse12546) {
        if (!__sync_val_compare_and_swap(&(myTP->TP12546_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse12546;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse12546;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (ist < iend) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse12546 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse12546 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP12546>(myTP, myTP->codeletsPerTP12546 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP12546Ptr[idx]));
#else
            place<TP12546>(idx, myTP, myTP->codeletsPerTP12546 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP12546Ptr[idx]));
#endif
        } else {
            if (myTP->TP12546Ptr[idx] != nullptr) {
                myTP->TP12546Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10896::_barrierCodelets12546::fire(void)
{
    TP10896* myTP = static_cast<TP10896*>(myTP_);
    myTP->TPParent->barrierCodelets10896[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets10896[0]));
}
TP10896::TP10896(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts10896(new int[this->numThreads]) /*VARIABLE*/
    , L2_darts10896(new int[this->numThreads]) /*VARIABLE*/
    , i_darts10896(new int[this->numThreads]) /*VARIABLE*/
    , iend1_darts10896(new int[this->numThreads]) /*VARIABLE*/
    , ist1_darts10896(new int[this->numThreads]) /*VARIABLE*/
    , j_darts10896(new int[this->numThreads]) /*VARIABLE*/
    , jend1_darts10896(new int[this->numThreads]) /*VARIABLE*/
    , jst1_darts10896(new int[this->numThreads]) /*VARIABLE*/
    , k_darts10896(new int[this->numThreads]) /*VARIABLE*/
    , m_darts10896(new int[this->numThreads]) /*VARIABLE*/
    , q_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , tmp_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u21_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u21i_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u21im1_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u21j_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u21jm1_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u21k_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u21km1_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u31_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u31i_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u31im1_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u31j_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u31jm1_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u31k_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u31km1_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u41_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u41i_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u41im1_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u41j_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u41jm1_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u41k_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u41km1_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u51i_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u51im1_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u51j_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u51jm1_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u51k_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , u51km1_darts10896(new double[this->numThreads]) /*VARIABLE*/
    , TP10911Ptr(new TP10911*[NUMTPS10911])
    , TP10911_alreadyLaunched(new size_t[NUMTPS10911])
    , numTPsSet10911(0)
    , numTPsReady10911(0)
    , TPsToUse10911(NUMTPS10911)
    , codeletsPerTP10911(this->numThreads / NUMTPS10911)
    , totalCodelets10911(this->TPsToUse10911 * this->codeletsPerTP10911)
    , TP10971Ptr(new TP10971*[NUMTPS10971])
    , TP10971_alreadyLaunched(new size_t[NUMTPS10971])
    , numTPsSet10971(0)
    , numTPsReady10971(0)
    , TPsToUse10971(NUMTPS10971)
    , codeletsPerTP10971(this->numThreads / NUMTPS10971)
    , totalCodelets10971(this->TPsToUse10971 * this->codeletsPerTP10971)
    , TP11120Ptr(new TP11120*[NUMTPS11120])
    , TP11120_alreadyLaunched(new size_t[NUMTPS11120])
    , numTPsSet11120(0)
    , numTPsReady11120(0)
    , TPsToUse11120(NUMTPS11120)
    , codeletsPerTP11120(this->numThreads / NUMTPS11120)
    , totalCodelets11120(this->TPsToUse11120 * this->codeletsPerTP11120)
    , TP11760Ptr(new TP11760*[NUMTPS11760])
    , TP11760_alreadyLaunched(new size_t[NUMTPS11760])
    , numTPsSet11760(0)
    , numTPsReady11760(0)
    , TPsToUse11760(NUMTPS11760)
    , codeletsPerTP11760(this->numThreads / NUMTPS11760)
    , totalCodelets11760(this->TPsToUse11760 * this->codeletsPerTP11760)
    , TP11909Ptr(new TP11909*[NUMTPS11909])
    , TP11909_alreadyLaunched(new size_t[NUMTPS11909])
    , numTPsSet11909(0)
    , numTPsReady11909(0)
    , TPsToUse11909(NUMTPS11909)
    , codeletsPerTP11909(this->numThreads / NUMTPS11909)
    , totalCodelets11909(this->TPsToUse11909 * this->codeletsPerTP11909)
    , TP12546Ptr(new TP12546*[NUMTPS12546])
    , TP12546_alreadyLaunched(new size_t[NUMTPS12546])
    , numTPsSet12546(0)
    , numTPsReady12546(0)
    , TPsToUse12546(NUMTPS12546)
    , codeletsPerTP12546(this->numThreads / NUMTPS12546)
    , totalCodelets12546(this->TPsToUse12546 * this->codeletsPerTP12546)
    , barrierCodelets10896(new _barrierCodelets10896[1])
    , checkInCodelets10911(new _checkInCodelets10911[this->numThreads])
    , barrierCodelets10911(new _barrierCodelets10911[1])
    , checkInCodelets10968(new _checkInCodelets10968[this->numThreads])
    , checkInCodelets10971(new _checkInCodelets10971[this->numThreads])
    , barrierCodelets10971(new _barrierCodelets10971[1])
    , checkInCodelets11120(new _checkInCodelets11120[this->numThreads])
    , barrierCodelets11120(new _barrierCodelets11120[1])
    , checkInCodelets11757(new _checkInCodelets11757[this->numThreads])
    , checkInCodelets11760(new _checkInCodelets11760[this->numThreads])
    , barrierCodelets11760(new _barrierCodelets11760[1])
    , checkInCodelets11909(new _checkInCodelets11909[this->numThreads])
    , barrierCodelets11909(new _barrierCodelets11909[1])
    , checkInCodelets12546(new _checkInCodelets12546[this->numThreads])
    , barrierCodelets12546(new _barrierCodelets12546[1])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets10896[0] = _barrierCodelets10896(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets12546[0] = _barrierCodelets12546(NUMTPS12546, NUMTPS12546, this, 0);
    barrierCodelets11909[0] = _barrierCodelets11909(NUMTPS11909, NUMTPS11909, this, 0);
    barrierCodelets11760[0] = _barrierCodelets11760(NUMTPS11760, NUMTPS11760, this, 0);
    barrierCodelets11120[0] = _barrierCodelets11120(NUMTPS11120, NUMTPS11120, this, 0);
    barrierCodelets10971[0] = _barrierCodelets10971(NUMTPS10971, NUMTPS10971, this, 0);
    barrierCodelets10911[0] = _barrierCodelets10911(NUMTPS10911, NUMTPS10911, this, 0);
    _checkInCodelets12546* checkInCodelets12546Ptr = (this->checkInCodelets12546);
    for (int i = 0; i < NUMTPS12546; i++) {
        TP12546Ptr[i] = nullptr;
        TP12546_alreadyLaunched[i] = 0;
    }
    _checkInCodelets11909* checkInCodelets11909Ptr = (this->checkInCodelets11909);
    for (int i = 0; i < NUMTPS11909; i++) {
        TP11909Ptr[i] = nullptr;
        TP11909_alreadyLaunched[i] = 0;
    }
    _checkInCodelets11760* checkInCodelets11760Ptr = (this->checkInCodelets11760);
    for (int i = 0; i < NUMTPS11760; i++) {
        TP11760Ptr[i] = nullptr;
        TP11760_alreadyLaunched[i] = 0;
    }
    _checkInCodelets11757* checkInCodelets11757Ptr = (this->checkInCodelets11757);
    _checkInCodelets11120* checkInCodelets11120Ptr = (this->checkInCodelets11120);
    for (int i = 0; i < NUMTPS11120; i++) {
        TP11120Ptr[i] = nullptr;
        TP11120_alreadyLaunched[i] = 0;
    }
    _checkInCodelets10971* checkInCodelets10971Ptr = (this->checkInCodelets10971);
    for (int i = 0; i < NUMTPS10971; i++) {
        TP10971Ptr[i] = nullptr;
        TP10971_alreadyLaunched[i] = 0;
    }
    _checkInCodelets10968* checkInCodelets10968Ptr = (this->checkInCodelets10968);
    _checkInCodelets10911* checkInCodelets10911Ptr = (this->checkInCodelets10911);
    for (int i = 0; i < NUMTPS10911; i++) {
        TP10911Ptr[i] = nullptr;
        TP10911_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets12546Ptr) = _checkInCodelets12546(1, 1, this, codeletCounter);
        checkInCodelets12546Ptr++;
        (*checkInCodelets11909Ptr) = _checkInCodelets11909(1, 1, this, codeletCounter);
        checkInCodelets11909Ptr++;
        (*checkInCodelets11760Ptr) = _checkInCodelets11760(1, 1, this, codeletCounter);
        checkInCodelets11760Ptr++;
        (*checkInCodelets11757Ptr) = _checkInCodelets11757(1, 1, this, codeletCounter);
        checkInCodelets11757Ptr++;
        (*checkInCodelets11120Ptr) = _checkInCodelets11120(1, 1, this, codeletCounter);
        checkInCodelets11120Ptr++;
        (*checkInCodelets10971Ptr) = _checkInCodelets10971(1, 1, this, codeletCounter);
        checkInCodelets10971Ptr++;
        (*checkInCodelets10968Ptr) = _checkInCodelets10968(1, 1, this, codeletCounter);
        checkInCodelets10968Ptr++;
        (*checkInCodelets10911Ptr) = _checkInCodelets10911(1, 1, this, codeletCounter);
        (*checkInCodelets10911Ptr).decDep();
        checkInCodelets10911Ptr++;
    }
}
TP10896::~TP10896()
{
    delete[] L1_darts10896;
    delete[] L2_darts10896;
    delete[] i_darts10896;
    delete[] iend1_darts10896;
    delete[] ist1_darts10896;
    delete[] j_darts10896;
    delete[] jend1_darts10896;
    delete[] jst1_darts10896;
    delete[] k_darts10896;
    delete[] m_darts10896;
    delete[] q_darts10896;
    delete[] tmp_darts10896;
    delete[] u21_darts10896;
    delete[] u21i_darts10896;
    delete[] u21im1_darts10896;
    delete[] u21j_darts10896;
    delete[] u21jm1_darts10896;
    delete[] u21k_darts10896;
    delete[] u21km1_darts10896;
    delete[] u31_darts10896;
    delete[] u31i_darts10896;
    delete[] u31im1_darts10896;
    delete[] u31j_darts10896;
    delete[] u31jm1_darts10896;
    delete[] u31k_darts10896;
    delete[] u31km1_darts10896;
    delete[] u41_darts10896;
    delete[] u41i_darts10896;
    delete[] u41im1_darts10896;
    delete[] u41j_darts10896;
    delete[] u41jm1_darts10896;
    delete[] u41k_darts10896;
    delete[] u41km1_darts10896;
    delete[] u51i_darts10896;
    delete[] u51im1_darts10896;
    delete[] u51j_darts10896;
    delete[] u51jm1_darts10896;
    delete[] u51k_darts10896;
    delete[] u51km1_darts10896;
    delete[] barrierCodelets10896;
    delete[] barrierCodelets12546;
    delete[] checkInCodelets12546;
    delete[] barrierCodelets11909;
    delete[] checkInCodelets11909;
    delete[] barrierCodelets11760;
    delete[] checkInCodelets11760;
    delete[] checkInCodelets11757;
    delete[] barrierCodelets11120;
    delete[] checkInCodelets11120;
    delete[] barrierCodelets10971;
    delete[] checkInCodelets10971;
    delete[] checkInCodelets10968;
    delete[] barrierCodelets10911;
    delete[] checkInCodelets10911;
}
/*TP10911: OMPForDirective*/
void TP10911::_barrierCodelets10911::fire(void)
{
    TP10911* myTP = static_cast<TP10911*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets10911[0].decDep();
}
bool TP10911::requestNewRangeIterations10911(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet10911 * codeletID;
        int tempEndRange = rangePerCodelet10911 * (codeletID + 1);
        if (remainderRange10911 != 0) {
            if (codeletID < (uint32_t)remainderRange10911) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange10911;
                tempEndRange += remainderRange10911;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration10911;
        tempEndRange = tempEndRange * 1 + minIteration10911;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration10911 < lastIteration10911) {
            (this->inputsTPParent->i_darts10911[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts10911[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration10911;
        }
    }
    return isThereNewIteration;
}
void TP10911::_checkInCodelets10912::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 10912: ForStmt*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    int* i = &(this->inputsTPParent->i_darts10911[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts10911[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts10911[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts10911[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations10911(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets10911[0].decDep();
        return;
    }
    for (int i_darts_counter_temp10911 = (*i); i_darts_counter_temp10911 <= endRange
         && i_darts_counter_temp10911 <= this->inputsTPParent->lastIteration10911;
         i_darts_counter_temp10911++) {
        {
            {
                /*Loop's init*/
                (*j) = 0;
                int j_darts_counter_temp10911 = (*j);
                for (; j_darts_counter_temp10911 <= ny - 1; j_darts_counter_temp10911++) {
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp10911 = (*k);
                        for (; k_darts_counter_temp10911 <= nz - 1; k_darts_counter_temp10911++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp10911 = (*m);
                                for (; m_darts_counter_temp10911 < 5; m_darts_counter_temp10911++) {
                                    rsd[(i_darts_counter_temp10911)][j_darts_counter_temp10911]
                                       [k_darts_counter_temp10911][m_darts_counter_temp10911]
                                        = -frct[(i_darts_counter_temp10911)]
                                               [j_darts_counter_temp10911]
                                               [k_darts_counter_temp10911]
                                               [m_darts_counter_temp10911];
                                }
                                (*m) = m_darts_counter_temp10911;
                            }
                        }
                        (*k) = k_darts_counter_temp10911;
                    }
                }
                (*j) = j_darts_counter_temp10911;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets10911[0].decDep();
}
TP10911::TP10911(int in_numThreads, int in_mainCodeletID, TP10896* in_TPParent,
    int in_initIteration, int in_lastIteration, TP10911** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts10911(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts10911(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts10911(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts10911(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration10911(in_initIteration)
    , lastIteration10911(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets10911(new _barrierCodelets10911[1])
    , checkInCodelets10912(new _checkInCodelets10912[this->numThreads])
{
    /*Initialize the loop parameters*/
    range10911 = abs(lastIteration10911 - initIteration10911) / 1;
    rangePerCodelet10911 = range10911 / numThreads;
    minIteration10911 = min<int>(lastIteration10911, initIteration10911);
    remainderRange10911 = range10911 % numThreads;
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets10911[0] = _barrierCodelets10911(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets10912* checkInCodelets10912Ptr = (this->checkInCodelets10912);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets10912);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets10912Ptr) = _checkInCodelets10912(2, 1, this, codeletCounter);
#else
        (*checkInCodelets10912Ptr) = _checkInCodelets10912(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets10912Ptr).decDep();
        checkInCodelets10912Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP10911::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets10912[localID].setID(codeletID);
    this->checkInCodelets10912[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets10912[localID + this->baseNumThreads * i]
            = _checkInCodelets10912(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets10912[localID + this->baseNumThreads * i]
            = _checkInCodelets10912(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets10912[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets10912[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP10911::~TP10911()
{
    delete[] barrierCodelets10911;
    delete[] checkInCodelets10912;
}
/*TP10971: OMPForDirective*/
void TP10971::_barrierCodelets10971::fire(void)
{
    TP10971* myTP = static_cast<TP10971*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets10971[0].decDep();
}
bool TP10971::requestNewRangeIterations10971(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet10971 * codeletID;
        int tempEndRange = rangePerCodelet10971 * (codeletID + 1);
        if (remainderRange10971 != 0) {
            if (codeletID < (uint32_t)remainderRange10971) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange10971;
                tempEndRange += remainderRange10971;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration10971;
        tempEndRange = tempEndRange * 1 + minIteration10971;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration10971 < lastIteration10971) {
            (this->inputsTPParent->i_darts10971[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts10971[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration10971;
        }
    }
    return isThereNewIteration;
}
void TP10971::_checkInCodelets10972::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L1_darts10971[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L1_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L2_darts10971[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->q_darts10971[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21_darts10971[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21_darts10896[this->getID()]);

    /*printing node 10972: ForStmt*/
    /*var: L1*/
    /*var: L2*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: q*/
    /*var: u21*/
    int* i = &(this->inputsTPParent->i_darts10971[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts10971[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts10971[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts10971[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** u21 = &(this->inputsTPParent->u21_darts10971[this->getLocalID()]);
    (void)u21 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations10971(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets10971[0].decDep();
        return;
    }
    for (int i_darts_counter_temp10971 = (*i); i_darts_counter_temp10971 <= endRange
         && i_darts_counter_temp10971 <= this->inputsTPParent->lastIteration10971;
         i_darts_counter_temp10971++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp10971 = (*j);
                for (; j_darts_counter_temp10971 <= jend; j_darts_counter_temp10971++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp10971 = (*k);
                        for (; k_darts_counter_temp10971 <= nz - 2; k_darts_counter_temp10971++) {
                            flux[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                [k_darts_counter_temp10971][0]
                                = u[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                   [k_darts_counter_temp10971][1];
                            (*(*u21)) = u[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                         [k_darts_counter_temp10971][1]
                                / u[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                   [k_darts_counter_temp10971][0];
                            (*(*q)) = 0.5
                                * (u[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                    [k_darts_counter_temp10971][1]
                                        * u[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                           [k_darts_counter_temp10971][1]
                                    + u[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                       [k_darts_counter_temp10971][2]
                                        * u[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                           [k_darts_counter_temp10971][2]
                                    + u[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                       [k_darts_counter_temp10971][3]
                                        * u[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                           [k_darts_counter_temp10971][3])
                                / u[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                   [k_darts_counter_temp10971][0];
                            flux[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                [k_darts_counter_temp10971][1]
                                = u[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                   [k_darts_counter_temp10971][1]
                                    * (*(*u21))
                                + 0.40000000000000002
                                    * (u[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                        [k_darts_counter_temp10971][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                [k_darts_counter_temp10971][2]
                                = u[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                   [k_darts_counter_temp10971][2]
                                * (*(*u21));
                            flux[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                [k_darts_counter_temp10971][3]
                                = u[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                   [k_darts_counter_temp10971][3]
                                * (*(*u21));
                            flux[(i_darts_counter_temp10971)][j_darts_counter_temp10971]
                                [k_darts_counter_temp10971][4]
                                = (1.3999999999999999
                                          * u[(i_darts_counter_temp10971)]
                                             [j_darts_counter_temp10971][k_darts_counter_temp10971]
                                             [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u21));
                        }
                        (*k) = k_darts_counter_temp10971;
                    }
                }
                (*j) = j_darts_counter_temp10971;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets10971[0].decDep();
}
TP10971::TP10971(int in_numThreads, int in_mainCodeletID, TP10896* in_TPParent,
    int in_initIteration, int in_lastIteration, TP10971** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts10971(new int*[this->numThreads])
    , L2_darts10971(new int*[this->numThreads])
    , i_darts10971(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts10971(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts10971(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts10971(new double*[this->numThreads])
    , u21_darts10971(new double*[this->numThreads])
    , initIteration10971(in_initIteration)
    , lastIteration10971(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets10971(new _barrierCodelets10971[1])
    , checkInCodelets10972(new _checkInCodelets10972[this->numThreads])
{
    /*Initialize the loop parameters*/
    range10971 = abs(lastIteration10971 - initIteration10971) / 1;
    rangePerCodelet10971 = range10971 / numThreads;
    minIteration10971 = min<int>(lastIteration10971, initIteration10971);
    remainderRange10971 = range10971 % numThreads;
    /*Initialize inputs and vars.*/
    this->L1_darts10971 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->L2_darts10971 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->q_darts10971
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21_darts10971
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets10971[0] = _barrierCodelets10971(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets10972* checkInCodelets10972Ptr = (this->checkInCodelets10972);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets10972);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets10972Ptr) = _checkInCodelets10972(2, 1, this, codeletCounter);
#else
        (*checkInCodelets10972Ptr) = _checkInCodelets10972(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets10972Ptr).decDep();
        checkInCodelets10972Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP10971::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets10972[localID].setID(codeletID);
    this->checkInCodelets10972[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets10972[localID + this->baseNumThreads * i]
            = _checkInCodelets10972(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets10972[localID + this->baseNumThreads * i]
            = _checkInCodelets10972(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets10972[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets10972[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP10971::~TP10971()
{
    delete[] L1_darts10971;
    delete[] L2_darts10971;
    delete[] q_darts10971;
    delete[] u21_darts10971;
    delete[] barrierCodelets10971;
    delete[] checkInCodelets10972;
}
/*TP11120: OMPForDirective*/
void TP11120::_barrierCodelets11120::fire(void)
{
    TP11120* myTP = static_cast<TP11120*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets11120[0].decDep();
}
bool TP11120::requestNewRangeIterations11120(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet11120 * codeletID;
        int tempEndRange = rangePerCodelet11120 * (codeletID + 1);
        if (remainderRange11120 != 0) {
            if (codeletID < (uint32_t)remainderRange11120) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange11120;
                tempEndRange += remainderRange11120;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration11120;
        tempEndRange = tempEndRange * 1 + minIteration11120;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration11120 < lastIteration11120) {
            (this->inputsTPParent->j_darts11120[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->j_darts11120[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration11120;
        }
    }
    return isThereNewIteration;
}
void TP11120::_checkInCodelets11121::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L2_darts11120[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->iend1_darts11120[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iend1_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ist1_darts11120[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->ist1_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts11120[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21i_darts11120[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21i_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21im1_darts11120[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21im1_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31i_darts11120[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31i_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31im1_darts11120[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31im1_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41i_darts11120[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41i_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41im1_darts11120[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41im1_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51i_darts11120[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51i_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51im1_darts11120[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51im1_darts10896[this->getID()]);

    /*printing node 11121: ForStmt*/
    /*var: L2*/
    /*var: i*/
    /*var: iend1*/
    /*var: ist1*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    /*var: tmp*/
    /*var: u21i*/
    /*var: u21im1*/
    /*var: u31i*/
    /*var: u31im1*/
    /*var: u41i*/
    /*var: u41im1*/
    /*var: u51i*/
    /*var: u51im1*/
    int** L2 = &(this->inputsTPParent->L2_darts11120[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts11120[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iend1 = &(this->inputsTPParent->iend1_darts11120[this->getLocalID()]);
    (void)iend1 /*OMP_SHARED_PRIVATE*/;
    int** ist1 = &(this->inputsTPParent->ist1_darts11120[this->getLocalID()]);
    (void)ist1 /*OMP_SHARED_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts11120[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts11120[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts11120[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts11120[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21i = &(this->inputsTPParent->u21i_darts11120[this->getLocalID()]);
    (void)u21i /*OMP_SHARED_PRIVATE*/;
    double** u21im1 = &(this->inputsTPParent->u21im1_darts11120[this->getLocalID()]);
    (void)u21im1 /*OMP_SHARED_PRIVATE*/;
    double** u31i = &(this->inputsTPParent->u31i_darts11120[this->getLocalID()]);
    (void)u31i /*OMP_SHARED_PRIVATE*/;
    double** u31im1 = &(this->inputsTPParent->u31im1_darts11120[this->getLocalID()]);
    (void)u31im1 /*OMP_SHARED_PRIVATE*/;
    double** u41i = &(this->inputsTPParent->u41i_darts11120[this->getLocalID()]);
    (void)u41i /*OMP_SHARED_PRIVATE*/;
    double** u41im1 = &(this->inputsTPParent->u41im1_darts11120[this->getLocalID()]);
    (void)u41im1 /*OMP_SHARED_PRIVATE*/;
    double** u51i = &(this->inputsTPParent->u51i_darts11120[this->getLocalID()]);
    (void)u51i /*OMP_SHARED_PRIVATE*/;
    double** u51im1 = &(this->inputsTPParent->u51im1_darts11120[this->getLocalID()]);
    (void)u51im1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations11120(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets11120[0].decDep();
        return;
    }
    for (int j_darts_counter_temp11120 = (*j); j_darts_counter_temp11120 <= endRange
         && j_darts_counter_temp11120 <= this->inputsTPParent->lastIteration11120;
         j_darts_counter_temp11120++) {
        {
            {
                /*Loop's init*/
                (*k) = 1;
                int k_darts_counter_temp11120 = (*k);
                for (; k_darts_counter_temp11120 <= nz - 2; k_darts_counter_temp11120++) {
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp11120 = (*i);
                        for (; i_darts_counter_temp11120 <= iend; i_darts_counter_temp11120++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp11120 = (*m);
                                for (; m_darts_counter_temp11120 < 5; m_darts_counter_temp11120++) {
                                    rsd[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                                       [k_darts_counter_temp11120][m_darts_counter_temp11120]
                                        = rsd[i_darts_counter_temp11120]
                                             [(j_darts_counter_temp11120)]
                                             [k_darts_counter_temp11120][m_darts_counter_temp11120]
                                        - tx2
                                            * (flux[i_darts_counter_temp11120 + 1]
                                                   [(j_darts_counter_temp11120)]
                                                   [k_darts_counter_temp11120]
                                                   [m_darts_counter_temp11120]
                                                - flux[i_darts_counter_temp11120 - 1]
                                                      [(j_darts_counter_temp11120)]
                                                      [k_darts_counter_temp11120]
                                                      [m_darts_counter_temp11120]);
                                }
                                (*m) = m_darts_counter_temp11120;
                            }
                        }
                        (*i) = i_darts_counter_temp11120;
                    }
                    (*(*L2)) = nx - 1;
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp11120 = (*i);
                        for (; i_darts_counter_temp11120 <= (*(*L2)); i_darts_counter_temp11120++) {
                            (*(*tmp)) = 1.
                                / u[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                                   [k_darts_counter_temp11120][0];
                            (*(*u21i)) = (*(*tmp))
                                * u[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                                   [k_darts_counter_temp11120][1];
                            (*(*u31i)) = (*(*tmp))
                                * u[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                                   [k_darts_counter_temp11120][2];
                            (*(*u41i)) = (*(*tmp))
                                * u[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                                   [k_darts_counter_temp11120][3];
                            (*(*u51i)) = (*(*tmp))
                                * u[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                                   [k_darts_counter_temp11120][4];
                            (*(*tmp)) = 1.
                                / u[i_darts_counter_temp11120 - 1][(j_darts_counter_temp11120)]
                                   [k_darts_counter_temp11120][0];
                            (*(*u21im1)) = (*(*tmp))
                                * u[i_darts_counter_temp11120 - 1][(j_darts_counter_temp11120)]
                                   [k_darts_counter_temp11120][1];
                            (*(*u31im1)) = (*(*tmp))
                                * u[i_darts_counter_temp11120 - 1][(j_darts_counter_temp11120)]
                                   [k_darts_counter_temp11120][2];
                            (*(*u41im1)) = (*(*tmp))
                                * u[i_darts_counter_temp11120 - 1][(j_darts_counter_temp11120)]
                                   [k_darts_counter_temp11120][3];
                            (*(*u51im1)) = (*(*tmp))
                                * u[i_darts_counter_temp11120 - 1][(j_darts_counter_temp11120)]
                                   [k_darts_counter_temp11120][4];
                            flux[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                                [k_darts_counter_temp11120][1]
                                = (4. / 3.) * tx3 * ((*(*u21i)) - (*(*u21im1)));
                            flux[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                                [k_darts_counter_temp11120][2]
                                = tx3 * ((*(*u31i)) - (*(*u31im1)));
                            flux[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                                [k_darts_counter_temp11120][3]
                                = tx3 * ((*(*u41i)) - (*(*u41im1)));
                            flux[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                                [k_darts_counter_temp11120][4]
                                = 0.5 * (1. - 1.3999999999999999 * 1.3999999999999999) * tx3
                                    * (((((*(*u21i))) * ((*(*u21i))))
                                           + (((*(*u31i))) * ((*(*u31i))))
                                           + (((*(*u41i))) * ((*(*u41i)))))
                                        - ((((*(*u21im1))) * ((*(*u21im1))))
                                            + (((*(*u31im1))) * ((*(*u31im1))))
                                            + (((*(*u41im1))) * ((*(*u41im1))))))
                                + (1. / 6.) * tx3
                                    * ((((*(*u21i))) * ((*(*u21i))))
                                        - (((*(*u21im1))) * ((*(*u21im1)))))
                                + 1.3999999999999999 * 1.3999999999999999 * tx3
                                    * ((*(*u51i)) - (*(*u51im1)));
                        }
                        (*i) = i_darts_counter_temp11120;
                    }
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp11120 = (*i);
                        for (; i_darts_counter_temp11120 <= iend; i_darts_counter_temp11120++) {
                            rsd[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                               [k_darts_counter_temp11120][0]
                                = rsd[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                                     [k_darts_counter_temp11120][0]
                                + dx1 * tx1
                                    * (u[i_darts_counter_temp11120 - 1][(j_darts_counter_temp11120)]
                                        [k_darts_counter_temp11120][0]
                                        - 2.
                                            * u[i_darts_counter_temp11120]
                                               [(j_darts_counter_temp11120)]
                                               [k_darts_counter_temp11120][0]
                                        + u[i_darts_counter_temp11120 + 1]
                                           [(j_darts_counter_temp11120)][k_darts_counter_temp11120]
                                           [0]);
                            rsd[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                               [k_darts_counter_temp11120][1]
                                = rsd[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                                     [k_darts_counter_temp11120][1]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp11120 + 1][(
                                           j_darts_counter_temp11120)][k_darts_counter_temp11120][1]
                                        - flux[i_darts_counter_temp11120]
                                              [(j_darts_counter_temp11120)]
                                              [k_darts_counter_temp11120][1])
                                + dx2 * tx1
                                    * (u[i_darts_counter_temp11120 - 1][(j_darts_counter_temp11120)]
                                        [k_darts_counter_temp11120][1]
                                        - 2.
                                            * u[i_darts_counter_temp11120]
                                               [(j_darts_counter_temp11120)]
                                               [k_darts_counter_temp11120][1]
                                        + u[i_darts_counter_temp11120 + 1]
                                           [(j_darts_counter_temp11120)][k_darts_counter_temp11120]
                                           [1]);
                            rsd[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                               [k_darts_counter_temp11120][2]
                                = rsd[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                                     [k_darts_counter_temp11120][2]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp11120 + 1][(
                                           j_darts_counter_temp11120)][k_darts_counter_temp11120][2]
                                        - flux[i_darts_counter_temp11120]
                                              [(j_darts_counter_temp11120)]
                                              [k_darts_counter_temp11120][2])
                                + dx3 * tx1
                                    * (u[i_darts_counter_temp11120 - 1][(j_darts_counter_temp11120)]
                                        [k_darts_counter_temp11120][2]
                                        - 2.
                                            * u[i_darts_counter_temp11120]
                                               [(j_darts_counter_temp11120)]
                                               [k_darts_counter_temp11120][2]
                                        + u[i_darts_counter_temp11120 + 1]
                                           [(j_darts_counter_temp11120)][k_darts_counter_temp11120]
                                           [2]);
                            rsd[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                               [k_darts_counter_temp11120][3]
                                = rsd[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                                     [k_darts_counter_temp11120][3]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp11120 + 1][(
                                           j_darts_counter_temp11120)][k_darts_counter_temp11120][3]
                                        - flux[i_darts_counter_temp11120]
                                              [(j_darts_counter_temp11120)]
                                              [k_darts_counter_temp11120][3])
                                + dx4 * tx1
                                    * (u[i_darts_counter_temp11120 - 1][(j_darts_counter_temp11120)]
                                        [k_darts_counter_temp11120][3]
                                        - 2.
                                            * u[i_darts_counter_temp11120]
                                               [(j_darts_counter_temp11120)]
                                               [k_darts_counter_temp11120][3]
                                        + u[i_darts_counter_temp11120 + 1]
                                           [(j_darts_counter_temp11120)][k_darts_counter_temp11120]
                                           [3]);
                            rsd[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                               [k_darts_counter_temp11120][4]
                                = rsd[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                                     [k_darts_counter_temp11120][4]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp11120 + 1][(
                                           j_darts_counter_temp11120)][k_darts_counter_temp11120][4]
                                        - flux[i_darts_counter_temp11120]
                                              [(j_darts_counter_temp11120)]
                                              [k_darts_counter_temp11120][4])
                                + dx5 * tx1
                                    * (u[i_darts_counter_temp11120 - 1][(j_darts_counter_temp11120)]
                                        [k_darts_counter_temp11120][4]
                                        - 2.
                                            * u[i_darts_counter_temp11120]
                                               [(j_darts_counter_temp11120)]
                                               [k_darts_counter_temp11120][4]
                                        + u[i_darts_counter_temp11120 + 1]
                                           [(j_darts_counter_temp11120)][k_darts_counter_temp11120]
                                           [4]);
                        }
                        (*i) = i_darts_counter_temp11120;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp11120 = (*m);
                        for (; m_darts_counter_temp11120 < 5; m_darts_counter_temp11120++) {
                            rsd[1][(j_darts_counter_temp11120)][k_darts_counter_temp11120]
                               [m_darts_counter_temp11120]
                                = rsd[1][(j_darts_counter_temp11120)][k_darts_counter_temp11120]
                                     [m_darts_counter_temp11120]
                                - dssp
                                    * (+5.
                                            * u[1][(j_darts_counter_temp11120)]
                                               [k_darts_counter_temp11120]
                                               [m_darts_counter_temp11120]
                                        - 4.
                                            * u[2][(j_darts_counter_temp11120)]
                                               [k_darts_counter_temp11120]
                                               [m_darts_counter_temp11120]
                                        + u[3][(j_darts_counter_temp11120)]
                                           [k_darts_counter_temp11120][m_darts_counter_temp11120]);
                            rsd[2][(j_darts_counter_temp11120)][k_darts_counter_temp11120]
                               [m_darts_counter_temp11120]
                                = rsd[2][(j_darts_counter_temp11120)][k_darts_counter_temp11120]
                                     [m_darts_counter_temp11120]
                                - dssp
                                    * (-4.
                                            * u[1][(j_darts_counter_temp11120)]
                                               [k_darts_counter_temp11120]
                                               [m_darts_counter_temp11120]
                                        + 6.
                                            * u[2][(j_darts_counter_temp11120)]
                                               [k_darts_counter_temp11120]
                                               [m_darts_counter_temp11120]
                                        - 4.
                                            * u[3][(j_darts_counter_temp11120)]
                                               [k_darts_counter_temp11120]
                                               [m_darts_counter_temp11120]
                                        + u[4][(j_darts_counter_temp11120)]
                                           [k_darts_counter_temp11120][m_darts_counter_temp11120]);
                        }
                        (*m) = m_darts_counter_temp11120;
                    }
                    (*(*ist1)) = 3;
                    (*(*iend1)) = nx - 4;
                    {
                        /*Loop's init*/
                        (*i) = (*(*ist1));
                        int i_darts_counter_temp11120 = (*i);
                        for (; i_darts_counter_temp11120 <= (*(*iend1));
                             i_darts_counter_temp11120++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp11120 = (*m);
                                for (; m_darts_counter_temp11120 < 5; m_darts_counter_temp11120++) {
                                    rsd[i_darts_counter_temp11120][(j_darts_counter_temp11120)]
                                       [k_darts_counter_temp11120][m_darts_counter_temp11120]
                                        = rsd[i_darts_counter_temp11120]
                                             [(j_darts_counter_temp11120)]
                                             [k_darts_counter_temp11120][m_darts_counter_temp11120]
                                        - dssp
                                            * (u[i_darts_counter_temp11120 - 2]
                                                [(j_darts_counter_temp11120)]
                                                [k_darts_counter_temp11120]
                                                [m_darts_counter_temp11120]
                                                - 4.
                                                    * u[i_darts_counter_temp11120 - 1]
                                                       [(j_darts_counter_temp11120)]
                                                       [k_darts_counter_temp11120]
                                                       [m_darts_counter_temp11120]
                                                + 6.
                                                    * u[i_darts_counter_temp11120]
                                                       [(j_darts_counter_temp11120)]
                                                       [k_darts_counter_temp11120]
                                                       [m_darts_counter_temp11120]
                                                - 4.
                                                    * u[i_darts_counter_temp11120 + 1]
                                                       [(j_darts_counter_temp11120)]
                                                       [k_darts_counter_temp11120]
                                                       [m_darts_counter_temp11120]
                                                + u[i_darts_counter_temp11120 + 2]
                                                   [(j_darts_counter_temp11120)]
                                                   [k_darts_counter_temp11120]
                                                   [m_darts_counter_temp11120]);
                                }
                                (*m) = m_darts_counter_temp11120;
                            }
                        }
                        (*i) = i_darts_counter_temp11120;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp11120 = (*m);
                        for (; m_darts_counter_temp11120 < 5; m_darts_counter_temp11120++) {
                            rsd[nx - 3][(j_darts_counter_temp11120)][k_darts_counter_temp11120]
                               [m_darts_counter_temp11120]
                                = rsd[nx - 3][(j_darts_counter_temp11120)]
                                     [k_darts_counter_temp11120][m_darts_counter_temp11120]
                                - dssp
                                    * (u[nx - 5][(j_darts_counter_temp11120)]
                                        [k_darts_counter_temp11120][m_darts_counter_temp11120]
                                        - 4.
                                            * u[nx - 4][(j_darts_counter_temp11120)]
                                               [k_darts_counter_temp11120]
                                               [m_darts_counter_temp11120]
                                        + 6.
                                            * u[nx - 3][(j_darts_counter_temp11120)]
                                               [k_darts_counter_temp11120]
                                               [m_darts_counter_temp11120]
                                        - 4.
                                            * u[nx - 2][(j_darts_counter_temp11120)]
                                               [k_darts_counter_temp11120]
                                               [m_darts_counter_temp11120]);
                            rsd[nx - 2][(j_darts_counter_temp11120)][k_darts_counter_temp11120]
                               [m_darts_counter_temp11120]
                                = rsd[nx - 2][(j_darts_counter_temp11120)]
                                     [k_darts_counter_temp11120][m_darts_counter_temp11120]
                                - dssp
                                    * (u[nx - 4][(j_darts_counter_temp11120)]
                                        [k_darts_counter_temp11120][m_darts_counter_temp11120]
                                        - 4.
                                            * u[nx - 3][(j_darts_counter_temp11120)]
                                               [k_darts_counter_temp11120]
                                               [m_darts_counter_temp11120]
                                        + 5.
                                            * u[nx - 2][(j_darts_counter_temp11120)]
                                               [k_darts_counter_temp11120]
                                               [m_darts_counter_temp11120]);
                        }
                        (*m) = m_darts_counter_temp11120;
                    }
                }
                (*k) = k_darts_counter_temp11120;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets11120[0].decDep();
}
TP11120::TP11120(int in_numThreads, int in_mainCodeletID, TP10896* in_TPParent,
    int in_initIteration, int in_lastIteration, TP11120** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L2_darts11120(new int*[this->numThreads])
    , i_darts11120(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iend1_darts11120(new int*[this->numThreads])
    , ist1_darts11120(new int*[this->numThreads])
    , j_darts11120(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts11120(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts11120(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts11120(new double*[this->numThreads])
    , u21i_darts11120(new double*[this->numThreads])
    , u21im1_darts11120(new double*[this->numThreads])
    , u31i_darts11120(new double*[this->numThreads])
    , u31im1_darts11120(new double*[this->numThreads])
    , u41i_darts11120(new double*[this->numThreads])
    , u41im1_darts11120(new double*[this->numThreads])
    , u51i_darts11120(new double*[this->numThreads])
    , u51im1_darts11120(new double*[this->numThreads])
    , initIteration11120(in_initIteration)
    , lastIteration11120(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets11120(new _barrierCodelets11120[1])
    , checkInCodelets11121(new _checkInCodelets11121[this->numThreads])
{
    /*Initialize the loop parameters*/
    range11120 = abs(lastIteration11120 - initIteration11120) / 1;
    rangePerCodelet11120 = range11120 / numThreads;
    minIteration11120 = min<int>(lastIteration11120, initIteration11120);
    remainderRange11120 = range11120 % numThreads;
    /*Initialize inputs and vars.*/
    this->L2_darts11120 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->iend1_darts11120 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ist1_darts11120 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts11120
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21i_darts11120
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21im1_darts11120
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31i_darts11120
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31im1_darts11120
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41i_darts11120
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41im1_darts11120
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51i_darts11120
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51im1_darts11120
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets11120[0] = _barrierCodelets11120(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets11121* checkInCodelets11121Ptr = (this->checkInCodelets11121);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets11121);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets11121Ptr) = _checkInCodelets11121(2, 1, this, codeletCounter);
#else
        (*checkInCodelets11121Ptr) = _checkInCodelets11121(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets11121Ptr).decDep();
        checkInCodelets11121Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP11120::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets11121[localID].setID(codeletID);
    this->checkInCodelets11121[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets11121[localID + this->baseNumThreads * i]
            = _checkInCodelets11121(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets11121[localID + this->baseNumThreads * i]
            = _checkInCodelets11121(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets11121[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets11121[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP11120::~TP11120()
{
    delete[] L2_darts11120;
    delete[] iend1_darts11120;
    delete[] ist1_darts11120;
    delete[] tmp_darts11120;
    delete[] u21i_darts11120;
    delete[] u21im1_darts11120;
    delete[] u31i_darts11120;
    delete[] u31im1_darts11120;
    delete[] u41i_darts11120;
    delete[] u41im1_darts11120;
    delete[] u51i_darts11120;
    delete[] u51im1_darts11120;
    delete[] barrierCodelets11120;
    delete[] checkInCodelets11121;
}
/*TP11760: OMPForDirective*/
void TP11760::_barrierCodelets11760::fire(void)
{
    TP11760* myTP = static_cast<TP11760*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets11760[0].decDep();
}
bool TP11760::requestNewRangeIterations11760(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet11760 * codeletID;
        int tempEndRange = rangePerCodelet11760 * (codeletID + 1);
        if (remainderRange11760 != 0) {
            if (codeletID < (uint32_t)remainderRange11760) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange11760;
                tempEndRange += remainderRange11760;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration11760;
        tempEndRange = tempEndRange * 1 + minIteration11760;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration11760 < lastIteration11760) {
            (this->inputsTPParent->i_darts11760[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts11760[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration11760;
        }
    }
    return isThereNewIteration;
}
void TP11760::_checkInCodelets11761::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L1_darts11760[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L1_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L2_darts11760[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->q_darts11760[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31_darts11760[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31_darts10896[this->getID()]);

    /*printing node 11761: ForStmt*/
    /*var: L1*/
    /*var: L2*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: q*/
    /*var: u31*/
    int** L1 = &(this->inputsTPParent->L1_darts11760[this->getLocalID()]);
    (void)L1 /*OMP_SHARED_PRIVATE*/;
    int** L2 = &(this->inputsTPParent->L2_darts11760[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts11760[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts11760[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts11760[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts11760[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** u31 = &(this->inputsTPParent->u31_darts11760[this->getLocalID()]);
    (void)u31 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations11760(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets11760[0].decDep();
        return;
    }
    for (int i_darts_counter_temp11760 = (*i); i_darts_counter_temp11760 <= endRange
         && i_darts_counter_temp11760 <= this->inputsTPParent->lastIteration11760;
         i_darts_counter_temp11760++) {
        {
            {
                /*Loop's init*/
                (*j) = (*(*L1));
                int j_darts_counter_temp11760 = (*j);
                for (; j_darts_counter_temp11760 <= (*(*L2)); j_darts_counter_temp11760++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp11760 = (*k);
                        for (; k_darts_counter_temp11760 <= nz - 2; k_darts_counter_temp11760++) {
                            flux[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                [k_darts_counter_temp11760][0]
                                = u[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                   [k_darts_counter_temp11760][2];
                            (*(*u31)) = u[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                         [k_darts_counter_temp11760][2]
                                / u[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                   [k_darts_counter_temp11760][0];
                            (*(*q)) = 0.5
                                * (u[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                    [k_darts_counter_temp11760][1]
                                        * u[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                           [k_darts_counter_temp11760][1]
                                    + u[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                       [k_darts_counter_temp11760][2]
                                        * u[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                           [k_darts_counter_temp11760][2]
                                    + u[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                       [k_darts_counter_temp11760][3]
                                        * u[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                           [k_darts_counter_temp11760][3])
                                / u[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                   [k_darts_counter_temp11760][0];
                            flux[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                [k_darts_counter_temp11760][1]
                                = u[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                   [k_darts_counter_temp11760][1]
                                * (*(*u31));
                            flux[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                [k_darts_counter_temp11760][2]
                                = u[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                   [k_darts_counter_temp11760][2]
                                    * (*(*u31))
                                + 0.40000000000000002
                                    * (u[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                        [k_darts_counter_temp11760][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                [k_darts_counter_temp11760][3]
                                = u[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                   [k_darts_counter_temp11760][3]
                                * (*(*u31));
                            flux[(i_darts_counter_temp11760)][j_darts_counter_temp11760]
                                [k_darts_counter_temp11760][4]
                                = (1.3999999999999999
                                          * u[(i_darts_counter_temp11760)]
                                             [j_darts_counter_temp11760][k_darts_counter_temp11760]
                                             [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u31));
                        }
                        (*k) = k_darts_counter_temp11760;
                    }
                }
                (*j) = j_darts_counter_temp11760;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets11760[0].decDep();
}
TP11760::TP11760(int in_numThreads, int in_mainCodeletID, TP10896* in_TPParent,
    int in_initIteration, int in_lastIteration, TP11760** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts11760(new int*[this->numThreads])
    , L2_darts11760(new int*[this->numThreads])
    , i_darts11760(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts11760(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts11760(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts11760(new double*[this->numThreads])
    , u31_darts11760(new double*[this->numThreads])
    , initIteration11760(in_initIteration)
    , lastIteration11760(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets11760(new _barrierCodelets11760[1])
    , checkInCodelets11761(new _checkInCodelets11761[this->numThreads])
{
    /*Initialize the loop parameters*/
    range11760 = abs(lastIteration11760 - initIteration11760) / 1;
    rangePerCodelet11760 = range11760 / numThreads;
    minIteration11760 = min<int>(lastIteration11760, initIteration11760);
    remainderRange11760 = range11760 % numThreads;
    /*Initialize inputs and vars.*/
    this->L1_darts11760 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->L2_darts11760 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->q_darts11760
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31_darts11760
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets11760[0] = _barrierCodelets11760(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets11761* checkInCodelets11761Ptr = (this->checkInCodelets11761);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets11761);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets11761Ptr) = _checkInCodelets11761(2, 1, this, codeletCounter);
#else
        (*checkInCodelets11761Ptr) = _checkInCodelets11761(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets11761Ptr).decDep();
        checkInCodelets11761Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP11760::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets11761[localID].setID(codeletID);
    this->checkInCodelets11761[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets11761[localID + this->baseNumThreads * i]
            = _checkInCodelets11761(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets11761[localID + this->baseNumThreads * i]
            = _checkInCodelets11761(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets11761[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets11761[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP11760::~TP11760()
{
    delete[] L1_darts11760;
    delete[] L2_darts11760;
    delete[] q_darts11760;
    delete[] u31_darts11760;
    delete[] barrierCodelets11760;
    delete[] checkInCodelets11761;
}
/*TP11909: OMPForDirective*/
void TP11909::_barrierCodelets11909::fire(void)
{
    TP11909* myTP = static_cast<TP11909*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets11909[0].decDep();
}
bool TP11909::requestNewRangeIterations11909(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet11909 * codeletID;
        int tempEndRange = rangePerCodelet11909 * (codeletID + 1);
        if (remainderRange11909 != 0) {
            if (codeletID < (uint32_t)remainderRange11909) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange11909;
                tempEndRange += remainderRange11909;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration11909;
        tempEndRange = tempEndRange * 1 + minIteration11909;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration11909 < lastIteration11909) {
            (this->inputsTPParent->i_darts11909[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts11909[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration11909;
        }
    }
    return isThereNewIteration;
}
void TP11909::_checkInCodelets11910::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L2_darts11909[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jend1_darts11909[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jend1_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jst1_darts11909[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jst1_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts11909[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21j_darts11909[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21j_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21jm1_darts11909[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21jm1_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31j_darts11909[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31j_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31jm1_darts11909[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31jm1_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41j_darts11909[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41j_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41jm1_darts11909[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41jm1_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51j_darts11909[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51j_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51jm1_darts11909[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51jm1_darts10896[this->getID()]);

    /*printing node 11910: ForStmt*/
    /*var: L2*/
    /*var: i*/
    /*var: j*/
    /*var: jend1*/
    /*var: jst1*/
    /*var: k*/
    /*var: m*/
    /*var: tmp*/
    /*var: u21j*/
    /*var: u21jm1*/
    /*var: u31j*/
    /*var: u31jm1*/
    /*var: u41j*/
    /*var: u41jm1*/
    /*var: u51j*/
    /*var: u51jm1*/
    int** L2 = &(this->inputsTPParent->L2_darts11909[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts11909[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts11909[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jend1 = &(this->inputsTPParent->jend1_darts11909[this->getLocalID()]);
    (void)jend1 /*OMP_SHARED_PRIVATE*/;
    int** jst1 = &(this->inputsTPParent->jst1_darts11909[this->getLocalID()]);
    (void)jst1 /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts11909[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts11909[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts11909[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21j = &(this->inputsTPParent->u21j_darts11909[this->getLocalID()]);
    (void)u21j /*OMP_SHARED_PRIVATE*/;
    double** u21jm1 = &(this->inputsTPParent->u21jm1_darts11909[this->getLocalID()]);
    (void)u21jm1 /*OMP_SHARED_PRIVATE*/;
    double** u31j = &(this->inputsTPParent->u31j_darts11909[this->getLocalID()]);
    (void)u31j /*OMP_SHARED_PRIVATE*/;
    double** u31jm1 = &(this->inputsTPParent->u31jm1_darts11909[this->getLocalID()]);
    (void)u31jm1 /*OMP_SHARED_PRIVATE*/;
    double** u41j = &(this->inputsTPParent->u41j_darts11909[this->getLocalID()]);
    (void)u41j /*OMP_SHARED_PRIVATE*/;
    double** u41jm1 = &(this->inputsTPParent->u41jm1_darts11909[this->getLocalID()]);
    (void)u41jm1 /*OMP_SHARED_PRIVATE*/;
    double** u51j = &(this->inputsTPParent->u51j_darts11909[this->getLocalID()]);
    (void)u51j /*OMP_SHARED_PRIVATE*/;
    double** u51jm1 = &(this->inputsTPParent->u51jm1_darts11909[this->getLocalID()]);
    (void)u51jm1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations11909(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets11909[0].decDep();
        return;
    }
    for (int i_darts_counter_temp11909 = (*i); i_darts_counter_temp11909 <= endRange
         && i_darts_counter_temp11909 <= this->inputsTPParent->lastIteration11909;
         i_darts_counter_temp11909++) {
        {
            {
                /*Loop's init*/
                (*k) = 1;
                int k_darts_counter_temp11909 = (*k);
                for (; k_darts_counter_temp11909 <= nz - 2; k_darts_counter_temp11909++) {
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp11909 = (*j);
                        for (; j_darts_counter_temp11909 <= jend; j_darts_counter_temp11909++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp11909 = (*m);
                                for (; m_darts_counter_temp11909 < 5; m_darts_counter_temp11909++) {
                                    rsd[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                                       [k_darts_counter_temp11909][m_darts_counter_temp11909]
                                        = rsd[(i_darts_counter_temp11909)]
                                             [j_darts_counter_temp11909][k_darts_counter_temp11909]
                                             [m_darts_counter_temp11909]
                                        - ty2
                                            * (flux[(i_darts_counter_temp11909)]
                                                   [j_darts_counter_temp11909 + 1]
                                                   [k_darts_counter_temp11909]
                                                   [m_darts_counter_temp11909]
                                                - flux[(i_darts_counter_temp11909)]
                                                      [j_darts_counter_temp11909 - 1]
                                                      [k_darts_counter_temp11909]
                                                      [m_darts_counter_temp11909]);
                                }
                                (*m) = m_darts_counter_temp11909;
                            }
                        }
                        (*j) = j_darts_counter_temp11909;
                    }
                    (*(*L2)) = ny - 1;
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp11909 = (*j);
                        for (; j_darts_counter_temp11909 <= (*(*L2)); j_darts_counter_temp11909++) {
                            (*(*tmp)) = 1.
                                / u[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                                   [k_darts_counter_temp11909][0];
                            (*(*u21j)) = (*(*tmp))
                                * u[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                                   [k_darts_counter_temp11909][1];
                            (*(*u31j)) = (*(*tmp))
                                * u[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                                   [k_darts_counter_temp11909][2];
                            (*(*u41j)) = (*(*tmp))
                                * u[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                                   [k_darts_counter_temp11909][3];
                            (*(*u51j)) = (*(*tmp))
                                * u[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                                   [k_darts_counter_temp11909][4];
                            (*(*tmp)) = 1.
                                / u[(i_darts_counter_temp11909)][j_darts_counter_temp11909 - 1]
                                   [k_darts_counter_temp11909][0];
                            (*(*u21jm1)) = (*(*tmp))
                                * u[(i_darts_counter_temp11909)][j_darts_counter_temp11909 - 1]
                                   [k_darts_counter_temp11909][1];
                            (*(*u31jm1)) = (*(*tmp))
                                * u[(i_darts_counter_temp11909)][j_darts_counter_temp11909 - 1]
                                   [k_darts_counter_temp11909][2];
                            (*(*u41jm1)) = (*(*tmp))
                                * u[(i_darts_counter_temp11909)][j_darts_counter_temp11909 - 1]
                                   [k_darts_counter_temp11909][3];
                            (*(*u51jm1)) = (*(*tmp))
                                * u[(i_darts_counter_temp11909)][j_darts_counter_temp11909 - 1]
                                   [k_darts_counter_temp11909][4];
                            flux[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                                [k_darts_counter_temp11909][1]
                                = ty3 * ((*(*u21j)) - (*(*u21jm1)));
                            flux[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                                [k_darts_counter_temp11909][2]
                                = (4. / 3.) * ty3 * ((*(*u31j)) - (*(*u31jm1)));
                            flux[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                                [k_darts_counter_temp11909][3]
                                = ty3 * ((*(*u41j)) - (*(*u41jm1)));
                            flux[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                                [k_darts_counter_temp11909][4]
                                = 0.5 * (1. - 1.3999999999999999 * 1.3999999999999999) * ty3
                                    * (((((*(*u21j))) * ((*(*u21j))))
                                           + (((*(*u31j))) * ((*(*u31j))))
                                           + (((*(*u41j))) * ((*(*u41j)))))
                                        - ((((*(*u21jm1))) * ((*(*u21jm1))))
                                            + (((*(*u31jm1))) * ((*(*u31jm1))))
                                            + (((*(*u41jm1))) * ((*(*u41jm1))))))
                                + (1. / 6.) * ty3
                                    * ((((*(*u31j))) * ((*(*u31j))))
                                        - (((*(*u31jm1))) * ((*(*u31jm1)))))
                                + 1.3999999999999999 * 1.3999999999999999 * ty3
                                    * ((*(*u51j)) - (*(*u51jm1)));
                        }
                        (*j) = j_darts_counter_temp11909;
                    }
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp11909 = (*j);
                        for (; j_darts_counter_temp11909 <= jend; j_darts_counter_temp11909++) {
                            rsd[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                               [k_darts_counter_temp11909][0]
                                = rsd[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                                     [k_darts_counter_temp11909][0]
                                + dy1 * ty1
                                    * (u[(i_darts_counter_temp11909)][j_darts_counter_temp11909 - 1]
                                        [k_darts_counter_temp11909][0]
                                        - 2.
                                            * u[(i_darts_counter_temp11909)]
                                               [j_darts_counter_temp11909]
                                               [k_darts_counter_temp11909][0]
                                        + u[(i_darts_counter_temp11909)][j_darts_counter_temp11909
                                            + 1][k_darts_counter_temp11909][0]);
                            rsd[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                               [k_darts_counter_temp11909][1]
                                = rsd[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                                     [k_darts_counter_temp11909][1]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp11909)][j_darts_counter_temp11909
                                           + 1][k_darts_counter_temp11909][1]
                                        - flux[(i_darts_counter_temp11909)]
                                              [j_darts_counter_temp11909][k_darts_counter_temp11909]
                                              [1])
                                + dy2 * ty1
                                    * (u[(i_darts_counter_temp11909)][j_darts_counter_temp11909 - 1]
                                        [k_darts_counter_temp11909][1]
                                        - 2.
                                            * u[(i_darts_counter_temp11909)]
                                               [j_darts_counter_temp11909]
                                               [k_darts_counter_temp11909][1]
                                        + u[(i_darts_counter_temp11909)][j_darts_counter_temp11909
                                            + 1][k_darts_counter_temp11909][1]);
                            rsd[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                               [k_darts_counter_temp11909][2]
                                = rsd[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                                     [k_darts_counter_temp11909][2]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp11909)][j_darts_counter_temp11909
                                           + 1][k_darts_counter_temp11909][2]
                                        - flux[(i_darts_counter_temp11909)]
                                              [j_darts_counter_temp11909][k_darts_counter_temp11909]
                                              [2])
                                + dy3 * ty1
                                    * (u[(i_darts_counter_temp11909)][j_darts_counter_temp11909 - 1]
                                        [k_darts_counter_temp11909][2]
                                        - 2.
                                            * u[(i_darts_counter_temp11909)]
                                               [j_darts_counter_temp11909]
                                               [k_darts_counter_temp11909][2]
                                        + u[(i_darts_counter_temp11909)][j_darts_counter_temp11909
                                            + 1][k_darts_counter_temp11909][2]);
                            rsd[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                               [k_darts_counter_temp11909][3]
                                = rsd[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                                     [k_darts_counter_temp11909][3]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp11909)][j_darts_counter_temp11909
                                           + 1][k_darts_counter_temp11909][3]
                                        - flux[(i_darts_counter_temp11909)]
                                              [j_darts_counter_temp11909][k_darts_counter_temp11909]
                                              [3])
                                + dy4 * ty1
                                    * (u[(i_darts_counter_temp11909)][j_darts_counter_temp11909 - 1]
                                        [k_darts_counter_temp11909][3]
                                        - 2.
                                            * u[(i_darts_counter_temp11909)]
                                               [j_darts_counter_temp11909]
                                               [k_darts_counter_temp11909][3]
                                        + u[(i_darts_counter_temp11909)][j_darts_counter_temp11909
                                            + 1][k_darts_counter_temp11909][3]);
                            rsd[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                               [k_darts_counter_temp11909][4]
                                = rsd[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                                     [k_darts_counter_temp11909][4]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp11909)][j_darts_counter_temp11909
                                           + 1][k_darts_counter_temp11909][4]
                                        - flux[(i_darts_counter_temp11909)]
                                              [j_darts_counter_temp11909][k_darts_counter_temp11909]
                                              [4])
                                + dy5 * ty1
                                    * (u[(i_darts_counter_temp11909)][j_darts_counter_temp11909 - 1]
                                        [k_darts_counter_temp11909][4]
                                        - 2.
                                            * u[(i_darts_counter_temp11909)]
                                               [j_darts_counter_temp11909]
                                               [k_darts_counter_temp11909][4]
                                        + u[(i_darts_counter_temp11909)][j_darts_counter_temp11909
                                            + 1][k_darts_counter_temp11909][4]);
                        }
                        (*j) = j_darts_counter_temp11909;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp11909 = (*m);
                        for (; m_darts_counter_temp11909 < 5; m_darts_counter_temp11909++) {
                            rsd[(i_darts_counter_temp11909)][1][k_darts_counter_temp11909]
                               [m_darts_counter_temp11909]
                                = rsd[(i_darts_counter_temp11909)][1][k_darts_counter_temp11909]
                                     [m_darts_counter_temp11909]
                                - dssp
                                    * (+5.
                                            * u[(i_darts_counter_temp11909)][1]
                                               [k_darts_counter_temp11909]
                                               [m_darts_counter_temp11909]
                                        - 4.
                                            * u[(i_darts_counter_temp11909)][2]
                                               [k_darts_counter_temp11909]
                                               [m_darts_counter_temp11909]
                                        + u[(i_darts_counter_temp11909)][3]
                                           [k_darts_counter_temp11909][m_darts_counter_temp11909]);
                            rsd[(i_darts_counter_temp11909)][2][k_darts_counter_temp11909]
                               [m_darts_counter_temp11909]
                                = rsd[(i_darts_counter_temp11909)][2][k_darts_counter_temp11909]
                                     [m_darts_counter_temp11909]
                                - dssp
                                    * (-4.
                                            * u[(i_darts_counter_temp11909)][1]
                                               [k_darts_counter_temp11909]
                                               [m_darts_counter_temp11909]
                                        + 6.
                                            * u[(i_darts_counter_temp11909)][2]
                                               [k_darts_counter_temp11909]
                                               [m_darts_counter_temp11909]
                                        - 4.
                                            * u[(i_darts_counter_temp11909)][3]
                                               [k_darts_counter_temp11909]
                                               [m_darts_counter_temp11909]
                                        + u[(i_darts_counter_temp11909)][4]
                                           [k_darts_counter_temp11909][m_darts_counter_temp11909]);
                        }
                        (*m) = m_darts_counter_temp11909;
                    }
                    (*(*jst1)) = 3;
                    (*(*jend1)) = ny - 4;
                    {
                        /*Loop's init*/
                        (*j) = (*(*jst1));
                        int j_darts_counter_temp11909 = (*j);
                        for (; j_darts_counter_temp11909 <= (*(*jend1));
                             j_darts_counter_temp11909++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp11909 = (*m);
                                for (; m_darts_counter_temp11909 < 5; m_darts_counter_temp11909++) {
                                    rsd[(i_darts_counter_temp11909)][j_darts_counter_temp11909]
                                       [k_darts_counter_temp11909][m_darts_counter_temp11909]
                                        = rsd[(i_darts_counter_temp11909)]
                                             [j_darts_counter_temp11909][k_darts_counter_temp11909]
                                             [m_darts_counter_temp11909]
                                        - dssp
                                            * (u[(i_darts_counter_temp11909)]
                                                [j_darts_counter_temp11909 - 2]
                                                [k_darts_counter_temp11909]
                                                [m_darts_counter_temp11909]
                                                - 4.
                                                    * u[(i_darts_counter_temp11909)]
                                                       [j_darts_counter_temp11909 - 1]
                                                       [k_darts_counter_temp11909]
                                                       [m_darts_counter_temp11909]
                                                + 6.
                                                    * u[(i_darts_counter_temp11909)]
                                                       [j_darts_counter_temp11909]
                                                       [k_darts_counter_temp11909]
                                                       [m_darts_counter_temp11909]
                                                - 4.
                                                    * u[(i_darts_counter_temp11909)]
                                                       [j_darts_counter_temp11909 + 1]
                                                       [k_darts_counter_temp11909]
                                                       [m_darts_counter_temp11909]
                                                + u[(i_darts_counter_temp11909)]
                                                   [j_darts_counter_temp11909 + 2]
                                                   [k_darts_counter_temp11909]
                                                   [m_darts_counter_temp11909]);
                                }
                                (*m) = m_darts_counter_temp11909;
                            }
                        }
                        (*j) = j_darts_counter_temp11909;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp11909 = (*m);
                        for (; m_darts_counter_temp11909 < 5; m_darts_counter_temp11909++) {
                            rsd[(i_darts_counter_temp11909)][ny - 3][k_darts_counter_temp11909]
                               [m_darts_counter_temp11909]
                                = rsd[(i_darts_counter_temp11909)][ny - 3]
                                     [k_darts_counter_temp11909][m_darts_counter_temp11909]
                                - dssp
                                    * (u[(i_darts_counter_temp11909)][ny - 5]
                                        [k_darts_counter_temp11909][m_darts_counter_temp11909]
                                        - 4.
                                            * u[(i_darts_counter_temp11909)][ny - 4]
                                               [k_darts_counter_temp11909]
                                               [m_darts_counter_temp11909]
                                        + 6.
                                            * u[(i_darts_counter_temp11909)][ny - 3]
                                               [k_darts_counter_temp11909]
                                               [m_darts_counter_temp11909]
                                        - 4.
                                            * u[(i_darts_counter_temp11909)][ny - 2]
                                               [k_darts_counter_temp11909]
                                               [m_darts_counter_temp11909]);
                            rsd[(i_darts_counter_temp11909)][ny - 2][k_darts_counter_temp11909]
                               [m_darts_counter_temp11909]
                                = rsd[(i_darts_counter_temp11909)][ny - 2]
                                     [k_darts_counter_temp11909][m_darts_counter_temp11909]
                                - dssp
                                    * (u[(i_darts_counter_temp11909)][ny - 4]
                                        [k_darts_counter_temp11909][m_darts_counter_temp11909]
                                        - 4.
                                            * u[(i_darts_counter_temp11909)][ny - 3]
                                               [k_darts_counter_temp11909]
                                               [m_darts_counter_temp11909]
                                        + 5.
                                            * u[(i_darts_counter_temp11909)][ny - 2]
                                               [k_darts_counter_temp11909]
                                               [m_darts_counter_temp11909]);
                        }
                        (*m) = m_darts_counter_temp11909;
                    }
                }
                (*k) = k_darts_counter_temp11909;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets11909[0].decDep();
}
TP11909::TP11909(int in_numThreads, int in_mainCodeletID, TP10896* in_TPParent,
    int in_initIteration, int in_lastIteration, TP11909** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L2_darts11909(new int*[this->numThreads])
    , i_darts11909(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts11909(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jend1_darts11909(new int*[this->numThreads])
    , jst1_darts11909(new int*[this->numThreads])
    , k_darts11909(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts11909(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts11909(new double*[this->numThreads])
    , u21j_darts11909(new double*[this->numThreads])
    , u21jm1_darts11909(new double*[this->numThreads])
    , u31j_darts11909(new double*[this->numThreads])
    , u31jm1_darts11909(new double*[this->numThreads])
    , u41j_darts11909(new double*[this->numThreads])
    , u41jm1_darts11909(new double*[this->numThreads])
    , u51j_darts11909(new double*[this->numThreads])
    , u51jm1_darts11909(new double*[this->numThreads])
    , initIteration11909(in_initIteration)
    , lastIteration11909(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets11909(new _barrierCodelets11909[1])
    , checkInCodelets11910(new _checkInCodelets11910[this->numThreads])
{
    /*Initialize the loop parameters*/
    range11909 = abs(lastIteration11909 - initIteration11909) / 1;
    rangePerCodelet11909 = range11909 / numThreads;
    minIteration11909 = min<int>(lastIteration11909, initIteration11909);
    remainderRange11909 = range11909 % numThreads;
    /*Initialize inputs and vars.*/
    this->L2_darts11909 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jend1_darts11909 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jst1_darts11909 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts11909
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21j_darts11909
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21jm1_darts11909
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31j_darts11909
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31jm1_darts11909
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41j_darts11909
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41jm1_darts11909
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51j_darts11909
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51jm1_darts11909
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets11909[0] = _barrierCodelets11909(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets11910* checkInCodelets11910Ptr = (this->checkInCodelets11910);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets11910);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets11910Ptr) = _checkInCodelets11910(2, 1, this, codeletCounter);
#else
        (*checkInCodelets11910Ptr) = _checkInCodelets11910(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets11910Ptr).decDep();
        checkInCodelets11910Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP11909::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets11910[localID].setID(codeletID);
    this->checkInCodelets11910[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets11910[localID + this->baseNumThreads * i]
            = _checkInCodelets11910(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets11910[localID + this->baseNumThreads * i]
            = _checkInCodelets11910(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets11910[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets11910[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP11909::~TP11909()
{
    delete[] L2_darts11909;
    delete[] jend1_darts11909;
    delete[] jst1_darts11909;
    delete[] tmp_darts11909;
    delete[] u21j_darts11909;
    delete[] u21jm1_darts11909;
    delete[] u31j_darts11909;
    delete[] u31jm1_darts11909;
    delete[] u41j_darts11909;
    delete[] u41jm1_darts11909;
    delete[] u51j_darts11909;
    delete[] u51jm1_darts11909;
    delete[] barrierCodelets11909;
    delete[] checkInCodelets11910;
}
/*TP12546: OMPForDirective*/
void TP12546::_barrierCodelets12546::fire(void)
{
    TP12546* myTP = static_cast<TP12546*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets12546[0].decDep();
}
bool TP12546::requestNewRangeIterations12546(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet12546 * codeletID;
        int tempEndRange = rangePerCodelet12546 * (codeletID + 1);
        if (remainderRange12546 != 0) {
            if (codeletID < (uint32_t)remainderRange12546) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange12546;
                tempEndRange += remainderRange12546;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration12546;
        tempEndRange = tempEndRange * 1 + minIteration12546;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration12546 < lastIteration12546) {
            (this->inputsTPParent->i_darts12546[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts12546[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration12546;
        }
    }
    return isThereNewIteration;
}
void TP12546::_checkInCodelets12547::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->q_darts12546[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts12546[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21k_darts12546[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21k_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21km1_darts12546[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21km1_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31k_darts12546[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31k_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31km1_darts12546[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31km1_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41_darts12546[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41k_darts12546[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41k_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41km1_darts12546[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41km1_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51k_darts12546[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51k_darts10896[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51km1_darts12546[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51km1_darts10896[this->getID()]);

    /*printing node 12547: ForStmt*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    /*var: q*/
    /*var: tmp*/
    /*var: u21k*/
    /*var: u21km1*/
    /*var: u31k*/
    /*var: u31km1*/
    /*var: u41*/
    /*var: u41k*/
    /*var: u41km1*/
    /*var: u51k*/
    /*var: u51km1*/
    int* i = &(this->inputsTPParent->i_darts12546[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts12546[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts12546[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts12546[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts12546[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts12546[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21k = &(this->inputsTPParent->u21k_darts12546[this->getLocalID()]);
    (void)u21k /*OMP_SHARED_PRIVATE*/;
    double** u21km1 = &(this->inputsTPParent->u21km1_darts12546[this->getLocalID()]);
    (void)u21km1 /*OMP_SHARED_PRIVATE*/;
    double** u31k = &(this->inputsTPParent->u31k_darts12546[this->getLocalID()]);
    (void)u31k /*OMP_SHARED_PRIVATE*/;
    double** u31km1 = &(this->inputsTPParent->u31km1_darts12546[this->getLocalID()]);
    (void)u31km1 /*OMP_SHARED_PRIVATE*/;
    double** u41 = &(this->inputsTPParent->u41_darts12546[this->getLocalID()]);
    (void)u41 /*OMP_SHARED_PRIVATE*/;
    double** u41k = &(this->inputsTPParent->u41k_darts12546[this->getLocalID()]);
    (void)u41k /*OMP_SHARED_PRIVATE*/;
    double** u41km1 = &(this->inputsTPParent->u41km1_darts12546[this->getLocalID()]);
    (void)u41km1 /*OMP_SHARED_PRIVATE*/;
    double** u51k = &(this->inputsTPParent->u51k_darts12546[this->getLocalID()]);
    (void)u51k /*OMP_SHARED_PRIVATE*/;
    double** u51km1 = &(this->inputsTPParent->u51km1_darts12546[this->getLocalID()]);
    (void)u51km1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations12546(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets12546[0].decDep();
        return;
    }
    for (int i_darts_counter_temp12546 = (*i); i_darts_counter_temp12546 <= endRange
         && i_darts_counter_temp12546 <= this->inputsTPParent->lastIteration12546;
         i_darts_counter_temp12546++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp12546 = (*j);
                for (; j_darts_counter_temp12546 <= jend; j_darts_counter_temp12546++) {
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp12546 = (*k);
                        for (; k_darts_counter_temp12546 <= nz - 1; k_darts_counter_temp12546++) {
                            flux[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                [k_darts_counter_temp12546][0]
                                = u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                   [k_darts_counter_temp12546][3];
                            (*(*u41)) = u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                         [k_darts_counter_temp12546][3]
                                / u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                   [k_darts_counter_temp12546][0];
                            (*(*q)) = 0.5
                                * (u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                    [k_darts_counter_temp12546][1]
                                        * u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                           [k_darts_counter_temp12546][1]
                                    + u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                       [k_darts_counter_temp12546][2]
                                        * u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                           [k_darts_counter_temp12546][2]
                                    + u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                       [k_darts_counter_temp12546][3]
                                        * u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                           [k_darts_counter_temp12546][3])
                                / u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                   [k_darts_counter_temp12546][0];
                            flux[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                [k_darts_counter_temp12546][1]
                                = u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                   [k_darts_counter_temp12546][1]
                                * (*(*u41));
                            flux[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                [k_darts_counter_temp12546][2]
                                = u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                   [k_darts_counter_temp12546][2]
                                * (*(*u41));
                            flux[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                [k_darts_counter_temp12546][3]
                                = u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                   [k_darts_counter_temp12546][3]
                                    * (*(*u41))
                                + 0.40000000000000002
                                    * (u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                        [k_darts_counter_temp12546][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                [k_darts_counter_temp12546][4]
                                = (1.3999999999999999
                                          * u[(i_darts_counter_temp12546)]
                                             [j_darts_counter_temp12546][k_darts_counter_temp12546]
                                             [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u41));
                        }
                        (*k) = k_darts_counter_temp12546;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp12546 = (*k);
                        for (; k_darts_counter_temp12546 <= nz - 2; k_darts_counter_temp12546++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp12546 = (*m);
                                for (; m_darts_counter_temp12546 < 5; m_darts_counter_temp12546++) {
                                    rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                       [k_darts_counter_temp12546][m_darts_counter_temp12546]
                                        = rsd[(i_darts_counter_temp12546)]
                                             [j_darts_counter_temp12546][k_darts_counter_temp12546]
                                             [m_darts_counter_temp12546]
                                        - tz2
                                            * (flux[(i_darts_counter_temp12546)]
                                                   [j_darts_counter_temp12546]
                                                   [k_darts_counter_temp12546 + 1]
                                                   [m_darts_counter_temp12546]
                                                - flux[k_darts_counter_temp12546 - 1]
                                                      [(i_darts_counter_temp12546)]
                                                      [j_darts_counter_temp12546]
                                                      [m_darts_counter_temp12546]);
                                }
                                (*m) = m_darts_counter_temp12546;
                            }
                        }
                        (*k) = k_darts_counter_temp12546;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp12546 = (*k);
                        for (; k_darts_counter_temp12546 <= nz - 1; k_darts_counter_temp12546++) {
                            (*(*tmp)) = 1.
                                / u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                   [k_darts_counter_temp12546][0];
                            (*(*u21k)) = (*(*tmp))
                                * u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                   [k_darts_counter_temp12546][1];
                            (*(*u31k)) = (*(*tmp))
                                * u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                   [k_darts_counter_temp12546][2];
                            (*(*u41k)) = (*(*tmp))
                                * u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                   [k_darts_counter_temp12546][3];
                            (*(*u51k)) = (*(*tmp))
                                * u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                   [k_darts_counter_temp12546][4];
                            (*(*tmp)) = 1.
                                / u[k_darts_counter_temp12546 - 1][(i_darts_counter_temp12546)]
                                   [j_darts_counter_temp12546][0];
                            (*(*u21km1)) = (*(*tmp))
                                * u[k_darts_counter_temp12546 - 1][(i_darts_counter_temp12546)]
                                   [j_darts_counter_temp12546][1];
                            (*(*u31km1)) = (*(*tmp))
                                * u[k_darts_counter_temp12546 - 1][(i_darts_counter_temp12546)]
                                   [j_darts_counter_temp12546][2];
                            (*(*u41km1)) = (*(*tmp))
                                * u[k_darts_counter_temp12546 - 1][(i_darts_counter_temp12546)]
                                   [j_darts_counter_temp12546][3];
                            (*(*u51km1)) = (*(*tmp))
                                * u[k_darts_counter_temp12546 - 1][(i_darts_counter_temp12546)]
                                   [j_darts_counter_temp12546][4];
                            flux[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                [k_darts_counter_temp12546][1]
                                = tz3 * ((*(*u21k)) - (*(*u21km1)));
                            flux[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                [k_darts_counter_temp12546][2]
                                = tz3 * ((*(*u31k)) - (*(*u31km1)));
                            flux[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                [k_darts_counter_temp12546][3]
                                = (4. / 3.) * tz3 * ((*(*u41k)) - (*(*u41km1)));
                            flux[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                [k_darts_counter_temp12546][4]
                                = 0.5 * (1. - 1.3999999999999999 * 1.3999999999999999) * tz3
                                    * (((((*(*u21k))) * ((*(*u21k))))
                                           + (((*(*u31k))) * ((*(*u31k))))
                                           + (((*(*u41k))) * ((*(*u41k)))))
                                        - ((((*(*u21km1))) * ((*(*u21km1))))
                                            + (((*(*u31km1))) * ((*(*u31km1))))
                                            + (((*(*u41km1))) * ((*(*u41km1))))))
                                + (1. / 6.) * tz3
                                    * ((((*(*u41k))) * ((*(*u41k))))
                                        - (((*(*u41km1))) * ((*(*u41km1)))))
                                + 1.3999999999999999 * 1.3999999999999999 * tz3
                                    * ((*(*u51k)) - (*(*u51km1)));
                        }
                        (*k) = k_darts_counter_temp12546;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp12546 = (*k);
                        for (; k_darts_counter_temp12546 <= nz - 2; k_darts_counter_temp12546++) {
                            rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                               [k_darts_counter_temp12546][0]
                                = rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                     [k_darts_counter_temp12546][0]
                                + dz1 * tz1
                                    * (u[k_darts_counter_temp12546 - 1][(i_darts_counter_temp12546)]
                                        [j_darts_counter_temp12546][0]
                                        - 2.
                                            * u[(i_darts_counter_temp12546)]
                                               [j_darts_counter_temp12546]
                                               [k_darts_counter_temp12546][0]
                                        + u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                           [k_darts_counter_temp12546 + 1][0]);
                            rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                               [k_darts_counter_temp12546][1]
                                = rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                     [k_darts_counter_temp12546][1]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                           [k_darts_counter_temp12546 + 1][1]
                                        - flux[(i_darts_counter_temp12546)]
                                              [j_darts_counter_temp12546][k_darts_counter_temp12546]
                                              [1])
                                + dz2 * tz1
                                    * (u[k_darts_counter_temp12546 - 1][(i_darts_counter_temp12546)]
                                        [j_darts_counter_temp12546][1]
                                        - 2.
                                            * u[(i_darts_counter_temp12546)]
                                               [j_darts_counter_temp12546]
                                               [k_darts_counter_temp12546][1]
                                        + u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                           [k_darts_counter_temp12546 + 1][1]);
                            rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                               [k_darts_counter_temp12546][2]
                                = rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                     [k_darts_counter_temp12546][2]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                           [k_darts_counter_temp12546 + 1][2]
                                        - flux[(i_darts_counter_temp12546)]
                                              [j_darts_counter_temp12546][k_darts_counter_temp12546]
                                              [2])
                                + dz3 * tz1
                                    * (u[k_darts_counter_temp12546 - 1][(i_darts_counter_temp12546)]
                                        [j_darts_counter_temp12546][2]
                                        - 2.
                                            * u[(i_darts_counter_temp12546)]
                                               [j_darts_counter_temp12546]
                                               [k_darts_counter_temp12546][2]
                                        + u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                           [k_darts_counter_temp12546 + 1][2]);
                            rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                               [k_darts_counter_temp12546][3]
                                = rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                     [k_darts_counter_temp12546][3]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                           [k_darts_counter_temp12546 + 1][3]
                                        - flux[(i_darts_counter_temp12546)]
                                              [j_darts_counter_temp12546][k_darts_counter_temp12546]
                                              [3])
                                + dz4 * tz1
                                    * (u[k_darts_counter_temp12546 - 1][(i_darts_counter_temp12546)]
                                        [j_darts_counter_temp12546][3]
                                        - 2.
                                            * u[(i_darts_counter_temp12546)]
                                               [j_darts_counter_temp12546]
                                               [k_darts_counter_temp12546][3]
                                        + u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                           [k_darts_counter_temp12546 + 1][3]);
                            rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                               [k_darts_counter_temp12546][4]
                                = rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                     [k_darts_counter_temp12546][4]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                           [k_darts_counter_temp12546 + 1][4]
                                        - flux[(i_darts_counter_temp12546)]
                                              [j_darts_counter_temp12546][k_darts_counter_temp12546]
                                              [4])
                                + dz5 * tz1
                                    * (u[k_darts_counter_temp12546 - 1][(i_darts_counter_temp12546)]
                                        [j_darts_counter_temp12546][4]
                                        - 2.
                                            * u[(i_darts_counter_temp12546)]
                                               [j_darts_counter_temp12546]
                                               [k_darts_counter_temp12546][4]
                                        + u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                           [k_darts_counter_temp12546 + 1][4]);
                        }
                        (*k) = k_darts_counter_temp12546;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp12546 = (*m);
                        for (; m_darts_counter_temp12546 < 5; m_darts_counter_temp12546++) {
                            rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546][1]
                               [m_darts_counter_temp12546]
                                = rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546][1]
                                     [m_darts_counter_temp12546]
                                - dssp
                                    * (+5.
                                            * u[(i_darts_counter_temp12546)]
                                               [j_darts_counter_temp12546][1]
                                               [m_darts_counter_temp12546]
                                        - 4.
                                            * u[(i_darts_counter_temp12546)]
                                               [j_darts_counter_temp12546][2]
                                               [m_darts_counter_temp12546]
                                        + u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                           [3][m_darts_counter_temp12546]);
                            rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546][2]
                               [m_darts_counter_temp12546]
                                = rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546][2]
                                     [m_darts_counter_temp12546]
                                - dssp
                                    * (-4.
                                            * u[(i_darts_counter_temp12546)]
                                               [j_darts_counter_temp12546][1]
                                               [m_darts_counter_temp12546]
                                        + 6.
                                            * u[(i_darts_counter_temp12546)]
                                               [j_darts_counter_temp12546][2]
                                               [m_darts_counter_temp12546]
                                        - 4.
                                            * u[(i_darts_counter_temp12546)]
                                               [j_darts_counter_temp12546][3]
                                               [m_darts_counter_temp12546]
                                        + u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                           [4][m_darts_counter_temp12546]);
                        }
                        (*m) = m_darts_counter_temp12546;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 3;
                        int k_darts_counter_temp12546 = (*k);
                        for (; k_darts_counter_temp12546 <= nz - 4; k_darts_counter_temp12546++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp12546 = (*m);
                                for (; m_darts_counter_temp12546 < 5; m_darts_counter_temp12546++) {
                                    rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                       [k_darts_counter_temp12546][m_darts_counter_temp12546]
                                        = rsd[(i_darts_counter_temp12546)]
                                             [j_darts_counter_temp12546][k_darts_counter_temp12546]
                                             [m_darts_counter_temp12546]
                                        - dssp
                                            * (u[(i_darts_counter_temp12546)]
                                                [j_darts_counter_temp12546]
                                                [k_darts_counter_temp12546 - 2]
                                                [m_darts_counter_temp12546]
                                                - 4.
                                                    * u[k_darts_counter_temp12546 - 1]
                                                       [(i_darts_counter_temp12546)]
                                                       [j_darts_counter_temp12546]
                                                       [m_darts_counter_temp12546]
                                                + 6.
                                                    * u[(i_darts_counter_temp12546)]
                                                       [j_darts_counter_temp12546]
                                                       [k_darts_counter_temp12546]
                                                       [m_darts_counter_temp12546]
                                                - 4.
                                                    * u[(i_darts_counter_temp12546)]
                                                       [j_darts_counter_temp12546]
                                                       [k_darts_counter_temp12546 + 1]
                                                       [m_darts_counter_temp12546]
                                                + u[(i_darts_counter_temp12546)]
                                                   [j_darts_counter_temp12546]
                                                   [k_darts_counter_temp12546 + 2]
                                                   [m_darts_counter_temp12546]);
                                }
                                (*m) = m_darts_counter_temp12546;
                            }
                        }
                        (*k) = k_darts_counter_temp12546;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp12546 = (*m);
                        for (; m_darts_counter_temp12546 < 5; m_darts_counter_temp12546++) {
                            rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546][nz - 3]
                               [m_darts_counter_temp12546]
                                = rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                     [nz - 3][m_darts_counter_temp12546]
                                - dssp
                                    * (u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                        [nz - 5][m_darts_counter_temp12546]
                                        - 4.
                                            * u[(i_darts_counter_temp12546)]
                                               [j_darts_counter_temp12546][nz - 4]
                                               [m_darts_counter_temp12546]
                                        + 6.
                                            * u[(i_darts_counter_temp12546)]
                                               [j_darts_counter_temp12546][nz - 3]
                                               [m_darts_counter_temp12546]
                                        - 4.
                                            * u[(i_darts_counter_temp12546)]
                                               [j_darts_counter_temp12546][nz - 2]
                                               [m_darts_counter_temp12546]);
                            rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546][nz - 2]
                               [m_darts_counter_temp12546]
                                = rsd[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                     [nz - 2][m_darts_counter_temp12546]
                                - dssp
                                    * (u[(i_darts_counter_temp12546)][j_darts_counter_temp12546]
                                        [nz - 4][m_darts_counter_temp12546]
                                        - 4.
                                            * u[(i_darts_counter_temp12546)]
                                               [j_darts_counter_temp12546][nz - 3]
                                               [m_darts_counter_temp12546]
                                        + 5.
                                            * u[(i_darts_counter_temp12546)]
                                               [j_darts_counter_temp12546][nz - 2]
                                               [m_darts_counter_temp12546]);
                        }
                        (*m) = m_darts_counter_temp12546;
                    }
                }
                (*j) = j_darts_counter_temp12546;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets12546[0].decDep();
}
TP12546::TP12546(int in_numThreads, int in_mainCodeletID, TP10896* in_TPParent,
    int in_initIteration, int in_lastIteration, TP12546** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts12546(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts12546(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts12546(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts12546(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts12546(new double*[this->numThreads])
    , tmp_darts12546(new double*[this->numThreads])
    , u21k_darts12546(new double*[this->numThreads])
    , u21km1_darts12546(new double*[this->numThreads])
    , u31k_darts12546(new double*[this->numThreads])
    , u31km1_darts12546(new double*[this->numThreads])
    , u41_darts12546(new double*[this->numThreads])
    , u41k_darts12546(new double*[this->numThreads])
    , u41km1_darts12546(new double*[this->numThreads])
    , u51k_darts12546(new double*[this->numThreads])
    , u51km1_darts12546(new double*[this->numThreads])
    , initIteration12546(in_initIteration)
    , lastIteration12546(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets12546(new _barrierCodelets12546[1])
    , checkInCodelets12547(new _checkInCodelets12547[this->numThreads])
{
    /*Initialize the loop parameters*/
    range12546 = abs(lastIteration12546 - initIteration12546) / 1;
    rangePerCodelet12546 = range12546 / numThreads;
    minIteration12546 = min<int>(lastIteration12546, initIteration12546);
    remainderRange12546 = range12546 % numThreads;
    /*Initialize inputs and vars.*/
    this->q_darts12546
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts12546
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21k_darts12546
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21km1_darts12546
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31k_darts12546
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31km1_darts12546
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41_darts12546
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41k_darts12546
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41km1_darts12546
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51k_darts12546
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51km1_darts12546
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets12546[0] = _barrierCodelets12546(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets12547* checkInCodelets12547Ptr = (this->checkInCodelets12547);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets12547);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets12547Ptr) = _checkInCodelets12547(2, 1, this, codeletCounter);
#else
        (*checkInCodelets12547Ptr) = _checkInCodelets12547(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets12547Ptr).decDep();
        checkInCodelets12547Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP12546::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets12547[localID].setID(codeletID);
    this->checkInCodelets12547[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets12547[localID + this->baseNumThreads * i]
            = _checkInCodelets12547(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets12547[localID + this->baseNumThreads * i]
            = _checkInCodelets12547(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets12547[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets12547[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP12546::~TP12546()
{
    delete[] q_darts12546;
    delete[] tmp_darts12546;
    delete[] u21k_darts12546;
    delete[] u21km1_darts12546;
    delete[] u31k_darts12546;
    delete[] u31km1_darts12546;
    delete[] u41_darts12546;
    delete[] u41k_darts12546;
    delete[] u41km1_darts12546;
    delete[] u51k_darts12546;
    delete[] u51km1_darts12546;
    delete[] barrierCodelets12546;
    delete[] checkInCodelets12547;
}
/*TP13297: OMPParallelDirective*/
void TP13297::_barrierCodelets13297::fire(void)
{
    TP13297* myTP = static_cast<TP13297*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP13297::_checkInCodelets13301::fire(void)
{
    /*region 13301 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13301;
    if (idx < myTP->TPsToUse13301) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13301_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(nx - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13301;
            int minIteration = min<int>(nx, 0);
            int remainderRange = range % myTP->TPsToUse13301;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (0 < nx) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse13301 - 1) {
                lastIteration = nx;
            }
#if USEINVOKE == 1
            invoke<TP13301>(myTP, myTP->codeletsPerTP13301 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13301Ptr[idx]));
#else
            place<TP13301>(idx, myTP, myTP->codeletsPerTP13301 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13301Ptr[idx]));
#endif
        } else {
            if (myTP->TP13301Ptr[idx] != nullptr) {
                myTP->TP13301Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13297::_barrierCodelets13301::fire(void)
{
    TP13297* myTP = static_cast<TP13297*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets13352[codeletsCounter].decDep();
        }
    }
}
void TP13297::_checkInCodelets13352::fire(void)
{
    /*region 13352 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13352;
    if (idx < myTP->TPsToUse13352) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13352_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(nx - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13352;
            int minIteration = min<int>(nx, 0);
            int remainderRange = range % myTP->TPsToUse13352;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (0 < nx) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse13352 - 1) {
                lastIteration = nx;
            }
#if USEINVOKE == 1
            invoke<TP13352>(myTP, myTP->codeletsPerTP13352 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13352Ptr[idx]));
#else
            place<TP13352>(idx, myTP, myTP->codeletsPerTP13352 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13352Ptr[idx]));
#endif
        } else {
            if (myTP->TP13352Ptr[idx] != nullptr) {
                myTP->TP13352Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13297::_barrierCodelets13352::fire(void)
{
    TP13297* myTP = static_cast<TP13297*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets13394[codeletsCounter].decDep();
        }
    }
}
void TP13297::_checkInCodelets13394::fire(void)
{
    /*region 13394 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13394;
    if (idx < myTP->TPsToUse13394) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13394_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(nx - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13394;
            int minIteration = min<int>(nx, 0);
            int remainderRange = range % myTP->TPsToUse13394;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (0 < nx) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse13394 - 1) {
                lastIteration = nx;
            }
#if USEINVOKE == 1
            invoke<TP13394>(myTP, myTP->codeletsPerTP13394 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13394Ptr[idx]));
#else
            place<TP13394>(idx, myTP, myTP->codeletsPerTP13394 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13394Ptr[idx]));
#endif
        } else {
            if (myTP->TP13394Ptr[idx] != nullptr) {
                myTP->TP13394Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13297::_barrierCodelets13394::fire(void)
{
    TP13297* myTP = static_cast<TP13297*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets13438[codeletsCounter].decDep();
        }
    }
}
void TP13297::_checkInCodelets13438::fire(void)
{
    /*region 13438 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13438;
    if (idx < myTP->TPsToUse13438) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13438_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(ny - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13438;
            int minIteration = min<int>(ny, 0);
            int remainderRange = range % myTP->TPsToUse13438;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (0 < ny) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse13438 - 1) {
                lastIteration = ny;
            }
#if USEINVOKE == 1
            invoke<TP13438>(myTP, myTP->codeletsPerTP13438 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13438Ptr[idx]));
#else
            place<TP13438>(idx, myTP, myTP->codeletsPerTP13438 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13438Ptr[idx]));
#endif
        } else {
            if (myTP->TP13438Ptr[idx] != nullptr) {
                myTP->TP13438Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13297::_barrierCodelets13438::fire(void)
{
    TP13297* myTP = static_cast<TP13297*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets13480[codeletsCounter].decDep();
        }
    }
}
void TP13297::_checkInCodelets13480::fire(void)
{
    /*region 13480 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13480;
    if (idx < myTP->TPsToUse13480) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13480_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(ny - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13480;
            int minIteration = min<int>(ny, 0);
            int remainderRange = range % myTP->TPsToUse13480;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (0 < ny) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse13480 - 1) {
                lastIteration = ny;
            }
#if USEINVOKE == 1
            invoke<TP13480>(myTP, myTP->codeletsPerTP13480 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13480Ptr[idx]));
#else
            place<TP13480>(idx, myTP, myTP->codeletsPerTP13480 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13480Ptr[idx]));
#endif
        } else {
            if (myTP->TP13480Ptr[idx] != nullptr) {
                myTP->TP13480Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13297::_barrierCodelets13480::fire(void)
{
    TP13297* myTP = static_cast<TP13297*>(myTP_);
    myTP->TPParent->barrierCodelets13297[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets13297[0]));
}
TP13297::TP13297(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13297(new int[this->numThreads]) /*VARIABLE*/
    , iglob_darts13297(new int[this->numThreads]) /*VARIABLE*/
    , j_darts13297(new int[this->numThreads]) /*VARIABLE*/
    , jglob_darts13297(new int[this->numThreads]) /*VARIABLE*/
    , k_darts13297(new int[this->numThreads]) /*VARIABLE*/
    , TP13301Ptr(new TP13301*[NUMTPS13301])
    , TP13301_alreadyLaunched(new size_t[NUMTPS13301])
    , numTPsSet13301(0)
    , numTPsReady13301(0)
    , TPsToUse13301(NUMTPS13301)
    , codeletsPerTP13301(this->numThreads / NUMTPS13301)
    , totalCodelets13301(this->TPsToUse13301 * this->codeletsPerTP13301)
    , TP13352Ptr(new TP13352*[NUMTPS13352])
    , TP13352_alreadyLaunched(new size_t[NUMTPS13352])
    , numTPsSet13352(0)
    , numTPsReady13352(0)
    , TPsToUse13352(NUMTPS13352)
    , codeletsPerTP13352(this->numThreads / NUMTPS13352)
    , totalCodelets13352(this->TPsToUse13352 * this->codeletsPerTP13352)
    , TP13394Ptr(new TP13394*[NUMTPS13394])
    , TP13394_alreadyLaunched(new size_t[NUMTPS13394])
    , numTPsSet13394(0)
    , numTPsReady13394(0)
    , TPsToUse13394(NUMTPS13394)
    , codeletsPerTP13394(this->numThreads / NUMTPS13394)
    , totalCodelets13394(this->TPsToUse13394 * this->codeletsPerTP13394)
    , TP13438Ptr(new TP13438*[NUMTPS13438])
    , TP13438_alreadyLaunched(new size_t[NUMTPS13438])
    , numTPsSet13438(0)
    , numTPsReady13438(0)
    , TPsToUse13438(NUMTPS13438)
    , codeletsPerTP13438(this->numThreads / NUMTPS13438)
    , totalCodelets13438(this->TPsToUse13438 * this->codeletsPerTP13438)
    , TP13480Ptr(new TP13480*[NUMTPS13480])
    , TP13480_alreadyLaunched(new size_t[NUMTPS13480])
    , numTPsSet13480(0)
    , numTPsReady13480(0)
    , TPsToUse13480(NUMTPS13480)
    , codeletsPerTP13480(this->numThreads / NUMTPS13480)
    , totalCodelets13480(this->TPsToUse13480 * this->codeletsPerTP13480)
    , barrierCodelets13297(new _barrierCodelets13297[1])
    , checkInCodelets13301(new _checkInCodelets13301[this->numThreads])
    , barrierCodelets13301(new _barrierCodelets13301[1])
    , checkInCodelets13352(new _checkInCodelets13352[this->numThreads])
    , barrierCodelets13352(new _barrierCodelets13352[1])
    , checkInCodelets13394(new _checkInCodelets13394[this->numThreads])
    , barrierCodelets13394(new _barrierCodelets13394[1])
    , checkInCodelets13438(new _checkInCodelets13438[this->numThreads])
    , barrierCodelets13438(new _barrierCodelets13438[1])
    , checkInCodelets13480(new _checkInCodelets13480[this->numThreads])
    , barrierCodelets13480(new _barrierCodelets13480[1])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets13297[0] = _barrierCodelets13297(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets13480[0] = _barrierCodelets13480(NUMTPS13480, NUMTPS13480, this, 0);
    barrierCodelets13438[0] = _barrierCodelets13438(NUMTPS13438, NUMTPS13438, this, 0);
    barrierCodelets13394[0] = _barrierCodelets13394(NUMTPS13394, NUMTPS13394, this, 0);
    barrierCodelets13352[0] = _barrierCodelets13352(NUMTPS13352, NUMTPS13352, this, 0);
    barrierCodelets13301[0] = _barrierCodelets13301(NUMTPS13301, NUMTPS13301, this, 0);
    _checkInCodelets13480* checkInCodelets13480Ptr = (this->checkInCodelets13480);
    for (int i = 0; i < NUMTPS13480; i++) {
        TP13480Ptr[i] = nullptr;
        TP13480_alreadyLaunched[i] = 0;
    }
    _checkInCodelets13438* checkInCodelets13438Ptr = (this->checkInCodelets13438);
    for (int i = 0; i < NUMTPS13438; i++) {
        TP13438Ptr[i] = nullptr;
        TP13438_alreadyLaunched[i] = 0;
    }
    _checkInCodelets13394* checkInCodelets13394Ptr = (this->checkInCodelets13394);
    for (int i = 0; i < NUMTPS13394; i++) {
        TP13394Ptr[i] = nullptr;
        TP13394_alreadyLaunched[i] = 0;
    }
    _checkInCodelets13352* checkInCodelets13352Ptr = (this->checkInCodelets13352);
    for (int i = 0; i < NUMTPS13352; i++) {
        TP13352Ptr[i] = nullptr;
        TP13352_alreadyLaunched[i] = 0;
    }
    _checkInCodelets13301* checkInCodelets13301Ptr = (this->checkInCodelets13301);
    for (int i = 0; i < NUMTPS13301; i++) {
        TP13301Ptr[i] = nullptr;
        TP13301_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets13480Ptr) = _checkInCodelets13480(1, 1, this, codeletCounter);
        checkInCodelets13480Ptr++;
        (*checkInCodelets13438Ptr) = _checkInCodelets13438(1, 1, this, codeletCounter);
        checkInCodelets13438Ptr++;
        (*checkInCodelets13394Ptr) = _checkInCodelets13394(1, 1, this, codeletCounter);
        checkInCodelets13394Ptr++;
        (*checkInCodelets13352Ptr) = _checkInCodelets13352(1, 1, this, codeletCounter);
        checkInCodelets13352Ptr++;
        (*checkInCodelets13301Ptr) = _checkInCodelets13301(1, 1, this, codeletCounter);
        (*checkInCodelets13301Ptr).decDep();
        checkInCodelets13301Ptr++;
    }
}
TP13297::~TP13297()
{
    delete[] i_darts13297;
    delete[] iglob_darts13297;
    delete[] j_darts13297;
    delete[] jglob_darts13297;
    delete[] k_darts13297;
    delete[] barrierCodelets13297;
    delete[] barrierCodelets13480;
    delete[] checkInCodelets13480;
    delete[] barrierCodelets13438;
    delete[] checkInCodelets13438;
    delete[] barrierCodelets13394;
    delete[] checkInCodelets13394;
    delete[] barrierCodelets13352;
    delete[] checkInCodelets13352;
    delete[] barrierCodelets13301;
    delete[] checkInCodelets13301;
}
/*TP13301: OMPForDirective*/
void TP13301::_barrierCodelets13301::fire(void)
{
    TP13301* myTP = static_cast<TP13301*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13301[0].decDep();
}
bool TP13301::requestNewRangeIterations13301(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13301 * codeletID;
        int tempEndRange = rangePerCodelet13301 * (codeletID + 1);
        if (remainderRange13301 != 0) {
            if (codeletID < (uint32_t)remainderRange13301) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13301;
                tempEndRange += remainderRange13301;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13301;
        tempEndRange = tempEndRange * 1 + minIteration13301;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13301 < lastIteration13301) {
            (this->inputsTPParent->i_darts13301[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts13301[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13301;
        }
    }
    return isThereNewIteration;
}
void TP13301::_checkInCodelets13302::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->iglob_darts13301[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts13297[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jglob_darts13301[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts13297[this->getID()]);

    /*printing node 13302: ForStmt*/
    /*var: i*/
    /*var: iglob*/
    /*var: j*/
    /*var: jglob*/
    int* i = &(this->inputsTPParent->i_darts13301[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iglob = &(this->inputsTPParent->iglob_darts13301[this->getLocalID()]);
    (void)iglob /*OMP_SHARED_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts13301[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jglob = &(this->inputsTPParent->jglob_darts13301[this->getLocalID()]);
    (void)jglob /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13301(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13301[0].decDep();
        return;
    }
    for (int i_darts_counter_temp13301 = (*i); i_darts_counter_temp13301 < endRange
         && i_darts_counter_temp13301 < this->inputsTPParent->lastIteration13301;
         i_darts_counter_temp13301++) {
        {
            (*(*iglob)) = (i_darts_counter_temp13301);
            {
                /*Loop's init*/
                (*j) = 0;
                int j_darts_counter_temp13301 = (*j);
                for (; j_darts_counter_temp13301 < ny; j_darts_counter_temp13301++) {
                    (*(*jglob)) = j_darts_counter_temp13301;
                    exact((*(*iglob)), (*(*jglob)), 0,
                        &u[(i_darts_counter_temp13301)][j_darts_counter_temp13301][0][0]);
                    exact((*(*iglob)), (*(*jglob)), nz - 1,
                        &u[(i_darts_counter_temp13301)][j_darts_counter_temp13301][nz - 1][0]);
                }
                (*j) = j_darts_counter_temp13301;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13301[0].decDep();
}
TP13301::TP13301(int in_numThreads, int in_mainCodeletID, TP13297* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13301** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13301(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iglob_darts13301(new int*[this->numThreads])
    , j_darts13301(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jglob_darts13301(new int*[this->numThreads])
    , initIteration13301(in_initIteration)
    , lastIteration13301(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13301(new _barrierCodelets13301[1])
    , checkInCodelets13302(new _checkInCodelets13302[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13301 = abs(lastIteration13301 - initIteration13301) / 1;
    rangePerCodelet13301 = range13301 / numThreads;
    minIteration13301 = min<int>(lastIteration13301, initIteration13301);
    remainderRange13301 = range13301 % numThreads;
    /*Initialize inputs and vars.*/
    this->iglob_darts13301 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jglob_darts13301 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets13301[0] = _barrierCodelets13301(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13302* checkInCodelets13302Ptr = (this->checkInCodelets13302);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13302);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13302Ptr) = _checkInCodelets13302(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13302Ptr) = _checkInCodelets13302(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13302Ptr).decDep();
        checkInCodelets13302Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13301::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13302[localID].setID(codeletID);
    this->checkInCodelets13302[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13302[localID + this->baseNumThreads * i]
            = _checkInCodelets13302(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13302[localID + this->baseNumThreads * i]
            = _checkInCodelets13302(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13302[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13302[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP13301::~TP13301()
{
    delete[] iglob_darts13301;
    delete[] jglob_darts13301;
    delete[] barrierCodelets13301;
    delete[] checkInCodelets13302;
}
/*TP13352: OMPForDirective*/
void TP13352::_barrierCodelets13352::fire(void)
{
    TP13352* myTP = static_cast<TP13352*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13352[0].decDep();
}
bool TP13352::requestNewRangeIterations13352(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13352 * codeletID;
        int tempEndRange = rangePerCodelet13352 * (codeletID + 1);
        if (remainderRange13352 != 0) {
            if (codeletID < (uint32_t)remainderRange13352) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13352;
                tempEndRange += remainderRange13352;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13352;
        tempEndRange = tempEndRange * 1 + minIteration13352;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13352 < lastIteration13352) {
            (this->inputsTPParent->i_darts13352[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts13352[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13352;
        }
    }
    return isThereNewIteration;
}
void TP13352::_checkInCodelets13353::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->iglob_darts13352[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts13297[this->getID()]);

    /*printing node 13353: ForStmt*/
    /*var: i*/
    /*var: iglob*/
    /*var: k*/
    int* i = &(this->inputsTPParent->i_darts13352[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iglob = &(this->inputsTPParent->iglob_darts13352[this->getLocalID()]);
    (void)iglob /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts13352[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13352(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13352[0].decDep();
        return;
    }
    for (int i_darts_counter_temp13352 = (*i); i_darts_counter_temp13352 < endRange
         && i_darts_counter_temp13352 < this->inputsTPParent->lastIteration13352;
         i_darts_counter_temp13352++) {
        {
            (*(*iglob)) = (i_darts_counter_temp13352);
            {
                /*Loop's init*/
                (*k) = 0;
                int k_darts_counter_temp13352 = (*k);
                for (; k_darts_counter_temp13352 < nz; k_darts_counter_temp13352++) {
                    exact((*(*iglob)), 0, k_darts_counter_temp13352,
                        &u[(i_darts_counter_temp13352)][0][k_darts_counter_temp13352][0]);
                }
                (*k) = k_darts_counter_temp13352;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13352[0].decDep();
}
TP13352::TP13352(int in_numThreads, int in_mainCodeletID, TP13297* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13352** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13352(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iglob_darts13352(new int*[this->numThreads])
    , k_darts13352(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration13352(in_initIteration)
    , lastIteration13352(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13352(new _barrierCodelets13352[1])
    , checkInCodelets13353(new _checkInCodelets13353[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13352 = abs(lastIteration13352 - initIteration13352) / 1;
    rangePerCodelet13352 = range13352 / numThreads;
    minIteration13352 = min<int>(lastIteration13352, initIteration13352);
    remainderRange13352 = range13352 % numThreads;
    /*Initialize inputs and vars.*/
    this->iglob_darts13352 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets13352[0] = _barrierCodelets13352(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13353* checkInCodelets13353Ptr = (this->checkInCodelets13353);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13353);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13353Ptr) = _checkInCodelets13353(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13353Ptr) = _checkInCodelets13353(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13353Ptr).decDep();
        checkInCodelets13353Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13352::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13353[localID].setID(codeletID);
    this->checkInCodelets13353[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13353[localID + this->baseNumThreads * i]
            = _checkInCodelets13353(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13353[localID + this->baseNumThreads * i]
            = _checkInCodelets13353(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13353[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13353[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP13352::~TP13352()
{
    delete[] iglob_darts13352;
    delete[] barrierCodelets13352;
    delete[] checkInCodelets13353;
}
/*TP13394: OMPForDirective*/
void TP13394::_barrierCodelets13394::fire(void)
{
    TP13394* myTP = static_cast<TP13394*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13394[0].decDep();
}
bool TP13394::requestNewRangeIterations13394(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13394 * codeletID;
        int tempEndRange = rangePerCodelet13394 * (codeletID + 1);
        if (remainderRange13394 != 0) {
            if (codeletID < (uint32_t)remainderRange13394) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13394;
                tempEndRange += remainderRange13394;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13394;
        tempEndRange = tempEndRange * 1 + minIteration13394;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13394 < lastIteration13394) {
            (this->inputsTPParent->i_darts13394[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts13394[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13394;
        }
    }
    return isThereNewIteration;
}
void TP13394::_checkInCodelets13395::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->iglob_darts13394[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts13297[this->getID()]);

    /*printing node 13395: ForStmt*/
    /*var: i*/
    /*var: iglob*/
    /*var: k*/
    int* i = &(this->inputsTPParent->i_darts13394[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iglob = &(this->inputsTPParent->iglob_darts13394[this->getLocalID()]);
    (void)iglob /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts13394[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13394(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13394[0].decDep();
        return;
    }
    for (int i_darts_counter_temp13394 = (*i); i_darts_counter_temp13394 < endRange
         && i_darts_counter_temp13394 < this->inputsTPParent->lastIteration13394;
         i_darts_counter_temp13394++) {
        {
            (*(*iglob)) = (i_darts_counter_temp13394);
            {
                /*Loop's init*/
                (*k) = 0;
                int k_darts_counter_temp13394 = (*k);
                for (; k_darts_counter_temp13394 < nz; k_darts_counter_temp13394++) {
                    exact((*(*iglob)), ny0 - 1, k_darts_counter_temp13394,
                        &u[(i_darts_counter_temp13394)][ny - 1][k_darts_counter_temp13394][0]);
                }
                (*k) = k_darts_counter_temp13394;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13394[0].decDep();
}
TP13394::TP13394(int in_numThreads, int in_mainCodeletID, TP13297* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13394** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13394(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iglob_darts13394(new int*[this->numThreads])
    , k_darts13394(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration13394(in_initIteration)
    , lastIteration13394(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13394(new _barrierCodelets13394[1])
    , checkInCodelets13395(new _checkInCodelets13395[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13394 = abs(lastIteration13394 - initIteration13394) / 1;
    rangePerCodelet13394 = range13394 / numThreads;
    minIteration13394 = min<int>(lastIteration13394, initIteration13394);
    remainderRange13394 = range13394 % numThreads;
    /*Initialize inputs and vars.*/
    this->iglob_darts13394 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets13394[0] = _barrierCodelets13394(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13395* checkInCodelets13395Ptr = (this->checkInCodelets13395);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13395);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13395Ptr) = _checkInCodelets13395(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13395Ptr) = _checkInCodelets13395(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13395Ptr).decDep();
        checkInCodelets13395Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13394::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13395[localID].setID(codeletID);
    this->checkInCodelets13395[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13395[localID + this->baseNumThreads * i]
            = _checkInCodelets13395(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13395[localID + this->baseNumThreads * i]
            = _checkInCodelets13395(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13395[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13395[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP13394::~TP13394()
{
    delete[] iglob_darts13394;
    delete[] barrierCodelets13394;
    delete[] checkInCodelets13395;
}
/*TP13438: OMPForDirective*/
void TP13438::_barrierCodelets13438::fire(void)
{
    TP13438* myTP = static_cast<TP13438*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13438[0].decDep();
}
bool TP13438::requestNewRangeIterations13438(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13438 * codeletID;
        int tempEndRange = rangePerCodelet13438 * (codeletID + 1);
        if (remainderRange13438 != 0) {
            if (codeletID < (uint32_t)remainderRange13438) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13438;
                tempEndRange += remainderRange13438;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13438;
        tempEndRange = tempEndRange * 1 + minIteration13438;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13438 < lastIteration13438) {
            (this->inputsTPParent->j_darts13438[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->j_darts13438[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13438;
        }
    }
    return isThereNewIteration;
}
void TP13438::_checkInCodelets13439::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jglob_darts13438[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts13297[this->getID()]);

    /*printing node 13439: ForStmt*/
    /*var: j*/
    /*var: jglob*/
    /*var: k*/
    int* j = &(this->inputsTPParent->j_darts13438[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jglob = &(this->inputsTPParent->jglob_darts13438[this->getLocalID()]);
    (void)jglob /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts13438[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13438(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13438[0].decDep();
        return;
    }
    for (int j_darts_counter_temp13438 = (*j); j_darts_counter_temp13438 < endRange
         && j_darts_counter_temp13438 < this->inputsTPParent->lastIteration13438;
         j_darts_counter_temp13438++) {
        {
            (*(*jglob)) = (j_darts_counter_temp13438);
            {
                /*Loop's init*/
                (*k) = 0;
                int k_darts_counter_temp13438 = (*k);
                for (; k_darts_counter_temp13438 < nz; k_darts_counter_temp13438++) {
                    exact(0, (*(*jglob)), k_darts_counter_temp13438,
                        &u[0][(j_darts_counter_temp13438)][k_darts_counter_temp13438][0]);
                }
                (*k) = k_darts_counter_temp13438;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13438[0].decDep();
}
TP13438::TP13438(int in_numThreads, int in_mainCodeletID, TP13297* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13438** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , j_darts13438(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jglob_darts13438(new int*[this->numThreads])
    , k_darts13438(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration13438(in_initIteration)
    , lastIteration13438(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13438(new _barrierCodelets13438[1])
    , checkInCodelets13439(new _checkInCodelets13439[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13438 = abs(lastIteration13438 - initIteration13438) / 1;
    rangePerCodelet13438 = range13438 / numThreads;
    minIteration13438 = min<int>(lastIteration13438, initIteration13438);
    remainderRange13438 = range13438 % numThreads;
    /*Initialize inputs and vars.*/
    this->jglob_darts13438 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets13438[0] = _barrierCodelets13438(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13439* checkInCodelets13439Ptr = (this->checkInCodelets13439);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13439);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13439Ptr) = _checkInCodelets13439(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13439Ptr) = _checkInCodelets13439(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13439Ptr).decDep();
        checkInCodelets13439Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13438::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13439[localID].setID(codeletID);
    this->checkInCodelets13439[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13439[localID + this->baseNumThreads * i]
            = _checkInCodelets13439(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13439[localID + this->baseNumThreads * i]
            = _checkInCodelets13439(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13439[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13439[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP13438::~TP13438()
{
    delete[] jglob_darts13438;
    delete[] barrierCodelets13438;
    delete[] checkInCodelets13439;
}
/*TP13480: OMPForDirective*/
void TP13480::_barrierCodelets13480::fire(void)
{
    TP13480* myTP = static_cast<TP13480*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13480[0].decDep();
}
bool TP13480::requestNewRangeIterations13480(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13480 * codeletID;
        int tempEndRange = rangePerCodelet13480 * (codeletID + 1);
        if (remainderRange13480 != 0) {
            if (codeletID < (uint32_t)remainderRange13480) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13480;
                tempEndRange += remainderRange13480;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13480;
        tempEndRange = tempEndRange * 1 + minIteration13480;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13480 < lastIteration13480) {
            (this->inputsTPParent->j_darts13480[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->j_darts13480[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13480;
        }
    }
    return isThereNewIteration;
}
void TP13480::_checkInCodelets13481::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jglob_darts13480[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts13297[this->getID()]);

    /*printing node 13481: ForStmt*/
    /*var: j*/
    /*var: jglob*/
    /*var: k*/
    int* j = &(this->inputsTPParent->j_darts13480[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jglob = &(this->inputsTPParent->jglob_darts13480[this->getLocalID()]);
    (void)jglob /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts13480[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13480(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13480[0].decDep();
        return;
    }
    for (int j_darts_counter_temp13480 = (*j); j_darts_counter_temp13480 < endRange
         && j_darts_counter_temp13480 < this->inputsTPParent->lastIteration13480;
         j_darts_counter_temp13480++) {
        {
            (*(*jglob)) = (j_darts_counter_temp13480);
            {
                /*Loop's init*/
                (*k) = 0;
                int k_darts_counter_temp13480 = (*k);
                for (; k_darts_counter_temp13480 < nz; k_darts_counter_temp13480++) {
                    exact(nx0 - 1, (*(*jglob)), k_darts_counter_temp13480,
                        &u[nx - 1][(j_darts_counter_temp13480)][k_darts_counter_temp13480][0]);
                }
                (*k) = k_darts_counter_temp13480;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13480[0].decDep();
}
TP13480::TP13480(int in_numThreads, int in_mainCodeletID, TP13297* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13480** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , j_darts13480(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jglob_darts13480(new int*[this->numThreads])
    , k_darts13480(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration13480(in_initIteration)
    , lastIteration13480(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13480(new _barrierCodelets13480[1])
    , checkInCodelets13481(new _checkInCodelets13481[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13480 = abs(lastIteration13480 - initIteration13480) / 1;
    rangePerCodelet13480 = range13480 / numThreads;
    minIteration13480 = min<int>(lastIteration13480, initIteration13480);
    remainderRange13480 = range13480 % numThreads;
    /*Initialize inputs and vars.*/
    this->jglob_darts13480 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets13480[0] = _barrierCodelets13480(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13481* checkInCodelets13481Ptr = (this->checkInCodelets13481);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13481);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13481Ptr) = _checkInCodelets13481(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13481Ptr) = _checkInCodelets13481(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13481Ptr).decDep();
        checkInCodelets13481Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13480::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13481[localID].setID(codeletID);
    this->checkInCodelets13481[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13481[localID + this->baseNumThreads * i]
            = _checkInCodelets13481(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13481[localID + this->baseNumThreads * i]
            = _checkInCodelets13481(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13481[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13481[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP13480::~TP13480()
{
    delete[] jglob_darts13480;
    delete[] barrierCodelets13480;
    delete[] checkInCodelets13481;
}
/*TP13864: OMPParallelDirective*/
void TP13864::_barrierCodelets13864::fire(void)
{
    TP13864* myTP = static_cast<TP13864*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP13864::_checkInCodelets13871::fire(void)
{
    /*Allocate each array var in the codelet */
    this->inputsTPParent->ue_1jk_outer13864_size = 5;
    this->inputsTPParent->ue_1jk_darts13864[this->getID()]
        = (double*)malloc(sizeof(double) * this->inputsTPParent->ue_1jk_outer13864_size);
    /*Allocate each array var in the codelet */
    this->inputsTPParent->ue_i1k_outer13864_size = 5;
    this->inputsTPParent->ue_i1k_darts13864[this->getID()]
        = (double*)malloc(sizeof(double) * this->inputsTPParent->ue_i1k_outer13864_size);
    /*Allocate each array var in the codelet */
    this->inputsTPParent->ue_ij1_outer13864_size = 5;
    this->inputsTPParent->ue_ij1_darts13864[this->getID()]
        = (double*)malloc(sizeof(double) * this->inputsTPParent->ue_ij1_outer13864_size);
    /*Allocate each array var in the codelet */
    this->inputsTPParent->ue_ijnz_outer13864_size = 5;
    this->inputsTPParent->ue_ijnz_darts13864[this->getID()]
        = (double*)malloc(sizeof(double) * this->inputsTPParent->ue_ijnz_outer13864_size);
    /*Allocate each array var in the codelet */
    this->inputsTPParent->ue_iny0k_outer13864_size = 5;
    this->inputsTPParent->ue_iny0k_darts13864[this->getID()]
        = (double*)malloc(sizeof(double) * this->inputsTPParent->ue_iny0k_outer13864_size);
    /*Allocate each array var in the codelet */
    this->inputsTPParent->ue_nx0jk_outer13864_size = 5;
    this->inputsTPParent->ue_nx0jk_darts13864[this->getID()]
        = (double*)malloc(sizeof(double) * this->inputsTPParent->ue_nx0jk_outer13864_size);
    /*region 13871 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13871;
    if (idx < myTP->TPsToUse13871) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13871_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(ny - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13871;
            int minIteration = min<int>(ny, 0);
            int remainderRange = range % myTP->TPsToUse13871;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (0 < ny) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse13871 - 1) {
                lastIteration = ny;
            }
#if USEINVOKE == 1
            invoke<TP13871>(myTP, myTP->codeletsPerTP13871 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13871Ptr[idx]));
#else
            place<TP13871>(idx, myTP, myTP->codeletsPerTP13871 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13871Ptr[idx]));
#endif
        } else {
            if (myTP->TP13871Ptr[idx] != nullptr) {
                myTP->TP13871Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13864::_barrierCodelets13871::fire(void)
{
    TP13864* myTP = static_cast<TP13864*>(myTP_);
    myTP->TPParent->barrierCodelets13864[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets13864[0]));
}
TP13864::TP13864(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , eta_darts13864(new double[this->numThreads]) /*VARIABLE*/
    , i_darts13864(new int[this->numThreads]) /*VARIABLE*/
    , iglob_darts13864(new int[this->numThreads]) /*VARIABLE*/
    , j_darts13864(new int[this->numThreads]) /*VARIABLE*/
    , jglob_darts13864(new int[this->numThreads]) /*VARIABLE*/
    , k_darts13864(new int[this->numThreads]) /*VARIABLE*/
    , m_darts13864(new int[this->numThreads]) /*VARIABLE*/
    , peta_darts13864(new double[this->numThreads]) /*VARIABLE*/
    , pxi_darts13864(new double[this->numThreads]) /*VARIABLE*/
    , pzeta_darts13864(new double[this->numThreads]) /*VARIABLE*/
    , xi_darts13864(new double[this->numThreads]) /*VARIABLE*/
    , zeta_darts13864(new double[this->numThreads]) /*VARIABLE*/
    , TP13871Ptr(new TP13871*[NUMTPS13871])
    , TP13871_alreadyLaunched(new size_t[NUMTPS13871])
    , numTPsSet13871(0)
    , numTPsReady13871(0)
    , TPsToUse13871(NUMTPS13871)
    , codeletsPerTP13871(this->numThreads / NUMTPS13871)
    , totalCodelets13871(this->TPsToUse13871 * this->codeletsPerTP13871)
    , barrierCodelets13864(new _barrierCodelets13864[1])
    , checkInCodelets13871(new _checkInCodelets13871[this->numThreads])
    , barrierCodelets13871(new _barrierCodelets13871[1])
{
    /*Initialize inputs and vars.*/
    this->ue_1jk_darts13864 = (double**)malloc(sizeof(double*) * this->numThreads); /*VARIABLE*/
    this->ue_i1k_darts13864 = (double**)malloc(sizeof(double*) * this->numThreads); /*VARIABLE*/
    this->ue_ij1_darts13864 = (double**)malloc(sizeof(double*) * this->numThreads); /*VARIABLE*/
    this->ue_ijnz_darts13864 = (double**)malloc(sizeof(double*) * this->numThreads); /*VARIABLE*/
    this->ue_iny0k_darts13864 = (double**)malloc(sizeof(double*) * this->numThreads); /*VARIABLE*/
    this->ue_nx0jk_darts13864 = (double**)malloc(sizeof(double*) * this->numThreads); /*VARIABLE*/
    /*Initialize Codelets*/
    barrierCodelets13864[0] = _barrierCodelets13864(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets13871[0] = _barrierCodelets13871(NUMTPS13871, NUMTPS13871, this, 0);
    _checkInCodelets13871* checkInCodelets13871Ptr = (this->checkInCodelets13871);
    for (int i = 0; i < NUMTPS13871; i++) {
        TP13871Ptr[i] = nullptr;
        TP13871_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets13871Ptr) = _checkInCodelets13871(1, 1, this, codeletCounter);
        (*checkInCodelets13871Ptr).decDep();
        checkInCodelets13871Ptr++;
    }
}
TP13864::~TP13864()
{
    delete[] eta_darts13864;
    delete[] i_darts13864;
    delete[] iglob_darts13864;
    delete[] j_darts13864;
    delete[] jglob_darts13864;
    delete[] k_darts13864;
    delete[] m_darts13864;
    delete[] peta_darts13864;
    delete[] pxi_darts13864;
    delete[] pzeta_darts13864;
    delete[] xi_darts13864;
    delete[] zeta_darts13864;
    delete[] barrierCodelets13864;
    delete[] barrierCodelets13871;
    delete[] checkInCodelets13871;
}
/*TP13871: OMPForDirective*/
void TP13871::_barrierCodelets13871::fire(void)
{
    TP13871* myTP = static_cast<TP13871*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13871[0].decDep();
}
bool TP13871::requestNewRangeIterations13871(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13871 * codeletID;
        int tempEndRange = rangePerCodelet13871 * (codeletID + 1);
        if (remainderRange13871 != 0) {
            if (codeletID < (uint32_t)remainderRange13871) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13871;
                tempEndRange += remainderRange13871;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13871;
        tempEndRange = tempEndRange * 1 + minIteration13871;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13871 < lastIteration13871) {
            (this->inputsTPParent->j_darts13871[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->j_darts13871[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13871;
        }
    }
    return isThereNewIteration;
}
void TP13871::_checkInCodelets13872::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->eta_darts13871[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->eta_darts13864[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->iglob_darts13871[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts13864[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jglob_darts13871[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts13864[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->peta_darts13871[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->peta_darts13864[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->pxi_darts13871[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->pxi_darts13864[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->pzeta_darts13871[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->pzeta_darts13864[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ue_1jk_darts13871[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->ue_1jk_darts13864[this->getID()][0]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ue_i1k_darts13871[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->ue_i1k_darts13864[this->getID()][0]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ue_ij1_darts13871[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->ue_ij1_darts13864[this->getID()][0]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ue_ijnz_darts13871[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->ue_ijnz_darts13864[this->getID()][0]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ue_iny0k_darts13871[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->ue_iny0k_darts13864[this->getID()][0]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ue_nx0jk_darts13871[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->ue_nx0jk_darts13864[this->getID()][0]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->xi_darts13871[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->xi_darts13864[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->zeta_darts13871[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->zeta_darts13864[this->getID()]);

    /*printing node 13872: ForStmt*/
    /*var: eta*/
    /*var: i*/
    /*var: iglob*/
    /*var: j*/
    /*var: jglob*/
    /*var: k*/
    /*var: m*/
    /*var: peta*/
    /*var: pxi*/
    /*var: pzeta*/
    /*var: ue_1jk*/
    /*var: ue_i1k*/
    /*var: ue_ij1*/
    /*var: ue_ijnz*/
    /*var: ue_iny0k*/
    /*var: ue_nx0jk*/
    /*var: xi*/
    /*var: zeta*/
    double** eta = &(this->inputsTPParent->eta_darts13871[this->getLocalID()]);
    (void)eta /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts13871[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iglob = &(this->inputsTPParent->iglob_darts13871[this->getLocalID()]);
    (void)iglob /*OMP_SHARED_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts13871[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jglob = &(this->inputsTPParent->jglob_darts13871[this->getLocalID()]);
    (void)jglob /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts13871[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts13871[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** peta = &(this->inputsTPParent->peta_darts13871[this->getLocalID()]);
    (void)peta /*OMP_SHARED_PRIVATE*/;
    double** pxi = &(this->inputsTPParent->pxi_darts13871[this->getLocalID()]);
    (void)pxi /*OMP_SHARED_PRIVATE*/;
    double** pzeta = &(this->inputsTPParent->pzeta_darts13871[this->getLocalID()]);
    (void)pzeta /*OMP_SHARED_PRIVATE*/;
    double** ue_1jk = &(this->inputsTPParent->ue_1jk_darts13871[this->getLocalID()]);
    (void)ue_1jk /*OMP_SHARED_PRIVATE*/;
    double** ue_i1k = &(this->inputsTPParent->ue_i1k_darts13871[this->getLocalID()]);
    (void)ue_i1k /*OMP_SHARED_PRIVATE*/;
    double** ue_ij1 = &(this->inputsTPParent->ue_ij1_darts13871[this->getLocalID()]);
    (void)ue_ij1 /*OMP_SHARED_PRIVATE*/;
    double** ue_ijnz = &(this->inputsTPParent->ue_ijnz_darts13871[this->getLocalID()]);
    (void)ue_ijnz /*OMP_SHARED_PRIVATE*/;
    double** ue_iny0k = &(this->inputsTPParent->ue_iny0k_darts13871[this->getLocalID()]);
    (void)ue_iny0k /*OMP_SHARED_PRIVATE*/;
    double** ue_nx0jk = &(this->inputsTPParent->ue_nx0jk_darts13871[this->getLocalID()]);
    (void)ue_nx0jk /*OMP_SHARED_PRIVATE*/;
    double** xi = &(this->inputsTPParent->xi_darts13871[this->getLocalID()]);
    (void)xi /*OMP_SHARED_PRIVATE*/;
    double** zeta = &(this->inputsTPParent->zeta_darts13871[this->getLocalID()]);
    (void)zeta /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13871(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13871[0].decDep();
        return;
    }
    for (int j_darts_counter_temp13871 = (*j); j_darts_counter_temp13871 < endRange
         && j_darts_counter_temp13871 < this->inputsTPParent->lastIteration13871;
         j_darts_counter_temp13871++) {
        {
            (*(*jglob)) = (j_darts_counter_temp13871);
            {
                /*Loop's init*/
                (*k) = 1;
                int k_darts_counter_temp13871 = (*k);
                for (; k_darts_counter_temp13871 < nz - 1; k_darts_counter_temp13871++) {
                    (*(*zeta)) = ((double)k_darts_counter_temp13871) / (nz - 1);
                    if ((*(*jglob)) != 0 && (*(*jglob)) != ny0 - 1) {
                        (*(*eta)) = ((double)((*(*jglob)))) / (ny0 - 1);
                        {
                            /*Loop's init*/
                            (*i) = 0;
                            int i_darts_counter_temp13871 = (*i);
                            for (; i_darts_counter_temp13871 < nx; i_darts_counter_temp13871++) {
                                (*(*iglob)) = i_darts_counter_temp13871;
                                if ((*(*iglob)) != 0 && (*(*iglob)) != nx0 - 1) {
                                    (*(*xi)) = ((double)((*(*iglob)))) / (nx0 - 1);
                                    exact(0, (*(*jglob)), k_darts_counter_temp13871, ((*ue_1jk)));
                                    exact(nx0 - 1, (*(*jglob)), k_darts_counter_temp13871,
                                        ((*ue_nx0jk)));
                                    exact((*(*iglob)), 0, k_darts_counter_temp13871, ((*ue_i1k)));
                                    exact((*(*iglob)), ny0 - 1, k_darts_counter_temp13871,
                                        ((*ue_iny0k)));
                                    exact((*(*iglob)), (*(*jglob)), 0, ((*ue_ij1)));
                                    exact((*(*iglob)), (*(*jglob)), nz - 1, ((*ue_ijnz)));
                                    {
                                        /*Loop's init*/
                                        (*m) = 0;
                                        int m_darts_counter_temp13871 = (*m);
                                        for (; m_darts_counter_temp13871 < 5;
                                             m_darts_counter_temp13871++) {
                                            (*(*pxi)) = (1. - (*(*xi)))
                                                    * ((*ue_1jk))[m_darts_counter_temp13871]
                                                + (*(*xi))
                                                    * ((*ue_nx0jk))[m_darts_counter_temp13871];
                                            (*(*peta)) = (1. - (*(*eta)))
                                                    * ((*ue_i1k))[m_darts_counter_temp13871]
                                                + (*(*eta))
                                                    * ((*ue_iny0k))[m_darts_counter_temp13871];
                                            (*(*pzeta)) = (1. - (*(*zeta)))
                                                    * ((*ue_ij1))[m_darts_counter_temp13871]
                                                + (*(*zeta))
                                                    * ((*ue_ijnz))[m_darts_counter_temp13871];
                                            u[i_darts_counter_temp13871]
                                             [(j_darts_counter_temp13871)]
                                             [k_darts_counter_temp13871][m_darts_counter_temp13871]
                                                = (*(*pxi)) + (*(*peta)) + (*(*pzeta))
                                                - (*(*pxi)) * (*(*peta)) - (*(*peta)) * (*(*pzeta))
                                                - (*(*pzeta)) * (*(*pxi))
                                                + (*(*pxi)) * (*(*peta)) * (*(*pzeta));
                                        }
                                        (*m) = m_darts_counter_temp13871;
                                    }
                                }
                            }
                            (*i) = i_darts_counter_temp13871;
                        }
                    }
                }
                (*k) = k_darts_counter_temp13871;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13871[0].decDep();
}
TP13871::TP13871(int in_numThreads, int in_mainCodeletID, TP13864* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13871** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , eta_darts13871(new double*[this->numThreads])
    , i_darts13871(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iglob_darts13871(new int*[this->numThreads])
    , j_darts13871(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jglob_darts13871(new int*[this->numThreads])
    , k_darts13871(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts13871(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , peta_darts13871(new double*[this->numThreads])
    , pxi_darts13871(new double*[this->numThreads])
    , pzeta_darts13871(new double*[this->numThreads])
    , ue_1jk_darts13871(new double*[this->numThreads][5])
    , ue_i1k_darts13871(new double*[this->numThreads][5])
    , ue_ij1_darts13871(new double*[this->numThreads][5])
    , ue_ijnz_darts13871(new double*[this->numThreads][5])
    , ue_iny0k_darts13871(new double*[this->numThreads][5])
    , ue_nx0jk_darts13871(new double*[this->numThreads][5])
    , xi_darts13871(new double*[this->numThreads])
    , zeta_darts13871(new double*[this->numThreads])
    , initIteration13871(in_initIteration)
    , lastIteration13871(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13871(new _barrierCodelets13871[1])
    , checkInCodelets13872(new _checkInCodelets13872[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13871 = abs(lastIteration13871 - initIteration13871) / 1;
    rangePerCodelet13871 = range13871 / numThreads;
    minIteration13871 = min<int>(lastIteration13871, initIteration13871);
    remainderRange13871 = range13871 % numThreads;
    /*Initialize inputs and vars.*/
    this->eta_darts13871
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->iglob_darts13871 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jglob_darts13871 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->peta_darts13871
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->pxi_darts13871
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->pzeta_darts13871
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ue_1jk_darts13871
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ue_i1k_darts13871
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ue_ij1_darts13871
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ue_ijnz_darts13871
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ue_iny0k_darts13871
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ue_nx0jk_darts13871
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->xi_darts13871
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->zeta_darts13871
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets13871[0] = _barrierCodelets13871(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13872* checkInCodelets13872Ptr = (this->checkInCodelets13872);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13872);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13872Ptr) = _checkInCodelets13872(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13872Ptr) = _checkInCodelets13872(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13872Ptr).decDep();
        checkInCodelets13872Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13871::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13872[localID].setID(codeletID);
    this->checkInCodelets13872[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13872[localID + this->baseNumThreads * i]
            = _checkInCodelets13872(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13872[localID + this->baseNumThreads * i]
            = _checkInCodelets13872(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13872[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13872[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP13871::~TP13871()
{
    delete[] eta_darts13871;
    delete[] iglob_darts13871;
    delete[] jglob_darts13871;
    delete[] peta_darts13871;
    delete[] pxi_darts13871;
    delete[] pzeta_darts13871;
    delete[] ue_1jk_darts13871;
    delete[] ue_i1k_darts13871;
    delete[] ue_ij1_darts13871;
    delete[] ue_ijnz_darts13871;
    delete[] ue_iny0k_darts13871;
    delete[] ue_nx0jk_darts13871;
    delete[] xi_darts13871;
    delete[] zeta_darts13871;
    delete[] barrierCodelets13871;
    delete[] checkInCodelets13872;
}
/*TP13997: OMPParallelDirective*/
void TP13997::_barrierCodelets13997::fire(void)
{
    TP13997* myTP = static_cast<TP13997*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP13997::_checkInCodelets13999::fire(void)
{
    /*region 13999 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13999;
    if (idx < myTP->TPsToUse13999) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13999_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(162 - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13999;
            int minIteration = min<int>(162, 0);
            int remainderRange = range % myTP->TPsToUse13999;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (0 < 162) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse13999 - 1) {
                lastIteration = 162;
            }
#if USEINVOKE == 1
            invoke<TP13999>(myTP, myTP->codeletsPerTP13999 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13999Ptr[idx]));
#else
            place<TP13999>(idx, myTP, myTP->codeletsPerTP13999 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13999Ptr[idx]));
#endif
        } else {
            if (myTP->TP13999Ptr[idx] != nullptr) {
                myTP->TP13999Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13997::_barrierCodelets13999::fire(void)
{
    TP13997* myTP = static_cast<TP13997*>(myTP_);
    myTP->TPParent->barrierCodelets13997[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets13997[0]));
}
TP13997::TP13997(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13997(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts13997(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts13997(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts13997(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , TP13999Ptr(new TP13999*[NUMTPS13999])
    , TP13999_alreadyLaunched(new size_t[NUMTPS13999])
    , numTPsSet13999(0)
    , numTPsReady13999(0)
    , TPsToUse13999(NUMTPS13999)
    , codeletsPerTP13999(this->numThreads / NUMTPS13999)
    , totalCodelets13999(this->TPsToUse13999 * this->codeletsPerTP13999)
    , barrierCodelets13997(new _barrierCodelets13997[1])
    , checkInCodelets13999(new _checkInCodelets13999[this->numThreads])
    , barrierCodelets13999(new _barrierCodelets13999[1])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets13997[0] = _barrierCodelets13997(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets13999[0] = _barrierCodelets13999(NUMTPS13999, NUMTPS13999, this, 0);
    _checkInCodelets13999* checkInCodelets13999Ptr = (this->checkInCodelets13999);
    for (int i = 0; i < NUMTPS13999; i++) {
        TP13999Ptr[i] = nullptr;
        TP13999_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets13999Ptr) = _checkInCodelets13999(1, 1, this, codeletCounter);
        (*checkInCodelets13999Ptr).decDep();
        checkInCodelets13999Ptr++;
    }
}
TP13997::~TP13997()
{
    delete[] barrierCodelets13997;
    delete[] barrierCodelets13999;
    delete[] checkInCodelets13999;
}
/*TP13999: OMPForDirective*/
void TP13999::_barrierCodelets13999::fire(void)
{
    TP13999* myTP = static_cast<TP13999*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13999[0].decDep();
}
bool TP13999::requestNewRangeIterations13999(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13999 * codeletID;
        int tempEndRange = rangePerCodelet13999 * (codeletID + 1);
        if (remainderRange13999 != 0) {
            if (codeletID < (uint32_t)remainderRange13999) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13999;
                tempEndRange += remainderRange13999;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13999;
        tempEndRange = tempEndRange * 1 + minIteration13999;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13999 < lastIteration13999) {
            (this->inputsTPParent->i_darts13999[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts13999[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13999;
        }
    }
    return isThereNewIteration;
}
void TP13999::_checkInCodelets14000::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 14000: ForStmt*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    int* i = &(this->inputsTPParent->i_darts13999[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts13999[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts13999[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts13999[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13999(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13999[0].decDep();
        return;
    }
    for (int i_darts_counter_temp13999 = (*i); i_darts_counter_temp13999 < endRange
         && i_darts_counter_temp13999 < this->inputsTPParent->lastIteration13999;
         i_darts_counter_temp13999++) {
        {
            {
                /*Loop's init*/
                (*j) = 0;
                int j_darts_counter_temp13999 = (*j);
                for (; j_darts_counter_temp13999 < 162; j_darts_counter_temp13999++) {
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp13999 = (*k);
                        for (; k_darts_counter_temp13999 < 5; k_darts_counter_temp13999++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp13999 = (*m);
                                for (; m_darts_counter_temp13999 < 5; m_darts_counter_temp13999++) {
                                    a[(i_darts_counter_temp13999)][j_darts_counter_temp13999]
                                     [k_darts_counter_temp13999][m_darts_counter_temp13999]
                                        = 0.;
                                    b[(i_darts_counter_temp13999)][j_darts_counter_temp13999]
                                     [k_darts_counter_temp13999][m_darts_counter_temp13999]
                                        = 0.;
                                    c[(i_darts_counter_temp13999)][j_darts_counter_temp13999]
                                     [k_darts_counter_temp13999][m_darts_counter_temp13999]
                                        = 0.;
                                    d[(i_darts_counter_temp13999)][j_darts_counter_temp13999]
                                     [k_darts_counter_temp13999][m_darts_counter_temp13999]
                                        = 0.;
                                }
                                (*m) = m_darts_counter_temp13999;
                            }
                        }
                        (*k) = k_darts_counter_temp13999;
                    }
                }
                (*j) = j_darts_counter_temp13999;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13999[0].decDep();
}
TP13999::TP13999(int in_numThreads, int in_mainCodeletID, TP13997* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13999** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13999(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts13999(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts13999(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts13999(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration13999(in_initIteration)
    , lastIteration13999(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13999(new _barrierCodelets13999[1])
    , checkInCodelets14000(new _checkInCodelets14000[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13999 = abs(lastIteration13999 - initIteration13999) / 1;
    rangePerCodelet13999 = range13999 / numThreads;
    minIteration13999 = min<int>(lastIteration13999, initIteration13999);
    remainderRange13999 = range13999 % numThreads;
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets13999[0] = _barrierCodelets13999(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets14000* checkInCodelets14000Ptr = (this->checkInCodelets14000);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets14000);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets14000Ptr) = _checkInCodelets14000(2, 1, this, codeletCounter);
#else
        (*checkInCodelets14000Ptr) = _checkInCodelets14000(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets14000Ptr).decDep();
        checkInCodelets14000Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13999::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets14000[localID].setID(codeletID);
    this->checkInCodelets14000[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets14000[localID + this->baseNumThreads * i]
            = _checkInCodelets14000(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets14000[localID + this->baseNumThreads * i]
            = _checkInCodelets14000(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets14000[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets14000[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP13999::~TP13999()
{
    delete[] barrierCodelets13999;
    delete[] checkInCodelets14000;
}
/*TP14096: OMPParallelDirective*/
void TP14096::_barrierCodelets14096::fire(void)
{
    TP14096* myTP = static_cast<TP14096*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP14096::_checkInCodelets14098::fire(void)
{
    /*region 14098 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP14098;
    if (idx < myTP->TPsToUse14098) {
        if (!__sync_val_compare_and_swap(&(myTP->TP14098_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse14098;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse14098;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (ist < iend) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse14098 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse14098 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP14098>(myTP, myTP->codeletsPerTP14098 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP14098Ptr[idx]));
#else
            place<TP14098>(idx, myTP, myTP->codeletsPerTP14098 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP14098Ptr[idx]));
#endif
        } else {
            if (myTP->TP14098Ptr[idx] != nullptr) {
                myTP->TP14098Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP14096::_barrierCodelets14098::fire(void)
{
    TP14096* myTP = static_cast<TP14096*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets14154[codeletsCounter].decDep();
        }
    }
}
void TP14096::_checkInCodelets14154::fire(void)
{

    /*printing node 14154: BinaryOperator*/
    (this->inputsTPParent->k_darts14096[this->getID()]) = 1;

    /*printing node 14155: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts14096[this->getID()]) <= nz - 2) {
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets14153[this->getID()].decDep();
        return;
    } else {
        /*Signal the codelet after the loop from the end condional node.*/
        /*Signaling next codelet region: 14157 nextRegion: 14161 */
        myTP->controlTPParent->barrierCodelets14161[0].decDep();
        return;
    }
}
void TP14096::_checkInCodelets14153::fire(void)
{

    /*printing node 14153: ForStmt*/
    bool haveToLaunch = __sync_bool_compare_and_swap(&(myTP->controlTPParent->TP14153_LoopCounter),
        myTP->controlTPParent->TP14153_LoopCounterPerThread[this->getID()],
        myTP->controlTPParent->TP14153_LoopCounterPerThread[this->getID()] + 1);
    unsigned int iterIdx = myTP->controlTPParent->TP14153_LoopCounterPerThread[this->getID()];
    if (haveToLaunch) {
        this->resetCodelet();
        myTP->controlTPParent->TP14153PtrVec.push_back(nullptr);
        myTP->controlTPParent->TP14153_LoopCounterPerThread[this->getID()] += 1;
        invoke<TP14153>(myTP, myTP->numThreads, this->getID(), myTP, myTP->inputsTPParent,
            &(myTP->controlTPParent->TP14153PtrVec.back()));
    } else {
        if (myTP->controlTPParent->TP14153PtrVec.size() == 0) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP14153PtrVec.size() < (iterIdx + 1)) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP14153PtrVec[iterIdx] == nullptr) {
            this->resetCodelet();
            this->decDep();
            return;
        } else {
            this->resetCodelet();
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14153PtrVec[iterIdx]->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14153PtrVec[iterIdx]->availableCodelets[this->getID()] = 1;
#endif
            myTP->controlTPParent->TP14153_LoopCounterPerThread[this->getID()] += 1;
        }
    }
}
void TP14096::_checkInCodelets14157::fire(void)
{

    /*printing node 14157: UnaryOperator*/
    (this->inputsTPParent->k_darts14096[this->getID()])++;

    /*printing node 14612: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts14096[this->getID()]) <= nz - 2) {
        this->resetCodelet();
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets14153[this->getID()].decDep();
        return;
    } else {
        /*Signal the codelet after the loop from the condtional node.*/
        /*Signaling next codelet region: 14157 nextRegion: 14161 */
        myTP->controlTPParent->barrierCodelets14161[0].decDep();
        return;
    }
}
void TP14096::_barrierCodelets14161::fire(void)
{
    TP14096* myTP = static_cast<TP14096*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets14163[codeletsCounter].decDep();
        }
    }
}
void TP14096::_checkInCodelets14163::fire(void)
{

    /*printing node 14163: BinaryOperator*/
    (this->inputsTPParent->k_darts14096[this->getID()]) = nz - 2;

    /*printing node 14165: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts14096[this->getID()]) >= 1) {
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets14162[this->getID()].decDep();
        return;
    } else {
        /*Signal the codelet after the loop from the end condional node.*/
        /*Signaling next codelet region: 14166 nextRegion: 14170 */
        myTP->controlTPParent->barrierCodelets14170[0].decDep();
        return;
    }
}
void TP14096::_checkInCodelets14162::fire(void)
{

    /*printing node 14162: ForStmt*/
    bool haveToLaunch = __sync_bool_compare_and_swap(&(myTP->controlTPParent->TP14162_LoopCounter),
        myTP->controlTPParent->TP14162_LoopCounterPerThread[this->getID()],
        myTP->controlTPParent->TP14162_LoopCounterPerThread[this->getID()] + 1);
    unsigned int iterIdx = myTP->controlTPParent->TP14162_LoopCounterPerThread[this->getID()];
    if (haveToLaunch) {
        this->resetCodelet();
        myTP->controlTPParent->TP14162PtrVec.push_back(nullptr);
        myTP->controlTPParent->TP14162_LoopCounterPerThread[this->getID()] += 1;
        invoke<TP14162>(myTP, myTP->numThreads, this->getID(), myTP, myTP->inputsTPParent,
            &(myTP->controlTPParent->TP14162PtrVec.back()));
    } else {
        if (myTP->controlTPParent->TP14162PtrVec.size() == 0) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP14162PtrVec.size() < (iterIdx + 1)) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP14162PtrVec[iterIdx] == nullptr) {
            this->resetCodelet();
            this->decDep();
            return;
        } else {
            this->resetCodelet();
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14162PtrVec[iterIdx]->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14162PtrVec[iterIdx]->availableCodelets[this->getID()] = 1;
#endif
            myTP->controlTPParent->TP14162_LoopCounterPerThread[this->getID()] += 1;
        }
    }
}
void TP14096::_checkInCodelets14166::fire(void)
{

    /*printing node 14166: UnaryOperator*/
    (this->inputsTPParent->k_darts14096[this->getID()])--;

    /*printing node 14613: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts14096[this->getID()]) >= 1) {
        this->resetCodelet();
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets14162[this->getID()].decDep();
        return;
    } else {
        /*Signal the codelet after the loop from the condtional node.*/
        /*Signaling next codelet region: 14166 nextRegion: 14170 */
        myTP->controlTPParent->barrierCodelets14170[0].decDep();
        return;
    }
}
void TP14096::_barrierCodelets14170::fire(void)
{
    TP14096* myTP = static_cast<TP14096*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets14171[codeletsCounter].decDep();
        }
    }
}
void TP14096::_checkInCodelets14171::fire(void)
{
    /*region 14171 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP14171;
    if (idx < myTP->TPsToUse14171) {
        if (!__sync_val_compare_and_swap(&(myTP->TP14171_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse14171;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse14171;
            int initIteration = rangePerCodelet * idx;
            int lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (ist < iend) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse14171 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse14171 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP14171>(myTP, myTP->codeletsPerTP14171 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(*(this->inputsTPParent->tmp_darts14096)),
                &(myTP->TP14171Ptr[idx]));
#else
            place<TP14171>(idx, myTP, myTP->codeletsPerTP14171 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(*(this->inputsTPParent->tmp_darts14096)),
                &(myTP->TP14171Ptr[idx]));
#endif
        } else {
            if (myTP->TP14171Ptr[idx] != nullptr) {
                myTP->TP14171Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP14096::_barrierCodelets14171::fire(void)
{
    TP14096* myTP = static_cast<TP14096*>(myTP_);
    myTP->TPParent->barrierCodelets14096[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets14096[0]));
}
TP14096::TP14096(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet, double* in_tmp,
    double* in_tv, int in_tv_outer14096_size)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts14096(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , istep_darts14096(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts14096(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts14096(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts14096(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts14096(in_tmp) /*OMP_SHARED - INPUT*/
    , tv_darts14096(in_tv) /*OMP_SHARED - INPUT*/
    , tv_outer14096_size(in_tv_outer14096_size)
    , TP14098Ptr(new TP14098*[NUMTPS14098])
    , TP14098_alreadyLaunched(new size_t[NUMTPS14098])
    , numTPsSet14098(0)
    , numTPsReady14098(0)
    , TPsToUse14098(NUMTPS14098)
    , codeletsPerTP14098(this->numThreads / NUMTPS14098)
    , totalCodelets14098(this->TPsToUse14098 * this->codeletsPerTP14098)
    , TP14153_LoopCounter(0)
    , TP14153_LoopCounterPerThread(new unsigned int[this->numThreads])
    , TP14162_LoopCounter(0)
    , TP14162_LoopCounterPerThread(new unsigned int[this->numThreads])
    , TP14171Ptr(new TP14171*[NUMTPS14171])
    , TP14171_alreadyLaunched(new size_t[NUMTPS14171])
    , numTPsSet14171(0)
    , numTPsReady14171(0)
    , TPsToUse14171(NUMTPS14171)
    , codeletsPerTP14171(this->numThreads / NUMTPS14171)
    , totalCodelets14171(this->TPsToUse14171 * this->codeletsPerTP14171)
    , barrierCodelets14096(new _barrierCodelets14096[1])
    , checkInCodelets14098(new _checkInCodelets14098[this->numThreads])
    , barrierCodelets14098(new _barrierCodelets14098[1])
    , checkInCodelets14154(new _checkInCodelets14154[this->numThreads])
    , checkInCodelets14153(new _checkInCodelets14153[this->numThreads])
    , checkInCodelets14157(new _checkInCodelets14157[this->numThreads])
    , barrierCodelets14161(new _barrierCodelets14161[1])
    , checkInCodelets14163(new _checkInCodelets14163[this->numThreads])
    , checkInCodelets14162(new _checkInCodelets14162[this->numThreads])
    , checkInCodelets14166(new _checkInCodelets14166[this->numThreads])
    , barrierCodelets14170(new _barrierCodelets14170[1])
    , checkInCodelets14171(new _checkInCodelets14171[this->numThreads])
    , barrierCodelets14171(new _barrierCodelets14171[1])
{
    memset((void*)TP14153_LoopCounterPerThread, 0, this->numThreads * sizeof(unsigned int));
    memset((void*)TP14162_LoopCounterPerThread, 0, this->numThreads * sizeof(unsigned int));
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets14096[0] = _barrierCodelets14096(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets14171[0] = _barrierCodelets14171(NUMTPS14171, NUMTPS14171, this, 0);
    barrierCodelets14170[0] = _barrierCodelets14170(this->numThreads, this->numThreads, this, 0);
    barrierCodelets14161[0] = _barrierCodelets14161(this->numThreads, this->numThreads, this, 0);
    barrierCodelets14098[0] = _barrierCodelets14098(NUMTPS14098, NUMTPS14098, this, 0);
    _checkInCodelets14171* checkInCodelets14171Ptr = (this->checkInCodelets14171);
    for (int i = 0; i < NUMTPS14171; i++) {
        TP14171Ptr[i] = nullptr;
        TP14171_alreadyLaunched[i] = 0;
    }
    _checkInCodelets14166* checkInCodelets14166Ptr = (this->checkInCodelets14166);
    _checkInCodelets14162* checkInCodelets14162Ptr = (this->checkInCodelets14162);
    _checkInCodelets14163* checkInCodelets14163Ptr = (this->checkInCodelets14163);
    _checkInCodelets14157* checkInCodelets14157Ptr = (this->checkInCodelets14157);
    _checkInCodelets14153* checkInCodelets14153Ptr = (this->checkInCodelets14153);
    _checkInCodelets14154* checkInCodelets14154Ptr = (this->checkInCodelets14154);
    _checkInCodelets14098* checkInCodelets14098Ptr = (this->checkInCodelets14098);
    for (int i = 0; i < NUMTPS14098; i++) {
        TP14098Ptr[i] = nullptr;
        TP14098_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets14171Ptr) = _checkInCodelets14171(1, 1, this, codeletCounter);
        checkInCodelets14171Ptr++;
        (*checkInCodelets14166Ptr) = _checkInCodelets14166(1, 1, this, codeletCounter);
        checkInCodelets14166Ptr++;
        (*checkInCodelets14162Ptr) = _checkInCodelets14162(1, 1, this, codeletCounter);
        checkInCodelets14162Ptr++;
        (*checkInCodelets14163Ptr) = _checkInCodelets14163(1, 1, this, codeletCounter);
        checkInCodelets14163Ptr++;
        (*checkInCodelets14157Ptr) = _checkInCodelets14157(1, 1, this, codeletCounter);
        checkInCodelets14157Ptr++;
        (*checkInCodelets14153Ptr) = _checkInCodelets14153(1, 1, this, codeletCounter);
        checkInCodelets14153Ptr++;
        (*checkInCodelets14154Ptr) = _checkInCodelets14154(1, 1, this, codeletCounter);
        checkInCodelets14154Ptr++;
        (*checkInCodelets14098Ptr) = _checkInCodelets14098(1, 1, this, codeletCounter);
        (*checkInCodelets14098Ptr).decDep();
        checkInCodelets14098Ptr++;
    }
}
TP14096::~TP14096()
{
    delete[] TP14153_LoopCounterPerThread;
    delete[] TP14162_LoopCounterPerThread;
    delete[] barrierCodelets14096;
    delete[] barrierCodelets14171;
    delete[] checkInCodelets14171;
    delete[] barrierCodelets14170;
    delete[] checkInCodelets14166;
    delete[] checkInCodelets14162;
    delete[] checkInCodelets14163;
    delete[] barrierCodelets14161;
    delete[] checkInCodelets14157;
    delete[] checkInCodelets14153;
    delete[] checkInCodelets14154;
    delete[] barrierCodelets14098;
    delete[] checkInCodelets14098;
}
/*TP14098: OMPForDirective*/
void TP14098::_barrierCodelets14098::fire(void)
{
    TP14098* myTP = static_cast<TP14098*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets14098[0].decDep();
}
bool TP14098::requestNewRangeIterations14098(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet14098 * codeletID;
        int tempEndRange = rangePerCodelet14098 * (codeletID + 1);
        if (remainderRange14098 != 0) {
            if (codeletID < (uint32_t)remainderRange14098) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange14098;
                tempEndRange += remainderRange14098;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration14098;
        tempEndRange = tempEndRange * 1 + minIteration14098;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration14098 < lastIteration14098) {
            (this->inputsTPParent->i_darts14098[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts14098[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration14098;
        }
    }
    return isThereNewIteration;
}
void TP14098::_checkInCodelets14099::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 14099: ForStmt*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    int* i = &(this->inputsTPParent->i_darts14098[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts14098[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts14098[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts14098[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations14098(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets14098[0].decDep();
        return;
    }
    for (int i_darts_counter_temp14098 = (*i); i_darts_counter_temp14098 <= endRange
         && i_darts_counter_temp14098 <= this->inputsTPParent->lastIteration14098;
         i_darts_counter_temp14098++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp14098 = (*j);
                for (; j_darts_counter_temp14098 <= jend; j_darts_counter_temp14098++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp14098 = (*k);
                        for (; k_darts_counter_temp14098 <= nz - 2; k_darts_counter_temp14098++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp14098 = (*m);
                                for (; m_darts_counter_temp14098 < 5; m_darts_counter_temp14098++) {
                                    rsd[(i_darts_counter_temp14098)][j_darts_counter_temp14098]
                                       [k_darts_counter_temp14098][m_darts_counter_temp14098]
                                        = dt
                                        * rsd[(i_darts_counter_temp14098)]
                                             [j_darts_counter_temp14098][k_darts_counter_temp14098]
                                             [m_darts_counter_temp14098];
                                }
                                (*m) = m_darts_counter_temp14098;
                            }
                        }
                        (*k) = k_darts_counter_temp14098;
                    }
                }
                (*j) = j_darts_counter_temp14098;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets14098[0].decDep();
}
TP14098::TP14098(int in_numThreads, int in_mainCodeletID, TP14096* in_TPParent,
    int in_initIteration, int in_lastIteration, TP14098** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts14098(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts14098(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts14098(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts14098(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration14098(in_initIteration)
    , lastIteration14098(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets14098(new _barrierCodelets14098[1])
    , checkInCodelets14099(new _checkInCodelets14099[this->numThreads])
{
    /*Initialize the loop parameters*/
    range14098 = abs(lastIteration14098 - initIteration14098) / 1;
    rangePerCodelet14098 = range14098 / numThreads;
    minIteration14098 = min<int>(lastIteration14098, initIteration14098);
    remainderRange14098 = range14098 % numThreads;
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets14098[0] = _barrierCodelets14098(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets14099* checkInCodelets14099Ptr = (this->checkInCodelets14099);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets14099);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets14099Ptr) = _checkInCodelets14099(2, 1, this, codeletCounter);
#else
        (*checkInCodelets14099Ptr) = _checkInCodelets14099(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets14099Ptr).decDep();
        checkInCodelets14099Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP14098::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets14099[localID].setID(codeletID);
    this->checkInCodelets14099[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets14099[localID + this->baseNumThreads * i]
            = _checkInCodelets14099(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets14099[localID + this->baseNumThreads * i]
            = _checkInCodelets14099(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets14099[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets14099[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP14098::~TP14098()
{
    delete[] barrierCodelets14098;
    delete[] checkInCodelets14099;
}
/*TP14153: ForStmt*/
void TP14153::_checkInCodelets14159::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif

    /*printing node 14159: CallExpr*/
    if (!__sync_val_compare_and_swap(&(myTP->controlTPParent->TP14159_alreadyLaunched), 0, 1)) {
        /*Make the function call*/
        invoke<TP_jacld>(myTP, myTP->numThreads, this->getID(),
            &(myTP->controlTPParent->checkInCodelets14160[this->getID()]),
            &(myTP->controlTPParent->TPParent->barrierCodelets14161[0]),
            &(myTP->controlTPParent->TP14159Ptr),
            (this->inputsTPParent->k_darts14096[this->getID()]));
    } else {
        if (myTP->controlTPParent->TP14159Ptr == nullptr) {
            myTP->add(this);
            return;
        } else {
            myTP->controlTPParent->TP14159Ptr->setNewInputs(
                (this->inputsTPParent->k_darts14096[this->getID()]), this->getID());
            myTP->controlTPParent->TP14159Ptr->nextCodeletsjacld[this->getID()]
                = &(myTP->controlTPParent->checkInCodelets14160[this->getID()]);
            myTP->controlTPParent->TP14159Ptr->nextSyncCodeletsjacld[this->getID()]
                = &(myTP->controlTPParent->TPParent->barrierCodelets14161[0]);
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14159Ptr->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14159Ptr->availableCodelets[this->getID()] = 1;
#endif
        }
    }
}
void TP14153::_checkInCodelets14160::fire(void)
{

    /*printing node 14160: CallExpr*/
    if (!__sync_val_compare_and_swap(&(myTP->controlTPParent->TP14160_alreadyLaunched), 0, 1)) {
        /*Make the function call*/
        invoke<TP_blts>(myTP, myTP->numThreads, this->getID(),
            &(myTP->controlTPParent->TPParent->checkInCodelets14157[this->getID()]),
            &(myTP->controlTPParent->TPParent->barrierCodelets14161[0]),
            &(myTP->controlTPParent->TP14160Ptr), nx, ny, nz,
            (this->inputsTPParent->k_darts14096[this->getID()]), omega, rsd, a, b, c, d, ist, iend,
            jst, jend, nx0, ny0);
    } else {
        if (myTP->controlTPParent->TP14160Ptr == nullptr) {
            myTP->add(this);
            return;
        } else {
            myTP->controlTPParent->TP14160Ptr->setNewInputs(nx, ny, nz,
                (this->inputsTPParent->k_darts14096[this->getID()]), omega, rsd, a, b, c, d, ist,
                iend, jst, jend, nx0, ny0, this->getID());
            myTP->controlTPParent->TP14160Ptr->nextCodeletsblts[this->getID()]
                = &(myTP->controlTPParent->TPParent->checkInCodelets14157[this->getID()]);
            myTP->controlTPParent->TP14160Ptr->nextSyncCodeletsblts[this->getID()]
                = &(myTP->controlTPParent->TPParent->barrierCodelets14161[0]);
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14160Ptr->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14160Ptr->availableCodelets[this->getID()] = 1;
#endif
        }
    }
}
TP14153::TP14153(int in_numThreads, int in_mainCodeletID, TP14096* in_TPParent,
    TP14096* in_inputsTPParent, TP14153** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(in_inputsTPParent)
    , ptrToThisTP(in_ptrToThisTP)
    , TP14159Ptr(nullptr)
    , TP14159_alreadyLaunched(0)
    , TP14160Ptr(nullptr)
    , TP14160_alreadyLaunched(0)
    , checkInCodelets14159(new _checkInCodelets14159[this->numThreads])
    , checkInCodelets14160(new _checkInCodelets14160[this->numThreads])
{
    /*Initialize Codelets*/
    _checkInCodelets14160* checkInCodelets14160Ptr = (this->checkInCodelets14160);
    _checkInCodelets14159* checkInCodelets14159Ptr = (this->checkInCodelets14159);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets14159);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets14160Ptr) = _checkInCodelets14160(1, 1, this, codeletCounter);
        checkInCodelets14160Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets14159Ptr) = _checkInCodelets14159(2, 1, this, codeletCounter);
#else
        (*checkInCodelets14159Ptr) = _checkInCodelets14159(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets14159Ptr).decDep();
        checkInCodelets14159Ptr++;
    }
    *(this->ptrToThisTP) = this;
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[this->getID()].decDep();
#else
    this->availableCodelets[this->getID()] = 1;
#endif
}
TP14153::~TP14153()
{
    delete[] checkInCodelets14160;
    delete[] checkInCodelets14159;
}
/*TP14162: ForStmt*/
void TP14162::_checkInCodelets14168::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif

    /*printing node 14168: CallExpr*/
    if (!__sync_val_compare_and_swap(&(myTP->controlTPParent->TP14168_alreadyLaunched), 0, 1)) {
        /*Make the function call*/
        invoke<TP_jacu>(myTP, myTP->numThreads, this->getID(),
            &(myTP->controlTPParent->checkInCodelets14169[this->getID()]),
            &(myTP->controlTPParent->TPParent->barrierCodelets14170[0]),
            &(myTP->controlTPParent->TP14168Ptr),
            (this->inputsTPParent->k_darts14096[this->getID()]));
    } else {
        if (myTP->controlTPParent->TP14168Ptr == nullptr) {
            myTP->add(this);
            return;
        } else {
            myTP->controlTPParent->TP14168Ptr->setNewInputs(
                (this->inputsTPParent->k_darts14096[this->getID()]), this->getID());
            myTP->controlTPParent->TP14168Ptr->nextCodeletsjacu[this->getID()]
                = &(myTP->controlTPParent->checkInCodelets14169[this->getID()]);
            myTP->controlTPParent->TP14168Ptr->nextSyncCodeletsjacu[this->getID()]
                = &(myTP->controlTPParent->TPParent->barrierCodelets14170[0]);
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14168Ptr->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14168Ptr->availableCodelets[this->getID()] = 1;
#endif
        }
    }
}
void TP14162::_checkInCodelets14169::fire(void)
{

    /*printing node 14169: CallExpr*/
    if (!__sync_val_compare_and_swap(&(myTP->controlTPParent->TP14169_alreadyLaunched), 0, 1)) {
        /*Make the function call*/
        invoke<TP_buts>(myTP, myTP->numThreads, this->getID(),
            &(myTP->controlTPParent->TPParent->checkInCodelets14166[this->getID()]),
            &(myTP->controlTPParent->TPParent->barrierCodelets14170[0]),
            &(myTP->controlTPParent->TP14169Ptr), nx, ny, nz,
            (this->inputsTPParent->k_darts14096[this->getID()]), omega, rsd,
            ((this->inputsTPParent->tv_darts14096)), d, a, b, c, ist, iend, jst, jend, nx0, ny0);
    } else {
        if (myTP->controlTPParent->TP14169Ptr == nullptr) {
            myTP->add(this);
            return;
        } else {
            myTP->controlTPParent->TP14169Ptr->setNewInputs(nx, ny, nz,
                (this->inputsTPParent->k_darts14096[this->getID()]), omega, rsd,
                ((this->inputsTPParent->tv_darts14096)), d, a, b, c, ist, iend, jst, jend, nx0, ny0,
                this->getID());
            myTP->controlTPParent->TP14169Ptr->nextCodeletsbuts[this->getID()]
                = &(myTP->controlTPParent->TPParent->checkInCodelets14166[this->getID()]);
            myTP->controlTPParent->TP14169Ptr->nextSyncCodeletsbuts[this->getID()]
                = &(myTP->controlTPParent->TPParent->barrierCodelets14170[0]);
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14169Ptr->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14169Ptr->availableCodelets[this->getID()] = 1;
#endif
        }
    }
}
TP14162::TP14162(int in_numThreads, int in_mainCodeletID, TP14096* in_TPParent,
    TP14096* in_inputsTPParent, TP14162** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(in_inputsTPParent)
    , ptrToThisTP(in_ptrToThisTP)
    , TP14168Ptr(nullptr)
    , TP14168_alreadyLaunched(0)
    , TP14169Ptr(nullptr)
    , TP14169_alreadyLaunched(0)
    , checkInCodelets14168(new _checkInCodelets14168[this->numThreads])
    , checkInCodelets14169(new _checkInCodelets14169[this->numThreads])
{
    /*Initialize Codelets*/
    _checkInCodelets14169* checkInCodelets14169Ptr = (this->checkInCodelets14169);
    _checkInCodelets14168* checkInCodelets14168Ptr = (this->checkInCodelets14168);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets14168);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets14169Ptr) = _checkInCodelets14169(1, 1, this, codeletCounter);
        checkInCodelets14169Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets14168Ptr) = _checkInCodelets14168(2, 1, this, codeletCounter);
#else
        (*checkInCodelets14168Ptr) = _checkInCodelets14168(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets14168Ptr).decDep();
        checkInCodelets14168Ptr++;
    }
    *(this->ptrToThisTP) = this;
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[this->getID()].decDep();
#else
    this->availableCodelets[this->getID()] = 1;
#endif
}
TP14162::~TP14162()
{
    delete[] checkInCodelets14169;
    delete[] checkInCodelets14168;
}
/*TP14171: OMPForDirective*/
void TP14171::_barrierCodelets14171::fire(void)
{
    TP14171* myTP = static_cast<TP14171*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets14171[0].decDep();
}
bool TP14171::requestNewRangeIterations14171(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet14171 * codeletID;
        int tempEndRange = rangePerCodelet14171 * (codeletID + 1);
        if (remainderRange14171 != 0) {
            if (codeletID < (uint32_t)remainderRange14171) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange14171;
                tempEndRange += remainderRange14171;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration14171;
        tempEndRange = tempEndRange * 1 + minIteration14171;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration14171 < lastIteration14171) {
            (this->inputsTPParent->i_darts14171[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts14171[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration14171;
        }
    }
    return isThereNewIteration;
}
void TP14171::_checkInCodelets14172::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 14172: ForStmt*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    /*var: tmp*/
    int* i = &(this->inputsTPParent->i_darts14171[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts14171[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts14171[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts14171[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double* tmp = (this->inputsTPParent->tmp_darts14171);
    (void)tmp /*OMP_SHARED*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations14171(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets14171[0].decDep();
        return;
    }
    for (int i_darts_counter_temp14171 = (*i); i_darts_counter_temp14171 <= endRange
         && i_darts_counter_temp14171 <= this->inputsTPParent->lastIteration14171;
         i_darts_counter_temp14171++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp14171 = (*j);
                for (; j_darts_counter_temp14171 <= jend; j_darts_counter_temp14171++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp14171 = (*k);
                        for (; k_darts_counter_temp14171 <= nz - 2; k_darts_counter_temp14171++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp14171 = (*m);
                                for (; m_darts_counter_temp14171 < 5; m_darts_counter_temp14171++) {
                                    u[(i_darts_counter_temp14171)][j_darts_counter_temp14171]
                                     [k_darts_counter_temp14171][m_darts_counter_temp14171]
                                        = u[(i_darts_counter_temp14171)][j_darts_counter_temp14171]
                                           [k_darts_counter_temp14171][m_darts_counter_temp14171]
                                        + (*(tmp))
                                            * rsd[(i_darts_counter_temp14171)]
                                                 [j_darts_counter_temp14171]
                                                 [k_darts_counter_temp14171]
                                                 [m_darts_counter_temp14171];
                                }
                                (*m) = m_darts_counter_temp14171;
                            }
                        }
                        (*k) = k_darts_counter_temp14171;
                    }
                }
                (*j) = j_darts_counter_temp14171;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets14171[0].decDep();
}
TP14171::TP14171(int in_numThreads, int in_mainCodeletID, TP14096* in_TPParent,
    int in_initIteration, int in_lastIteration, double* in_tmp, TP14171** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts14171(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts14171(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts14171(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts14171(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts14171(in_tmp) /*OMP_SHARED - INPUT*/
    , initIteration14171(in_initIteration)
    , lastIteration14171(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets14171(new _barrierCodelets14171[1])
    , checkInCodelets14172(new _checkInCodelets14172[this->numThreads])
{
    /*Initialize the loop parameters*/
    range14171 = abs(lastIteration14171 - initIteration14171) / 1;
    rangePerCodelet14171 = range14171 / numThreads;
    minIteration14171 = min<int>(lastIteration14171, initIteration14171);
    remainderRange14171 = range14171 % numThreads;
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets14171[0] = _barrierCodelets14171(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets14172* checkInCodelets14172Ptr = (this->checkInCodelets14172);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets14172);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets14172Ptr) = _checkInCodelets14172(2, 1, this, codeletCounter);
#else
        (*checkInCodelets14172Ptr) = _checkInCodelets14172(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets14172Ptr).decDep();
        checkInCodelets14172Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP14171::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets14172[localID].setID(codeletID);
    this->checkInCodelets14172[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets14172[localID + this->baseNumThreads * i]
            = _checkInCodelets14172(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets14172[localID + this->baseNumThreads * i]
            = _checkInCodelets14172(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets14172[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets14172[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP14171::~TP14171()
{
    delete[] barrierCodelets14171;
    delete[] checkInCodelets14172;
}
