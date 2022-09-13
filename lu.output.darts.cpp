#include "lu.output.darts.h"
using namespace darts;
using namespace std;
std::mutex TP9984mutex;
static boolean flag[13];
static void verify(double xcr[5], double xce[5], double xci, char* class_is, boolean* verified);
static void ssor();
static void setcoeff();
static void setbv();
static void pintgr();
static void domain();
static void exact(int i, int j, int k, double u000ijk[5]);
static void erhs();
static void error();
static void setiv();
static void rhs();
static void l2norm(int nx0, int ny0, int nz0, int ist, int iend, int jst, int jend, double sum[5]);
static void read_input();
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
    /*CompoundStmt:175*/
    char class_is;
    boolean verified;
    double mflops;
    int nthreads = 1;
    read_input();
    printf("read_input complete...\n");
    domain();
    printf("domain complete...\n");
    setcoeff();
    printf("setcoeff complete...\n");
    setbv();
    printf("setbv complete...\n");
    setiv();
    printf("setiv complete...\n");
    erhs();
    printf("erhs complete...\n");
    if (affinMaskRes) {
        myDARTSRuntime->run(launch<TP192>(
            ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet, (int*)&((nthreads))));
    }
    ssor();
    printf("ssor complete...\n");
    error();
    printf("error complete...\n");
    pintgr();
    printf("pintgr complete...\n");
    verify(rsdnm, errnm, frc, &class_is, &verified);
    mflops = (double)itmax
        * (1984.77 * (double)nx0 * (double)ny0 * (double)nz0
            - 10923.299999999999
                * (((double)(nx0 + ny0 + nz0) / 3.) * ((double)(nx0 + ny0 + nz0) / 3.))
            + 27770.900000000001 * (double)(nx0 + ny0 + nz0) / 3. - 144010.)
        / (maxtime * 1.0E+6);
    printf("verify complete...\n");
    c_print_results("LU", class_is, nx0, ny0, nz0, itmax, nthreads, maxtime, mflops,
        "          floating point", verified, "3.0 structured", "09 Sep 2022", "gcc", "gcc",
        "-lm -fopenmp", "-I../common", "-O3 -fopenmp", "(none)", "(none)");
}
/*Function: domain, ID: 3*/
static void domain()
{
    /*domain:3*/
    /*CompoundStmt:2178*/
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
    if (nx > 12 || ny > 12 || nz > 12) {
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
    /*CompoundStmt:2198*/
    if (affinMaskRes) {
        myDARTSRuntime->run(
            launch<TP2199>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
    }
}
/*Function: error, ID: 5*/
static void error()
{
    /*error:5*/
    /*CompoundStmt:4725*/
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
    /*CompoundStmt:4789*/
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
static void l2norm(int nx0, int ny0, int nz0, int ist, int iend, int jst, int jend, double sum[5])
{
    /*l2norm:9*/
    /*CompoundStmt:9870*/
    if (affinMaskRes) {
        myDARTSRuntime->run(launch<TP9871>(ompNumThreads * DARTS_CODELETS_MULT, 0,
            RuntimeFinalCodelet, (int*)&((iend)), (int*)&((ist)), (int*)&((jend)), (int*)&((jst)),
            (int*)&((nx0)), (int*)&((ny0)), (int*)&((nz0)), (double**)&((sum))));
    }
}
/*Function: pintgr, ID: 10*/
static void pintgr()
{
    /*pintgr:10*/
    /*CompoundStmt:10013*/
    int i, j, k;
    int ibeg, ifin, ifin1;
    int jbeg, jfin, jfin1;
    int iglob, iglob1, iglob2;
    int jglob, jglob1, jglob2;
    double phi1[14][14];
    double phi2[14][14];
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
    for (i = 0; i <= 12 + 1; i++) {
        for (k = 0; k <= 12 + 1; k++) {
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
    for (i = 0; i <= 12 + 1; i++) {
        for (k = 0; k <= 12 + 1; k++) {
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
    for (i = 0; i <= 12 + 1; i++) {
        for (k = 0; k <= 12 + 1; k++) {
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
    /*CompoundStmt:10643*/
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
        inorm = 50;
        itmax = 50;
        dt = 0.5;
        omega = 1.2;
        tolrsd[0] = 1.0E-8;
        tolrsd[1] = 1.0E-8;
        tolrsd[2] = 1.0E-8;
        tolrsd[3] = 1.0E-8;
        tolrsd[4] = 1.0E-8;
        nx0 = 12;
        ny0 = 12;
        nz0 = 12;
    }
    if (nx0 < 4 || ny0 < 4 || nz0 < 4) {
        printf("     PROBLEM SIZE IS TOO SMALL - \n     SET EACH OF NX, NY AND NZ AT LEAST EQUAL "
               "TO 5\n");
        exit(1);
    }
    if (nx0 > 12 || ny0 > 12 || nz0 > 12) {
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
    /*CompoundStmt:10787*/
    if (affinMaskRes) {
        myDARTSRuntime->run(
            launch<TP10788>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
    }
}
/*Function: setbv, ID: 13*/
static void setbv()
{
    /*setbv:13*/
    /*CompoundStmt:13188*/
    if (affinMaskRes) {
        myDARTSRuntime->run(
            launch<TP13189>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
    }
}
/*Function: setcoeff, ID: 14*/
static void setcoeff()
{
    /*setcoeff:14*/
    /*CompoundStmt:13416*/
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
    /*CompoundStmt:13755*/
    if (affinMaskRes) {
        myDARTSRuntime->run(
            launch<TP13756>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
    }
}
/*Function: ssor, ID: 16*/
static void ssor()
{
    /*ssor:16*/
    /*CompoundStmt:13877*/
    int i, j, k, m;
    int istep;
    double tmp;
    double delunm[5];
    tmp = 1. / (omega * (2. - omega));
    if (affinMaskRes) {
        myDARTSRuntime->run(
            launch<TP13888>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
    }
    rhs();
    l2norm(nx0, ny0, nz0, ist, iend, jst, jend, rsdnm);
    timer_clear(1);
    timer_start(1);
    for (istep = 1; istep <= itmax; istep++) {
        /*CompoundStmt:13981*/
        if (istep % 20 == 0 || istep == itmax || istep == 1) {
            //#pragma omp master
            printf(" Time step %4d\n", istep);
        }
        if (affinMaskRes) {
            myDARTSRuntime->run(launch<TP13987>(
                ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet, (double*)&((tmp))));
        }
        if (istep % inorm == 0) {
            l2norm(nx0, ny0, nz0, ist, iend, jst, jend, delunm);
        }
        rhs();
        if ((istep % inorm == 0) || (istep == itmax)) {
            l2norm(nx0, ny0, nz0, ist, iend, jst, jend, rsdnm);
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
static void verify(double xcr[5], double xce[5], double xci, char* class_is, boolean* verified)
{
    /*verify:17*/
    /*CompoundStmt:14138*/
    double xcrref[5], xceref[5], xciref, xcrdif[5], xcedif[5], xcidif, epsilon, dtref;
    int m;
    epsilon = 1.0E-8;
    *class_is = 'U';
    *verified = 1;
    for (m = 0; m < 5; m++) {
        xcrref[m] = 1.;
        xceref[m] = 1.;
    }
    xciref = 1.;
    if (nx0 == 12 && ny0 == 12 && nz0 == 12 && itmax == 50) {
        *class_is = 'S';
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
        *class_is = 'W';
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
        *class_is = 'A';
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
        *class_is = 'B';
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
        *class_is = 'C';
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
    if (*class_is != 'U') {
        printf("\n Verification being performed for class_is %1c\n", *class_is);
        printf(" Accuracy setting for epsilon = %20.13e\n", epsilon);
        if (fabs(dt - dtref) > epsilon) {
            *verified = 0;
            *class_is = 'U';
            printf(" DT does not match the reference value of %15.8e\n", dtref);
        }
    } else {
        printf(" Unknown class_is\n");
    }
    if (*class_is != 'U') {
        printf(" Comparison of RMS-norms of residual\n");
    } else {
        printf(" RMS-norms of residual\n");
    }
    for (m = 0; m < 5; m++) {
        if (*class_is == 'U') {
            printf("          %2d  %20.13e\n", m, xcr[m]);
        } else if (xcrdif[m] > epsilon) {
            *verified = 0;
            printf(" FAILURE: %2d  %20.13e%20.13e%20.13e\n", m, xcr[m], xcrref[m], xcrdif[m]);
        } else {
            printf("          %2d  %20.13e%20.13e%20.13e\n", m, xcr[m], xcrref[m], xcrdif[m]);
        }
    }
    if (*class_is != 'U') {
        printf(" Comparison of RMS-norms of solution error\n");
    } else {
        printf(" RMS-norms of solution error\n");
    }
    for (m = 0; m < 5; m++) {
        if (*class_is == 'U') {
            printf("          %2d  %20.13e\n", m, xce[m]);
        } else if (xcedif[m] > epsilon) {
            *verified = 0;
            printf(" FAILURE: %2d  %20.13e%20.13e%20.13e\n", m, xce[m], xceref[m], xcedif[m]);
        } else {
            printf("          %2d  %20.13e%20.13e%20.13e\n", m, xce[m], xceref[m], xcedif[m]);
        }
    }
    if (*class_is != 'U') {
        printf(" Comparison of surface integral\n");
    } else {
        printf(" Surface integral\n");
    }
    if (*class_is == 'U') {
        printf("              %20.13e\n", xci);
    } else if (xcidif > epsilon) {
        *verified = 0;
        printf(" FAILURE:     %20.13e%20.13e%20.13e\n", xci, xciref, xcidif);
    } else {
        printf("              %20.13e%20.13e%20.13e\n", xci, xciref, xcidif);
    }
    if (*class_is == 'U') {
        printf(" No reference values provided\n");
        printf(" No verification performed\n");
    } else if (*verified) {
        printf(" Verification Successful\n");
    } else {
        printf(" Verification failed\n");
    }
}
/*TP192: OMPParallelDirective*/
void TP192::_barrierCodelets192::fire(void)
{
    TP192* myTP = static_cast<TP192*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP192::_checkInCodelets194::fire(void)
{
    if (this->getID() == 0) {
        /*Init the vars for this region*/
        /*Initialize the vars of the inlined region*/
        this->inputsTPParent->nthreads_darts194
            = (this->inputsTPParent->nthreads_darts192) /*OMP_SHARED - VAR INLINED*/;

        /*printing node 195: BinaryOperator*/
        (*(this->inputsTPParent->nthreads_darts194)) = omp_get_num_threads();
        /*Signaling next codelet from last stmt in the codelet*/
        /*Find and signal the next codelet*/
        myTP->controlTPParent->TPParent->barrierCodelets192[0].decDep();
    } else {
        /*Find and signal the next codelet*/
        myTP->TPParent->barrierCodelets192[0].decDep();
    }
}
TP192::TP192(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet, int* in_nthreads)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , nthreads_darts192(in_nthreads) /*OMP_SHARED - INPUT*/
    , TP194_alreadyLaunched(0)
    , barrierCodelets192(new _barrierCodelets192[1])
    , checkInCodelets194(new _checkInCodelets194[this->numThreads])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets192[0] = _barrierCodelets192(ompNumThreads, ompNumThreads, this, 0);
    _checkInCodelets194* checkInCodelets194Ptr = (this->checkInCodelets194);
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets194Ptr) = _checkInCodelets194(1, 1, this, codeletCounter);
        (*checkInCodelets194Ptr).decDep();
        checkInCodelets194Ptr++;
    }
}
TP192::~TP192()
{
    delete[] barrierCodelets192;
    delete[] checkInCodelets194;
}
/*TP1: TP_blts*/
void TP1::_checkInCodelets248::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*region 248 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP248;
    if (idx < myTP->TPsToUse248) {
        if (!__sync_val_compare_and_swap(&(myTP->TP248_alreadyLaunched[idx]), 0, 1)) {
            int range = abs((this->inputsTPParent->iend_darts1[this->getID()])
                            - (this->inputsTPParent->ist_darts1[this->getID()]))
                / 1;
            int rangePerCodelet = range / myTP->TPsToUse248;
            int minIteration = min<int>((this->inputsTPParent->iend_darts1[this->getID()]),
                (this->inputsTPParent->ist_darts1[this->getID()]));
            int remainderRange = range % myTP->TPsToUse248;
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
            if ((this->inputsTPParent->ist_darts1[this->getID()])
                < (this->inputsTPParent->iend_darts1[this->getID()])) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse248 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse248 - 1) {
                lastIteration = (this->inputsTPParent->iend_darts1[this->getID()]);
            }
#if USEINVOKE == 1
            invoke<TP248>(myTP, myTP->codeletsPerTP248 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP248Ptr[idx]));
#else
            place<TP248>(idx, myTP, myTP->codeletsPerTP248 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP248Ptr[idx]));
#endif
        } else {
            if (myTP->TP248Ptr[idx] != nullptr) {
                myTP->TP248Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP1::_barrierCodelets248::fire(void)
{
    TP1* myTP = static_cast<TP1*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets352[codeletsCounter].decDep();
        }
    }
}
void TP1::_checkInCodelets352::fire(void)
{
    if (this->getID() == 0) {
        /*Init the vars for this region*/
        /*Initialize the vars of the inlined region*/
        this->inputsTPParent->iend_darts352
            = &(this->inputsTPParent
                    ->iend_darts1[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;
        this->inputsTPParent->ist_darts352
            = &(this->inputsTPParent
                    ->ist_darts1[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;
        this->inputsTPParent->jend_darts352
            = &(this->inputsTPParent
                    ->jend_darts1[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;
        this->inputsTPParent->jst_darts352
            = &(this->inputsTPParent
                    ->jst_darts1[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;
        this->inputsTPParent->k_darts352
            = &(this->inputsTPParent
                    ->k_darts1[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;
        this->inputsTPParent->omega_darts352
            = &(this->inputsTPParent
                    ->omega_darts1[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;
        this->inputsTPParent->tmp_darts352
            = &(this->inputsTPParent
                    ->tmp_darts1[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;
        this->inputsTPParent->tmp1_darts352
            = &(this->inputsTPParent
                    ->tmp1_darts1[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;

        /*printing node 353: ForStmt*/
        {
            /*Loop's init*/
            (this->inputsTPParent->i_darts352) = (*(this->inputsTPParent->ist_darts352));
            int i_darts_counter_temp352 = (this->inputsTPParent->i_darts352);
            for (; i_darts_counter_temp352 <= (*(this->inputsTPParent->iend_darts352));
                 i_darts_counter_temp352++) {
                if (i_darts_counter_temp352 != (*(this->inputsTPParent->ist_darts352))) {
                    while (flag[i_darts_counter_temp352 - 1] == 0) { }
                }
                if (i_darts_counter_temp352 != (*(this->inputsTPParent->iend_darts352))) {
                    while (flag[i_darts_counter_temp352] == 1) { }
                }
                {
                    /*Loop's init*/
                    (this->inputsTPParent->j_darts352) = (*(this->inputsTPParent->jst_darts352));
                    int j_darts_counter_temp352 = (this->inputsTPParent->j_darts352);
                    for (; j_darts_counter_temp352 <= (*(this->inputsTPParent->jend_darts352));
                         j_darts_counter_temp352++) {
                        {
                            /*Loop's init*/
                            (this->inputsTPParent->m_darts352) = 0;
                            int m_darts_counter_temp352 = (this->inputsTPParent->m_darts352);
                            for (; m_darts_counter_temp352 < 5; m_darts_counter_temp352++) {
                                rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                   [(*(this->inputsTPParent->k_darts352))][m_darts_counter_temp352]
                                    = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                         [(*(this->inputsTPParent->k_darts352))]
                                         [m_darts_counter_temp352]
                                    - (*(this->inputsTPParent->omega_darts352))
                                        * (b[i_darts_counter_temp352][j_darts_counter_temp352]
                                            [m_darts_counter_temp352][0]
                                                * rsd[(*(this->inputsTPParent->k_darts352))]
                                                     [i_darts_counter_temp352]
                                                     [j_darts_counter_temp352 - 1][0]
                                            + c[i_darts_counter_temp352][j_darts_counter_temp352]
                                               [m_darts_counter_temp352][0]
                                                * rsd[(*(this->inputsTPParent->k_darts352))]
                                                     [i_darts_counter_temp352 - 1]
                                                     [j_darts_counter_temp352][0]
                                            + b[i_darts_counter_temp352][j_darts_counter_temp352]
                                               [m_darts_counter_temp352][1]
                                                * rsd[(*(this->inputsTPParent->k_darts352))]
                                                     [i_darts_counter_temp352]
                                                     [j_darts_counter_temp352 - 1][1]
                                            + c[i_darts_counter_temp352][j_darts_counter_temp352]
                                               [m_darts_counter_temp352][1]
                                                * rsd[(*(this->inputsTPParent->k_darts352))]
                                                     [i_darts_counter_temp352 - 1]
                                                     [j_darts_counter_temp352][1]
                                            + b[i_darts_counter_temp352][j_darts_counter_temp352]
                                               [m_darts_counter_temp352][2]
                                                * rsd[(*(this->inputsTPParent->k_darts352))]
                                                     [i_darts_counter_temp352]
                                                     [j_darts_counter_temp352 - 1][2]
                                            + c[i_darts_counter_temp352][j_darts_counter_temp352]
                                               [m_darts_counter_temp352][2]
                                                * rsd[(*(this->inputsTPParent->k_darts352))]
                                                     [i_darts_counter_temp352 - 1]
                                                     [j_darts_counter_temp352][2]
                                            + b[i_darts_counter_temp352][j_darts_counter_temp352]
                                               [m_darts_counter_temp352][3]
                                                * rsd[(*(this->inputsTPParent->k_darts352))]
                                                     [i_darts_counter_temp352]
                                                     [j_darts_counter_temp352 - 1][3]
                                            + c[i_darts_counter_temp352][j_darts_counter_temp352]
                                               [m_darts_counter_temp352][3]
                                                * rsd[(*(this->inputsTPParent->k_darts352))]
                                                     [i_darts_counter_temp352 - 1]
                                                     [j_darts_counter_temp352][3]
                                            + b[i_darts_counter_temp352][j_darts_counter_temp352]
                                               [m_darts_counter_temp352][4]
                                                * rsd[(*(this->inputsTPParent->k_darts352))]
                                                     [i_darts_counter_temp352]
                                                     [j_darts_counter_temp352 - 1][4]
                                            + c[i_darts_counter_temp352][j_darts_counter_temp352]
                                               [m_darts_counter_temp352][4]
                                                * rsd[(*(this->inputsTPParent->k_darts352))]
                                                     [i_darts_counter_temp352 - 1]
                                                     [j_darts_counter_temp352][4]);
                            }
                            (this->inputsTPParent->m_darts352) = m_darts_counter_temp352;
                        }
                        {
                            /*Loop's init*/
                            (this->inputsTPParent->m_darts352) = 0;
                            int m_darts_counter_temp352 = (this->inputsTPParent->m_darts352);
                            for (; m_darts_counter_temp352 < 5; m_darts_counter_temp352++) {
                                tmat[m_darts_counter_temp352][0]
                                    = d[i_darts_counter_temp352][j_darts_counter_temp352]
                                       [m_darts_counter_temp352][0];
                                tmat[m_darts_counter_temp352][1]
                                    = d[i_darts_counter_temp352][j_darts_counter_temp352]
                                       [m_darts_counter_temp352][1];
                                tmat[m_darts_counter_temp352][2]
                                    = d[i_darts_counter_temp352][j_darts_counter_temp352]
                                       [m_darts_counter_temp352][2];
                                tmat[m_darts_counter_temp352][3]
                                    = d[i_darts_counter_temp352][j_darts_counter_temp352]
                                       [m_darts_counter_temp352][3];
                                tmat[m_darts_counter_temp352][4]
                                    = d[i_darts_counter_temp352][j_darts_counter_temp352]
                                       [m_darts_counter_temp352][4];
                            }
                            (this->inputsTPParent->m_darts352) = m_darts_counter_temp352;
                        }
                        (*(this->inputsTPParent->tmp1_darts352)) = 1. / tmat[0][0];
                        (*(this->inputsTPParent->tmp_darts352))
                            = (*(this->inputsTPParent->tmp1_darts352)) * tmat[1][0];
                        tmat[1][1]
                            = tmat[1][1] - (*(this->inputsTPParent->tmp_darts352)) * tmat[0][1];
                        tmat[1][2]
                            = tmat[1][2] - (*(this->inputsTPParent->tmp_darts352)) * tmat[0][2];
                        tmat[1][3]
                            = tmat[1][3] - (*(this->inputsTPParent->tmp_darts352)) * tmat[0][3];
                        tmat[1][4]
                            = tmat[1][4] - (*(this->inputsTPParent->tmp_darts352)) * tmat[0][4];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][1]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][1]
                            - rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][0]
                                * (*(this->inputsTPParent->tmp_darts352));
                        (*(this->inputsTPParent->tmp_darts352))
                            = (*(this->inputsTPParent->tmp1_darts352)) * tmat[2][0];
                        tmat[2][1]
                            = tmat[2][1] - (*(this->inputsTPParent->tmp_darts352)) * tmat[0][1];
                        tmat[2][2]
                            = tmat[2][2] - (*(this->inputsTPParent->tmp_darts352)) * tmat[0][2];
                        tmat[2][3]
                            = tmat[2][3] - (*(this->inputsTPParent->tmp_darts352)) * tmat[0][3];
                        tmat[2][4]
                            = tmat[2][4] - (*(this->inputsTPParent->tmp_darts352)) * tmat[0][4];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][2]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][2]
                            - rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][0]
                                * (*(this->inputsTPParent->tmp_darts352));
                        (*(this->inputsTPParent->tmp_darts352))
                            = (*(this->inputsTPParent->tmp1_darts352)) * tmat[3][0];
                        tmat[3][1]
                            = tmat[3][1] - (*(this->inputsTPParent->tmp_darts352)) * tmat[0][1];
                        tmat[3][2]
                            = tmat[3][2] - (*(this->inputsTPParent->tmp_darts352)) * tmat[0][2];
                        tmat[3][3]
                            = tmat[3][3] - (*(this->inputsTPParent->tmp_darts352)) * tmat[0][3];
                        tmat[3][4]
                            = tmat[3][4] - (*(this->inputsTPParent->tmp_darts352)) * tmat[0][4];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][3]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][3]
                            - rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][0]
                                * (*(this->inputsTPParent->tmp_darts352));
                        (*(this->inputsTPParent->tmp_darts352))
                            = (*(this->inputsTPParent->tmp1_darts352)) * tmat[4][0];
                        tmat[4][1]
                            = tmat[4][1] - (*(this->inputsTPParent->tmp_darts352)) * tmat[0][1];
                        tmat[4][2]
                            = tmat[4][2] - (*(this->inputsTPParent->tmp_darts352)) * tmat[0][2];
                        tmat[4][3]
                            = tmat[4][3] - (*(this->inputsTPParent->tmp_darts352)) * tmat[0][3];
                        tmat[4][4]
                            = tmat[4][4] - (*(this->inputsTPParent->tmp_darts352)) * tmat[0][4];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][4]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][4]
                            - rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][0]
                                * (*(this->inputsTPParent->tmp_darts352));
                        (*(this->inputsTPParent->tmp1_darts352)) = 1. / tmat[1][1];
                        (*(this->inputsTPParent->tmp_darts352))
                            = (*(this->inputsTPParent->tmp1_darts352)) * tmat[2][1];
                        tmat[2][2]
                            = tmat[2][2] - (*(this->inputsTPParent->tmp_darts352)) * tmat[1][2];
                        tmat[2][3]
                            = tmat[2][3] - (*(this->inputsTPParent->tmp_darts352)) * tmat[1][3];
                        tmat[2][4]
                            = tmat[2][4] - (*(this->inputsTPParent->tmp_darts352)) * tmat[1][4];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][2]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][2]
                            - rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][1]
                                * (*(this->inputsTPParent->tmp_darts352));
                        (*(this->inputsTPParent->tmp_darts352))
                            = (*(this->inputsTPParent->tmp1_darts352)) * tmat[3][1];
                        tmat[3][2]
                            = tmat[3][2] - (*(this->inputsTPParent->tmp_darts352)) * tmat[1][2];
                        tmat[3][3]
                            = tmat[3][3] - (*(this->inputsTPParent->tmp_darts352)) * tmat[1][3];
                        tmat[3][4]
                            = tmat[3][4] - (*(this->inputsTPParent->tmp_darts352)) * tmat[1][4];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][3]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][3]
                            - rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][1]
                                * (*(this->inputsTPParent->tmp_darts352));
                        (*(this->inputsTPParent->tmp_darts352))
                            = (*(this->inputsTPParent->tmp1_darts352)) * tmat[4][1];
                        tmat[4][2]
                            = tmat[4][2] - (*(this->inputsTPParent->tmp_darts352)) * tmat[1][2];
                        tmat[4][3]
                            = tmat[4][3] - (*(this->inputsTPParent->tmp_darts352)) * tmat[1][3];
                        tmat[4][4]
                            = tmat[4][4] - (*(this->inputsTPParent->tmp_darts352)) * tmat[1][4];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][4]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][4]
                            - rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][1]
                                * (*(this->inputsTPParent->tmp_darts352));
                        (*(this->inputsTPParent->tmp1_darts352)) = 1. / tmat[2][2];
                        (*(this->inputsTPParent->tmp_darts352))
                            = (*(this->inputsTPParent->tmp1_darts352)) * tmat[3][2];
                        tmat[3][3]
                            = tmat[3][3] - (*(this->inputsTPParent->tmp_darts352)) * tmat[2][3];
                        tmat[3][4]
                            = tmat[3][4] - (*(this->inputsTPParent->tmp_darts352)) * tmat[2][4];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][3]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][3]
                            - rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][2]
                                * (*(this->inputsTPParent->tmp_darts352));
                        (*(this->inputsTPParent->tmp_darts352))
                            = (*(this->inputsTPParent->tmp1_darts352)) * tmat[4][2];
                        tmat[4][3]
                            = tmat[4][3] - (*(this->inputsTPParent->tmp_darts352)) * tmat[2][3];
                        tmat[4][4]
                            = tmat[4][4] - (*(this->inputsTPParent->tmp_darts352)) * tmat[2][4];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][4]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][4]
                            - rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][2]
                                * (*(this->inputsTPParent->tmp_darts352));
                        (*(this->inputsTPParent->tmp1_darts352)) = 1. / tmat[3][3];
                        (*(this->inputsTPParent->tmp_darts352))
                            = (*(this->inputsTPParent->tmp1_darts352)) * tmat[4][3];
                        tmat[4][4]
                            = tmat[4][4] - (*(this->inputsTPParent->tmp_darts352)) * tmat[3][4];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][4]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][4]
                            - rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][3]
                                * (*(this->inputsTPParent->tmp_darts352));
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][4]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][4]
                            / tmat[4][4];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][3]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][3]
                            - tmat[3][4]
                                * rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                     [(*(this->inputsTPParent->k_darts352))][4];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][3]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][3]
                            / tmat[3][3];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][2]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][2]
                            - tmat[2][3]
                                * rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                     [(*(this->inputsTPParent->k_darts352))][3]
                            - tmat[2][4]
                                * rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                     [(*(this->inputsTPParent->k_darts352))][4];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][2]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][2]
                            / tmat[2][2];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][1]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][1]
                            - tmat[1][2]
                                * rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                     [(*(this->inputsTPParent->k_darts352))][2]
                            - tmat[1][3]
                                * rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                     [(*(this->inputsTPParent->k_darts352))][3]
                            - tmat[1][4]
                                * rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                     [(*(this->inputsTPParent->k_darts352))][4];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][1]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][1]
                            / tmat[1][1];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][0]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][0]
                            - tmat[0][1]
                                * rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                     [(*(this->inputsTPParent->k_darts352))][1]
                            - tmat[0][2]
                                * rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                     [(*(this->inputsTPParent->k_darts352))][2]
                            - tmat[0][3]
                                * rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                     [(*(this->inputsTPParent->k_darts352))][3]
                            - tmat[0][4]
                                * rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                     [(*(this->inputsTPParent->k_darts352))][4];
                        rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                           [(*(this->inputsTPParent->k_darts352))][0]
                            = rsd[i_darts_counter_temp352][j_darts_counter_temp352]
                                 [(*(this->inputsTPParent->k_darts352))][0]
                            / tmat[0][0];
                    }
                    (this->inputsTPParent->j_darts352) = j_darts_counter_temp352;
                }
                if (i_darts_counter_temp352 != (*(this->inputsTPParent->ist_darts352))) {
                    flag[i_darts_counter_temp352 - 1] = 0;
                }
                if (i_darts_counter_temp352 != (*(this->inputsTPParent->iend_darts352))) {
                    flag[i_darts_counter_temp352] = 1;
                }
            }
            (this->inputsTPParent->i_darts352) = i_darts_counter_temp352;
        }
        /*Signaling next codelet from last stmt in the codelet*/
        /*Find and signal the next codelet*/

        myTP->controlTPParent->nextCodeletsblts[this->getID()]->decDep();
    } else {
        /*Find and signal the next codelet*/

        myTP->nextCodeletsblts[this->getID()]->decDep();
    }
}
TP1::TP1(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
    darts::Codelet* in_mainSyncCodelet, TP1** in_ptrToThisFunctionTP, int in_nx, int in_ny,
    int in_nz, int in_k, double in_omega, int in_ist, int in_iend, int in_jst, int in_jend,
    int in_nx0, int in_ny0)
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
    , ist_darts1(new int[this->numThreads])
    , iend_darts1(new int[this->numThreads])
    , jst_darts1(new int[this->numThreads])
    , jend_darts1(new int[this->numThreads])
    , nx0_darts1(new int[this->numThreads])
    , ny0_darts1(new int[this->numThreads])
    , i_darts1(new int[this->numThreads])
    , j_darts1(new int[this->numThreads])
    , m_darts1(new int[this->numThreads])
    , tmp_darts1(new double[this->numThreads])
    , tmp1_darts1(new double[this->numThreads])
    , TP248Ptr(new TP248*[NUMTPS248])
    , TP248_alreadyLaunched(new size_t[NUMTPS248])
    , numTPsSet248(0)
    , numTPsReady248(0)
    , TPsToUse248(NUMTPS248)
    , codeletsPerTP248(this->numThreads / NUMTPS248)
    , totalCodelets248(this->TPsToUse248 * this->codeletsPerTP248)
    , TP352_alreadyLaunched(0)
    , checkInCodelets248(new _checkInCodelets248[this->numThreads])
    , barrierCodelets248(new _barrierCodelets248[1])
    , checkInCodelets352(new _checkInCodelets352[this->numThreads])
{
    barrierCodelets248[0] = _barrierCodelets248(NUMTPS248, NUMTPS248, this, 0);
    _checkInCodelets352* checkInCodelets352Ptr = (this->checkInCodelets352);
    _checkInCodelets248* checkInCodelets248Ptr = (this->checkInCodelets248);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets248);
#endif
    for (int i = 0; i < NUMTPS248; i++) {
        TP248Ptr[i] = nullptr;
        TP248_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++) {
        (*checkInCodelets352Ptr) = _checkInCodelets352(1, 1, this, codeletCounter);
        checkInCodelets352Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets248Ptr) = _checkInCodelets248(2, 1, this, codeletCounter);
#else
        (*checkInCodelets248Ptr) = _checkInCodelets248(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets248Ptr).decDep();
        checkInCodelets248Ptr++;
    }
    if (this->numThreads == 1) {
        this->nextCodeletsblts[0] = in_mainNextCodelet;
        this->nextSyncCodeletsblts[0] = in_mainSyncCodelet;
        this->nx_darts1[0] = in_nx;
        this->ny_darts1[0] = in_ny;
        this->nz_darts1[0] = in_nz;
        this->k_darts1[0] = in_k;
        this->omega_darts1[0] = in_omega;
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
    delete[] checkInCodelets352;
    delete[] barrierCodelets248;
    delete[] checkInCodelets248;
    delete[] nextCodeletsblts;
    delete[] nextSyncCodeletsblts;
    delete[] nx_darts1;
    delete[] ny_darts1;
    delete[] nz_darts1;
    delete[] k_darts1;
    delete[] omega_darts1;
    delete[] ist_darts1;
    delete[] iend_darts1;
    delete[] jst_darts1;
    delete[] jend_darts1;
    delete[] nx0_darts1;
    delete[] ny0_darts1;
    delete[] i_darts1;
    delete[] j_darts1;
    delete[] m_darts1;
    delete[] tmp_darts1;
    delete[] tmp1_darts1;
}
void TP1::setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, int in_ist,
    int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID)
{
    this->nx_darts1[codeletID] = in_nx;
    this->ny_darts1[codeletID] = in_ny;
    this->nz_darts1[codeletID] = in_nz;
    this->k_darts1[codeletID] = in_k;
    this->omega_darts1[codeletID] = in_omega;
    this->ist_darts1[codeletID] = in_ist;
    this->iend_darts1[codeletID] = in_iend;
    this->jst_darts1[codeletID] = in_jst;
    this->jend_darts1[codeletID] = in_jend;
    this->nx0_darts1[codeletID] = in_nx0;
    this->ny0_darts1[codeletID] = in_ny0;
}
/*TP248: OMPForDirective*/
void TP248::_barrierCodelets248::fire(void)
{
    TP248* myTP = static_cast<TP248*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets248[0].decDep();
}
bool TP248::requestNewRangeIterations248(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 1*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet248 * codeletID;
        int tempEndRange = rangePerCodelet248 * (codeletID + 1);
        if (remainderRange248 != 0) {
            if (codeletID < (uint32_t)remainderRange248) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange248;
                tempEndRange += remainderRange248;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration248;
        tempEndRange = tempEndRange * 1 + minIteration248;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration248 < lastIteration248) {
            (this->inputsTPParent->i_darts248[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts248[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration248;
        }
    }
    return isThereNewIteration;
}
void TP248::_checkInCodelets249::fire(void)
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
    this->inputsTPParent->iend_darts248[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iend_darts1[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ist_darts248[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->ist_darts1[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jend_darts248[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jend_darts1[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jst_darts248[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jst_darts1[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->k_darts248[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->k_darts1[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->omega_darts248[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->omega_darts1[this->getID()]);

    /*printing node 249: ForStmt*/
    /*var: i*/
    /*var: iend*/
    /*var: ist*/
    /*var: j*/
    /*var: jend*/
    /*var: jst*/
    /*var: k*/
    /*var: m*/
    /*var: omega*/
    int* i = &(this->inputsTPParent->i_darts248[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts248[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jend = &(this->inputsTPParent->jend_darts248[this->getLocalID()]);
    (void)jend /*OMP_SHARED_PRIVATE*/;
    int** jst = &(this->inputsTPParent->jst_darts248[this->getLocalID()]);
    (void)jst /*OMP_SHARED_PRIVATE*/;
    int** k = &(this->inputsTPParent->k_darts248[this->getLocalID()]);
    (void)k /*OMP_SHARED_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts248[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** omega = &(this->inputsTPParent->omega_darts248[this->getLocalID()]);
    (void)omega /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations248(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets248[0].decDep();
        return;
    }
    for (int i_darts_counter_temp248 = (*i); i_darts_counter_temp248 <= endRange
         && i_darts_counter_temp248 <= this->inputsTPParent->lastIteration248;
         i_darts_counter_temp248++) {
        {
            {
                /*Loop's init*/
                (*j) = (*(*jst));
                int j_darts_counter_temp248 = (*j);
                for (; j_darts_counter_temp248 <= (*(*jend)); j_darts_counter_temp248++) {
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp248 = (*m);
                        for (; m_darts_counter_temp248 < 5; m_darts_counter_temp248++) {
                            rsd[(i_darts_counter_temp248)][j_darts_counter_temp248][(*(*k))]
                               [m_darts_counter_temp248]
                                = rsd[(i_darts_counter_temp248)][j_darts_counter_temp248][(*(*k))]
                                     [m_darts_counter_temp248]
                                - (*(*omega))
                                    * (a[(i_darts_counter_temp248)][j_darts_counter_temp248]
                                        [m_darts_counter_temp248][0]
                                            * rsd[(*(*k)) - 1][(i_darts_counter_temp248)]
                                                 [j_darts_counter_temp248][0]
                                        + a[(i_darts_counter_temp248)][j_darts_counter_temp248]
                                           [m_darts_counter_temp248][1]
                                            * rsd[(*(*k)) - 1][(i_darts_counter_temp248)]
                                                 [j_darts_counter_temp248][1]
                                        + a[(i_darts_counter_temp248)][j_darts_counter_temp248]
                                           [m_darts_counter_temp248][2]
                                            * rsd[(*(*k)) - 1][(i_darts_counter_temp248)]
                                                 [j_darts_counter_temp248][2]
                                        + a[(i_darts_counter_temp248)][j_darts_counter_temp248]
                                           [m_darts_counter_temp248][3]
                                            * rsd[(*(*k)) - 1][(i_darts_counter_temp248)]
                                                 [j_darts_counter_temp248][3]
                                        + a[(i_darts_counter_temp248)][j_darts_counter_temp248]
                                           [m_darts_counter_temp248][4]
                                            * rsd[(*(*k)) - 1][(i_darts_counter_temp248)]
                                                 [j_darts_counter_temp248][4]);
                        }
                        (*m) = m_darts_counter_temp248;
                    }
                }
                (*j) = j_darts_counter_temp248;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets248[0].decDep();
}
TP248::TP248(int in_numThreads, int in_mainCodeletID, TP1* in_TPParent, int in_initIteration,
    int in_lastIteration, TP248** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts248(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iend_darts248(new int*[this->numThreads])
    , ist_darts248(new int*[this->numThreads])
    , j_darts248(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jend_darts248(new int*[this->numThreads])
    , jst_darts248(new int*[this->numThreads])
    , k_darts248(new int*[this->numThreads])
    , m_darts248(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , omega_darts248(new double*[this->numThreads])
    , initIteration248(in_initIteration)
    , lastIteration248(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets248(new _barrierCodelets248[1])
    , checkInCodelets249(new _checkInCodelets249[this->numThreads])
{
    /*Initialize the loop parameters*/
    range248 = abs(lastIteration248 - initIteration248) / 1;
    rangePerCodelet248 = range248 / numThreads;
    minIteration248 = min<int>(lastIteration248, initIteration248);
    remainderRange248 = range248 % numThreads;
    /*Initialize inputs and vars.*/
    this->iend_darts248 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ist_darts248 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jend_darts248 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jst_darts248 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->k_darts248 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->omega_darts248
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets248[0] = _barrierCodelets248(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets249* checkInCodelets249Ptr = (this->checkInCodelets249);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets249);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets249Ptr) = _checkInCodelets249(2, 1, this, codeletCounter);
#else
        (*checkInCodelets249Ptr) = _checkInCodelets249(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets249Ptr).decDep();
        checkInCodelets249Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP248::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets249[localID].setID(codeletID);
    this->checkInCodelets249[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets249[localID + this->baseNumThreads * i]
            = _checkInCodelets249(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets249[localID + this->baseNumThreads * i]
            = _checkInCodelets249(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets249[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets249[localID + this->baseNumThreads * i].decDep();
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
TP248::~TP248()
{
    delete[] iend_darts248;
    delete[] ist_darts248;
    delete[] jend_darts248;
    delete[] jst_darts248;
    delete[] k_darts248;
    delete[] omega_darts248;
    delete[] barrierCodelets248;
    delete[] checkInCodelets249;
}
/*TP2: TP_buts*/
void TP2::_checkInCodelets1213::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*region 1213 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP1213;
    if (idx < myTP->TPsToUse1213) {
        if (!__sync_val_compare_and_swap(&(myTP->TP1213_alreadyLaunched[idx]), 0, 1)) {
            int range = abs((this->inputsTPParent->ist_darts2[this->getID()])
                            - (this->inputsTPParent->iend_darts2[this->getID()]))
                / 1;
            int rangePerCodelet = range / myTP->TPsToUse1213;
            int minIteration = min<int>((this->inputsTPParent->ist_darts2[this->getID()]),
                (this->inputsTPParent->iend_darts2[this->getID()]));
            int remainderRange = range % myTP->TPsToUse1213;
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
            if ((this->inputsTPParent->iend_darts2[this->getID()])
                < (this->inputsTPParent->ist_darts2[this->getID()])) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == 0) {
                lastIteration = lastIteration - 1;
            }
            if (idx == myTP->TPsToUse1213 - 1) {
                lastIteration = (this->inputsTPParent->ist_darts2[this->getID()]);
            }
#if USEINVOKE == 1
            invoke<TP1213>(myTP, myTP->codeletsPerTP1213 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP1213Ptr[idx]));
#else
            place<TP1213>(idx, myTP, myTP->codeletsPerTP1213 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP1213Ptr[idx]));
#endif
        } else {
            if (myTP->TP1213Ptr[idx] != nullptr) {
                myTP->TP1213Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2::_barrierCodelets1213::fire(void)
{
    TP2* myTP = static_cast<TP2*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets1315[codeletsCounter].decDep();
        }
    }
}
void TP2::_checkInCodelets1315::fire(void)
{
    if (this->getID() == 0) {
        /*Init the vars for this region*/
        /*Initialize the vars of the inlined region*/
        this->inputsTPParent->iend_darts1315
            = &(this->inputsTPParent
                    ->iend_darts2[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;
        this->inputsTPParent->ist_darts1315
            = &(this->inputsTPParent
                    ->ist_darts2[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;
        this->inputsTPParent->jend_darts1315
            = &(this->inputsTPParent
                    ->jend_darts2[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;
        this->inputsTPParent->jst_darts1315
            = &(this->inputsTPParent
                    ->jst_darts2[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;
        this->inputsTPParent->k_darts1315
            = &(this->inputsTPParent
                    ->k_darts2[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;
        this->inputsTPParent->omega_darts1315
            = &(this->inputsTPParent
                    ->omega_darts2[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;
        this->inputsTPParent->tmp_darts1315
            = &(this->inputsTPParent
                    ->tmp_darts2[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;
        this->inputsTPParent->tmp1_darts1315
            = &(this->inputsTPParent
                    ->tmp1_darts2[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;

        /*printing node 1316: ForStmt*/
        {
            /*Loop's init*/
            (this->inputsTPParent->i_darts1315) = (*(this->inputsTPParent->iend_darts1315));
            int i_darts_counter_temp1315 = (this->inputsTPParent->i_darts1315);
            for (; i_darts_counter_temp1315 >= (*(this->inputsTPParent->ist_darts1315));
                 i_darts_counter_temp1315--) {
                if (i_darts_counter_temp1315 != (*(this->inputsTPParent->iend_darts1315))) {
                    while (flag[i_darts_counter_temp1315 + 1] == 0) { }
                }
                if (i_darts_counter_temp1315 != (*(this->inputsTPParent->ist_darts1315))) {
                    while (flag[i_darts_counter_temp1315] == 1) { }
                }
                {
                    /*Loop's init*/
                    (this->inputsTPParent->j_darts1315) = (*(this->inputsTPParent->jend_darts1315));
                    int j_darts_counter_temp1315 = (this->inputsTPParent->j_darts1315);
                    for (; j_darts_counter_temp1315 >= (*(this->inputsTPParent->jst_darts1315));
                         j_darts_counter_temp1315--) {
                        {
                            /*Loop's init*/
                            (this->inputsTPParent->m_darts1315) = 0;
                            int m_darts_counter_temp1315 = (this->inputsTPParent->m_darts1315);
                            for (; m_darts_counter_temp1315 < 5; m_darts_counter_temp1315++) {
                                tv[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                  [m_darts_counter_temp1315]
                                    = tv[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                        [m_darts_counter_temp1315]
                                    + (*(this->inputsTPParent->omega_darts1315))
                                        * (b[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                            [m_darts_counter_temp1315][0]
                                                * rsd[i_darts_counter_temp1315]
                                                     [j_darts_counter_temp1315 + 1]
                                                     [(*(this->inputsTPParent->k_darts1315))][0]
                                            + a[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                               [m_darts_counter_temp1315][0]
                                                * rsd[i_darts_counter_temp1315 + 1]
                                                     [j_darts_counter_temp1315]
                                                     [(*(this->inputsTPParent->k_darts1315))][0]
                                            + b[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                               [m_darts_counter_temp1315][1]
                                                * rsd[i_darts_counter_temp1315]
                                                     [j_darts_counter_temp1315 + 1]
                                                     [(*(this->inputsTPParent->k_darts1315))][1]
                                            + a[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                               [m_darts_counter_temp1315][1]
                                                * rsd[i_darts_counter_temp1315 + 1]
                                                     [j_darts_counter_temp1315]
                                                     [(*(this->inputsTPParent->k_darts1315))][1]
                                            + b[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                               [m_darts_counter_temp1315][2]
                                                * rsd[i_darts_counter_temp1315]
                                                     [j_darts_counter_temp1315 + 1]
                                                     [(*(this->inputsTPParent->k_darts1315))][2]
                                            + a[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                               [m_darts_counter_temp1315][2]
                                                * rsd[i_darts_counter_temp1315 + 1]
                                                     [j_darts_counter_temp1315]
                                                     [(*(this->inputsTPParent->k_darts1315))][2]
                                            + b[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                               [m_darts_counter_temp1315][3]
                                                * rsd[i_darts_counter_temp1315]
                                                     [j_darts_counter_temp1315 + 1]
                                                     [(*(this->inputsTPParent->k_darts1315))][3]
                                            + a[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                               [m_darts_counter_temp1315][3]
                                                * rsd[i_darts_counter_temp1315 + 1]
                                                     [j_darts_counter_temp1315]
                                                     [(*(this->inputsTPParent->k_darts1315))][3]
                                            + b[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                               [m_darts_counter_temp1315][4]
                                                * rsd[i_darts_counter_temp1315]
                                                     [j_darts_counter_temp1315 + 1]
                                                     [(*(this->inputsTPParent->k_darts1315))][4]
                                            + a[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                               [m_darts_counter_temp1315][4]
                                                * rsd[i_darts_counter_temp1315 + 1]
                                                     [j_darts_counter_temp1315]
                                                     [(*(this->inputsTPParent->k_darts1315))][4]);
                            }
                            (this->inputsTPParent->m_darts1315) = m_darts_counter_temp1315;
                        }
                        {
                            /*Loop's init*/
                            (this->inputsTPParent->m_darts1315) = 0;
                            int m_darts_counter_temp1315 = (this->inputsTPParent->m_darts1315);
                            for (; m_darts_counter_temp1315 < 5; m_darts_counter_temp1315++) {
                                tmat[m_darts_counter_temp1315][0]
                                    = d[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                       [m_darts_counter_temp1315][0];
                                tmat[m_darts_counter_temp1315][1]
                                    = d[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                       [m_darts_counter_temp1315][1];
                                tmat[m_darts_counter_temp1315][2]
                                    = d[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                       [m_darts_counter_temp1315][2];
                                tmat[m_darts_counter_temp1315][3]
                                    = d[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                       [m_darts_counter_temp1315][3];
                                tmat[m_darts_counter_temp1315][4]
                                    = d[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                       [m_darts_counter_temp1315][4];
                            }
                            (this->inputsTPParent->m_darts1315) = m_darts_counter_temp1315;
                        }
                        (*(this->inputsTPParent->tmp1_darts1315)) = 1. / tmat[0][0];
                        (*(this->inputsTPParent->tmp_darts1315))
                            = (*(this->inputsTPParent->tmp1_darts1315)) * tmat[1][0];
                        tmat[1][1]
                            = tmat[1][1] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[0][1];
                        tmat[1][2]
                            = tmat[1][2] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[0][2];
                        tmat[1][3]
                            = tmat[1][3] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[0][3];
                        tmat[1][4]
                            = tmat[1][4] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[0][4];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][1]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][1]
                            - tv[i_darts_counter_temp1315][j_darts_counter_temp1315][0]
                                * (*(this->inputsTPParent->tmp_darts1315));
                        (*(this->inputsTPParent->tmp_darts1315))
                            = (*(this->inputsTPParent->tmp1_darts1315)) * tmat[2][0];
                        tmat[2][1]
                            = tmat[2][1] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[0][1];
                        tmat[2][2]
                            = tmat[2][2] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[0][2];
                        tmat[2][3]
                            = tmat[2][3] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[0][3];
                        tmat[2][4]
                            = tmat[2][4] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[0][4];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][2]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][2]
                            - tv[i_darts_counter_temp1315][j_darts_counter_temp1315][0]
                                * (*(this->inputsTPParent->tmp_darts1315));
                        (*(this->inputsTPParent->tmp_darts1315))
                            = (*(this->inputsTPParent->tmp1_darts1315)) * tmat[3][0];
                        tmat[3][1]
                            = tmat[3][1] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[0][1];
                        tmat[3][2]
                            = tmat[3][2] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[0][2];
                        tmat[3][3]
                            = tmat[3][3] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[0][3];
                        tmat[3][4]
                            = tmat[3][4] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[0][4];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][3]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][3]
                            - tv[i_darts_counter_temp1315][j_darts_counter_temp1315][0]
                                * (*(this->inputsTPParent->tmp_darts1315));
                        (*(this->inputsTPParent->tmp_darts1315))
                            = (*(this->inputsTPParent->tmp1_darts1315)) * tmat[4][0];
                        tmat[4][1]
                            = tmat[4][1] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[0][1];
                        tmat[4][2]
                            = tmat[4][2] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[0][2];
                        tmat[4][3]
                            = tmat[4][3] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[0][3];
                        tmat[4][4]
                            = tmat[4][4] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[0][4];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][4]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][4]
                            - tv[i_darts_counter_temp1315][j_darts_counter_temp1315][0]
                                * (*(this->inputsTPParent->tmp_darts1315));
                        (*(this->inputsTPParent->tmp1_darts1315)) = 1. / tmat[1][1];
                        (*(this->inputsTPParent->tmp_darts1315))
                            = (*(this->inputsTPParent->tmp1_darts1315)) * tmat[2][1];
                        tmat[2][2]
                            = tmat[2][2] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[1][2];
                        tmat[2][3]
                            = tmat[2][3] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[1][3];
                        tmat[2][4]
                            = tmat[2][4] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[1][4];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][2]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][2]
                            - tv[i_darts_counter_temp1315][j_darts_counter_temp1315][1]
                                * (*(this->inputsTPParent->tmp_darts1315));
                        (*(this->inputsTPParent->tmp_darts1315))
                            = (*(this->inputsTPParent->tmp1_darts1315)) * tmat[3][1];
                        tmat[3][2]
                            = tmat[3][2] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[1][2];
                        tmat[3][3]
                            = tmat[3][3] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[1][3];
                        tmat[3][4]
                            = tmat[3][4] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[1][4];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][3]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][3]
                            - tv[i_darts_counter_temp1315][j_darts_counter_temp1315][1]
                                * (*(this->inputsTPParent->tmp_darts1315));
                        (*(this->inputsTPParent->tmp_darts1315))
                            = (*(this->inputsTPParent->tmp1_darts1315)) * tmat[4][1];
                        tmat[4][2]
                            = tmat[4][2] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[1][2];
                        tmat[4][3]
                            = tmat[4][3] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[1][3];
                        tmat[4][4]
                            = tmat[4][4] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[1][4];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][4]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][4]
                            - tv[i_darts_counter_temp1315][j_darts_counter_temp1315][1]
                                * (*(this->inputsTPParent->tmp_darts1315));
                        (*(this->inputsTPParent->tmp1_darts1315)) = 1. / tmat[2][2];
                        (*(this->inputsTPParent->tmp_darts1315))
                            = (*(this->inputsTPParent->tmp1_darts1315)) * tmat[3][2];
                        tmat[3][3]
                            = tmat[3][3] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[2][3];
                        tmat[3][4]
                            = tmat[3][4] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[2][4];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][3]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][3]
                            - tv[i_darts_counter_temp1315][j_darts_counter_temp1315][2]
                                * (*(this->inputsTPParent->tmp_darts1315));
                        (*(this->inputsTPParent->tmp_darts1315))
                            = (*(this->inputsTPParent->tmp1_darts1315)) * tmat[4][2];
                        tmat[4][3]
                            = tmat[4][3] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[2][3];
                        tmat[4][4]
                            = tmat[4][4] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[2][4];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][4]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][4]
                            - tv[i_darts_counter_temp1315][j_darts_counter_temp1315][2]
                                * (*(this->inputsTPParent->tmp_darts1315));
                        (*(this->inputsTPParent->tmp1_darts1315)) = 1. / tmat[3][3];
                        (*(this->inputsTPParent->tmp_darts1315))
                            = (*(this->inputsTPParent->tmp1_darts1315)) * tmat[4][3];
                        tmat[4][4]
                            = tmat[4][4] - (*(this->inputsTPParent->tmp_darts1315)) * tmat[3][4];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][4]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][4]
                            - tv[i_darts_counter_temp1315][j_darts_counter_temp1315][3]
                                * (*(this->inputsTPParent->tmp_darts1315));
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][4]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][4]
                            / tmat[4][4];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][3]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][3]
                            - tmat[3][4]
                                * tv[i_darts_counter_temp1315][j_darts_counter_temp1315][4];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][3]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][3]
                            / tmat[3][3];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][2]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][2]
                            - tmat[2][3] * tv[i_darts_counter_temp1315][j_darts_counter_temp1315][3]
                            - tmat[2][4]
                                * tv[i_darts_counter_temp1315][j_darts_counter_temp1315][4];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][2]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][2]
                            / tmat[2][2];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][1]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][1]
                            - tmat[1][2] * tv[i_darts_counter_temp1315][j_darts_counter_temp1315][2]
                            - tmat[1][3] * tv[i_darts_counter_temp1315][j_darts_counter_temp1315][3]
                            - tmat[1][4]
                                * tv[i_darts_counter_temp1315][j_darts_counter_temp1315][4];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][1]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][1]
                            / tmat[1][1];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][0]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][0]
                            - tmat[0][1] * tv[i_darts_counter_temp1315][j_darts_counter_temp1315][1]
                            - tmat[0][2] * tv[i_darts_counter_temp1315][j_darts_counter_temp1315][2]
                            - tmat[0][3] * tv[i_darts_counter_temp1315][j_darts_counter_temp1315][3]
                            - tmat[0][4]
                                * tv[i_darts_counter_temp1315][j_darts_counter_temp1315][4];
                        tv[i_darts_counter_temp1315][j_darts_counter_temp1315][0]
                            = tv[i_darts_counter_temp1315][j_darts_counter_temp1315][0]
                            / tmat[0][0];
                        rsd[i_darts_counter_temp1315][j_darts_counter_temp1315]
                           [(*(this->inputsTPParent->k_darts1315))][0]
                            = rsd[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                 [(*(this->inputsTPParent->k_darts1315))][0]
                            - tv[i_darts_counter_temp1315][j_darts_counter_temp1315][0];
                        rsd[i_darts_counter_temp1315][j_darts_counter_temp1315]
                           [(*(this->inputsTPParent->k_darts1315))][1]
                            = rsd[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                 [(*(this->inputsTPParent->k_darts1315))][1]
                            - tv[i_darts_counter_temp1315][j_darts_counter_temp1315][1];
                        rsd[i_darts_counter_temp1315][j_darts_counter_temp1315]
                           [(*(this->inputsTPParent->k_darts1315))][2]
                            = rsd[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                 [(*(this->inputsTPParent->k_darts1315))][2]
                            - tv[i_darts_counter_temp1315][j_darts_counter_temp1315][2];
                        rsd[i_darts_counter_temp1315][j_darts_counter_temp1315]
                           [(*(this->inputsTPParent->k_darts1315))][3]
                            = rsd[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                 [(*(this->inputsTPParent->k_darts1315))][3]
                            - tv[i_darts_counter_temp1315][j_darts_counter_temp1315][3];
                        rsd[i_darts_counter_temp1315][j_darts_counter_temp1315]
                           [(*(this->inputsTPParent->k_darts1315))][4]
                            = rsd[i_darts_counter_temp1315][j_darts_counter_temp1315]
                                 [(*(this->inputsTPParent->k_darts1315))][4]
                            - tv[i_darts_counter_temp1315][j_darts_counter_temp1315][4];
                    }
                    (this->inputsTPParent->j_darts1315) = j_darts_counter_temp1315;
                }
                if (i_darts_counter_temp1315 != (*(this->inputsTPParent->iend_darts1315))) {
                    flag[i_darts_counter_temp1315 + 1] = 0;
                }
                if (i_darts_counter_temp1315 != (*(this->inputsTPParent->ist_darts1315))) {
                    flag[i_darts_counter_temp1315] = 1;
                }
            }
            (this->inputsTPParent->i_darts1315) = i_darts_counter_temp1315;
        }
        /*Signaling next codelet from last stmt in the codelet*/
        /*Find and signal the next codelet*/

        myTP->controlTPParent->nextCodeletsbuts[this->getID()]->decDep();
    } else {
        /*Find and signal the next codelet*/

        myTP->nextCodeletsbuts[this->getID()]->decDep();
    }
}
TP2::TP2(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_mainNextCodelet,
    darts::Codelet* in_mainSyncCodelet, TP2** in_ptrToThisFunctionTP, int in_nx, int in_ny,
    int in_nz, int in_k, double in_omega, int in_ist, int in_iend, int in_jst, int in_jend,
    int in_nx0, int in_ny0)
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
    , ist_darts2(new int[this->numThreads])
    , iend_darts2(new int[this->numThreads])
    , jst_darts2(new int[this->numThreads])
    , jend_darts2(new int[this->numThreads])
    , nx0_darts2(new int[this->numThreads])
    , ny0_darts2(new int[this->numThreads])
    , i_darts2(new int[this->numThreads])
    , j_darts2(new int[this->numThreads])
    , m_darts2(new int[this->numThreads])
    , tmp_darts2(new double[this->numThreads])
    , tmp1_darts2(new double[this->numThreads])
    , TP1213Ptr(new TP1213*[NUMTPS1213])
    , TP1213_alreadyLaunched(new size_t[NUMTPS1213])
    , numTPsSet1213(0)
    , numTPsReady1213(0)
    , TPsToUse1213(NUMTPS1213)
    , codeletsPerTP1213(this->numThreads / NUMTPS1213)
    , totalCodelets1213(this->TPsToUse1213 * this->codeletsPerTP1213)
    , TP1315_alreadyLaunched(0)
    , checkInCodelets1213(new _checkInCodelets1213[this->numThreads])
    , barrierCodelets1213(new _barrierCodelets1213[1])
    , checkInCodelets1315(new _checkInCodelets1315[this->numThreads])
{
    barrierCodelets1213[0] = _barrierCodelets1213(NUMTPS1213, NUMTPS1213, this, 0);
    _checkInCodelets1315* checkInCodelets1315Ptr = (this->checkInCodelets1315);
    _checkInCodelets1213* checkInCodelets1213Ptr = (this->checkInCodelets1213);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets1213);
#endif
    for (int i = 0; i < NUMTPS1213; i++) {
        TP1213Ptr[i] = nullptr;
        TP1213_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++) {
        (*checkInCodelets1315Ptr) = _checkInCodelets1315(1, 1, this, codeletCounter);
        checkInCodelets1315Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets1213Ptr) = _checkInCodelets1213(2, 1, this, codeletCounter);
#else
        (*checkInCodelets1213Ptr) = _checkInCodelets1213(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets1213Ptr).decDep();
        checkInCodelets1213Ptr++;
    }
    if (this->numThreads == 1) {
        this->nextCodeletsbuts[0] = in_mainNextCodelet;
        this->nextSyncCodeletsbuts[0] = in_mainSyncCodelet;
        this->nx_darts2[0] = in_nx;
        this->ny_darts2[0] = in_ny;
        this->nz_darts2[0] = in_nz;
        this->k_darts2[0] = in_k;
        this->omega_darts2[0] = in_omega;
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
    delete[] checkInCodelets1315;
    delete[] barrierCodelets1213;
    delete[] checkInCodelets1213;
    delete[] nextCodeletsbuts;
    delete[] nextSyncCodeletsbuts;
    delete[] nx_darts2;
    delete[] ny_darts2;
    delete[] nz_darts2;
    delete[] k_darts2;
    delete[] omega_darts2;
    delete[] ist_darts2;
    delete[] iend_darts2;
    delete[] jst_darts2;
    delete[] jend_darts2;
    delete[] nx0_darts2;
    delete[] ny0_darts2;
    delete[] i_darts2;
    delete[] j_darts2;
    delete[] m_darts2;
    delete[] tmp_darts2;
    delete[] tmp1_darts2;
}
void TP2::setNewInputs(int in_nx, int in_ny, int in_nz, int in_k, double in_omega, int in_ist,
    int in_iend, int in_jst, int in_jend, int in_nx0, int in_ny0, size_t codeletID)
{
    this->nx_darts2[codeletID] = in_nx;
    this->ny_darts2[codeletID] = in_ny;
    this->nz_darts2[codeletID] = in_nz;
    this->k_darts2[codeletID] = in_k;
    this->omega_darts2[codeletID] = in_omega;
    this->ist_darts2[codeletID] = in_ist;
    this->iend_darts2[codeletID] = in_iend;
    this->jst_darts2[codeletID] = in_jst;
    this->jend_darts2[codeletID] = in_jend;
    this->nx0_darts2[codeletID] = in_nx0;
    this->ny0_darts2[codeletID] = in_ny0;
}
/*TP1213: OMPForDirective*/
void TP1213::_barrierCodelets1213::fire(void)
{
    TP1213* myTP = static_cast<TP1213*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets1213[0].decDep();
}
bool TP1213::requestNewRangeIterations1213(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 1*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet1213 * codeletID;
        int tempEndRange = rangePerCodelet1213 * (codeletID + 1);
        if (remainderRange1213 != 0) {
            if (codeletID < (uint32_t)remainderRange1213) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange1213;
                tempEndRange += remainderRange1213;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration1213;
        tempEndRange = tempEndRange * 1 + minIteration1213;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration1213 < lastIteration1213) {
            (this->inputsTPParent->i_darts1213[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts1213[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == 0) {
            *endRange = *endRange - 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration1213;
        }
    }
    return isThereNewIteration;
}
void TP1213::_checkInCodelets1214::fire(void)
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
    this->inputsTPParent->iend_darts1213[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iend_darts2[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ist_darts1213[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->ist_darts2[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jend_darts1213[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jend_darts2[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jst_darts1213[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jst_darts2[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->k_darts1213[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->k_darts2[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->omega_darts1213[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->omega_darts2[this->getID()]);

    /*printing node 1214: ForStmt*/
    /*var: i*/
    /*var: iend*/
    /*var: ist*/
    /*var: j*/
    /*var: jend*/
    /*var: jst*/
    /*var: k*/
    /*var: m*/
    /*var: omega*/
    int* i = &(this->inputsTPParent->i_darts1213[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts1213[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jend = &(this->inputsTPParent->jend_darts1213[this->getLocalID()]);
    (void)jend /*OMP_SHARED_PRIVATE*/;
    int** jst = &(this->inputsTPParent->jst_darts1213[this->getLocalID()]);
    (void)jst /*OMP_SHARED_PRIVATE*/;
    int** k = &(this->inputsTPParent->k_darts1213[this->getLocalID()]);
    (void)k /*OMP_SHARED_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts1213[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** omega = &(this->inputsTPParent->omega_darts1213[this->getLocalID()]);
    (void)omega /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations1213(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets1213[0].decDep();
        return;
    }
    for (int i_darts_counter_temp1213 = (*i); i_darts_counter_temp1213 >= endRange
         && i_darts_counter_temp1213 >= this->inputsTPParent->lastIteration1213;
         i_darts_counter_temp1213--) {
        {
            {
                /*Loop's init*/
                (*j) = (*(*jend));
                int j_darts_counter_temp1213 = (*j);
                for (; j_darts_counter_temp1213 >= (*(*jst)); j_darts_counter_temp1213--) {
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp1213 = (*m);
                        for (; m_darts_counter_temp1213 < 5; m_darts_counter_temp1213++) {
                            tv[(i_darts_counter_temp1213)][j_darts_counter_temp1213]
                              [m_darts_counter_temp1213]
                                = (*(*omega))
                                * (c[(i_darts_counter_temp1213)][j_darts_counter_temp1213]
                                    [m_darts_counter_temp1213][0]
                                        * rsd[(i_darts_counter_temp1213)][j_darts_counter_temp1213]
                                             [(*(*k)) + 1][0]
                                    + c[(i_darts_counter_temp1213)][j_darts_counter_temp1213]
                                       [m_darts_counter_temp1213][1]
                                        * rsd[(i_darts_counter_temp1213)][j_darts_counter_temp1213]
                                             [(*(*k)) + 1][1]
                                    + c[(i_darts_counter_temp1213)][j_darts_counter_temp1213]
                                       [m_darts_counter_temp1213][2]
                                        * rsd[(i_darts_counter_temp1213)][j_darts_counter_temp1213]
                                             [(*(*k)) + 1][2]
                                    + c[(i_darts_counter_temp1213)][j_darts_counter_temp1213]
                                       [m_darts_counter_temp1213][3]
                                        * rsd[(i_darts_counter_temp1213)][j_darts_counter_temp1213]
                                             [(*(*k)) + 1][3]
                                    + c[(i_darts_counter_temp1213)][j_darts_counter_temp1213]
                                       [m_darts_counter_temp1213][4]
                                        * rsd[(i_darts_counter_temp1213)][j_darts_counter_temp1213]
                                             [(*(*k)) + 1][4]);
                        }
                        (*m) = m_darts_counter_temp1213;
                    }
                }
                (*j) = j_darts_counter_temp1213;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets1213[0].decDep();
}
TP1213::TP1213(int in_numThreads, int in_mainCodeletID, TP2* in_TPParent, int in_initIteration,
    int in_lastIteration, TP1213** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts1213(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iend_darts1213(new int*[this->numThreads])
    , ist_darts1213(new int*[this->numThreads])
    , j_darts1213(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jend_darts1213(new int*[this->numThreads])
    , jst_darts1213(new int*[this->numThreads])
    , k_darts1213(new int*[this->numThreads])
    , m_darts1213(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , omega_darts1213(new double*[this->numThreads])
    , initIteration1213(in_initIteration)
    , lastIteration1213(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets1213(new _barrierCodelets1213[1])
    , checkInCodelets1214(new _checkInCodelets1214[this->numThreads])
{
    /*Initialize the loop parameters*/
    range1213 = abs(lastIteration1213 - initIteration1213) / 1;
    rangePerCodelet1213 = range1213 / numThreads;
    minIteration1213 = min<int>(lastIteration1213, initIteration1213);
    remainderRange1213 = range1213 % numThreads;
    /*Initialize inputs and vars.*/
    this->iend_darts1213 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ist_darts1213 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jend_darts1213 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jst_darts1213 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->k_darts1213 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->omega_darts1213
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets1213[0] = _barrierCodelets1213(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets1214* checkInCodelets1214Ptr = (this->checkInCodelets1214);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets1214);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets1214Ptr) = _checkInCodelets1214(2, 1, this, codeletCounter);
#else
        (*checkInCodelets1214Ptr) = _checkInCodelets1214(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets1214Ptr).decDep();
        checkInCodelets1214Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP1213::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets1214[localID].setID(codeletID);
    this->checkInCodelets1214[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets1214[localID + this->baseNumThreads * i]
            = _checkInCodelets1214(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets1214[localID + this->baseNumThreads * i]
            = _checkInCodelets1214(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets1214[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets1214[localID + this->baseNumThreads * i].decDep();
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
TP1213::~TP1213()
{
    delete[] iend_darts1213;
    delete[] ist_darts1213;
    delete[] jend_darts1213;
    delete[] jst_darts1213;
    delete[] k_darts1213;
    delete[] omega_darts1213;
    delete[] barrierCodelets1213;
    delete[] checkInCodelets1214;
}
/*TP2199: OMPParallelDirective*/
void TP2199::_barrierCodelets2199::fire(void)
{
    TP2199* myTP = static_cast<TP2199*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP2199::_checkInCodelets2201::fire(void)
{
    /*Init the vars for this region*/

    /*printing node 2201: DeclStmt*/

    /*printing node 2202: DeclStmt*/

    /*printing node 2203: DeclStmt*/

    /*printing node 2204: DeclStmt*/

    /*printing node 2205: DeclStmt*/

    /*printing node 2206: DeclStmt*/

    /*printing node 2207: DeclStmt*/

    /*printing node 2208: DeclStmt*/

    /*printing node 2209: DeclStmt*/

    /*printing node 2210: DeclStmt*/

    /*printing node 2211: DeclStmt*/

    /*printing node 2212: DeclStmt*/

    /*printing node 2213: DeclStmt*/

    /*printing node 2214: DeclStmt*/

    /*printing node 2215: DeclStmt*/

    /*printing node 2216: DeclStmt*/

    /*printing node 2217: BinaryOperator*/
    (this->inputsTPParent->dsspm_darts2199[this->getID()]) = dssp;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 2201 nextRegion: 2218 */
    myTP->controlTPParent->checkInCodelets2218[this->getID()].decDep();
}
void TP2199::_checkInCodelets2218::fire(void)
{
    /*region 2218 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP2218;
    if (idx < myTP->TPsToUse2218) {
        if (!__sync_val_compare_and_swap(&(myTP->TP2218_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(nx - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse2218;
            int minIteration = min<int>(nx, 0);
            int remainderRange = range % myTP->TPsToUse2218;
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
            if (idx == myTP->TPsToUse2218 - 1) {
                lastIteration = nx;
            }
#if USEINVOKE == 1
            invoke<TP2218>(myTP, myTP->codeletsPerTP2218 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP2218Ptr[idx]));
#else
            place<TP2218>(idx, myTP, myTP->codeletsPerTP2218 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP2218Ptr[idx]));
#endif
        } else {
            if (myTP->TP2218Ptr[idx] != nullptr) {
                myTP->TP2218Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2199::_barrierCodelets2218::fire(void)
{
    TP2199* myTP = static_cast<TP2199*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets2269[codeletsCounter].decDep();
        }
    }
}
void TP2199::_checkInCodelets2269::fire(void)
{
    /*region 2269 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP2269;
    if (idx < myTP->TPsToUse2269) {
        if (!__sync_val_compare_and_swap(&(myTP->TP2269_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(nx - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse2269;
            int minIteration = min<int>(nx, 0);
            int remainderRange = range % myTP->TPsToUse2269;
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
            if (idx == myTP->TPsToUse2269 - 1) {
                lastIteration = nx;
            }
#if USEINVOKE == 1
            invoke<TP2269>(myTP, myTP->codeletsPerTP2269 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP2269Ptr[idx]));
#else
            place<TP2269>(idx, myTP, myTP->codeletsPerTP2269 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP2269Ptr[idx]));
#endif
        } else {
            if (myTP->TP2269Ptr[idx] != nullptr) {
                myTP->TP2269Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2199::_barrierCodelets2269::fire(void)
{
    TP2199* myTP = static_cast<TP2199*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets2401[codeletsCounter].decDep();
        }
    }
}
void TP2199::_checkInCodelets2401::fire(void)
{

    /*printing node 2401: BinaryOperator*/
    (this->inputsTPParent->L1_darts2199[this->getID()]) = 0;

    /*printing node 2402: BinaryOperator*/
    (this->inputsTPParent->L2_darts2199[this->getID()]) = nx - 1;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 2401 nextRegion: 2404 */
    myTP->controlTPParent->checkInCodelets2404[this->getID()].decDep();
}
void TP2199::_checkInCodelets2404::fire(void)
{
    /*region 2404 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP2404;
    if (idx < myTP->TPsToUse2404) {
        if (!__sync_val_compare_and_swap(&(myTP->TP2404_alreadyLaunched[idx]), 0, 1)) {
            int range = abs((this->inputsTPParent->L2_darts2199[this->getID()])
                            - (this->inputsTPParent->L1_darts2199[this->getID()]))
                / 1;
            int rangePerCodelet = range / myTP->TPsToUse2404;
            int minIteration = min<int>((this->inputsTPParent->L2_darts2199[this->getID()]),
                (this->inputsTPParent->L1_darts2199[this->getID()]));
            int remainderRange = range % myTP->TPsToUse2404;
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
            if ((this->inputsTPParent->L1_darts2199[this->getID()])
                < (this->inputsTPParent->L2_darts2199[this->getID()])) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse2404 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse2404 - 1) {
                lastIteration = (this->inputsTPParent->L2_darts2199[this->getID()]);
            }
#if USEINVOKE == 1
            invoke<TP2404>(myTP, myTP->codeletsPerTP2404 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP2404Ptr[idx]));
#else
            place<TP2404>(idx, myTP, myTP->codeletsPerTP2404 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP2404Ptr[idx]));
#endif
        } else {
            if (myTP->TP2404Ptr[idx] != nullptr) {
                myTP->TP2404Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2199::_barrierCodelets2404::fire(void)
{
    TP2199* myTP = static_cast<TP2199*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets2553[codeletsCounter].decDep();
        }
    }
}
void TP2199::_checkInCodelets2553::fire(void)
{
    /*region 2553 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP2553;
    if (idx < myTP->TPsToUse2553) {
        if (!__sync_val_compare_and_swap(&(myTP->TP2553_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(jend - jst) / 1;
            int rangePerCodelet = range / myTP->TPsToUse2553;
            int minIteration = min<int>(jend, jst);
            int remainderRange = range % myTP->TPsToUse2553;
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
            if (idx == myTP->TPsToUse2553 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse2553 - 1) {
                lastIteration = jend;
            }
#if USEINVOKE == 1
            invoke<TP2553>(myTP, myTP->codeletsPerTP2553 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP2553Ptr[idx]));
#else
            place<TP2553>(idx, myTP, myTP->codeletsPerTP2553 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP2553Ptr[idx]));
#endif
        } else {
            if (myTP->TP2553Ptr[idx] != nullptr) {
                myTP->TP2553Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2199::_barrierCodelets2553::fire(void)
{
    TP2199* myTP = static_cast<TP2199*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets3188[codeletsCounter].decDep();
        }
    }
}
void TP2199::_checkInCodelets3188::fire(void)
{

    /*printing node 3188: BinaryOperator*/
    (this->inputsTPParent->L1_darts2199[this->getID()]) = 0;

    /*printing node 3189: BinaryOperator*/
    (this->inputsTPParent->L2_darts2199[this->getID()]) = ny - 1;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 3188 nextRegion: 3191 */
    myTP->controlTPParent->checkInCodelets3191[this->getID()].decDep();
}
void TP2199::_checkInCodelets3191::fire(void)
{
    /*region 3191 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP3191;
    if (idx < myTP->TPsToUse3191) {
        if (!__sync_val_compare_and_swap(&(myTP->TP3191_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse3191;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse3191;
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
            if (idx == myTP->TPsToUse3191 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse3191 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP3191>(myTP, myTP->codeletsPerTP3191 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP3191Ptr[idx]));
#else
            place<TP3191>(idx, myTP, myTP->codeletsPerTP3191 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP3191Ptr[idx]));
#endif
        } else {
            if (myTP->TP3191Ptr[idx] != nullptr) {
                myTP->TP3191Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2199::_barrierCodelets3191::fire(void)
{
    TP2199* myTP = static_cast<TP2199*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets3340[codeletsCounter].decDep();
        }
    }
}
void TP2199::_checkInCodelets3340::fire(void)
{
    /*region 3340 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP3340;
    if (idx < myTP->TPsToUse3340) {
        if (!__sync_val_compare_and_swap(&(myTP->TP3340_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse3340;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse3340;
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
            if (idx == myTP->TPsToUse3340 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse3340 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP3340>(myTP, myTP->codeletsPerTP3340 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP3340Ptr[idx]));
#else
            place<TP3340>(idx, myTP, myTP->codeletsPerTP3340 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP3340Ptr[idx]));
#endif
        } else {
            if (myTP->TP3340Ptr[idx] != nullptr) {
                myTP->TP3340Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2199::_barrierCodelets3340::fire(void)
{
    TP2199* myTP = static_cast<TP2199*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets3975[codeletsCounter].decDep();
        }
    }
}
void TP2199::_checkInCodelets3975::fire(void)
{
    /*region 3975 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP3975;
    if (idx < myTP->TPsToUse3975) {
        if (!__sync_val_compare_and_swap(&(myTP->TP3975_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse3975;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse3975;
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
            if (idx == myTP->TPsToUse3975 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse3975 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP3975>(myTP, myTP->codeletsPerTP3975 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP3975Ptr[idx]));
#else
            place<TP3975>(idx, myTP, myTP->codeletsPerTP3975 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP3975Ptr[idx]));
#endif
        } else {
            if (myTP->TP3975Ptr[idx] != nullptr) {
                myTP->TP3975Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2199::_barrierCodelets3975::fire(void)
{
    TP2199* myTP = static_cast<TP2199*>(myTP_);
    myTP->TPParent->barrierCodelets2199[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets2199[0]));
}
TP2199::TP2199(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts2199(new int[this->numThreads]) /*VARIABLE*/
    , L2_darts2199(new int[this->numThreads]) /*VARIABLE*/
    , dsspm_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , eta_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , i_darts2199(new int[this->numThreads]) /*VARIABLE*/
    , iend1_darts2199(new int[this->numThreads]) /*VARIABLE*/
    , iglob_darts2199(new int[this->numThreads]) /*VARIABLE*/
    , ist1_darts2199(new int[this->numThreads]) /*VARIABLE*/
    , j_darts2199(new int[this->numThreads]) /*VARIABLE*/
    , jend1_darts2199(new int[this->numThreads]) /*VARIABLE*/
    , jglob_darts2199(new int[this->numThreads]) /*VARIABLE*/
    , jst1_darts2199(new int[this->numThreads]) /*VARIABLE*/
    , k_darts2199(new int[this->numThreads]) /*VARIABLE*/
    , m_darts2199(new int[this->numThreads]) /*VARIABLE*/
    , q_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , tmp_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u21_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u21i_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u21im1_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u21j_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u21jm1_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u21k_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u21km1_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u31_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u31i_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u31im1_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u31j_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u31jm1_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u31k_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u31km1_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u41_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u41i_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u41im1_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u41j_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u41jm1_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u41k_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u41km1_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u51i_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u51im1_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u51j_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u51jm1_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u51k_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , u51km1_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , xi_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , zeta_darts2199(new double[this->numThreads]) /*VARIABLE*/
    , TP2218Ptr(new TP2218*[NUMTPS2218])
    , TP2218_alreadyLaunched(new size_t[NUMTPS2218])
    , numTPsSet2218(0)
    , numTPsReady2218(0)
    , TPsToUse2218(NUMTPS2218)
    , codeletsPerTP2218(this->numThreads / NUMTPS2218)
    , totalCodelets2218(this->TPsToUse2218 * this->codeletsPerTP2218)
    , TP2269Ptr(new TP2269*[NUMTPS2269])
    , TP2269_alreadyLaunched(new size_t[NUMTPS2269])
    , numTPsSet2269(0)
    , numTPsReady2269(0)
    , TPsToUse2269(NUMTPS2269)
    , codeletsPerTP2269(this->numThreads / NUMTPS2269)
    , totalCodelets2269(this->TPsToUse2269 * this->codeletsPerTP2269)
    , TP2404Ptr(new TP2404*[NUMTPS2404])
    , TP2404_alreadyLaunched(new size_t[NUMTPS2404])
    , numTPsSet2404(0)
    , numTPsReady2404(0)
    , TPsToUse2404(NUMTPS2404)
    , codeletsPerTP2404(this->numThreads / NUMTPS2404)
    , totalCodelets2404(this->TPsToUse2404 * this->codeletsPerTP2404)
    , TP2553Ptr(new TP2553*[NUMTPS2553])
    , TP2553_alreadyLaunched(new size_t[NUMTPS2553])
    , numTPsSet2553(0)
    , numTPsReady2553(0)
    , TPsToUse2553(NUMTPS2553)
    , codeletsPerTP2553(this->numThreads / NUMTPS2553)
    , totalCodelets2553(this->TPsToUse2553 * this->codeletsPerTP2553)
    , TP3191Ptr(new TP3191*[NUMTPS3191])
    , TP3191_alreadyLaunched(new size_t[NUMTPS3191])
    , numTPsSet3191(0)
    , numTPsReady3191(0)
    , TPsToUse3191(NUMTPS3191)
    , codeletsPerTP3191(this->numThreads / NUMTPS3191)
    , totalCodelets3191(this->TPsToUse3191 * this->codeletsPerTP3191)
    , TP3340Ptr(new TP3340*[NUMTPS3340])
    , TP3340_alreadyLaunched(new size_t[NUMTPS3340])
    , numTPsSet3340(0)
    , numTPsReady3340(0)
    , TPsToUse3340(NUMTPS3340)
    , codeletsPerTP3340(this->numThreads / NUMTPS3340)
    , totalCodelets3340(this->TPsToUse3340 * this->codeletsPerTP3340)
    , TP3975Ptr(new TP3975*[NUMTPS3975])
    , TP3975_alreadyLaunched(new size_t[NUMTPS3975])
    , numTPsSet3975(0)
    , numTPsReady3975(0)
    , TPsToUse3975(NUMTPS3975)
    , codeletsPerTP3975(this->numThreads / NUMTPS3975)
    , totalCodelets3975(this->TPsToUse3975 * this->codeletsPerTP3975)
    , barrierCodelets2199(new _barrierCodelets2199[1])
    , checkInCodelets2201(new _checkInCodelets2201[this->numThreads])
    , checkInCodelets2218(new _checkInCodelets2218[this->numThreads])
    , barrierCodelets2218(new _barrierCodelets2218[1])
    , checkInCodelets2269(new _checkInCodelets2269[this->numThreads])
    , barrierCodelets2269(new _barrierCodelets2269[1])
    , checkInCodelets2401(new _checkInCodelets2401[this->numThreads])
    , checkInCodelets2404(new _checkInCodelets2404[this->numThreads])
    , barrierCodelets2404(new _barrierCodelets2404[1])
    , checkInCodelets2553(new _checkInCodelets2553[this->numThreads])
    , barrierCodelets2553(new _barrierCodelets2553[1])
    , checkInCodelets3188(new _checkInCodelets3188[this->numThreads])
    , checkInCodelets3191(new _checkInCodelets3191[this->numThreads])
    , barrierCodelets3191(new _barrierCodelets3191[1])
    , checkInCodelets3340(new _checkInCodelets3340[this->numThreads])
    , barrierCodelets3340(new _barrierCodelets3340[1])
    , checkInCodelets3975(new _checkInCodelets3975[this->numThreads])
    , barrierCodelets3975(new _barrierCodelets3975[1])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets2199[0] = _barrierCodelets2199(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets3975[0] = _barrierCodelets3975(NUMTPS3975, NUMTPS3975, this, 0);
    barrierCodelets3340[0] = _barrierCodelets3340(NUMTPS3340, NUMTPS3340, this, 0);
    barrierCodelets3191[0] = _barrierCodelets3191(NUMTPS3191, NUMTPS3191, this, 0);
    barrierCodelets2553[0] = _barrierCodelets2553(NUMTPS2553, NUMTPS2553, this, 0);
    barrierCodelets2404[0] = _barrierCodelets2404(NUMTPS2404, NUMTPS2404, this, 0);
    barrierCodelets2269[0] = _barrierCodelets2269(NUMTPS2269, NUMTPS2269, this, 0);
    barrierCodelets2218[0] = _barrierCodelets2218(NUMTPS2218, NUMTPS2218, this, 0);
    _checkInCodelets3975* checkInCodelets3975Ptr = (this->checkInCodelets3975);
    for (int i = 0; i < NUMTPS3975; i++) {
        TP3975Ptr[i] = nullptr;
        TP3975_alreadyLaunched[i] = 0;
    }
    _checkInCodelets3340* checkInCodelets3340Ptr = (this->checkInCodelets3340);
    for (int i = 0; i < NUMTPS3340; i++) {
        TP3340Ptr[i] = nullptr;
        TP3340_alreadyLaunched[i] = 0;
    }
    _checkInCodelets3191* checkInCodelets3191Ptr = (this->checkInCodelets3191);
    for (int i = 0; i < NUMTPS3191; i++) {
        TP3191Ptr[i] = nullptr;
        TP3191_alreadyLaunched[i] = 0;
    }
    _checkInCodelets3188* checkInCodelets3188Ptr = (this->checkInCodelets3188);
    _checkInCodelets2553* checkInCodelets2553Ptr = (this->checkInCodelets2553);
    for (int i = 0; i < NUMTPS2553; i++) {
        TP2553Ptr[i] = nullptr;
        TP2553_alreadyLaunched[i] = 0;
    }
    _checkInCodelets2404* checkInCodelets2404Ptr = (this->checkInCodelets2404);
    for (int i = 0; i < NUMTPS2404; i++) {
        TP2404Ptr[i] = nullptr;
        TP2404_alreadyLaunched[i] = 0;
    }
    _checkInCodelets2401* checkInCodelets2401Ptr = (this->checkInCodelets2401);
    _checkInCodelets2269* checkInCodelets2269Ptr = (this->checkInCodelets2269);
    for (int i = 0; i < NUMTPS2269; i++) {
        TP2269Ptr[i] = nullptr;
        TP2269_alreadyLaunched[i] = 0;
    }
    _checkInCodelets2218* checkInCodelets2218Ptr = (this->checkInCodelets2218);
    for (int i = 0; i < NUMTPS2218; i++) {
        TP2218Ptr[i] = nullptr;
        TP2218_alreadyLaunched[i] = 0;
    }
    _checkInCodelets2201* checkInCodelets2201Ptr = (this->checkInCodelets2201);
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets3975Ptr) = _checkInCodelets3975(1, 1, this, codeletCounter);
        checkInCodelets3975Ptr++;
        (*checkInCodelets3340Ptr) = _checkInCodelets3340(1, 1, this, codeletCounter);
        checkInCodelets3340Ptr++;
        (*checkInCodelets3191Ptr) = _checkInCodelets3191(1, 1, this, codeletCounter);
        checkInCodelets3191Ptr++;
        (*checkInCodelets3188Ptr) = _checkInCodelets3188(1, 1, this, codeletCounter);
        checkInCodelets3188Ptr++;
        (*checkInCodelets2553Ptr) = _checkInCodelets2553(1, 1, this, codeletCounter);
        checkInCodelets2553Ptr++;
        (*checkInCodelets2404Ptr) = _checkInCodelets2404(1, 1, this, codeletCounter);
        checkInCodelets2404Ptr++;
        (*checkInCodelets2401Ptr) = _checkInCodelets2401(1, 1, this, codeletCounter);
        checkInCodelets2401Ptr++;
        (*checkInCodelets2269Ptr) = _checkInCodelets2269(1, 1, this, codeletCounter);
        checkInCodelets2269Ptr++;
        (*checkInCodelets2218Ptr) = _checkInCodelets2218(1, 1, this, codeletCounter);
        checkInCodelets2218Ptr++;
        (*checkInCodelets2201Ptr) = _checkInCodelets2201(1, 1, this, codeletCounter);
        (*checkInCodelets2201Ptr).decDep();
        checkInCodelets2201Ptr++;
    }
}
TP2199::~TP2199()
{
    delete[] L1_darts2199;
    delete[] L2_darts2199;
    delete[] dsspm_darts2199;
    delete[] eta_darts2199;
    delete[] i_darts2199;
    delete[] iend1_darts2199;
    delete[] iglob_darts2199;
    delete[] ist1_darts2199;
    delete[] j_darts2199;
    delete[] jend1_darts2199;
    delete[] jglob_darts2199;
    delete[] jst1_darts2199;
    delete[] k_darts2199;
    delete[] m_darts2199;
    delete[] q_darts2199;
    delete[] tmp_darts2199;
    delete[] u21_darts2199;
    delete[] u21i_darts2199;
    delete[] u21im1_darts2199;
    delete[] u21j_darts2199;
    delete[] u21jm1_darts2199;
    delete[] u21k_darts2199;
    delete[] u21km1_darts2199;
    delete[] u31_darts2199;
    delete[] u31i_darts2199;
    delete[] u31im1_darts2199;
    delete[] u31j_darts2199;
    delete[] u31jm1_darts2199;
    delete[] u31k_darts2199;
    delete[] u31km1_darts2199;
    delete[] u41_darts2199;
    delete[] u41i_darts2199;
    delete[] u41im1_darts2199;
    delete[] u41j_darts2199;
    delete[] u41jm1_darts2199;
    delete[] u41k_darts2199;
    delete[] u41km1_darts2199;
    delete[] u51i_darts2199;
    delete[] u51im1_darts2199;
    delete[] u51j_darts2199;
    delete[] u51jm1_darts2199;
    delete[] u51k_darts2199;
    delete[] u51km1_darts2199;
    delete[] xi_darts2199;
    delete[] zeta_darts2199;
    delete[] barrierCodelets2199;
    delete[] barrierCodelets3975;
    delete[] checkInCodelets3975;
    delete[] barrierCodelets3340;
    delete[] checkInCodelets3340;
    delete[] barrierCodelets3191;
    delete[] checkInCodelets3191;
    delete[] checkInCodelets3188;
    delete[] barrierCodelets2553;
    delete[] checkInCodelets2553;
    delete[] barrierCodelets2404;
    delete[] checkInCodelets2404;
    delete[] checkInCodelets2401;
    delete[] barrierCodelets2269;
    delete[] checkInCodelets2269;
    delete[] barrierCodelets2218;
    delete[] checkInCodelets2218;
    delete[] checkInCodelets2201;
}
/*TP2218: OMPForDirective*/
void TP2218::_barrierCodelets2218::fire(void)
{
    TP2218* myTP = static_cast<TP2218*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets2218[0].decDep();
}
bool TP2218::requestNewRangeIterations2218(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet2218 * codeletID;
        int tempEndRange = rangePerCodelet2218 * (codeletID + 1);
        if (remainderRange2218 != 0) {
            if (codeletID < (uint32_t)remainderRange2218) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange2218;
                tempEndRange += remainderRange2218;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration2218;
        tempEndRange = tempEndRange * 1 + minIteration2218;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration2218 < lastIteration2218) {
            (this->inputsTPParent->i_darts2218[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts2218[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration2218;
        }
    }
    return isThereNewIteration;
}
void TP2218::_checkInCodelets2219::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 2219: ForStmt*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    int* i = &(this->inputsTPParent->i_darts2218[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts2218[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts2218[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts2218[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations2218(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets2218[0].decDep();
        return;
    }
    for (int i_darts_counter_temp2218 = (*i); i_darts_counter_temp2218 < endRange
         && i_darts_counter_temp2218 < this->inputsTPParent->lastIteration2218;
         i_darts_counter_temp2218++) {
        {
            {
                /*Loop's init*/
                (*j) = 0;
                int j_darts_counter_temp2218 = (*j);
                for (; j_darts_counter_temp2218 < ny; j_darts_counter_temp2218++) {
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp2218 = (*k);
                        for (; k_darts_counter_temp2218 < nz; k_darts_counter_temp2218++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp2218 = (*m);
                                for (; m_darts_counter_temp2218 < 5; m_darts_counter_temp2218++) {
                                    frct[(i_darts_counter_temp2218)][j_darts_counter_temp2218]
                                        [k_darts_counter_temp2218][m_darts_counter_temp2218]
                                        = 0.;
                                }
                                (*m) = m_darts_counter_temp2218;
                            }
                        }
                        (*k) = k_darts_counter_temp2218;
                    }
                }
                (*j) = j_darts_counter_temp2218;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets2218[0].decDep();
}
TP2218::TP2218(int in_numThreads, int in_mainCodeletID, TP2199* in_TPParent, int in_initIteration,
    int in_lastIteration, TP2218** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts2218(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts2218(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts2218(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts2218(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration2218(in_initIteration)
    , lastIteration2218(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets2218(new _barrierCodelets2218[1])
    , checkInCodelets2219(new _checkInCodelets2219[this->numThreads])
{
    /*Initialize the loop parameters*/
    range2218 = abs(lastIteration2218 - initIteration2218) / 1;
    rangePerCodelet2218 = range2218 / numThreads;
    minIteration2218 = min<int>(lastIteration2218, initIteration2218);
    remainderRange2218 = range2218 % numThreads;
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets2218[0] = _barrierCodelets2218(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets2219* checkInCodelets2219Ptr = (this->checkInCodelets2219);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets2219);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets2219Ptr) = _checkInCodelets2219(2, 1, this, codeletCounter);
#else
        (*checkInCodelets2219Ptr) = _checkInCodelets2219(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets2219Ptr).decDep();
        checkInCodelets2219Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP2218::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets2219[localID].setID(codeletID);
    this->checkInCodelets2219[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets2219[localID + this->baseNumThreads * i]
            = _checkInCodelets2219(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets2219[localID + this->baseNumThreads * i]
            = _checkInCodelets2219(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets2219[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets2219[localID + this->baseNumThreads * i].decDep();
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
TP2218::~TP2218()
{
    delete[] barrierCodelets2218;
    delete[] checkInCodelets2219;
}
/*TP2269: OMPForDirective*/
void TP2269::_barrierCodelets2269::fire(void)
{
    TP2269* myTP = static_cast<TP2269*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets2269[0].decDep();
}
bool TP2269::requestNewRangeIterations2269(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet2269 * codeletID;
        int tempEndRange = rangePerCodelet2269 * (codeletID + 1);
        if (remainderRange2269 != 0) {
            if (codeletID < (uint32_t)remainderRange2269) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange2269;
                tempEndRange += remainderRange2269;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration2269;
        tempEndRange = tempEndRange * 1 + minIteration2269;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration2269 < lastIteration2269) {
            (this->inputsTPParent->i_darts2269[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts2269[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration2269;
        }
    }
    return isThereNewIteration;
}
void TP2269::_checkInCodelets2270::fire(void)
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
    this->inputsTPParent->eta_darts2269[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->eta_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->iglob_darts2269[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jglob_darts2269[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->xi_darts2269[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->xi_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->zeta_darts2269[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->zeta_darts2199[this->getID()]);

    /*printing node 2270: ForStmt*/
    /*var: eta*/
    /*var: i*/
    /*var: iglob*/
    /*var: j*/
    /*var: jglob*/
    /*var: k*/
    /*var: m*/
    /*var: xi*/
    /*var: zeta*/
    double** eta = &(this->inputsTPParent->eta_darts2269[this->getLocalID()]);
    (void)eta /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts2269[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iglob = &(this->inputsTPParent->iglob_darts2269[this->getLocalID()]);
    (void)iglob /*OMP_SHARED_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts2269[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jglob = &(this->inputsTPParent->jglob_darts2269[this->getLocalID()]);
    (void)jglob /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts2269[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts2269[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** xi = &(this->inputsTPParent->xi_darts2269[this->getLocalID()]);
    (void)xi /*OMP_SHARED_PRIVATE*/;
    double** zeta = &(this->inputsTPParent->zeta_darts2269[this->getLocalID()]);
    (void)zeta /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations2269(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets2269[0].decDep();
        return;
    }
    for (int i_darts_counter_temp2269 = (*i); i_darts_counter_temp2269 < endRange
         && i_darts_counter_temp2269 < this->inputsTPParent->lastIteration2269;
         i_darts_counter_temp2269++) {
        {
            (*(*iglob)) = (i_darts_counter_temp2269);
            (*(*xi)) = ((double)((*(*iglob)))) / (nx0 - 1);
            {
                /*Loop's init*/
                (*j) = 0;
                int j_darts_counter_temp2269 = (*j);
                for (; j_darts_counter_temp2269 < ny; j_darts_counter_temp2269++) {
                    (*(*jglob)) = j_darts_counter_temp2269;
                    (*(*eta)) = ((double)((*(*jglob)))) / (ny0 - 1);
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp2269 = (*k);
                        for (; k_darts_counter_temp2269 < nz; k_darts_counter_temp2269++) {
                            (*(*zeta)) = ((double)(k_darts_counter_temp2269)) / (nz - 1);
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp2269 = (*m);
                                for (; m_darts_counter_temp2269 < 5; m_darts_counter_temp2269++) {
                                    rsd[(i_darts_counter_temp2269)][j_darts_counter_temp2269]
                                       [k_darts_counter_temp2269][m_darts_counter_temp2269]
                                        = ce[m_darts_counter_temp2269][0]
                                        + ce[m_darts_counter_temp2269][1] * (*(*xi))
                                        + ce[m_darts_counter_temp2269][2] * (*(*eta))
                                        + ce[m_darts_counter_temp2269][3] * (*(*zeta))
                                        + ce[m_darts_counter_temp2269][4] * (*(*xi)) * (*(*xi))
                                        + ce[m_darts_counter_temp2269][5] * (*(*eta)) * (*(*eta))
                                        + ce[m_darts_counter_temp2269][6] * (*(*zeta)) * (*(*zeta))
                                        + ce[m_darts_counter_temp2269][7] * (*(*xi)) * (*(*xi))
                                            * (*(*xi))
                                        + ce[m_darts_counter_temp2269][8] * (*(*eta)) * (*(*eta))
                                            * (*(*eta))
                                        + ce[m_darts_counter_temp2269][9] * (*(*zeta)) * (*(*zeta))
                                            * (*(*zeta))
                                        + ce[m_darts_counter_temp2269][10] * (*(*xi)) * (*(*xi))
                                            * (*(*xi)) * (*(*xi))
                                        + ce[m_darts_counter_temp2269][11] * (*(*eta)) * (*(*eta))
                                            * (*(*eta)) * (*(*eta))
                                        + ce[m_darts_counter_temp2269][12] * (*(*zeta)) * (*(*zeta))
                                            * (*(*zeta)) * (*(*zeta));
                                }
                                (*m) = m_darts_counter_temp2269;
                            }
                        }
                        (*k) = k_darts_counter_temp2269;
                    }
                }
                (*j) = j_darts_counter_temp2269;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets2269[0].decDep();
}
TP2269::TP2269(int in_numThreads, int in_mainCodeletID, TP2199* in_TPParent, int in_initIteration,
    int in_lastIteration, TP2269** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , eta_darts2269(new double*[this->numThreads])
    , i_darts2269(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iglob_darts2269(new int*[this->numThreads])
    , j_darts2269(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jglob_darts2269(new int*[this->numThreads])
    , k_darts2269(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts2269(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , xi_darts2269(new double*[this->numThreads])
    , zeta_darts2269(new double*[this->numThreads])
    , initIteration2269(in_initIteration)
    , lastIteration2269(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets2269(new _barrierCodelets2269[1])
    , checkInCodelets2270(new _checkInCodelets2270[this->numThreads])
{
    /*Initialize the loop parameters*/
    range2269 = abs(lastIteration2269 - initIteration2269) / 1;
    rangePerCodelet2269 = range2269 / numThreads;
    minIteration2269 = min<int>(lastIteration2269, initIteration2269);
    remainderRange2269 = range2269 % numThreads;
    /*Initialize inputs and vars.*/
    this->eta_darts2269
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->iglob_darts2269 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jglob_darts2269 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->xi_darts2269
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->zeta_darts2269
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets2269[0] = _barrierCodelets2269(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets2270* checkInCodelets2270Ptr = (this->checkInCodelets2270);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets2270);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets2270Ptr) = _checkInCodelets2270(2, 1, this, codeletCounter);
#else
        (*checkInCodelets2270Ptr) = _checkInCodelets2270(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets2270Ptr).decDep();
        checkInCodelets2270Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP2269::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets2270[localID].setID(codeletID);
    this->checkInCodelets2270[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets2270[localID + this->baseNumThreads * i]
            = _checkInCodelets2270(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets2270[localID + this->baseNumThreads * i]
            = _checkInCodelets2270(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets2270[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets2270[localID + this->baseNumThreads * i].decDep();
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
TP2269::~TP2269()
{
    delete[] eta_darts2269;
    delete[] iglob_darts2269;
    delete[] jglob_darts2269;
    delete[] xi_darts2269;
    delete[] zeta_darts2269;
    delete[] barrierCodelets2269;
    delete[] checkInCodelets2270;
}
/*TP2404: OMPForDirective*/
void TP2404::_barrierCodelets2404::fire(void)
{
    TP2404* myTP = static_cast<TP2404*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets2404[0].decDep();
}
bool TP2404::requestNewRangeIterations2404(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet2404 * codeletID;
        int tempEndRange = rangePerCodelet2404 * (codeletID + 1);
        if (remainderRange2404 != 0) {
            if (codeletID < (uint32_t)remainderRange2404) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange2404;
                tempEndRange += remainderRange2404;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration2404;
        tempEndRange = tempEndRange * 1 + minIteration2404;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration2404 < lastIteration2404) {
            (this->inputsTPParent->i_darts2404[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts2404[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration2404;
        }
    }
    return isThereNewIteration;
}
void TP2404::_checkInCodelets2405::fire(void)
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
    this->inputsTPParent->L1_darts2404[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L1_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L2_darts2404[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->q_darts2404[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21_darts2404[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21_darts2199[this->getID()]);

    /*printing node 2405: ForStmt*/
    /*var: L1*/
    /*var: L2*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: q*/
    /*var: u21*/
    int* i = &(this->inputsTPParent->i_darts2404[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts2404[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts2404[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts2404[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** u21 = &(this->inputsTPParent->u21_darts2404[this->getLocalID()]);
    (void)u21 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations2404(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets2404[0].decDep();
        return;
    }
    for (int i_darts_counter_temp2404 = (*i); i_darts_counter_temp2404 <= endRange
         && i_darts_counter_temp2404 <= this->inputsTPParent->lastIteration2404;
         i_darts_counter_temp2404++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp2404 = (*j);
                for (; j_darts_counter_temp2404 <= jend; j_darts_counter_temp2404++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp2404 = (*k);
                        for (; k_darts_counter_temp2404 < nz - 1; k_darts_counter_temp2404++) {
                            flux[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                [k_darts_counter_temp2404][0]
                                = rsd[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                     [k_darts_counter_temp2404][1];
                            (*(*u21)) = rsd[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                           [k_darts_counter_temp2404][1]
                                / rsd[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                     [k_darts_counter_temp2404][0];
                            (*(*q)) = 0.5
                                * (rsd[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                      [k_darts_counter_temp2404][1]
                                        * rsd[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                             [k_darts_counter_temp2404][1]
                                    + rsd[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                         [k_darts_counter_temp2404][2]
                                        * rsd[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                             [k_darts_counter_temp2404][2]
                                    + rsd[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                         [k_darts_counter_temp2404][3]
                                        * rsd[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                             [k_darts_counter_temp2404][3])
                                / rsd[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                     [k_darts_counter_temp2404][0];
                            flux[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                [k_darts_counter_temp2404][1]
                                = rsd[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                     [k_darts_counter_temp2404][1]
                                    * (*(*u21))
                                + 0.40000000000000002
                                    * (rsd[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                          [k_darts_counter_temp2404][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                [k_darts_counter_temp2404][2]
                                = rsd[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                     [k_darts_counter_temp2404][2]
                                * (*(*u21));
                            flux[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                [k_darts_counter_temp2404][3]
                                = rsd[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                     [k_darts_counter_temp2404][3]
                                * (*(*u21));
                            flux[(i_darts_counter_temp2404)][j_darts_counter_temp2404]
                                [k_darts_counter_temp2404][4]
                                = (1.3999999999999999
                                          * rsd[(i_darts_counter_temp2404)]
                                               [j_darts_counter_temp2404][k_darts_counter_temp2404]
                                               [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u21));
                        }
                        (*k) = k_darts_counter_temp2404;
                    }
                }
                (*j) = j_darts_counter_temp2404;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets2404[0].decDep();
}
TP2404::TP2404(int in_numThreads, int in_mainCodeletID, TP2199* in_TPParent, int in_initIteration,
    int in_lastIteration, TP2404** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts2404(new int*[this->numThreads])
    , L2_darts2404(new int*[this->numThreads])
    , i_darts2404(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts2404(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts2404(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts2404(new double*[this->numThreads])
    , u21_darts2404(new double*[this->numThreads])
    , initIteration2404(in_initIteration)
    , lastIteration2404(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets2404(new _barrierCodelets2404[1])
    , checkInCodelets2405(new _checkInCodelets2405[this->numThreads])
{
    /*Initialize the loop parameters*/
    range2404 = abs(lastIteration2404 - initIteration2404) / 1;
    rangePerCodelet2404 = range2404 / numThreads;
    minIteration2404 = min<int>(lastIteration2404, initIteration2404);
    remainderRange2404 = range2404 % numThreads;
    /*Initialize inputs and vars.*/
    this->L1_darts2404 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->L2_darts2404 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->q_darts2404 = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21_darts2404
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets2404[0] = _barrierCodelets2404(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets2405* checkInCodelets2405Ptr = (this->checkInCodelets2405);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets2405);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets2405Ptr) = _checkInCodelets2405(2, 1, this, codeletCounter);
#else
        (*checkInCodelets2405Ptr) = _checkInCodelets2405(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets2405Ptr).decDep();
        checkInCodelets2405Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP2404::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets2405[localID].setID(codeletID);
    this->checkInCodelets2405[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets2405[localID + this->baseNumThreads * i]
            = _checkInCodelets2405(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets2405[localID + this->baseNumThreads * i]
            = _checkInCodelets2405(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets2405[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets2405[localID + this->baseNumThreads * i].decDep();
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
TP2404::~TP2404()
{
    delete[] L1_darts2404;
    delete[] L2_darts2404;
    delete[] q_darts2404;
    delete[] u21_darts2404;
    delete[] barrierCodelets2404;
    delete[] checkInCodelets2405;
}
/*TP2553: OMPForDirective*/
void TP2553::_barrierCodelets2553::fire(void)
{
    TP2553* myTP = static_cast<TP2553*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets2553[0].decDep();
}
bool TP2553::requestNewRangeIterations2553(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet2553 * codeletID;
        int tempEndRange = rangePerCodelet2553 * (codeletID + 1);
        if (remainderRange2553 != 0) {
            if (codeletID < (uint32_t)remainderRange2553) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange2553;
                tempEndRange += remainderRange2553;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration2553;
        tempEndRange = tempEndRange * 1 + minIteration2553;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration2553 < lastIteration2553) {
            (this->inputsTPParent->j_darts2553[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->j_darts2553[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration2553;
        }
    }
    return isThereNewIteration;
}
void TP2553::_checkInCodelets2554::fire(void)
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
    this->inputsTPParent->L2_darts2553[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->dsspm_darts2553[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->dsspm_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->iend1_darts2553[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iend1_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ist1_darts2553[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->ist1_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts2553[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21i_darts2553[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21i_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21im1_darts2553[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21im1_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31i_darts2553[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31i_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31im1_darts2553[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31im1_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41i_darts2553[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41i_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41im1_darts2553[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41im1_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51i_darts2553[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51i_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51im1_darts2553[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51im1_darts2199[this->getID()]);

    /*printing node 2554: ForStmt*/
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
    int** L2 = &(this->inputsTPParent->L2_darts2553[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    double** dsspm = &(this->inputsTPParent->dsspm_darts2553[this->getLocalID()]);
    (void)dsspm /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts2553[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iend1 = &(this->inputsTPParent->iend1_darts2553[this->getLocalID()]);
    (void)iend1 /*OMP_SHARED_PRIVATE*/;
    int** ist1 = &(this->inputsTPParent->ist1_darts2553[this->getLocalID()]);
    (void)ist1 /*OMP_SHARED_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts2553[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts2553[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts2553[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts2553[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21i = &(this->inputsTPParent->u21i_darts2553[this->getLocalID()]);
    (void)u21i /*OMP_SHARED_PRIVATE*/;
    double** u21im1 = &(this->inputsTPParent->u21im1_darts2553[this->getLocalID()]);
    (void)u21im1 /*OMP_SHARED_PRIVATE*/;
    double** u31i = &(this->inputsTPParent->u31i_darts2553[this->getLocalID()]);
    (void)u31i /*OMP_SHARED_PRIVATE*/;
    double** u31im1 = &(this->inputsTPParent->u31im1_darts2553[this->getLocalID()]);
    (void)u31im1 /*OMP_SHARED_PRIVATE*/;
    double** u41i = &(this->inputsTPParent->u41i_darts2553[this->getLocalID()]);
    (void)u41i /*OMP_SHARED_PRIVATE*/;
    double** u41im1 = &(this->inputsTPParent->u41im1_darts2553[this->getLocalID()]);
    (void)u41im1 /*OMP_SHARED_PRIVATE*/;
    double** u51i = &(this->inputsTPParent->u51i_darts2553[this->getLocalID()]);
    (void)u51i /*OMP_SHARED_PRIVATE*/;
    double** u51im1 = &(this->inputsTPParent->u51im1_darts2553[this->getLocalID()]);
    (void)u51im1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations2553(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets2553[0].decDep();
        return;
    }
    for (int j_darts_counter_temp2553 = (*j); j_darts_counter_temp2553 <= endRange
         && j_darts_counter_temp2553 <= this->inputsTPParent->lastIteration2553;
         j_darts_counter_temp2553++) {
        {
            {
                /*Loop's init*/
                (*k) = 1;
                int k_darts_counter_temp2553 = (*k);
                for (; k_darts_counter_temp2553 <= nz - 2; k_darts_counter_temp2553++) {
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp2553 = (*i);
                        for (; i_darts_counter_temp2553 <= iend; i_darts_counter_temp2553++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp2553 = (*m);
                                for (; m_darts_counter_temp2553 < 5; m_darts_counter_temp2553++) {
                                    frct[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                        [k_darts_counter_temp2553][m_darts_counter_temp2553]
                                        = frct[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                              [k_darts_counter_temp2553][m_darts_counter_temp2553]
                                        - tx2
                                            * (flux[i_darts_counter_temp2553 + 1]
                                                   [(j_darts_counter_temp2553)]
                                                   [k_darts_counter_temp2553]
                                                   [m_darts_counter_temp2553]
                                                - flux[i_darts_counter_temp2553 - 1]
                                                      [(j_darts_counter_temp2553)]
                                                      [k_darts_counter_temp2553]
                                                      [m_darts_counter_temp2553]);
                                }
                                (*m) = m_darts_counter_temp2553;
                            }
                        }
                        (*i) = i_darts_counter_temp2553;
                    }
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp2553 = (*i);
                        for (; i_darts_counter_temp2553 <= (*(*L2)); i_darts_counter_temp2553++) {
                            (*(*tmp)) = 1.
                                / rsd[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                     [k_darts_counter_temp2553][0];
                            (*(*u21i)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                     [k_darts_counter_temp2553][1];
                            (*(*u31i)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                     [k_darts_counter_temp2553][2];
                            (*(*u41i)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                     [k_darts_counter_temp2553][3];
                            (*(*u51i)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                     [k_darts_counter_temp2553][4];
                            (*(*tmp)) = 1.
                                / rsd[i_darts_counter_temp2553 - 1][(j_darts_counter_temp2553)]
                                     [k_darts_counter_temp2553][0];
                            (*(*u21im1)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2553 - 1][(j_darts_counter_temp2553)]
                                     [k_darts_counter_temp2553][1];
                            (*(*u31im1)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2553 - 1][(j_darts_counter_temp2553)]
                                     [k_darts_counter_temp2553][2];
                            (*(*u41im1)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2553 - 1][(j_darts_counter_temp2553)]
                                     [k_darts_counter_temp2553][3];
                            (*(*u51im1)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2553 - 1][(j_darts_counter_temp2553)]
                                     [k_darts_counter_temp2553][4];
                            flux[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                [k_darts_counter_temp2553][1]
                                = (4. / 3.) * tx3 * ((*(*u21i)) - (*(*u21im1)));
                            flux[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                [k_darts_counter_temp2553][2]
                                = tx3 * ((*(*u31i)) - (*(*u31im1)));
                            flux[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                [k_darts_counter_temp2553][3]
                                = tx3 * ((*(*u41i)) - (*(*u41im1)));
                            flux[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                [k_darts_counter_temp2553][4]
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
                        (*i) = i_darts_counter_temp2553;
                    }
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp2553 = (*i);
                        for (; i_darts_counter_temp2553 <= iend; i_darts_counter_temp2553++) {
                            frct[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                [k_darts_counter_temp2553][0]
                                = frct[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                      [k_darts_counter_temp2553][0]
                                + dx1 * tx1
                                    * (rsd[i_darts_counter_temp2553 - 1][(j_darts_counter_temp2553)]
                                          [k_darts_counter_temp2553][0]
                                        - 2.
                                            * rsd[i_darts_counter_temp2553]
                                                 [(j_darts_counter_temp2553)]
                                                 [k_darts_counter_temp2553][0]
                                        + rsd[i_darts_counter_temp2553 + 1]
                                             [(j_darts_counter_temp2553)][k_darts_counter_temp2553]
                                             [0]);
                            frct[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                [k_darts_counter_temp2553][1]
                                = frct[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                      [k_darts_counter_temp2553][1]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp2553 + 1]
                                           [(j_darts_counter_temp2553)][k_darts_counter_temp2553][1]
                                        - flux[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                              [k_darts_counter_temp2553][1])
                                + dx2 * tx1
                                    * (rsd[i_darts_counter_temp2553 - 1][(j_darts_counter_temp2553)]
                                          [k_darts_counter_temp2553][1]
                                        - 2.
                                            * rsd[i_darts_counter_temp2553]
                                                 [(j_darts_counter_temp2553)]
                                                 [k_darts_counter_temp2553][1]
                                        + rsd[i_darts_counter_temp2553 + 1]
                                             [(j_darts_counter_temp2553)][k_darts_counter_temp2553]
                                             [1]);
                            frct[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                [k_darts_counter_temp2553][2]
                                = frct[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                      [k_darts_counter_temp2553][2]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp2553 + 1]
                                           [(j_darts_counter_temp2553)][k_darts_counter_temp2553][2]
                                        - flux[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                              [k_darts_counter_temp2553][2])
                                + dx3 * tx1
                                    * (rsd[i_darts_counter_temp2553 - 1][(j_darts_counter_temp2553)]
                                          [k_darts_counter_temp2553][2]
                                        - 2.
                                            * rsd[i_darts_counter_temp2553]
                                                 [(j_darts_counter_temp2553)]
                                                 [k_darts_counter_temp2553][2]
                                        + rsd[i_darts_counter_temp2553 + 1]
                                             [(j_darts_counter_temp2553)][k_darts_counter_temp2553]
                                             [2]);
                            frct[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                [k_darts_counter_temp2553][3]
                                = frct[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                      [k_darts_counter_temp2553][3]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp2553 + 1]
                                           [(j_darts_counter_temp2553)][k_darts_counter_temp2553][3]
                                        - flux[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                              [k_darts_counter_temp2553][3])
                                + dx4 * tx1
                                    * (rsd[i_darts_counter_temp2553 - 1][(j_darts_counter_temp2553)]
                                          [k_darts_counter_temp2553][3]
                                        - 2.
                                            * rsd[i_darts_counter_temp2553]
                                                 [(j_darts_counter_temp2553)]
                                                 [k_darts_counter_temp2553][3]
                                        + rsd[i_darts_counter_temp2553 + 1]
                                             [(j_darts_counter_temp2553)][k_darts_counter_temp2553]
                                             [3]);
                            frct[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                [k_darts_counter_temp2553][4]
                                = frct[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                      [k_darts_counter_temp2553][4]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp2553 + 1]
                                           [(j_darts_counter_temp2553)][k_darts_counter_temp2553][4]
                                        - flux[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                              [k_darts_counter_temp2553][4])
                                + dx5 * tx1
                                    * (rsd[i_darts_counter_temp2553 - 1][(j_darts_counter_temp2553)]
                                          [k_darts_counter_temp2553][4]
                                        - 2.
                                            * rsd[i_darts_counter_temp2553]
                                                 [(j_darts_counter_temp2553)]
                                                 [k_darts_counter_temp2553][4]
                                        + rsd[i_darts_counter_temp2553 + 1]
                                             [(j_darts_counter_temp2553)][k_darts_counter_temp2553]
                                             [4]);
                        }
                        (*i) = i_darts_counter_temp2553;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp2553 = (*m);
                        for (; m_darts_counter_temp2553 < 5; m_darts_counter_temp2553++) {
                            frct[1][(j_darts_counter_temp2553)][k_darts_counter_temp2553]
                                [m_darts_counter_temp2553]
                                = frct[1][(j_darts_counter_temp2553)][k_darts_counter_temp2553]
                                      [m_darts_counter_temp2553]
                                - (*(*dsspm))
                                    * (+5.
                                            * rsd[1][(j_darts_counter_temp2553)]
                                                 [k_darts_counter_temp2553]
                                                 [m_darts_counter_temp2553]
                                        - 4.
                                            * rsd[2][(j_darts_counter_temp2553)]
                                                 [k_darts_counter_temp2553]
                                                 [m_darts_counter_temp2553]
                                        + rsd[3][(j_darts_counter_temp2553)]
                                             [k_darts_counter_temp2553][m_darts_counter_temp2553]);
                            frct[2][(j_darts_counter_temp2553)][k_darts_counter_temp2553]
                                [m_darts_counter_temp2553]
                                = frct[2][(j_darts_counter_temp2553)][k_darts_counter_temp2553]
                                      [m_darts_counter_temp2553]
                                - (*(*dsspm))
                                    * (-4.
                                            * rsd[1][(j_darts_counter_temp2553)]
                                                 [k_darts_counter_temp2553]
                                                 [m_darts_counter_temp2553]
                                        + 6.
                                            * rsd[2][(j_darts_counter_temp2553)]
                                                 [k_darts_counter_temp2553]
                                                 [m_darts_counter_temp2553]
                                        - 4.
                                            * rsd[3][(j_darts_counter_temp2553)]
                                                 [k_darts_counter_temp2553]
                                                 [m_darts_counter_temp2553]
                                        + rsd[4][(j_darts_counter_temp2553)]
                                             [k_darts_counter_temp2553][m_darts_counter_temp2553]);
                        }
                        (*m) = m_darts_counter_temp2553;
                    }
                    (*(*ist1)) = 3;
                    (*(*iend1)) = nx - 4;
                    {
                        /*Loop's init*/
                        (*i) = (*(*ist1));
                        int i_darts_counter_temp2553 = (*i);
                        for (; i_darts_counter_temp2553 <= (*(*iend1));
                             i_darts_counter_temp2553++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp2553 = (*m);
                                for (; m_darts_counter_temp2553 < 5; m_darts_counter_temp2553++) {
                                    frct[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                        [k_darts_counter_temp2553][m_darts_counter_temp2553]
                                        = frct[i_darts_counter_temp2553][(j_darts_counter_temp2553)]
                                              [k_darts_counter_temp2553][m_darts_counter_temp2553]
                                        - (*(*dsspm))
                                            * (rsd[i_darts_counter_temp2553 - 2]
                                                  [(j_darts_counter_temp2553)]
                                                  [k_darts_counter_temp2553]
                                                  [m_darts_counter_temp2553]
                                                - 4.
                                                    * rsd[i_darts_counter_temp2553 - 1]
                                                         [(j_darts_counter_temp2553)]
                                                         [k_darts_counter_temp2553]
                                                         [m_darts_counter_temp2553]
                                                + 6.
                                                    * rsd[i_darts_counter_temp2553]
                                                         [(j_darts_counter_temp2553)]
                                                         [k_darts_counter_temp2553]
                                                         [m_darts_counter_temp2553]
                                                - 4.
                                                    * rsd[i_darts_counter_temp2553 + 1]
                                                         [(j_darts_counter_temp2553)]
                                                         [k_darts_counter_temp2553]
                                                         [m_darts_counter_temp2553]
                                                + rsd[i_darts_counter_temp2553 + 2]
                                                     [(j_darts_counter_temp2553)]
                                                     [k_darts_counter_temp2553]
                                                     [m_darts_counter_temp2553]);
                                }
                                (*m) = m_darts_counter_temp2553;
                            }
                        }
                        (*i) = i_darts_counter_temp2553;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp2553 = (*m);
                        for (; m_darts_counter_temp2553 < 5; m_darts_counter_temp2553++) {
                            frct[nx - 3][(j_darts_counter_temp2553)][k_darts_counter_temp2553]
                                [m_darts_counter_temp2553]
                                = frct[nx - 3][(j_darts_counter_temp2553)][k_darts_counter_temp2553]
                                      [m_darts_counter_temp2553]
                                - (*(*dsspm))
                                    * (rsd[nx - 5][(j_darts_counter_temp2553)]
                                          [k_darts_counter_temp2553][m_darts_counter_temp2553]
                                        - 4.
                                            * rsd[nx - 4][(j_darts_counter_temp2553)]
                                                 [k_darts_counter_temp2553]
                                                 [m_darts_counter_temp2553]
                                        + 6.
                                            * rsd[nx - 3][(j_darts_counter_temp2553)]
                                                 [k_darts_counter_temp2553]
                                                 [m_darts_counter_temp2553]
                                        - 4.
                                            * rsd[nx - 2][(j_darts_counter_temp2553)]
                                                 [k_darts_counter_temp2553]
                                                 [m_darts_counter_temp2553]);
                            frct[nx - 2][(j_darts_counter_temp2553)][k_darts_counter_temp2553]
                                [m_darts_counter_temp2553]
                                = frct[nx - 2][(j_darts_counter_temp2553)][k_darts_counter_temp2553]
                                      [m_darts_counter_temp2553]
                                - (*(*dsspm))
                                    * (rsd[nx - 4][(j_darts_counter_temp2553)]
                                          [k_darts_counter_temp2553][m_darts_counter_temp2553]
                                        - 4.
                                            * rsd[nx - 3][(j_darts_counter_temp2553)]
                                                 [k_darts_counter_temp2553]
                                                 [m_darts_counter_temp2553]
                                        + 5.
                                            * rsd[nx - 2][(j_darts_counter_temp2553)]
                                                 [k_darts_counter_temp2553]
                                                 [m_darts_counter_temp2553]);
                        }
                        (*m) = m_darts_counter_temp2553;
                    }
                }
                (*k) = k_darts_counter_temp2553;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets2553[0].decDep();
}
TP2553::TP2553(int in_numThreads, int in_mainCodeletID, TP2199* in_TPParent, int in_initIteration,
    int in_lastIteration, TP2553** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L2_darts2553(new int*[this->numThreads])
    , dsspm_darts2553(new double*[this->numThreads])
    , i_darts2553(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iend1_darts2553(new int*[this->numThreads])
    , ist1_darts2553(new int*[this->numThreads])
    , j_darts2553(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts2553(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts2553(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts2553(new double*[this->numThreads])
    , u21i_darts2553(new double*[this->numThreads])
    , u21im1_darts2553(new double*[this->numThreads])
    , u31i_darts2553(new double*[this->numThreads])
    , u31im1_darts2553(new double*[this->numThreads])
    , u41i_darts2553(new double*[this->numThreads])
    , u41im1_darts2553(new double*[this->numThreads])
    , u51i_darts2553(new double*[this->numThreads])
    , u51im1_darts2553(new double*[this->numThreads])
    , initIteration2553(in_initIteration)
    , lastIteration2553(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets2553(new _barrierCodelets2553[1])
    , checkInCodelets2554(new _checkInCodelets2554[this->numThreads])
{
    /*Initialize the loop parameters*/
    range2553 = abs(lastIteration2553 - initIteration2553) / 1;
    rangePerCodelet2553 = range2553 / numThreads;
    minIteration2553 = min<int>(lastIteration2553, initIteration2553);
    remainderRange2553 = range2553 % numThreads;
    /*Initialize inputs and vars.*/
    this->L2_darts2553 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->dsspm_darts2553
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->iend1_darts2553 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ist1_darts2553 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts2553
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21i_darts2553
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21im1_darts2553
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31i_darts2553
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31im1_darts2553
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41i_darts2553
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41im1_darts2553
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51i_darts2553
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51im1_darts2553
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets2553[0] = _barrierCodelets2553(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets2554* checkInCodelets2554Ptr = (this->checkInCodelets2554);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets2554);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets2554Ptr) = _checkInCodelets2554(2, 1, this, codeletCounter);
#else
        (*checkInCodelets2554Ptr) = _checkInCodelets2554(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets2554Ptr).decDep();
        checkInCodelets2554Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP2553::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets2554[localID].setID(codeletID);
    this->checkInCodelets2554[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets2554[localID + this->baseNumThreads * i]
            = _checkInCodelets2554(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets2554[localID + this->baseNumThreads * i]
            = _checkInCodelets2554(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets2554[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets2554[localID + this->baseNumThreads * i].decDep();
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
TP2553::~TP2553()
{
    delete[] L2_darts2553;
    delete[] dsspm_darts2553;
    delete[] iend1_darts2553;
    delete[] ist1_darts2553;
    delete[] tmp_darts2553;
    delete[] u21i_darts2553;
    delete[] u21im1_darts2553;
    delete[] u31i_darts2553;
    delete[] u31im1_darts2553;
    delete[] u41i_darts2553;
    delete[] u41im1_darts2553;
    delete[] u51i_darts2553;
    delete[] u51im1_darts2553;
    delete[] barrierCodelets2553;
    delete[] checkInCodelets2554;
}
/*TP3191: OMPForDirective*/
void TP3191::_barrierCodelets3191::fire(void)
{
    TP3191* myTP = static_cast<TP3191*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets3191[0].decDep();
}
bool TP3191::requestNewRangeIterations3191(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet3191 * codeletID;
        int tempEndRange = rangePerCodelet3191 * (codeletID + 1);
        if (remainderRange3191 != 0) {
            if (codeletID < (uint32_t)remainderRange3191) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange3191;
                tempEndRange += remainderRange3191;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration3191;
        tempEndRange = tempEndRange * 1 + minIteration3191;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration3191 < lastIteration3191) {
            (this->inputsTPParent->i_darts3191[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts3191[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration3191;
        }
    }
    return isThereNewIteration;
}
void TP3191::_checkInCodelets3192::fire(void)
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
    this->inputsTPParent->L1_darts3191[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L1_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L2_darts3191[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->q_darts3191[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31_darts3191[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31_darts2199[this->getID()]);

    /*printing node 3192: ForStmt*/
    /*var: L1*/
    /*var: L2*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: q*/
    /*var: u31*/
    int** L1 = &(this->inputsTPParent->L1_darts3191[this->getLocalID()]);
    (void)L1 /*OMP_SHARED_PRIVATE*/;
    int** L2 = &(this->inputsTPParent->L2_darts3191[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts3191[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts3191[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts3191[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts3191[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** u31 = &(this->inputsTPParent->u31_darts3191[this->getLocalID()]);
    (void)u31 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations3191(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets3191[0].decDep();
        return;
    }
    for (int i_darts_counter_temp3191 = (*i); i_darts_counter_temp3191 <= endRange
         && i_darts_counter_temp3191 <= this->inputsTPParent->lastIteration3191;
         i_darts_counter_temp3191++) {
        {
            {
                /*Loop's init*/
                (*j) = (*(*L1));
                int j_darts_counter_temp3191 = (*j);
                for (; j_darts_counter_temp3191 <= (*(*L2)); j_darts_counter_temp3191++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp3191 = (*k);
                        for (; k_darts_counter_temp3191 <= nz - 2; k_darts_counter_temp3191++) {
                            flux[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                [k_darts_counter_temp3191][0]
                                = rsd[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                     [k_darts_counter_temp3191][2];
                            (*(*u31)) = rsd[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                           [k_darts_counter_temp3191][2]
                                / rsd[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                     [k_darts_counter_temp3191][0];
                            (*(*q)) = 0.5
                                * (rsd[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                      [k_darts_counter_temp3191][1]
                                        * rsd[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                             [k_darts_counter_temp3191][1]
                                    + rsd[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                         [k_darts_counter_temp3191][2]
                                        * rsd[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                             [k_darts_counter_temp3191][2]
                                    + rsd[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                         [k_darts_counter_temp3191][3]
                                        * rsd[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                             [k_darts_counter_temp3191][3])
                                / rsd[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                     [k_darts_counter_temp3191][0];
                            flux[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                [k_darts_counter_temp3191][1]
                                = rsd[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                     [k_darts_counter_temp3191][1]
                                * (*(*u31));
                            flux[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                [k_darts_counter_temp3191][2]
                                = rsd[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                     [k_darts_counter_temp3191][2]
                                    * (*(*u31))
                                + 0.40000000000000002
                                    * (rsd[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                          [k_darts_counter_temp3191][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                [k_darts_counter_temp3191][3]
                                = rsd[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                     [k_darts_counter_temp3191][3]
                                * (*(*u31));
                            flux[(i_darts_counter_temp3191)][j_darts_counter_temp3191]
                                [k_darts_counter_temp3191][4]
                                = (1.3999999999999999
                                          * rsd[(i_darts_counter_temp3191)]
                                               [j_darts_counter_temp3191][k_darts_counter_temp3191]
                                               [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u31));
                        }
                        (*k) = k_darts_counter_temp3191;
                    }
                }
                (*j) = j_darts_counter_temp3191;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets3191[0].decDep();
}
TP3191::TP3191(int in_numThreads, int in_mainCodeletID, TP2199* in_TPParent, int in_initIteration,
    int in_lastIteration, TP3191** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts3191(new int*[this->numThreads])
    , L2_darts3191(new int*[this->numThreads])
    , i_darts3191(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts3191(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts3191(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts3191(new double*[this->numThreads])
    , u31_darts3191(new double*[this->numThreads])
    , initIteration3191(in_initIteration)
    , lastIteration3191(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets3191(new _barrierCodelets3191[1])
    , checkInCodelets3192(new _checkInCodelets3192[this->numThreads])
{
    /*Initialize the loop parameters*/
    range3191 = abs(lastIteration3191 - initIteration3191) / 1;
    rangePerCodelet3191 = range3191 / numThreads;
    minIteration3191 = min<int>(lastIteration3191, initIteration3191);
    remainderRange3191 = range3191 % numThreads;
    /*Initialize inputs and vars.*/
    this->L1_darts3191 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->L2_darts3191 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->q_darts3191 = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31_darts3191
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets3191[0] = _barrierCodelets3191(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets3192* checkInCodelets3192Ptr = (this->checkInCodelets3192);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets3192);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets3192Ptr) = _checkInCodelets3192(2, 1, this, codeletCounter);
#else
        (*checkInCodelets3192Ptr) = _checkInCodelets3192(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets3192Ptr).decDep();
        checkInCodelets3192Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP3191::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets3192[localID].setID(codeletID);
    this->checkInCodelets3192[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets3192[localID + this->baseNumThreads * i]
            = _checkInCodelets3192(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets3192[localID + this->baseNumThreads * i]
            = _checkInCodelets3192(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets3192[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets3192[localID + this->baseNumThreads * i].decDep();
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
TP3191::~TP3191()
{
    delete[] L1_darts3191;
    delete[] L2_darts3191;
    delete[] q_darts3191;
    delete[] u31_darts3191;
    delete[] barrierCodelets3191;
    delete[] checkInCodelets3192;
}
/*TP3340: OMPForDirective*/
void TP3340::_barrierCodelets3340::fire(void)
{
    TP3340* myTP = static_cast<TP3340*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets3340[0].decDep();
}
bool TP3340::requestNewRangeIterations3340(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet3340 * codeletID;
        int tempEndRange = rangePerCodelet3340 * (codeletID + 1);
        if (remainderRange3340 != 0) {
            if (codeletID < (uint32_t)remainderRange3340) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange3340;
                tempEndRange += remainderRange3340;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration3340;
        tempEndRange = tempEndRange * 1 + minIteration3340;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration3340 < lastIteration3340) {
            (this->inputsTPParent->i_darts3340[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts3340[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration3340;
        }
    }
    return isThereNewIteration;
}
void TP3340::_checkInCodelets3341::fire(void)
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
    this->inputsTPParent->L2_darts3340[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->dsspm_darts3340[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->dsspm_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jend1_darts3340[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jend1_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jst1_darts3340[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jst1_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts3340[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21j_darts3340[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21j_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21jm1_darts3340[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21jm1_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31j_darts3340[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31j_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31jm1_darts3340[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31jm1_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41j_darts3340[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41j_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41jm1_darts3340[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41jm1_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51j_darts3340[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51j_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51jm1_darts3340[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51jm1_darts2199[this->getID()]);

    /*printing node 3341: ForStmt*/
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
    int** L2 = &(this->inputsTPParent->L2_darts3340[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    double** dsspm = &(this->inputsTPParent->dsspm_darts3340[this->getLocalID()]);
    (void)dsspm /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts3340[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts3340[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jend1 = &(this->inputsTPParent->jend1_darts3340[this->getLocalID()]);
    (void)jend1 /*OMP_SHARED_PRIVATE*/;
    int** jst1 = &(this->inputsTPParent->jst1_darts3340[this->getLocalID()]);
    (void)jst1 /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts3340[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts3340[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts3340[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21j = &(this->inputsTPParent->u21j_darts3340[this->getLocalID()]);
    (void)u21j /*OMP_SHARED_PRIVATE*/;
    double** u21jm1 = &(this->inputsTPParent->u21jm1_darts3340[this->getLocalID()]);
    (void)u21jm1 /*OMP_SHARED_PRIVATE*/;
    double** u31j = &(this->inputsTPParent->u31j_darts3340[this->getLocalID()]);
    (void)u31j /*OMP_SHARED_PRIVATE*/;
    double** u31jm1 = &(this->inputsTPParent->u31jm1_darts3340[this->getLocalID()]);
    (void)u31jm1 /*OMP_SHARED_PRIVATE*/;
    double** u41j = &(this->inputsTPParent->u41j_darts3340[this->getLocalID()]);
    (void)u41j /*OMP_SHARED_PRIVATE*/;
    double** u41jm1 = &(this->inputsTPParent->u41jm1_darts3340[this->getLocalID()]);
    (void)u41jm1 /*OMP_SHARED_PRIVATE*/;
    double** u51j = &(this->inputsTPParent->u51j_darts3340[this->getLocalID()]);
    (void)u51j /*OMP_SHARED_PRIVATE*/;
    double** u51jm1 = &(this->inputsTPParent->u51jm1_darts3340[this->getLocalID()]);
    (void)u51jm1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations3340(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets3340[0].decDep();
        return;
    }
    for (int i_darts_counter_temp3340 = (*i); i_darts_counter_temp3340 <= endRange
         && i_darts_counter_temp3340 <= this->inputsTPParent->lastIteration3340;
         i_darts_counter_temp3340++) {
        {
            {
                /*Loop's init*/
                (*k) = 1;
                int k_darts_counter_temp3340 = (*k);
                for (; k_darts_counter_temp3340 <= nz - 2; k_darts_counter_temp3340++) {
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp3340 = (*j);
                        for (; j_darts_counter_temp3340 <= jend; j_darts_counter_temp3340++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp3340 = (*m);
                                for (; m_darts_counter_temp3340 < 5; m_darts_counter_temp3340++) {
                                    frct[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                        [k_darts_counter_temp3340][m_darts_counter_temp3340]
                                        = frct[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                              [k_darts_counter_temp3340][m_darts_counter_temp3340]
                                        - ty2
                                            * (flux[(i_darts_counter_temp3340)]
                                                   [j_darts_counter_temp3340 + 1]
                                                   [k_darts_counter_temp3340]
                                                   [m_darts_counter_temp3340]
                                                - flux[(i_darts_counter_temp3340)]
                                                      [j_darts_counter_temp3340 - 1]
                                                      [k_darts_counter_temp3340]
                                                      [m_darts_counter_temp3340]);
                                }
                                (*m) = m_darts_counter_temp3340;
                            }
                        }
                        (*j) = j_darts_counter_temp3340;
                    }
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp3340 = (*j);
                        for (; j_darts_counter_temp3340 <= (*(*L2)); j_darts_counter_temp3340++) {
                            (*(*tmp)) = 1.
                                / rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                     [k_darts_counter_temp3340][0];
                            (*(*u21j)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                     [k_darts_counter_temp3340][1];
                            (*(*u31j)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                     [k_darts_counter_temp3340][2];
                            (*(*u41j)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                     [k_darts_counter_temp3340][3];
                            (*(*u51j)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                     [k_darts_counter_temp3340][4];
                            (*(*tmp)) = 1.
                                / rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340 - 1]
                                     [k_darts_counter_temp3340][0];
                            (*(*u21jm1)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340 - 1]
                                     [k_darts_counter_temp3340][1];
                            (*(*u31jm1)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340 - 1]
                                     [k_darts_counter_temp3340][2];
                            (*(*u41jm1)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340 - 1]
                                     [k_darts_counter_temp3340][3];
                            (*(*u51jm1)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340 - 1]
                                     [k_darts_counter_temp3340][4];
                            flux[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                [k_darts_counter_temp3340][1]
                                = ty3 * ((*(*u21j)) - (*(*u21jm1)));
                            flux[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                [k_darts_counter_temp3340][2]
                                = (4. / 3.) * ty3 * ((*(*u31j)) - (*(*u31jm1)));
                            flux[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                [k_darts_counter_temp3340][3]
                                = ty3 * ((*(*u41j)) - (*(*u41jm1)));
                            flux[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                [k_darts_counter_temp3340][4]
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
                        (*j) = j_darts_counter_temp3340;
                    }
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp3340 = (*j);
                        for (; j_darts_counter_temp3340 <= jend; j_darts_counter_temp3340++) {
                            frct[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                [k_darts_counter_temp3340][0]
                                = frct[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                      [k_darts_counter_temp3340][0]
                                + dy1 * ty1
                                    * (rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340 - 1]
                                          [k_darts_counter_temp3340][0]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3340)]
                                                 [j_darts_counter_temp3340]
                                                 [k_darts_counter_temp3340][0]
                                        + rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340
                                            + 1][k_darts_counter_temp3340][0]);
                            frct[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                [k_darts_counter_temp3340][1]
                                = frct[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                      [k_darts_counter_temp3340][1]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp3340)][j_darts_counter_temp3340
                                           + 1][k_darts_counter_temp3340][1]
                                        - flux[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                              [k_darts_counter_temp3340][1])
                                + dy2 * ty1
                                    * (rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340 - 1]
                                          [k_darts_counter_temp3340][1]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3340)]
                                                 [j_darts_counter_temp3340]
                                                 [k_darts_counter_temp3340][1]
                                        + rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340
                                            + 1][k_darts_counter_temp3340][1]);
                            frct[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                [k_darts_counter_temp3340][2]
                                = frct[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                      [k_darts_counter_temp3340][2]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp3340)][j_darts_counter_temp3340
                                           + 1][k_darts_counter_temp3340][2]
                                        - flux[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                              [k_darts_counter_temp3340][2])
                                + dy3 * ty1
                                    * (rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340 - 1]
                                          [k_darts_counter_temp3340][2]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3340)]
                                                 [j_darts_counter_temp3340]
                                                 [k_darts_counter_temp3340][2]
                                        + rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340
                                            + 1][k_darts_counter_temp3340][2]);
                            frct[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                [k_darts_counter_temp3340][3]
                                = frct[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                      [k_darts_counter_temp3340][3]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp3340)][j_darts_counter_temp3340
                                           + 1][k_darts_counter_temp3340][3]
                                        - flux[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                              [k_darts_counter_temp3340][3])
                                + dy4 * ty1
                                    * (rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340 - 1]
                                          [k_darts_counter_temp3340][3]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3340)]
                                                 [j_darts_counter_temp3340]
                                                 [k_darts_counter_temp3340][3]
                                        + rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340
                                            + 1][k_darts_counter_temp3340][3]);
                            frct[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                [k_darts_counter_temp3340][4]
                                = frct[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                      [k_darts_counter_temp3340][4]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp3340)][j_darts_counter_temp3340
                                           + 1][k_darts_counter_temp3340][4]
                                        - flux[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                              [k_darts_counter_temp3340][4])
                                + dy5 * ty1
                                    * (rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340 - 1]
                                          [k_darts_counter_temp3340][4]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3340)]
                                                 [j_darts_counter_temp3340]
                                                 [k_darts_counter_temp3340][4]
                                        + rsd[(i_darts_counter_temp3340)][j_darts_counter_temp3340
                                            + 1][k_darts_counter_temp3340][4]);
                        }
                        (*j) = j_darts_counter_temp3340;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp3340 = (*m);
                        for (; m_darts_counter_temp3340 < 5; m_darts_counter_temp3340++) {
                            frct[(i_darts_counter_temp3340)][1][k_darts_counter_temp3340]
                                [m_darts_counter_temp3340]
                                = frct[(i_darts_counter_temp3340)][1][k_darts_counter_temp3340]
                                      [m_darts_counter_temp3340]
                                - (*(*dsspm))
                                    * (+5.
                                            * rsd[(i_darts_counter_temp3340)][1]
                                                 [k_darts_counter_temp3340]
                                                 [m_darts_counter_temp3340]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3340)][2]
                                                 [k_darts_counter_temp3340]
                                                 [m_darts_counter_temp3340]
                                        + rsd[(i_darts_counter_temp3340)][3]
                                             [k_darts_counter_temp3340][m_darts_counter_temp3340]);
                            frct[(i_darts_counter_temp3340)][2][k_darts_counter_temp3340]
                                [m_darts_counter_temp3340]
                                = frct[(i_darts_counter_temp3340)][2][k_darts_counter_temp3340]
                                      [m_darts_counter_temp3340]
                                - (*(*dsspm))
                                    * (-4.
                                            * rsd[(i_darts_counter_temp3340)][1]
                                                 [k_darts_counter_temp3340]
                                                 [m_darts_counter_temp3340]
                                        + 6.
                                            * rsd[(i_darts_counter_temp3340)][2]
                                                 [k_darts_counter_temp3340]
                                                 [m_darts_counter_temp3340]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3340)][3]
                                                 [k_darts_counter_temp3340]
                                                 [m_darts_counter_temp3340]
                                        + rsd[(i_darts_counter_temp3340)][4]
                                             [k_darts_counter_temp3340][m_darts_counter_temp3340]);
                        }
                        (*m) = m_darts_counter_temp3340;
                    }
                    (*(*jst1)) = 3;
                    (*(*jend1)) = ny - 4;
                    {
                        /*Loop's init*/
                        (*j) = (*(*jst1));
                        int j_darts_counter_temp3340 = (*j);
                        for (; j_darts_counter_temp3340 <= (*(*jend1));
                             j_darts_counter_temp3340++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp3340 = (*m);
                                for (; m_darts_counter_temp3340 < 5; m_darts_counter_temp3340++) {
                                    frct[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                        [k_darts_counter_temp3340][m_darts_counter_temp3340]
                                        = frct[(i_darts_counter_temp3340)][j_darts_counter_temp3340]
                                              [k_darts_counter_temp3340][m_darts_counter_temp3340]
                                        - (*(*dsspm))
                                            * (rsd[(i_darts_counter_temp3340)]
                                                  [j_darts_counter_temp3340 - 2]
                                                  [k_darts_counter_temp3340]
                                                  [m_darts_counter_temp3340]
                                                - 4.
                                                    * rsd[(i_darts_counter_temp3340)]
                                                         [j_darts_counter_temp3340 - 1]
                                                         [k_darts_counter_temp3340]
                                                         [m_darts_counter_temp3340]
                                                + 6.
                                                    * rsd[(i_darts_counter_temp3340)]
                                                         [j_darts_counter_temp3340]
                                                         [k_darts_counter_temp3340]
                                                         [m_darts_counter_temp3340]
                                                - 4.
                                                    * rsd[(i_darts_counter_temp3340)]
                                                         [j_darts_counter_temp3340 + 1]
                                                         [k_darts_counter_temp3340]
                                                         [m_darts_counter_temp3340]
                                                + rsd[(i_darts_counter_temp3340)]
                                                     [j_darts_counter_temp3340 + 2]
                                                     [k_darts_counter_temp3340]
                                                     [m_darts_counter_temp3340]);
                                }
                                (*m) = m_darts_counter_temp3340;
                            }
                        }
                        (*j) = j_darts_counter_temp3340;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp3340 = (*m);
                        for (; m_darts_counter_temp3340 < 5; m_darts_counter_temp3340++) {
                            frct[(i_darts_counter_temp3340)][ny - 3][k_darts_counter_temp3340]
                                [m_darts_counter_temp3340]
                                = frct[(i_darts_counter_temp3340)][ny - 3][k_darts_counter_temp3340]
                                      [m_darts_counter_temp3340]
                                - (*(*dsspm))
                                    * (rsd[(i_darts_counter_temp3340)][ny - 5]
                                          [k_darts_counter_temp3340][m_darts_counter_temp3340]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3340)][ny - 4]
                                                 [k_darts_counter_temp3340]
                                                 [m_darts_counter_temp3340]
                                        + 6.
                                            * rsd[(i_darts_counter_temp3340)][ny - 3]
                                                 [k_darts_counter_temp3340]
                                                 [m_darts_counter_temp3340]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3340)][ny - 2]
                                                 [k_darts_counter_temp3340]
                                                 [m_darts_counter_temp3340]);
                            frct[(i_darts_counter_temp3340)][ny - 2][k_darts_counter_temp3340]
                                [m_darts_counter_temp3340]
                                = frct[(i_darts_counter_temp3340)][ny - 2][k_darts_counter_temp3340]
                                      [m_darts_counter_temp3340]
                                - (*(*dsspm))
                                    * (rsd[(i_darts_counter_temp3340)][ny - 4]
                                          [k_darts_counter_temp3340][m_darts_counter_temp3340]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3340)][ny - 3]
                                                 [k_darts_counter_temp3340]
                                                 [m_darts_counter_temp3340]
                                        + 5.
                                            * rsd[(i_darts_counter_temp3340)][ny - 2]
                                                 [k_darts_counter_temp3340]
                                                 [m_darts_counter_temp3340]);
                        }
                        (*m) = m_darts_counter_temp3340;
                    }
                }
                (*k) = k_darts_counter_temp3340;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets3340[0].decDep();
}
TP3340::TP3340(int in_numThreads, int in_mainCodeletID, TP2199* in_TPParent, int in_initIteration,
    int in_lastIteration, TP3340** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L2_darts3340(new int*[this->numThreads])
    , dsspm_darts3340(new double*[this->numThreads])
    , i_darts3340(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts3340(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jend1_darts3340(new int*[this->numThreads])
    , jst1_darts3340(new int*[this->numThreads])
    , k_darts3340(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts3340(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts3340(new double*[this->numThreads])
    , u21j_darts3340(new double*[this->numThreads])
    , u21jm1_darts3340(new double*[this->numThreads])
    , u31j_darts3340(new double*[this->numThreads])
    , u31jm1_darts3340(new double*[this->numThreads])
    , u41j_darts3340(new double*[this->numThreads])
    , u41jm1_darts3340(new double*[this->numThreads])
    , u51j_darts3340(new double*[this->numThreads])
    , u51jm1_darts3340(new double*[this->numThreads])
    , initIteration3340(in_initIteration)
    , lastIteration3340(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets3340(new _barrierCodelets3340[1])
    , checkInCodelets3341(new _checkInCodelets3341[this->numThreads])
{
    /*Initialize the loop parameters*/
    range3340 = abs(lastIteration3340 - initIteration3340) / 1;
    rangePerCodelet3340 = range3340 / numThreads;
    minIteration3340 = min<int>(lastIteration3340, initIteration3340);
    remainderRange3340 = range3340 % numThreads;
    /*Initialize inputs and vars.*/
    this->L2_darts3340 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->dsspm_darts3340
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jend1_darts3340 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jst1_darts3340 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts3340
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21j_darts3340
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21jm1_darts3340
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31j_darts3340
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31jm1_darts3340
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41j_darts3340
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41jm1_darts3340
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51j_darts3340
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51jm1_darts3340
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets3340[0] = _barrierCodelets3340(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets3341* checkInCodelets3341Ptr = (this->checkInCodelets3341);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets3341);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets3341Ptr) = _checkInCodelets3341(2, 1, this, codeletCounter);
#else
        (*checkInCodelets3341Ptr) = _checkInCodelets3341(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets3341Ptr).decDep();
        checkInCodelets3341Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP3340::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets3341[localID].setID(codeletID);
    this->checkInCodelets3341[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets3341[localID + this->baseNumThreads * i]
            = _checkInCodelets3341(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets3341[localID + this->baseNumThreads * i]
            = _checkInCodelets3341(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets3341[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets3341[localID + this->baseNumThreads * i].decDep();
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
TP3340::~TP3340()
{
    delete[] L2_darts3340;
    delete[] dsspm_darts3340;
    delete[] jend1_darts3340;
    delete[] jst1_darts3340;
    delete[] tmp_darts3340;
    delete[] u21j_darts3340;
    delete[] u21jm1_darts3340;
    delete[] u31j_darts3340;
    delete[] u31jm1_darts3340;
    delete[] u41j_darts3340;
    delete[] u41jm1_darts3340;
    delete[] u51j_darts3340;
    delete[] u51jm1_darts3340;
    delete[] barrierCodelets3340;
    delete[] checkInCodelets3341;
}
/*TP3975: OMPForDirective*/
void TP3975::_barrierCodelets3975::fire(void)
{
    TP3975* myTP = static_cast<TP3975*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets3975[0].decDep();
}
bool TP3975::requestNewRangeIterations3975(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet3975 * codeletID;
        int tempEndRange = rangePerCodelet3975 * (codeletID + 1);
        if (remainderRange3975 != 0) {
            if (codeletID < (uint32_t)remainderRange3975) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange3975;
                tempEndRange += remainderRange3975;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration3975;
        tempEndRange = tempEndRange * 1 + minIteration3975;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration3975 < lastIteration3975) {
            (this->inputsTPParent->i_darts3975[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts3975[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration3975;
        }
    }
    return isThereNewIteration;
}
void TP3975::_checkInCodelets3976::fire(void)
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
    this->inputsTPParent->dsspm_darts3975[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->dsspm_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->q_darts3975[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts3975[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21k_darts3975[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21k_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21km1_darts3975[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21km1_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31k_darts3975[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31k_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31km1_darts3975[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31km1_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41_darts3975[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41k_darts3975[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41k_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41km1_darts3975[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41km1_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51k_darts3975[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51k_darts2199[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51km1_darts3975[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51km1_darts2199[this->getID()]);

    /*printing node 3976: ForStmt*/
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
    double** dsspm = &(this->inputsTPParent->dsspm_darts3975[this->getLocalID()]);
    (void)dsspm /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts3975[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts3975[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts3975[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts3975[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts3975[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts3975[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21k = &(this->inputsTPParent->u21k_darts3975[this->getLocalID()]);
    (void)u21k /*OMP_SHARED_PRIVATE*/;
    double** u21km1 = &(this->inputsTPParent->u21km1_darts3975[this->getLocalID()]);
    (void)u21km1 /*OMP_SHARED_PRIVATE*/;
    double** u31k = &(this->inputsTPParent->u31k_darts3975[this->getLocalID()]);
    (void)u31k /*OMP_SHARED_PRIVATE*/;
    double** u31km1 = &(this->inputsTPParent->u31km1_darts3975[this->getLocalID()]);
    (void)u31km1 /*OMP_SHARED_PRIVATE*/;
    double** u41 = &(this->inputsTPParent->u41_darts3975[this->getLocalID()]);
    (void)u41 /*OMP_SHARED_PRIVATE*/;
    double** u41k = &(this->inputsTPParent->u41k_darts3975[this->getLocalID()]);
    (void)u41k /*OMP_SHARED_PRIVATE*/;
    double** u41km1 = &(this->inputsTPParent->u41km1_darts3975[this->getLocalID()]);
    (void)u41km1 /*OMP_SHARED_PRIVATE*/;
    double** u51k = &(this->inputsTPParent->u51k_darts3975[this->getLocalID()]);
    (void)u51k /*OMP_SHARED_PRIVATE*/;
    double** u51km1 = &(this->inputsTPParent->u51km1_darts3975[this->getLocalID()]);
    (void)u51km1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations3975(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets3975[0].decDep();
        return;
    }
    for (int i_darts_counter_temp3975 = (*i); i_darts_counter_temp3975 <= endRange
         && i_darts_counter_temp3975 <= this->inputsTPParent->lastIteration3975;
         i_darts_counter_temp3975++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp3975 = (*j);
                for (; j_darts_counter_temp3975 <= jend; j_darts_counter_temp3975++) {
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp3975 = (*k);
                        for (; k_darts_counter_temp3975 <= nz - 1; k_darts_counter_temp3975++) {
                            flux[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                [k_darts_counter_temp3975][0]
                                = rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                     [k_darts_counter_temp3975][3];
                            (*(*u41)) = rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                           [k_darts_counter_temp3975][3]
                                / rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                     [k_darts_counter_temp3975][0];
                            (*(*q)) = 0.5
                                * (rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                      [k_darts_counter_temp3975][1]
                                        * rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                             [k_darts_counter_temp3975][1]
                                    + rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                         [k_darts_counter_temp3975][2]
                                        * rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                             [k_darts_counter_temp3975][2]
                                    + rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                         [k_darts_counter_temp3975][3]
                                        * rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                             [k_darts_counter_temp3975][3])
                                / rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                     [k_darts_counter_temp3975][0];
                            flux[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                [k_darts_counter_temp3975][1]
                                = rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                     [k_darts_counter_temp3975][1]
                                * (*(*u41));
                            flux[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                [k_darts_counter_temp3975][2]
                                = rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                     [k_darts_counter_temp3975][2]
                                * (*(*u41));
                            flux[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                [k_darts_counter_temp3975][3]
                                = rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                     [k_darts_counter_temp3975][3]
                                    * (*(*u41))
                                + 0.40000000000000002
                                    * (rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                          [k_darts_counter_temp3975][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                [k_darts_counter_temp3975][4]
                                = (1.3999999999999999
                                          * rsd[(i_darts_counter_temp3975)]
                                               [j_darts_counter_temp3975][k_darts_counter_temp3975]
                                               [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u41));
                        }
                        (*k) = k_darts_counter_temp3975;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp3975 = (*k);
                        for (; k_darts_counter_temp3975 <= nz - 2; k_darts_counter_temp3975++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp3975 = (*m);
                                for (; m_darts_counter_temp3975 < 5; m_darts_counter_temp3975++) {
                                    frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                        [k_darts_counter_temp3975][m_darts_counter_temp3975]
                                        = frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                              [k_darts_counter_temp3975][m_darts_counter_temp3975]
                                        - tz2
                                            * (flux[(i_darts_counter_temp3975)]
                                                   [j_darts_counter_temp3975]
                                                   [k_darts_counter_temp3975 + 1]
                                                   [m_darts_counter_temp3975]
                                                - flux[k_darts_counter_temp3975 - 1]
                                                      [(i_darts_counter_temp3975)]
                                                      [j_darts_counter_temp3975]
                                                      [m_darts_counter_temp3975]);
                                }
                                (*m) = m_darts_counter_temp3975;
                            }
                        }
                        (*k) = k_darts_counter_temp3975;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp3975 = (*k);
                        for (; k_darts_counter_temp3975 <= nz - 1; k_darts_counter_temp3975++) {
                            (*(*tmp)) = 1.
                                / rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                     [k_darts_counter_temp3975][0];
                            (*(*u21k)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                     [k_darts_counter_temp3975][1];
                            (*(*u31k)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                     [k_darts_counter_temp3975][2];
                            (*(*u41k)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                     [k_darts_counter_temp3975][3];
                            (*(*u51k)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                     [k_darts_counter_temp3975][4];
                            (*(*tmp)) = 1.
                                / rsd[k_darts_counter_temp3975 - 1][(i_darts_counter_temp3975)]
                                     [j_darts_counter_temp3975][0];
                            (*(*u21km1)) = (*(*tmp))
                                * rsd[k_darts_counter_temp3975 - 1][(i_darts_counter_temp3975)]
                                     [j_darts_counter_temp3975][1];
                            (*(*u31km1)) = (*(*tmp))
                                * rsd[k_darts_counter_temp3975 - 1][(i_darts_counter_temp3975)]
                                     [j_darts_counter_temp3975][2];
                            (*(*u41km1)) = (*(*tmp))
                                * rsd[k_darts_counter_temp3975 - 1][(i_darts_counter_temp3975)]
                                     [j_darts_counter_temp3975][3];
                            (*(*u51km1)) = (*(*tmp))
                                * rsd[k_darts_counter_temp3975 - 1][(i_darts_counter_temp3975)]
                                     [j_darts_counter_temp3975][4];
                            flux[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                [k_darts_counter_temp3975][1]
                                = tz3 * ((*(*u21k)) - (*(*u21km1)));
                            flux[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                [k_darts_counter_temp3975][2]
                                = tz3 * ((*(*u31k)) - (*(*u31km1)));
                            flux[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                [k_darts_counter_temp3975][3]
                                = (4. / 3.) * tz3 * ((*(*u41k)) - (*(*u41km1)));
                            flux[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                [k_darts_counter_temp3975][4]
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
                        (*k) = k_darts_counter_temp3975;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp3975 = (*k);
                        for (; k_darts_counter_temp3975 <= nz - 2; k_darts_counter_temp3975++) {
                            frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                [k_darts_counter_temp3975][0]
                                = frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                      [k_darts_counter_temp3975][0]
                                + dz1 * tz1
                                    * (rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                          [k_darts_counter_temp3975 + 1][0]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3975)]
                                                 [j_darts_counter_temp3975]
                                                 [k_darts_counter_temp3975][0]
                                        + rsd[k_darts_counter_temp3975 - 1]
                                             [(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                             [0]);
                            frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                [k_darts_counter_temp3975][1]
                                = frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                      [k_darts_counter_temp3975][1]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                           [k_darts_counter_temp3975 + 1][1]
                                        - flux[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                              [k_darts_counter_temp3975][1])
                                + dz2 * tz1
                                    * (rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                          [k_darts_counter_temp3975 + 1][1]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3975)]
                                                 [j_darts_counter_temp3975]
                                                 [k_darts_counter_temp3975][1]
                                        + rsd[k_darts_counter_temp3975 - 1]
                                             [(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                             [1]);
                            frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                [k_darts_counter_temp3975][2]
                                = frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                      [k_darts_counter_temp3975][2]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                           [k_darts_counter_temp3975 + 1][2]
                                        - flux[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                              [k_darts_counter_temp3975][2])
                                + dz3 * tz1
                                    * (rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                          [k_darts_counter_temp3975 + 1][2]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3975)]
                                                 [j_darts_counter_temp3975]
                                                 [k_darts_counter_temp3975][2]
                                        + rsd[k_darts_counter_temp3975 - 1]
                                             [(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                             [2]);
                            frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                [k_darts_counter_temp3975][3]
                                = frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                      [k_darts_counter_temp3975][3]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                           [k_darts_counter_temp3975 + 1][3]
                                        - flux[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                              [k_darts_counter_temp3975][3])
                                + dz4 * tz1
                                    * (rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                          [k_darts_counter_temp3975 + 1][3]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3975)]
                                                 [j_darts_counter_temp3975]
                                                 [k_darts_counter_temp3975][3]
                                        + rsd[k_darts_counter_temp3975 - 1]
                                             [(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                             [3]);
                            frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                [k_darts_counter_temp3975][4]
                                = frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                      [k_darts_counter_temp3975][4]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                           [k_darts_counter_temp3975 + 1][4]
                                        - flux[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                              [k_darts_counter_temp3975][4])
                                + dz5 * tz1
                                    * (rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                          [k_darts_counter_temp3975 + 1][4]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3975)]
                                                 [j_darts_counter_temp3975]
                                                 [k_darts_counter_temp3975][4]
                                        + rsd[k_darts_counter_temp3975 - 1]
                                             [(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                             [4]);
                        }
                        (*k) = k_darts_counter_temp3975;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp3975 = (*m);
                        for (; m_darts_counter_temp3975 < 5; m_darts_counter_temp3975++) {
                            frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975][1]
                                [m_darts_counter_temp3975]
                                = frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975][1]
                                      [m_darts_counter_temp3975]
                                - (*(*dsspm))
                                    * (+5.
                                            * rsd[(i_darts_counter_temp3975)]
                                                 [j_darts_counter_temp3975][1]
                                                 [m_darts_counter_temp3975]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3975)]
                                                 [j_darts_counter_temp3975][2]
                                                 [m_darts_counter_temp3975]
                                        + rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                             [3][m_darts_counter_temp3975]);
                            frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975][2]
                                [m_darts_counter_temp3975]
                                = frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975][2]
                                      [m_darts_counter_temp3975]
                                - (*(*dsspm))
                                    * (-4.
                                            * rsd[(i_darts_counter_temp3975)]
                                                 [j_darts_counter_temp3975][1]
                                                 [m_darts_counter_temp3975]
                                        + 6.
                                            * rsd[(i_darts_counter_temp3975)]
                                                 [j_darts_counter_temp3975][2]
                                                 [m_darts_counter_temp3975]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3975)]
                                                 [j_darts_counter_temp3975][3]
                                                 [m_darts_counter_temp3975]
                                        + rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                             [4][m_darts_counter_temp3975]);
                        }
                        (*m) = m_darts_counter_temp3975;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 3;
                        int k_darts_counter_temp3975 = (*k);
                        for (; k_darts_counter_temp3975 <= nz - 4; k_darts_counter_temp3975++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp3975 = (*m);
                                for (; m_darts_counter_temp3975 < 5; m_darts_counter_temp3975++) {
                                    frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                        [k_darts_counter_temp3975][m_darts_counter_temp3975]
                                        = frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                              [k_darts_counter_temp3975][m_darts_counter_temp3975]
                                        - (*(*dsspm))
                                            * (rsd[(i_darts_counter_temp3975)]
                                                  [j_darts_counter_temp3975]
                                                  [k_darts_counter_temp3975 - 2]
                                                  [m_darts_counter_temp3975]
                                                - 4.
                                                    * rsd[k_darts_counter_temp3975 - 1]
                                                         [(i_darts_counter_temp3975)]
                                                         [j_darts_counter_temp3975]
                                                         [m_darts_counter_temp3975]
                                                + 6.
                                                    * rsd[(i_darts_counter_temp3975)]
                                                         [j_darts_counter_temp3975]
                                                         [k_darts_counter_temp3975]
                                                         [m_darts_counter_temp3975]
                                                - 4.
                                                    * rsd[(i_darts_counter_temp3975)]
                                                         [j_darts_counter_temp3975]
                                                         [k_darts_counter_temp3975 + 1]
                                                         [m_darts_counter_temp3975]
                                                + rsd[(i_darts_counter_temp3975)]
                                                     [j_darts_counter_temp3975]
                                                     [k_darts_counter_temp3975 + 2]
                                                     [m_darts_counter_temp3975]);
                                }
                                (*m) = m_darts_counter_temp3975;
                            }
                        }
                        (*k) = k_darts_counter_temp3975;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp3975 = (*m);
                        for (; m_darts_counter_temp3975 < 5; m_darts_counter_temp3975++) {
                            frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975][nz - 3]
                                [m_darts_counter_temp3975]
                                = frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975][nz - 3]
                                      [m_darts_counter_temp3975]
                                - (*(*dsspm))
                                    * (rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                          [nz - 5][m_darts_counter_temp3975]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3975)]
                                                 [j_darts_counter_temp3975][nz - 4]
                                                 [m_darts_counter_temp3975]
                                        + 6.
                                            * rsd[(i_darts_counter_temp3975)]
                                                 [j_darts_counter_temp3975][nz - 3]
                                                 [m_darts_counter_temp3975]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3975)]
                                                 [j_darts_counter_temp3975][nz - 2]
                                                 [m_darts_counter_temp3975]);
                            frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975][nz - 2]
                                [m_darts_counter_temp3975]
                                = frct[(i_darts_counter_temp3975)][j_darts_counter_temp3975][nz - 2]
                                      [m_darts_counter_temp3975]
                                - (*(*dsspm))
                                    * (rsd[(i_darts_counter_temp3975)][j_darts_counter_temp3975]
                                          [nz - 4][m_darts_counter_temp3975]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3975)]
                                                 [j_darts_counter_temp3975][nz - 3]
                                                 [m_darts_counter_temp3975]
                                        + 5.
                                            * rsd[(i_darts_counter_temp3975)]
                                                 [j_darts_counter_temp3975][nz - 2]
                                                 [m_darts_counter_temp3975]);
                        }
                        (*m) = m_darts_counter_temp3975;
                    }
                }
                (*j) = j_darts_counter_temp3975;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets3975[0].decDep();
}
TP3975::TP3975(int in_numThreads, int in_mainCodeletID, TP2199* in_TPParent, int in_initIteration,
    int in_lastIteration, TP3975** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , dsspm_darts3975(new double*[this->numThreads])
    , i_darts3975(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts3975(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts3975(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts3975(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts3975(new double*[this->numThreads])
    , tmp_darts3975(new double*[this->numThreads])
    , u21k_darts3975(new double*[this->numThreads])
    , u21km1_darts3975(new double*[this->numThreads])
    , u31k_darts3975(new double*[this->numThreads])
    , u31km1_darts3975(new double*[this->numThreads])
    , u41_darts3975(new double*[this->numThreads])
    , u41k_darts3975(new double*[this->numThreads])
    , u41km1_darts3975(new double*[this->numThreads])
    , u51k_darts3975(new double*[this->numThreads])
    , u51km1_darts3975(new double*[this->numThreads])
    , initIteration3975(in_initIteration)
    , lastIteration3975(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets3975(new _barrierCodelets3975[1])
    , checkInCodelets3976(new _checkInCodelets3976[this->numThreads])
{
    /*Initialize the loop parameters*/
    range3975 = abs(lastIteration3975 - initIteration3975) / 1;
    rangePerCodelet3975 = range3975 / numThreads;
    minIteration3975 = min<int>(lastIteration3975, initIteration3975);
    remainderRange3975 = range3975 % numThreads;
    /*Initialize inputs and vars.*/
    this->dsspm_darts3975
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->q_darts3975 = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts3975
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21k_darts3975
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21km1_darts3975
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31k_darts3975
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31km1_darts3975
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41_darts3975
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41k_darts3975
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41km1_darts3975
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51k_darts3975
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51km1_darts3975
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets3975[0] = _barrierCodelets3975(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets3976* checkInCodelets3976Ptr = (this->checkInCodelets3976);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets3976);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets3976Ptr) = _checkInCodelets3976(2, 1, this, codeletCounter);
#else
        (*checkInCodelets3976Ptr) = _checkInCodelets3976(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets3976Ptr).decDep();
        checkInCodelets3976Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP3975::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets3976[localID].setID(codeletID);
    this->checkInCodelets3976[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets3976[localID + this->baseNumThreads * i]
            = _checkInCodelets3976(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets3976[localID + this->baseNumThreads * i]
            = _checkInCodelets3976(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets3976[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets3976[localID + this->baseNumThreads * i].decDep();
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
TP3975::~TP3975()
{
    delete[] dsspm_darts3975;
    delete[] q_darts3975;
    delete[] tmp_darts3975;
    delete[] u21k_darts3975;
    delete[] u21km1_darts3975;
    delete[] u31k_darts3975;
    delete[] u31km1_darts3975;
    delete[] u41_darts3975;
    delete[] u41k_darts3975;
    delete[] u41km1_darts3975;
    delete[] u51k_darts3975;
    delete[] u51km1_darts3975;
    delete[] barrierCodelets3975;
    delete[] checkInCodelets3976;
}
/*TP7: TP_jacld*/
void TP7::_checkInCodelets4880::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 4880: DeclStmt*/

    /*printing node 4881: DeclStmt*/

    /*printing node 4882: DeclStmt*/

    /*printing node 4883: DeclStmt*/

    /*printing node 4884: DeclStmt*/

    /*printing node 4885: BinaryOperator*/
    (this->inputsTPParent->r43_darts7[this->getID()]) = (4. / 3.);

    /*printing node 4889: BinaryOperator*/
    (this->inputsTPParent->c1345_darts7[this->getID()])
        = 1.3999999999999999 * 0.10000000000000001 * 1. * 1.3999999999999999;

    /*printing node 4897: BinaryOperator*/
    (this->inputsTPParent->c34_darts7[this->getID()]) = 0.10000000000000001 * 1.;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 4880 nextRegion: 4901 */
    myTP->controlTPParent->checkInCodelets4901[this->getID()].decDep();
}
void TP7::_checkInCodelets4901::fire(void)
{
    /*region 4901 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP4901;
    if (idx < myTP->TPsToUse4901) {
        if (!__sync_val_compare_and_swap(&(myTP->TP4901_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse4901;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse4901;
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
            if (idx == myTP->TPsToUse4901 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse4901 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP4901>(myTP, myTP->codeletsPerTP4901 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP4901Ptr[idx]));
#else
            place<TP4901>(idx, myTP, myTP->codeletsPerTP4901 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP4901Ptr[idx]));
#endif
        } else {
            if (myTP->TP4901Ptr[idx] != nullptr) {
                myTP->TP4901Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP7::_barrierCodelets4901::fire(void)
{
    TP7* myTP = static_cast<TP7*>(myTP_);

    for (size_t codeletsCounter = 0; codeletsCounter < (size_t)myTP->numThreads;
         codeletsCounter++) {
        myTP->nextCodeletsjacld[codeletsCounter]->decDep();
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
    , TP4901Ptr(new TP4901*[NUMTPS4901])
    , TP4901_alreadyLaunched(new size_t[NUMTPS4901])
    , numTPsSet4901(0)
    , numTPsReady4901(0)
    , TPsToUse4901(NUMTPS4901)
    , codeletsPerTP4901(this->numThreads / NUMTPS4901)
    , totalCodelets4901(this->TPsToUse4901 * this->codeletsPerTP4901)
    , checkInCodelets4880(new _checkInCodelets4880[this->numThreads])
    , checkInCodelets4901(new _checkInCodelets4901[this->numThreads])
    , barrierCodelets4901(new _barrierCodelets4901[1])
{
    barrierCodelets4901[0] = _barrierCodelets4901(NUMTPS4901, NUMTPS4901, this, 0);
    _checkInCodelets4901* checkInCodelets4901Ptr = (this->checkInCodelets4901);
    for (int i = 0; i < NUMTPS4901; i++) {
        TP4901Ptr[i] = nullptr;
        TP4901_alreadyLaunched[i] = 0;
    }
    _checkInCodelets4880* checkInCodelets4880Ptr = (this->checkInCodelets4880);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets4880);
#endif
    for (size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++) {
        (*checkInCodelets4901Ptr) = _checkInCodelets4901(1, 1, this, codeletCounter);
        checkInCodelets4901Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets4880Ptr) = _checkInCodelets4880(2, 1, this, codeletCounter);
#else
        (*checkInCodelets4880Ptr) = _checkInCodelets4880(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets4880Ptr).decDep();
        checkInCodelets4880Ptr++;
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
    delete[] barrierCodelets4901;
    delete[] checkInCodelets4901;
    delete[] checkInCodelets4880;
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
/*TP4901: OMPForDirective*/
void TP4901::_barrierCodelets4901::fire(void)
{
    TP4901* myTP = static_cast<TP4901*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets4901[0].decDep();
}
bool TP4901::requestNewRangeIterations4901(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 1*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet4901 * codeletID;
        int tempEndRange = rangePerCodelet4901 * (codeletID + 1);
        if (remainderRange4901 != 0) {
            if (codeletID < (uint32_t)remainderRange4901) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange4901;
                tempEndRange += remainderRange4901;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration4901;
        tempEndRange = tempEndRange * 1 + minIteration4901;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration4901 < lastIteration4901) {
            (this->inputsTPParent->i_darts4901[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts4901[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration4901;
        }
    }
    return isThereNewIteration;
}
void TP4901::_checkInCodelets4902::fire(void)
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
    this->inputsTPParent->c1345_darts4901[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->c1345_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->c34_darts4901[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->c34_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->k_darts4901[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->k_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->r43_darts4901[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->r43_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp1_darts4901[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp1_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp2_darts4901[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp2_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp3_darts4901[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp3_darts7[this->getID()]);

    /*printing node 4902: ForStmt*/
    /*var: c1345*/
    /*var: c34*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: r43*/
    /*var: tmp1*/
    /*var: tmp2*/
    /*var: tmp3*/
    double** c1345 = &(this->inputsTPParent->c1345_darts4901[this->getLocalID()]);
    (void)c1345 /*OMP_SHARED_PRIVATE*/;
    double** c34 = &(this->inputsTPParent->c34_darts4901[this->getLocalID()]);
    (void)c34 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts4901[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts4901[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** k = &(this->inputsTPParent->k_darts4901[this->getLocalID()]);
    (void)k /*OMP_SHARED_PRIVATE*/;
    double** r43 = &(this->inputsTPParent->r43_darts4901[this->getLocalID()]);
    (void)r43 /*OMP_SHARED_PRIVATE*/;
    double** tmp1 = &(this->inputsTPParent->tmp1_darts4901[this->getLocalID()]);
    (void)tmp1 /*OMP_SHARED_PRIVATE*/;
    double** tmp2 = &(this->inputsTPParent->tmp2_darts4901[this->getLocalID()]);
    (void)tmp2 /*OMP_SHARED_PRIVATE*/;
    double** tmp3 = &(this->inputsTPParent->tmp3_darts4901[this->getLocalID()]);
    (void)tmp3 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations4901(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets4901[0].decDep();
        return;
    }
    for (int i_darts_counter_temp4901 = (*i); i_darts_counter_temp4901 <= endRange
         && i_darts_counter_temp4901 <= this->inputsTPParent->lastIteration4901;
         i_darts_counter_temp4901++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp4901 = (*j);
                for (; j_darts_counter_temp4901 <= jend; j_darts_counter_temp4901++) {
                    (*(*tmp1))
                        = 1. / u[(i_darts_counter_temp4901)][j_darts_counter_temp4901][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][0]
                        = 1. + dt * 2. * (tx1 * dx1 + ty1 * dy1 + tz1 * dz1);
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][1] = 0.;
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][2] = 0.;
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][3] = 0.;
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][4] = 0.;
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][0] = dt * 2.
                        * (tx1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                       [(*(*k))][1])
                            + ty1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                       [(*(*k))][1])
                            + tz1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                       [(*(*k))][1]));
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][1] = 1.
                        + dt * 2.
                            * (tx1 * (*(*r43)) * (*(*c34)) * (*(*tmp1))
                                + ty1 * (*(*c34)) * (*(*tmp1)) + tz1 * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx2 + ty1 * dy2 + tz1 * dz2);
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][2] = 0.;
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][3] = 0.;
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][4] = 0.;
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][0] = dt * 2.
                        * (tx1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                       [(*(*k))][2])
                            + ty1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                       [(*(*k))][2])
                            + tz1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                       [(*(*k))][2]));
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][1] = 0.;
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][2] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c34)) * (*(*tmp1))
                                + ty1 * (*(*r43)) * (*(*c34)) * (*(*tmp1))
                                + tz1 * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx3 + ty1 * dy3 + tz1 * dz3);
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][3] = 0.;
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][4] = 0.;
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][0] = dt * 2.
                        * (tx1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                       [(*(*k))][3])
                            + ty1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                       [(*(*k))][3])
                            + tz1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                       [(*(*k))][3]));
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][1] = 0.;
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][2] = 0.;
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][3] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c34)) * (*(*tmp1)) + ty1 * (*(*c34)) * (*(*tmp1))
                                + tz1 * (*(*r43)) * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx4 + ty1 * dy4 + tz1 * dz4);
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][4] = 0.;
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][0] = dt * 2.
                        * (tx1
                                * (-((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp4901)]
                                                [j_darts_counter_temp4901][(*(*k))][1])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp4901)]
                                                [j_darts_counter_temp4901][(*(*k))][2])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp4901)]
                                                [j_darts_counter_temp4901][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                           [(*(*k))][4])
                            + ty1
                                * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp4901)]
                                                [j_darts_counter_temp4901][(*(*k))][1])))
                                    - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp4901)]
                                                [j_darts_counter_temp4901][(*(*k))][2])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp4901)]
                                                [j_darts_counter_temp4901][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                           [(*(*k))][4])
                            + tz1
                                * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp4901)]
                                                [j_darts_counter_temp4901][(*(*k))][1])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp4901)]
                                                [j_darts_counter_temp4901][(*(*k))][2])))
                                    - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp4901)]
                                                [j_darts_counter_temp4901][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                           [(*(*k))][4]));
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][1] = dt * 2.
                        * (tx1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901][(*(*k))]
                                   [1]
                            + ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901][(*(*k))]
                                   [1]
                            + tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901][(*(*k))]
                                   [1]);
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][2] = dt * 2.
                        * (tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901][(*(*k))]
                                   [2]
                            + ty1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901][(*(*k))]
                                   [2]
                            + tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901][(*(*k))]
                                   [2]);
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][3] = dt * 2.
                        * (tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901][(*(*k))]
                                   [3]
                            + ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901][(*(*k))]
                                   [3]
                            + tz1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901][(*(*k))]
                                   [3]);
                    d[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][4] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c1345)) * (*(*tmp1)) + ty1 * (*(*c1345)) * (*(*tmp1))
                                + tz1 * (*(*c1345)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx5 + ty1 * dy5 + tz1 * dz5);
                    (*(*tmp1)) = 1.
                        / u[(*(*k)) - 1][(i_darts_counter_temp4901)][j_darts_counter_temp4901][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][0] = -dt * tz1 * dz1;
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][1] = 0.;
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][2] = 0.;
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][3] = -dt * tz2;
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][4] = 0.;
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][0] = -dt * tz2
                            * (-(u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                  [j_darts_counter_temp4901][1]
                                   * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                      [j_darts_counter_temp4901][3])
                                * (*(*tmp2)))
                        - dt * tz1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                   [j_darts_counter_temp4901][1]);
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][1] = -dt * tz2
                            * (u[(*(*k)) - 1][(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * (*(*c34)) * (*(*tmp1)) - dt * tz1 * dz2;
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][2] = 0.;
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][3] = -dt * tz2
                        * (u[(*(*k)) - 1][(i_darts_counter_temp4901)][j_darts_counter_temp4901][1]
                            * (*(*tmp1)));
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][4] = 0.;
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][0] = -dt * tz2
                            * (-(u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                  [j_darts_counter_temp4901][2]
                                   * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                      [j_darts_counter_temp4901][3])
                                * (*(*tmp2)))
                        - dt * tz1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                   [j_darts_counter_temp4901][2]);
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][1] = 0.;
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][2] = -dt * tz2
                            * (u[(*(*k)) - 1][(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * ((*(*c34)) * (*(*tmp1))) - dt * tz1 * dz3;
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][3] = -dt * tz2
                        * (u[(*(*k)) - 1][(i_darts_counter_temp4901)][j_darts_counter_temp4901][2]
                            * (*(*tmp1)));
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][4] = 0.;
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][0] = -dt * tz2
                            * (-(u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                  [j_darts_counter_temp4901][3]
                                   * (*(*tmp1)))
                                    * (u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                        [j_darts_counter_temp4901][3]
                                        * (*(*tmp1)))
                                + 0.5 * 0.40000000000000002
                                    * ((u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                         [j_darts_counter_temp4901][1]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901][1]
                                           + u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                              [j_darts_counter_temp4901][2]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901][2]
                                           + u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                              [j_darts_counter_temp4901][3]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901][3])
                                        * (*(*tmp2))))
                        - dt * tz1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                   [j_darts_counter_temp4901][3]);
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][1] = -dt * tz2
                        * (-0.40000000000000002
                            * (u[(*(*k)) - 1][(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                [1]
                                * (*(*tmp1))));
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][2] = -dt * tz2
                        * (-0.40000000000000002
                            * (u[(*(*k)) - 1][(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                [2]
                                * (*(*tmp1))));
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][3] = -dt * tz2
                            * (2. - 0.40000000000000002)
                            * (u[(*(*k)) - 1][(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * tz1 * dz4;
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][4]
                        = -dt * tz2 * 0.40000000000000002;
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][0] = -dt * tz2
                            * ((0.40000000000000002
                                       * (u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                           [j_darts_counter_temp4901][1]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901][1]
                                           + u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                              [j_darts_counter_temp4901][2]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901][2]
                                           + u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                              [j_darts_counter_temp4901][3]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                           [j_darts_counter_temp4901][4]
                                           * (*(*tmp1))))
                                * (u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                    [j_darts_counter_temp4901][3]
                                    * (*(*tmp1))))
                        - dt * tz1
                            * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                        [j_darts_counter_temp4901][1]
                                        * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                           [j_darts_counter_temp4901][1])
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                        [j_darts_counter_temp4901][2]
                                        * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                           [j_darts_counter_temp4901][2])
                                - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                        [j_darts_counter_temp4901][3]
                                        * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                           [j_darts_counter_temp4901][3])
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                       [j_darts_counter_temp4901][4]);
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][1] = -dt * tz2
                            * (-0.40000000000000002
                                * (u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                    [j_darts_counter_temp4901][1]
                                    * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                       [j_darts_counter_temp4901][3])
                                * (*(*tmp2)))
                        - dt * tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(*(*k)) - 1][(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                               [1];
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][2] = -dt * tz2
                            * (-0.40000000000000002
                                * (u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                    [j_darts_counter_temp4901][2]
                                    * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                       [j_darts_counter_temp4901][3])
                                * (*(*tmp2)))
                        - dt * tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(*(*k)) - 1][(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                               [2];
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][3] = -dt * tz2
                            * (1.3999999999999999
                                    * (u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                        [j_darts_counter_temp4901][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                         [j_darts_counter_temp4901][1]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901][1]
                                           + u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                              [j_darts_counter_temp4901][2]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901][2]
                                           + 3.
                                               * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901][3]
                                               * u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901][3])
                                        * (*(*tmp2))))
                        - dt * tz1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(*(*k)) - 1][(i_darts_counter_temp4901)][j_darts_counter_temp4901]
                               [3];
                    a[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][4] = -dt * tz2
                            * (1.3999999999999999
                                * (u[(*(*k)) - 1][(i_darts_counter_temp4901)]
                                    [j_darts_counter_temp4901][3]
                                    * (*(*tmp1))))
                        - dt * tz1 * (*(*c1345)) * (*(*tmp1)) - dt * tz1 * dz5;
                    (*(*tmp1)) = 1.
                        / u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][0] = -dt * ty1 * dy1;
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][1] = 0.;
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][2] = -dt * ty2;
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][3] = 0.;
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][4] = 0.;
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][0] = -dt * ty2
                            * (-(u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                  [(*(*k))][1]
                                   * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                      [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * ty1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                   [(*(*k))][1]);
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][1] = -dt * ty2
                            * (u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1][(*(*k))]
                                [2]
                                * (*(*tmp1)))
                        - dt * ty1 * ((*(*c34)) * (*(*tmp1))) - dt * ty1 * dy2;
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][2] = -dt * ty2
                        * (u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1][(*(*k))][1]
                            * (*(*tmp1)));
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][3] = 0.;
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][4] = 0.;
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][0] = -dt * ty2
                            * (-(u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                  [(*(*k))][2]
                                   * (*(*tmp1)))
                                    * (u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                        [(*(*k))][2]
                                        * (*(*tmp1)))
                                + 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                         [(*(*k))][1]
                                               * u[(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901 - 1][(*(*k))][1]
                                           + u[(i_darts_counter_temp4901)]
                                              [j_darts_counter_temp4901 - 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901 - 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp4901)]
                                              [j_darts_counter_temp4901 - 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901 - 1][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * ty1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                   [(*(*k))][2]);
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][1] = -dt * ty2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1][(*(*k))]
                                [1]
                                * (*(*tmp1))));
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][2] = -dt * ty2
                            * ((2. - 0.40000000000000002)
                                * (u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * ty1 * dy3;
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][3] = -dt * ty2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1][(*(*k))]
                                [3]
                                * (*(*tmp1))));
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][4]
                        = -dt * ty2 * 0.40000000000000002;
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][0] = -dt * ty2
                            * (-(u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                  [(*(*k))][2]
                                   * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                      [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * ty1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                   [(*(*k))][3]);
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][1] = 0.;
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][2] = -dt * ty2
                        * (u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1][(*(*k))][3]
                            * (*(*tmp1)));
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][3] = -dt * ty2
                            * (u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1][(*(*k))]
                                [2]
                                * (*(*tmp1)))
                        - dt * ty1 * ((*(*c34)) * (*(*tmp1))) - dt * ty1 * dy4;
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][4] = 0.;
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][0] = -dt * ty2
                            * ((0.40000000000000002
                                       * (u[(i_darts_counter_temp4901)]
                                           [j_darts_counter_temp4901 - 1][(*(*k))][1]
                                               * u[(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901 - 1][(*(*k))][1]
                                           + u[(i_darts_counter_temp4901)]
                                              [j_darts_counter_temp4901 - 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901 - 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp4901)]
                                              [j_darts_counter_temp4901 - 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901 - 1][(*(*k))][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(i_darts_counter_temp4901)]
                                           [j_darts_counter_temp4901 - 1][(*(*k))][4]
                                           * (*(*tmp1))))
                                * (u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1
                            * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                          [(*(*k))][1])
                                        * (u[(i_darts_counter_temp4901)]
                                            [j_darts_counter_temp4901 - 1][(*(*k))][1])))
                                - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                          [(*(*k))][2])
                                        * (u[(i_darts_counter_temp4901)]
                                            [j_darts_counter_temp4901 - 1][(*(*k))][2])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                          [(*(*k))][3])
                                        * (u[(i_darts_counter_temp4901)]
                                            [j_darts_counter_temp4901 - 1][(*(*k))][3])))
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                       [(*(*k))][4]);
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][1] = -dt * ty2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                    [(*(*k))][1]
                                    * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                       [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1][(*(*k))]
                               [1];
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][2] = -dt * ty2
                            * (1.3999999999999999
                                    * (u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                        [(*(*k))][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                         [(*(*k))][1]
                                               * u[(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901 - 1][(*(*k))][1]
                                           + 3.
                                               * u[(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901 - 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901 - 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp4901)]
                                              [j_darts_counter_temp4901 - 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp4901)]
                                                  [j_darts_counter_temp4901 - 1][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * ty1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1][(*(*k))]
                               [2];
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][3] = -dt * ty2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                    [(*(*k))][2]
                                    * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                       [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1][(*(*k))]
                               [3];
                    b[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][4] = -dt * ty2
                            * (1.3999999999999999
                                * (u[(i_darts_counter_temp4901)][j_darts_counter_temp4901 - 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1 * (*(*c1345)) * (*(*tmp1)) - dt * ty1 * dy5;
                    (*(*tmp1)) = 1.
                        / u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][0] = -dt * tx1 * dx1;
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][1] = -dt * tx2;
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][2] = 0.;
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][3] = 0.;
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][0][4] = 0.;
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][0] = -dt * tx2
                            * (-(u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901][(*(*k))]
                                  [1]
                                   * (*(*tmp1)))
                                    * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                        [(*(*k))][1]
                                        * (*(*tmp1)))
                                + 0.40000000000000002 * 0.5
                                    * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                        [(*(*k))][1]
                                            * u[(i_darts_counter_temp4901)-1]
                                               [j_darts_counter_temp4901][(*(*k))][1]
                                        + u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                           [(*(*k))][2]
                                            * u[(i_darts_counter_temp4901)-1]
                                               [j_darts_counter_temp4901][(*(*k))][2]
                                        + u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                           [(*(*k))][3]
                                            * u[(i_darts_counter_temp4901)-1]
                                               [j_darts_counter_temp4901][(*(*k))][3])
                                    * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901][(*(*k))]
                                   [1]);
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][1] = -dt * tx2
                            * ((2. - 0.40000000000000002)
                                * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * tx1 * dx2;
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][2] = -dt * tx2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901][(*(*k))][2]
                                * (*(*tmp1))));
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][3] = -dt * tx2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901][(*(*k))][3]
                                * (*(*tmp1))));
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][1][4]
                        = -dt * tx2 * 0.40000000000000002;
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][0] = -dt * tx2
                            * (-(u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901][(*(*k))]
                                  [1]
                                   * u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                      [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901][(*(*k))]
                                   [2]);
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][1] = -dt * tx2
                        * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901][(*(*k))][2]
                            * (*(*tmp1)));
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][2] = -dt * tx2
                            * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901][(*(*k))][1]
                                * (*(*tmp1)))
                        - dt * tx1 * ((*(*c34)) * (*(*tmp1))) - dt * tx1 * dx3;
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][3] = 0.;
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][2][4] = 0.;
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][0] = -dt * tx2
                            * (-(u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901][(*(*k))]
                                  [1]
                                   * u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                      [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901][(*(*k))]
                                   [3]);
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][1] = -dt * tx2
                        * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901][(*(*k))][3]
                            * (*(*tmp1)));
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][2] = 0.;
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][3] = -dt * tx2
                            * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901][(*(*k))][1]
                                * (*(*tmp1)))
                        - dt * tx1 * ((*(*c34)) * (*(*tmp1))) - dt * tx1 * dx4;
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][3][4] = 0.;
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][0] = -dt * tx2
                            * ((0.40000000000000002
                                       * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                           [(*(*k))][1]
                                               * u[(i_darts_counter_temp4901)-1]
                                                  [j_darts_counter_temp4901][(*(*k))][1]
                                           + u[(i_darts_counter_temp4901)-1]
                                              [j_darts_counter_temp4901][(*(*k))][2]
                                               * u[(i_darts_counter_temp4901)-1]
                                                  [j_darts_counter_temp4901][(*(*k))][2]
                                           + u[(i_darts_counter_temp4901)-1]
                                              [j_darts_counter_temp4901][(*(*k))][3]
                                               * u[(i_darts_counter_temp4901)-1]
                                                  [j_darts_counter_temp4901][(*(*k))][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                           [(*(*k))][4]
                                           * (*(*tmp1))))
                                * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1
                            * (-((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                          [(*(*k))][1])
                                        * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                            [(*(*k))][1])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                          [(*(*k))][2])
                                        * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                            [(*(*k))][2])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                          [(*(*k))][3])
                                        * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                            [(*(*k))][3])))
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                       [(*(*k))][4]);
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][1] = -dt * tx2
                            * (1.3999999999999999
                                    * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                        [(*(*k))][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((3.
                                               * u[(i_darts_counter_temp4901)-1]
                                                  [j_darts_counter_temp4901][(*(*k))][1]
                                               * u[(i_darts_counter_temp4901)-1]
                                                  [j_darts_counter_temp4901][(*(*k))][1]
                                           + u[(i_darts_counter_temp4901)-1]
                                              [j_darts_counter_temp4901][(*(*k))][2]
                                               * u[(i_darts_counter_temp4901)-1]
                                                  [j_darts_counter_temp4901][(*(*k))][2]
                                           + u[(i_darts_counter_temp4901)-1]
                                              [j_darts_counter_temp4901][(*(*k))][3]
                                               * u[(i_darts_counter_temp4901)-1]
                                                  [j_darts_counter_temp4901][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * tx1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901][(*(*k))][1];
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][2] = -dt * tx2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                    [(*(*k))][2]
                                    * u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                       [(*(*k))][1])
                                * (*(*tmp2)))
                        - dt * tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901][(*(*k))][2];
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][3] = -dt * tx2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                    [(*(*k))][3]
                                    * u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                       [(*(*k))][1])
                                * (*(*tmp2)))
                        - dt * tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901][(*(*k))][3];
                    c[(i_darts_counter_temp4901)][j_darts_counter_temp4901][4][4] = -dt * tx2
                            * (1.3999999999999999
                                * (u[(i_darts_counter_temp4901)-1][j_darts_counter_temp4901]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1 * (*(*c1345)) * (*(*tmp1)) - dt * tx1 * dx5;
                }
                (*j) = j_darts_counter_temp4901;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets4901[0].decDep();
}
TP4901::TP4901(int in_numThreads, int in_mainCodeletID, TP7* in_TPParent, int in_initIteration,
    int in_lastIteration, TP4901** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , c1345_darts4901(new double*[this->numThreads])
    , c34_darts4901(new double*[this->numThreads])
    , i_darts4901(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts4901(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts4901(new int*[this->numThreads])
    , r43_darts4901(new double*[this->numThreads])
    , tmp1_darts4901(new double*[this->numThreads])
    , tmp2_darts4901(new double*[this->numThreads])
    , tmp3_darts4901(new double*[this->numThreads])
    , initIteration4901(in_initIteration)
    , lastIteration4901(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets4901(new _barrierCodelets4901[1])
    , checkInCodelets4902(new _checkInCodelets4902[this->numThreads])
{
    /*Initialize the loop parameters*/
    range4901 = abs(lastIteration4901 - initIteration4901) / 1;
    rangePerCodelet4901 = range4901 / numThreads;
    minIteration4901 = min<int>(lastIteration4901, initIteration4901);
    remainderRange4901 = range4901 % numThreads;
    /*Initialize inputs and vars.*/
    this->c1345_darts4901
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->c34_darts4901
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->k_darts4901 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->r43_darts4901
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp1_darts4901
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp2_darts4901
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp3_darts4901
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets4901[0] = _barrierCodelets4901(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets4902* checkInCodelets4902Ptr = (this->checkInCodelets4902);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets4902);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets4902Ptr) = _checkInCodelets4902(2, 1, this, codeletCounter);
#else
        (*checkInCodelets4902Ptr) = _checkInCodelets4902(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets4902Ptr).decDep();
        checkInCodelets4902Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP4901::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets4902[localID].setID(codeletID);
    this->checkInCodelets4902[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets4902[localID + this->baseNumThreads * i]
            = _checkInCodelets4902(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets4902[localID + this->baseNumThreads * i]
            = _checkInCodelets4902(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets4902[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets4902[localID + this->baseNumThreads * i].decDep();
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
TP4901::~TP4901()
{
    delete[] c1345_darts4901;
    delete[] c34_darts4901;
    delete[] k_darts4901;
    delete[] r43_darts4901;
    delete[] tmp1_darts4901;
    delete[] tmp2_darts4901;
    delete[] tmp3_darts4901;
    delete[] barrierCodelets4901;
    delete[] checkInCodelets4902;
}
/*TP8: TP_jacu*/
void TP8::_checkInCodelets7399::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 7399: DeclStmt*/

    /*printing node 7400: DeclStmt*/

    /*printing node 7401: DeclStmt*/

    /*printing node 7402: DeclStmt*/

    /*printing node 7403: DeclStmt*/

    /*printing node 7404: BinaryOperator*/
    (this->inputsTPParent->r43_darts8[this->getID()]) = (4. / 3.);

    /*printing node 7408: BinaryOperator*/
    (this->inputsTPParent->c1345_darts8[this->getID()])
        = 1.3999999999999999 * 0.10000000000000001 * 1. * 1.3999999999999999;

    /*printing node 7416: BinaryOperator*/
    (this->inputsTPParent->c34_darts8[this->getID()]) = 0.10000000000000001 * 1.;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 7399 nextRegion: 7420 */
    myTP->controlTPParent->checkInCodelets7420[this->getID()].decDep();
}
void TP8::_checkInCodelets7420::fire(void)
{
    /*region 7420 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP7420;
    if (idx < myTP->TPsToUse7420) {
        if (!__sync_val_compare_and_swap(&(myTP->TP7420_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(ist - iend) / 1;
            int rangePerCodelet = range / myTP->TPsToUse7420;
            int minIteration = min<int>(ist, iend);
            int remainderRange = range % myTP->TPsToUse7420;
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
            if (idx == myTP->TPsToUse7420 - 1) {
                lastIteration = ist;
            }
#if USEINVOKE == 1
            invoke<TP7420>(myTP, myTP->codeletsPerTP7420 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP7420Ptr[idx]));
#else
            place<TP7420>(idx, myTP, myTP->codeletsPerTP7420 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP7420Ptr[idx]));
#endif
        } else {
            if (myTP->TP7420Ptr[idx] != nullptr) {
                myTP->TP7420Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP8::_barrierCodelets7420::fire(void)
{
    TP8* myTP = static_cast<TP8*>(myTP_);

    for (size_t codeletsCounter = 0; codeletsCounter < (size_t)myTP->numThreads;
         codeletsCounter++) {
        myTP->nextCodeletsjacu[codeletsCounter]->decDep();
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
    , TP7420Ptr(new TP7420*[NUMTPS7420])
    , TP7420_alreadyLaunched(new size_t[NUMTPS7420])
    , numTPsSet7420(0)
    , numTPsReady7420(0)
    , TPsToUse7420(NUMTPS7420)
    , codeletsPerTP7420(this->numThreads / NUMTPS7420)
    , totalCodelets7420(this->TPsToUse7420 * this->codeletsPerTP7420)
    , checkInCodelets7399(new _checkInCodelets7399[this->numThreads])
    , checkInCodelets7420(new _checkInCodelets7420[this->numThreads])
    , barrierCodelets7420(new _barrierCodelets7420[1])
{
    barrierCodelets7420[0] = _barrierCodelets7420(NUMTPS7420, NUMTPS7420, this, 0);
    _checkInCodelets7420* checkInCodelets7420Ptr = (this->checkInCodelets7420);
    for (int i = 0; i < NUMTPS7420; i++) {
        TP7420Ptr[i] = nullptr;
        TP7420_alreadyLaunched[i] = 0;
    }
    _checkInCodelets7399* checkInCodelets7399Ptr = (this->checkInCodelets7399);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets7399);
#endif
    for (size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++) {
        (*checkInCodelets7420Ptr) = _checkInCodelets7420(1, 1, this, codeletCounter);
        checkInCodelets7420Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets7399Ptr) = _checkInCodelets7399(2, 1, this, codeletCounter);
#else
        (*checkInCodelets7399Ptr) = _checkInCodelets7399(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets7399Ptr).decDep();
        checkInCodelets7399Ptr++;
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
    delete[] barrierCodelets7420;
    delete[] checkInCodelets7420;
    delete[] checkInCodelets7399;
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
/*TP7420: OMPForDirective*/
void TP7420::_barrierCodelets7420::fire(void)
{
    TP7420* myTP = static_cast<TP7420*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets7420[0].decDep();
}
bool TP7420::requestNewRangeIterations7420(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 1*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet7420 * codeletID;
        int tempEndRange = rangePerCodelet7420 * (codeletID + 1);
        if (remainderRange7420 != 0) {
            if (codeletID < (uint32_t)remainderRange7420) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange7420;
                tempEndRange += remainderRange7420;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration7420;
        tempEndRange = tempEndRange * 1 + minIteration7420;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration7420 < lastIteration7420) {
            (this->inputsTPParent->i_darts7420[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts7420[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == 0) {
            *endRange = *endRange - 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration7420;
        }
    }
    return isThereNewIteration;
}
void TP7420::_checkInCodelets7421::fire(void)
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
    this->inputsTPParent->c1345_darts7420[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->c1345_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->c34_darts7420[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->c34_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->k_darts7420[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->k_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->r43_darts7420[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->r43_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp1_darts7420[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp1_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp2_darts7420[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp2_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp3_darts7420[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp3_darts8[this->getID()]);

    /*printing node 7421: ForStmt*/
    /*var: c1345*/
    /*var: c34*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: r43*/
    /*var: tmp1*/
    /*var: tmp2*/
    /*var: tmp3*/
    double** c1345 = &(this->inputsTPParent->c1345_darts7420[this->getLocalID()]);
    (void)c1345 /*OMP_SHARED_PRIVATE*/;
    double** c34 = &(this->inputsTPParent->c34_darts7420[this->getLocalID()]);
    (void)c34 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts7420[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts7420[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** k = &(this->inputsTPParent->k_darts7420[this->getLocalID()]);
    (void)k /*OMP_SHARED_PRIVATE*/;
    double** r43 = &(this->inputsTPParent->r43_darts7420[this->getLocalID()]);
    (void)r43 /*OMP_SHARED_PRIVATE*/;
    double** tmp1 = &(this->inputsTPParent->tmp1_darts7420[this->getLocalID()]);
    (void)tmp1 /*OMP_SHARED_PRIVATE*/;
    double** tmp2 = &(this->inputsTPParent->tmp2_darts7420[this->getLocalID()]);
    (void)tmp2 /*OMP_SHARED_PRIVATE*/;
    double** tmp3 = &(this->inputsTPParent->tmp3_darts7420[this->getLocalID()]);
    (void)tmp3 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations7420(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets7420[0].decDep();
        return;
    }
    for (int i_darts_counter_temp7420 = (*i); i_darts_counter_temp7420 >= endRange
         && i_darts_counter_temp7420 >= this->inputsTPParent->lastIteration7420;
         i_darts_counter_temp7420--) {
        {
            {
                /*Loop's init*/
                (*j) = jend;
                int j_darts_counter_temp7420 = (*j);
                for (; j_darts_counter_temp7420 >= jst; j_darts_counter_temp7420--) {
                    (*(*tmp1))
                        = 1. / u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][0]
                        = 1. + dt * 2. * (tx1 * dx1 + ty1 * dy1 + tz1 * dz1);
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][1] = 0.;
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][2] = 0.;
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][3] = 0.;
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][4] = 0.;
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][0] = dt * 2.
                        * (tx1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                       [(*(*k))][1])
                            + ty1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                       [(*(*k))][1])
                            + tz1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                       [(*(*k))][1]));
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][1] = 1.
                        + dt * 2.
                            * (tx1 * (*(*r43)) * (*(*c34)) * (*(*tmp1))
                                + ty1 * (*(*c34)) * (*(*tmp1)) + tz1 * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx2 + ty1 * dy2 + tz1 * dz2);
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][2] = 0.;
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][3] = 0.;
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][4] = 0.;
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][0] = dt * 2.
                        * (tx1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                       [(*(*k))][2])
                            + ty1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                       [(*(*k))][2])
                            + tz1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                       [(*(*k))][2]));
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][1] = 0.;
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][2] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c34)) * (*(*tmp1))
                                + ty1 * (*(*r43)) * (*(*c34)) * (*(*tmp1))
                                + tz1 * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx3 + ty1 * dy3 + tz1 * dz3);
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][3] = 0.;
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][4] = 0.;
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][0] = dt * 2.
                        * (tx1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                       [(*(*k))][3])
                            + ty1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                       [(*(*k))][3])
                            + tz1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                       [(*(*k))][3]));
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][1] = 0.;
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][2] = 0.;
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][3] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c34)) * (*(*tmp1)) + ty1 * (*(*c34)) * (*(*tmp1))
                                + tz1 * (*(*r43)) * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx4 + ty1 * dy4 + tz1 * dz4);
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][4] = 0.;
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][0] = dt * 2.
                        * (tx1
                                * (-((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp7420)]
                                                [j_darts_counter_temp7420][(*(*k))][1])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp7420)]
                                                [j_darts_counter_temp7420][(*(*k))][2])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp7420)]
                                                [j_darts_counter_temp7420][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                           [(*(*k))][4])
                            + ty1
                                * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp7420)]
                                                [j_darts_counter_temp7420][(*(*k))][1])))
                                    - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp7420)]
                                                [j_darts_counter_temp7420][(*(*k))][2])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp7420)]
                                                [j_darts_counter_temp7420][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                           [(*(*k))][4])
                            + tz1
                                * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp7420)]
                                                [j_darts_counter_temp7420][(*(*k))][1])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp7420)]
                                                [j_darts_counter_temp7420][(*(*k))][2])))
                                    - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp7420)]
                                                [j_darts_counter_temp7420][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                           [(*(*k))][4]));
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][1] = dt * 2.
                        * (tx1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k))]
                                   [1]
                            + ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k))]
                                   [1]
                            + tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k))]
                                   [1]);
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][2] = dt * 2.
                        * (tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k))]
                                   [2]
                            + ty1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k))]
                                   [2]
                            + tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k))]
                                   [2]);
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][3] = dt * 2.
                        * (tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k))]
                                   [3]
                            + ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k))]
                                   [3]
                            + tz1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k))]
                                   [3]);
                    d[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][4] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c1345)) * (*(*tmp1)) + ty1 * (*(*c1345)) * (*(*tmp1))
                                + tz1 * (*(*c1345)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx5 + ty1 * dy5 + tz1 * dz5);
                    (*(*tmp1)) = 1.
                        / u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][0] = -dt * tx1 * dx1;
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][1] = dt * tx2;
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][2] = 0.;
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][3] = 0.;
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][4] = 0.;
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][0] = dt * tx2
                            * (-(u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                  [(*(*k))][1]
                                   * (*(*tmp1)))
                                    * (u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                        [(*(*k))][1]
                                        * (*(*tmp1)))
                                + 0.40000000000000002 * 0.5
                                    * (u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                        [(*(*k))][1]
                                            * u[(i_darts_counter_temp7420) + 1]
                                               [j_darts_counter_temp7420][(*(*k))][1]
                                        + u[(i_darts_counter_temp7420) + 1]
                                           [j_darts_counter_temp7420][(*(*k))][2]
                                            * u[(i_darts_counter_temp7420) + 1]
                                               [j_darts_counter_temp7420][(*(*k))][2]
                                        + u[(i_darts_counter_temp7420) + 1]
                                           [j_darts_counter_temp7420][(*(*k))][3]
                                            * u[(i_darts_counter_temp7420) + 1]
                                               [j_darts_counter_temp7420][(*(*k))][3])
                                    * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                   [(*(*k))][1]);
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][1] = dt * tx2
                            * ((2. - 0.40000000000000002)
                                * (u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * tx1 * dx2;
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][2] = dt * tx2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420][(*(*k))]
                                [2]
                                * (*(*tmp1))));
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][3] = dt * tx2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420][(*(*k))]
                                [3]
                                * (*(*tmp1))));
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][4]
                        = dt * tx2 * 0.40000000000000002;
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][0] = dt * tx2
                            * (-(u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                  [(*(*k))][1]
                                   * u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                      [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                   [(*(*k))][2]);
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][1] = dt * tx2
                        * (u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420][(*(*k))][2]
                            * (*(*tmp1)));
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][2] = dt * tx2
                            * (u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420][(*(*k))]
                                [1]
                                * (*(*tmp1)))
                        - dt * tx1 * ((*(*c34)) * (*(*tmp1))) - dt * tx1 * dx3;
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][3] = 0.;
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][4] = 0.;
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][0] = dt * tx2
                            * (-(u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                  [(*(*k))][1]
                                   * u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                      [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                   [(*(*k))][3]);
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][1] = dt * tx2
                        * (u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420][(*(*k))][3]
                            * (*(*tmp1)));
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][2] = 0.;
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][3] = dt * tx2
                            * (u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420][(*(*k))]
                                [1]
                                * (*(*tmp1)))
                        - dt * tx1 * ((*(*c34)) * (*(*tmp1))) - dt * tx1 * dx4;
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][4] = 0.;
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][0] = dt * tx2
                            * ((0.40000000000000002
                                       * (u[(i_darts_counter_temp7420) + 1]
                                           [j_darts_counter_temp7420][(*(*k))][1]
                                               * u[(i_darts_counter_temp7420) + 1]
                                                  [j_darts_counter_temp7420][(*(*k))][1]
                                           + u[(i_darts_counter_temp7420) + 1]
                                              [j_darts_counter_temp7420][(*(*k))][2]
                                               * u[(i_darts_counter_temp7420) + 1]
                                                  [j_darts_counter_temp7420][(*(*k))][2]
                                           + u[(i_darts_counter_temp7420) + 1]
                                              [j_darts_counter_temp7420][(*(*k))][3]
                                               * u[(i_darts_counter_temp7420) + 1]
                                                  [j_darts_counter_temp7420][(*(*k))][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(i_darts_counter_temp7420) + 1]
                                           [j_darts_counter_temp7420][(*(*k))][4]
                                           * (*(*tmp1))))
                                * (u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1
                            * (-((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                          [(*(*k))][1])
                                        * (u[(i_darts_counter_temp7420) + 1]
                                            [j_darts_counter_temp7420][(*(*k))][1])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                          [(*(*k))][2])
                                        * (u[(i_darts_counter_temp7420) + 1]
                                            [j_darts_counter_temp7420][(*(*k))][2])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                          [(*(*k))][3])
                                        * (u[(i_darts_counter_temp7420) + 1]
                                            [j_darts_counter_temp7420][(*(*k))][3])))
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                       [(*(*k))][4]);
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][1] = dt * tx2
                            * (1.3999999999999999
                                    * (u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                        [(*(*k))][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((3.
                                               * u[(i_darts_counter_temp7420) + 1]
                                                  [j_darts_counter_temp7420][(*(*k))][1]
                                               * u[(i_darts_counter_temp7420) + 1]
                                                  [j_darts_counter_temp7420][(*(*k))][1]
                                           + u[(i_darts_counter_temp7420) + 1]
                                              [j_darts_counter_temp7420][(*(*k))][2]
                                               * u[(i_darts_counter_temp7420) + 1]
                                                  [j_darts_counter_temp7420][(*(*k))][2]
                                           + u[(i_darts_counter_temp7420) + 1]
                                              [j_darts_counter_temp7420][(*(*k))][3]
                                               * u[(i_darts_counter_temp7420) + 1]
                                                  [j_darts_counter_temp7420][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * tx1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420][(*(*k))]
                               [1];
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][2] = dt * tx2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                    [(*(*k))][2]
                                    * u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                       [(*(*k))][1])
                                * (*(*tmp2)))
                        - dt * tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420][(*(*k))]
                               [2];
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][3] = dt * tx2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                    [(*(*k))][3]
                                    * u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                       [(*(*k))][1])
                                * (*(*tmp2)))
                        - dt * tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420][(*(*k))]
                               [3];
                    a[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][4] = dt * tx2
                            * (1.3999999999999999
                                * (u[(i_darts_counter_temp7420) + 1][j_darts_counter_temp7420]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1 * (*(*c1345)) * (*(*tmp1)) - dt * tx1 * dx5;
                    (*(*tmp1)) = 1.
                        / u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][0] = -dt * ty1 * dy1;
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][1] = 0.;
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][2] = dt * ty2;
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][3] = 0.;
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][4] = 0.;
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][0] = dt * ty2
                            * (-(u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                  [(*(*k))][1]
                                   * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                      [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * ty1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                   [(*(*k))][1]);
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][1] = dt * ty2
                            * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1][(*(*k))]
                                [2]
                                * (*(*tmp1)))
                        - dt * ty1 * ((*(*c34)) * (*(*tmp1))) - dt * ty1 * dy2;
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][2] = dt * ty2
                        * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1][(*(*k))][1]
                            * (*(*tmp1)));
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][3] = 0.;
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][4] = 0.;
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][0] = dt * ty2
                            * (-(u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                  [(*(*k))][2]
                                   * (*(*tmp1)))
                                    * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                        [(*(*k))][2]
                                        * (*(*tmp1)))
                                + 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                         [(*(*k))][1]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420 + 1][(*(*k))][1]
                                           + u[(i_darts_counter_temp7420)]
                                              [j_darts_counter_temp7420 + 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420 + 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp7420)]
                                              [j_darts_counter_temp7420 + 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420 + 1][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * ty1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                   [(*(*k))][2]);
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][1] = dt * ty2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1][(*(*k))]
                                [1]
                                * (*(*tmp1))));
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][2] = dt * ty2
                            * ((2. - 0.40000000000000002)
                                * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * ty1 * dy3;
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][3] = dt * ty2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1][(*(*k))]
                                [3]
                                * (*(*tmp1))));
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][4]
                        = dt * ty2 * 0.40000000000000002;
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][0] = dt * ty2
                            * (-(u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                  [(*(*k))][2]
                                   * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                      [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * ty1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                   [(*(*k))][3]);
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][1] = 0.;
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][2] = dt * ty2
                        * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1][(*(*k))][3]
                            * (*(*tmp1)));
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][3] = dt * ty2
                            * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1][(*(*k))]
                                [2]
                                * (*(*tmp1)))
                        - dt * ty1 * ((*(*c34)) * (*(*tmp1))) - dt * ty1 * dy4;
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][4] = 0.;
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][0] = dt * ty2
                            * ((0.40000000000000002
                                       * (u[(i_darts_counter_temp7420)]
                                           [j_darts_counter_temp7420 + 1][(*(*k))][1]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420 + 1][(*(*k))][1]
                                           + u[(i_darts_counter_temp7420)]
                                              [j_darts_counter_temp7420 + 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420 + 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp7420)]
                                              [j_darts_counter_temp7420 + 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420 + 1][(*(*k))][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(i_darts_counter_temp7420)]
                                           [j_darts_counter_temp7420 + 1][(*(*k))][4]
                                           * (*(*tmp1))))
                                * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1
                            * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                          [(*(*k))][1])
                                        * (u[(i_darts_counter_temp7420)]
                                            [j_darts_counter_temp7420 + 1][(*(*k))][1])))
                                - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                          [(*(*k))][2])
                                        * (u[(i_darts_counter_temp7420)]
                                            [j_darts_counter_temp7420 + 1][(*(*k))][2])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                          [(*(*k))][3])
                                        * (u[(i_darts_counter_temp7420)]
                                            [j_darts_counter_temp7420 + 1][(*(*k))][3])))
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                       [(*(*k))][4]);
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][1] = dt * ty2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                    [(*(*k))][1]
                                    * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                       [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1][(*(*k))]
                               [1];
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][2] = dt * ty2
                            * (1.3999999999999999
                                    * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                        [(*(*k))][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                         [(*(*k))][1]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420 + 1][(*(*k))][1]
                                           + 3.
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420 + 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420 + 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp7420)]
                                              [j_darts_counter_temp7420 + 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420 + 1][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * ty1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1][(*(*k))]
                               [2];
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][3] = dt * ty2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                    [(*(*k))][2]
                                    * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                       [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1][(*(*k))]
                               [3];
                    b[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][4] = dt * ty2
                            * (1.3999999999999999
                                * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420 + 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1 * (*(*c1345)) * (*(*tmp1)) - dt * ty1 * dy5;
                    (*(*tmp1)) = 1.
                        / u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k)) + 1][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][0] = -dt * tz1 * dz1;
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][1] = 0.;
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][2] = 0.;
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][3] = dt * tz2;
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][0][4] = 0.;
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][0] = dt * tz2
                            * (-(u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                  [(*(*k)) + 1][1]
                                   * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                      [(*(*k)) + 1][3])
                                * (*(*tmp2)))
                        - dt * tz1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                   [(*(*k)) + 1][1]);
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][1] = dt * tz2
                            * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k)) + 1]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * (*(*c34)) * (*(*tmp1)) - dt * tz1 * dz2;
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][2] = 0.;
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][3] = dt * tz2
                        * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k)) + 1][1]
                            * (*(*tmp1)));
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][1][4] = 0.;
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][0] = dt * tz2
                            * (-(u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                  [(*(*k)) + 1][2]
                                   * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                      [(*(*k)) + 1][3])
                                * (*(*tmp2)))
                        - dt * tz1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                   [(*(*k)) + 1][2]);
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][1] = 0.;
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][2] = dt * tz2
                            * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k)) + 1]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * ((*(*c34)) * (*(*tmp1))) - dt * tz1 * dz3;
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][3] = dt * tz2
                        * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k)) + 1][2]
                            * (*(*tmp1)));
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][2][4] = 0.;
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][0] = dt * tz2
                            * (-(u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                  [(*(*k)) + 1][3]
                                   * (*(*tmp1)))
                                    * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                        [(*(*k)) + 1][3]
                                        * (*(*tmp1)))
                                + 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                         [(*(*k)) + 1][1]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420][(*(*k)) + 1][1]
                                           + u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                              [(*(*k)) + 1][2]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420][(*(*k)) + 1][2]
                                           + u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                              [(*(*k)) + 1][3]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420][(*(*k)) + 1][3])
                                        * (*(*tmp2))))
                        - dt * tz1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                   [(*(*k)) + 1][3]);
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][1] = dt * tz2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k)) + 1]
                                [1]
                                * (*(*tmp1))));
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][2] = dt * tz2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k)) + 1]
                                [2]
                                * (*(*tmp1))));
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][3] = dt * tz2
                            * (2. - 0.40000000000000002)
                            * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k)) + 1]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * tz1 * dz4;
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][3][4]
                        = dt * tz2 * 0.40000000000000002;
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][0] = dt * tz2
                            * ((0.40000000000000002
                                       * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                           [(*(*k)) + 1][1]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420][(*(*k)) + 1][1]
                                           + u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                              [(*(*k)) + 1][2]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420][(*(*k)) + 1][2]
                                           + u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                              [(*(*k)) + 1][3]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420][(*(*k)) + 1][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                           [(*(*k)) + 1][4]
                                           * (*(*tmp1))))
                                * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                    [(*(*k)) + 1][3]
                                    * (*(*tmp1))))
                        - dt * tz1
                            * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                          [(*(*k)) + 1][1])
                                        * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                            [(*(*k)) + 1][1])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                          [(*(*k)) + 1][2])
                                        * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                            [(*(*k)) + 1][2])))
                                - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                          [(*(*k)) + 1][3])
                                        * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                            [(*(*k)) + 1][3])))
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                       [(*(*k)) + 1][4]);
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][1] = dt * tz2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                    [(*(*k)) + 1][1]
                                    * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                       [(*(*k)) + 1][3])
                                * (*(*tmp2)))
                        - dt * tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k)) + 1]
                               [1];
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][2] = dt * tz2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                    [(*(*k)) + 1][2]
                                    * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                       [(*(*k)) + 1][3])
                                * (*(*tmp2)))
                        - dt * tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k)) + 1]
                               [2];
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][3] = dt * tz2
                            * (1.3999999999999999
                                    * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                        [(*(*k)) + 1][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                         [(*(*k)) + 1][1]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420][(*(*k)) + 1][1]
                                           + u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                              [(*(*k)) + 1][2]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420][(*(*k)) + 1][2]
                                           + 3.
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420][(*(*k)) + 1][3]
                                               * u[(i_darts_counter_temp7420)]
                                                  [j_darts_counter_temp7420][(*(*k)) + 1][3])
                                        * (*(*tmp2))))
                        - dt * tz1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7420)][j_darts_counter_temp7420][(*(*k)) + 1]
                               [3];
                    c[(i_darts_counter_temp7420)][j_darts_counter_temp7420][4][4] = dt * tz2
                            * (1.3999999999999999
                                * (u[(i_darts_counter_temp7420)][j_darts_counter_temp7420]
                                    [(*(*k)) + 1][3]
                                    * (*(*tmp1))))
                        - dt * tz1 * (*(*c1345)) * (*(*tmp1)) - dt * tz1 * dz5;
                }
                (*j) = j_darts_counter_temp7420;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets7420[0].decDep();
}
TP7420::TP7420(int in_numThreads, int in_mainCodeletID, TP8* in_TPParent, int in_initIteration,
    int in_lastIteration, TP7420** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , c1345_darts7420(new double*[this->numThreads])
    , c34_darts7420(new double*[this->numThreads])
    , i_darts7420(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts7420(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts7420(new int*[this->numThreads])
    , r43_darts7420(new double*[this->numThreads])
    , tmp1_darts7420(new double*[this->numThreads])
    , tmp2_darts7420(new double*[this->numThreads])
    , tmp3_darts7420(new double*[this->numThreads])
    , initIteration7420(in_initIteration)
    , lastIteration7420(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets7420(new _barrierCodelets7420[1])
    , checkInCodelets7421(new _checkInCodelets7421[this->numThreads])
{
    /*Initialize the loop parameters*/
    range7420 = abs(lastIteration7420 - initIteration7420) / 1;
    rangePerCodelet7420 = range7420 / numThreads;
    minIteration7420 = min<int>(lastIteration7420, initIteration7420);
    remainderRange7420 = range7420 % numThreads;
    /*Initialize inputs and vars.*/
    this->c1345_darts7420
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->c34_darts7420
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->k_darts7420 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->r43_darts7420
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp1_darts7420
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp2_darts7420
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp3_darts7420
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets7420[0] = _barrierCodelets7420(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets7421* checkInCodelets7421Ptr = (this->checkInCodelets7421);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets7421);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets7421Ptr) = _checkInCodelets7421(2, 1, this, codeletCounter);
#else
        (*checkInCodelets7421Ptr) = _checkInCodelets7421(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets7421Ptr).decDep();
        checkInCodelets7421Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP7420::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets7421[localID].setID(codeletID);
    this->checkInCodelets7421[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets7421[localID + this->baseNumThreads * i]
            = _checkInCodelets7421(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets7421[localID + this->baseNumThreads * i]
            = _checkInCodelets7421(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets7421[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets7421[localID + this->baseNumThreads * i].decDep();
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
TP7420::~TP7420()
{
    delete[] c1345_darts7420;
    delete[] c34_darts7420;
    delete[] k_darts7420;
    delete[] r43_darts7420;
    delete[] tmp1_darts7420;
    delete[] tmp2_darts7420;
    delete[] tmp3_darts7420;
    delete[] barrierCodelets7420;
    delete[] checkInCodelets7421;
}
/*TP9871: OMPParallelDirective*/
void TP9871::_barrierCodelets9871::fire(void)
{
    TP9871* myTP = static_cast<TP9871*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP9871::_checkInCodelets9873::fire(void)
{
    /*Init the vars for this region*/

    /*printing node 9873: DeclStmt*/

    /*printing node 9874: DeclStmt*/
    this->inputsTPParent->sum0_darts9871[this->getID()] = 0.;
    this->inputsTPParent->sum1_darts9871[this->getID()] = 0.;
    this->inputsTPParent->sum2_darts9871[this->getID()] = 0.;
    this->inputsTPParent->sum3_darts9871[this->getID()] = 0.;
    this->inputsTPParent->sum4_darts9871[this->getID()] = 0.;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 9873 nextRegion: 9880 */
    myTP->controlTPParent->checkInCodelets9880[this->getID()].decDep();
}
void TP9871::_checkInCodelets9880::fire(void)
{
    /*Select the thread executing OMPSingleDirective 9880*/
    if (!__sync_val_compare_and_swap(&(myTP->TP9880_alreadyLaunched), 0, 1)) {
        /*Init the vars for this region*/
        /*Initialize the vars of the inlined region*/
        this->inputsTPParent->sum_darts9880
            = (this->inputsTPParent->sum_darts9871) /*OMP_SHARED - VAR INLINED*/;

        /*printing node 9881: ForStmt*/
        {
            /*Loop's init*/
            (this->inputsTPParent->m_darts9880) = 0;
            int m_darts_counter_temp9880 = (this->inputsTPParent->m_darts9880);
            for (; m_darts_counter_temp9880 < 5; m_darts_counter_temp9880++) {
                (*(this->inputsTPParent->sum_darts9880))[m_darts_counter_temp9880] = 0.;
            }
            (this->inputsTPParent->m_darts9880) = m_darts_counter_temp9880;
        }
        /*Signaling next codelet from last stmt in the codelet*/
        /*Signaling omp region's barrier*/
        myTP->controlTPParent->barrierCodelets9880[0].decDep();
    } else {
        /*Signaling omp region's barrier*/
        myTP->barrierCodelets9880[0].decDep();
    }
}
void TP9871::_barrierCodelets9880::fire(void)
{
    TP9871* myTP = static_cast<TP9871*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets9889[codeletsCounter].decDep();
        }
    }
}
void TP9871::_checkInCodelets9889::fire(void)
{
    /*region 9889 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP9889;
    if (idx < myTP->TPsToUse9889) {
        if (!__sync_val_compare_and_swap(&(myTP->TP9889_alreadyLaunched[idx]), 0, 1)) {
            int range = abs((*(this->inputsTPParent->iend_darts9871))
                            - (*(this->inputsTPParent->ist_darts9871)))
                / 1;
            int rangePerCodelet = range / myTP->TPsToUse9889;
            int minIteration = min<int>((*(this->inputsTPParent->iend_darts9871)),
                (*(this->inputsTPParent->ist_darts9871)));
            int remainderRange = range % myTP->TPsToUse9889;
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
            if ((*(this->inputsTPParent->ist_darts9871))
                < (*(this->inputsTPParent->iend_darts9871))) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse9889 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse9889 - 1) {
                lastIteration = (*(this->inputsTPParent->iend_darts9871));
            }
#if USEINVOKE == 1
            invoke<TP9889>(myTP, myTP->codeletsPerTP9889 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(*(this->inputsTPParent->iend_darts9871)),
                &(*(this->inputsTPParent->ist_darts9871)),
                &(*(this->inputsTPParent->jend_darts9871)),
                &(*(this->inputsTPParent->jst_darts9871)),
                &(*(this->inputsTPParent->nz0_darts9871)), &(myTP->TP9889Ptr[idx]));
#else
            place<TP9889>(idx, myTP, myTP->codeletsPerTP9889 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(*(this->inputsTPParent->iend_darts9871)),
                &(*(this->inputsTPParent->ist_darts9871)),
                &(*(this->inputsTPParent->jend_darts9871)),
                &(*(this->inputsTPParent->jst_darts9871)),
                &(*(this->inputsTPParent->nz0_darts9871)), &(myTP->TP9889Ptr[idx]));
#endif
        } else {
            if (myTP->TP9889Ptr[idx] != nullptr) {
                myTP->TP9889Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    } else {
        /*Signaling next codelet region: 9889 nextRegion: 9986 */
        myTP->controlTPParent->checkInCodelets9986[this->getID()].decDep();
    }
}
void TP9871::_checkInCodelets9986::fire(void)
{

    /*printing node 9986: CompoundAssignOperator*/
    TP9984mutex.lock();
    (*(this->inputsTPParent->sum_darts9871))[0]
        += (this->inputsTPParent->sum0_darts9871[this->getID()]);

    /*printing node 9988: CompoundAssignOperator*/
    (*(this->inputsTPParent->sum_darts9871))[1]
        += (this->inputsTPParent->sum1_darts9871[this->getID()]);

    /*printing node 9990: CompoundAssignOperator*/
    (*(this->inputsTPParent->sum_darts9871))[2]
        += (this->inputsTPParent->sum2_darts9871[this->getID()]);

    /*printing node 9992: CompoundAssignOperator*/
    (*(this->inputsTPParent->sum_darts9871))[3]
        += (this->inputsTPParent->sum3_darts9871[this->getID()]);

    /*printing node 9994: CompoundAssignOperator*/
    (*(this->inputsTPParent->sum_darts9871))[4]
        += (this->inputsTPParent->sum4_darts9871[this->getID()]);
    TP9984mutex.unlock();
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 9986 nextRegion: 9996 */
    myTP->controlTPParent->barrierCodelets9996[0].decDep();
}
void TP9871::_barrierCodelets9996::fire(void)
{
    TP9871* myTP = static_cast<TP9871*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets9997[codeletsCounter].decDep();
        }
    }
}
void TP9871::_checkInCodelets9997::fire(void)
{
    /*Select the thread executing OMPSingleDirective 9997*/
    if (!__sync_val_compare_and_swap(&(myTP->TP9997_alreadyLaunched), 0, 1)) {
        /*Init the vars for this region*/
        /*Initialize the vars of the inlined region*/
        this->inputsTPParent->nx0_darts9997
            = (this->inputsTPParent->nx0_darts9871) /*OMP_SHARED - VAR INLINED*/;
        this->inputsTPParent->ny0_darts9997
            = (this->inputsTPParent->ny0_darts9871) /*OMP_SHARED - VAR INLINED*/;
        this->inputsTPParent->nz0_darts9997
            = (this->inputsTPParent->nz0_darts9871) /*OMP_SHARED - VAR INLINED*/;
        this->inputsTPParent->sum_darts9997
            = (this->inputsTPParent->sum_darts9871) /*OMP_SHARED - VAR INLINED*/;

        /*printing node 9998: ForStmt*/
        {
            /*Loop's init*/
            (this->inputsTPParent->m_darts9997) = 0;
            int m_darts_counter_temp9997 = (this->inputsTPParent->m_darts9997);
            for (; m_darts_counter_temp9997 < 5; m_darts_counter_temp9997++) {
                (*(this->inputsTPParent->sum_darts9997))[m_darts_counter_temp9997]
                    = sqrt((*(this->inputsTPParent->sum_darts9997))[m_darts_counter_temp9997]
                        / (((*(this->inputsTPParent->nx0_darts9997)) - 2)
                            * ((*(this->inputsTPParent->ny0_darts9997)) - 2)
                            * ((*(this->inputsTPParent->nz0_darts9997)) - 2)));
            }
            (this->inputsTPParent->m_darts9997) = m_darts_counter_temp9997;
        }
        /*Signaling next codelet from last stmt in the codelet*/
        /*Signaling omp region's barrier*/
        myTP->controlTPParent->barrierCodelets9997[0].decDep();
    } else {
        /*Signaling omp region's barrier*/
        myTP->barrierCodelets9997[0].decDep();
    }
}
void TP9871::_barrierCodelets9997::fire(void)
{
    TP9871* myTP = static_cast<TP9871*>(myTP_);
    myTP->TPParent->barrierCodelets9871[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets9871[0]));
}
TP9871::TP9871(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet, int* in_iend,
    int* in_ist, int* in_jend, int* in_jst, int* in_nx0, int* in_ny0, int* in_nz0, double** in_sum)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , iend_darts9871(in_iend) /*OMP_SHARED - INPUT*/
    , ist_darts9871(in_ist) /*OMP_SHARED - INPUT*/
    , jend_darts9871(in_jend) /*OMP_SHARED - INPUT*/
    , jst_darts9871(in_jst) /*OMP_SHARED - INPUT*/
    , nx0_darts9871(in_nx0) /*OMP_SHARED - INPUT*/
    , ny0_darts9871(in_ny0) /*OMP_SHARED - INPUT*/
    , nz0_darts9871(in_nz0) /*OMP_SHARED - INPUT*/
    , sum_darts9871(in_sum) /*OMP_SHARED - INPUT*/
    , i_darts9871(new int[this->numThreads]) /*VARIABLE*/
    , j_darts9871(new int[this->numThreads]) /*VARIABLE*/
    , k_darts9871(new int[this->numThreads]) /*VARIABLE*/
    , m_darts9871(new int[this->numThreads]) /*VARIABLE*/
    , sum0_darts9871(new double[this->numThreads]) /*VARIABLE*/
    , sum1_darts9871(new double[this->numThreads]) /*VARIABLE*/
    , sum2_darts9871(new double[this->numThreads]) /*VARIABLE*/
    , sum3_darts9871(new double[this->numThreads]) /*VARIABLE*/
    , sum4_darts9871(new double[this->numThreads]) /*VARIABLE*/
    , TP9880_alreadyLaunched(0)
    , TP9889Ptr(new TP9889*[NUMTPS9889])
    , TP9889_alreadyLaunched(new size_t[NUMTPS9889])
    , numTPsSet9889(0)
    , numTPsReady9889(0)
    , TPsToUse9889(NUMTPS9889)
    , codeletsPerTP9889(this->numThreads / NUMTPS9889)
    , totalCodelets9889(this->TPsToUse9889 * this->codeletsPerTP9889)
    , TP9997_alreadyLaunched(0)
    , barrierCodelets9871(new _barrierCodelets9871[1])
    , checkInCodelets9873(new _checkInCodelets9873[this->numThreads])
    , checkInCodelets9880(new _checkInCodelets9880[this->numThreads])
    , barrierCodelets9880(new _barrierCodelets9880[1])
    , checkInCodelets9889(new _checkInCodelets9889[this->numThreads])
    , checkInCodelets9986(new _checkInCodelets9986[this->numThreads])
    , barrierCodelets9996(new _barrierCodelets9996[1])
    , checkInCodelets9997(new _checkInCodelets9997[this->numThreads])
    , barrierCodelets9997(new _barrierCodelets9997[1])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets9871[0] = _barrierCodelets9871(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets9997[0] = _barrierCodelets9997(this->numThreads, this->numThreads, this, 0);
    barrierCodelets9996[0] = _barrierCodelets9996(this->numThreads, this->numThreads, this, 0);
    barrierCodelets9880[0] = _barrierCodelets9880(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets9997* checkInCodelets9997Ptr = (this->checkInCodelets9997);
    _checkInCodelets9986* checkInCodelets9986Ptr = (this->checkInCodelets9986);
    _checkInCodelets9889* checkInCodelets9889Ptr = (this->checkInCodelets9889);
    for (int i = 0; i < NUMTPS9889; i++) {
        TP9889Ptr[i] = nullptr;
        TP9889_alreadyLaunched[i] = 0;
    }
    _checkInCodelets9880* checkInCodelets9880Ptr = (this->checkInCodelets9880);
    _checkInCodelets9873* checkInCodelets9873Ptr = (this->checkInCodelets9873);
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets9997Ptr) = _checkInCodelets9997(1, 1, this, codeletCounter);
        checkInCodelets9997Ptr++;
        (*checkInCodelets9986Ptr) = _checkInCodelets9986(1, 1, this, codeletCounter);
        checkInCodelets9986Ptr++;
        (*checkInCodelets9889Ptr) = _checkInCodelets9889(1, 1, this, codeletCounter);
        checkInCodelets9889Ptr++;
        (*checkInCodelets9880Ptr) = _checkInCodelets9880(1, 1, this, codeletCounter);
        checkInCodelets9880Ptr++;
        (*checkInCodelets9873Ptr) = _checkInCodelets9873(1, 1, this, codeletCounter);
        (*checkInCodelets9873Ptr).decDep();
        checkInCodelets9873Ptr++;
    }
}
TP9871::~TP9871()
{
    delete[] i_darts9871;
    delete[] j_darts9871;
    delete[] k_darts9871;
    delete[] m_darts9871;
    delete[] sum0_darts9871;
    delete[] sum1_darts9871;
    delete[] sum2_darts9871;
    delete[] sum3_darts9871;
    delete[] sum4_darts9871;
    delete[] barrierCodelets9871;
    delete[] barrierCodelets9997;
    delete[] checkInCodelets9997;
    delete[] barrierCodelets9996;
    delete[] checkInCodelets9986;
    delete[] checkInCodelets9889;
    delete[] barrierCodelets9880;
    delete[] checkInCodelets9880;
    delete[] checkInCodelets9873;
}
/*TP9889: OMPForDirective*/
bool TP9889::requestNewRangeIterations9889(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet9889 * codeletID;
        int tempEndRange = rangePerCodelet9889 * (codeletID + 1);
        if (remainderRange9889 != 0) {
            if (codeletID < (uint32_t)remainderRange9889) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange9889;
                tempEndRange += remainderRange9889;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration9889;
        tempEndRange = tempEndRange * 1 + minIteration9889;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration9889 < lastIteration9889) {
            (this->inputsTPParent->i_darts9889[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts9889[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration9889;
        }
    }
    return isThereNewIteration;
}
void TP9889::_checkInCodelets9890::fire(void)
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
    this->inputsTPParent->sum0_darts9889[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->sum0_darts9871[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->sum1_darts9889[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->sum1_darts9871[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->sum2_darts9889[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->sum2_darts9871[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->sum3_darts9889[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->sum3_darts9871[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->sum4_darts9889[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->sum4_darts9871[this->getID()]);

    /*printing node 9890: ForStmt*/
    /*var: i*/
    /*var: iend*/
    /*var: ist*/
    /*var: j*/
    /*var: jend*/
    /*var: jst*/
    /*var: k*/
    /*var: nz0*/
    /*var: sum0*/
    /*var: sum1*/
    /*var: sum2*/
    /*var: sum3*/
    /*var: sum4*/
    int* i = &(this->inputsTPParent->i_darts9889[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts9889[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* jend = (this->inputsTPParent->jend_darts9889);
    (void)jend /*OMP_SHARED*/;
    int* jst = (this->inputsTPParent->jst_darts9889);
    (void)jst /*OMP_SHARED*/;
    int* k = &(this->inputsTPParent->k_darts9889[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* nz0 = (this->inputsTPParent->nz0_darts9889);
    (void)nz0 /*OMP_SHARED*/;
    double** sum0 = &(this->inputsTPParent->sum0_darts9889[this->getLocalID()]);
    (void)sum0 /*OMP_SHARED_PRIVATE*/;
    double** sum1 = &(this->inputsTPParent->sum1_darts9889[this->getLocalID()]);
    (void)sum1 /*OMP_SHARED_PRIVATE*/;
    double** sum2 = &(this->inputsTPParent->sum2_darts9889[this->getLocalID()]);
    (void)sum2 /*OMP_SHARED_PRIVATE*/;
    double** sum3 = &(this->inputsTPParent->sum3_darts9889[this->getLocalID()]);
    (void)sum3 /*OMP_SHARED_PRIVATE*/;
    double** sum4 = &(this->inputsTPParent->sum4_darts9889[this->getLocalID()]);
    (void)sum4 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations9889(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Find and signal the next codelet*/
        myTP->controlTPParent->TPParent->checkInCodelets9986[this->getID()].decDep();
        return;
    }
    for (int i_darts_counter_temp9889 = (*i); i_darts_counter_temp9889 <= endRange
         && i_darts_counter_temp9889 <= this->inputsTPParent->lastIteration9889;
         i_darts_counter_temp9889++) {
        {
            {
                /*Loop's init*/
                (*j) = (*(jst));
                int j_darts_counter_temp9889 = (*j);
                for (; j_darts_counter_temp9889 <= (*(jend)); j_darts_counter_temp9889++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp9889 = (*k);
                        for (; k_darts_counter_temp9889 <= (*(nz0)) - 2;
                             k_darts_counter_temp9889++) {
                            (*(*sum0)) = (*(*sum0))
                                + rsd[(i_darts_counter_temp9889)][j_darts_counter_temp9889]
                                     [k_darts_counter_temp9889][0]
                                    * rsd[(i_darts_counter_temp9889)][j_darts_counter_temp9889]
                                         [k_darts_counter_temp9889][0];
                            (*(*sum1)) = (*(*sum1))
                                + rsd[(i_darts_counter_temp9889)][j_darts_counter_temp9889]
                                     [k_darts_counter_temp9889][1]
                                    * rsd[(i_darts_counter_temp9889)][j_darts_counter_temp9889]
                                         [k_darts_counter_temp9889][1];
                            (*(*sum2)) = (*(*sum2))
                                + rsd[(i_darts_counter_temp9889)][j_darts_counter_temp9889]
                                     [k_darts_counter_temp9889][2]
                                    * rsd[(i_darts_counter_temp9889)][j_darts_counter_temp9889]
                                         [k_darts_counter_temp9889][2];
                            (*(*sum3)) = (*(*sum3))
                                + rsd[(i_darts_counter_temp9889)][j_darts_counter_temp9889]
                                     [k_darts_counter_temp9889][3]
                                    * rsd[(i_darts_counter_temp9889)][j_darts_counter_temp9889]
                                         [k_darts_counter_temp9889][3];
                            (*(*sum4)) = (*(*sum4))
                                + rsd[(i_darts_counter_temp9889)][j_darts_counter_temp9889]
                                     [k_darts_counter_temp9889][4]
                                    * rsd[(i_darts_counter_temp9889)][j_darts_counter_temp9889]
                                         [k_darts_counter_temp9889][4];
                        }
                        (*k) = k_darts_counter_temp9889;
                    }
                }
                (*j) = j_darts_counter_temp9889;
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
    myTP->controlTPParent->TPParent->checkInCodelets9986[this->getID()].decDep();
}
TP9889::TP9889(int in_numThreads, int in_mainCodeletID, TP9871* in_TPParent, int in_initIteration,
    int in_lastIteration, int* in_iend, int* in_ist, int* in_jend, int* in_jst, int* in_nz0,
    TP9889** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts9889(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iend_darts9889(in_iend) /*OMP_SHARED - INPUT*/
    , ist_darts9889(in_ist) /*OMP_SHARED - INPUT*/
    , j_darts9889(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jend_darts9889(in_jend) /*OMP_SHARED - INPUT*/
    , jst_darts9889(in_jst) /*OMP_SHARED - INPUT*/
    , k_darts9889(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , nz0_darts9889(in_nz0) /*OMP_SHARED - INPUT*/
    , sum0_darts9889(new double*[this->numThreads])
    , sum1_darts9889(new double*[this->numThreads])
    , sum2_darts9889(new double*[this->numThreads])
    , sum3_darts9889(new double*[this->numThreads])
    , sum4_darts9889(new double*[this->numThreads])
    , initIteration9889(in_initIteration)
    , lastIteration9889(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , signalNextReady(new int[baseNumThreads])
    , checkInCodelets9890(new _checkInCodelets9890[this->numThreads])
{
    /*Initialize the loop parameters*/
    range9889 = abs(lastIteration9889 - initIteration9889) / 1;
    rangePerCodelet9889 = range9889 / numThreads;
    minIteration9889 = min<int>(lastIteration9889, initIteration9889);
    remainderRange9889 = range9889 % numThreads;
    /*Initialize inputs and vars.*/
    this->sum0_darts9889
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->sum1_darts9889
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->sum2_darts9889
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->sum3_darts9889
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->sum4_darts9889
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    _checkInCodelets9890* checkInCodelets9890Ptr = (this->checkInCodelets9890);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets9890);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
        this->signalNextReady[codeletCounter] = 0;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets9890Ptr) = _checkInCodelets9890(2, 1, this, codeletCounter);
#else
        (*checkInCodelets9890Ptr) = _checkInCodelets9890(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets9890Ptr).decDep();
        checkInCodelets9890Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP9889::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets9890[localID].setID(codeletID);
    this->checkInCodelets9890[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets9890[localID + this->baseNumThreads * i]
            = _checkInCodelets9890(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets9890[localID + this->baseNumThreads * i]
            = _checkInCodelets9890(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets9890[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets9890[localID + this->baseNumThreads * i].decDep();
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
TP9889::~TP9889()
{
    delete[] sum0_darts9889;
    delete[] sum1_darts9889;
    delete[] sum2_darts9889;
    delete[] sum3_darts9889;
    delete[] sum4_darts9889;
    delete[] checkInCodelets9890;
}
/*TP10788: OMPParallelDirective*/
void TP10788::_barrierCodelets10788::fire(void)
{
    TP10788* myTP = static_cast<TP10788*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP10788::_checkInCodelets10803::fire(void)
{
    /*region 10803 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP10803;
    if (idx < myTP->TPsToUse10803) {
        if (!__sync_val_compare_and_swap(&(myTP->TP10803_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(nx - 1 - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse10803;
            int minIteration = min<int>(nx - 1, 0);
            int remainderRange = range % myTP->TPsToUse10803;
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
            if (idx == myTP->TPsToUse10803 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse10803 - 1) {
                lastIteration = nx - 1;
            }
#if USEINVOKE == 1
            invoke<TP10803>(myTP, myTP->codeletsPerTP10803 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP10803Ptr[idx]));
#else
            place<TP10803>(idx, myTP, myTP->codeletsPerTP10803 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP10803Ptr[idx]));
#endif
        } else {
            if (myTP->TP10803Ptr[idx] != nullptr) {
                myTP->TP10803Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10788::_barrierCodelets10803::fire(void)
{
    TP10788* myTP = static_cast<TP10788*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets10860[codeletsCounter].decDep();
        }
    }
}
void TP10788::_checkInCodelets10860::fire(void)
{

    /*printing node 10860: BinaryOperator*/
    (this->inputsTPParent->L1_darts10788[this->getID()]) = 0;

    /*printing node 10861: BinaryOperator*/
    (this->inputsTPParent->L2_darts10788[this->getID()]) = nx - 1;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 10860 nextRegion: 10863 */
    myTP->controlTPParent->checkInCodelets10863[this->getID()].decDep();
}
void TP10788::_checkInCodelets10863::fire(void)
{
    /*region 10863 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP10863;
    if (idx < myTP->TPsToUse10863) {
        if (!__sync_val_compare_and_swap(&(myTP->TP10863_alreadyLaunched[idx]), 0, 1)) {
            int range = abs((this->inputsTPParent->L2_darts10788[this->getID()])
                            - (this->inputsTPParent->L1_darts10788[this->getID()]))
                / 1;
            int rangePerCodelet = range / myTP->TPsToUse10863;
            int minIteration = min<int>((this->inputsTPParent->L2_darts10788[this->getID()]),
                (this->inputsTPParent->L1_darts10788[this->getID()]));
            int remainderRange = range % myTP->TPsToUse10863;
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
            if ((this->inputsTPParent->L1_darts10788[this->getID()])
                < (this->inputsTPParent->L2_darts10788[this->getID()])) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse10863 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse10863 - 1) {
                lastIteration = (this->inputsTPParent->L2_darts10788[this->getID()]);
            }
#if USEINVOKE == 1
            invoke<TP10863>(myTP, myTP->codeletsPerTP10863 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP10863Ptr[idx]));
#else
            place<TP10863>(idx, myTP, myTP->codeletsPerTP10863 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP10863Ptr[idx]));
#endif
        } else {
            if (myTP->TP10863Ptr[idx] != nullptr) {
                myTP->TP10863Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10788::_barrierCodelets10863::fire(void)
{
    TP10788* myTP = static_cast<TP10788*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets11012[codeletsCounter].decDep();
        }
    }
}
void TP10788::_checkInCodelets11012::fire(void)
{
    /*region 11012 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP11012;
    if (idx < myTP->TPsToUse11012) {
        if (!__sync_val_compare_and_swap(&(myTP->TP11012_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(jend - jst) / 1;
            int rangePerCodelet = range / myTP->TPsToUse11012;
            int minIteration = min<int>(jend, jst);
            int remainderRange = range % myTP->TPsToUse11012;
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
            if (idx == myTP->TPsToUse11012 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse11012 - 1) {
                lastIteration = jend;
            }
#if USEINVOKE == 1
            invoke<TP11012>(myTP, myTP->codeletsPerTP11012 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11012Ptr[idx]));
#else
            place<TP11012>(idx, myTP, myTP->codeletsPerTP11012 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11012Ptr[idx]));
#endif
        } else {
            if (myTP->TP11012Ptr[idx] != nullptr) {
                myTP->TP11012Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10788::_barrierCodelets11012::fire(void)
{
    TP10788* myTP = static_cast<TP10788*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets11649[codeletsCounter].decDep();
        }
    }
}
void TP10788::_checkInCodelets11649::fire(void)
{

    /*printing node 11649: BinaryOperator*/
    (this->inputsTPParent->L1_darts10788[this->getID()]) = 0;

    /*printing node 11650: BinaryOperator*/
    (this->inputsTPParent->L2_darts10788[this->getID()]) = ny - 1;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 11649 nextRegion: 11652 */
    myTP->controlTPParent->checkInCodelets11652[this->getID()].decDep();
}
void TP10788::_checkInCodelets11652::fire(void)
{
    /*region 11652 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP11652;
    if (idx < myTP->TPsToUse11652) {
        if (!__sync_val_compare_and_swap(&(myTP->TP11652_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse11652;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse11652;
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
            if (idx == myTP->TPsToUse11652 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse11652 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP11652>(myTP, myTP->codeletsPerTP11652 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11652Ptr[idx]));
#else
            place<TP11652>(idx, myTP, myTP->codeletsPerTP11652 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11652Ptr[idx]));
#endif
        } else {
            if (myTP->TP11652Ptr[idx] != nullptr) {
                myTP->TP11652Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10788::_barrierCodelets11652::fire(void)
{
    TP10788* myTP = static_cast<TP10788*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets11801[codeletsCounter].decDep();
        }
    }
}
void TP10788::_checkInCodelets11801::fire(void)
{
    /*region 11801 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP11801;
    if (idx < myTP->TPsToUse11801) {
        if (!__sync_val_compare_and_swap(&(myTP->TP11801_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse11801;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse11801;
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
            if (idx == myTP->TPsToUse11801 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse11801 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP11801>(myTP, myTP->codeletsPerTP11801 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11801Ptr[idx]));
#else
            place<TP11801>(idx, myTP, myTP->codeletsPerTP11801 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11801Ptr[idx]));
#endif
        } else {
            if (myTP->TP11801Ptr[idx] != nullptr) {
                myTP->TP11801Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10788::_barrierCodelets11801::fire(void)
{
    TP10788* myTP = static_cast<TP10788*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets12438[codeletsCounter].decDep();
        }
    }
}
void TP10788::_checkInCodelets12438::fire(void)
{
    /*region 12438 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP12438;
    if (idx < myTP->TPsToUse12438) {
        if (!__sync_val_compare_and_swap(&(myTP->TP12438_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse12438;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse12438;
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
            if (idx == myTP->TPsToUse12438 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse12438 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP12438>(myTP, myTP->codeletsPerTP12438 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP12438Ptr[idx]));
#else
            place<TP12438>(idx, myTP, myTP->codeletsPerTP12438 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP12438Ptr[idx]));
#endif
        } else {
            if (myTP->TP12438Ptr[idx] != nullptr) {
                myTP->TP12438Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10788::_barrierCodelets12438::fire(void)
{
    TP10788* myTP = static_cast<TP10788*>(myTP_);
    myTP->TPParent->barrierCodelets10788[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets10788[0]));
}
TP10788::TP10788(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts10788(new int[this->numThreads]) /*VARIABLE*/
    , L2_darts10788(new int[this->numThreads]) /*VARIABLE*/
    , i_darts10788(new int[this->numThreads]) /*VARIABLE*/
    , iend1_darts10788(new int[this->numThreads]) /*VARIABLE*/
    , ist1_darts10788(new int[this->numThreads]) /*VARIABLE*/
    , j_darts10788(new int[this->numThreads]) /*VARIABLE*/
    , jend1_darts10788(new int[this->numThreads]) /*VARIABLE*/
    , jst1_darts10788(new int[this->numThreads]) /*VARIABLE*/
    , k_darts10788(new int[this->numThreads]) /*VARIABLE*/
    , m_darts10788(new int[this->numThreads]) /*VARIABLE*/
    , q_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , tmp_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u21_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u21i_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u21im1_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u21j_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u21jm1_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u21k_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u21km1_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u31_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u31i_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u31im1_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u31j_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u31jm1_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u31k_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u31km1_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u41_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u41i_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u41im1_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u41j_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u41jm1_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u41k_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u41km1_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u51i_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u51im1_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u51j_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u51jm1_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u51k_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , u51km1_darts10788(new double[this->numThreads]) /*VARIABLE*/
    , TP10803Ptr(new TP10803*[NUMTPS10803])
    , TP10803_alreadyLaunched(new size_t[NUMTPS10803])
    , numTPsSet10803(0)
    , numTPsReady10803(0)
    , TPsToUse10803(NUMTPS10803)
    , codeletsPerTP10803(this->numThreads / NUMTPS10803)
    , totalCodelets10803(this->TPsToUse10803 * this->codeletsPerTP10803)
    , TP10863Ptr(new TP10863*[NUMTPS10863])
    , TP10863_alreadyLaunched(new size_t[NUMTPS10863])
    , numTPsSet10863(0)
    , numTPsReady10863(0)
    , TPsToUse10863(NUMTPS10863)
    , codeletsPerTP10863(this->numThreads / NUMTPS10863)
    , totalCodelets10863(this->TPsToUse10863 * this->codeletsPerTP10863)
    , TP11012Ptr(new TP11012*[NUMTPS11012])
    , TP11012_alreadyLaunched(new size_t[NUMTPS11012])
    , numTPsSet11012(0)
    , numTPsReady11012(0)
    , TPsToUse11012(NUMTPS11012)
    , codeletsPerTP11012(this->numThreads / NUMTPS11012)
    , totalCodelets11012(this->TPsToUse11012 * this->codeletsPerTP11012)
    , TP11652Ptr(new TP11652*[NUMTPS11652])
    , TP11652_alreadyLaunched(new size_t[NUMTPS11652])
    , numTPsSet11652(0)
    , numTPsReady11652(0)
    , TPsToUse11652(NUMTPS11652)
    , codeletsPerTP11652(this->numThreads / NUMTPS11652)
    , totalCodelets11652(this->TPsToUse11652 * this->codeletsPerTP11652)
    , TP11801Ptr(new TP11801*[NUMTPS11801])
    , TP11801_alreadyLaunched(new size_t[NUMTPS11801])
    , numTPsSet11801(0)
    , numTPsReady11801(0)
    , TPsToUse11801(NUMTPS11801)
    , codeletsPerTP11801(this->numThreads / NUMTPS11801)
    , totalCodelets11801(this->TPsToUse11801 * this->codeletsPerTP11801)
    , TP12438Ptr(new TP12438*[NUMTPS12438])
    , TP12438_alreadyLaunched(new size_t[NUMTPS12438])
    , numTPsSet12438(0)
    , numTPsReady12438(0)
    , TPsToUse12438(NUMTPS12438)
    , codeletsPerTP12438(this->numThreads / NUMTPS12438)
    , totalCodelets12438(this->TPsToUse12438 * this->codeletsPerTP12438)
    , barrierCodelets10788(new _barrierCodelets10788[1])
    , checkInCodelets10803(new _checkInCodelets10803[this->numThreads])
    , barrierCodelets10803(new _barrierCodelets10803[1])
    , checkInCodelets10860(new _checkInCodelets10860[this->numThreads])
    , checkInCodelets10863(new _checkInCodelets10863[this->numThreads])
    , barrierCodelets10863(new _barrierCodelets10863[1])
    , checkInCodelets11012(new _checkInCodelets11012[this->numThreads])
    , barrierCodelets11012(new _barrierCodelets11012[1])
    , checkInCodelets11649(new _checkInCodelets11649[this->numThreads])
    , checkInCodelets11652(new _checkInCodelets11652[this->numThreads])
    , barrierCodelets11652(new _barrierCodelets11652[1])
    , checkInCodelets11801(new _checkInCodelets11801[this->numThreads])
    , barrierCodelets11801(new _barrierCodelets11801[1])
    , checkInCodelets12438(new _checkInCodelets12438[this->numThreads])
    , barrierCodelets12438(new _barrierCodelets12438[1])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets10788[0] = _barrierCodelets10788(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets12438[0] = _barrierCodelets12438(NUMTPS12438, NUMTPS12438, this, 0);
    barrierCodelets11801[0] = _barrierCodelets11801(NUMTPS11801, NUMTPS11801, this, 0);
    barrierCodelets11652[0] = _barrierCodelets11652(NUMTPS11652, NUMTPS11652, this, 0);
    barrierCodelets11012[0] = _barrierCodelets11012(NUMTPS11012, NUMTPS11012, this, 0);
    barrierCodelets10863[0] = _barrierCodelets10863(NUMTPS10863, NUMTPS10863, this, 0);
    barrierCodelets10803[0] = _barrierCodelets10803(NUMTPS10803, NUMTPS10803, this, 0);
    _checkInCodelets12438* checkInCodelets12438Ptr = (this->checkInCodelets12438);
    for (int i = 0; i < NUMTPS12438; i++) {
        TP12438Ptr[i] = nullptr;
        TP12438_alreadyLaunched[i] = 0;
    }
    _checkInCodelets11801* checkInCodelets11801Ptr = (this->checkInCodelets11801);
    for (int i = 0; i < NUMTPS11801; i++) {
        TP11801Ptr[i] = nullptr;
        TP11801_alreadyLaunched[i] = 0;
    }
    _checkInCodelets11652* checkInCodelets11652Ptr = (this->checkInCodelets11652);
    for (int i = 0; i < NUMTPS11652; i++) {
        TP11652Ptr[i] = nullptr;
        TP11652_alreadyLaunched[i] = 0;
    }
    _checkInCodelets11649* checkInCodelets11649Ptr = (this->checkInCodelets11649);
    _checkInCodelets11012* checkInCodelets11012Ptr = (this->checkInCodelets11012);
    for (int i = 0; i < NUMTPS11012; i++) {
        TP11012Ptr[i] = nullptr;
        TP11012_alreadyLaunched[i] = 0;
    }
    _checkInCodelets10863* checkInCodelets10863Ptr = (this->checkInCodelets10863);
    for (int i = 0; i < NUMTPS10863; i++) {
        TP10863Ptr[i] = nullptr;
        TP10863_alreadyLaunched[i] = 0;
    }
    _checkInCodelets10860* checkInCodelets10860Ptr = (this->checkInCodelets10860);
    _checkInCodelets10803* checkInCodelets10803Ptr = (this->checkInCodelets10803);
    for (int i = 0; i < NUMTPS10803; i++) {
        TP10803Ptr[i] = nullptr;
        TP10803_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets12438Ptr) = _checkInCodelets12438(1, 1, this, codeletCounter);
        checkInCodelets12438Ptr++;
        (*checkInCodelets11801Ptr) = _checkInCodelets11801(1, 1, this, codeletCounter);
        checkInCodelets11801Ptr++;
        (*checkInCodelets11652Ptr) = _checkInCodelets11652(1, 1, this, codeletCounter);
        checkInCodelets11652Ptr++;
        (*checkInCodelets11649Ptr) = _checkInCodelets11649(1, 1, this, codeletCounter);
        checkInCodelets11649Ptr++;
        (*checkInCodelets11012Ptr) = _checkInCodelets11012(1, 1, this, codeletCounter);
        checkInCodelets11012Ptr++;
        (*checkInCodelets10863Ptr) = _checkInCodelets10863(1, 1, this, codeletCounter);
        checkInCodelets10863Ptr++;
        (*checkInCodelets10860Ptr) = _checkInCodelets10860(1, 1, this, codeletCounter);
        checkInCodelets10860Ptr++;
        (*checkInCodelets10803Ptr) = _checkInCodelets10803(1, 1, this, codeletCounter);
        (*checkInCodelets10803Ptr).decDep();
        checkInCodelets10803Ptr++;
    }
}
TP10788::~TP10788()
{
    delete[] L1_darts10788;
    delete[] L2_darts10788;
    delete[] i_darts10788;
    delete[] iend1_darts10788;
    delete[] ist1_darts10788;
    delete[] j_darts10788;
    delete[] jend1_darts10788;
    delete[] jst1_darts10788;
    delete[] k_darts10788;
    delete[] m_darts10788;
    delete[] q_darts10788;
    delete[] tmp_darts10788;
    delete[] u21_darts10788;
    delete[] u21i_darts10788;
    delete[] u21im1_darts10788;
    delete[] u21j_darts10788;
    delete[] u21jm1_darts10788;
    delete[] u21k_darts10788;
    delete[] u21km1_darts10788;
    delete[] u31_darts10788;
    delete[] u31i_darts10788;
    delete[] u31im1_darts10788;
    delete[] u31j_darts10788;
    delete[] u31jm1_darts10788;
    delete[] u31k_darts10788;
    delete[] u31km1_darts10788;
    delete[] u41_darts10788;
    delete[] u41i_darts10788;
    delete[] u41im1_darts10788;
    delete[] u41j_darts10788;
    delete[] u41jm1_darts10788;
    delete[] u41k_darts10788;
    delete[] u41km1_darts10788;
    delete[] u51i_darts10788;
    delete[] u51im1_darts10788;
    delete[] u51j_darts10788;
    delete[] u51jm1_darts10788;
    delete[] u51k_darts10788;
    delete[] u51km1_darts10788;
    delete[] barrierCodelets10788;
    delete[] barrierCodelets12438;
    delete[] checkInCodelets12438;
    delete[] barrierCodelets11801;
    delete[] checkInCodelets11801;
    delete[] barrierCodelets11652;
    delete[] checkInCodelets11652;
    delete[] checkInCodelets11649;
    delete[] barrierCodelets11012;
    delete[] checkInCodelets11012;
    delete[] barrierCodelets10863;
    delete[] checkInCodelets10863;
    delete[] checkInCodelets10860;
    delete[] barrierCodelets10803;
    delete[] checkInCodelets10803;
}
/*TP10803: OMPForDirective*/
void TP10803::_barrierCodelets10803::fire(void)
{
    TP10803* myTP = static_cast<TP10803*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets10803[0].decDep();
}
bool TP10803::requestNewRangeIterations10803(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet10803 * codeletID;
        int tempEndRange = rangePerCodelet10803 * (codeletID + 1);
        if (remainderRange10803 != 0) {
            if (codeletID < (uint32_t)remainderRange10803) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange10803;
                tempEndRange += remainderRange10803;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration10803;
        tempEndRange = tempEndRange * 1 + minIteration10803;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration10803 < lastIteration10803) {
            (this->inputsTPParent->i_darts10803[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts10803[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration10803;
        }
    }
    return isThereNewIteration;
}
void TP10803::_checkInCodelets10804::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 10804: ForStmt*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    int* i = &(this->inputsTPParent->i_darts10803[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts10803[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts10803[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts10803[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations10803(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets10803[0].decDep();
        return;
    }
    for (int i_darts_counter_temp10803 = (*i); i_darts_counter_temp10803 <= endRange
         && i_darts_counter_temp10803 <= this->inputsTPParent->lastIteration10803;
         i_darts_counter_temp10803++) {
        {
            {
                /*Loop's init*/
                (*j) = 0;
                int j_darts_counter_temp10803 = (*j);
                for (; j_darts_counter_temp10803 <= ny - 1; j_darts_counter_temp10803++) {
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp10803 = (*k);
                        for (; k_darts_counter_temp10803 <= nz - 1; k_darts_counter_temp10803++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp10803 = (*m);
                                for (; m_darts_counter_temp10803 < 5; m_darts_counter_temp10803++) {
                                    rsd[(i_darts_counter_temp10803)][j_darts_counter_temp10803]
                                       [k_darts_counter_temp10803][m_darts_counter_temp10803]
                                        = -frct[(i_darts_counter_temp10803)]
                                               [j_darts_counter_temp10803]
                                               [k_darts_counter_temp10803]
                                               [m_darts_counter_temp10803];
                                }
                                (*m) = m_darts_counter_temp10803;
                            }
                        }
                        (*k) = k_darts_counter_temp10803;
                    }
                }
                (*j) = j_darts_counter_temp10803;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets10803[0].decDep();
}
TP10803::TP10803(int in_numThreads, int in_mainCodeletID, TP10788* in_TPParent,
    int in_initIteration, int in_lastIteration, TP10803** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts10803(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts10803(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts10803(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts10803(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration10803(in_initIteration)
    , lastIteration10803(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets10803(new _barrierCodelets10803[1])
    , checkInCodelets10804(new _checkInCodelets10804[this->numThreads])
{
    /*Initialize the loop parameters*/
    range10803 = abs(lastIteration10803 - initIteration10803) / 1;
    rangePerCodelet10803 = range10803 / numThreads;
    minIteration10803 = min<int>(lastIteration10803, initIteration10803);
    remainderRange10803 = range10803 % numThreads;
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets10803[0] = _barrierCodelets10803(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets10804* checkInCodelets10804Ptr = (this->checkInCodelets10804);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets10804);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets10804Ptr) = _checkInCodelets10804(2, 1, this, codeletCounter);
#else
        (*checkInCodelets10804Ptr) = _checkInCodelets10804(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets10804Ptr).decDep();
        checkInCodelets10804Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP10803::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets10804[localID].setID(codeletID);
    this->checkInCodelets10804[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets10804[localID + this->baseNumThreads * i]
            = _checkInCodelets10804(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets10804[localID + this->baseNumThreads * i]
            = _checkInCodelets10804(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets10804[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets10804[localID + this->baseNumThreads * i].decDep();
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
TP10803::~TP10803()
{
    delete[] barrierCodelets10803;
    delete[] checkInCodelets10804;
}
/*TP10863: OMPForDirective*/
void TP10863::_barrierCodelets10863::fire(void)
{
    TP10863* myTP = static_cast<TP10863*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets10863[0].decDep();
}
bool TP10863::requestNewRangeIterations10863(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet10863 * codeletID;
        int tempEndRange = rangePerCodelet10863 * (codeletID + 1);
        if (remainderRange10863 != 0) {
            if (codeletID < (uint32_t)remainderRange10863) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange10863;
                tempEndRange += remainderRange10863;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration10863;
        tempEndRange = tempEndRange * 1 + minIteration10863;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration10863 < lastIteration10863) {
            (this->inputsTPParent->i_darts10863[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts10863[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration10863;
        }
    }
    return isThereNewIteration;
}
void TP10863::_checkInCodelets10864::fire(void)
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
    this->inputsTPParent->L1_darts10863[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L1_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L2_darts10863[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->q_darts10863[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21_darts10863[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21_darts10788[this->getID()]);

    /*printing node 10864: ForStmt*/
    /*var: L1*/
    /*var: L2*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: q*/
    /*var: u21*/
    int* i = &(this->inputsTPParent->i_darts10863[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts10863[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts10863[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts10863[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** u21 = &(this->inputsTPParent->u21_darts10863[this->getLocalID()]);
    (void)u21 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations10863(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets10863[0].decDep();
        return;
    }
    for (int i_darts_counter_temp10863 = (*i); i_darts_counter_temp10863 <= endRange
         && i_darts_counter_temp10863 <= this->inputsTPParent->lastIteration10863;
         i_darts_counter_temp10863++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp10863 = (*j);
                for (; j_darts_counter_temp10863 <= jend; j_darts_counter_temp10863++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp10863 = (*k);
                        for (; k_darts_counter_temp10863 <= nz - 2; k_darts_counter_temp10863++) {
                            flux[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                [k_darts_counter_temp10863][0]
                                = u[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                   [k_darts_counter_temp10863][1];
                            (*(*u21)) = u[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                         [k_darts_counter_temp10863][1]
                                / u[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                   [k_darts_counter_temp10863][0];
                            (*(*q)) = 0.5
                                * (u[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                    [k_darts_counter_temp10863][1]
                                        * u[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                           [k_darts_counter_temp10863][1]
                                    + u[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                       [k_darts_counter_temp10863][2]
                                        * u[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                           [k_darts_counter_temp10863][2]
                                    + u[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                       [k_darts_counter_temp10863][3]
                                        * u[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                           [k_darts_counter_temp10863][3])
                                / u[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                   [k_darts_counter_temp10863][0];
                            flux[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                [k_darts_counter_temp10863][1]
                                = u[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                   [k_darts_counter_temp10863][1]
                                    * (*(*u21))
                                + 0.40000000000000002
                                    * (u[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                        [k_darts_counter_temp10863][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                [k_darts_counter_temp10863][2]
                                = u[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                   [k_darts_counter_temp10863][2]
                                * (*(*u21));
                            flux[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                [k_darts_counter_temp10863][3]
                                = u[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                   [k_darts_counter_temp10863][3]
                                * (*(*u21));
                            flux[(i_darts_counter_temp10863)][j_darts_counter_temp10863]
                                [k_darts_counter_temp10863][4]
                                = (1.3999999999999999
                                          * u[(i_darts_counter_temp10863)]
                                             [j_darts_counter_temp10863][k_darts_counter_temp10863]
                                             [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u21));
                        }
                        (*k) = k_darts_counter_temp10863;
                    }
                }
                (*j) = j_darts_counter_temp10863;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets10863[0].decDep();
}
TP10863::TP10863(int in_numThreads, int in_mainCodeletID, TP10788* in_TPParent,
    int in_initIteration, int in_lastIteration, TP10863** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts10863(new int*[this->numThreads])
    , L2_darts10863(new int*[this->numThreads])
    , i_darts10863(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts10863(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts10863(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts10863(new double*[this->numThreads])
    , u21_darts10863(new double*[this->numThreads])
    , initIteration10863(in_initIteration)
    , lastIteration10863(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets10863(new _barrierCodelets10863[1])
    , checkInCodelets10864(new _checkInCodelets10864[this->numThreads])
{
    /*Initialize the loop parameters*/
    range10863 = abs(lastIteration10863 - initIteration10863) / 1;
    rangePerCodelet10863 = range10863 / numThreads;
    minIteration10863 = min<int>(lastIteration10863, initIteration10863);
    remainderRange10863 = range10863 % numThreads;
    /*Initialize inputs and vars.*/
    this->L1_darts10863 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->L2_darts10863 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->q_darts10863
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21_darts10863
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets10863[0] = _barrierCodelets10863(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets10864* checkInCodelets10864Ptr = (this->checkInCodelets10864);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets10864);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets10864Ptr) = _checkInCodelets10864(2, 1, this, codeletCounter);
#else
        (*checkInCodelets10864Ptr) = _checkInCodelets10864(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets10864Ptr).decDep();
        checkInCodelets10864Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP10863::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets10864[localID].setID(codeletID);
    this->checkInCodelets10864[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets10864[localID + this->baseNumThreads * i]
            = _checkInCodelets10864(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets10864[localID + this->baseNumThreads * i]
            = _checkInCodelets10864(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets10864[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets10864[localID + this->baseNumThreads * i].decDep();
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
TP10863::~TP10863()
{
    delete[] L1_darts10863;
    delete[] L2_darts10863;
    delete[] q_darts10863;
    delete[] u21_darts10863;
    delete[] barrierCodelets10863;
    delete[] checkInCodelets10864;
}
/*TP11012: OMPForDirective*/
void TP11012::_barrierCodelets11012::fire(void)
{
    TP11012* myTP = static_cast<TP11012*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets11012[0].decDep();
}
bool TP11012::requestNewRangeIterations11012(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet11012 * codeletID;
        int tempEndRange = rangePerCodelet11012 * (codeletID + 1);
        if (remainderRange11012 != 0) {
            if (codeletID < (uint32_t)remainderRange11012) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange11012;
                tempEndRange += remainderRange11012;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration11012;
        tempEndRange = tempEndRange * 1 + minIteration11012;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration11012 < lastIteration11012) {
            (this->inputsTPParent->j_darts11012[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->j_darts11012[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration11012;
        }
    }
    return isThereNewIteration;
}
void TP11012::_checkInCodelets11013::fire(void)
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
    this->inputsTPParent->L2_darts11012[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->iend1_darts11012[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iend1_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ist1_darts11012[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->ist1_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts11012[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21i_darts11012[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21i_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21im1_darts11012[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21im1_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31i_darts11012[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31i_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31im1_darts11012[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31im1_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41i_darts11012[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41i_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41im1_darts11012[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41im1_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51i_darts11012[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51i_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51im1_darts11012[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51im1_darts10788[this->getID()]);

    /*printing node 11013: ForStmt*/
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
    int** L2 = &(this->inputsTPParent->L2_darts11012[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts11012[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iend1 = &(this->inputsTPParent->iend1_darts11012[this->getLocalID()]);
    (void)iend1 /*OMP_SHARED_PRIVATE*/;
    int** ist1 = &(this->inputsTPParent->ist1_darts11012[this->getLocalID()]);
    (void)ist1 /*OMP_SHARED_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts11012[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts11012[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts11012[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts11012[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21i = &(this->inputsTPParent->u21i_darts11012[this->getLocalID()]);
    (void)u21i /*OMP_SHARED_PRIVATE*/;
    double** u21im1 = &(this->inputsTPParent->u21im1_darts11012[this->getLocalID()]);
    (void)u21im1 /*OMP_SHARED_PRIVATE*/;
    double** u31i = &(this->inputsTPParent->u31i_darts11012[this->getLocalID()]);
    (void)u31i /*OMP_SHARED_PRIVATE*/;
    double** u31im1 = &(this->inputsTPParent->u31im1_darts11012[this->getLocalID()]);
    (void)u31im1 /*OMP_SHARED_PRIVATE*/;
    double** u41i = &(this->inputsTPParent->u41i_darts11012[this->getLocalID()]);
    (void)u41i /*OMP_SHARED_PRIVATE*/;
    double** u41im1 = &(this->inputsTPParent->u41im1_darts11012[this->getLocalID()]);
    (void)u41im1 /*OMP_SHARED_PRIVATE*/;
    double** u51i = &(this->inputsTPParent->u51i_darts11012[this->getLocalID()]);
    (void)u51i /*OMP_SHARED_PRIVATE*/;
    double** u51im1 = &(this->inputsTPParent->u51im1_darts11012[this->getLocalID()]);
    (void)u51im1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations11012(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets11012[0].decDep();
        return;
    }
    for (int j_darts_counter_temp11012 = (*j); j_darts_counter_temp11012 <= endRange
         && j_darts_counter_temp11012 <= this->inputsTPParent->lastIteration11012;
         j_darts_counter_temp11012++) {
        {
            {
                /*Loop's init*/
                (*k) = 1;
                int k_darts_counter_temp11012 = (*k);
                for (; k_darts_counter_temp11012 <= nz - 2; k_darts_counter_temp11012++) {
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp11012 = (*i);
                        for (; i_darts_counter_temp11012 <= iend; i_darts_counter_temp11012++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp11012 = (*m);
                                for (; m_darts_counter_temp11012 < 5; m_darts_counter_temp11012++) {
                                    rsd[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                                       [k_darts_counter_temp11012][m_darts_counter_temp11012]
                                        = rsd[i_darts_counter_temp11012]
                                             [(j_darts_counter_temp11012)]
                                             [k_darts_counter_temp11012][m_darts_counter_temp11012]
                                        - tx2
                                            * (flux[i_darts_counter_temp11012 + 1]
                                                   [(j_darts_counter_temp11012)]
                                                   [k_darts_counter_temp11012]
                                                   [m_darts_counter_temp11012]
                                                - flux[i_darts_counter_temp11012 - 1]
                                                      [(j_darts_counter_temp11012)]
                                                      [k_darts_counter_temp11012]
                                                      [m_darts_counter_temp11012]);
                                }
                                (*m) = m_darts_counter_temp11012;
                            }
                        }
                        (*i) = i_darts_counter_temp11012;
                    }
                    (*(*L2)) = nx - 1;
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp11012 = (*i);
                        for (; i_darts_counter_temp11012 <= (*(*L2)); i_darts_counter_temp11012++) {
                            (*(*tmp)) = 1.
                                / u[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                                   [k_darts_counter_temp11012][0];
                            (*(*u21i)) = (*(*tmp))
                                * u[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                                   [k_darts_counter_temp11012][1];
                            (*(*u31i)) = (*(*tmp))
                                * u[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                                   [k_darts_counter_temp11012][2];
                            (*(*u41i)) = (*(*tmp))
                                * u[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                                   [k_darts_counter_temp11012][3];
                            (*(*u51i)) = (*(*tmp))
                                * u[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                                   [k_darts_counter_temp11012][4];
                            (*(*tmp)) = 1.
                                / u[i_darts_counter_temp11012 - 1][(j_darts_counter_temp11012)]
                                   [k_darts_counter_temp11012][0];
                            (*(*u21im1)) = (*(*tmp))
                                * u[i_darts_counter_temp11012 - 1][(j_darts_counter_temp11012)]
                                   [k_darts_counter_temp11012][1];
                            (*(*u31im1)) = (*(*tmp))
                                * u[i_darts_counter_temp11012 - 1][(j_darts_counter_temp11012)]
                                   [k_darts_counter_temp11012][2];
                            (*(*u41im1)) = (*(*tmp))
                                * u[i_darts_counter_temp11012 - 1][(j_darts_counter_temp11012)]
                                   [k_darts_counter_temp11012][3];
                            (*(*u51im1)) = (*(*tmp))
                                * u[i_darts_counter_temp11012 - 1][(j_darts_counter_temp11012)]
                                   [k_darts_counter_temp11012][4];
                            flux[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                                [k_darts_counter_temp11012][1]
                                = (4. / 3.) * tx3 * ((*(*u21i)) - (*(*u21im1)));
                            flux[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                                [k_darts_counter_temp11012][2]
                                = tx3 * ((*(*u31i)) - (*(*u31im1)));
                            flux[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                                [k_darts_counter_temp11012][3]
                                = tx3 * ((*(*u41i)) - (*(*u41im1)));
                            flux[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                                [k_darts_counter_temp11012][4]
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
                        (*i) = i_darts_counter_temp11012;
                    }
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp11012 = (*i);
                        for (; i_darts_counter_temp11012 <= iend; i_darts_counter_temp11012++) {
                            rsd[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                               [k_darts_counter_temp11012][0]
                                = rsd[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                                     [k_darts_counter_temp11012][0]
                                + dx1 * tx1
                                    * (u[i_darts_counter_temp11012 - 1][(j_darts_counter_temp11012)]
                                        [k_darts_counter_temp11012][0]
                                        - 2.
                                            * u[i_darts_counter_temp11012]
                                               [(j_darts_counter_temp11012)]
                                               [k_darts_counter_temp11012][0]
                                        + u[i_darts_counter_temp11012 + 1]
                                           [(j_darts_counter_temp11012)][k_darts_counter_temp11012]
                                           [0]);
                            rsd[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                               [k_darts_counter_temp11012][1]
                                = rsd[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                                     [k_darts_counter_temp11012][1]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp11012 + 1][(
                                           j_darts_counter_temp11012)][k_darts_counter_temp11012][1]
                                        - flux[i_darts_counter_temp11012]
                                              [(j_darts_counter_temp11012)]
                                              [k_darts_counter_temp11012][1])
                                + dx2 * tx1
                                    * (u[i_darts_counter_temp11012 - 1][(j_darts_counter_temp11012)]
                                        [k_darts_counter_temp11012][1]
                                        - 2.
                                            * u[i_darts_counter_temp11012]
                                               [(j_darts_counter_temp11012)]
                                               [k_darts_counter_temp11012][1]
                                        + u[i_darts_counter_temp11012 + 1]
                                           [(j_darts_counter_temp11012)][k_darts_counter_temp11012]
                                           [1]);
                            rsd[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                               [k_darts_counter_temp11012][2]
                                = rsd[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                                     [k_darts_counter_temp11012][2]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp11012 + 1][(
                                           j_darts_counter_temp11012)][k_darts_counter_temp11012][2]
                                        - flux[i_darts_counter_temp11012]
                                              [(j_darts_counter_temp11012)]
                                              [k_darts_counter_temp11012][2])
                                + dx3 * tx1
                                    * (u[i_darts_counter_temp11012 - 1][(j_darts_counter_temp11012)]
                                        [k_darts_counter_temp11012][2]
                                        - 2.
                                            * u[i_darts_counter_temp11012]
                                               [(j_darts_counter_temp11012)]
                                               [k_darts_counter_temp11012][2]
                                        + u[i_darts_counter_temp11012 + 1]
                                           [(j_darts_counter_temp11012)][k_darts_counter_temp11012]
                                           [2]);
                            rsd[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                               [k_darts_counter_temp11012][3]
                                = rsd[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                                     [k_darts_counter_temp11012][3]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp11012 + 1][(
                                           j_darts_counter_temp11012)][k_darts_counter_temp11012][3]
                                        - flux[i_darts_counter_temp11012]
                                              [(j_darts_counter_temp11012)]
                                              [k_darts_counter_temp11012][3])
                                + dx4 * tx1
                                    * (u[i_darts_counter_temp11012 - 1][(j_darts_counter_temp11012)]
                                        [k_darts_counter_temp11012][3]
                                        - 2.
                                            * u[i_darts_counter_temp11012]
                                               [(j_darts_counter_temp11012)]
                                               [k_darts_counter_temp11012][3]
                                        + u[i_darts_counter_temp11012 + 1]
                                           [(j_darts_counter_temp11012)][k_darts_counter_temp11012]
                                           [3]);
                            rsd[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                               [k_darts_counter_temp11012][4]
                                = rsd[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                                     [k_darts_counter_temp11012][4]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp11012 + 1][(
                                           j_darts_counter_temp11012)][k_darts_counter_temp11012][4]
                                        - flux[i_darts_counter_temp11012]
                                              [(j_darts_counter_temp11012)]
                                              [k_darts_counter_temp11012][4])
                                + dx5 * tx1
                                    * (u[i_darts_counter_temp11012 - 1][(j_darts_counter_temp11012)]
                                        [k_darts_counter_temp11012][4]
                                        - 2.
                                            * u[i_darts_counter_temp11012]
                                               [(j_darts_counter_temp11012)]
                                               [k_darts_counter_temp11012][4]
                                        + u[i_darts_counter_temp11012 + 1]
                                           [(j_darts_counter_temp11012)][k_darts_counter_temp11012]
                                           [4]);
                        }
                        (*i) = i_darts_counter_temp11012;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp11012 = (*m);
                        for (; m_darts_counter_temp11012 < 5; m_darts_counter_temp11012++) {
                            rsd[1][(j_darts_counter_temp11012)][k_darts_counter_temp11012]
                               [m_darts_counter_temp11012]
                                = rsd[1][(j_darts_counter_temp11012)][k_darts_counter_temp11012]
                                     [m_darts_counter_temp11012]
                                - dssp
                                    * (+5.
                                            * u[1][(j_darts_counter_temp11012)]
                                               [k_darts_counter_temp11012]
                                               [m_darts_counter_temp11012]
                                        - 4.
                                            * u[2][(j_darts_counter_temp11012)]
                                               [k_darts_counter_temp11012]
                                               [m_darts_counter_temp11012]
                                        + u[3][(j_darts_counter_temp11012)]
                                           [k_darts_counter_temp11012][m_darts_counter_temp11012]);
                            rsd[2][(j_darts_counter_temp11012)][k_darts_counter_temp11012]
                               [m_darts_counter_temp11012]
                                = rsd[2][(j_darts_counter_temp11012)][k_darts_counter_temp11012]
                                     [m_darts_counter_temp11012]
                                - dssp
                                    * (-4.
                                            * u[1][(j_darts_counter_temp11012)]
                                               [k_darts_counter_temp11012]
                                               [m_darts_counter_temp11012]
                                        + 6.
                                            * u[2][(j_darts_counter_temp11012)]
                                               [k_darts_counter_temp11012]
                                               [m_darts_counter_temp11012]
                                        - 4.
                                            * u[3][(j_darts_counter_temp11012)]
                                               [k_darts_counter_temp11012]
                                               [m_darts_counter_temp11012]
                                        + u[4][(j_darts_counter_temp11012)]
                                           [k_darts_counter_temp11012][m_darts_counter_temp11012]);
                        }
                        (*m) = m_darts_counter_temp11012;
                    }
                    (*(*ist1)) = 3;
                    (*(*iend1)) = nx - 4;
                    {
                        /*Loop's init*/
                        (*i) = (*(*ist1));
                        int i_darts_counter_temp11012 = (*i);
                        for (; i_darts_counter_temp11012 <= (*(*iend1));
                             i_darts_counter_temp11012++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp11012 = (*m);
                                for (; m_darts_counter_temp11012 < 5; m_darts_counter_temp11012++) {
                                    rsd[i_darts_counter_temp11012][(j_darts_counter_temp11012)]
                                       [k_darts_counter_temp11012][m_darts_counter_temp11012]
                                        = rsd[i_darts_counter_temp11012]
                                             [(j_darts_counter_temp11012)]
                                             [k_darts_counter_temp11012][m_darts_counter_temp11012]
                                        - dssp
                                            * (u[i_darts_counter_temp11012 - 2]
                                                [(j_darts_counter_temp11012)]
                                                [k_darts_counter_temp11012]
                                                [m_darts_counter_temp11012]
                                                - 4.
                                                    * u[i_darts_counter_temp11012 - 1]
                                                       [(j_darts_counter_temp11012)]
                                                       [k_darts_counter_temp11012]
                                                       [m_darts_counter_temp11012]
                                                + 6.
                                                    * u[i_darts_counter_temp11012]
                                                       [(j_darts_counter_temp11012)]
                                                       [k_darts_counter_temp11012]
                                                       [m_darts_counter_temp11012]
                                                - 4.
                                                    * u[i_darts_counter_temp11012 + 1]
                                                       [(j_darts_counter_temp11012)]
                                                       [k_darts_counter_temp11012]
                                                       [m_darts_counter_temp11012]
                                                + u[i_darts_counter_temp11012 + 2]
                                                   [(j_darts_counter_temp11012)]
                                                   [k_darts_counter_temp11012]
                                                   [m_darts_counter_temp11012]);
                                }
                                (*m) = m_darts_counter_temp11012;
                            }
                        }
                        (*i) = i_darts_counter_temp11012;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp11012 = (*m);
                        for (; m_darts_counter_temp11012 < 5; m_darts_counter_temp11012++) {
                            rsd[nx - 3][(j_darts_counter_temp11012)][k_darts_counter_temp11012]
                               [m_darts_counter_temp11012]
                                = rsd[nx - 3][(j_darts_counter_temp11012)]
                                     [k_darts_counter_temp11012][m_darts_counter_temp11012]
                                - dssp
                                    * (u[nx - 5][(j_darts_counter_temp11012)]
                                        [k_darts_counter_temp11012][m_darts_counter_temp11012]
                                        - 4.
                                            * u[nx - 4][(j_darts_counter_temp11012)]
                                               [k_darts_counter_temp11012]
                                               [m_darts_counter_temp11012]
                                        + 6.
                                            * u[nx - 3][(j_darts_counter_temp11012)]
                                               [k_darts_counter_temp11012]
                                               [m_darts_counter_temp11012]
                                        - 4.
                                            * u[nx - 2][(j_darts_counter_temp11012)]
                                               [k_darts_counter_temp11012]
                                               [m_darts_counter_temp11012]);
                            rsd[nx - 2][(j_darts_counter_temp11012)][k_darts_counter_temp11012]
                               [m_darts_counter_temp11012]
                                = rsd[nx - 2][(j_darts_counter_temp11012)]
                                     [k_darts_counter_temp11012][m_darts_counter_temp11012]
                                - dssp
                                    * (u[nx - 4][(j_darts_counter_temp11012)]
                                        [k_darts_counter_temp11012][m_darts_counter_temp11012]
                                        - 4.
                                            * u[nx - 3][(j_darts_counter_temp11012)]
                                               [k_darts_counter_temp11012]
                                               [m_darts_counter_temp11012]
                                        + 5.
                                            * u[nx - 2][(j_darts_counter_temp11012)]
                                               [k_darts_counter_temp11012]
                                               [m_darts_counter_temp11012]);
                        }
                        (*m) = m_darts_counter_temp11012;
                    }
                }
                (*k) = k_darts_counter_temp11012;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets11012[0].decDep();
}
TP11012::TP11012(int in_numThreads, int in_mainCodeletID, TP10788* in_TPParent,
    int in_initIteration, int in_lastIteration, TP11012** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L2_darts11012(new int*[this->numThreads])
    , i_darts11012(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iend1_darts11012(new int*[this->numThreads])
    , ist1_darts11012(new int*[this->numThreads])
    , j_darts11012(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts11012(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts11012(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts11012(new double*[this->numThreads])
    , u21i_darts11012(new double*[this->numThreads])
    , u21im1_darts11012(new double*[this->numThreads])
    , u31i_darts11012(new double*[this->numThreads])
    , u31im1_darts11012(new double*[this->numThreads])
    , u41i_darts11012(new double*[this->numThreads])
    , u41im1_darts11012(new double*[this->numThreads])
    , u51i_darts11012(new double*[this->numThreads])
    , u51im1_darts11012(new double*[this->numThreads])
    , initIteration11012(in_initIteration)
    , lastIteration11012(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets11012(new _barrierCodelets11012[1])
    , checkInCodelets11013(new _checkInCodelets11013[this->numThreads])
{
    /*Initialize the loop parameters*/
    range11012 = abs(lastIteration11012 - initIteration11012) / 1;
    rangePerCodelet11012 = range11012 / numThreads;
    minIteration11012 = min<int>(lastIteration11012, initIteration11012);
    remainderRange11012 = range11012 % numThreads;
    /*Initialize inputs and vars.*/
    this->L2_darts11012 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->iend1_darts11012 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ist1_darts11012 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts11012
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21i_darts11012
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21im1_darts11012
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31i_darts11012
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31im1_darts11012
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41i_darts11012
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41im1_darts11012
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51i_darts11012
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51im1_darts11012
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets11012[0] = _barrierCodelets11012(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets11013* checkInCodelets11013Ptr = (this->checkInCodelets11013);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets11013);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets11013Ptr) = _checkInCodelets11013(2, 1, this, codeletCounter);
#else
        (*checkInCodelets11013Ptr) = _checkInCodelets11013(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets11013Ptr).decDep();
        checkInCodelets11013Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP11012::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets11013[localID].setID(codeletID);
    this->checkInCodelets11013[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets11013[localID + this->baseNumThreads * i]
            = _checkInCodelets11013(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets11013[localID + this->baseNumThreads * i]
            = _checkInCodelets11013(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets11013[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets11013[localID + this->baseNumThreads * i].decDep();
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
TP11012::~TP11012()
{
    delete[] L2_darts11012;
    delete[] iend1_darts11012;
    delete[] ist1_darts11012;
    delete[] tmp_darts11012;
    delete[] u21i_darts11012;
    delete[] u21im1_darts11012;
    delete[] u31i_darts11012;
    delete[] u31im1_darts11012;
    delete[] u41i_darts11012;
    delete[] u41im1_darts11012;
    delete[] u51i_darts11012;
    delete[] u51im1_darts11012;
    delete[] barrierCodelets11012;
    delete[] checkInCodelets11013;
}
/*TP11652: OMPForDirective*/
void TP11652::_barrierCodelets11652::fire(void)
{
    TP11652* myTP = static_cast<TP11652*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets11652[0].decDep();
}
bool TP11652::requestNewRangeIterations11652(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet11652 * codeletID;
        int tempEndRange = rangePerCodelet11652 * (codeletID + 1);
        if (remainderRange11652 != 0) {
            if (codeletID < (uint32_t)remainderRange11652) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange11652;
                tempEndRange += remainderRange11652;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration11652;
        tempEndRange = tempEndRange * 1 + minIteration11652;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration11652 < lastIteration11652) {
            (this->inputsTPParent->i_darts11652[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts11652[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration11652;
        }
    }
    return isThereNewIteration;
}
void TP11652::_checkInCodelets11653::fire(void)
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
    this->inputsTPParent->L1_darts11652[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L1_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L2_darts11652[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->q_darts11652[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31_darts11652[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31_darts10788[this->getID()]);

    /*printing node 11653: ForStmt*/
    /*var: L1*/
    /*var: L2*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: q*/
    /*var: u31*/
    int** L1 = &(this->inputsTPParent->L1_darts11652[this->getLocalID()]);
    (void)L1 /*OMP_SHARED_PRIVATE*/;
    int** L2 = &(this->inputsTPParent->L2_darts11652[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts11652[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts11652[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts11652[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts11652[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** u31 = &(this->inputsTPParent->u31_darts11652[this->getLocalID()]);
    (void)u31 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations11652(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets11652[0].decDep();
        return;
    }
    for (int i_darts_counter_temp11652 = (*i); i_darts_counter_temp11652 <= endRange
         && i_darts_counter_temp11652 <= this->inputsTPParent->lastIteration11652;
         i_darts_counter_temp11652++) {
        {
            {
                /*Loop's init*/
                (*j) = (*(*L1));
                int j_darts_counter_temp11652 = (*j);
                for (; j_darts_counter_temp11652 <= (*(*L2)); j_darts_counter_temp11652++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp11652 = (*k);
                        for (; k_darts_counter_temp11652 <= nz - 2; k_darts_counter_temp11652++) {
                            flux[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                [k_darts_counter_temp11652][0]
                                = u[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                   [k_darts_counter_temp11652][2];
                            (*(*u31)) = u[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                         [k_darts_counter_temp11652][2]
                                / u[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                   [k_darts_counter_temp11652][0];
                            (*(*q)) = 0.5
                                * (u[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                    [k_darts_counter_temp11652][1]
                                        * u[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                           [k_darts_counter_temp11652][1]
                                    + u[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                       [k_darts_counter_temp11652][2]
                                        * u[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                           [k_darts_counter_temp11652][2]
                                    + u[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                       [k_darts_counter_temp11652][3]
                                        * u[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                           [k_darts_counter_temp11652][3])
                                / u[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                   [k_darts_counter_temp11652][0];
                            flux[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                [k_darts_counter_temp11652][1]
                                = u[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                   [k_darts_counter_temp11652][1]
                                * (*(*u31));
                            flux[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                [k_darts_counter_temp11652][2]
                                = u[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                   [k_darts_counter_temp11652][2]
                                    * (*(*u31))
                                + 0.40000000000000002
                                    * (u[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                        [k_darts_counter_temp11652][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                [k_darts_counter_temp11652][3]
                                = u[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                   [k_darts_counter_temp11652][3]
                                * (*(*u31));
                            flux[(i_darts_counter_temp11652)][j_darts_counter_temp11652]
                                [k_darts_counter_temp11652][4]
                                = (1.3999999999999999
                                          * u[(i_darts_counter_temp11652)]
                                             [j_darts_counter_temp11652][k_darts_counter_temp11652]
                                             [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u31));
                        }
                        (*k) = k_darts_counter_temp11652;
                    }
                }
                (*j) = j_darts_counter_temp11652;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets11652[0].decDep();
}
TP11652::TP11652(int in_numThreads, int in_mainCodeletID, TP10788* in_TPParent,
    int in_initIteration, int in_lastIteration, TP11652** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts11652(new int*[this->numThreads])
    , L2_darts11652(new int*[this->numThreads])
    , i_darts11652(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts11652(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts11652(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts11652(new double*[this->numThreads])
    , u31_darts11652(new double*[this->numThreads])
    , initIteration11652(in_initIteration)
    , lastIteration11652(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets11652(new _barrierCodelets11652[1])
    , checkInCodelets11653(new _checkInCodelets11653[this->numThreads])
{
    /*Initialize the loop parameters*/
    range11652 = abs(lastIteration11652 - initIteration11652) / 1;
    rangePerCodelet11652 = range11652 / numThreads;
    minIteration11652 = min<int>(lastIteration11652, initIteration11652);
    remainderRange11652 = range11652 % numThreads;
    /*Initialize inputs and vars.*/
    this->L1_darts11652 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->L2_darts11652 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->q_darts11652
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31_darts11652
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets11652[0] = _barrierCodelets11652(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets11653* checkInCodelets11653Ptr = (this->checkInCodelets11653);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets11653);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets11653Ptr) = _checkInCodelets11653(2, 1, this, codeletCounter);
#else
        (*checkInCodelets11653Ptr) = _checkInCodelets11653(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets11653Ptr).decDep();
        checkInCodelets11653Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP11652::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets11653[localID].setID(codeletID);
    this->checkInCodelets11653[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets11653[localID + this->baseNumThreads * i]
            = _checkInCodelets11653(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets11653[localID + this->baseNumThreads * i]
            = _checkInCodelets11653(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets11653[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets11653[localID + this->baseNumThreads * i].decDep();
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
TP11652::~TP11652()
{
    delete[] L1_darts11652;
    delete[] L2_darts11652;
    delete[] q_darts11652;
    delete[] u31_darts11652;
    delete[] barrierCodelets11652;
    delete[] checkInCodelets11653;
}
/*TP11801: OMPForDirective*/
void TP11801::_barrierCodelets11801::fire(void)
{
    TP11801* myTP = static_cast<TP11801*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets11801[0].decDep();
}
bool TP11801::requestNewRangeIterations11801(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet11801 * codeletID;
        int tempEndRange = rangePerCodelet11801 * (codeletID + 1);
        if (remainderRange11801 != 0) {
            if (codeletID < (uint32_t)remainderRange11801) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange11801;
                tempEndRange += remainderRange11801;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration11801;
        tempEndRange = tempEndRange * 1 + minIteration11801;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration11801 < lastIteration11801) {
            (this->inputsTPParent->i_darts11801[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts11801[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration11801;
        }
    }
    return isThereNewIteration;
}
void TP11801::_checkInCodelets11802::fire(void)
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
    this->inputsTPParent->L2_darts11801[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jend1_darts11801[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jend1_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jst1_darts11801[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jst1_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts11801[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21j_darts11801[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21j_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21jm1_darts11801[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21jm1_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31j_darts11801[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31j_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31jm1_darts11801[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31jm1_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41j_darts11801[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41j_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41jm1_darts11801[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41jm1_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51j_darts11801[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51j_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51jm1_darts11801[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51jm1_darts10788[this->getID()]);

    /*printing node 11802: ForStmt*/
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
    int** L2 = &(this->inputsTPParent->L2_darts11801[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts11801[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts11801[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jend1 = &(this->inputsTPParent->jend1_darts11801[this->getLocalID()]);
    (void)jend1 /*OMP_SHARED_PRIVATE*/;
    int** jst1 = &(this->inputsTPParent->jst1_darts11801[this->getLocalID()]);
    (void)jst1 /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts11801[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts11801[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts11801[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21j = &(this->inputsTPParent->u21j_darts11801[this->getLocalID()]);
    (void)u21j /*OMP_SHARED_PRIVATE*/;
    double** u21jm1 = &(this->inputsTPParent->u21jm1_darts11801[this->getLocalID()]);
    (void)u21jm1 /*OMP_SHARED_PRIVATE*/;
    double** u31j = &(this->inputsTPParent->u31j_darts11801[this->getLocalID()]);
    (void)u31j /*OMP_SHARED_PRIVATE*/;
    double** u31jm1 = &(this->inputsTPParent->u31jm1_darts11801[this->getLocalID()]);
    (void)u31jm1 /*OMP_SHARED_PRIVATE*/;
    double** u41j = &(this->inputsTPParent->u41j_darts11801[this->getLocalID()]);
    (void)u41j /*OMP_SHARED_PRIVATE*/;
    double** u41jm1 = &(this->inputsTPParent->u41jm1_darts11801[this->getLocalID()]);
    (void)u41jm1 /*OMP_SHARED_PRIVATE*/;
    double** u51j = &(this->inputsTPParent->u51j_darts11801[this->getLocalID()]);
    (void)u51j /*OMP_SHARED_PRIVATE*/;
    double** u51jm1 = &(this->inputsTPParent->u51jm1_darts11801[this->getLocalID()]);
    (void)u51jm1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations11801(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets11801[0].decDep();
        return;
    }
    for (int i_darts_counter_temp11801 = (*i); i_darts_counter_temp11801 <= endRange
         && i_darts_counter_temp11801 <= this->inputsTPParent->lastIteration11801;
         i_darts_counter_temp11801++) {
        {
            {
                /*Loop's init*/
                (*k) = 1;
                int k_darts_counter_temp11801 = (*k);
                for (; k_darts_counter_temp11801 <= nz - 2; k_darts_counter_temp11801++) {
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp11801 = (*j);
                        for (; j_darts_counter_temp11801 <= jend; j_darts_counter_temp11801++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp11801 = (*m);
                                for (; m_darts_counter_temp11801 < 5; m_darts_counter_temp11801++) {
                                    rsd[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                                       [k_darts_counter_temp11801][m_darts_counter_temp11801]
                                        = rsd[(i_darts_counter_temp11801)]
                                             [j_darts_counter_temp11801][k_darts_counter_temp11801]
                                             [m_darts_counter_temp11801]
                                        - ty2
                                            * (flux[(i_darts_counter_temp11801)]
                                                   [j_darts_counter_temp11801 + 1]
                                                   [k_darts_counter_temp11801]
                                                   [m_darts_counter_temp11801]
                                                - flux[(i_darts_counter_temp11801)]
                                                      [j_darts_counter_temp11801 - 1]
                                                      [k_darts_counter_temp11801]
                                                      [m_darts_counter_temp11801]);
                                }
                                (*m) = m_darts_counter_temp11801;
                            }
                        }
                        (*j) = j_darts_counter_temp11801;
                    }
                    (*(*L2)) = ny - 1;
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp11801 = (*j);
                        for (; j_darts_counter_temp11801 <= (*(*L2)); j_darts_counter_temp11801++) {
                            (*(*tmp)) = 1.
                                / u[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                                   [k_darts_counter_temp11801][0];
                            (*(*u21j)) = (*(*tmp))
                                * u[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                                   [k_darts_counter_temp11801][1];
                            (*(*u31j)) = (*(*tmp))
                                * u[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                                   [k_darts_counter_temp11801][2];
                            (*(*u41j)) = (*(*tmp))
                                * u[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                                   [k_darts_counter_temp11801][3];
                            (*(*u51j)) = (*(*tmp))
                                * u[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                                   [k_darts_counter_temp11801][4];
                            (*(*tmp)) = 1.
                                / u[(i_darts_counter_temp11801)][j_darts_counter_temp11801 - 1]
                                   [k_darts_counter_temp11801][0];
                            (*(*u21jm1)) = (*(*tmp))
                                * u[(i_darts_counter_temp11801)][j_darts_counter_temp11801 - 1]
                                   [k_darts_counter_temp11801][1];
                            (*(*u31jm1)) = (*(*tmp))
                                * u[(i_darts_counter_temp11801)][j_darts_counter_temp11801 - 1]
                                   [k_darts_counter_temp11801][2];
                            (*(*u41jm1)) = (*(*tmp))
                                * u[(i_darts_counter_temp11801)][j_darts_counter_temp11801 - 1]
                                   [k_darts_counter_temp11801][3];
                            (*(*u51jm1)) = (*(*tmp))
                                * u[(i_darts_counter_temp11801)][j_darts_counter_temp11801 - 1]
                                   [k_darts_counter_temp11801][4];
                            flux[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                                [k_darts_counter_temp11801][1]
                                = ty3 * ((*(*u21j)) - (*(*u21jm1)));
                            flux[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                                [k_darts_counter_temp11801][2]
                                = (4. / 3.) * ty3 * ((*(*u31j)) - (*(*u31jm1)));
                            flux[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                                [k_darts_counter_temp11801][3]
                                = ty3 * ((*(*u41j)) - (*(*u41jm1)));
                            flux[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                                [k_darts_counter_temp11801][4]
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
                        (*j) = j_darts_counter_temp11801;
                    }
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp11801 = (*j);
                        for (; j_darts_counter_temp11801 <= jend; j_darts_counter_temp11801++) {
                            rsd[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                               [k_darts_counter_temp11801][0]
                                = rsd[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                                     [k_darts_counter_temp11801][0]
                                + dy1 * ty1
                                    * (u[(i_darts_counter_temp11801)][j_darts_counter_temp11801 - 1]
                                        [k_darts_counter_temp11801][0]
                                        - 2.
                                            * u[(i_darts_counter_temp11801)]
                                               [j_darts_counter_temp11801]
                                               [k_darts_counter_temp11801][0]
                                        + u[(i_darts_counter_temp11801)][j_darts_counter_temp11801
                                            + 1][k_darts_counter_temp11801][0]);
                            rsd[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                               [k_darts_counter_temp11801][1]
                                = rsd[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                                     [k_darts_counter_temp11801][1]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp11801)][j_darts_counter_temp11801
                                           + 1][k_darts_counter_temp11801][1]
                                        - flux[(i_darts_counter_temp11801)]
                                              [j_darts_counter_temp11801][k_darts_counter_temp11801]
                                              [1])
                                + dy2 * ty1
                                    * (u[(i_darts_counter_temp11801)][j_darts_counter_temp11801 - 1]
                                        [k_darts_counter_temp11801][1]
                                        - 2.
                                            * u[(i_darts_counter_temp11801)]
                                               [j_darts_counter_temp11801]
                                               [k_darts_counter_temp11801][1]
                                        + u[(i_darts_counter_temp11801)][j_darts_counter_temp11801
                                            + 1][k_darts_counter_temp11801][1]);
                            rsd[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                               [k_darts_counter_temp11801][2]
                                = rsd[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                                     [k_darts_counter_temp11801][2]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp11801)][j_darts_counter_temp11801
                                           + 1][k_darts_counter_temp11801][2]
                                        - flux[(i_darts_counter_temp11801)]
                                              [j_darts_counter_temp11801][k_darts_counter_temp11801]
                                              [2])
                                + dy3 * ty1
                                    * (u[(i_darts_counter_temp11801)][j_darts_counter_temp11801 - 1]
                                        [k_darts_counter_temp11801][2]
                                        - 2.
                                            * u[(i_darts_counter_temp11801)]
                                               [j_darts_counter_temp11801]
                                               [k_darts_counter_temp11801][2]
                                        + u[(i_darts_counter_temp11801)][j_darts_counter_temp11801
                                            + 1][k_darts_counter_temp11801][2]);
                            rsd[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                               [k_darts_counter_temp11801][3]
                                = rsd[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                                     [k_darts_counter_temp11801][3]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp11801)][j_darts_counter_temp11801
                                           + 1][k_darts_counter_temp11801][3]
                                        - flux[(i_darts_counter_temp11801)]
                                              [j_darts_counter_temp11801][k_darts_counter_temp11801]
                                              [3])
                                + dy4 * ty1
                                    * (u[(i_darts_counter_temp11801)][j_darts_counter_temp11801 - 1]
                                        [k_darts_counter_temp11801][3]
                                        - 2.
                                            * u[(i_darts_counter_temp11801)]
                                               [j_darts_counter_temp11801]
                                               [k_darts_counter_temp11801][3]
                                        + u[(i_darts_counter_temp11801)][j_darts_counter_temp11801
                                            + 1][k_darts_counter_temp11801][3]);
                            rsd[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                               [k_darts_counter_temp11801][4]
                                = rsd[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                                     [k_darts_counter_temp11801][4]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp11801)][j_darts_counter_temp11801
                                           + 1][k_darts_counter_temp11801][4]
                                        - flux[(i_darts_counter_temp11801)]
                                              [j_darts_counter_temp11801][k_darts_counter_temp11801]
                                              [4])
                                + dy5 * ty1
                                    * (u[(i_darts_counter_temp11801)][j_darts_counter_temp11801 - 1]
                                        [k_darts_counter_temp11801][4]
                                        - 2.
                                            * u[(i_darts_counter_temp11801)]
                                               [j_darts_counter_temp11801]
                                               [k_darts_counter_temp11801][4]
                                        + u[(i_darts_counter_temp11801)][j_darts_counter_temp11801
                                            + 1][k_darts_counter_temp11801][4]);
                        }
                        (*j) = j_darts_counter_temp11801;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp11801 = (*m);
                        for (; m_darts_counter_temp11801 < 5; m_darts_counter_temp11801++) {
                            rsd[(i_darts_counter_temp11801)][1][k_darts_counter_temp11801]
                               [m_darts_counter_temp11801]
                                = rsd[(i_darts_counter_temp11801)][1][k_darts_counter_temp11801]
                                     [m_darts_counter_temp11801]
                                - dssp
                                    * (+5.
                                            * u[(i_darts_counter_temp11801)][1]
                                               [k_darts_counter_temp11801]
                                               [m_darts_counter_temp11801]
                                        - 4.
                                            * u[(i_darts_counter_temp11801)][2]
                                               [k_darts_counter_temp11801]
                                               [m_darts_counter_temp11801]
                                        + u[(i_darts_counter_temp11801)][3]
                                           [k_darts_counter_temp11801][m_darts_counter_temp11801]);
                            rsd[(i_darts_counter_temp11801)][2][k_darts_counter_temp11801]
                               [m_darts_counter_temp11801]
                                = rsd[(i_darts_counter_temp11801)][2][k_darts_counter_temp11801]
                                     [m_darts_counter_temp11801]
                                - dssp
                                    * (-4.
                                            * u[(i_darts_counter_temp11801)][1]
                                               [k_darts_counter_temp11801]
                                               [m_darts_counter_temp11801]
                                        + 6.
                                            * u[(i_darts_counter_temp11801)][2]
                                               [k_darts_counter_temp11801]
                                               [m_darts_counter_temp11801]
                                        - 4.
                                            * u[(i_darts_counter_temp11801)][3]
                                               [k_darts_counter_temp11801]
                                               [m_darts_counter_temp11801]
                                        + u[(i_darts_counter_temp11801)][4]
                                           [k_darts_counter_temp11801][m_darts_counter_temp11801]);
                        }
                        (*m) = m_darts_counter_temp11801;
                    }
                    (*(*jst1)) = 3;
                    (*(*jend1)) = ny - 4;
                    {
                        /*Loop's init*/
                        (*j) = (*(*jst1));
                        int j_darts_counter_temp11801 = (*j);
                        for (; j_darts_counter_temp11801 <= (*(*jend1));
                             j_darts_counter_temp11801++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp11801 = (*m);
                                for (; m_darts_counter_temp11801 < 5; m_darts_counter_temp11801++) {
                                    rsd[(i_darts_counter_temp11801)][j_darts_counter_temp11801]
                                       [k_darts_counter_temp11801][m_darts_counter_temp11801]
                                        = rsd[(i_darts_counter_temp11801)]
                                             [j_darts_counter_temp11801][k_darts_counter_temp11801]
                                             [m_darts_counter_temp11801]
                                        - dssp
                                            * (u[(i_darts_counter_temp11801)]
                                                [j_darts_counter_temp11801 - 2]
                                                [k_darts_counter_temp11801]
                                                [m_darts_counter_temp11801]
                                                - 4.
                                                    * u[(i_darts_counter_temp11801)]
                                                       [j_darts_counter_temp11801 - 1]
                                                       [k_darts_counter_temp11801]
                                                       [m_darts_counter_temp11801]
                                                + 6.
                                                    * u[(i_darts_counter_temp11801)]
                                                       [j_darts_counter_temp11801]
                                                       [k_darts_counter_temp11801]
                                                       [m_darts_counter_temp11801]
                                                - 4.
                                                    * u[(i_darts_counter_temp11801)]
                                                       [j_darts_counter_temp11801 + 1]
                                                       [k_darts_counter_temp11801]
                                                       [m_darts_counter_temp11801]
                                                + u[(i_darts_counter_temp11801)]
                                                   [j_darts_counter_temp11801 + 2]
                                                   [k_darts_counter_temp11801]
                                                   [m_darts_counter_temp11801]);
                                }
                                (*m) = m_darts_counter_temp11801;
                            }
                        }
                        (*j) = j_darts_counter_temp11801;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp11801 = (*m);
                        for (; m_darts_counter_temp11801 < 5; m_darts_counter_temp11801++) {
                            rsd[(i_darts_counter_temp11801)][ny - 3][k_darts_counter_temp11801]
                               [m_darts_counter_temp11801]
                                = rsd[(i_darts_counter_temp11801)][ny - 3]
                                     [k_darts_counter_temp11801][m_darts_counter_temp11801]
                                - dssp
                                    * (u[(i_darts_counter_temp11801)][ny - 5]
                                        [k_darts_counter_temp11801][m_darts_counter_temp11801]
                                        - 4.
                                            * u[(i_darts_counter_temp11801)][ny - 4]
                                               [k_darts_counter_temp11801]
                                               [m_darts_counter_temp11801]
                                        + 6.
                                            * u[(i_darts_counter_temp11801)][ny - 3]
                                               [k_darts_counter_temp11801]
                                               [m_darts_counter_temp11801]
                                        - 4.
                                            * u[(i_darts_counter_temp11801)][ny - 2]
                                               [k_darts_counter_temp11801]
                                               [m_darts_counter_temp11801]);
                            rsd[(i_darts_counter_temp11801)][ny - 2][k_darts_counter_temp11801]
                               [m_darts_counter_temp11801]
                                = rsd[(i_darts_counter_temp11801)][ny - 2]
                                     [k_darts_counter_temp11801][m_darts_counter_temp11801]
                                - dssp
                                    * (u[(i_darts_counter_temp11801)][ny - 4]
                                        [k_darts_counter_temp11801][m_darts_counter_temp11801]
                                        - 4.
                                            * u[(i_darts_counter_temp11801)][ny - 3]
                                               [k_darts_counter_temp11801]
                                               [m_darts_counter_temp11801]
                                        + 5.
                                            * u[(i_darts_counter_temp11801)][ny - 2]
                                               [k_darts_counter_temp11801]
                                               [m_darts_counter_temp11801]);
                        }
                        (*m) = m_darts_counter_temp11801;
                    }
                }
                (*k) = k_darts_counter_temp11801;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets11801[0].decDep();
}
TP11801::TP11801(int in_numThreads, int in_mainCodeletID, TP10788* in_TPParent,
    int in_initIteration, int in_lastIteration, TP11801** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L2_darts11801(new int*[this->numThreads])
    , i_darts11801(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts11801(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jend1_darts11801(new int*[this->numThreads])
    , jst1_darts11801(new int*[this->numThreads])
    , k_darts11801(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts11801(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts11801(new double*[this->numThreads])
    , u21j_darts11801(new double*[this->numThreads])
    , u21jm1_darts11801(new double*[this->numThreads])
    , u31j_darts11801(new double*[this->numThreads])
    , u31jm1_darts11801(new double*[this->numThreads])
    , u41j_darts11801(new double*[this->numThreads])
    , u41jm1_darts11801(new double*[this->numThreads])
    , u51j_darts11801(new double*[this->numThreads])
    , u51jm1_darts11801(new double*[this->numThreads])
    , initIteration11801(in_initIteration)
    , lastIteration11801(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets11801(new _barrierCodelets11801[1])
    , checkInCodelets11802(new _checkInCodelets11802[this->numThreads])
{
    /*Initialize the loop parameters*/
    range11801 = abs(lastIteration11801 - initIteration11801) / 1;
    rangePerCodelet11801 = range11801 / numThreads;
    minIteration11801 = min<int>(lastIteration11801, initIteration11801);
    remainderRange11801 = range11801 % numThreads;
    /*Initialize inputs and vars.*/
    this->L2_darts11801 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jend1_darts11801 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jst1_darts11801 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts11801
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21j_darts11801
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21jm1_darts11801
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31j_darts11801
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31jm1_darts11801
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41j_darts11801
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41jm1_darts11801
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51j_darts11801
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51jm1_darts11801
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets11801[0] = _barrierCodelets11801(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets11802* checkInCodelets11802Ptr = (this->checkInCodelets11802);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets11802);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets11802Ptr) = _checkInCodelets11802(2, 1, this, codeletCounter);
#else
        (*checkInCodelets11802Ptr) = _checkInCodelets11802(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets11802Ptr).decDep();
        checkInCodelets11802Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP11801::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets11802[localID].setID(codeletID);
    this->checkInCodelets11802[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets11802[localID + this->baseNumThreads * i]
            = _checkInCodelets11802(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets11802[localID + this->baseNumThreads * i]
            = _checkInCodelets11802(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets11802[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets11802[localID + this->baseNumThreads * i].decDep();
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
TP11801::~TP11801()
{
    delete[] L2_darts11801;
    delete[] jend1_darts11801;
    delete[] jst1_darts11801;
    delete[] tmp_darts11801;
    delete[] u21j_darts11801;
    delete[] u21jm1_darts11801;
    delete[] u31j_darts11801;
    delete[] u31jm1_darts11801;
    delete[] u41j_darts11801;
    delete[] u41jm1_darts11801;
    delete[] u51j_darts11801;
    delete[] u51jm1_darts11801;
    delete[] barrierCodelets11801;
    delete[] checkInCodelets11802;
}
/*TP12438: OMPForDirective*/
void TP12438::_barrierCodelets12438::fire(void)
{
    TP12438* myTP = static_cast<TP12438*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets12438[0].decDep();
}
bool TP12438::requestNewRangeIterations12438(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet12438 * codeletID;
        int tempEndRange = rangePerCodelet12438 * (codeletID + 1);
        if (remainderRange12438 != 0) {
            if (codeletID < (uint32_t)remainderRange12438) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange12438;
                tempEndRange += remainderRange12438;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration12438;
        tempEndRange = tempEndRange * 1 + minIteration12438;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration12438 < lastIteration12438) {
            (this->inputsTPParent->i_darts12438[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts12438[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration12438;
        }
    }
    return isThereNewIteration;
}
void TP12438::_checkInCodelets12439::fire(void)
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
    this->inputsTPParent->q_darts12438[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts12438[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21k_darts12438[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21k_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21km1_darts12438[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21km1_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31k_darts12438[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31k_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31km1_darts12438[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31km1_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41_darts12438[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41k_darts12438[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41k_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41km1_darts12438[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41km1_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51k_darts12438[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51k_darts10788[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51km1_darts12438[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51km1_darts10788[this->getID()]);

    /*printing node 12439: ForStmt*/
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
    int* i = &(this->inputsTPParent->i_darts12438[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts12438[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts12438[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts12438[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts12438[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts12438[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21k = &(this->inputsTPParent->u21k_darts12438[this->getLocalID()]);
    (void)u21k /*OMP_SHARED_PRIVATE*/;
    double** u21km1 = &(this->inputsTPParent->u21km1_darts12438[this->getLocalID()]);
    (void)u21km1 /*OMP_SHARED_PRIVATE*/;
    double** u31k = &(this->inputsTPParent->u31k_darts12438[this->getLocalID()]);
    (void)u31k /*OMP_SHARED_PRIVATE*/;
    double** u31km1 = &(this->inputsTPParent->u31km1_darts12438[this->getLocalID()]);
    (void)u31km1 /*OMP_SHARED_PRIVATE*/;
    double** u41 = &(this->inputsTPParent->u41_darts12438[this->getLocalID()]);
    (void)u41 /*OMP_SHARED_PRIVATE*/;
    double** u41k = &(this->inputsTPParent->u41k_darts12438[this->getLocalID()]);
    (void)u41k /*OMP_SHARED_PRIVATE*/;
    double** u41km1 = &(this->inputsTPParent->u41km1_darts12438[this->getLocalID()]);
    (void)u41km1 /*OMP_SHARED_PRIVATE*/;
    double** u51k = &(this->inputsTPParent->u51k_darts12438[this->getLocalID()]);
    (void)u51k /*OMP_SHARED_PRIVATE*/;
    double** u51km1 = &(this->inputsTPParent->u51km1_darts12438[this->getLocalID()]);
    (void)u51km1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations12438(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets12438[0].decDep();
        return;
    }
    for (int i_darts_counter_temp12438 = (*i); i_darts_counter_temp12438 <= endRange
         && i_darts_counter_temp12438 <= this->inputsTPParent->lastIteration12438;
         i_darts_counter_temp12438++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp12438 = (*j);
                for (; j_darts_counter_temp12438 <= jend; j_darts_counter_temp12438++) {
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp12438 = (*k);
                        for (; k_darts_counter_temp12438 <= nz - 1; k_darts_counter_temp12438++) {
                            flux[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                [k_darts_counter_temp12438][0]
                                = u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                   [k_darts_counter_temp12438][3];
                            (*(*u41)) = u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                         [k_darts_counter_temp12438][3]
                                / u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                   [k_darts_counter_temp12438][0];
                            (*(*q)) = 0.5
                                * (u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                    [k_darts_counter_temp12438][1]
                                        * u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                           [k_darts_counter_temp12438][1]
                                    + u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                       [k_darts_counter_temp12438][2]
                                        * u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                           [k_darts_counter_temp12438][2]
                                    + u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                       [k_darts_counter_temp12438][3]
                                        * u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                           [k_darts_counter_temp12438][3])
                                / u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                   [k_darts_counter_temp12438][0];
                            flux[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                [k_darts_counter_temp12438][1]
                                = u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                   [k_darts_counter_temp12438][1]
                                * (*(*u41));
                            flux[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                [k_darts_counter_temp12438][2]
                                = u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                   [k_darts_counter_temp12438][2]
                                * (*(*u41));
                            flux[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                [k_darts_counter_temp12438][3]
                                = u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                   [k_darts_counter_temp12438][3]
                                    * (*(*u41))
                                + 0.40000000000000002
                                    * (u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                        [k_darts_counter_temp12438][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                [k_darts_counter_temp12438][4]
                                = (1.3999999999999999
                                          * u[(i_darts_counter_temp12438)]
                                             [j_darts_counter_temp12438][k_darts_counter_temp12438]
                                             [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u41));
                        }
                        (*k) = k_darts_counter_temp12438;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp12438 = (*k);
                        for (; k_darts_counter_temp12438 <= nz - 2; k_darts_counter_temp12438++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp12438 = (*m);
                                for (; m_darts_counter_temp12438 < 5; m_darts_counter_temp12438++) {
                                    rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                       [k_darts_counter_temp12438][m_darts_counter_temp12438]
                                        = rsd[(i_darts_counter_temp12438)]
                                             [j_darts_counter_temp12438][k_darts_counter_temp12438]
                                             [m_darts_counter_temp12438]
                                        - tz2
                                            * (flux[(i_darts_counter_temp12438)]
                                                   [j_darts_counter_temp12438]
                                                   [k_darts_counter_temp12438 + 1]
                                                   [m_darts_counter_temp12438]
                                                - flux[k_darts_counter_temp12438 - 1]
                                                      [(i_darts_counter_temp12438)]
                                                      [j_darts_counter_temp12438]
                                                      [m_darts_counter_temp12438]);
                                }
                                (*m) = m_darts_counter_temp12438;
                            }
                        }
                        (*k) = k_darts_counter_temp12438;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp12438 = (*k);
                        for (; k_darts_counter_temp12438 <= nz - 1; k_darts_counter_temp12438++) {
                            (*(*tmp)) = 1.
                                / u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                   [k_darts_counter_temp12438][0];
                            (*(*u21k)) = (*(*tmp))
                                * u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                   [k_darts_counter_temp12438][1];
                            (*(*u31k)) = (*(*tmp))
                                * u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                   [k_darts_counter_temp12438][2];
                            (*(*u41k)) = (*(*tmp))
                                * u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                   [k_darts_counter_temp12438][3];
                            (*(*u51k)) = (*(*tmp))
                                * u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                   [k_darts_counter_temp12438][4];
                            (*(*tmp)) = 1.
                                / u[k_darts_counter_temp12438 - 1][(i_darts_counter_temp12438)]
                                   [j_darts_counter_temp12438][0];
                            (*(*u21km1)) = (*(*tmp))
                                * u[k_darts_counter_temp12438 - 1][(i_darts_counter_temp12438)]
                                   [j_darts_counter_temp12438][1];
                            (*(*u31km1)) = (*(*tmp))
                                * u[k_darts_counter_temp12438 - 1][(i_darts_counter_temp12438)]
                                   [j_darts_counter_temp12438][2];
                            (*(*u41km1)) = (*(*tmp))
                                * u[k_darts_counter_temp12438 - 1][(i_darts_counter_temp12438)]
                                   [j_darts_counter_temp12438][3];
                            (*(*u51km1)) = (*(*tmp))
                                * u[k_darts_counter_temp12438 - 1][(i_darts_counter_temp12438)]
                                   [j_darts_counter_temp12438][4];
                            flux[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                [k_darts_counter_temp12438][1]
                                = tz3 * ((*(*u21k)) - (*(*u21km1)));
                            flux[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                [k_darts_counter_temp12438][2]
                                = tz3 * ((*(*u31k)) - (*(*u31km1)));
                            flux[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                [k_darts_counter_temp12438][3]
                                = (4. / 3.) * tz3 * ((*(*u41k)) - (*(*u41km1)));
                            flux[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                [k_darts_counter_temp12438][4]
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
                        (*k) = k_darts_counter_temp12438;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp12438 = (*k);
                        for (; k_darts_counter_temp12438 <= nz - 2; k_darts_counter_temp12438++) {
                            rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                               [k_darts_counter_temp12438][0]
                                = rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                     [k_darts_counter_temp12438][0]
                                + dz1 * tz1
                                    * (u[k_darts_counter_temp12438 - 1][(i_darts_counter_temp12438)]
                                        [j_darts_counter_temp12438][0]
                                        - 2.
                                            * u[(i_darts_counter_temp12438)]
                                               [j_darts_counter_temp12438]
                                               [k_darts_counter_temp12438][0]
                                        + u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                           [k_darts_counter_temp12438 + 1][0]);
                            rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                               [k_darts_counter_temp12438][1]
                                = rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                     [k_darts_counter_temp12438][1]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                           [k_darts_counter_temp12438 + 1][1]
                                        - flux[(i_darts_counter_temp12438)]
                                              [j_darts_counter_temp12438][k_darts_counter_temp12438]
                                              [1])
                                + dz2 * tz1
                                    * (u[k_darts_counter_temp12438 - 1][(i_darts_counter_temp12438)]
                                        [j_darts_counter_temp12438][1]
                                        - 2.
                                            * u[(i_darts_counter_temp12438)]
                                               [j_darts_counter_temp12438]
                                               [k_darts_counter_temp12438][1]
                                        + u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                           [k_darts_counter_temp12438 + 1][1]);
                            rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                               [k_darts_counter_temp12438][2]
                                = rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                     [k_darts_counter_temp12438][2]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                           [k_darts_counter_temp12438 + 1][2]
                                        - flux[(i_darts_counter_temp12438)]
                                              [j_darts_counter_temp12438][k_darts_counter_temp12438]
                                              [2])
                                + dz3 * tz1
                                    * (u[k_darts_counter_temp12438 - 1][(i_darts_counter_temp12438)]
                                        [j_darts_counter_temp12438][2]
                                        - 2.
                                            * u[(i_darts_counter_temp12438)]
                                               [j_darts_counter_temp12438]
                                               [k_darts_counter_temp12438][2]
                                        + u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                           [k_darts_counter_temp12438 + 1][2]);
                            rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                               [k_darts_counter_temp12438][3]
                                = rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                     [k_darts_counter_temp12438][3]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                           [k_darts_counter_temp12438 + 1][3]
                                        - flux[(i_darts_counter_temp12438)]
                                              [j_darts_counter_temp12438][k_darts_counter_temp12438]
                                              [3])
                                + dz4 * tz1
                                    * (u[k_darts_counter_temp12438 - 1][(i_darts_counter_temp12438)]
                                        [j_darts_counter_temp12438][3]
                                        - 2.
                                            * u[(i_darts_counter_temp12438)]
                                               [j_darts_counter_temp12438]
                                               [k_darts_counter_temp12438][3]
                                        + u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                           [k_darts_counter_temp12438 + 1][3]);
                            rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                               [k_darts_counter_temp12438][4]
                                = rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                     [k_darts_counter_temp12438][4]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                           [k_darts_counter_temp12438 + 1][4]
                                        - flux[(i_darts_counter_temp12438)]
                                              [j_darts_counter_temp12438][k_darts_counter_temp12438]
                                              [4])
                                + dz5 * tz1
                                    * (u[k_darts_counter_temp12438 - 1][(i_darts_counter_temp12438)]
                                        [j_darts_counter_temp12438][4]
                                        - 2.
                                            * u[(i_darts_counter_temp12438)]
                                               [j_darts_counter_temp12438]
                                               [k_darts_counter_temp12438][4]
                                        + u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                           [k_darts_counter_temp12438 + 1][4]);
                        }
                        (*k) = k_darts_counter_temp12438;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp12438 = (*m);
                        for (; m_darts_counter_temp12438 < 5; m_darts_counter_temp12438++) {
                            rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438][1]
                               [m_darts_counter_temp12438]
                                = rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438][1]
                                     [m_darts_counter_temp12438]
                                - dssp
                                    * (+5.
                                            * u[(i_darts_counter_temp12438)]
                                               [j_darts_counter_temp12438][1]
                                               [m_darts_counter_temp12438]
                                        - 4.
                                            * u[(i_darts_counter_temp12438)]
                                               [j_darts_counter_temp12438][2]
                                               [m_darts_counter_temp12438]
                                        + u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                           [3][m_darts_counter_temp12438]);
                            rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438][2]
                               [m_darts_counter_temp12438]
                                = rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438][2]
                                     [m_darts_counter_temp12438]
                                - dssp
                                    * (-4.
                                            * u[(i_darts_counter_temp12438)]
                                               [j_darts_counter_temp12438][1]
                                               [m_darts_counter_temp12438]
                                        + 6.
                                            * u[(i_darts_counter_temp12438)]
                                               [j_darts_counter_temp12438][2]
                                               [m_darts_counter_temp12438]
                                        - 4.
                                            * u[(i_darts_counter_temp12438)]
                                               [j_darts_counter_temp12438][3]
                                               [m_darts_counter_temp12438]
                                        + u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                           [4][m_darts_counter_temp12438]);
                        }
                        (*m) = m_darts_counter_temp12438;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 3;
                        int k_darts_counter_temp12438 = (*k);
                        for (; k_darts_counter_temp12438 <= nz - 4; k_darts_counter_temp12438++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp12438 = (*m);
                                for (; m_darts_counter_temp12438 < 5; m_darts_counter_temp12438++) {
                                    rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                       [k_darts_counter_temp12438][m_darts_counter_temp12438]
                                        = rsd[(i_darts_counter_temp12438)]
                                             [j_darts_counter_temp12438][k_darts_counter_temp12438]
                                             [m_darts_counter_temp12438]
                                        - dssp
                                            * (u[(i_darts_counter_temp12438)]
                                                [j_darts_counter_temp12438]
                                                [k_darts_counter_temp12438 - 2]
                                                [m_darts_counter_temp12438]
                                                - 4.
                                                    * u[k_darts_counter_temp12438 - 1]
                                                       [(i_darts_counter_temp12438)]
                                                       [j_darts_counter_temp12438]
                                                       [m_darts_counter_temp12438]
                                                + 6.
                                                    * u[(i_darts_counter_temp12438)]
                                                       [j_darts_counter_temp12438]
                                                       [k_darts_counter_temp12438]
                                                       [m_darts_counter_temp12438]
                                                - 4.
                                                    * u[(i_darts_counter_temp12438)]
                                                       [j_darts_counter_temp12438]
                                                       [k_darts_counter_temp12438 + 1]
                                                       [m_darts_counter_temp12438]
                                                + u[(i_darts_counter_temp12438)]
                                                   [j_darts_counter_temp12438]
                                                   [k_darts_counter_temp12438 + 2]
                                                   [m_darts_counter_temp12438]);
                                }
                                (*m) = m_darts_counter_temp12438;
                            }
                        }
                        (*k) = k_darts_counter_temp12438;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp12438 = (*m);
                        for (; m_darts_counter_temp12438 < 5; m_darts_counter_temp12438++) {
                            rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438][nz - 3]
                               [m_darts_counter_temp12438]
                                = rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                     [nz - 3][m_darts_counter_temp12438]
                                - dssp
                                    * (u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                        [nz - 5][m_darts_counter_temp12438]
                                        - 4.
                                            * u[(i_darts_counter_temp12438)]
                                               [j_darts_counter_temp12438][nz - 4]
                                               [m_darts_counter_temp12438]
                                        + 6.
                                            * u[(i_darts_counter_temp12438)]
                                               [j_darts_counter_temp12438][nz - 3]
                                               [m_darts_counter_temp12438]
                                        - 4.
                                            * u[(i_darts_counter_temp12438)]
                                               [j_darts_counter_temp12438][nz - 2]
                                               [m_darts_counter_temp12438]);
                            rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438][nz - 2]
                               [m_darts_counter_temp12438]
                                = rsd[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                     [nz - 2][m_darts_counter_temp12438]
                                - dssp
                                    * (u[(i_darts_counter_temp12438)][j_darts_counter_temp12438]
                                        [nz - 4][m_darts_counter_temp12438]
                                        - 4.
                                            * u[(i_darts_counter_temp12438)]
                                               [j_darts_counter_temp12438][nz - 3]
                                               [m_darts_counter_temp12438]
                                        + 5.
                                            * u[(i_darts_counter_temp12438)]
                                               [j_darts_counter_temp12438][nz - 2]
                                               [m_darts_counter_temp12438]);
                        }
                        (*m) = m_darts_counter_temp12438;
                    }
                }
                (*j) = j_darts_counter_temp12438;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets12438[0].decDep();
}
TP12438::TP12438(int in_numThreads, int in_mainCodeletID, TP10788* in_TPParent,
    int in_initIteration, int in_lastIteration, TP12438** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts12438(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts12438(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts12438(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts12438(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts12438(new double*[this->numThreads])
    , tmp_darts12438(new double*[this->numThreads])
    , u21k_darts12438(new double*[this->numThreads])
    , u21km1_darts12438(new double*[this->numThreads])
    , u31k_darts12438(new double*[this->numThreads])
    , u31km1_darts12438(new double*[this->numThreads])
    , u41_darts12438(new double*[this->numThreads])
    , u41k_darts12438(new double*[this->numThreads])
    , u41km1_darts12438(new double*[this->numThreads])
    , u51k_darts12438(new double*[this->numThreads])
    , u51km1_darts12438(new double*[this->numThreads])
    , initIteration12438(in_initIteration)
    , lastIteration12438(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets12438(new _barrierCodelets12438[1])
    , checkInCodelets12439(new _checkInCodelets12439[this->numThreads])
{
    /*Initialize the loop parameters*/
    range12438 = abs(lastIteration12438 - initIteration12438) / 1;
    rangePerCodelet12438 = range12438 / numThreads;
    minIteration12438 = min<int>(lastIteration12438, initIteration12438);
    remainderRange12438 = range12438 % numThreads;
    /*Initialize inputs and vars.*/
    this->q_darts12438
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts12438
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21k_darts12438
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21km1_darts12438
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31k_darts12438
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31km1_darts12438
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41_darts12438
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41k_darts12438
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41km1_darts12438
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51k_darts12438
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51km1_darts12438
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets12438[0] = _barrierCodelets12438(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets12439* checkInCodelets12439Ptr = (this->checkInCodelets12439);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets12439);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets12439Ptr) = _checkInCodelets12439(2, 1, this, codeletCounter);
#else
        (*checkInCodelets12439Ptr) = _checkInCodelets12439(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets12439Ptr).decDep();
        checkInCodelets12439Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP12438::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets12439[localID].setID(codeletID);
    this->checkInCodelets12439[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets12439[localID + this->baseNumThreads * i]
            = _checkInCodelets12439(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets12439[localID + this->baseNumThreads * i]
            = _checkInCodelets12439(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets12439[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets12439[localID + this->baseNumThreads * i].decDep();
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
TP12438::~TP12438()
{
    delete[] q_darts12438;
    delete[] tmp_darts12438;
    delete[] u21k_darts12438;
    delete[] u21km1_darts12438;
    delete[] u31k_darts12438;
    delete[] u31km1_darts12438;
    delete[] u41_darts12438;
    delete[] u41k_darts12438;
    delete[] u41km1_darts12438;
    delete[] u51k_darts12438;
    delete[] u51km1_darts12438;
    delete[] barrierCodelets12438;
    delete[] checkInCodelets12439;
}
/*TP13189: OMPParallelDirective*/
void TP13189::_barrierCodelets13189::fire(void)
{
    TP13189* myTP = static_cast<TP13189*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP13189::_checkInCodelets13193::fire(void)
{
    /*region 13193 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13193;
    if (idx < myTP->TPsToUse13193) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13193_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(nx - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13193;
            int minIteration = min<int>(nx, 0);
            int remainderRange = range % myTP->TPsToUse13193;
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
            if (idx == myTP->TPsToUse13193 - 1) {
                lastIteration = nx;
            }
#if USEINVOKE == 1
            invoke<TP13193>(myTP, myTP->codeletsPerTP13193 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13193Ptr[idx]));
#else
            place<TP13193>(idx, myTP, myTP->codeletsPerTP13193 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13193Ptr[idx]));
#endif
        } else {
            if (myTP->TP13193Ptr[idx] != nullptr) {
                myTP->TP13193Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13189::_barrierCodelets13193::fire(void)
{
    TP13189* myTP = static_cast<TP13189*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets13244[codeletsCounter].decDep();
        }
    }
}
void TP13189::_checkInCodelets13244::fire(void)
{
    /*region 13244 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13244;
    if (idx < myTP->TPsToUse13244) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13244_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(nx - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13244;
            int minIteration = min<int>(nx, 0);
            int remainderRange = range % myTP->TPsToUse13244;
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
            if (idx == myTP->TPsToUse13244 - 1) {
                lastIteration = nx;
            }
#if USEINVOKE == 1
            invoke<TP13244>(myTP, myTP->codeletsPerTP13244 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13244Ptr[idx]));
#else
            place<TP13244>(idx, myTP, myTP->codeletsPerTP13244 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13244Ptr[idx]));
#endif
        } else {
            if (myTP->TP13244Ptr[idx] != nullptr) {
                myTP->TP13244Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13189::_barrierCodelets13244::fire(void)
{
    TP13189* myTP = static_cast<TP13189*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets13286[codeletsCounter].decDep();
        }
    }
}
void TP13189::_checkInCodelets13286::fire(void)
{
    /*region 13286 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13286;
    if (idx < myTP->TPsToUse13286) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13286_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(nx - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13286;
            int minIteration = min<int>(nx, 0);
            int remainderRange = range % myTP->TPsToUse13286;
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
            if (idx == myTP->TPsToUse13286 - 1) {
                lastIteration = nx;
            }
#if USEINVOKE == 1
            invoke<TP13286>(myTP, myTP->codeletsPerTP13286 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13286Ptr[idx]));
#else
            place<TP13286>(idx, myTP, myTP->codeletsPerTP13286 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13286Ptr[idx]));
#endif
        } else {
            if (myTP->TP13286Ptr[idx] != nullptr) {
                myTP->TP13286Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13189::_barrierCodelets13286::fire(void)
{
    TP13189* myTP = static_cast<TP13189*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets13330[codeletsCounter].decDep();
        }
    }
}
void TP13189::_checkInCodelets13330::fire(void)
{
    /*region 13330 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13330;
    if (idx < myTP->TPsToUse13330) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13330_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(ny - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13330;
            int minIteration = min<int>(ny, 0);
            int remainderRange = range % myTP->TPsToUse13330;
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
            if (idx == myTP->TPsToUse13330 - 1) {
                lastIteration = ny;
            }
#if USEINVOKE == 1
            invoke<TP13330>(myTP, myTP->codeletsPerTP13330 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13330Ptr[idx]));
#else
            place<TP13330>(idx, myTP, myTP->codeletsPerTP13330 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13330Ptr[idx]));
#endif
        } else {
            if (myTP->TP13330Ptr[idx] != nullptr) {
                myTP->TP13330Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13189::_barrierCodelets13330::fire(void)
{
    TP13189* myTP = static_cast<TP13189*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets13372[codeletsCounter].decDep();
        }
    }
}
void TP13189::_checkInCodelets13372::fire(void)
{
    /*region 13372 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13372;
    if (idx < myTP->TPsToUse13372) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13372_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(ny - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13372;
            int minIteration = min<int>(ny, 0);
            int remainderRange = range % myTP->TPsToUse13372;
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
            if (idx == myTP->TPsToUse13372 - 1) {
                lastIteration = ny;
            }
#if USEINVOKE == 1
            invoke<TP13372>(myTP, myTP->codeletsPerTP13372 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13372Ptr[idx]));
#else
            place<TP13372>(idx, myTP, myTP->codeletsPerTP13372 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13372Ptr[idx]));
#endif
        } else {
            if (myTP->TP13372Ptr[idx] != nullptr) {
                myTP->TP13372Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13189::_barrierCodelets13372::fire(void)
{
    TP13189* myTP = static_cast<TP13189*>(myTP_);
    myTP->TPParent->barrierCodelets13189[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets13189[0]));
}
TP13189::TP13189(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13189(new int[this->numThreads]) /*VARIABLE*/
    , iglob_darts13189(new int[this->numThreads]) /*VARIABLE*/
    , j_darts13189(new int[this->numThreads]) /*VARIABLE*/
    , jglob_darts13189(new int[this->numThreads]) /*VARIABLE*/
    , k_darts13189(new int[this->numThreads]) /*VARIABLE*/
    , TP13193Ptr(new TP13193*[NUMTPS13193])
    , TP13193_alreadyLaunched(new size_t[NUMTPS13193])
    , numTPsSet13193(0)
    , numTPsReady13193(0)
    , TPsToUse13193(NUMTPS13193)
    , codeletsPerTP13193(this->numThreads / NUMTPS13193)
    , totalCodelets13193(this->TPsToUse13193 * this->codeletsPerTP13193)
    , TP13244Ptr(new TP13244*[NUMTPS13244])
    , TP13244_alreadyLaunched(new size_t[NUMTPS13244])
    , numTPsSet13244(0)
    , numTPsReady13244(0)
    , TPsToUse13244(NUMTPS13244)
    , codeletsPerTP13244(this->numThreads / NUMTPS13244)
    , totalCodelets13244(this->TPsToUse13244 * this->codeletsPerTP13244)
    , TP13286Ptr(new TP13286*[NUMTPS13286])
    , TP13286_alreadyLaunched(new size_t[NUMTPS13286])
    , numTPsSet13286(0)
    , numTPsReady13286(0)
    , TPsToUse13286(NUMTPS13286)
    , codeletsPerTP13286(this->numThreads / NUMTPS13286)
    , totalCodelets13286(this->TPsToUse13286 * this->codeletsPerTP13286)
    , TP13330Ptr(new TP13330*[NUMTPS13330])
    , TP13330_alreadyLaunched(new size_t[NUMTPS13330])
    , numTPsSet13330(0)
    , numTPsReady13330(0)
    , TPsToUse13330(NUMTPS13330)
    , codeletsPerTP13330(this->numThreads / NUMTPS13330)
    , totalCodelets13330(this->TPsToUse13330 * this->codeletsPerTP13330)
    , TP13372Ptr(new TP13372*[NUMTPS13372])
    , TP13372_alreadyLaunched(new size_t[NUMTPS13372])
    , numTPsSet13372(0)
    , numTPsReady13372(0)
    , TPsToUse13372(NUMTPS13372)
    , codeletsPerTP13372(this->numThreads / NUMTPS13372)
    , totalCodelets13372(this->TPsToUse13372 * this->codeletsPerTP13372)
    , barrierCodelets13189(new _barrierCodelets13189[1])
    , checkInCodelets13193(new _checkInCodelets13193[this->numThreads])
    , barrierCodelets13193(new _barrierCodelets13193[1])
    , checkInCodelets13244(new _checkInCodelets13244[this->numThreads])
    , barrierCodelets13244(new _barrierCodelets13244[1])
    , checkInCodelets13286(new _checkInCodelets13286[this->numThreads])
    , barrierCodelets13286(new _barrierCodelets13286[1])
    , checkInCodelets13330(new _checkInCodelets13330[this->numThreads])
    , barrierCodelets13330(new _barrierCodelets13330[1])
    , checkInCodelets13372(new _checkInCodelets13372[this->numThreads])
    , barrierCodelets13372(new _barrierCodelets13372[1])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets13189[0] = _barrierCodelets13189(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets13372[0] = _barrierCodelets13372(NUMTPS13372, NUMTPS13372, this, 0);
    barrierCodelets13330[0] = _barrierCodelets13330(NUMTPS13330, NUMTPS13330, this, 0);
    barrierCodelets13286[0] = _barrierCodelets13286(NUMTPS13286, NUMTPS13286, this, 0);
    barrierCodelets13244[0] = _barrierCodelets13244(NUMTPS13244, NUMTPS13244, this, 0);
    barrierCodelets13193[0] = _barrierCodelets13193(NUMTPS13193, NUMTPS13193, this, 0);
    _checkInCodelets13372* checkInCodelets13372Ptr = (this->checkInCodelets13372);
    for (int i = 0; i < NUMTPS13372; i++) {
        TP13372Ptr[i] = nullptr;
        TP13372_alreadyLaunched[i] = 0;
    }
    _checkInCodelets13330* checkInCodelets13330Ptr = (this->checkInCodelets13330);
    for (int i = 0; i < NUMTPS13330; i++) {
        TP13330Ptr[i] = nullptr;
        TP13330_alreadyLaunched[i] = 0;
    }
    _checkInCodelets13286* checkInCodelets13286Ptr = (this->checkInCodelets13286);
    for (int i = 0; i < NUMTPS13286; i++) {
        TP13286Ptr[i] = nullptr;
        TP13286_alreadyLaunched[i] = 0;
    }
    _checkInCodelets13244* checkInCodelets13244Ptr = (this->checkInCodelets13244);
    for (int i = 0; i < NUMTPS13244; i++) {
        TP13244Ptr[i] = nullptr;
        TP13244_alreadyLaunched[i] = 0;
    }
    _checkInCodelets13193* checkInCodelets13193Ptr = (this->checkInCodelets13193);
    for (int i = 0; i < NUMTPS13193; i++) {
        TP13193Ptr[i] = nullptr;
        TP13193_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets13372Ptr) = _checkInCodelets13372(1, 1, this, codeletCounter);
        checkInCodelets13372Ptr++;
        (*checkInCodelets13330Ptr) = _checkInCodelets13330(1, 1, this, codeletCounter);
        checkInCodelets13330Ptr++;
        (*checkInCodelets13286Ptr) = _checkInCodelets13286(1, 1, this, codeletCounter);
        checkInCodelets13286Ptr++;
        (*checkInCodelets13244Ptr) = _checkInCodelets13244(1, 1, this, codeletCounter);
        checkInCodelets13244Ptr++;
        (*checkInCodelets13193Ptr) = _checkInCodelets13193(1, 1, this, codeletCounter);
        (*checkInCodelets13193Ptr).decDep();
        checkInCodelets13193Ptr++;
    }
}
TP13189::~TP13189()
{
    delete[] i_darts13189;
    delete[] iglob_darts13189;
    delete[] j_darts13189;
    delete[] jglob_darts13189;
    delete[] k_darts13189;
    delete[] barrierCodelets13189;
    delete[] barrierCodelets13372;
    delete[] checkInCodelets13372;
    delete[] barrierCodelets13330;
    delete[] checkInCodelets13330;
    delete[] barrierCodelets13286;
    delete[] checkInCodelets13286;
    delete[] barrierCodelets13244;
    delete[] checkInCodelets13244;
    delete[] barrierCodelets13193;
    delete[] checkInCodelets13193;
}
/*TP13193: OMPForDirective*/
void TP13193::_barrierCodelets13193::fire(void)
{
    TP13193* myTP = static_cast<TP13193*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13193[0].decDep();
}
bool TP13193::requestNewRangeIterations13193(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13193 * codeletID;
        int tempEndRange = rangePerCodelet13193 * (codeletID + 1);
        if (remainderRange13193 != 0) {
            if (codeletID < (uint32_t)remainderRange13193) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13193;
                tempEndRange += remainderRange13193;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13193;
        tempEndRange = tempEndRange * 1 + minIteration13193;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13193 < lastIteration13193) {
            (this->inputsTPParent->i_darts13193[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts13193[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13193;
        }
    }
    return isThereNewIteration;
}
void TP13193::_checkInCodelets13194::fire(void)
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
    this->inputsTPParent->iglob_darts13193[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts13189[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jglob_darts13193[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts13189[this->getID()]);

    /*printing node 13194: ForStmt*/
    /*var: i*/
    /*var: iglob*/
    /*var: j*/
    /*var: jglob*/
    int* i = &(this->inputsTPParent->i_darts13193[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iglob = &(this->inputsTPParent->iglob_darts13193[this->getLocalID()]);
    (void)iglob /*OMP_SHARED_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts13193[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jglob = &(this->inputsTPParent->jglob_darts13193[this->getLocalID()]);
    (void)jglob /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13193(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13193[0].decDep();
        return;
    }
    for (int i_darts_counter_temp13193 = (*i); i_darts_counter_temp13193 < endRange
         && i_darts_counter_temp13193 < this->inputsTPParent->lastIteration13193;
         i_darts_counter_temp13193++) {
        {
            (*(*iglob)) = (i_darts_counter_temp13193);
            {
                /*Loop's init*/
                (*j) = 0;
                int j_darts_counter_temp13193 = (*j);
                for (; j_darts_counter_temp13193 < ny; j_darts_counter_temp13193++) {
                    (*(*jglob)) = j_darts_counter_temp13193;
                    exact((*(*iglob)), (*(*jglob)), 0,
                        &u[(i_darts_counter_temp13193)][j_darts_counter_temp13193][0][0]);
                    exact((*(*iglob)), (*(*jglob)), nz - 1,
                        &u[(i_darts_counter_temp13193)][j_darts_counter_temp13193][nz - 1][0]);
                }
                (*j) = j_darts_counter_temp13193;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13193[0].decDep();
}
TP13193::TP13193(int in_numThreads, int in_mainCodeletID, TP13189* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13193** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13193(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iglob_darts13193(new int*[this->numThreads])
    , j_darts13193(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jglob_darts13193(new int*[this->numThreads])
    , initIteration13193(in_initIteration)
    , lastIteration13193(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13193(new _barrierCodelets13193[1])
    , checkInCodelets13194(new _checkInCodelets13194[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13193 = abs(lastIteration13193 - initIteration13193) / 1;
    rangePerCodelet13193 = range13193 / numThreads;
    minIteration13193 = min<int>(lastIteration13193, initIteration13193);
    remainderRange13193 = range13193 % numThreads;
    /*Initialize inputs and vars.*/
    this->iglob_darts13193 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jglob_darts13193 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets13193[0] = _barrierCodelets13193(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13194* checkInCodelets13194Ptr = (this->checkInCodelets13194);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13194);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13194Ptr) = _checkInCodelets13194(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13194Ptr) = _checkInCodelets13194(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13194Ptr).decDep();
        checkInCodelets13194Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13193::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13194[localID].setID(codeletID);
    this->checkInCodelets13194[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13194[localID + this->baseNumThreads * i]
            = _checkInCodelets13194(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13194[localID + this->baseNumThreads * i]
            = _checkInCodelets13194(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13194[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13194[localID + this->baseNumThreads * i].decDep();
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
TP13193::~TP13193()
{
    delete[] iglob_darts13193;
    delete[] jglob_darts13193;
    delete[] barrierCodelets13193;
    delete[] checkInCodelets13194;
}
/*TP13244: OMPForDirective*/
void TP13244::_barrierCodelets13244::fire(void)
{
    TP13244* myTP = static_cast<TP13244*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13244[0].decDep();
}
bool TP13244::requestNewRangeIterations13244(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13244 * codeletID;
        int tempEndRange = rangePerCodelet13244 * (codeletID + 1);
        if (remainderRange13244 != 0) {
            if (codeletID < (uint32_t)remainderRange13244) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13244;
                tempEndRange += remainderRange13244;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13244;
        tempEndRange = tempEndRange * 1 + minIteration13244;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13244 < lastIteration13244) {
            (this->inputsTPParent->i_darts13244[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts13244[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13244;
        }
    }
    return isThereNewIteration;
}
void TP13244::_checkInCodelets13245::fire(void)
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
    this->inputsTPParent->iglob_darts13244[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts13189[this->getID()]);

    /*printing node 13245: ForStmt*/
    /*var: i*/
    /*var: iglob*/
    /*var: k*/
    int* i = &(this->inputsTPParent->i_darts13244[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iglob = &(this->inputsTPParent->iglob_darts13244[this->getLocalID()]);
    (void)iglob /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts13244[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13244(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13244[0].decDep();
        return;
    }
    for (int i_darts_counter_temp13244 = (*i); i_darts_counter_temp13244 < endRange
         && i_darts_counter_temp13244 < this->inputsTPParent->lastIteration13244;
         i_darts_counter_temp13244++) {
        {
            (*(*iglob)) = (i_darts_counter_temp13244);
            {
                /*Loop's init*/
                (*k) = 0;
                int k_darts_counter_temp13244 = (*k);
                for (; k_darts_counter_temp13244 < nz; k_darts_counter_temp13244++) {
                    exact((*(*iglob)), 0, k_darts_counter_temp13244,
                        &u[(i_darts_counter_temp13244)][0][k_darts_counter_temp13244][0]);
                }
                (*k) = k_darts_counter_temp13244;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13244[0].decDep();
}
TP13244::TP13244(int in_numThreads, int in_mainCodeletID, TP13189* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13244** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13244(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iglob_darts13244(new int*[this->numThreads])
    , k_darts13244(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration13244(in_initIteration)
    , lastIteration13244(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13244(new _barrierCodelets13244[1])
    , checkInCodelets13245(new _checkInCodelets13245[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13244 = abs(lastIteration13244 - initIteration13244) / 1;
    rangePerCodelet13244 = range13244 / numThreads;
    minIteration13244 = min<int>(lastIteration13244, initIteration13244);
    remainderRange13244 = range13244 % numThreads;
    /*Initialize inputs and vars.*/
    this->iglob_darts13244 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets13244[0] = _barrierCodelets13244(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13245* checkInCodelets13245Ptr = (this->checkInCodelets13245);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13245);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13245Ptr) = _checkInCodelets13245(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13245Ptr) = _checkInCodelets13245(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13245Ptr).decDep();
        checkInCodelets13245Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13244::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13245[localID].setID(codeletID);
    this->checkInCodelets13245[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13245[localID + this->baseNumThreads * i]
            = _checkInCodelets13245(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13245[localID + this->baseNumThreads * i]
            = _checkInCodelets13245(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13245[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13245[localID + this->baseNumThreads * i].decDep();
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
TP13244::~TP13244()
{
    delete[] iglob_darts13244;
    delete[] barrierCodelets13244;
    delete[] checkInCodelets13245;
}
/*TP13286: OMPForDirective*/
void TP13286::_barrierCodelets13286::fire(void)
{
    TP13286* myTP = static_cast<TP13286*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13286[0].decDep();
}
bool TP13286::requestNewRangeIterations13286(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13286 * codeletID;
        int tempEndRange = rangePerCodelet13286 * (codeletID + 1);
        if (remainderRange13286 != 0) {
            if (codeletID < (uint32_t)remainderRange13286) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13286;
                tempEndRange += remainderRange13286;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13286;
        tempEndRange = tempEndRange * 1 + minIteration13286;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13286 < lastIteration13286) {
            (this->inputsTPParent->i_darts13286[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts13286[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13286;
        }
    }
    return isThereNewIteration;
}
void TP13286::_checkInCodelets13287::fire(void)
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
    this->inputsTPParent->iglob_darts13286[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts13189[this->getID()]);

    /*printing node 13287: ForStmt*/
    /*var: i*/
    /*var: iglob*/
    /*var: k*/
    int* i = &(this->inputsTPParent->i_darts13286[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iglob = &(this->inputsTPParent->iglob_darts13286[this->getLocalID()]);
    (void)iglob /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts13286[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13286(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13286[0].decDep();
        return;
    }
    for (int i_darts_counter_temp13286 = (*i); i_darts_counter_temp13286 < endRange
         && i_darts_counter_temp13286 < this->inputsTPParent->lastIteration13286;
         i_darts_counter_temp13286++) {
        {
            (*(*iglob)) = (i_darts_counter_temp13286);
            {
                /*Loop's init*/
                (*k) = 0;
                int k_darts_counter_temp13286 = (*k);
                for (; k_darts_counter_temp13286 < nz; k_darts_counter_temp13286++) {
                    exact((*(*iglob)), ny0 - 1, k_darts_counter_temp13286,
                        &u[(i_darts_counter_temp13286)][ny - 1][k_darts_counter_temp13286][0]);
                }
                (*k) = k_darts_counter_temp13286;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13286[0].decDep();
}
TP13286::TP13286(int in_numThreads, int in_mainCodeletID, TP13189* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13286** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13286(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iglob_darts13286(new int*[this->numThreads])
    , k_darts13286(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration13286(in_initIteration)
    , lastIteration13286(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13286(new _barrierCodelets13286[1])
    , checkInCodelets13287(new _checkInCodelets13287[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13286 = abs(lastIteration13286 - initIteration13286) / 1;
    rangePerCodelet13286 = range13286 / numThreads;
    minIteration13286 = min<int>(lastIteration13286, initIteration13286);
    remainderRange13286 = range13286 % numThreads;
    /*Initialize inputs and vars.*/
    this->iglob_darts13286 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets13286[0] = _barrierCodelets13286(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13287* checkInCodelets13287Ptr = (this->checkInCodelets13287);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13287);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13287Ptr) = _checkInCodelets13287(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13287Ptr) = _checkInCodelets13287(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13287Ptr).decDep();
        checkInCodelets13287Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13286::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13287[localID].setID(codeletID);
    this->checkInCodelets13287[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13287[localID + this->baseNumThreads * i]
            = _checkInCodelets13287(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13287[localID + this->baseNumThreads * i]
            = _checkInCodelets13287(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13287[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13287[localID + this->baseNumThreads * i].decDep();
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
TP13286::~TP13286()
{
    delete[] iglob_darts13286;
    delete[] barrierCodelets13286;
    delete[] checkInCodelets13287;
}
/*TP13330: OMPForDirective*/
void TP13330::_barrierCodelets13330::fire(void)
{
    TP13330* myTP = static_cast<TP13330*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13330[0].decDep();
}
bool TP13330::requestNewRangeIterations13330(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13330 * codeletID;
        int tempEndRange = rangePerCodelet13330 * (codeletID + 1);
        if (remainderRange13330 != 0) {
            if (codeletID < (uint32_t)remainderRange13330) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13330;
                tempEndRange += remainderRange13330;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13330;
        tempEndRange = tempEndRange * 1 + minIteration13330;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13330 < lastIteration13330) {
            (this->inputsTPParent->j_darts13330[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->j_darts13330[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13330;
        }
    }
    return isThereNewIteration;
}
void TP13330::_checkInCodelets13331::fire(void)
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
    this->inputsTPParent->jglob_darts13330[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts13189[this->getID()]);

    /*printing node 13331: ForStmt*/
    /*var: j*/
    /*var: jglob*/
    /*var: k*/
    int* j = &(this->inputsTPParent->j_darts13330[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jglob = &(this->inputsTPParent->jglob_darts13330[this->getLocalID()]);
    (void)jglob /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts13330[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13330(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13330[0].decDep();
        return;
    }
    for (int j_darts_counter_temp13330 = (*j); j_darts_counter_temp13330 < endRange
         && j_darts_counter_temp13330 < this->inputsTPParent->lastIteration13330;
         j_darts_counter_temp13330++) {
        {
            (*(*jglob)) = (j_darts_counter_temp13330);
            {
                /*Loop's init*/
                (*k) = 0;
                int k_darts_counter_temp13330 = (*k);
                for (; k_darts_counter_temp13330 < nz; k_darts_counter_temp13330++) {
                    exact(0, (*(*jglob)), k_darts_counter_temp13330,
                        &u[0][(j_darts_counter_temp13330)][k_darts_counter_temp13330][0]);
                }
                (*k) = k_darts_counter_temp13330;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13330[0].decDep();
}
TP13330::TP13330(int in_numThreads, int in_mainCodeletID, TP13189* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13330** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , j_darts13330(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jglob_darts13330(new int*[this->numThreads])
    , k_darts13330(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration13330(in_initIteration)
    , lastIteration13330(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13330(new _barrierCodelets13330[1])
    , checkInCodelets13331(new _checkInCodelets13331[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13330 = abs(lastIteration13330 - initIteration13330) / 1;
    rangePerCodelet13330 = range13330 / numThreads;
    minIteration13330 = min<int>(lastIteration13330, initIteration13330);
    remainderRange13330 = range13330 % numThreads;
    /*Initialize inputs and vars.*/
    this->jglob_darts13330 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets13330[0] = _barrierCodelets13330(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13331* checkInCodelets13331Ptr = (this->checkInCodelets13331);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13331);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13331Ptr) = _checkInCodelets13331(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13331Ptr) = _checkInCodelets13331(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13331Ptr).decDep();
        checkInCodelets13331Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13330::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13331[localID].setID(codeletID);
    this->checkInCodelets13331[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13331[localID + this->baseNumThreads * i]
            = _checkInCodelets13331(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13331[localID + this->baseNumThreads * i]
            = _checkInCodelets13331(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13331[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13331[localID + this->baseNumThreads * i].decDep();
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
TP13330::~TP13330()
{
    delete[] jglob_darts13330;
    delete[] barrierCodelets13330;
    delete[] checkInCodelets13331;
}
/*TP13372: OMPForDirective*/
void TP13372::_barrierCodelets13372::fire(void)
{
    TP13372* myTP = static_cast<TP13372*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13372[0].decDep();
}
bool TP13372::requestNewRangeIterations13372(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13372 * codeletID;
        int tempEndRange = rangePerCodelet13372 * (codeletID + 1);
        if (remainderRange13372 != 0) {
            if (codeletID < (uint32_t)remainderRange13372) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13372;
                tempEndRange += remainderRange13372;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13372;
        tempEndRange = tempEndRange * 1 + minIteration13372;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13372 < lastIteration13372) {
            (this->inputsTPParent->j_darts13372[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->j_darts13372[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13372;
        }
    }
    return isThereNewIteration;
}
void TP13372::_checkInCodelets13373::fire(void)
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
    this->inputsTPParent->jglob_darts13372[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts13189[this->getID()]);

    /*printing node 13373: ForStmt*/
    /*var: j*/
    /*var: jglob*/
    /*var: k*/
    int* j = &(this->inputsTPParent->j_darts13372[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jglob = &(this->inputsTPParent->jglob_darts13372[this->getLocalID()]);
    (void)jglob /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts13372[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13372(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13372[0].decDep();
        return;
    }
    for (int j_darts_counter_temp13372 = (*j); j_darts_counter_temp13372 < endRange
         && j_darts_counter_temp13372 < this->inputsTPParent->lastIteration13372;
         j_darts_counter_temp13372++) {
        {
            (*(*jglob)) = (j_darts_counter_temp13372);
            {
                /*Loop's init*/
                (*k) = 0;
                int k_darts_counter_temp13372 = (*k);
                for (; k_darts_counter_temp13372 < nz; k_darts_counter_temp13372++) {
                    exact(nx0 - 1, (*(*jglob)), k_darts_counter_temp13372,
                        &u[nx - 1][(j_darts_counter_temp13372)][k_darts_counter_temp13372][0]);
                }
                (*k) = k_darts_counter_temp13372;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13372[0].decDep();
}
TP13372::TP13372(int in_numThreads, int in_mainCodeletID, TP13189* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13372** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , j_darts13372(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jglob_darts13372(new int*[this->numThreads])
    , k_darts13372(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration13372(in_initIteration)
    , lastIteration13372(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13372(new _barrierCodelets13372[1])
    , checkInCodelets13373(new _checkInCodelets13373[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13372 = abs(lastIteration13372 - initIteration13372) / 1;
    rangePerCodelet13372 = range13372 / numThreads;
    minIteration13372 = min<int>(lastIteration13372, initIteration13372);
    remainderRange13372 = range13372 % numThreads;
    /*Initialize inputs and vars.*/
    this->jglob_darts13372 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets13372[0] = _barrierCodelets13372(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13373* checkInCodelets13373Ptr = (this->checkInCodelets13373);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13373);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13373Ptr) = _checkInCodelets13373(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13373Ptr) = _checkInCodelets13373(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13373Ptr).decDep();
        checkInCodelets13373Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13372::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13373[localID].setID(codeletID);
    this->checkInCodelets13373[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13373[localID + this->baseNumThreads * i]
            = _checkInCodelets13373(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13373[localID + this->baseNumThreads * i]
            = _checkInCodelets13373(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13373[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13373[localID + this->baseNumThreads * i].decDep();
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
TP13372::~TP13372()
{
    delete[] jglob_darts13372;
    delete[] barrierCodelets13372;
    delete[] checkInCodelets13373;
}
/*TP13756: OMPParallelDirective*/
void TP13756::_barrierCodelets13756::fire(void)
{
    TP13756* myTP = static_cast<TP13756*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP13756::_checkInCodelets13762::fire(void)
{
    /*region 13762 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13762;
    if (idx < myTP->TPsToUse13762) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13762_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(ny - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13762;
            int minIteration = min<int>(ny, 0);
            int remainderRange = range % myTP->TPsToUse13762;
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
            if (idx == myTP->TPsToUse13762 - 1) {
                lastIteration = ny;
            }
#if USEINVOKE == 1
            invoke<TP13762>(myTP, myTP->codeletsPerTP13762 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13762Ptr[idx]));
#else
            place<TP13762>(idx, myTP, myTP->codeletsPerTP13762 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13762Ptr[idx]));
#endif
        } else {
            if (myTP->TP13762Ptr[idx] != nullptr) {
                myTP->TP13762Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13756::_barrierCodelets13762::fire(void)
{
    TP13756* myTP = static_cast<TP13756*>(myTP_);
    myTP->TPParent->barrierCodelets13756[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets13756[0]));
}
TP13756::TP13756(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , eta_darts13756(new double[this->numThreads]) /*VARIABLE*/
    , i_darts13756(new int[this->numThreads]) /*VARIABLE*/
    , iglob_darts13756(new int[this->numThreads]) /*VARIABLE*/
    , j_darts13756(new int[this->numThreads]) /*VARIABLE*/
    , jglob_darts13756(new int[this->numThreads]) /*VARIABLE*/
    , k_darts13756(new int[this->numThreads]) /*VARIABLE*/
    , m_darts13756(new int[this->numThreads]) /*VARIABLE*/
    , peta_darts13756(new double[this->numThreads]) /*VARIABLE*/
    , pxi_darts13756(new double[this->numThreads]) /*VARIABLE*/
    , pzeta_darts13756(new double[this->numThreads]) /*VARIABLE*/
    , xi_darts13756(new double[this->numThreads]) /*VARIABLE*/
    , zeta_darts13756(new double[this->numThreads]) /*VARIABLE*/
    , TP13762Ptr(new TP13762*[NUMTPS13762])
    , TP13762_alreadyLaunched(new size_t[NUMTPS13762])
    , numTPsSet13762(0)
    , numTPsReady13762(0)
    , TPsToUse13762(NUMTPS13762)
    , codeletsPerTP13762(this->numThreads / NUMTPS13762)
    , totalCodelets13762(this->TPsToUse13762 * this->codeletsPerTP13762)
    , barrierCodelets13756(new _barrierCodelets13756[1])
    , checkInCodelets13762(new _checkInCodelets13762[this->numThreads])
    , barrierCodelets13762(new _barrierCodelets13762[1])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets13756[0] = _barrierCodelets13756(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets13762[0] = _barrierCodelets13762(NUMTPS13762, NUMTPS13762, this, 0);
    _checkInCodelets13762* checkInCodelets13762Ptr = (this->checkInCodelets13762);
    for (int i = 0; i < NUMTPS13762; i++) {
        TP13762Ptr[i] = nullptr;
        TP13762_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets13762Ptr) = _checkInCodelets13762(1, 1, this, codeletCounter);
        (*checkInCodelets13762Ptr).decDep();
        checkInCodelets13762Ptr++;
    }
}
TP13756::~TP13756()
{
    delete[] eta_darts13756;
    delete[] i_darts13756;
    delete[] iglob_darts13756;
    delete[] j_darts13756;
    delete[] jglob_darts13756;
    delete[] k_darts13756;
    delete[] m_darts13756;
    delete[] peta_darts13756;
    delete[] pxi_darts13756;
    delete[] pzeta_darts13756;
    delete[] xi_darts13756;
    delete[] zeta_darts13756;
    delete[] barrierCodelets13756;
    delete[] barrierCodelets13762;
    delete[] checkInCodelets13762;
}
/*TP13762: OMPForDirective*/
void TP13762::_barrierCodelets13762::fire(void)
{
    TP13762* myTP = static_cast<TP13762*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13762[0].decDep();
}
bool TP13762::requestNewRangeIterations13762(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13762 * codeletID;
        int tempEndRange = rangePerCodelet13762 * (codeletID + 1);
        if (remainderRange13762 != 0) {
            if (codeletID < (uint32_t)remainderRange13762) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13762;
                tempEndRange += remainderRange13762;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13762;
        tempEndRange = tempEndRange * 1 + minIteration13762;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13762 < lastIteration13762) {
            (this->inputsTPParent->j_darts13762[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->j_darts13762[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13762;
        }
    }
    return isThereNewIteration;
}
void TP13762::_checkInCodelets13763::fire(void)
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
    this->inputsTPParent->eta_darts13762[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->eta_darts13756[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->iglob_darts13762[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts13756[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jglob_darts13762[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts13756[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->peta_darts13762[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->peta_darts13756[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->pxi_darts13762[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->pxi_darts13756[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->pzeta_darts13762[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->pzeta_darts13756[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->xi_darts13762[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->xi_darts13756[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->zeta_darts13762[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->zeta_darts13756[this->getID()]);

    /*printing node 13763: ForStmt*/
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
    /*var: xi*/
    /*var: zeta*/
    double** eta = &(this->inputsTPParent->eta_darts13762[this->getLocalID()]);
    (void)eta /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts13762[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iglob = &(this->inputsTPParent->iglob_darts13762[this->getLocalID()]);
    (void)iglob /*OMP_SHARED_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts13762[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jglob = &(this->inputsTPParent->jglob_darts13762[this->getLocalID()]);
    (void)jglob /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts13762[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts13762[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** peta = &(this->inputsTPParent->peta_darts13762[this->getLocalID()]);
    (void)peta /*OMP_SHARED_PRIVATE*/;
    double** pxi = &(this->inputsTPParent->pxi_darts13762[this->getLocalID()]);
    (void)pxi /*OMP_SHARED_PRIVATE*/;
    double** pzeta = &(this->inputsTPParent->pzeta_darts13762[this->getLocalID()]);
    (void)pzeta /*OMP_SHARED_PRIVATE*/;
    double** xi = &(this->inputsTPParent->xi_darts13762[this->getLocalID()]);
    (void)xi /*OMP_SHARED_PRIVATE*/;
    double** zeta = &(this->inputsTPParent->zeta_darts13762[this->getLocalID()]);
    (void)zeta /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13762(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13762[0].decDep();
        return;
    }
    for (int j_darts_counter_temp13762 = (*j); j_darts_counter_temp13762 < endRange
         && j_darts_counter_temp13762 < this->inputsTPParent->lastIteration13762;
         j_darts_counter_temp13762++) {
        {
            (*(*jglob)) = (j_darts_counter_temp13762);
            {
                /*Loop's init*/
                (*k) = 1;
                int k_darts_counter_temp13762 = (*k);
                for (; k_darts_counter_temp13762 < nz - 1; k_darts_counter_temp13762++) {
                    (*(*zeta)) = ((double)k_darts_counter_temp13762) / (nz - 1);
                    if ((*(*jglob)) != 0 && (*(*jglob)) != ny0 - 1) {
                        (*(*eta)) = ((double)((*(*jglob)))) / (ny0 - 1);
                        {
                            /*Loop's init*/
                            (*i) = 0;
                            int i_darts_counter_temp13762 = (*i);
                            for (; i_darts_counter_temp13762 < nx; i_darts_counter_temp13762++) {
                                (*(*iglob)) = i_darts_counter_temp13762;
                                if ((*(*iglob)) != 0 && (*(*iglob)) != nx0 - 1) {
                                    (*(*xi)) = ((double)((*(*iglob)))) / (nx0 - 1);
                                    exact(0, (*(*jglob)), k_darts_counter_temp13762, ue_1jk);
                                    exact(
                                        nx0 - 1, (*(*jglob)), k_darts_counter_temp13762, ue_nx0jk);
                                    exact((*(*iglob)), 0, k_darts_counter_temp13762, ue_i1k);
                                    exact(
                                        (*(*iglob)), ny0 - 1, k_darts_counter_temp13762, ue_iny0k);
                                    exact((*(*iglob)), (*(*jglob)), 0, ue_ij1);
                                    exact((*(*iglob)), (*(*jglob)), nz - 1, ue_ijnz);
                                    {
                                        /*Loop's init*/
                                        (*m) = 0;
                                        int m_darts_counter_temp13762 = (*m);
                                        for (; m_darts_counter_temp13762 < 5;
                                             m_darts_counter_temp13762++) {
                                            (*(*pxi)) = (1. - (*(*xi)))
                                                    * ue_1jk[m_darts_counter_temp13762]
                                                + (*(*xi)) * ue_nx0jk[m_darts_counter_temp13762];
                                            (*(*peta)) = (1. - (*(*eta)))
                                                    * ue_i1k[m_darts_counter_temp13762]
                                                + (*(*eta)) * ue_iny0k[m_darts_counter_temp13762];
                                            (*(*pzeta)) = (1. - (*(*zeta)))
                                                    * ue_ij1[m_darts_counter_temp13762]
                                                + (*(*zeta)) * ue_ijnz[m_darts_counter_temp13762];
                                            u[i_darts_counter_temp13762]
                                             [(j_darts_counter_temp13762)]
                                             [k_darts_counter_temp13762][m_darts_counter_temp13762]
                                                = (*(*pxi)) + (*(*peta)) + (*(*pzeta))
                                                - (*(*pxi)) * (*(*peta)) - (*(*peta)) * (*(*pzeta))
                                                - (*(*pzeta)) * (*(*pxi))
                                                + (*(*pxi)) * (*(*peta)) * (*(*pzeta));
                                        }
                                        (*m) = m_darts_counter_temp13762;
                                    }
                                }
                            }
                            (*i) = i_darts_counter_temp13762;
                        }
                    }
                }
                (*k) = k_darts_counter_temp13762;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13762[0].decDep();
}
TP13762::TP13762(int in_numThreads, int in_mainCodeletID, TP13756* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13762** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , eta_darts13762(new double*[this->numThreads])
    , i_darts13762(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iglob_darts13762(new int*[this->numThreads])
    , j_darts13762(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jglob_darts13762(new int*[this->numThreads])
    , k_darts13762(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts13762(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , peta_darts13762(new double*[this->numThreads])
    , pxi_darts13762(new double*[this->numThreads])
    , pzeta_darts13762(new double*[this->numThreads])
    , xi_darts13762(new double*[this->numThreads])
    , zeta_darts13762(new double*[this->numThreads])
    , initIteration13762(in_initIteration)
    , lastIteration13762(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13762(new _barrierCodelets13762[1])
    , checkInCodelets13763(new _checkInCodelets13763[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13762 = abs(lastIteration13762 - initIteration13762) / 1;
    rangePerCodelet13762 = range13762 / numThreads;
    minIteration13762 = min<int>(lastIteration13762, initIteration13762);
    remainderRange13762 = range13762 % numThreads;
    /*Initialize inputs and vars.*/
    this->eta_darts13762
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->iglob_darts13762 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jglob_darts13762 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->peta_darts13762
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->pxi_darts13762
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->pzeta_darts13762
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->xi_darts13762
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->zeta_darts13762
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets13762[0] = _barrierCodelets13762(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13763* checkInCodelets13763Ptr = (this->checkInCodelets13763);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13763);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13763Ptr) = _checkInCodelets13763(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13763Ptr) = _checkInCodelets13763(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13763Ptr).decDep();
        checkInCodelets13763Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13762::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13763[localID].setID(codeletID);
    this->checkInCodelets13763[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13763[localID + this->baseNumThreads * i]
            = _checkInCodelets13763(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13763[localID + this->baseNumThreads * i]
            = _checkInCodelets13763(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13763[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13763[localID + this->baseNumThreads * i].decDep();
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
TP13762::~TP13762()
{
    delete[] eta_darts13762;
    delete[] iglob_darts13762;
    delete[] jglob_darts13762;
    delete[] peta_darts13762;
    delete[] pxi_darts13762;
    delete[] pzeta_darts13762;
    delete[] xi_darts13762;
    delete[] zeta_darts13762;
    delete[] barrierCodelets13762;
    delete[] checkInCodelets13763;
}
/*TP13888: OMPParallelDirective*/
void TP13888::_barrierCodelets13888::fire(void)
{
    TP13888* myTP = static_cast<TP13888*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP13888::_checkInCodelets13890::fire(void)
{
    /*region 13890 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13890;
    if (idx < myTP->TPsToUse13890) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13890_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(12 - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13890;
            int minIteration = min<int>(12, 0);
            int remainderRange = range % myTP->TPsToUse13890;
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
            if (0 < 12) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse13890 - 1) {
                lastIteration = 12;
            }
#if USEINVOKE == 1
            invoke<TP13890>(myTP, myTP->codeletsPerTP13890 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13890Ptr[idx]));
#else
            place<TP13890>(idx, myTP, myTP->codeletsPerTP13890 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13890Ptr[idx]));
#endif
        } else {
            if (myTP->TP13890Ptr[idx] != nullptr) {
                myTP->TP13890Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13888::_barrierCodelets13890::fire(void)
{
    TP13888* myTP = static_cast<TP13888*>(myTP_);
    myTP->TPParent->barrierCodelets13888[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets13888[0]));
}
TP13888::TP13888(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13888(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts13888(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts13888(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts13888(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , TP13890Ptr(new TP13890*[NUMTPS13890])
    , TP13890_alreadyLaunched(new size_t[NUMTPS13890])
    , numTPsSet13890(0)
    , numTPsReady13890(0)
    , TPsToUse13890(NUMTPS13890)
    , codeletsPerTP13890(this->numThreads / NUMTPS13890)
    , totalCodelets13890(this->TPsToUse13890 * this->codeletsPerTP13890)
    , barrierCodelets13888(new _barrierCodelets13888[1])
    , checkInCodelets13890(new _checkInCodelets13890[this->numThreads])
    , barrierCodelets13890(new _barrierCodelets13890[1])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets13888[0] = _barrierCodelets13888(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets13890[0] = _barrierCodelets13890(NUMTPS13890, NUMTPS13890, this, 0);
    _checkInCodelets13890* checkInCodelets13890Ptr = (this->checkInCodelets13890);
    for (int i = 0; i < NUMTPS13890; i++) {
        TP13890Ptr[i] = nullptr;
        TP13890_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets13890Ptr) = _checkInCodelets13890(1, 1, this, codeletCounter);
        (*checkInCodelets13890Ptr).decDep();
        checkInCodelets13890Ptr++;
    }
}
TP13888::~TP13888()
{
    delete[] barrierCodelets13888;
    delete[] barrierCodelets13890;
    delete[] checkInCodelets13890;
}
/*TP13890: OMPForDirective*/
void TP13890::_barrierCodelets13890::fire(void)
{
    TP13890* myTP = static_cast<TP13890*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13890[0].decDep();
}
bool TP13890::requestNewRangeIterations13890(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13890 * codeletID;
        int tempEndRange = rangePerCodelet13890 * (codeletID + 1);
        if (remainderRange13890 != 0) {
            if (codeletID < (uint32_t)remainderRange13890) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13890;
                tempEndRange += remainderRange13890;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13890;
        tempEndRange = tempEndRange * 1 + minIteration13890;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13890 < lastIteration13890) {
            (this->inputsTPParent->i_darts13890[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts13890[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13890;
        }
    }
    return isThereNewIteration;
}
void TP13890::_checkInCodelets13891::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 13891: ForStmt*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    int* i = &(this->inputsTPParent->i_darts13890[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts13890[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts13890[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts13890[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13890(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13890[0].decDep();
        return;
    }
    for (int i_darts_counter_temp13890 = (*i); i_darts_counter_temp13890 < endRange
         && i_darts_counter_temp13890 < this->inputsTPParent->lastIteration13890;
         i_darts_counter_temp13890++) {
        {
            {
                /*Loop's init*/
                (*j) = 0;
                int j_darts_counter_temp13890 = (*j);
                for (; j_darts_counter_temp13890 < 12; j_darts_counter_temp13890++) {
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp13890 = (*k);
                        for (; k_darts_counter_temp13890 < 5; k_darts_counter_temp13890++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp13890 = (*m);
                                for (; m_darts_counter_temp13890 < 5; m_darts_counter_temp13890++) {
                                    a[(i_darts_counter_temp13890)][j_darts_counter_temp13890]
                                     [k_darts_counter_temp13890][m_darts_counter_temp13890]
                                        = 0.;
                                    b[(i_darts_counter_temp13890)][j_darts_counter_temp13890]
                                     [k_darts_counter_temp13890][m_darts_counter_temp13890]
                                        = 0.;
                                    c[(i_darts_counter_temp13890)][j_darts_counter_temp13890]
                                     [k_darts_counter_temp13890][m_darts_counter_temp13890]
                                        = 0.;
                                    d[(i_darts_counter_temp13890)][j_darts_counter_temp13890]
                                     [k_darts_counter_temp13890][m_darts_counter_temp13890]
                                        = 0.;
                                }
                                (*m) = m_darts_counter_temp13890;
                            }
                        }
                        (*k) = k_darts_counter_temp13890;
                    }
                }
                (*j) = j_darts_counter_temp13890;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13890[0].decDep();
}
TP13890::TP13890(int in_numThreads, int in_mainCodeletID, TP13888* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13890** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13890(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts13890(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts13890(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts13890(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration13890(in_initIteration)
    , lastIteration13890(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13890(new _barrierCodelets13890[1])
    , checkInCodelets13891(new _checkInCodelets13891[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13890 = abs(lastIteration13890 - initIteration13890) / 1;
    rangePerCodelet13890 = range13890 / numThreads;
    minIteration13890 = min<int>(lastIteration13890, initIteration13890);
    remainderRange13890 = range13890 % numThreads;
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets13890[0] = _barrierCodelets13890(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13891* checkInCodelets13891Ptr = (this->checkInCodelets13891);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13891);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13891Ptr) = _checkInCodelets13891(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13891Ptr) = _checkInCodelets13891(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13891Ptr).decDep();
        checkInCodelets13891Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13890::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13891[localID].setID(codeletID);
    this->checkInCodelets13891[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13891[localID + this->baseNumThreads * i]
            = _checkInCodelets13891(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13891[localID + this->baseNumThreads * i]
            = _checkInCodelets13891(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13891[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13891[localID + this->baseNumThreads * i].decDep();
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
TP13890::~TP13890()
{
    delete[] barrierCodelets13890;
    delete[] checkInCodelets13891;
}
/*TP13987: OMPParallelDirective*/
void TP13987::_barrierCodelets13987::fire(void)
{
    TP13987* myTP = static_cast<TP13987*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP13987::_checkInCodelets13989::fire(void)
{
    /*region 13989 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13989;
    if (idx < myTP->TPsToUse13989) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13989_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13989;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse13989;
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
            if (idx == myTP->TPsToUse13989 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse13989 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP13989>(myTP, myTP->codeletsPerTP13989 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13989Ptr[idx]));
#else
            place<TP13989>(idx, myTP, myTP->codeletsPerTP13989 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13989Ptr[idx]));
#endif
        } else {
            if (myTP->TP13989Ptr[idx] != nullptr) {
                myTP->TP13989Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13987::_barrierCodelets13989::fire(void)
{
    TP13987* myTP = static_cast<TP13987*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets14045[codeletsCounter].decDep();
        }
    }
}
void TP13987::_checkInCodelets14045::fire(void)
{

    /*printing node 14045: BinaryOperator*/
    (this->inputsTPParent->k_darts13987[this->getID()]) = 1;

    /*printing node 14046: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts13987[this->getID()]) <= nz - 2) {
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets14044[this->getID()].decDep();
        return;
    } else {
        /*Signal the codelet after the loop from the end condional node.*/
        /*Signaling next codelet region: 14048 nextRegion: 14052 */
        myTP->controlTPParent->barrierCodelets14052[0].decDep();
        return;
    }
}
void TP13987::_checkInCodelets14044::fire(void)
{

    /*printing node 14044: ForStmt*/
    bool haveToLaunch = __sync_bool_compare_and_swap(&(myTP->controlTPParent->TP14044_LoopCounter),
        myTP->controlTPParent->TP14044_LoopCounterPerThread[this->getID()],
        myTP->controlTPParent->TP14044_LoopCounterPerThread[this->getID()] + 1);
    unsigned int iterIdx = myTP->controlTPParent->TP14044_LoopCounterPerThread[this->getID()];
    if (haveToLaunch) {
        this->resetCodelet();
        myTP->controlTPParent->TP14044PtrVec.push_back(nullptr);
        myTP->controlTPParent->TP14044_LoopCounterPerThread[this->getID()] += 1;
        invoke<TP14044>(myTP, myTP->numThreads, this->getID(), myTP, myTP->inputsTPParent,
            &(myTP->controlTPParent->TP14044PtrVec.back()));
    } else {
        if (myTP->controlTPParent->TP14044PtrVec.size() == 0) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP14044PtrVec.size() < (iterIdx + 1)) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP14044PtrVec[iterIdx] == nullptr) {
            this->resetCodelet();
            this->decDep();
            return;
        } else {
            this->resetCodelet();
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14044PtrVec[iterIdx]->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14044PtrVec[iterIdx]->availableCodelets[this->getID()] = 1;
#endif
            myTP->controlTPParent->TP14044_LoopCounterPerThread[this->getID()] += 1;
        }
    }
}
void TP13987::_checkInCodelets14048::fire(void)
{

    /*printing node 14048: UnaryOperator*/
    (this->inputsTPParent->k_darts13987[this->getID()])++;

    /*printing node 14503: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts13987[this->getID()]) <= nz - 2) {
        this->resetCodelet();
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets14044[this->getID()].decDep();
        return;
    } else {
        /*Signal the codelet after the loop from the condtional node.*/
        /*Signaling next codelet region: 14048 nextRegion: 14052 */
        myTP->controlTPParent->barrierCodelets14052[0].decDep();
        return;
    }
}
void TP13987::_barrierCodelets14052::fire(void)
{
    TP13987* myTP = static_cast<TP13987*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets14054[codeletsCounter].decDep();
        }
    }
}
void TP13987::_checkInCodelets14054::fire(void)
{

    /*printing node 14054: BinaryOperator*/
    (this->inputsTPParent->k_darts13987[this->getID()]) = nz - 2;

    /*printing node 14056: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts13987[this->getID()]) >= 1) {
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets14053[this->getID()].decDep();
        return;
    } else {
        /*Signal the codelet after the loop from the end condional node.*/
        /*Signaling next codelet region: 14057 nextRegion: 14061 */
        myTP->controlTPParent->barrierCodelets14061[0].decDep();
        return;
    }
}
void TP13987::_checkInCodelets14053::fire(void)
{

    /*printing node 14053: ForStmt*/
    bool haveToLaunch = __sync_bool_compare_and_swap(&(myTP->controlTPParent->TP14053_LoopCounter),
        myTP->controlTPParent->TP14053_LoopCounterPerThread[this->getID()],
        myTP->controlTPParent->TP14053_LoopCounterPerThread[this->getID()] + 1);
    unsigned int iterIdx = myTP->controlTPParent->TP14053_LoopCounterPerThread[this->getID()];
    if (haveToLaunch) {
        this->resetCodelet();
        myTP->controlTPParent->TP14053PtrVec.push_back(nullptr);
        myTP->controlTPParent->TP14053_LoopCounterPerThread[this->getID()] += 1;
        invoke<TP14053>(myTP, myTP->numThreads, this->getID(), myTP, myTP->inputsTPParent,
            &(myTP->controlTPParent->TP14053PtrVec.back()));
    } else {
        if (myTP->controlTPParent->TP14053PtrVec.size() == 0) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP14053PtrVec.size() < (iterIdx + 1)) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP14053PtrVec[iterIdx] == nullptr) {
            this->resetCodelet();
            this->decDep();
            return;
        } else {
            this->resetCodelet();
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14053PtrVec[iterIdx]->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14053PtrVec[iterIdx]->availableCodelets[this->getID()] = 1;
#endif
            myTP->controlTPParent->TP14053_LoopCounterPerThread[this->getID()] += 1;
        }
    }
}
void TP13987::_checkInCodelets14057::fire(void)
{

    /*printing node 14057: UnaryOperator*/
    (this->inputsTPParent->k_darts13987[this->getID()])--;

    /*printing node 14504: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts13987[this->getID()]) >= 1) {
        this->resetCodelet();
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets14053[this->getID()].decDep();
        return;
    } else {
        /*Signal the codelet after the loop from the condtional node.*/
        /*Signaling next codelet region: 14057 nextRegion: 14061 */
        myTP->controlTPParent->barrierCodelets14061[0].decDep();
        return;
    }
}
void TP13987::_barrierCodelets14061::fire(void)
{
    TP13987* myTP = static_cast<TP13987*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets14062[codeletsCounter].decDep();
        }
    }
}
void TP13987::_checkInCodelets14062::fire(void)
{
    /*region 14062 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP14062;
    if (idx < myTP->TPsToUse14062) {
        if (!__sync_val_compare_and_swap(&(myTP->TP14062_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse14062;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse14062;
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
            if (idx == myTP->TPsToUse14062 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse14062 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP14062>(myTP, myTP->codeletsPerTP14062 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(*(this->inputsTPParent->tmp_darts13987)),
                &(myTP->TP14062Ptr[idx]));
#else
            place<TP14062>(idx, myTP, myTP->codeletsPerTP14062 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(*(this->inputsTPParent->tmp_darts13987)),
                &(myTP->TP14062Ptr[idx]));
#endif
        } else {
            if (myTP->TP14062Ptr[idx] != nullptr) {
                myTP->TP14062Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13987::_barrierCodelets14062::fire(void)
{
    TP13987* myTP = static_cast<TP13987*>(myTP_);
    myTP->TPParent->barrierCodelets13987[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets13987[0]));
}
TP13987::TP13987(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet, double* in_tmp)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13987(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , istep_darts13987(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts13987(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts13987(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts13987(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts13987(in_tmp) /*OMP_SHARED - INPUT*/
    , TP13989Ptr(new TP13989*[NUMTPS13989])
    , TP13989_alreadyLaunched(new size_t[NUMTPS13989])
    , numTPsSet13989(0)
    , numTPsReady13989(0)
    , TPsToUse13989(NUMTPS13989)
    , codeletsPerTP13989(this->numThreads / NUMTPS13989)
    , totalCodelets13989(this->TPsToUse13989 * this->codeletsPerTP13989)
    , TP14044_LoopCounter(0)
    , TP14044_LoopCounterPerThread(new unsigned int[this->numThreads])
    , TP14053_LoopCounter(0)
    , TP14053_LoopCounterPerThread(new unsigned int[this->numThreads])
    , TP14062Ptr(new TP14062*[NUMTPS14062])
    , TP14062_alreadyLaunched(new size_t[NUMTPS14062])
    , numTPsSet14062(0)
    , numTPsReady14062(0)
    , TPsToUse14062(NUMTPS14062)
    , codeletsPerTP14062(this->numThreads / NUMTPS14062)
    , totalCodelets14062(this->TPsToUse14062 * this->codeletsPerTP14062)
    , barrierCodelets13987(new _barrierCodelets13987[1])
    , checkInCodelets13989(new _checkInCodelets13989[this->numThreads])
    , barrierCodelets13989(new _barrierCodelets13989[1])
    , checkInCodelets14045(new _checkInCodelets14045[this->numThreads])
    , checkInCodelets14044(new _checkInCodelets14044[this->numThreads])
    , checkInCodelets14048(new _checkInCodelets14048[this->numThreads])
    , barrierCodelets14052(new _barrierCodelets14052[1])
    , checkInCodelets14054(new _checkInCodelets14054[this->numThreads])
    , checkInCodelets14053(new _checkInCodelets14053[this->numThreads])
    , checkInCodelets14057(new _checkInCodelets14057[this->numThreads])
    , barrierCodelets14061(new _barrierCodelets14061[1])
    , checkInCodelets14062(new _checkInCodelets14062[this->numThreads])
    , barrierCodelets14062(new _barrierCodelets14062[1])
{
    memset((void*)TP14044_LoopCounterPerThread, 0, this->numThreads * sizeof(unsigned int));
    memset((void*)TP14053_LoopCounterPerThread, 0, this->numThreads * sizeof(unsigned int));
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets13987[0] = _barrierCodelets13987(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets14062[0] = _barrierCodelets14062(NUMTPS14062, NUMTPS14062, this, 0);
    barrierCodelets14061[0] = _barrierCodelets14061(this->numThreads, this->numThreads, this, 0);
    barrierCodelets14052[0] = _barrierCodelets14052(this->numThreads, this->numThreads, this, 0);
    barrierCodelets13989[0] = _barrierCodelets13989(NUMTPS13989, NUMTPS13989, this, 0);
    _checkInCodelets14062* checkInCodelets14062Ptr = (this->checkInCodelets14062);
    for (int i = 0; i < NUMTPS14062; i++) {
        TP14062Ptr[i] = nullptr;
        TP14062_alreadyLaunched[i] = 0;
    }
    _checkInCodelets14057* checkInCodelets14057Ptr = (this->checkInCodelets14057);
    _checkInCodelets14053* checkInCodelets14053Ptr = (this->checkInCodelets14053);
    _checkInCodelets14054* checkInCodelets14054Ptr = (this->checkInCodelets14054);
    _checkInCodelets14048* checkInCodelets14048Ptr = (this->checkInCodelets14048);
    _checkInCodelets14044* checkInCodelets14044Ptr = (this->checkInCodelets14044);
    _checkInCodelets14045* checkInCodelets14045Ptr = (this->checkInCodelets14045);
    _checkInCodelets13989* checkInCodelets13989Ptr = (this->checkInCodelets13989);
    for (int i = 0; i < NUMTPS13989; i++) {
        TP13989Ptr[i] = nullptr;
        TP13989_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets14062Ptr) = _checkInCodelets14062(1, 1, this, codeletCounter);
        checkInCodelets14062Ptr++;
        (*checkInCodelets14057Ptr) = _checkInCodelets14057(1, 1, this, codeletCounter);
        checkInCodelets14057Ptr++;
        (*checkInCodelets14053Ptr) = _checkInCodelets14053(1, 1, this, codeletCounter);
        checkInCodelets14053Ptr++;
        (*checkInCodelets14054Ptr) = _checkInCodelets14054(1, 1, this, codeletCounter);
        checkInCodelets14054Ptr++;
        (*checkInCodelets14048Ptr) = _checkInCodelets14048(1, 1, this, codeletCounter);
        checkInCodelets14048Ptr++;
        (*checkInCodelets14044Ptr) = _checkInCodelets14044(1, 1, this, codeletCounter);
        checkInCodelets14044Ptr++;
        (*checkInCodelets14045Ptr) = _checkInCodelets14045(1, 1, this, codeletCounter);
        checkInCodelets14045Ptr++;
        (*checkInCodelets13989Ptr) = _checkInCodelets13989(1, 1, this, codeletCounter);
        (*checkInCodelets13989Ptr).decDep();
        checkInCodelets13989Ptr++;
    }
}
TP13987::~TP13987()
{
    delete[] TP14044_LoopCounterPerThread;
    delete[] TP14053_LoopCounterPerThread;
    delete[] barrierCodelets13987;
    delete[] barrierCodelets14062;
    delete[] checkInCodelets14062;
    delete[] barrierCodelets14061;
    delete[] checkInCodelets14057;
    delete[] checkInCodelets14053;
    delete[] checkInCodelets14054;
    delete[] barrierCodelets14052;
    delete[] checkInCodelets14048;
    delete[] checkInCodelets14044;
    delete[] checkInCodelets14045;
    delete[] barrierCodelets13989;
    delete[] checkInCodelets13989;
}
/*TP13989: OMPForDirective*/
void TP13989::_barrierCodelets13989::fire(void)
{
    TP13989* myTP = static_cast<TP13989*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13989[0].decDep();
}
bool TP13989::requestNewRangeIterations13989(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13989 * codeletID;
        int tempEndRange = rangePerCodelet13989 * (codeletID + 1);
        if (remainderRange13989 != 0) {
            if (codeletID < (uint32_t)remainderRange13989) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13989;
                tempEndRange += remainderRange13989;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13989;
        tempEndRange = tempEndRange * 1 + minIteration13989;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13989 < lastIteration13989) {
            (this->inputsTPParent->i_darts13989[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts13989[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13989;
        }
    }
    return isThereNewIteration;
}
void TP13989::_checkInCodelets13990::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 13990: ForStmt*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    int* i = &(this->inputsTPParent->i_darts13989[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts13989[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts13989[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts13989[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13989(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13989[0].decDep();
        return;
    }
    for (int i_darts_counter_temp13989 = (*i); i_darts_counter_temp13989 <= endRange
         && i_darts_counter_temp13989 <= this->inputsTPParent->lastIteration13989;
         i_darts_counter_temp13989++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp13989 = (*j);
                for (; j_darts_counter_temp13989 <= jend; j_darts_counter_temp13989++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp13989 = (*k);
                        for (; k_darts_counter_temp13989 <= nz - 2; k_darts_counter_temp13989++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp13989 = (*m);
                                for (; m_darts_counter_temp13989 < 5; m_darts_counter_temp13989++) {
                                    rsd[(i_darts_counter_temp13989)][j_darts_counter_temp13989]
                                       [k_darts_counter_temp13989][m_darts_counter_temp13989]
                                        = dt
                                        * rsd[(i_darts_counter_temp13989)]
                                             [j_darts_counter_temp13989][k_darts_counter_temp13989]
                                             [m_darts_counter_temp13989];
                                }
                                (*m) = m_darts_counter_temp13989;
                            }
                        }
                        (*k) = k_darts_counter_temp13989;
                    }
                }
                (*j) = j_darts_counter_temp13989;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13989[0].decDep();
}
TP13989::TP13989(int in_numThreads, int in_mainCodeletID, TP13987* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13989** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13989(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts13989(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts13989(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts13989(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration13989(in_initIteration)
    , lastIteration13989(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13989(new _barrierCodelets13989[1])
    , checkInCodelets13990(new _checkInCodelets13990[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13989 = abs(lastIteration13989 - initIteration13989) / 1;
    rangePerCodelet13989 = range13989 / numThreads;
    minIteration13989 = min<int>(lastIteration13989, initIteration13989);
    remainderRange13989 = range13989 % numThreads;
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets13989[0] = _barrierCodelets13989(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13990* checkInCodelets13990Ptr = (this->checkInCodelets13990);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13990);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13990Ptr) = _checkInCodelets13990(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13990Ptr) = _checkInCodelets13990(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13990Ptr).decDep();
        checkInCodelets13990Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13989::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13990[localID].setID(codeletID);
    this->checkInCodelets13990[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13990[localID + this->baseNumThreads * i]
            = _checkInCodelets13990(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13990[localID + this->baseNumThreads * i]
            = _checkInCodelets13990(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13990[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13990[localID + this->baseNumThreads * i].decDep();
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
TP13989::~TP13989()
{
    delete[] barrierCodelets13989;
    delete[] checkInCodelets13990;
}
/*TP14044: ForStmt*/
void TP14044::_checkInCodelets14050::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif

    /*printing node 14050: CallExpr*/
    if (!__sync_val_compare_and_swap(&(myTP->controlTPParent->TP14050_alreadyLaunched), 0, 1)) {
        /*Make the function call*/
        invoke<TP_jacld>(myTP, myTP->numThreads, this->getID(),
            &(myTP->controlTPParent->checkInCodelets14051[this->getID()]),
            &(myTP->controlTPParent->TPParent->barrierCodelets14052[0]),
            &(myTP->controlTPParent->TP14050Ptr),
            (this->inputsTPParent->k_darts13987[this->getID()]));
    } else {
        if (myTP->controlTPParent->TP14050Ptr == nullptr) {
            myTP->add(this);
            return;
        } else {
            myTP->controlTPParent->TP14050Ptr->setNewInputs(
                (this->inputsTPParent->k_darts13987[this->getID()]), this->getID());
            myTP->controlTPParent->TP14050Ptr->nextCodeletsjacld[this->getID()]
                = &(myTP->controlTPParent->checkInCodelets14051[this->getID()]);
            myTP->controlTPParent->TP14050Ptr->nextSyncCodeletsjacld[this->getID()]
                = &(myTP->controlTPParent->TPParent->barrierCodelets14052[0]);
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14050Ptr->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14050Ptr->availableCodelets[this->getID()] = 1;
#endif
        }
    }
}
void TP14044::_checkInCodelets14051::fire(void)
{

    /*printing node 14051: CallExpr*/
    if (!__sync_val_compare_and_swap(&(myTP->controlTPParent->TP14051_alreadyLaunched), 0, 1)) {
        /*Make the function call*/
        invoke<TP_blts>(myTP, myTP->numThreads, this->getID(),
            &(myTP->controlTPParent->TPParent->checkInCodelets14048[this->getID()]),
            &(myTP->controlTPParent->TPParent->barrierCodelets14052[0]),
            &(myTP->controlTPParent->TP14051Ptr), nx, ny, nz,
            (this->inputsTPParent->k_darts13987[this->getID()]), omega, ist, iend, jst, jend, nx0,
            ny0);
    } else {
        if (myTP->controlTPParent->TP14051Ptr == nullptr) {
            myTP->add(this);
            return;
        } else {
            myTP->controlTPParent->TP14051Ptr->setNewInputs(nx, ny, nz,
                (this->inputsTPParent->k_darts13987[this->getID()]), omega, ist, iend, jst, jend,
                nx0, ny0, this->getID());
            myTP->controlTPParent->TP14051Ptr->nextCodeletsblts[this->getID()]
                = &(myTP->controlTPParent->TPParent->checkInCodelets14048[this->getID()]);
            myTP->controlTPParent->TP14051Ptr->nextSyncCodeletsblts[this->getID()]
                = &(myTP->controlTPParent->TPParent->barrierCodelets14052[0]);
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14051Ptr->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14051Ptr->availableCodelets[this->getID()] = 1;
#endif
        }
    }
}
TP14044::TP14044(int in_numThreads, int in_mainCodeletID, TP13987* in_TPParent,
    TP13987* in_inputsTPParent, TP14044** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(in_inputsTPParent)
    , ptrToThisTP(in_ptrToThisTP)
    , TP14050Ptr(nullptr)
    , TP14050_alreadyLaunched(0)
    , TP14051Ptr(nullptr)
    , TP14051_alreadyLaunched(0)
    , checkInCodelets14050(new _checkInCodelets14050[this->numThreads])
    , checkInCodelets14051(new _checkInCodelets14051[this->numThreads])
{
    /*Initialize Codelets*/
    _checkInCodelets14051* checkInCodelets14051Ptr = (this->checkInCodelets14051);
    _checkInCodelets14050* checkInCodelets14050Ptr = (this->checkInCodelets14050);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets14050);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets14051Ptr) = _checkInCodelets14051(1, 1, this, codeletCounter);
        checkInCodelets14051Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets14050Ptr) = _checkInCodelets14050(2, 1, this, codeletCounter);
#else
        (*checkInCodelets14050Ptr) = _checkInCodelets14050(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets14050Ptr).decDep();
        checkInCodelets14050Ptr++;
    }
    *(this->ptrToThisTP) = this;
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[this->getID()].decDep();
#else
    this->availableCodelets[this->getID()] = 1;
#endif
}
TP14044::~TP14044()
{
    delete[] checkInCodelets14051;
    delete[] checkInCodelets14050;
}
/*TP14053: ForStmt*/
void TP14053::_checkInCodelets14059::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif

    /*printing node 14059: CallExpr*/
    if (!__sync_val_compare_and_swap(&(myTP->controlTPParent->TP14059_alreadyLaunched), 0, 1)) {
        /*Make the function call*/
        invoke<TP_jacu>(myTP, myTP->numThreads, this->getID(),
            &(myTP->controlTPParent->checkInCodelets14060[this->getID()]),
            &(myTP->controlTPParent->TPParent->barrierCodelets14061[0]),
            &(myTP->controlTPParent->TP14059Ptr),
            (this->inputsTPParent->k_darts13987[this->getID()]));
    } else {
        if (myTP->controlTPParent->TP14059Ptr == nullptr) {
            myTP->add(this);
            return;
        } else {
            myTP->controlTPParent->TP14059Ptr->setNewInputs(
                (this->inputsTPParent->k_darts13987[this->getID()]), this->getID());
            myTP->controlTPParent->TP14059Ptr->nextCodeletsjacu[this->getID()]
                = &(myTP->controlTPParent->checkInCodelets14060[this->getID()]);
            myTP->controlTPParent->TP14059Ptr->nextSyncCodeletsjacu[this->getID()]
                = &(myTP->controlTPParent->TPParent->barrierCodelets14061[0]);
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14059Ptr->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14059Ptr->availableCodelets[this->getID()] = 1;
#endif
        }
    }
}
void TP14053::_checkInCodelets14060::fire(void)
{

    /*printing node 14060: CallExpr*/
    if (!__sync_val_compare_and_swap(&(myTP->controlTPParent->TP14060_alreadyLaunched), 0, 1)) {
        /*Make the function call*/
        invoke<TP_buts>(myTP, myTP->numThreads, this->getID(),
            &(myTP->controlTPParent->TPParent->checkInCodelets14057[this->getID()]),
            &(myTP->controlTPParent->TPParent->barrierCodelets14061[0]),
            &(myTP->controlTPParent->TP14060Ptr), nx, ny, nz,
            (this->inputsTPParent->k_darts13987[this->getID()]), omega, ist, iend, jst, jend, nx0,
            ny0);
    } else {
        if (myTP->controlTPParent->TP14060Ptr == nullptr) {
            myTP->add(this);
            return;
        } else {
            myTP->controlTPParent->TP14060Ptr->setNewInputs(nx, ny, nz,
                (this->inputsTPParent->k_darts13987[this->getID()]), omega, ist, iend, jst, jend,
                nx0, ny0, this->getID());
            myTP->controlTPParent->TP14060Ptr->nextCodeletsbuts[this->getID()]
                = &(myTP->controlTPParent->TPParent->checkInCodelets14057[this->getID()]);
            myTP->controlTPParent->TP14060Ptr->nextSyncCodeletsbuts[this->getID()]
                = &(myTP->controlTPParent->TPParent->barrierCodelets14061[0]);
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14060Ptr->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14060Ptr->availableCodelets[this->getID()] = 1;
#endif
        }
    }
}
TP14053::TP14053(int in_numThreads, int in_mainCodeletID, TP13987* in_TPParent,
    TP13987* in_inputsTPParent, TP14053** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(in_inputsTPParent)
    , ptrToThisTP(in_ptrToThisTP)
    , TP14059Ptr(nullptr)
    , TP14059_alreadyLaunched(0)
    , TP14060Ptr(nullptr)
    , TP14060_alreadyLaunched(0)
    , checkInCodelets14059(new _checkInCodelets14059[this->numThreads])
    , checkInCodelets14060(new _checkInCodelets14060[this->numThreads])
{
    /*Initialize Codelets*/
    _checkInCodelets14060* checkInCodelets14060Ptr = (this->checkInCodelets14060);
    _checkInCodelets14059* checkInCodelets14059Ptr = (this->checkInCodelets14059);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets14059);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets14060Ptr) = _checkInCodelets14060(1, 1, this, codeletCounter);
        checkInCodelets14060Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets14059Ptr) = _checkInCodelets14059(2, 1, this, codeletCounter);
#else
        (*checkInCodelets14059Ptr) = _checkInCodelets14059(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets14059Ptr).decDep();
        checkInCodelets14059Ptr++;
    }
    *(this->ptrToThisTP) = this;
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[this->getID()].decDep();
#else
    this->availableCodelets[this->getID()] = 1;
#endif
}
TP14053::~TP14053()
{
    delete[] checkInCodelets14060;
    delete[] checkInCodelets14059;
}
/*TP14062: OMPForDirective*/
void TP14062::_barrierCodelets14062::fire(void)
{
    TP14062* myTP = static_cast<TP14062*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets14062[0].decDep();
}
bool TP14062::requestNewRangeIterations14062(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet14062 * codeletID;
        int tempEndRange = rangePerCodelet14062 * (codeletID + 1);
        if (remainderRange14062 != 0) {
            if (codeletID < (uint32_t)remainderRange14062) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange14062;
                tempEndRange += remainderRange14062;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration14062;
        tempEndRange = tempEndRange * 1 + minIteration14062;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration14062 < lastIteration14062) {
            (this->inputsTPParent->i_darts14062[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts14062[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration14062;
        }
    }
    return isThereNewIteration;
}
void TP14062::_checkInCodelets14063::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 14063: ForStmt*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    /*var: tmp*/
    int* i = &(this->inputsTPParent->i_darts14062[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts14062[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts14062[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts14062[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double* tmp = (this->inputsTPParent->tmp_darts14062);
    (void)tmp /*OMP_SHARED*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations14062(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets14062[0].decDep();
        return;
    }
    for (int i_darts_counter_temp14062 = (*i); i_darts_counter_temp14062 <= endRange
         && i_darts_counter_temp14062 <= this->inputsTPParent->lastIteration14062;
         i_darts_counter_temp14062++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp14062 = (*j);
                for (; j_darts_counter_temp14062 <= jend; j_darts_counter_temp14062++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp14062 = (*k);
                        for (; k_darts_counter_temp14062 <= nz - 2; k_darts_counter_temp14062++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp14062 = (*m);
                                for (; m_darts_counter_temp14062 < 5; m_darts_counter_temp14062++) {
                                    u[(i_darts_counter_temp14062)][j_darts_counter_temp14062]
                                     [k_darts_counter_temp14062][m_darts_counter_temp14062]
                                        = u[(i_darts_counter_temp14062)][j_darts_counter_temp14062]
                                           [k_darts_counter_temp14062][m_darts_counter_temp14062]
                                        + (*(tmp))
                                            * rsd[(i_darts_counter_temp14062)]
                                                 [j_darts_counter_temp14062]
                                                 [k_darts_counter_temp14062]
                                                 [m_darts_counter_temp14062];
                                }
                                (*m) = m_darts_counter_temp14062;
                            }
                        }
                        (*k) = k_darts_counter_temp14062;
                    }
                }
                (*j) = j_darts_counter_temp14062;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets14062[0].decDep();
}
TP14062::TP14062(int in_numThreads, int in_mainCodeletID, TP13987* in_TPParent,
    int in_initIteration, int in_lastIteration, double* in_tmp, TP14062** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts14062(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts14062(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts14062(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts14062(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts14062(in_tmp) /*OMP_SHARED - INPUT*/
    , initIteration14062(in_initIteration)
    , lastIteration14062(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets14062(new _barrierCodelets14062[1])
    , checkInCodelets14063(new _checkInCodelets14063[this->numThreads])
{
    /*Initialize the loop parameters*/
    range14062 = abs(lastIteration14062 - initIteration14062) / 1;
    rangePerCodelet14062 = range14062 / numThreads;
    minIteration14062 = min<int>(lastIteration14062, initIteration14062);
    remainderRange14062 = range14062 % numThreads;
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets14062[0] = _barrierCodelets14062(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets14063* checkInCodelets14063Ptr = (this->checkInCodelets14063);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets14063);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets14063Ptr) = _checkInCodelets14063(2, 1, this, codeletCounter);
#else
        (*checkInCodelets14063Ptr) = _checkInCodelets14063(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets14063Ptr).decDep();
        checkInCodelets14063Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP14062::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets14063[localID].setID(codeletID);
    this->checkInCodelets14063[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets14063[localID + this->baseNumThreads * i]
            = _checkInCodelets14063(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets14063[localID + this->baseNumThreads * i]
            = _checkInCodelets14063(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets14063[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets14063[localID + this->baseNumThreads * i].decDep();
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
TP14062::~TP14062()
{
    delete[] barrierCodelets14062;
    delete[] checkInCodelets14063;
}
