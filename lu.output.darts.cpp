#include "lu.output.darts.h"
using namespace darts;
using namespace std;
std::mutex TP10026mutex;
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
    domain();
    setcoeff();
    setbv();
    setiv();
    erhs();
    if (affinMaskRes) {
        myDARTSRuntime->run(launch<TP186>(
            ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet, (int*)&((nthreads))));
    }
    ssor();
    error();
    pintgr();
    verify(rsdnm, errnm, frc, &class_is, &verified);
    mflops = (double)itmax
        * (1984.77 * (double)nx0 * (double)ny0 * (double)nz0
            - 10923.299999999999
                * (((double)(nx0 + ny0 + nz0) / 3.) * ((double)(nx0 + ny0 + nz0) / 3.))
            + 27770.900000000001 * (double)(nx0 + ny0 + nz0) / 3. - 144010.)
        / (maxtime * 1.0E+6);
    c_print_results("LU", class_is, nx0, ny0, nz0, itmax, nthreads, maxtime, mflops,
        "          floating point", verified, "3.0 structured", "09 Sep 2022", "gcc", "gcc",
        "-lm -fopenmp", "-I../common", "-O3 -fopenmp", "(none)", "(none)");
}
/*Function: domain, ID: 3*/
static void domain()
{
    /*domain:3*/
    /*CompoundStmt:2220*/
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
    /*CompoundStmt:2240*/
    if (affinMaskRes) {
        myDARTSRuntime->run(
            launch<TP2241>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
    }
}
/*Function: error, ID: 5*/
static void error()
{
    /*error:5*/
    /*CompoundStmt:4767*/
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
    /*CompoundStmt:4831*/
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
    /*CompoundStmt:9912*/
    if (affinMaskRes) {
        myDARTSRuntime->run(launch<TP9913>(ompNumThreads * DARTS_CODELETS_MULT, 0,
            RuntimeFinalCodelet, (int*)&((iend)), (int*)&((ist)), (int*)&((jend)), (int*)&((jst)),
            (int*)&((nx0)), (int*)&((ny0)), (int*)&((nz0)), (double**)&((sum))));
    }
}
/*Function: pintgr, ID: 10*/
static void pintgr()
{
    /*pintgr:10*/
    /*CompoundStmt:10055*/
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
    /*CompoundStmt:10685*/
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
    /*CompoundStmt:10829*/
    if (affinMaskRes) {
        myDARTSRuntime->run(
            launch<TP10830>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
    }
}
/*Function: setbv, ID: 13*/
static void setbv()
{
    /*setbv:13*/
    /*CompoundStmt:13230*/
    if (affinMaskRes) {
        myDARTSRuntime->run(
            launch<TP13231>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
    }
}
/*Function: setcoeff, ID: 14*/
static void setcoeff()
{
    /*setcoeff:14*/
    /*CompoundStmt:13458*/
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
    /*CompoundStmt:13797*/
    int i, j, k, m;
    int iglob, jglob;
    double xi, eta, zeta;
    double pxi, peta, pzeta;
    for (j = 0; j < ny; j++) {
        jglob = j;
        for (k = 1; k < nz - 1; k++) {
            zeta = ((double)k) / (nz - 1);
            if (jglob != 0 && jglob != ny0 - 1) {
                eta = ((double)(jglob)) / (ny0 - 1);
                for (i = 0; i < nx; i++) {
                    iglob = i;
                    if (iglob != 0 && iglob != nx0 - 1) {
                        xi = ((double)(iglob)) / (nx0 - 1);
                        exact(0, jglob, k, ue_1jk);
                        exact(nx0 - 1, jglob, k, ue_nx0jk);
                        exact(iglob, 0, k, ue_i1k);
                        exact(iglob, ny0 - 1, k, ue_iny0k);
                        exact(iglob, jglob, 0, ue_ij1);
                        exact(iglob, jglob, nz - 1, ue_ijnz);
                        for (m = 0; m < 5; m++) {
                            pxi = (1. - xi) * ue_1jk[m] + xi * ue_nx0jk[m];
                            peta = (1. - eta) * ue_i1k[m] + eta * ue_iny0k[m];
                            pzeta = (1. - zeta) * ue_ij1[m] + zeta * ue_ijnz[m];
                            u[i][j][k][m] = pxi + peta + pzeta - pxi * peta - peta * pzeta
                                - pzeta * pxi + pxi * peta * pzeta;
                        }
                    }
                }
            }
        }
    }
}
/*Function: ssor, ID: 16*/
static void ssor()
{
    /*ssor:16*/
    /*CompoundStmt:13892*/
    int i, j, k, m;
    int istep;
    double tmp;
    double delunm[5];
    tmp = 1. / (omega * (2. - omega));
    if (affinMaskRes) {
        myDARTSRuntime->run(
            launch<TP13903>(ompNumThreads * DARTS_CODELETS_MULT, 0, RuntimeFinalCodelet));
    }
    rhs();
    l2norm(nx0, ny0, nz0, ist, iend, jst, jend, rsdnm);
    timer_clear(1);
    timer_start(1);
    for (istep = 1; istep <= itmax; istep++) {
        /*CompoundStmt:13996*/
        if (istep % 20 == 0 || istep == itmax || istep == 1) {
            //#pragma omp master
            printf(" Time step %4d\n", istep);
        }
        if (affinMaskRes) {
            myDARTSRuntime->run(launch<TP14002>(
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
    /*CompoundStmt:14157*/
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
/*TP186: OMPParallelDirective*/
void TP186::_barrierCodelets186::fire(void)
{
    TP186* myTP = static_cast<TP186*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP186::_checkInCodelets188::fire(void)
{
    if (this->getID() == 0) {
        /*Init the vars for this region*/
        /*Initialize the vars of the inlined region*/
        this->inputsTPParent->nthreads_darts188
            = (this->inputsTPParent->nthreads_darts186) /*OMP_SHARED - VAR INLINED*/;

        /*printing node 189: BinaryOperator*/
        (*(this->inputsTPParent->nthreads_darts188)) = omp_get_num_threads();
        /*Signaling next codelet from last stmt in the codelet*/
        /*Find and signal the next codelet*/
        myTP->controlTPParent->TPParent->barrierCodelets186[0].decDep();
    } else {
        /*Find and signal the next codelet*/
        myTP->TPParent->barrierCodelets186[0].decDep();
    }
}
TP186::TP186(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet, int* in_nthreads)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , nthreads_darts186(in_nthreads) /*OMP_SHARED - INPUT*/
    , TP188_alreadyLaunched(0)
    , barrierCodelets186(new _barrierCodelets186[1])
    , checkInCodelets188(new _checkInCodelets188[this->numThreads])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets186[0] = _barrierCodelets186(ompNumThreads, ompNumThreads, this, 0);
    _checkInCodelets188* checkInCodelets188Ptr = (this->checkInCodelets188);
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets188Ptr) = _checkInCodelets188(1, 1, this, codeletCounter);
        (*checkInCodelets188Ptr).decDep();
        checkInCodelets188Ptr++;
    }
}
TP186::~TP186()
{
    delete[] barrierCodelets186;
    delete[] checkInCodelets188;
}
/*TP1: TP_blts*/
void TP1::_checkInCodelets238::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*region 238 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP238;
    if (idx < myTP->TPsToUse238) {
        if (!__sync_val_compare_and_swap(&(myTP->TP238_alreadyLaunched[idx]), 0, 1)) {
            int range = abs((this->inputsTPParent->iend_darts1[this->getID()])
                            - (this->inputsTPParent->ist_darts1[this->getID()]))
                / 1;
            int rangePerCodelet = range / myTP->TPsToUse238;
            int minIteration = min<int>((this->inputsTPParent->iend_darts1[this->getID()]),
                (this->inputsTPParent->ist_darts1[this->getID()]));
            int remainderRange = range % myTP->TPsToUse238;
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
            if (idx == myTP->TPsToUse238 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse238 - 1) {
                lastIteration = (this->inputsTPParent->iend_darts1[this->getID()]);
            }
#if USEINVOKE == 1
            invoke<TP238>(myTP, myTP->codeletsPerTP238 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP238Ptr[idx]));
#else
            place<TP238>(idx, myTP, myTP->codeletsPerTP238 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP238Ptr[idx]));
#endif
        } else {
            if (myTP->TP238Ptr[idx] != nullptr) {
                myTP->TP238Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP1::_barrierCodelets238::fire(void)
{
    TP1* myTP = static_cast<TP1*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets342[codeletsCounter].decDep();
        }
    }
}
void TP1::_checkInCodelets342::fire(void)
{
    if (this->getID() == 0) {
        /*Init the vars for this region*/

        /*printing node 343: CallExpr*/
        printf("part1\n");
        /*Signaling next codelet from last stmt in the codelet*/
        /*Signaling next codelet region: 342 nextRegion: 344 */
        myTP->controlTPParent->checkInCodelets344[this->getID()].decDep();
    } else {
        /*Signaling next codelet region: 342 nextRegion: 344 */
        myTP->checkInCodelets344[this->getID()].decDep();
    }
}
void TP1::_checkInCodelets344::fire(void)
{
    /*region 344 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP344;
    if (idx < myTP->TPsToUse344) {
        if (!__sync_val_compare_and_swap(&(myTP->TP344_alreadyLaunched[idx]), 0, 1)) {
            int range = abs((this->inputsTPParent->iend_darts1[this->getID()])
                            - (this->inputsTPParent->ist_darts1[this->getID()]))
                / 1;
            int rangePerCodelet = range / myTP->TPsToUse344;
            int minIteration = min<int>((this->inputsTPParent->iend_darts1[this->getID()]),
                (this->inputsTPParent->ist_darts1[this->getID()]));
            int remainderRange = range % myTP->TPsToUse344;
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
            if (idx == myTP->TPsToUse344 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse344 - 1) {
                lastIteration = (this->inputsTPParent->iend_darts1[this->getID()]);
            }
#if USEINVOKE == 1
            invoke<TP344>(myTP, myTP->codeletsPerTP344 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP344Ptr[idx]));
#else
            place<TP344>(idx, myTP, myTP->codeletsPerTP344 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP344Ptr[idx]));
#endif
        } else {
            if (myTP->TP344Ptr[idx] != nullptr) {
                myTP->TP344Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP1::_barrierCodelets344::fire(void)
{
    TP1* myTP = static_cast<TP1*>(myTP_);

    for (size_t codeletsCounter = 0; codeletsCounter < (size_t)myTP->numThreads;
         codeletsCounter++) {
        myTP->nextCodeletsblts[codeletsCounter]->decDep();
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
    , TP238Ptr(new TP238*[NUMTPS238])
    , TP238_alreadyLaunched(new size_t[NUMTPS238])
    , numTPsSet238(0)
    , numTPsReady238(0)
    , TPsToUse238(NUMTPS238)
    , codeletsPerTP238(this->numThreads / NUMTPS238)
    , totalCodelets238(this->TPsToUse238 * this->codeletsPerTP238)
    , TP342_alreadyLaunched(0)
    , TP344Ptr(new TP344*[NUMTPS344])
    , TP344_alreadyLaunched(new size_t[NUMTPS344])
    , numTPsSet344(0)
    , numTPsReady344(0)
    , TPsToUse344(NUMTPS344)
    , codeletsPerTP344(this->numThreads / NUMTPS344)
    , totalCodelets344(this->TPsToUse344 * this->codeletsPerTP344)
    , checkInCodelets238(new _checkInCodelets238[this->numThreads])
    , barrierCodelets238(new _barrierCodelets238[1])
    , checkInCodelets342(new _checkInCodelets342[this->numThreads])
    , checkInCodelets344(new _checkInCodelets344[this->numThreads])
    , barrierCodelets344(new _barrierCodelets344[1])
{
    barrierCodelets344[0] = _barrierCodelets344(NUMTPS344, NUMTPS344, this, 0);
    barrierCodelets238[0] = _barrierCodelets238(NUMTPS238, NUMTPS238, this, 0);
    _checkInCodelets344* checkInCodelets344Ptr = (this->checkInCodelets344);
    for (int i = 0; i < NUMTPS344; i++) {
        TP344Ptr[i] = nullptr;
        TP344_alreadyLaunched[i] = 0;
    }
    _checkInCodelets342* checkInCodelets342Ptr = (this->checkInCodelets342);
    _checkInCodelets238* checkInCodelets238Ptr = (this->checkInCodelets238);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets238);
#endif
    for (int i = 0; i < NUMTPS238; i++) {
        TP238Ptr[i] = nullptr;
        TP238_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++) {
        (*checkInCodelets344Ptr) = _checkInCodelets344(1, 1, this, codeletCounter);
        checkInCodelets344Ptr++;
        (*checkInCodelets342Ptr) = _checkInCodelets342(1, 1, this, codeletCounter);
        checkInCodelets342Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets238Ptr) = _checkInCodelets238(2, 1, this, codeletCounter);
#else
        (*checkInCodelets238Ptr) = _checkInCodelets238(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets238Ptr).decDep();
        checkInCodelets238Ptr++;
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
    delete[] barrierCodelets344;
    delete[] checkInCodelets344;
    delete[] checkInCodelets342;
    delete[] barrierCodelets238;
    delete[] checkInCodelets238;
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
/*TP238: OMPForDirective*/
void TP238::_barrierCodelets238::fire(void)
{
    TP238* myTP = static_cast<TP238*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets238[0].decDep();
}
bool TP238::requestNewRangeIterations238(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 1*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet238 * codeletID;
        int tempEndRange = rangePerCodelet238 * (codeletID + 1);
        if (remainderRange238 != 0) {
            if (codeletID < (uint32_t)remainderRange238) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange238;
                tempEndRange += remainderRange238;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration238;
        tempEndRange = tempEndRange * 1 + minIteration238;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration238 < lastIteration238) {
            (this->inputsTPParent->i_darts238[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts238[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration238;
        }
    }
    return isThereNewIteration;
}
void TP238::_checkInCodelets239::fire(void)
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
    this->inputsTPParent->iend_darts238[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iend_darts1[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ist_darts238[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->ist_darts1[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jend_darts238[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jend_darts1[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jst_darts238[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jst_darts1[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->k_darts238[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->k_darts1[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->omega_darts238[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->omega_darts1[this->getID()]);

    /*printing node 239: ForStmt*/
    /*var: i*/
    /*var: iend*/
    /*var: ist*/
    /*var: j*/
    /*var: jend*/
    /*var: jst*/
    /*var: k*/
    /*var: m*/
    /*var: omega*/
    int* i = &(this->inputsTPParent->i_darts238[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts238[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jend = &(this->inputsTPParent->jend_darts238[this->getLocalID()]);
    (void)jend /*OMP_SHARED_PRIVATE*/;
    int** jst = &(this->inputsTPParent->jst_darts238[this->getLocalID()]);
    (void)jst /*OMP_SHARED_PRIVATE*/;
    int** k = &(this->inputsTPParent->k_darts238[this->getLocalID()]);
    (void)k /*OMP_SHARED_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts238[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** omega = &(this->inputsTPParent->omega_darts238[this->getLocalID()]);
    (void)omega /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations238(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets238[0].decDep();
        return;
    }
    for (int i_darts_counter_temp238 = (*i); i_darts_counter_temp238 <= endRange
         && i_darts_counter_temp238 <= this->inputsTPParent->lastIteration238;
         i_darts_counter_temp238++) {
        {
            {
                /*Loop's init*/
                (*j) = (*(*jst));
                int j_darts_counter_temp238 = (*j);
                for (; j_darts_counter_temp238 <= (*(*jend)); j_darts_counter_temp238++) {
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp238 = (*m);
                        for (; m_darts_counter_temp238 < 5; m_darts_counter_temp238++) {
                            rsd[(i_darts_counter_temp238)][j_darts_counter_temp238][(*(*k))]
                               [m_darts_counter_temp238]
                                = rsd[(i_darts_counter_temp238)][j_darts_counter_temp238][(*(*k))]
                                     [m_darts_counter_temp238]
                                - (*(*omega))
                                    * (a[(i_darts_counter_temp238)][j_darts_counter_temp238]
                                        [m_darts_counter_temp238][0]
                                            * rsd[(i_darts_counter_temp238)]
                                                 [j_darts_counter_temp238][(*(*k)) - 1][0]
                                        + a[(i_darts_counter_temp238)][j_darts_counter_temp238]
                                           [m_darts_counter_temp238][1]
                                            * rsd[(i_darts_counter_temp238)]
                                                 [j_darts_counter_temp238][(*(*k)) - 1][1]
                                        + a[(i_darts_counter_temp238)][j_darts_counter_temp238]
                                           [m_darts_counter_temp238][2]
                                            * rsd[(i_darts_counter_temp238)]
                                                 [j_darts_counter_temp238][(*(*k)) - 1][2]
                                        + a[(i_darts_counter_temp238)][j_darts_counter_temp238]
                                           [m_darts_counter_temp238][3]
                                            * rsd[(i_darts_counter_temp238)]
                                                 [j_darts_counter_temp238][(*(*k)) - 1][3]
                                        + a[(i_darts_counter_temp238)][j_darts_counter_temp238]
                                           [m_darts_counter_temp238][4]
                                            * rsd[(i_darts_counter_temp238)]
                                                 [j_darts_counter_temp238][(*(*k)) - 1][4]);
                        }
                        (*m) = m_darts_counter_temp238;
                    }
                }
                (*j) = j_darts_counter_temp238;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets238[0].decDep();
}
TP238::TP238(int in_numThreads, int in_mainCodeletID, TP1* in_TPParent, int in_initIteration,
    int in_lastIteration, TP238** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts238(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iend_darts238(new int*[this->numThreads])
    , ist_darts238(new int*[this->numThreads])
    , j_darts238(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jend_darts238(new int*[this->numThreads])
    , jst_darts238(new int*[this->numThreads])
    , k_darts238(new int*[this->numThreads])
    , m_darts238(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , omega_darts238(new double*[this->numThreads])
    , initIteration238(in_initIteration)
    , lastIteration238(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets238(new _barrierCodelets238[1])
    , checkInCodelets239(new _checkInCodelets239[this->numThreads])
{
    /*Initialize the loop parameters*/
    range238 = abs(lastIteration238 - initIteration238) / 1;
    rangePerCodelet238 = range238 / numThreads;
    minIteration238 = min<int>(lastIteration238, initIteration238);
    remainderRange238 = range238 % numThreads;
    /*Initialize inputs and vars.*/
    this->iend_darts238 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ist_darts238 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jend_darts238 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jst_darts238 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->k_darts238 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->omega_darts238
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets238[0] = _barrierCodelets238(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets239* checkInCodelets239Ptr = (this->checkInCodelets239);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets239);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets239Ptr) = _checkInCodelets239(2, 1, this, codeletCounter);
#else
        (*checkInCodelets239Ptr) = _checkInCodelets239(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets239Ptr).decDep();
        checkInCodelets239Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP238::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets239[localID].setID(codeletID);
    this->checkInCodelets239[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets239[localID + this->baseNumThreads * i]
            = _checkInCodelets239(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets239[localID + this->baseNumThreads * i]
            = _checkInCodelets239(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets239[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets239[localID + this->baseNumThreads * i].decDep();
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
TP238::~TP238()
{
    delete[] iend_darts238;
    delete[] ist_darts238;
    delete[] jend_darts238;
    delete[] jst_darts238;
    delete[] k_darts238;
    delete[] omega_darts238;
    delete[] barrierCodelets238;
    delete[] checkInCodelets239;
}
/*TP344: OMPForDirective*/
void TP344::_barrierCodelets344::fire(void)
{
    TP344* myTP = static_cast<TP344*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets344[0].decDep();
}
bool TP344::requestNewRangeIterations344(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 1*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet344 * codeletID;
        int tempEndRange = rangePerCodelet344 * (codeletID + 1);
        if (remainderRange344 != 0) {
            if (codeletID < (uint32_t)remainderRange344) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange344;
                tempEndRange += remainderRange344;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration344;
        tempEndRange = tempEndRange * 1 + minIteration344;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration344 < lastIteration344) {
            (this->inputsTPParent->i_darts344[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts344[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration344;
        }
    }
    return isThereNewIteration;
}
void TP344::_checkInCodelets345::fire(void)
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
    this->inputsTPParent->iend_darts344[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iend_darts1[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ist_darts344[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->ist_darts1[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jend_darts344[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jend_darts1[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jst_darts344[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jst_darts1[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->k_darts344[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->k_darts1[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->omega_darts344[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->omega_darts1[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts344[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts1[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp1_darts344[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp1_darts1[this->getID()]);

    /*printing node 345: ForStmt*/
    /*var: i*/
    /*var: iend*/
    /*var: ist*/
    /*var: j*/
    /*var: jend*/
    /*var: jst*/
    /*var: k*/
    /*var: m*/
    /*var: omega*/
    /*var: tmp*/
    /*var: tmp1*/
    int* i = &(this->inputsTPParent->i_darts344[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iend = &(this->inputsTPParent->iend_darts344[this->getLocalID()]);
    (void)iend /*OMP_SHARED_PRIVATE*/;
    int** ist = &(this->inputsTPParent->ist_darts344[this->getLocalID()]);
    (void)ist /*OMP_SHARED_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts344[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jend = &(this->inputsTPParent->jend_darts344[this->getLocalID()]);
    (void)jend /*OMP_SHARED_PRIVATE*/;
    int** jst = &(this->inputsTPParent->jst_darts344[this->getLocalID()]);
    (void)jst /*OMP_SHARED_PRIVATE*/;
    int** k = &(this->inputsTPParent->k_darts344[this->getLocalID()]);
    (void)k /*OMP_SHARED_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts344[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** omega = &(this->inputsTPParent->omega_darts344[this->getLocalID()]);
    (void)omega /*OMP_SHARED_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts344[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** tmp1 = &(this->inputsTPParent->tmp1_darts344[this->getLocalID()]);
    (void)tmp1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations344(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets344[0].decDep();
        return;
    }
    for (int i_darts_counter_temp344 = (*i); i_darts_counter_temp344 <= endRange
         && i_darts_counter_temp344 <= this->inputsTPParent->lastIteration344;
         i_darts_counter_temp344++) {
        {
            if ((i_darts_counter_temp344) != (*(*ist))) {
                while (flag[(i_darts_counter_temp344)-1] == 0) { }
            }
            if ((i_darts_counter_temp344) != (*(*iend))) {
                while (flag[(i_darts_counter_temp344)] == 1) { }
            }
            {
                /*Loop's init*/
                (*j) = (*(*jst));
                int j_darts_counter_temp344 = (*j);
                for (; j_darts_counter_temp344 <= (*(*jend)); j_darts_counter_temp344++) {
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp344 = (*m);
                        for (; m_darts_counter_temp344 < 5; m_darts_counter_temp344++) {
                            rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))]
                               [m_darts_counter_temp344]
                                = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))]
                                     [m_darts_counter_temp344]
                                - (*(*omega))
                                    * (b[(i_darts_counter_temp344)][j_darts_counter_temp344]
                                        [m_darts_counter_temp344][0]
                                            * rsd[(i_darts_counter_temp344)]
                                                 [j_darts_counter_temp344 - 1][(*(*k))][0]
                                        + c[(i_darts_counter_temp344)][j_darts_counter_temp344]
                                           [m_darts_counter_temp344][0]
                                            * rsd[(i_darts_counter_temp344)-1]
                                                 [j_darts_counter_temp344][(*(*k))][0]
                                        + b[(i_darts_counter_temp344)][j_darts_counter_temp344]
                                           [m_darts_counter_temp344][1]
                                            * rsd[(i_darts_counter_temp344)]
                                                 [j_darts_counter_temp344 - 1][(*(*k))][1]
                                        + c[(i_darts_counter_temp344)][j_darts_counter_temp344]
                                           [m_darts_counter_temp344][1]
                                            * rsd[(i_darts_counter_temp344)-1]
                                                 [j_darts_counter_temp344][(*(*k))][1]
                                        + b[(i_darts_counter_temp344)][j_darts_counter_temp344]
                                           [m_darts_counter_temp344][2]
                                            * rsd[(i_darts_counter_temp344)]
                                                 [j_darts_counter_temp344 - 1][(*(*k))][2]
                                        + c[(i_darts_counter_temp344)][j_darts_counter_temp344]
                                           [m_darts_counter_temp344][2]
                                            * rsd[(i_darts_counter_temp344)-1]
                                                 [j_darts_counter_temp344][(*(*k))][2]
                                        + b[(i_darts_counter_temp344)][j_darts_counter_temp344]
                                           [m_darts_counter_temp344][3]
                                            * rsd[(i_darts_counter_temp344)]
                                                 [j_darts_counter_temp344 - 1][(*(*k))][3]
                                        + c[(i_darts_counter_temp344)][j_darts_counter_temp344]
                                           [m_darts_counter_temp344][3]
                                            * rsd[(i_darts_counter_temp344)-1]
                                                 [j_darts_counter_temp344][(*(*k))][3]
                                        + b[(i_darts_counter_temp344)][j_darts_counter_temp344]
                                           [m_darts_counter_temp344][4]
                                            * rsd[(i_darts_counter_temp344)]
                                                 [j_darts_counter_temp344 - 1][(*(*k))][4]
                                        + c[(i_darts_counter_temp344)][j_darts_counter_temp344]
                                           [m_darts_counter_temp344][4]
                                            * rsd[(i_darts_counter_temp344)-1]
                                                 [j_darts_counter_temp344][(*(*k))][4]);
                        }
                        (*m) = m_darts_counter_temp344;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp344 = (*m);
                        for (; m_darts_counter_temp344 < 5; m_darts_counter_temp344++) {
                            tmat[m_darts_counter_temp344][0]
                                = d[(i_darts_counter_temp344)][j_darts_counter_temp344]
                                   [m_darts_counter_temp344][0];
                            tmat[m_darts_counter_temp344][1]
                                = d[(i_darts_counter_temp344)][j_darts_counter_temp344]
                                   [m_darts_counter_temp344][1];
                            tmat[m_darts_counter_temp344][2]
                                = d[(i_darts_counter_temp344)][j_darts_counter_temp344]
                                   [m_darts_counter_temp344][2];
                            tmat[m_darts_counter_temp344][3]
                                = d[(i_darts_counter_temp344)][j_darts_counter_temp344]
                                   [m_darts_counter_temp344][3];
                            tmat[m_darts_counter_temp344][4]
                                = d[(i_darts_counter_temp344)][j_darts_counter_temp344]
                                   [m_darts_counter_temp344][4];
                        }
                        (*m) = m_darts_counter_temp344;
                    }
                    (*(*tmp1)) = 1. / tmat[0][0];
                    (*(*tmp)) = (*(*tmp1)) * tmat[1][0];
                    tmat[1][1] = tmat[1][1] - (*(*tmp)) * tmat[0][1];
                    tmat[1][2] = tmat[1][2] - (*(*tmp)) * tmat[0][2];
                    tmat[1][3] = tmat[1][3] - (*(*tmp)) * tmat[0][3];
                    tmat[1][4] = tmat[1][4] - (*(*tmp)) * tmat[0][4];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][1]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][1]
                        - rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][0]
                            * (*(*tmp));
                    (*(*tmp)) = (*(*tmp1)) * tmat[2][0];
                    tmat[2][1] = tmat[2][1] - (*(*tmp)) * tmat[0][1];
                    tmat[2][2] = tmat[2][2] - (*(*tmp)) * tmat[0][2];
                    tmat[2][3] = tmat[2][3] - (*(*tmp)) * tmat[0][3];
                    tmat[2][4] = tmat[2][4] - (*(*tmp)) * tmat[0][4];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][2]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][2]
                        - rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][0]
                            * (*(*tmp));
                    (*(*tmp)) = (*(*tmp1)) * tmat[3][0];
                    tmat[3][1] = tmat[3][1] - (*(*tmp)) * tmat[0][1];
                    tmat[3][2] = tmat[3][2] - (*(*tmp)) * tmat[0][2];
                    tmat[3][3] = tmat[3][3] - (*(*tmp)) * tmat[0][3];
                    tmat[3][4] = tmat[3][4] - (*(*tmp)) * tmat[0][4];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][3]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][3]
                        - rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][0]
                            * (*(*tmp));
                    (*(*tmp)) = (*(*tmp1)) * tmat[4][0];
                    tmat[4][1] = tmat[4][1] - (*(*tmp)) * tmat[0][1];
                    tmat[4][2] = tmat[4][2] - (*(*tmp)) * tmat[0][2];
                    tmat[4][3] = tmat[4][3] - (*(*tmp)) * tmat[0][3];
                    tmat[4][4] = tmat[4][4] - (*(*tmp)) * tmat[0][4];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][4]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][4]
                        - rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][0]
                            * (*(*tmp));
                    (*(*tmp1)) = 1. / tmat[1][1];
                    (*(*tmp)) = (*(*tmp1)) * tmat[2][1];
                    tmat[2][2] = tmat[2][2] - (*(*tmp)) * tmat[1][2];
                    tmat[2][3] = tmat[2][3] - (*(*tmp)) * tmat[1][3];
                    tmat[2][4] = tmat[2][4] - (*(*tmp)) * tmat[1][4];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][2]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][2]
                        - rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][1]
                            * (*(*tmp));
                    (*(*tmp)) = (*(*tmp1)) * tmat[3][1];
                    tmat[3][2] = tmat[3][2] - (*(*tmp)) * tmat[1][2];
                    tmat[3][3] = tmat[3][3] - (*(*tmp)) * tmat[1][3];
                    tmat[3][4] = tmat[3][4] - (*(*tmp)) * tmat[1][4];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][3]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][3]
                        - rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][1]
                            * (*(*tmp));
                    (*(*tmp)) = (*(*tmp1)) * tmat[4][1];
                    tmat[4][2] = tmat[4][2] - (*(*tmp)) * tmat[1][2];
                    tmat[4][3] = tmat[4][3] - (*(*tmp)) * tmat[1][3];
                    tmat[4][4] = tmat[4][4] - (*(*tmp)) * tmat[1][4];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][4]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][4]
                        - rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][1]
                            * (*(*tmp));
                    (*(*tmp1)) = 1. / tmat[2][2];
                    (*(*tmp)) = (*(*tmp1)) * tmat[3][2];
                    tmat[3][3] = tmat[3][3] - (*(*tmp)) * tmat[2][3];
                    tmat[3][4] = tmat[3][4] - (*(*tmp)) * tmat[2][4];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][3]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][3]
                        - rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][2]
                            * (*(*tmp));
                    (*(*tmp)) = (*(*tmp1)) * tmat[4][2];
                    tmat[4][3] = tmat[4][3] - (*(*tmp)) * tmat[2][3];
                    tmat[4][4] = tmat[4][4] - (*(*tmp)) * tmat[2][4];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][4]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][4]
                        - rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][2]
                            * (*(*tmp));
                    (*(*tmp1)) = 1. / tmat[3][3];
                    (*(*tmp)) = (*(*tmp1)) * tmat[4][3];
                    tmat[4][4] = tmat[4][4] - (*(*tmp)) * tmat[3][4];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][4]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][4]
                        - rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][3]
                            * (*(*tmp));
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][4]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][4]
                        / tmat[4][4];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][3]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][3]
                        - tmat[3][4]
                            * rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][4];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][3]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][3]
                        / tmat[3][3];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][2]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][2]
                        - tmat[2][3]
                            * rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][3]
                        - tmat[2][4]
                            * rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][4];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][2]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][2]
                        / tmat[2][2];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][1]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][1]
                        - tmat[1][2]
                            * rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][2]
                        - tmat[1][3]
                            * rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][3]
                        - tmat[1][4]
                            * rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][4];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][1]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][1]
                        / tmat[1][1];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][0]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][0]
                        - tmat[0][1]
                            * rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][1]
                        - tmat[0][2]
                            * rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][2]
                        - tmat[0][3]
                            * rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][3]
                        - tmat[0][4]
                            * rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][4];
                    rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][0]
                        = rsd[(i_darts_counter_temp344)][j_darts_counter_temp344][(*(*k))][0]
                        / tmat[0][0];
                }
                (*j) = j_darts_counter_temp344;
            }
            if ((i_darts_counter_temp344) != (*(*ist))) {
                flag[(i_darts_counter_temp344)-1] = 0;
            }
            if ((i_darts_counter_temp344) != (*(*iend))) {
                flag[(i_darts_counter_temp344)] = 1;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets344[0].decDep();
}
TP344::TP344(int in_numThreads, int in_mainCodeletID, TP1* in_TPParent, int in_initIteration,
    int in_lastIteration, TP344** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts344(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iend_darts344(new int*[this->numThreads])
    , ist_darts344(new int*[this->numThreads])
    , j_darts344(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jend_darts344(new int*[this->numThreads])
    , jst_darts344(new int*[this->numThreads])
    , k_darts344(new int*[this->numThreads])
    , m_darts344(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , omega_darts344(new double*[this->numThreads])
    , tmp_darts344(new double*[this->numThreads])
    , tmp1_darts344(new double*[this->numThreads])
    , initIteration344(in_initIteration)
    , lastIteration344(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets344(new _barrierCodelets344[1])
    , checkInCodelets345(new _checkInCodelets345[this->numThreads])
{
    /*Initialize the loop parameters*/
    range344 = abs(lastIteration344 - initIteration344) / 1;
    rangePerCodelet344 = range344 / numThreads;
    minIteration344 = min<int>(lastIteration344, initIteration344);
    remainderRange344 = range344 % numThreads;
    /*Initialize inputs and vars.*/
    this->iend_darts344 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ist_darts344 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jend_darts344 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jst_darts344 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->k_darts344 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->omega_darts344
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts344
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp1_darts344
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets344[0] = _barrierCodelets344(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets345* checkInCodelets345Ptr = (this->checkInCodelets345);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets345);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets345Ptr) = _checkInCodelets345(2, 1, this, codeletCounter);
#else
        (*checkInCodelets345Ptr) = _checkInCodelets345(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets345Ptr).decDep();
        checkInCodelets345Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP344::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets345[localID].setID(codeletID);
    this->checkInCodelets345[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets345[localID + this->baseNumThreads * i]
            = _checkInCodelets345(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets345[localID + this->baseNumThreads * i]
            = _checkInCodelets345(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets345[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets345[localID + this->baseNumThreads * i].decDep();
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
TP344::~TP344()
{
    delete[] iend_darts344;
    delete[] ist_darts344;
    delete[] jend_darts344;
    delete[] jst_darts344;
    delete[] k_darts344;
    delete[] omega_darts344;
    delete[] tmp_darts344;
    delete[] tmp1_darts344;
    delete[] barrierCodelets344;
    delete[] checkInCodelets345;
}
/*TP2: TP_buts*/
void TP2::_checkInCodelets1228::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*region 1228 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP1228;
    if (idx < myTP->TPsToUse1228) {
        if (!__sync_val_compare_and_swap(&(myTP->TP1228_alreadyLaunched[idx]), 0, 1)) {
            int range = abs((this->inputsTPParent->ist_darts2[this->getID()])
                            - (this->inputsTPParent->iend_darts2[this->getID()]))
                / 1;
            int rangePerCodelet = range / myTP->TPsToUse1228;
            int minIteration = min<int>((this->inputsTPParent->ist_darts2[this->getID()]),
                (this->inputsTPParent->iend_darts2[this->getID()]));
            int remainderRange = range % myTP->TPsToUse1228;
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
            if (idx == myTP->TPsToUse1228 - 1) {
                lastIteration = (this->inputsTPParent->ist_darts2[this->getID()]);
            }
#if USEINVOKE == 1
            invoke<TP1228>(myTP, myTP->codeletsPerTP1228 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP1228Ptr[idx]));
#else
            place<TP1228>(idx, myTP, myTP->codeletsPerTP1228 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP1228Ptr[idx]));
#endif
        } else {
            if (myTP->TP1228Ptr[idx] != nullptr) {
                myTP->TP1228Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    } else {
        /*Signaling next codelet region: 1228 nextRegion: 1330 */
        myTP->controlTPParent->checkInCodelets1330[this->getID()].decDep();
    }
}
void TP2::_checkInCodelets1330::fire(void)
{
    /*region 1330 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP1330;
    if (idx < myTP->TPsToUse1330) {
        if (!__sync_val_compare_and_swap(&(myTP->TP1330_alreadyLaunched[idx]), 0, 1)) {
            int range = abs((this->inputsTPParent->ist_darts2[this->getID()])
                            - (this->inputsTPParent->iend_darts2[this->getID()]))
                / 1;
            int rangePerCodelet = range / myTP->TPsToUse1330;
            int minIteration = min<int>((this->inputsTPParent->ist_darts2[this->getID()]),
                (this->inputsTPParent->iend_darts2[this->getID()]));
            int remainderRange = range % myTP->TPsToUse1330;
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
            if (idx == myTP->TPsToUse1330 - 1) {
                lastIteration = (this->inputsTPParent->ist_darts2[this->getID()]);
            }
#if USEINVOKE == 1
            invoke<TP1330>(myTP, myTP->codeletsPerTP1330 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP1330Ptr[idx]));
#else
            place<TP1330>(idx, myTP, myTP->codeletsPerTP1330 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP1330Ptr[idx]));
#endif
        } else {
            if (myTP->TP1330Ptr[idx] != nullptr) {
                myTP->TP1330Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    } else {
        /*Find and signal the next codelet*/

        myTP->controlTPParent->nextCodeletsbuts[this->getID()]->decDep();
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
    , TP1228Ptr(new TP1228*[NUMTPS1228])
    , TP1228_alreadyLaunched(new size_t[NUMTPS1228])
    , numTPsSet1228(0)
    , numTPsReady1228(0)
    , TPsToUse1228(NUMTPS1228)
    , codeletsPerTP1228(this->numThreads / NUMTPS1228)
    , totalCodelets1228(this->TPsToUse1228 * this->codeletsPerTP1228)
    , TP1330Ptr(new TP1330*[NUMTPS1330])
    , TP1330_alreadyLaunched(new size_t[NUMTPS1330])
    , numTPsSet1330(0)
    , numTPsReady1330(0)
    , TPsToUse1330(NUMTPS1330)
    , codeletsPerTP1330(this->numThreads / NUMTPS1330)
    , totalCodelets1330(this->TPsToUse1330 * this->codeletsPerTP1330)
    , checkInCodelets1228(new _checkInCodelets1228[this->numThreads])
    , checkInCodelets1330(new _checkInCodelets1330[this->numThreads])
{
    _checkInCodelets1330* checkInCodelets1330Ptr = (this->checkInCodelets1330);
    for (int i = 0; i < NUMTPS1330; i++) {
        TP1330Ptr[i] = nullptr;
        TP1330_alreadyLaunched[i] = 0;
    }
    _checkInCodelets1228* checkInCodelets1228Ptr = (this->checkInCodelets1228);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets1228);
#endif
    for (int i = 0; i < NUMTPS1228; i++) {
        TP1228Ptr[i] = nullptr;
        TP1228_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++) {
        (*checkInCodelets1330Ptr) = _checkInCodelets1330(1, 1, this, codeletCounter);
        checkInCodelets1330Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets1228Ptr) = _checkInCodelets1228(2, 1, this, codeletCounter);
#else
        (*checkInCodelets1228Ptr) = _checkInCodelets1228(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets1228Ptr).decDep();
        checkInCodelets1228Ptr++;
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
    delete[] checkInCodelets1330;
    delete[] checkInCodelets1228;
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
/*TP1228: OMPForDirective*/
bool TP1228::requestNewRangeIterations1228(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 1*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet1228 * codeletID;
        int tempEndRange = rangePerCodelet1228 * (codeletID + 1);
        if (remainderRange1228 != 0) {
            if (codeletID < (uint32_t)remainderRange1228) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange1228;
                tempEndRange += remainderRange1228;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration1228;
        tempEndRange = tempEndRange * 1 + minIteration1228;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration1228 < lastIteration1228) {
            (this->inputsTPParent->i_darts1228[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts1228[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == 0) {
            *endRange = *endRange - 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration1228;
        }
    }
    return isThereNewIteration;
}
void TP1228::_checkInCodelets1229::fire(void)
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
    this->inputsTPParent->iend_darts1228[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iend_darts2[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ist_darts1228[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->ist_darts2[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jend_darts1228[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jend_darts2[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jst_darts1228[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jst_darts2[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->k_darts1228[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->k_darts2[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->omega_darts1228[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->omega_darts2[this->getID()]);

    /*printing node 1229: ForStmt*/
    /*var: i*/
    /*var: iend*/
    /*var: ist*/
    /*var: j*/
    /*var: jend*/
    /*var: jst*/
    /*var: k*/
    /*var: m*/
    /*var: omega*/
    int* i = &(this->inputsTPParent->i_darts1228[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts1228[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jend = &(this->inputsTPParent->jend_darts1228[this->getLocalID()]);
    (void)jend /*OMP_SHARED_PRIVATE*/;
    int** jst = &(this->inputsTPParent->jst_darts1228[this->getLocalID()]);
    (void)jst /*OMP_SHARED_PRIVATE*/;
    int** k = &(this->inputsTPParent->k_darts1228[this->getLocalID()]);
    (void)k /*OMP_SHARED_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts1228[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** omega = &(this->inputsTPParent->omega_darts1228[this->getLocalID()]);
    (void)omega /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations1228(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Find and signal the next codelet*/
        myTP->controlTPParent->TPParent->checkInCodelets1330[this->getID()].decDep();
        return;
    }
    for (int i_darts_counter_temp1228 = (*i); i_darts_counter_temp1228 >= endRange
         && i_darts_counter_temp1228 >= this->inputsTPParent->lastIteration1228;
         i_darts_counter_temp1228--) {
        {
            {
                /*Loop's init*/
                (*j) = (*(*jend));
                int j_darts_counter_temp1228 = (*j);
                for (; j_darts_counter_temp1228 >= (*(*jst)); j_darts_counter_temp1228--) {
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp1228 = (*m);
                        for (; m_darts_counter_temp1228 < 5; m_darts_counter_temp1228++) {
                            tv[(i_darts_counter_temp1228)][j_darts_counter_temp1228]
                              [m_darts_counter_temp1228]
                                = (*(*omega))
                                * (c[(i_darts_counter_temp1228)][j_darts_counter_temp1228]
                                    [m_darts_counter_temp1228][0]
                                        * rsd[(i_darts_counter_temp1228)][j_darts_counter_temp1228]
                                             [(*(*k)) + 1][0]
                                    + c[(i_darts_counter_temp1228)][j_darts_counter_temp1228]
                                       [m_darts_counter_temp1228][1]
                                        * rsd[(i_darts_counter_temp1228)][j_darts_counter_temp1228]
                                             [(*(*k)) + 1][1]
                                    + c[(i_darts_counter_temp1228)][j_darts_counter_temp1228]
                                       [m_darts_counter_temp1228][2]
                                        * rsd[(i_darts_counter_temp1228)][j_darts_counter_temp1228]
                                             [(*(*k)) + 1][2]
                                    + c[(i_darts_counter_temp1228)][j_darts_counter_temp1228]
                                       [m_darts_counter_temp1228][3]
                                        * rsd[(i_darts_counter_temp1228)][j_darts_counter_temp1228]
                                             [(*(*k)) + 1][3]
                                    + c[(i_darts_counter_temp1228)][j_darts_counter_temp1228]
                                       [m_darts_counter_temp1228][4]
                                        * rsd[(i_darts_counter_temp1228)][j_darts_counter_temp1228]
                                             [(*(*k)) + 1][4]);
                        }
                        (*m) = m_darts_counter_temp1228;
                    }
                }
                (*j) = j_darts_counter_temp1228;
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
    myTP->controlTPParent->TPParent->checkInCodelets1330[this->getID()].decDep();
}
TP1228::TP1228(int in_numThreads, int in_mainCodeletID, TP2* in_TPParent, int in_initIteration,
    int in_lastIteration, TP1228** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts1228(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iend_darts1228(new int*[this->numThreads])
    , ist_darts1228(new int*[this->numThreads])
    , j_darts1228(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jend_darts1228(new int*[this->numThreads])
    , jst_darts1228(new int*[this->numThreads])
    , k_darts1228(new int*[this->numThreads])
    , m_darts1228(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , omega_darts1228(new double*[this->numThreads])
    , initIteration1228(in_initIteration)
    , lastIteration1228(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , signalNextReady(new int[baseNumThreads])
    , checkInCodelets1229(new _checkInCodelets1229[this->numThreads])
{
    /*Initialize the loop parameters*/
    range1228 = abs(lastIteration1228 - initIteration1228) / 1;
    rangePerCodelet1228 = range1228 / numThreads;
    minIteration1228 = min<int>(lastIteration1228, initIteration1228);
    remainderRange1228 = range1228 % numThreads;
    /*Initialize inputs and vars.*/
    this->iend_darts1228 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ist_darts1228 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jend_darts1228 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jst_darts1228 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->k_darts1228 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->omega_darts1228
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    _checkInCodelets1229* checkInCodelets1229Ptr = (this->checkInCodelets1229);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets1229);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
        this->signalNextReady[codeletCounter] = 0;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets1229Ptr) = _checkInCodelets1229(2, 1, this, codeletCounter);
#else
        (*checkInCodelets1229Ptr) = _checkInCodelets1229(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets1229Ptr).decDep();
        checkInCodelets1229Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP1228::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets1229[localID].setID(codeletID);
    this->checkInCodelets1229[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets1229[localID + this->baseNumThreads * i]
            = _checkInCodelets1229(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets1229[localID + this->baseNumThreads * i]
            = _checkInCodelets1229(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets1229[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets1229[localID + this->baseNumThreads * i].decDep();
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
TP1228::~TP1228()
{
    delete[] iend_darts1228;
    delete[] ist_darts1228;
    delete[] jend_darts1228;
    delete[] jst_darts1228;
    delete[] k_darts1228;
    delete[] omega_darts1228;
    delete[] checkInCodelets1229;
}
/*TP1330: OMPForDirective*/
bool TP1330::requestNewRangeIterations1330(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 1*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet1330 * codeletID;
        int tempEndRange = rangePerCodelet1330 * (codeletID + 1);
        if (remainderRange1330 != 0) {
            if (codeletID < (uint32_t)remainderRange1330) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange1330;
                tempEndRange += remainderRange1330;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration1330;
        tempEndRange = tempEndRange * 1 + minIteration1330;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration1330 < lastIteration1330) {
            (this->inputsTPParent->i_darts1330[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts1330[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == 0) {
            *endRange = *endRange - 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration1330;
        }
    }
    return isThereNewIteration;
}
void TP1330::_checkInCodelets1331::fire(void)
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
    this->inputsTPParent->iend_darts1330[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iend_darts2[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ist_darts1330[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->ist_darts2[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jend_darts1330[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jend_darts2[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jst_darts1330[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jst_darts2[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->k_darts1330[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->k_darts2[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->omega_darts1330[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->omega_darts2[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts1330[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts2[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp1_darts1330[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp1_darts2[this->getID()]);

    /*printing node 1331: ForStmt*/
    /*var: i*/
    /*var: iend*/
    /*var: ist*/
    /*var: j*/
    /*var: jend*/
    /*var: jst*/
    /*var: k*/
    /*var: m*/
    /*var: omega*/
    /*var: tmp*/
    /*var: tmp1*/
    int* i = &(this->inputsTPParent->i_darts1330[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iend = &(this->inputsTPParent->iend_darts1330[this->getLocalID()]);
    (void)iend /*OMP_SHARED_PRIVATE*/;
    int** ist = &(this->inputsTPParent->ist_darts1330[this->getLocalID()]);
    (void)ist /*OMP_SHARED_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts1330[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jend = &(this->inputsTPParent->jend_darts1330[this->getLocalID()]);
    (void)jend /*OMP_SHARED_PRIVATE*/;
    int** jst = &(this->inputsTPParent->jst_darts1330[this->getLocalID()]);
    (void)jst /*OMP_SHARED_PRIVATE*/;
    int** k = &(this->inputsTPParent->k_darts1330[this->getLocalID()]);
    (void)k /*OMP_SHARED_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts1330[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** omega = &(this->inputsTPParent->omega_darts1330[this->getLocalID()]);
    (void)omega /*OMP_SHARED_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts1330[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** tmp1 = &(this->inputsTPParent->tmp1_darts1330[this->getLocalID()]);
    (void)tmp1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations1330(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Find and signal the next codelet*/
        return;
    }
    for (int i_darts_counter_temp1330 = (*i); i_darts_counter_temp1330 >= endRange
         && i_darts_counter_temp1330 >= this->inputsTPParent->lastIteration1330;
         i_darts_counter_temp1330--) {
        {
            if ((i_darts_counter_temp1330) != (*(*iend))) {
                while (flag[(i_darts_counter_temp1330) + 1] == 0) { }
            }
            if ((i_darts_counter_temp1330) != (*(*ist))) {
                while (flag[(i_darts_counter_temp1330)] == 1) { }
            }
            {
                /*Loop's init*/
                (*j) = (*(*jend));
                int j_darts_counter_temp1330 = (*j);
                for (; j_darts_counter_temp1330 >= (*(*jst)); j_darts_counter_temp1330--) {
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp1330 = (*m);
                        for (; m_darts_counter_temp1330 < 5; m_darts_counter_temp1330++) {
                            tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330]
                              [m_darts_counter_temp1330]
                                = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330]
                                    [m_darts_counter_temp1330]
                                + (*(*omega))
                                    * (b[(i_darts_counter_temp1330)][j_darts_counter_temp1330]
                                        [m_darts_counter_temp1330][0]
                                            * rsd[(i_darts_counter_temp1330)]
                                                 [j_darts_counter_temp1330 + 1][(*(*k))][0]
                                        + a[(i_darts_counter_temp1330)][j_darts_counter_temp1330]
                                           [m_darts_counter_temp1330][0]
                                            * rsd[(i_darts_counter_temp1330) + 1]
                                                 [j_darts_counter_temp1330][(*(*k))][0]
                                        + b[(i_darts_counter_temp1330)][j_darts_counter_temp1330]
                                           [m_darts_counter_temp1330][1]
                                            * rsd[(i_darts_counter_temp1330)]
                                                 [j_darts_counter_temp1330 + 1][(*(*k))][1]
                                        + a[(i_darts_counter_temp1330)][j_darts_counter_temp1330]
                                           [m_darts_counter_temp1330][1]
                                            * rsd[(i_darts_counter_temp1330) + 1]
                                                 [j_darts_counter_temp1330][(*(*k))][1]
                                        + b[(i_darts_counter_temp1330)][j_darts_counter_temp1330]
                                           [m_darts_counter_temp1330][2]
                                            * rsd[(i_darts_counter_temp1330)]
                                                 [j_darts_counter_temp1330 + 1][(*(*k))][2]
                                        + a[(i_darts_counter_temp1330)][j_darts_counter_temp1330]
                                           [m_darts_counter_temp1330][2]
                                            * rsd[(i_darts_counter_temp1330) + 1]
                                                 [j_darts_counter_temp1330][(*(*k))][2]
                                        + b[(i_darts_counter_temp1330)][j_darts_counter_temp1330]
                                           [m_darts_counter_temp1330][3]
                                            * rsd[(i_darts_counter_temp1330)]
                                                 [j_darts_counter_temp1330 + 1][(*(*k))][3]
                                        + a[(i_darts_counter_temp1330)][j_darts_counter_temp1330]
                                           [m_darts_counter_temp1330][3]
                                            * rsd[(i_darts_counter_temp1330) + 1]
                                                 [j_darts_counter_temp1330][(*(*k))][3]
                                        + b[(i_darts_counter_temp1330)][j_darts_counter_temp1330]
                                           [m_darts_counter_temp1330][4]
                                            * rsd[(i_darts_counter_temp1330)]
                                                 [j_darts_counter_temp1330 + 1][(*(*k))][4]
                                        + a[(i_darts_counter_temp1330)][j_darts_counter_temp1330]
                                           [m_darts_counter_temp1330][4]
                                            * rsd[(i_darts_counter_temp1330) + 1]
                                                 [j_darts_counter_temp1330][(*(*k))][4]);
                        }
                        (*m) = m_darts_counter_temp1330;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp1330 = (*m);
                        for (; m_darts_counter_temp1330 < 5; m_darts_counter_temp1330++) {
                            tmat[m_darts_counter_temp1330][0]
                                = d[(i_darts_counter_temp1330)][j_darts_counter_temp1330]
                                   [m_darts_counter_temp1330][0];
                            tmat[m_darts_counter_temp1330][1]
                                = d[(i_darts_counter_temp1330)][j_darts_counter_temp1330]
                                   [m_darts_counter_temp1330][1];
                            tmat[m_darts_counter_temp1330][2]
                                = d[(i_darts_counter_temp1330)][j_darts_counter_temp1330]
                                   [m_darts_counter_temp1330][2];
                            tmat[m_darts_counter_temp1330][3]
                                = d[(i_darts_counter_temp1330)][j_darts_counter_temp1330]
                                   [m_darts_counter_temp1330][3];
                            tmat[m_darts_counter_temp1330][4]
                                = d[(i_darts_counter_temp1330)][j_darts_counter_temp1330]
                                   [m_darts_counter_temp1330][4];
                        }
                        (*m) = m_darts_counter_temp1330;
                    }
                    (*(*tmp1)) = 1. / tmat[0][0];
                    (*(*tmp)) = (*(*tmp1)) * tmat[1][0];
                    tmat[1][1] = tmat[1][1] - (*(*tmp)) * tmat[0][1];
                    tmat[1][2] = tmat[1][2] - (*(*tmp)) * tmat[0][2];
                    tmat[1][3] = tmat[1][3] - (*(*tmp)) * tmat[0][3];
                    tmat[1][4] = tmat[1][4] - (*(*tmp)) * tmat[0][4];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][1]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][1]
                        - tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][0] * (*(*tmp));
                    (*(*tmp)) = (*(*tmp1)) * tmat[2][0];
                    tmat[2][1] = tmat[2][1] - (*(*tmp)) * tmat[0][1];
                    tmat[2][2] = tmat[2][2] - (*(*tmp)) * tmat[0][2];
                    tmat[2][3] = tmat[2][3] - (*(*tmp)) * tmat[0][3];
                    tmat[2][4] = tmat[2][4] - (*(*tmp)) * tmat[0][4];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][2]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][2]
                        - tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][0] * (*(*tmp));
                    (*(*tmp)) = (*(*tmp1)) * tmat[3][0];
                    tmat[3][1] = tmat[3][1] - (*(*tmp)) * tmat[0][1];
                    tmat[3][2] = tmat[3][2] - (*(*tmp)) * tmat[0][2];
                    tmat[3][3] = tmat[3][3] - (*(*tmp)) * tmat[0][3];
                    tmat[3][4] = tmat[3][4] - (*(*tmp)) * tmat[0][4];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][3]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][3]
                        - tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][0] * (*(*tmp));
                    (*(*tmp)) = (*(*tmp1)) * tmat[4][0];
                    tmat[4][1] = tmat[4][1] - (*(*tmp)) * tmat[0][1];
                    tmat[4][2] = tmat[4][2] - (*(*tmp)) * tmat[0][2];
                    tmat[4][3] = tmat[4][3] - (*(*tmp)) * tmat[0][3];
                    tmat[4][4] = tmat[4][4] - (*(*tmp)) * tmat[0][4];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][4]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][4]
                        - tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][0] * (*(*tmp));
                    (*(*tmp1)) = 1. / tmat[1][1];
                    (*(*tmp)) = (*(*tmp1)) * tmat[2][1];
                    tmat[2][2] = tmat[2][2] - (*(*tmp)) * tmat[1][2];
                    tmat[2][3] = tmat[2][3] - (*(*tmp)) * tmat[1][3];
                    tmat[2][4] = tmat[2][4] - (*(*tmp)) * tmat[1][4];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][2]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][2]
                        - tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][1] * (*(*tmp));
                    (*(*tmp)) = (*(*tmp1)) * tmat[3][1];
                    tmat[3][2] = tmat[3][2] - (*(*tmp)) * tmat[1][2];
                    tmat[3][3] = tmat[3][3] - (*(*tmp)) * tmat[1][3];
                    tmat[3][4] = tmat[3][4] - (*(*tmp)) * tmat[1][4];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][3]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][3]
                        - tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][1] * (*(*tmp));
                    (*(*tmp)) = (*(*tmp1)) * tmat[4][1];
                    tmat[4][2] = tmat[4][2] - (*(*tmp)) * tmat[1][2];
                    tmat[4][3] = tmat[4][3] - (*(*tmp)) * tmat[1][3];
                    tmat[4][4] = tmat[4][4] - (*(*tmp)) * tmat[1][4];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][4]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][4]
                        - tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][1] * (*(*tmp));
                    (*(*tmp1)) = 1. / tmat[2][2];
                    (*(*tmp)) = (*(*tmp1)) * tmat[3][2];
                    tmat[3][3] = tmat[3][3] - (*(*tmp)) * tmat[2][3];
                    tmat[3][4] = tmat[3][4] - (*(*tmp)) * tmat[2][4];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][3]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][3]
                        - tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][2] * (*(*tmp));
                    (*(*tmp)) = (*(*tmp1)) * tmat[4][2];
                    tmat[4][3] = tmat[4][3] - (*(*tmp)) * tmat[2][3];
                    tmat[4][4] = tmat[4][4] - (*(*tmp)) * tmat[2][4];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][4]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][4]
                        - tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][2] * (*(*tmp));
                    (*(*tmp1)) = 1. / tmat[3][3];
                    (*(*tmp)) = (*(*tmp1)) * tmat[4][3];
                    tmat[4][4] = tmat[4][4] - (*(*tmp)) * tmat[3][4];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][4]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][4]
                        - tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][3] * (*(*tmp));
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][4]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][4] / tmat[4][4];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][3]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][3]
                        - tmat[3][4] * tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][4];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][3]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][3] / tmat[3][3];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][2]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][2]
                        - tmat[2][3] * tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][3]
                        - tmat[2][4] * tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][4];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][2]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][2] / tmat[2][2];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][1]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][1]
                        - tmat[1][2] * tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][2]
                        - tmat[1][3] * tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][3]
                        - tmat[1][4] * tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][4];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][1]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][1] / tmat[1][1];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][0]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][0]
                        - tmat[0][1] * tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][1]
                        - tmat[0][2] * tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][2]
                        - tmat[0][3] * tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][3]
                        - tmat[0][4] * tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][4];
                    tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][0]
                        = tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][0] / tmat[0][0];
                    rsd[(i_darts_counter_temp1330)][j_darts_counter_temp1330][(*(*k))][0]
                        = rsd[(i_darts_counter_temp1330)][j_darts_counter_temp1330][(*(*k))][0]
                        - tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][0];
                    rsd[(i_darts_counter_temp1330)][j_darts_counter_temp1330][(*(*k))][1]
                        = rsd[(i_darts_counter_temp1330)][j_darts_counter_temp1330][(*(*k))][1]
                        - tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][1];
                    rsd[(i_darts_counter_temp1330)][j_darts_counter_temp1330][(*(*k))][2]
                        = rsd[(i_darts_counter_temp1330)][j_darts_counter_temp1330][(*(*k))][2]
                        - tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][2];
                    rsd[(i_darts_counter_temp1330)][j_darts_counter_temp1330][(*(*k))][3]
                        = rsd[(i_darts_counter_temp1330)][j_darts_counter_temp1330][(*(*k))][3]
                        - tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][3];
                    rsd[(i_darts_counter_temp1330)][j_darts_counter_temp1330][(*(*k))][4]
                        = rsd[(i_darts_counter_temp1330)][j_darts_counter_temp1330][(*(*k))][4]
                        - tv[(i_darts_counter_temp1330)][j_darts_counter_temp1330][4];
                }
                (*j) = j_darts_counter_temp1330;
            }
            if ((i_darts_counter_temp1330) != (*(*iend))) {
                flag[(i_darts_counter_temp1330) + 1] = 0;
            }
            if ((i_darts_counter_temp1330) != (*(*ist))) {
                flag[(i_darts_counter_temp1330)] = 1;
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
TP1330::TP1330(int in_numThreads, int in_mainCodeletID, TP2* in_TPParent, int in_initIteration,
    int in_lastIteration, TP1330** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts1330(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iend_darts1330(new int*[this->numThreads])
    , ist_darts1330(new int*[this->numThreads])
    , j_darts1330(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jend_darts1330(new int*[this->numThreads])
    , jst_darts1330(new int*[this->numThreads])
    , k_darts1330(new int*[this->numThreads])
    , m_darts1330(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , omega_darts1330(new double*[this->numThreads])
    , tmp_darts1330(new double*[this->numThreads])
    , tmp1_darts1330(new double*[this->numThreads])
    , initIteration1330(in_initIteration)
    , lastIteration1330(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , signalNextReady(new int[baseNumThreads])
    , checkInCodelets1331(new _checkInCodelets1331[this->numThreads])
{
    /*Initialize the loop parameters*/
    range1330 = abs(lastIteration1330 - initIteration1330) / 1;
    rangePerCodelet1330 = range1330 / numThreads;
    minIteration1330 = min<int>(lastIteration1330, initIteration1330);
    remainderRange1330 = range1330 % numThreads;
    /*Initialize inputs and vars.*/
    this->iend_darts1330 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ist_darts1330 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jend_darts1330 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jst_darts1330 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->k_darts1330 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->omega_darts1330
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts1330
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp1_darts1330
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    _checkInCodelets1331* checkInCodelets1331Ptr = (this->checkInCodelets1331);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets1331);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
        this->signalNextReady[codeletCounter] = 0;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets1331Ptr) = _checkInCodelets1331(2, 1, this, codeletCounter);
#else
        (*checkInCodelets1331Ptr) = _checkInCodelets1331(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets1331Ptr).decDep();
        checkInCodelets1331Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP1330::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets1331[localID].setID(codeletID);
    this->checkInCodelets1331[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets1331[localID + this->baseNumThreads * i]
            = _checkInCodelets1331(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets1331[localID + this->baseNumThreads * i]
            = _checkInCodelets1331(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets1331[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets1331[localID + this->baseNumThreads * i].decDep();
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
TP1330::~TP1330()
{
    delete[] iend_darts1330;
    delete[] ist_darts1330;
    delete[] jend_darts1330;
    delete[] jst_darts1330;
    delete[] k_darts1330;
    delete[] omega_darts1330;
    delete[] tmp_darts1330;
    delete[] tmp1_darts1330;
    delete[] checkInCodelets1331;
}
/*TP2241: OMPParallelDirective*/
void TP2241::_barrierCodelets2241::fire(void)
{
    TP2241* myTP = static_cast<TP2241*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP2241::_checkInCodelets2243::fire(void)
{
    /*Init the vars for this region*/

    /*printing node 2243: DeclStmt*/

    /*printing node 2244: DeclStmt*/

    /*printing node 2245: DeclStmt*/

    /*printing node 2246: DeclStmt*/

    /*printing node 2247: DeclStmt*/

    /*printing node 2248: DeclStmt*/

    /*printing node 2249: DeclStmt*/

    /*printing node 2250: DeclStmt*/

    /*printing node 2251: DeclStmt*/

    /*printing node 2252: DeclStmt*/

    /*printing node 2253: DeclStmt*/

    /*printing node 2254: DeclStmt*/

    /*printing node 2255: DeclStmt*/

    /*printing node 2256: DeclStmt*/

    /*printing node 2257: DeclStmt*/

    /*printing node 2258: DeclStmt*/

    /*printing node 2259: BinaryOperator*/
    (this->inputsTPParent->dsspm_darts2241[this->getID()]) = dssp;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 2243 nextRegion: 2260 */
    myTP->controlTPParent->checkInCodelets2260[this->getID()].decDep();
}
void TP2241::_checkInCodelets2260::fire(void)
{
    /*region 2260 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP2260;
    if (idx < myTP->TPsToUse2260) {
        if (!__sync_val_compare_and_swap(&(myTP->TP2260_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(nx - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse2260;
            int minIteration = min<int>(nx, 0);
            int remainderRange = range % myTP->TPsToUse2260;
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
            if (idx == myTP->TPsToUse2260 - 1) {
                lastIteration = nx;
            }
#if USEINVOKE == 1
            invoke<TP2260>(myTP, myTP->codeletsPerTP2260 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP2260Ptr[idx]));
#else
            place<TP2260>(idx, myTP, myTP->codeletsPerTP2260 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP2260Ptr[idx]));
#endif
        } else {
            if (myTP->TP2260Ptr[idx] != nullptr) {
                myTP->TP2260Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2241::_barrierCodelets2260::fire(void)
{
    TP2241* myTP = static_cast<TP2241*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets2311[codeletsCounter].decDep();
        }
    }
}
void TP2241::_checkInCodelets2311::fire(void)
{
    /*region 2311 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP2311;
    if (idx < myTP->TPsToUse2311) {
        if (!__sync_val_compare_and_swap(&(myTP->TP2311_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(nx - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse2311;
            int minIteration = min<int>(nx, 0);
            int remainderRange = range % myTP->TPsToUse2311;
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
            if (idx == myTP->TPsToUse2311 - 1) {
                lastIteration = nx;
            }
#if USEINVOKE == 1
            invoke<TP2311>(myTP, myTP->codeletsPerTP2311 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP2311Ptr[idx]));
#else
            place<TP2311>(idx, myTP, myTP->codeletsPerTP2311 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP2311Ptr[idx]));
#endif
        } else {
            if (myTP->TP2311Ptr[idx] != nullptr) {
                myTP->TP2311Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2241::_barrierCodelets2311::fire(void)
{
    TP2241* myTP = static_cast<TP2241*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets2443[codeletsCounter].decDep();
        }
    }
}
void TP2241::_checkInCodelets2443::fire(void)
{

    /*printing node 2443: BinaryOperator*/
    (this->inputsTPParent->L1_darts2241[this->getID()]) = 0;

    /*printing node 2444: BinaryOperator*/
    (this->inputsTPParent->L2_darts2241[this->getID()]) = nx - 1;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 2443 nextRegion: 2446 */
    myTP->controlTPParent->checkInCodelets2446[this->getID()].decDep();
}
void TP2241::_checkInCodelets2446::fire(void)
{
    /*region 2446 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP2446;
    if (idx < myTP->TPsToUse2446) {
        if (!__sync_val_compare_and_swap(&(myTP->TP2446_alreadyLaunched[idx]), 0, 1)) {
            int range = abs((this->inputsTPParent->L2_darts2241[this->getID()])
                            - (this->inputsTPParent->L1_darts2241[this->getID()]))
                / 1;
            int rangePerCodelet = range / myTP->TPsToUse2446;
            int minIteration = min<int>((this->inputsTPParent->L2_darts2241[this->getID()]),
                (this->inputsTPParent->L1_darts2241[this->getID()]));
            int remainderRange = range % myTP->TPsToUse2446;
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
            if ((this->inputsTPParent->L1_darts2241[this->getID()])
                < (this->inputsTPParent->L2_darts2241[this->getID()])) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse2446 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse2446 - 1) {
                lastIteration = (this->inputsTPParent->L2_darts2241[this->getID()]);
            }
#if USEINVOKE == 1
            invoke<TP2446>(myTP, myTP->codeletsPerTP2446 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP2446Ptr[idx]));
#else
            place<TP2446>(idx, myTP, myTP->codeletsPerTP2446 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP2446Ptr[idx]));
#endif
        } else {
            if (myTP->TP2446Ptr[idx] != nullptr) {
                myTP->TP2446Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2241::_barrierCodelets2446::fire(void)
{
    TP2241* myTP = static_cast<TP2241*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets2595[codeletsCounter].decDep();
        }
    }
}
void TP2241::_checkInCodelets2595::fire(void)
{
    /*region 2595 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP2595;
    if (idx < myTP->TPsToUse2595) {
        if (!__sync_val_compare_and_swap(&(myTP->TP2595_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(jend - jst) / 1;
            int rangePerCodelet = range / myTP->TPsToUse2595;
            int minIteration = min<int>(jend, jst);
            int remainderRange = range % myTP->TPsToUse2595;
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
            if (idx == myTP->TPsToUse2595 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse2595 - 1) {
                lastIteration = jend;
            }
#if USEINVOKE == 1
            invoke<TP2595>(myTP, myTP->codeletsPerTP2595 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP2595Ptr[idx]));
#else
            place<TP2595>(idx, myTP, myTP->codeletsPerTP2595 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP2595Ptr[idx]));
#endif
        } else {
            if (myTP->TP2595Ptr[idx] != nullptr) {
                myTP->TP2595Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2241::_barrierCodelets2595::fire(void)
{
    TP2241* myTP = static_cast<TP2241*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets3230[codeletsCounter].decDep();
        }
    }
}
void TP2241::_checkInCodelets3230::fire(void)
{

    /*printing node 3230: BinaryOperator*/
    (this->inputsTPParent->L1_darts2241[this->getID()]) = 0;

    /*printing node 3231: BinaryOperator*/
    (this->inputsTPParent->L2_darts2241[this->getID()]) = ny - 1;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 3230 nextRegion: 3233 */
    myTP->controlTPParent->checkInCodelets3233[this->getID()].decDep();
}
void TP2241::_checkInCodelets3233::fire(void)
{
    /*region 3233 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP3233;
    if (idx < myTP->TPsToUse3233) {
        if (!__sync_val_compare_and_swap(&(myTP->TP3233_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse3233;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse3233;
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
            if (idx == myTP->TPsToUse3233 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse3233 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP3233>(myTP, myTP->codeletsPerTP3233 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP3233Ptr[idx]));
#else
            place<TP3233>(idx, myTP, myTP->codeletsPerTP3233 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP3233Ptr[idx]));
#endif
        } else {
            if (myTP->TP3233Ptr[idx] != nullptr) {
                myTP->TP3233Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2241::_barrierCodelets3233::fire(void)
{
    TP2241* myTP = static_cast<TP2241*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets3382[codeletsCounter].decDep();
        }
    }
}
void TP2241::_checkInCodelets3382::fire(void)
{
    /*region 3382 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP3382;
    if (idx < myTP->TPsToUse3382) {
        if (!__sync_val_compare_and_swap(&(myTP->TP3382_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse3382;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse3382;
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
            if (idx == myTP->TPsToUse3382 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse3382 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP3382>(myTP, myTP->codeletsPerTP3382 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP3382Ptr[idx]));
#else
            place<TP3382>(idx, myTP, myTP->codeletsPerTP3382 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP3382Ptr[idx]));
#endif
        } else {
            if (myTP->TP3382Ptr[idx] != nullptr) {
                myTP->TP3382Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2241::_barrierCodelets3382::fire(void)
{
    TP2241* myTP = static_cast<TP2241*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets4017[codeletsCounter].decDep();
        }
    }
}
void TP2241::_checkInCodelets4017::fire(void)
{
    /*region 4017 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP4017;
    if (idx < myTP->TPsToUse4017) {
        if (!__sync_val_compare_and_swap(&(myTP->TP4017_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse4017;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse4017;
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
            if (idx == myTP->TPsToUse4017 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse4017 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP4017>(myTP, myTP->codeletsPerTP4017 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP4017Ptr[idx]));
#else
            place<TP4017>(idx, myTP, myTP->codeletsPerTP4017 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP4017Ptr[idx]));
#endif
        } else {
            if (myTP->TP4017Ptr[idx] != nullptr) {
                myTP->TP4017Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP2241::_barrierCodelets4017::fire(void)
{
    TP2241* myTP = static_cast<TP2241*>(myTP_);
    myTP->TPParent->barrierCodelets2241[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets2241[0]));
}
TP2241::TP2241(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts2241(new int[this->numThreads]) /*VARIABLE*/
    , L2_darts2241(new int[this->numThreads]) /*VARIABLE*/
    , dsspm_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , eta_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , i_darts2241(new int[this->numThreads]) /*VARIABLE*/
    , iend1_darts2241(new int[this->numThreads]) /*VARIABLE*/
    , iglob_darts2241(new int[this->numThreads]) /*VARIABLE*/
    , ist1_darts2241(new int[this->numThreads]) /*VARIABLE*/
    , j_darts2241(new int[this->numThreads]) /*VARIABLE*/
    , jend1_darts2241(new int[this->numThreads]) /*VARIABLE*/
    , jglob_darts2241(new int[this->numThreads]) /*VARIABLE*/
    , jst1_darts2241(new int[this->numThreads]) /*VARIABLE*/
    , k_darts2241(new int[this->numThreads]) /*VARIABLE*/
    , m_darts2241(new int[this->numThreads]) /*VARIABLE*/
    , q_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , tmp_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u21_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u21i_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u21im1_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u21j_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u21jm1_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u21k_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u21km1_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u31_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u31i_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u31im1_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u31j_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u31jm1_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u31k_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u31km1_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u41_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u41i_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u41im1_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u41j_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u41jm1_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u41k_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u41km1_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u51i_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u51im1_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u51j_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u51jm1_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u51k_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , u51km1_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , xi_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , zeta_darts2241(new double[this->numThreads]) /*VARIABLE*/
    , TP2260Ptr(new TP2260*[NUMTPS2260])
    , TP2260_alreadyLaunched(new size_t[NUMTPS2260])
    , numTPsSet2260(0)
    , numTPsReady2260(0)
    , TPsToUse2260(NUMTPS2260)
    , codeletsPerTP2260(this->numThreads / NUMTPS2260)
    , totalCodelets2260(this->TPsToUse2260 * this->codeletsPerTP2260)
    , TP2311Ptr(new TP2311*[NUMTPS2311])
    , TP2311_alreadyLaunched(new size_t[NUMTPS2311])
    , numTPsSet2311(0)
    , numTPsReady2311(0)
    , TPsToUse2311(NUMTPS2311)
    , codeletsPerTP2311(this->numThreads / NUMTPS2311)
    , totalCodelets2311(this->TPsToUse2311 * this->codeletsPerTP2311)
    , TP2446Ptr(new TP2446*[NUMTPS2446])
    , TP2446_alreadyLaunched(new size_t[NUMTPS2446])
    , numTPsSet2446(0)
    , numTPsReady2446(0)
    , TPsToUse2446(NUMTPS2446)
    , codeletsPerTP2446(this->numThreads / NUMTPS2446)
    , totalCodelets2446(this->TPsToUse2446 * this->codeletsPerTP2446)
    , TP2595Ptr(new TP2595*[NUMTPS2595])
    , TP2595_alreadyLaunched(new size_t[NUMTPS2595])
    , numTPsSet2595(0)
    , numTPsReady2595(0)
    , TPsToUse2595(NUMTPS2595)
    , codeletsPerTP2595(this->numThreads / NUMTPS2595)
    , totalCodelets2595(this->TPsToUse2595 * this->codeletsPerTP2595)
    , TP3233Ptr(new TP3233*[NUMTPS3233])
    , TP3233_alreadyLaunched(new size_t[NUMTPS3233])
    , numTPsSet3233(0)
    , numTPsReady3233(0)
    , TPsToUse3233(NUMTPS3233)
    , codeletsPerTP3233(this->numThreads / NUMTPS3233)
    , totalCodelets3233(this->TPsToUse3233 * this->codeletsPerTP3233)
    , TP3382Ptr(new TP3382*[NUMTPS3382])
    , TP3382_alreadyLaunched(new size_t[NUMTPS3382])
    , numTPsSet3382(0)
    , numTPsReady3382(0)
    , TPsToUse3382(NUMTPS3382)
    , codeletsPerTP3382(this->numThreads / NUMTPS3382)
    , totalCodelets3382(this->TPsToUse3382 * this->codeletsPerTP3382)
    , TP4017Ptr(new TP4017*[NUMTPS4017])
    , TP4017_alreadyLaunched(new size_t[NUMTPS4017])
    , numTPsSet4017(0)
    , numTPsReady4017(0)
    , TPsToUse4017(NUMTPS4017)
    , codeletsPerTP4017(this->numThreads / NUMTPS4017)
    , totalCodelets4017(this->TPsToUse4017 * this->codeletsPerTP4017)
    , barrierCodelets2241(new _barrierCodelets2241[1])
    , checkInCodelets2243(new _checkInCodelets2243[this->numThreads])
    , checkInCodelets2260(new _checkInCodelets2260[this->numThreads])
    , barrierCodelets2260(new _barrierCodelets2260[1])
    , checkInCodelets2311(new _checkInCodelets2311[this->numThreads])
    , barrierCodelets2311(new _barrierCodelets2311[1])
    , checkInCodelets2443(new _checkInCodelets2443[this->numThreads])
    , checkInCodelets2446(new _checkInCodelets2446[this->numThreads])
    , barrierCodelets2446(new _barrierCodelets2446[1])
    , checkInCodelets2595(new _checkInCodelets2595[this->numThreads])
    , barrierCodelets2595(new _barrierCodelets2595[1])
    , checkInCodelets3230(new _checkInCodelets3230[this->numThreads])
    , checkInCodelets3233(new _checkInCodelets3233[this->numThreads])
    , barrierCodelets3233(new _barrierCodelets3233[1])
    , checkInCodelets3382(new _checkInCodelets3382[this->numThreads])
    , barrierCodelets3382(new _barrierCodelets3382[1])
    , checkInCodelets4017(new _checkInCodelets4017[this->numThreads])
    , barrierCodelets4017(new _barrierCodelets4017[1])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets2241[0] = _barrierCodelets2241(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets4017[0] = _barrierCodelets4017(NUMTPS4017, NUMTPS4017, this, 0);
    barrierCodelets3382[0] = _barrierCodelets3382(NUMTPS3382, NUMTPS3382, this, 0);
    barrierCodelets3233[0] = _barrierCodelets3233(NUMTPS3233, NUMTPS3233, this, 0);
    barrierCodelets2595[0] = _barrierCodelets2595(NUMTPS2595, NUMTPS2595, this, 0);
    barrierCodelets2446[0] = _barrierCodelets2446(NUMTPS2446, NUMTPS2446, this, 0);
    barrierCodelets2311[0] = _barrierCodelets2311(NUMTPS2311, NUMTPS2311, this, 0);
    barrierCodelets2260[0] = _barrierCodelets2260(NUMTPS2260, NUMTPS2260, this, 0);
    _checkInCodelets4017* checkInCodelets4017Ptr = (this->checkInCodelets4017);
    for (int i = 0; i < NUMTPS4017; i++) {
        TP4017Ptr[i] = nullptr;
        TP4017_alreadyLaunched[i] = 0;
    }
    _checkInCodelets3382* checkInCodelets3382Ptr = (this->checkInCodelets3382);
    for (int i = 0; i < NUMTPS3382; i++) {
        TP3382Ptr[i] = nullptr;
        TP3382_alreadyLaunched[i] = 0;
    }
    _checkInCodelets3233* checkInCodelets3233Ptr = (this->checkInCodelets3233);
    for (int i = 0; i < NUMTPS3233; i++) {
        TP3233Ptr[i] = nullptr;
        TP3233_alreadyLaunched[i] = 0;
    }
    _checkInCodelets3230* checkInCodelets3230Ptr = (this->checkInCodelets3230);
    _checkInCodelets2595* checkInCodelets2595Ptr = (this->checkInCodelets2595);
    for (int i = 0; i < NUMTPS2595; i++) {
        TP2595Ptr[i] = nullptr;
        TP2595_alreadyLaunched[i] = 0;
    }
    _checkInCodelets2446* checkInCodelets2446Ptr = (this->checkInCodelets2446);
    for (int i = 0; i < NUMTPS2446; i++) {
        TP2446Ptr[i] = nullptr;
        TP2446_alreadyLaunched[i] = 0;
    }
    _checkInCodelets2443* checkInCodelets2443Ptr = (this->checkInCodelets2443);
    _checkInCodelets2311* checkInCodelets2311Ptr = (this->checkInCodelets2311);
    for (int i = 0; i < NUMTPS2311; i++) {
        TP2311Ptr[i] = nullptr;
        TP2311_alreadyLaunched[i] = 0;
    }
    _checkInCodelets2260* checkInCodelets2260Ptr = (this->checkInCodelets2260);
    for (int i = 0; i < NUMTPS2260; i++) {
        TP2260Ptr[i] = nullptr;
        TP2260_alreadyLaunched[i] = 0;
    }
    _checkInCodelets2243* checkInCodelets2243Ptr = (this->checkInCodelets2243);
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets4017Ptr) = _checkInCodelets4017(1, 1, this, codeletCounter);
        checkInCodelets4017Ptr++;
        (*checkInCodelets3382Ptr) = _checkInCodelets3382(1, 1, this, codeletCounter);
        checkInCodelets3382Ptr++;
        (*checkInCodelets3233Ptr) = _checkInCodelets3233(1, 1, this, codeletCounter);
        checkInCodelets3233Ptr++;
        (*checkInCodelets3230Ptr) = _checkInCodelets3230(1, 1, this, codeletCounter);
        checkInCodelets3230Ptr++;
        (*checkInCodelets2595Ptr) = _checkInCodelets2595(1, 1, this, codeletCounter);
        checkInCodelets2595Ptr++;
        (*checkInCodelets2446Ptr) = _checkInCodelets2446(1, 1, this, codeletCounter);
        checkInCodelets2446Ptr++;
        (*checkInCodelets2443Ptr) = _checkInCodelets2443(1, 1, this, codeletCounter);
        checkInCodelets2443Ptr++;
        (*checkInCodelets2311Ptr) = _checkInCodelets2311(1, 1, this, codeletCounter);
        checkInCodelets2311Ptr++;
        (*checkInCodelets2260Ptr) = _checkInCodelets2260(1, 1, this, codeletCounter);
        checkInCodelets2260Ptr++;
        (*checkInCodelets2243Ptr) = _checkInCodelets2243(1, 1, this, codeletCounter);
        (*checkInCodelets2243Ptr).decDep();
        checkInCodelets2243Ptr++;
    }
}
TP2241::~TP2241()
{
    delete[] L1_darts2241;
    delete[] L2_darts2241;
    delete[] dsspm_darts2241;
    delete[] eta_darts2241;
    delete[] i_darts2241;
    delete[] iend1_darts2241;
    delete[] iglob_darts2241;
    delete[] ist1_darts2241;
    delete[] j_darts2241;
    delete[] jend1_darts2241;
    delete[] jglob_darts2241;
    delete[] jst1_darts2241;
    delete[] k_darts2241;
    delete[] m_darts2241;
    delete[] q_darts2241;
    delete[] tmp_darts2241;
    delete[] u21_darts2241;
    delete[] u21i_darts2241;
    delete[] u21im1_darts2241;
    delete[] u21j_darts2241;
    delete[] u21jm1_darts2241;
    delete[] u21k_darts2241;
    delete[] u21km1_darts2241;
    delete[] u31_darts2241;
    delete[] u31i_darts2241;
    delete[] u31im1_darts2241;
    delete[] u31j_darts2241;
    delete[] u31jm1_darts2241;
    delete[] u31k_darts2241;
    delete[] u31km1_darts2241;
    delete[] u41_darts2241;
    delete[] u41i_darts2241;
    delete[] u41im1_darts2241;
    delete[] u41j_darts2241;
    delete[] u41jm1_darts2241;
    delete[] u41k_darts2241;
    delete[] u41km1_darts2241;
    delete[] u51i_darts2241;
    delete[] u51im1_darts2241;
    delete[] u51j_darts2241;
    delete[] u51jm1_darts2241;
    delete[] u51k_darts2241;
    delete[] u51km1_darts2241;
    delete[] xi_darts2241;
    delete[] zeta_darts2241;
    delete[] barrierCodelets2241;
    delete[] barrierCodelets4017;
    delete[] checkInCodelets4017;
    delete[] barrierCodelets3382;
    delete[] checkInCodelets3382;
    delete[] barrierCodelets3233;
    delete[] checkInCodelets3233;
    delete[] checkInCodelets3230;
    delete[] barrierCodelets2595;
    delete[] checkInCodelets2595;
    delete[] barrierCodelets2446;
    delete[] checkInCodelets2446;
    delete[] checkInCodelets2443;
    delete[] barrierCodelets2311;
    delete[] checkInCodelets2311;
    delete[] barrierCodelets2260;
    delete[] checkInCodelets2260;
    delete[] checkInCodelets2243;
}
/*TP2260: OMPForDirective*/
void TP2260::_barrierCodelets2260::fire(void)
{
    TP2260* myTP = static_cast<TP2260*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets2260[0].decDep();
}
bool TP2260::requestNewRangeIterations2260(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet2260 * codeletID;
        int tempEndRange = rangePerCodelet2260 * (codeletID + 1);
        if (remainderRange2260 != 0) {
            if (codeletID < (uint32_t)remainderRange2260) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange2260;
                tempEndRange += remainderRange2260;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration2260;
        tempEndRange = tempEndRange * 1 + minIteration2260;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration2260 < lastIteration2260) {
            (this->inputsTPParent->i_darts2260[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts2260[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration2260;
        }
    }
    return isThereNewIteration;
}
void TP2260::_checkInCodelets2261::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 2261: ForStmt*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    int* i = &(this->inputsTPParent->i_darts2260[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts2260[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts2260[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts2260[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations2260(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets2260[0].decDep();
        return;
    }
    for (int i_darts_counter_temp2260 = (*i); i_darts_counter_temp2260 < endRange
         && i_darts_counter_temp2260 < this->inputsTPParent->lastIteration2260;
         i_darts_counter_temp2260++) {
        {
            {
                /*Loop's init*/
                (*j) = 0;
                int j_darts_counter_temp2260 = (*j);
                for (; j_darts_counter_temp2260 < ny; j_darts_counter_temp2260++) {
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp2260 = (*k);
                        for (; k_darts_counter_temp2260 < nz; k_darts_counter_temp2260++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp2260 = (*m);
                                for (; m_darts_counter_temp2260 < 5; m_darts_counter_temp2260++) {
                                    frct[(i_darts_counter_temp2260)][j_darts_counter_temp2260]
                                        [k_darts_counter_temp2260][m_darts_counter_temp2260]
                                        = 0.;
                                }
                                (*m) = m_darts_counter_temp2260;
                            }
                        }
                        (*k) = k_darts_counter_temp2260;
                    }
                }
                (*j) = j_darts_counter_temp2260;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets2260[0].decDep();
}
TP2260::TP2260(int in_numThreads, int in_mainCodeletID, TP2241* in_TPParent, int in_initIteration,
    int in_lastIteration, TP2260** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts2260(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts2260(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts2260(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts2260(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration2260(in_initIteration)
    , lastIteration2260(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets2260(new _barrierCodelets2260[1])
    , checkInCodelets2261(new _checkInCodelets2261[this->numThreads])
{
    /*Initialize the loop parameters*/
    range2260 = abs(lastIteration2260 - initIteration2260) / 1;
    rangePerCodelet2260 = range2260 / numThreads;
    minIteration2260 = min<int>(lastIteration2260, initIteration2260);
    remainderRange2260 = range2260 % numThreads;
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets2260[0] = _barrierCodelets2260(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets2261* checkInCodelets2261Ptr = (this->checkInCodelets2261);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets2261);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets2261Ptr) = _checkInCodelets2261(2, 1, this, codeletCounter);
#else
        (*checkInCodelets2261Ptr) = _checkInCodelets2261(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets2261Ptr).decDep();
        checkInCodelets2261Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP2260::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets2261[localID].setID(codeletID);
    this->checkInCodelets2261[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets2261[localID + this->baseNumThreads * i]
            = _checkInCodelets2261(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets2261[localID + this->baseNumThreads * i]
            = _checkInCodelets2261(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets2261[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets2261[localID + this->baseNumThreads * i].decDep();
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
TP2260::~TP2260()
{
    delete[] barrierCodelets2260;
    delete[] checkInCodelets2261;
}
/*TP2311: OMPForDirective*/
void TP2311::_barrierCodelets2311::fire(void)
{
    TP2311* myTP = static_cast<TP2311*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets2311[0].decDep();
}
bool TP2311::requestNewRangeIterations2311(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet2311 * codeletID;
        int tempEndRange = rangePerCodelet2311 * (codeletID + 1);
        if (remainderRange2311 != 0) {
            if (codeletID < (uint32_t)remainderRange2311) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange2311;
                tempEndRange += remainderRange2311;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration2311;
        tempEndRange = tempEndRange * 1 + minIteration2311;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration2311 < lastIteration2311) {
            (this->inputsTPParent->i_darts2311[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts2311[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration2311;
        }
    }
    return isThereNewIteration;
}
void TP2311::_checkInCodelets2312::fire(void)
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
    this->inputsTPParent->eta_darts2311[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->eta_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->iglob_darts2311[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jglob_darts2311[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->xi_darts2311[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->xi_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->zeta_darts2311[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->zeta_darts2241[this->getID()]);

    /*printing node 2312: ForStmt*/
    /*var: eta*/
    /*var: i*/
    /*var: iglob*/
    /*var: j*/
    /*var: jglob*/
    /*var: k*/
    /*var: m*/
    /*var: xi*/
    /*var: zeta*/
    double** eta = &(this->inputsTPParent->eta_darts2311[this->getLocalID()]);
    (void)eta /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts2311[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iglob = &(this->inputsTPParent->iglob_darts2311[this->getLocalID()]);
    (void)iglob /*OMP_SHARED_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts2311[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jglob = &(this->inputsTPParent->jglob_darts2311[this->getLocalID()]);
    (void)jglob /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts2311[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts2311[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** xi = &(this->inputsTPParent->xi_darts2311[this->getLocalID()]);
    (void)xi /*OMP_SHARED_PRIVATE*/;
    double** zeta = &(this->inputsTPParent->zeta_darts2311[this->getLocalID()]);
    (void)zeta /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations2311(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets2311[0].decDep();
        return;
    }
    for (int i_darts_counter_temp2311 = (*i); i_darts_counter_temp2311 < endRange
         && i_darts_counter_temp2311 < this->inputsTPParent->lastIteration2311;
         i_darts_counter_temp2311++) {
        {
            (*(*iglob)) = (i_darts_counter_temp2311);
            (*(*xi)) = ((double)((*(*iglob)))) / (nx0 - 1);
            {
                /*Loop's init*/
                (*j) = 0;
                int j_darts_counter_temp2311 = (*j);
                for (; j_darts_counter_temp2311 < ny; j_darts_counter_temp2311++) {
                    (*(*jglob)) = j_darts_counter_temp2311;
                    (*(*eta)) = ((double)((*(*jglob)))) / (ny0 - 1);
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp2311 = (*k);
                        for (; k_darts_counter_temp2311 < nz; k_darts_counter_temp2311++) {
                            (*(*zeta)) = ((double)(k_darts_counter_temp2311)) / (nz - 1);
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp2311 = (*m);
                                for (; m_darts_counter_temp2311 < 5; m_darts_counter_temp2311++) {
                                    rsd[(i_darts_counter_temp2311)][j_darts_counter_temp2311]
                                       [k_darts_counter_temp2311][m_darts_counter_temp2311]
                                        = ce[m_darts_counter_temp2311][0]
                                        + ce[m_darts_counter_temp2311][1] * (*(*xi))
                                        + ce[m_darts_counter_temp2311][2] * (*(*eta))
                                        + ce[m_darts_counter_temp2311][3] * (*(*zeta))
                                        + ce[m_darts_counter_temp2311][4] * (*(*xi)) * (*(*xi))
                                        + ce[m_darts_counter_temp2311][5] * (*(*eta)) * (*(*eta))
                                        + ce[m_darts_counter_temp2311][6] * (*(*zeta)) * (*(*zeta))
                                        + ce[m_darts_counter_temp2311][7] * (*(*xi)) * (*(*xi))
                                            * (*(*xi))
                                        + ce[m_darts_counter_temp2311][8] * (*(*eta)) * (*(*eta))
                                            * (*(*eta))
                                        + ce[m_darts_counter_temp2311][9] * (*(*zeta)) * (*(*zeta))
                                            * (*(*zeta))
                                        + ce[m_darts_counter_temp2311][10] * (*(*xi)) * (*(*xi))
                                            * (*(*xi)) * (*(*xi))
                                        + ce[m_darts_counter_temp2311][11] * (*(*eta)) * (*(*eta))
                                            * (*(*eta)) * (*(*eta))
                                        + ce[m_darts_counter_temp2311][12] * (*(*zeta)) * (*(*zeta))
                                            * (*(*zeta)) * (*(*zeta));
                                }
                                (*m) = m_darts_counter_temp2311;
                            }
                        }
                        (*k) = k_darts_counter_temp2311;
                    }
                }
                (*j) = j_darts_counter_temp2311;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets2311[0].decDep();
}
TP2311::TP2311(int in_numThreads, int in_mainCodeletID, TP2241* in_TPParent, int in_initIteration,
    int in_lastIteration, TP2311** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , eta_darts2311(new double*[this->numThreads])
    , i_darts2311(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iglob_darts2311(new int*[this->numThreads])
    , j_darts2311(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jglob_darts2311(new int*[this->numThreads])
    , k_darts2311(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts2311(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , xi_darts2311(new double*[this->numThreads])
    , zeta_darts2311(new double*[this->numThreads])
    , initIteration2311(in_initIteration)
    , lastIteration2311(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets2311(new _barrierCodelets2311[1])
    , checkInCodelets2312(new _checkInCodelets2312[this->numThreads])
{
    /*Initialize the loop parameters*/
    range2311 = abs(lastIteration2311 - initIteration2311) / 1;
    rangePerCodelet2311 = range2311 / numThreads;
    minIteration2311 = min<int>(lastIteration2311, initIteration2311);
    remainderRange2311 = range2311 % numThreads;
    /*Initialize inputs and vars.*/
    this->eta_darts2311
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->iglob_darts2311 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jglob_darts2311 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->xi_darts2311
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->zeta_darts2311
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets2311[0] = _barrierCodelets2311(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets2312* checkInCodelets2312Ptr = (this->checkInCodelets2312);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets2312);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets2312Ptr) = _checkInCodelets2312(2, 1, this, codeletCounter);
#else
        (*checkInCodelets2312Ptr) = _checkInCodelets2312(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets2312Ptr).decDep();
        checkInCodelets2312Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP2311::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets2312[localID].setID(codeletID);
    this->checkInCodelets2312[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets2312[localID + this->baseNumThreads * i]
            = _checkInCodelets2312(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets2312[localID + this->baseNumThreads * i]
            = _checkInCodelets2312(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets2312[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets2312[localID + this->baseNumThreads * i].decDep();
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
TP2311::~TP2311()
{
    delete[] eta_darts2311;
    delete[] iglob_darts2311;
    delete[] jglob_darts2311;
    delete[] xi_darts2311;
    delete[] zeta_darts2311;
    delete[] barrierCodelets2311;
    delete[] checkInCodelets2312;
}
/*TP2446: OMPForDirective*/
void TP2446::_barrierCodelets2446::fire(void)
{
    TP2446* myTP = static_cast<TP2446*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets2446[0].decDep();
}
bool TP2446::requestNewRangeIterations2446(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet2446 * codeletID;
        int tempEndRange = rangePerCodelet2446 * (codeletID + 1);
        if (remainderRange2446 != 0) {
            if (codeletID < (uint32_t)remainderRange2446) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange2446;
                tempEndRange += remainderRange2446;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration2446;
        tempEndRange = tempEndRange * 1 + minIteration2446;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration2446 < lastIteration2446) {
            (this->inputsTPParent->i_darts2446[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts2446[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration2446;
        }
    }
    return isThereNewIteration;
}
void TP2446::_checkInCodelets2447::fire(void)
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
    this->inputsTPParent->L1_darts2446[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L1_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L2_darts2446[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->q_darts2446[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21_darts2446[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21_darts2241[this->getID()]);

    /*printing node 2447: ForStmt*/
    /*var: L1*/
    /*var: L2*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: q*/
    /*var: u21*/
    int* i = &(this->inputsTPParent->i_darts2446[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts2446[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts2446[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts2446[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** u21 = &(this->inputsTPParent->u21_darts2446[this->getLocalID()]);
    (void)u21 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations2446(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets2446[0].decDep();
        return;
    }
    for (int i_darts_counter_temp2446 = (*i); i_darts_counter_temp2446 <= endRange
         && i_darts_counter_temp2446 <= this->inputsTPParent->lastIteration2446;
         i_darts_counter_temp2446++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp2446 = (*j);
                for (; j_darts_counter_temp2446 <= jend; j_darts_counter_temp2446++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp2446 = (*k);
                        for (; k_darts_counter_temp2446 < nz - 1; k_darts_counter_temp2446++) {
                            flux[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                [k_darts_counter_temp2446][0]
                                = rsd[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                     [k_darts_counter_temp2446][1];
                            (*(*u21)) = rsd[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                           [k_darts_counter_temp2446][1]
                                / rsd[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                     [k_darts_counter_temp2446][0];
                            (*(*q)) = 0.5
                                * (rsd[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                      [k_darts_counter_temp2446][1]
                                        * rsd[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                             [k_darts_counter_temp2446][1]
                                    + rsd[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                         [k_darts_counter_temp2446][2]
                                        * rsd[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                             [k_darts_counter_temp2446][2]
                                    + rsd[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                         [k_darts_counter_temp2446][3]
                                        * rsd[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                             [k_darts_counter_temp2446][3])
                                / rsd[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                     [k_darts_counter_temp2446][0];
                            flux[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                [k_darts_counter_temp2446][1]
                                = rsd[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                     [k_darts_counter_temp2446][1]
                                    * (*(*u21))
                                + 0.40000000000000002
                                    * (rsd[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                          [k_darts_counter_temp2446][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                [k_darts_counter_temp2446][2]
                                = rsd[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                     [k_darts_counter_temp2446][2]
                                * (*(*u21));
                            flux[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                [k_darts_counter_temp2446][3]
                                = rsd[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                     [k_darts_counter_temp2446][3]
                                * (*(*u21));
                            flux[(i_darts_counter_temp2446)][j_darts_counter_temp2446]
                                [k_darts_counter_temp2446][4]
                                = (1.3999999999999999
                                          * rsd[(i_darts_counter_temp2446)]
                                               [j_darts_counter_temp2446][k_darts_counter_temp2446]
                                               [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u21));
                        }
                        (*k) = k_darts_counter_temp2446;
                    }
                }
                (*j) = j_darts_counter_temp2446;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets2446[0].decDep();
}
TP2446::TP2446(int in_numThreads, int in_mainCodeletID, TP2241* in_TPParent, int in_initIteration,
    int in_lastIteration, TP2446** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts2446(new int*[this->numThreads])
    , L2_darts2446(new int*[this->numThreads])
    , i_darts2446(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts2446(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts2446(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts2446(new double*[this->numThreads])
    , u21_darts2446(new double*[this->numThreads])
    , initIteration2446(in_initIteration)
    , lastIteration2446(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets2446(new _barrierCodelets2446[1])
    , checkInCodelets2447(new _checkInCodelets2447[this->numThreads])
{
    /*Initialize the loop parameters*/
    range2446 = abs(lastIteration2446 - initIteration2446) / 1;
    rangePerCodelet2446 = range2446 / numThreads;
    minIteration2446 = min<int>(lastIteration2446, initIteration2446);
    remainderRange2446 = range2446 % numThreads;
    /*Initialize inputs and vars.*/
    this->L1_darts2446 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->L2_darts2446 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->q_darts2446 = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21_darts2446
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets2446[0] = _barrierCodelets2446(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets2447* checkInCodelets2447Ptr = (this->checkInCodelets2447);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets2447);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets2447Ptr) = _checkInCodelets2447(2, 1, this, codeletCounter);
#else
        (*checkInCodelets2447Ptr) = _checkInCodelets2447(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets2447Ptr).decDep();
        checkInCodelets2447Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP2446::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets2447[localID].setID(codeletID);
    this->checkInCodelets2447[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets2447[localID + this->baseNumThreads * i]
            = _checkInCodelets2447(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets2447[localID + this->baseNumThreads * i]
            = _checkInCodelets2447(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets2447[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets2447[localID + this->baseNumThreads * i].decDep();
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
TP2446::~TP2446()
{
    delete[] L1_darts2446;
    delete[] L2_darts2446;
    delete[] q_darts2446;
    delete[] u21_darts2446;
    delete[] barrierCodelets2446;
    delete[] checkInCodelets2447;
}
/*TP2595: OMPForDirective*/
void TP2595::_barrierCodelets2595::fire(void)
{
    TP2595* myTP = static_cast<TP2595*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets2595[0].decDep();
}
bool TP2595::requestNewRangeIterations2595(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet2595 * codeletID;
        int tempEndRange = rangePerCodelet2595 * (codeletID + 1);
        if (remainderRange2595 != 0) {
            if (codeletID < (uint32_t)remainderRange2595) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange2595;
                tempEndRange += remainderRange2595;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration2595;
        tempEndRange = tempEndRange * 1 + minIteration2595;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration2595 < lastIteration2595) {
            (this->inputsTPParent->j_darts2595[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->j_darts2595[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration2595;
        }
    }
    return isThereNewIteration;
}
void TP2595::_checkInCodelets2596::fire(void)
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
    this->inputsTPParent->L2_darts2595[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->dsspm_darts2595[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->dsspm_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->iend1_darts2595[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iend1_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ist1_darts2595[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->ist1_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts2595[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21i_darts2595[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21i_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21im1_darts2595[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21im1_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31i_darts2595[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31i_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31im1_darts2595[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31im1_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41i_darts2595[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41i_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41im1_darts2595[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41im1_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51i_darts2595[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51i_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51im1_darts2595[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51im1_darts2241[this->getID()]);

    /*printing node 2596: ForStmt*/
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
    int** L2 = &(this->inputsTPParent->L2_darts2595[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    double** dsspm = &(this->inputsTPParent->dsspm_darts2595[this->getLocalID()]);
    (void)dsspm /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts2595[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iend1 = &(this->inputsTPParent->iend1_darts2595[this->getLocalID()]);
    (void)iend1 /*OMP_SHARED_PRIVATE*/;
    int** ist1 = &(this->inputsTPParent->ist1_darts2595[this->getLocalID()]);
    (void)ist1 /*OMP_SHARED_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts2595[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts2595[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts2595[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts2595[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21i = &(this->inputsTPParent->u21i_darts2595[this->getLocalID()]);
    (void)u21i /*OMP_SHARED_PRIVATE*/;
    double** u21im1 = &(this->inputsTPParent->u21im1_darts2595[this->getLocalID()]);
    (void)u21im1 /*OMP_SHARED_PRIVATE*/;
    double** u31i = &(this->inputsTPParent->u31i_darts2595[this->getLocalID()]);
    (void)u31i /*OMP_SHARED_PRIVATE*/;
    double** u31im1 = &(this->inputsTPParent->u31im1_darts2595[this->getLocalID()]);
    (void)u31im1 /*OMP_SHARED_PRIVATE*/;
    double** u41i = &(this->inputsTPParent->u41i_darts2595[this->getLocalID()]);
    (void)u41i /*OMP_SHARED_PRIVATE*/;
    double** u41im1 = &(this->inputsTPParent->u41im1_darts2595[this->getLocalID()]);
    (void)u41im1 /*OMP_SHARED_PRIVATE*/;
    double** u51i = &(this->inputsTPParent->u51i_darts2595[this->getLocalID()]);
    (void)u51i /*OMP_SHARED_PRIVATE*/;
    double** u51im1 = &(this->inputsTPParent->u51im1_darts2595[this->getLocalID()]);
    (void)u51im1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations2595(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets2595[0].decDep();
        return;
    }
    for (int j_darts_counter_temp2595 = (*j); j_darts_counter_temp2595 <= endRange
         && j_darts_counter_temp2595 <= this->inputsTPParent->lastIteration2595;
         j_darts_counter_temp2595++) {
        {
            {
                /*Loop's init*/
                (*k) = 1;
                int k_darts_counter_temp2595 = (*k);
                for (; k_darts_counter_temp2595 <= nz - 2; k_darts_counter_temp2595++) {
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp2595 = (*i);
                        for (; i_darts_counter_temp2595 <= iend; i_darts_counter_temp2595++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp2595 = (*m);
                                for (; m_darts_counter_temp2595 < 5; m_darts_counter_temp2595++) {
                                    frct[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                        [k_darts_counter_temp2595][m_darts_counter_temp2595]
                                        = frct[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                              [k_darts_counter_temp2595][m_darts_counter_temp2595]
                                        - tx2
                                            * (flux[i_darts_counter_temp2595 + 1]
                                                   [(j_darts_counter_temp2595)]
                                                   [k_darts_counter_temp2595]
                                                   [m_darts_counter_temp2595]
                                                - flux[i_darts_counter_temp2595 - 1]
                                                      [(j_darts_counter_temp2595)]
                                                      [k_darts_counter_temp2595]
                                                      [m_darts_counter_temp2595]);
                                }
                                (*m) = m_darts_counter_temp2595;
                            }
                        }
                        (*i) = i_darts_counter_temp2595;
                    }
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp2595 = (*i);
                        for (; i_darts_counter_temp2595 <= (*(*L2)); i_darts_counter_temp2595++) {
                            (*(*tmp)) = 1.
                                / rsd[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                     [k_darts_counter_temp2595][0];
                            (*(*u21i)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                     [k_darts_counter_temp2595][1];
                            (*(*u31i)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                     [k_darts_counter_temp2595][2];
                            (*(*u41i)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                     [k_darts_counter_temp2595][3];
                            (*(*u51i)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                     [k_darts_counter_temp2595][4];
                            (*(*tmp)) = 1.
                                / rsd[i_darts_counter_temp2595 - 1][(j_darts_counter_temp2595)]
                                     [k_darts_counter_temp2595][0];
                            (*(*u21im1)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2595 - 1][(j_darts_counter_temp2595)]
                                     [k_darts_counter_temp2595][1];
                            (*(*u31im1)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2595 - 1][(j_darts_counter_temp2595)]
                                     [k_darts_counter_temp2595][2];
                            (*(*u41im1)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2595 - 1][(j_darts_counter_temp2595)]
                                     [k_darts_counter_temp2595][3];
                            (*(*u51im1)) = (*(*tmp))
                                * rsd[i_darts_counter_temp2595 - 1][(j_darts_counter_temp2595)]
                                     [k_darts_counter_temp2595][4];
                            flux[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                [k_darts_counter_temp2595][1]
                                = (4. / 3.) * tx3 * ((*(*u21i)) - (*(*u21im1)));
                            flux[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                [k_darts_counter_temp2595][2]
                                = tx3 * ((*(*u31i)) - (*(*u31im1)));
                            flux[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                [k_darts_counter_temp2595][3]
                                = tx3 * ((*(*u41i)) - (*(*u41im1)));
                            flux[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                [k_darts_counter_temp2595][4]
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
                        (*i) = i_darts_counter_temp2595;
                    }
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp2595 = (*i);
                        for (; i_darts_counter_temp2595 <= iend; i_darts_counter_temp2595++) {
                            frct[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                [k_darts_counter_temp2595][0]
                                = frct[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                      [k_darts_counter_temp2595][0]
                                + dx1 * tx1
                                    * (rsd[i_darts_counter_temp2595 - 1][(j_darts_counter_temp2595)]
                                          [k_darts_counter_temp2595][0]
                                        - 2.
                                            * rsd[i_darts_counter_temp2595]
                                                 [(j_darts_counter_temp2595)]
                                                 [k_darts_counter_temp2595][0]
                                        + rsd[i_darts_counter_temp2595 + 1]
                                             [(j_darts_counter_temp2595)][k_darts_counter_temp2595]
                                             [0]);
                            frct[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                [k_darts_counter_temp2595][1]
                                = frct[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                      [k_darts_counter_temp2595][1]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp2595 + 1]
                                           [(j_darts_counter_temp2595)][k_darts_counter_temp2595][1]
                                        - flux[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                              [k_darts_counter_temp2595][1])
                                + dx2 * tx1
                                    * (rsd[i_darts_counter_temp2595 - 1][(j_darts_counter_temp2595)]
                                          [k_darts_counter_temp2595][1]
                                        - 2.
                                            * rsd[i_darts_counter_temp2595]
                                                 [(j_darts_counter_temp2595)]
                                                 [k_darts_counter_temp2595][1]
                                        + rsd[i_darts_counter_temp2595 + 1]
                                             [(j_darts_counter_temp2595)][k_darts_counter_temp2595]
                                             [1]);
                            frct[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                [k_darts_counter_temp2595][2]
                                = frct[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                      [k_darts_counter_temp2595][2]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp2595 + 1]
                                           [(j_darts_counter_temp2595)][k_darts_counter_temp2595][2]
                                        - flux[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                              [k_darts_counter_temp2595][2])
                                + dx3 * tx1
                                    * (rsd[i_darts_counter_temp2595 - 1][(j_darts_counter_temp2595)]
                                          [k_darts_counter_temp2595][2]
                                        - 2.
                                            * rsd[i_darts_counter_temp2595]
                                                 [(j_darts_counter_temp2595)]
                                                 [k_darts_counter_temp2595][2]
                                        + rsd[i_darts_counter_temp2595 + 1]
                                             [(j_darts_counter_temp2595)][k_darts_counter_temp2595]
                                             [2]);
                            frct[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                [k_darts_counter_temp2595][3]
                                = frct[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                      [k_darts_counter_temp2595][3]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp2595 + 1]
                                           [(j_darts_counter_temp2595)][k_darts_counter_temp2595][3]
                                        - flux[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                              [k_darts_counter_temp2595][3])
                                + dx4 * tx1
                                    * (rsd[i_darts_counter_temp2595 - 1][(j_darts_counter_temp2595)]
                                          [k_darts_counter_temp2595][3]
                                        - 2.
                                            * rsd[i_darts_counter_temp2595]
                                                 [(j_darts_counter_temp2595)]
                                                 [k_darts_counter_temp2595][3]
                                        + rsd[i_darts_counter_temp2595 + 1]
                                             [(j_darts_counter_temp2595)][k_darts_counter_temp2595]
                                             [3]);
                            frct[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                [k_darts_counter_temp2595][4]
                                = frct[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                      [k_darts_counter_temp2595][4]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp2595 + 1]
                                           [(j_darts_counter_temp2595)][k_darts_counter_temp2595][4]
                                        - flux[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                              [k_darts_counter_temp2595][4])
                                + dx5 * tx1
                                    * (rsd[i_darts_counter_temp2595 - 1][(j_darts_counter_temp2595)]
                                          [k_darts_counter_temp2595][4]
                                        - 2.
                                            * rsd[i_darts_counter_temp2595]
                                                 [(j_darts_counter_temp2595)]
                                                 [k_darts_counter_temp2595][4]
                                        + rsd[i_darts_counter_temp2595 + 1]
                                             [(j_darts_counter_temp2595)][k_darts_counter_temp2595]
                                             [4]);
                        }
                        (*i) = i_darts_counter_temp2595;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp2595 = (*m);
                        for (; m_darts_counter_temp2595 < 5; m_darts_counter_temp2595++) {
                            frct[1][(j_darts_counter_temp2595)][k_darts_counter_temp2595]
                                [m_darts_counter_temp2595]
                                = frct[1][(j_darts_counter_temp2595)][k_darts_counter_temp2595]
                                      [m_darts_counter_temp2595]
                                - (*(*dsspm))
                                    * (+5.
                                            * rsd[1][(j_darts_counter_temp2595)]
                                                 [k_darts_counter_temp2595]
                                                 [m_darts_counter_temp2595]
                                        - 4.
                                            * rsd[2][(j_darts_counter_temp2595)]
                                                 [k_darts_counter_temp2595]
                                                 [m_darts_counter_temp2595]
                                        + rsd[3][(j_darts_counter_temp2595)]
                                             [k_darts_counter_temp2595][m_darts_counter_temp2595]);
                            frct[2][(j_darts_counter_temp2595)][k_darts_counter_temp2595]
                                [m_darts_counter_temp2595]
                                = frct[2][(j_darts_counter_temp2595)][k_darts_counter_temp2595]
                                      [m_darts_counter_temp2595]
                                - (*(*dsspm))
                                    * (-4.
                                            * rsd[1][(j_darts_counter_temp2595)]
                                                 [k_darts_counter_temp2595]
                                                 [m_darts_counter_temp2595]
                                        + 6.
                                            * rsd[2][(j_darts_counter_temp2595)]
                                                 [k_darts_counter_temp2595]
                                                 [m_darts_counter_temp2595]
                                        - 4.
                                            * rsd[3][(j_darts_counter_temp2595)]
                                                 [k_darts_counter_temp2595]
                                                 [m_darts_counter_temp2595]
                                        + rsd[4][(j_darts_counter_temp2595)]
                                             [k_darts_counter_temp2595][m_darts_counter_temp2595]);
                        }
                        (*m) = m_darts_counter_temp2595;
                    }
                    (*(*ist1)) = 3;
                    (*(*iend1)) = nx - 4;
                    {
                        /*Loop's init*/
                        (*i) = (*(*ist1));
                        int i_darts_counter_temp2595 = (*i);
                        for (; i_darts_counter_temp2595 <= (*(*iend1));
                             i_darts_counter_temp2595++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp2595 = (*m);
                                for (; m_darts_counter_temp2595 < 5; m_darts_counter_temp2595++) {
                                    frct[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                        [k_darts_counter_temp2595][m_darts_counter_temp2595]
                                        = frct[i_darts_counter_temp2595][(j_darts_counter_temp2595)]
                                              [k_darts_counter_temp2595][m_darts_counter_temp2595]
                                        - (*(*dsspm))
                                            * (rsd[i_darts_counter_temp2595 - 2]
                                                  [(j_darts_counter_temp2595)]
                                                  [k_darts_counter_temp2595]
                                                  [m_darts_counter_temp2595]
                                                - 4.
                                                    * rsd[i_darts_counter_temp2595 - 1]
                                                         [(j_darts_counter_temp2595)]
                                                         [k_darts_counter_temp2595]
                                                         [m_darts_counter_temp2595]
                                                + 6.
                                                    * rsd[i_darts_counter_temp2595]
                                                         [(j_darts_counter_temp2595)]
                                                         [k_darts_counter_temp2595]
                                                         [m_darts_counter_temp2595]
                                                - 4.
                                                    * rsd[i_darts_counter_temp2595 + 1]
                                                         [(j_darts_counter_temp2595)]
                                                         [k_darts_counter_temp2595]
                                                         [m_darts_counter_temp2595]
                                                + rsd[i_darts_counter_temp2595 + 2]
                                                     [(j_darts_counter_temp2595)]
                                                     [k_darts_counter_temp2595]
                                                     [m_darts_counter_temp2595]);
                                }
                                (*m) = m_darts_counter_temp2595;
                            }
                        }
                        (*i) = i_darts_counter_temp2595;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp2595 = (*m);
                        for (; m_darts_counter_temp2595 < 5; m_darts_counter_temp2595++) {
                            frct[nx - 3][(j_darts_counter_temp2595)][k_darts_counter_temp2595]
                                [m_darts_counter_temp2595]
                                = frct[nx - 3][(j_darts_counter_temp2595)][k_darts_counter_temp2595]
                                      [m_darts_counter_temp2595]
                                - (*(*dsspm))
                                    * (rsd[nx - 5][(j_darts_counter_temp2595)]
                                          [k_darts_counter_temp2595][m_darts_counter_temp2595]
                                        - 4.
                                            * rsd[nx - 4][(j_darts_counter_temp2595)]
                                                 [k_darts_counter_temp2595]
                                                 [m_darts_counter_temp2595]
                                        + 6.
                                            * rsd[nx - 3][(j_darts_counter_temp2595)]
                                                 [k_darts_counter_temp2595]
                                                 [m_darts_counter_temp2595]
                                        - 4.
                                            * rsd[nx - 2][(j_darts_counter_temp2595)]
                                                 [k_darts_counter_temp2595]
                                                 [m_darts_counter_temp2595]);
                            frct[nx - 2][(j_darts_counter_temp2595)][k_darts_counter_temp2595]
                                [m_darts_counter_temp2595]
                                = frct[nx - 2][(j_darts_counter_temp2595)][k_darts_counter_temp2595]
                                      [m_darts_counter_temp2595]
                                - (*(*dsspm))
                                    * (rsd[nx - 4][(j_darts_counter_temp2595)]
                                          [k_darts_counter_temp2595][m_darts_counter_temp2595]
                                        - 4.
                                            * rsd[nx - 3][(j_darts_counter_temp2595)]
                                                 [k_darts_counter_temp2595]
                                                 [m_darts_counter_temp2595]
                                        + 5.
                                            * rsd[nx - 2][(j_darts_counter_temp2595)]
                                                 [k_darts_counter_temp2595]
                                                 [m_darts_counter_temp2595]);
                        }
                        (*m) = m_darts_counter_temp2595;
                    }
                }
                (*k) = k_darts_counter_temp2595;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets2595[0].decDep();
}
TP2595::TP2595(int in_numThreads, int in_mainCodeletID, TP2241* in_TPParent, int in_initIteration,
    int in_lastIteration, TP2595** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L2_darts2595(new int*[this->numThreads])
    , dsspm_darts2595(new double*[this->numThreads])
    , i_darts2595(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iend1_darts2595(new int*[this->numThreads])
    , ist1_darts2595(new int*[this->numThreads])
    , j_darts2595(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts2595(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts2595(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts2595(new double*[this->numThreads])
    , u21i_darts2595(new double*[this->numThreads])
    , u21im1_darts2595(new double*[this->numThreads])
    , u31i_darts2595(new double*[this->numThreads])
    , u31im1_darts2595(new double*[this->numThreads])
    , u41i_darts2595(new double*[this->numThreads])
    , u41im1_darts2595(new double*[this->numThreads])
    , u51i_darts2595(new double*[this->numThreads])
    , u51im1_darts2595(new double*[this->numThreads])
    , initIteration2595(in_initIteration)
    , lastIteration2595(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets2595(new _barrierCodelets2595[1])
    , checkInCodelets2596(new _checkInCodelets2596[this->numThreads])
{
    /*Initialize the loop parameters*/
    range2595 = abs(lastIteration2595 - initIteration2595) / 1;
    rangePerCodelet2595 = range2595 / numThreads;
    minIteration2595 = min<int>(lastIteration2595, initIteration2595);
    remainderRange2595 = range2595 % numThreads;
    /*Initialize inputs and vars.*/
    this->L2_darts2595 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->dsspm_darts2595
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->iend1_darts2595 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ist1_darts2595 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts2595
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21i_darts2595
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21im1_darts2595
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31i_darts2595
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31im1_darts2595
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41i_darts2595
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41im1_darts2595
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51i_darts2595
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51im1_darts2595
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets2595[0] = _barrierCodelets2595(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets2596* checkInCodelets2596Ptr = (this->checkInCodelets2596);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets2596);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets2596Ptr) = _checkInCodelets2596(2, 1, this, codeletCounter);
#else
        (*checkInCodelets2596Ptr) = _checkInCodelets2596(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets2596Ptr).decDep();
        checkInCodelets2596Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP2595::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets2596[localID].setID(codeletID);
    this->checkInCodelets2596[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets2596[localID + this->baseNumThreads * i]
            = _checkInCodelets2596(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets2596[localID + this->baseNumThreads * i]
            = _checkInCodelets2596(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets2596[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets2596[localID + this->baseNumThreads * i].decDep();
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
TP2595::~TP2595()
{
    delete[] L2_darts2595;
    delete[] dsspm_darts2595;
    delete[] iend1_darts2595;
    delete[] ist1_darts2595;
    delete[] tmp_darts2595;
    delete[] u21i_darts2595;
    delete[] u21im1_darts2595;
    delete[] u31i_darts2595;
    delete[] u31im1_darts2595;
    delete[] u41i_darts2595;
    delete[] u41im1_darts2595;
    delete[] u51i_darts2595;
    delete[] u51im1_darts2595;
    delete[] barrierCodelets2595;
    delete[] checkInCodelets2596;
}
/*TP3233: OMPForDirective*/
void TP3233::_barrierCodelets3233::fire(void)
{
    TP3233* myTP = static_cast<TP3233*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets3233[0].decDep();
}
bool TP3233::requestNewRangeIterations3233(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet3233 * codeletID;
        int tempEndRange = rangePerCodelet3233 * (codeletID + 1);
        if (remainderRange3233 != 0) {
            if (codeletID < (uint32_t)remainderRange3233) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange3233;
                tempEndRange += remainderRange3233;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration3233;
        tempEndRange = tempEndRange * 1 + minIteration3233;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration3233 < lastIteration3233) {
            (this->inputsTPParent->i_darts3233[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts3233[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration3233;
        }
    }
    return isThereNewIteration;
}
void TP3233::_checkInCodelets3234::fire(void)
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
    this->inputsTPParent->L1_darts3233[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L1_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L2_darts3233[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->q_darts3233[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31_darts3233[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31_darts2241[this->getID()]);

    /*printing node 3234: ForStmt*/
    /*var: L1*/
    /*var: L2*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: q*/
    /*var: u31*/
    int** L1 = &(this->inputsTPParent->L1_darts3233[this->getLocalID()]);
    (void)L1 /*OMP_SHARED_PRIVATE*/;
    int** L2 = &(this->inputsTPParent->L2_darts3233[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts3233[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts3233[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts3233[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts3233[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** u31 = &(this->inputsTPParent->u31_darts3233[this->getLocalID()]);
    (void)u31 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations3233(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets3233[0].decDep();
        return;
    }
    for (int i_darts_counter_temp3233 = (*i); i_darts_counter_temp3233 <= endRange
         && i_darts_counter_temp3233 <= this->inputsTPParent->lastIteration3233;
         i_darts_counter_temp3233++) {
        {
            {
                /*Loop's init*/
                (*j) = (*(*L1));
                int j_darts_counter_temp3233 = (*j);
                for (; j_darts_counter_temp3233 <= (*(*L2)); j_darts_counter_temp3233++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp3233 = (*k);
                        for (; k_darts_counter_temp3233 <= nz - 2; k_darts_counter_temp3233++) {
                            flux[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                [k_darts_counter_temp3233][0]
                                = rsd[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                     [k_darts_counter_temp3233][2];
                            (*(*u31)) = rsd[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                           [k_darts_counter_temp3233][2]
                                / rsd[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                     [k_darts_counter_temp3233][0];
                            (*(*q)) = 0.5
                                * (rsd[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                      [k_darts_counter_temp3233][1]
                                        * rsd[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                             [k_darts_counter_temp3233][1]
                                    + rsd[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                         [k_darts_counter_temp3233][2]
                                        * rsd[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                             [k_darts_counter_temp3233][2]
                                    + rsd[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                         [k_darts_counter_temp3233][3]
                                        * rsd[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                             [k_darts_counter_temp3233][3])
                                / rsd[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                     [k_darts_counter_temp3233][0];
                            flux[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                [k_darts_counter_temp3233][1]
                                = rsd[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                     [k_darts_counter_temp3233][1]
                                * (*(*u31));
                            flux[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                [k_darts_counter_temp3233][2]
                                = rsd[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                     [k_darts_counter_temp3233][2]
                                    * (*(*u31))
                                + 0.40000000000000002
                                    * (rsd[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                          [k_darts_counter_temp3233][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                [k_darts_counter_temp3233][3]
                                = rsd[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                     [k_darts_counter_temp3233][3]
                                * (*(*u31));
                            flux[(i_darts_counter_temp3233)][j_darts_counter_temp3233]
                                [k_darts_counter_temp3233][4]
                                = (1.3999999999999999
                                          * rsd[(i_darts_counter_temp3233)]
                                               [j_darts_counter_temp3233][k_darts_counter_temp3233]
                                               [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u31));
                        }
                        (*k) = k_darts_counter_temp3233;
                    }
                }
                (*j) = j_darts_counter_temp3233;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets3233[0].decDep();
}
TP3233::TP3233(int in_numThreads, int in_mainCodeletID, TP2241* in_TPParent, int in_initIteration,
    int in_lastIteration, TP3233** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts3233(new int*[this->numThreads])
    , L2_darts3233(new int*[this->numThreads])
    , i_darts3233(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts3233(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts3233(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts3233(new double*[this->numThreads])
    , u31_darts3233(new double*[this->numThreads])
    , initIteration3233(in_initIteration)
    , lastIteration3233(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets3233(new _barrierCodelets3233[1])
    , checkInCodelets3234(new _checkInCodelets3234[this->numThreads])
{
    /*Initialize the loop parameters*/
    range3233 = abs(lastIteration3233 - initIteration3233) / 1;
    rangePerCodelet3233 = range3233 / numThreads;
    minIteration3233 = min<int>(lastIteration3233, initIteration3233);
    remainderRange3233 = range3233 % numThreads;
    /*Initialize inputs and vars.*/
    this->L1_darts3233 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->L2_darts3233 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->q_darts3233 = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31_darts3233
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets3233[0] = _barrierCodelets3233(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets3234* checkInCodelets3234Ptr = (this->checkInCodelets3234);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets3234);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets3234Ptr) = _checkInCodelets3234(2, 1, this, codeletCounter);
#else
        (*checkInCodelets3234Ptr) = _checkInCodelets3234(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets3234Ptr).decDep();
        checkInCodelets3234Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP3233::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets3234[localID].setID(codeletID);
    this->checkInCodelets3234[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets3234[localID + this->baseNumThreads * i]
            = _checkInCodelets3234(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets3234[localID + this->baseNumThreads * i]
            = _checkInCodelets3234(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets3234[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets3234[localID + this->baseNumThreads * i].decDep();
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
TP3233::~TP3233()
{
    delete[] L1_darts3233;
    delete[] L2_darts3233;
    delete[] q_darts3233;
    delete[] u31_darts3233;
    delete[] barrierCodelets3233;
    delete[] checkInCodelets3234;
}
/*TP3382: OMPForDirective*/
void TP3382::_barrierCodelets3382::fire(void)
{
    TP3382* myTP = static_cast<TP3382*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets3382[0].decDep();
}
bool TP3382::requestNewRangeIterations3382(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet3382 * codeletID;
        int tempEndRange = rangePerCodelet3382 * (codeletID + 1);
        if (remainderRange3382 != 0) {
            if (codeletID < (uint32_t)remainderRange3382) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange3382;
                tempEndRange += remainderRange3382;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration3382;
        tempEndRange = tempEndRange * 1 + minIteration3382;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration3382 < lastIteration3382) {
            (this->inputsTPParent->i_darts3382[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts3382[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration3382;
        }
    }
    return isThereNewIteration;
}
void TP3382::_checkInCodelets3383::fire(void)
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
    this->inputsTPParent->L2_darts3382[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->dsspm_darts3382[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->dsspm_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jend1_darts3382[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jend1_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jst1_darts3382[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jst1_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts3382[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21j_darts3382[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21j_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21jm1_darts3382[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21jm1_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31j_darts3382[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31j_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31jm1_darts3382[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31jm1_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41j_darts3382[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41j_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41jm1_darts3382[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41jm1_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51j_darts3382[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51j_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51jm1_darts3382[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51jm1_darts2241[this->getID()]);

    /*printing node 3383: ForStmt*/
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
    int** L2 = &(this->inputsTPParent->L2_darts3382[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    double** dsspm = &(this->inputsTPParent->dsspm_darts3382[this->getLocalID()]);
    (void)dsspm /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts3382[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts3382[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jend1 = &(this->inputsTPParent->jend1_darts3382[this->getLocalID()]);
    (void)jend1 /*OMP_SHARED_PRIVATE*/;
    int** jst1 = &(this->inputsTPParent->jst1_darts3382[this->getLocalID()]);
    (void)jst1 /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts3382[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts3382[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts3382[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21j = &(this->inputsTPParent->u21j_darts3382[this->getLocalID()]);
    (void)u21j /*OMP_SHARED_PRIVATE*/;
    double** u21jm1 = &(this->inputsTPParent->u21jm1_darts3382[this->getLocalID()]);
    (void)u21jm1 /*OMP_SHARED_PRIVATE*/;
    double** u31j = &(this->inputsTPParent->u31j_darts3382[this->getLocalID()]);
    (void)u31j /*OMP_SHARED_PRIVATE*/;
    double** u31jm1 = &(this->inputsTPParent->u31jm1_darts3382[this->getLocalID()]);
    (void)u31jm1 /*OMP_SHARED_PRIVATE*/;
    double** u41j = &(this->inputsTPParent->u41j_darts3382[this->getLocalID()]);
    (void)u41j /*OMP_SHARED_PRIVATE*/;
    double** u41jm1 = &(this->inputsTPParent->u41jm1_darts3382[this->getLocalID()]);
    (void)u41jm1 /*OMP_SHARED_PRIVATE*/;
    double** u51j = &(this->inputsTPParent->u51j_darts3382[this->getLocalID()]);
    (void)u51j /*OMP_SHARED_PRIVATE*/;
    double** u51jm1 = &(this->inputsTPParent->u51jm1_darts3382[this->getLocalID()]);
    (void)u51jm1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations3382(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets3382[0].decDep();
        return;
    }
    for (int i_darts_counter_temp3382 = (*i); i_darts_counter_temp3382 <= endRange
         && i_darts_counter_temp3382 <= this->inputsTPParent->lastIteration3382;
         i_darts_counter_temp3382++) {
        {
            {
                /*Loop's init*/
                (*k) = 1;
                int k_darts_counter_temp3382 = (*k);
                for (; k_darts_counter_temp3382 <= nz - 2; k_darts_counter_temp3382++) {
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp3382 = (*j);
                        for (; j_darts_counter_temp3382 <= jend; j_darts_counter_temp3382++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp3382 = (*m);
                                for (; m_darts_counter_temp3382 < 5; m_darts_counter_temp3382++) {
                                    frct[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                        [k_darts_counter_temp3382][m_darts_counter_temp3382]
                                        = frct[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                              [k_darts_counter_temp3382][m_darts_counter_temp3382]
                                        - ty2
                                            * (flux[(i_darts_counter_temp3382)]
                                                   [j_darts_counter_temp3382 + 1]
                                                   [k_darts_counter_temp3382]
                                                   [m_darts_counter_temp3382]
                                                - flux[(i_darts_counter_temp3382)]
                                                      [j_darts_counter_temp3382 - 1]
                                                      [k_darts_counter_temp3382]
                                                      [m_darts_counter_temp3382]);
                                }
                                (*m) = m_darts_counter_temp3382;
                            }
                        }
                        (*j) = j_darts_counter_temp3382;
                    }
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp3382 = (*j);
                        for (; j_darts_counter_temp3382 <= (*(*L2)); j_darts_counter_temp3382++) {
                            (*(*tmp)) = 1.
                                / rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                     [k_darts_counter_temp3382][0];
                            (*(*u21j)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                     [k_darts_counter_temp3382][1];
                            (*(*u31j)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                     [k_darts_counter_temp3382][2];
                            (*(*u41j)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                     [k_darts_counter_temp3382][3];
                            (*(*u51j)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                     [k_darts_counter_temp3382][4];
                            (*(*tmp)) = 1.
                                / rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382 - 1]
                                     [k_darts_counter_temp3382][0];
                            (*(*u21jm1)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382 - 1]
                                     [k_darts_counter_temp3382][1];
                            (*(*u31jm1)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382 - 1]
                                     [k_darts_counter_temp3382][2];
                            (*(*u41jm1)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382 - 1]
                                     [k_darts_counter_temp3382][3];
                            (*(*u51jm1)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382 - 1]
                                     [k_darts_counter_temp3382][4];
                            flux[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                [k_darts_counter_temp3382][1]
                                = ty3 * ((*(*u21j)) - (*(*u21jm1)));
                            flux[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                [k_darts_counter_temp3382][2]
                                = (4. / 3.) * ty3 * ((*(*u31j)) - (*(*u31jm1)));
                            flux[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                [k_darts_counter_temp3382][3]
                                = ty3 * ((*(*u41j)) - (*(*u41jm1)));
                            flux[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                [k_darts_counter_temp3382][4]
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
                        (*j) = j_darts_counter_temp3382;
                    }
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp3382 = (*j);
                        for (; j_darts_counter_temp3382 <= jend; j_darts_counter_temp3382++) {
                            frct[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                [k_darts_counter_temp3382][0]
                                = frct[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                      [k_darts_counter_temp3382][0]
                                + dy1 * ty1
                                    * (rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382 - 1]
                                          [k_darts_counter_temp3382][0]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3382)]
                                                 [j_darts_counter_temp3382]
                                                 [k_darts_counter_temp3382][0]
                                        + rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382
                                            + 1][k_darts_counter_temp3382][0]);
                            frct[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                [k_darts_counter_temp3382][1]
                                = frct[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                      [k_darts_counter_temp3382][1]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp3382)][j_darts_counter_temp3382
                                           + 1][k_darts_counter_temp3382][1]
                                        - flux[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                              [k_darts_counter_temp3382][1])
                                + dy2 * ty1
                                    * (rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382 - 1]
                                          [k_darts_counter_temp3382][1]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3382)]
                                                 [j_darts_counter_temp3382]
                                                 [k_darts_counter_temp3382][1]
                                        + rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382
                                            + 1][k_darts_counter_temp3382][1]);
                            frct[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                [k_darts_counter_temp3382][2]
                                = frct[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                      [k_darts_counter_temp3382][2]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp3382)][j_darts_counter_temp3382
                                           + 1][k_darts_counter_temp3382][2]
                                        - flux[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                              [k_darts_counter_temp3382][2])
                                + dy3 * ty1
                                    * (rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382 - 1]
                                          [k_darts_counter_temp3382][2]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3382)]
                                                 [j_darts_counter_temp3382]
                                                 [k_darts_counter_temp3382][2]
                                        + rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382
                                            + 1][k_darts_counter_temp3382][2]);
                            frct[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                [k_darts_counter_temp3382][3]
                                = frct[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                      [k_darts_counter_temp3382][3]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp3382)][j_darts_counter_temp3382
                                           + 1][k_darts_counter_temp3382][3]
                                        - flux[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                              [k_darts_counter_temp3382][3])
                                + dy4 * ty1
                                    * (rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382 - 1]
                                          [k_darts_counter_temp3382][3]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3382)]
                                                 [j_darts_counter_temp3382]
                                                 [k_darts_counter_temp3382][3]
                                        + rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382
                                            + 1][k_darts_counter_temp3382][3]);
                            frct[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                [k_darts_counter_temp3382][4]
                                = frct[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                      [k_darts_counter_temp3382][4]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp3382)][j_darts_counter_temp3382
                                           + 1][k_darts_counter_temp3382][4]
                                        - flux[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                              [k_darts_counter_temp3382][4])
                                + dy5 * ty1
                                    * (rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382 - 1]
                                          [k_darts_counter_temp3382][4]
                                        - 2.
                                            * rsd[(i_darts_counter_temp3382)]
                                                 [j_darts_counter_temp3382]
                                                 [k_darts_counter_temp3382][4]
                                        + rsd[(i_darts_counter_temp3382)][j_darts_counter_temp3382
                                            + 1][k_darts_counter_temp3382][4]);
                        }
                        (*j) = j_darts_counter_temp3382;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp3382 = (*m);
                        for (; m_darts_counter_temp3382 < 5; m_darts_counter_temp3382++) {
                            frct[(i_darts_counter_temp3382)][1][k_darts_counter_temp3382]
                                [m_darts_counter_temp3382]
                                = frct[(i_darts_counter_temp3382)][1][k_darts_counter_temp3382]
                                      [m_darts_counter_temp3382]
                                - (*(*dsspm))
                                    * (+5.
                                            * rsd[(i_darts_counter_temp3382)][1]
                                                 [k_darts_counter_temp3382]
                                                 [m_darts_counter_temp3382]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3382)][2]
                                                 [k_darts_counter_temp3382]
                                                 [m_darts_counter_temp3382]
                                        + rsd[(i_darts_counter_temp3382)][3]
                                             [k_darts_counter_temp3382][m_darts_counter_temp3382]);
                            frct[(i_darts_counter_temp3382)][2][k_darts_counter_temp3382]
                                [m_darts_counter_temp3382]
                                = frct[(i_darts_counter_temp3382)][2][k_darts_counter_temp3382]
                                      [m_darts_counter_temp3382]
                                - (*(*dsspm))
                                    * (-4.
                                            * rsd[(i_darts_counter_temp3382)][1]
                                                 [k_darts_counter_temp3382]
                                                 [m_darts_counter_temp3382]
                                        + 6.
                                            * rsd[(i_darts_counter_temp3382)][2]
                                                 [k_darts_counter_temp3382]
                                                 [m_darts_counter_temp3382]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3382)][3]
                                                 [k_darts_counter_temp3382]
                                                 [m_darts_counter_temp3382]
                                        + rsd[(i_darts_counter_temp3382)][4]
                                             [k_darts_counter_temp3382][m_darts_counter_temp3382]);
                        }
                        (*m) = m_darts_counter_temp3382;
                    }
                    (*(*jst1)) = 3;
                    (*(*jend1)) = ny - 4;
                    {
                        /*Loop's init*/
                        (*j) = (*(*jst1));
                        int j_darts_counter_temp3382 = (*j);
                        for (; j_darts_counter_temp3382 <= (*(*jend1));
                             j_darts_counter_temp3382++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp3382 = (*m);
                                for (; m_darts_counter_temp3382 < 5; m_darts_counter_temp3382++) {
                                    frct[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                        [k_darts_counter_temp3382][m_darts_counter_temp3382]
                                        = frct[(i_darts_counter_temp3382)][j_darts_counter_temp3382]
                                              [k_darts_counter_temp3382][m_darts_counter_temp3382]
                                        - (*(*dsspm))
                                            * (rsd[(i_darts_counter_temp3382)]
                                                  [j_darts_counter_temp3382 - 2]
                                                  [k_darts_counter_temp3382]
                                                  [m_darts_counter_temp3382]
                                                - 4.
                                                    * rsd[(i_darts_counter_temp3382)]
                                                         [j_darts_counter_temp3382 - 1]
                                                         [k_darts_counter_temp3382]
                                                         [m_darts_counter_temp3382]
                                                + 6.
                                                    * rsd[(i_darts_counter_temp3382)]
                                                         [j_darts_counter_temp3382]
                                                         [k_darts_counter_temp3382]
                                                         [m_darts_counter_temp3382]
                                                - 4.
                                                    * rsd[(i_darts_counter_temp3382)]
                                                         [j_darts_counter_temp3382 + 1]
                                                         [k_darts_counter_temp3382]
                                                         [m_darts_counter_temp3382]
                                                + rsd[(i_darts_counter_temp3382)]
                                                     [j_darts_counter_temp3382 + 2]
                                                     [k_darts_counter_temp3382]
                                                     [m_darts_counter_temp3382]);
                                }
                                (*m) = m_darts_counter_temp3382;
                            }
                        }
                        (*j) = j_darts_counter_temp3382;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp3382 = (*m);
                        for (; m_darts_counter_temp3382 < 5; m_darts_counter_temp3382++) {
                            frct[(i_darts_counter_temp3382)][ny - 3][k_darts_counter_temp3382]
                                [m_darts_counter_temp3382]
                                = frct[(i_darts_counter_temp3382)][ny - 3][k_darts_counter_temp3382]
                                      [m_darts_counter_temp3382]
                                - (*(*dsspm))
                                    * (rsd[(i_darts_counter_temp3382)][ny - 5]
                                          [k_darts_counter_temp3382][m_darts_counter_temp3382]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3382)][ny - 4]
                                                 [k_darts_counter_temp3382]
                                                 [m_darts_counter_temp3382]
                                        + 6.
                                            * rsd[(i_darts_counter_temp3382)][ny - 3]
                                                 [k_darts_counter_temp3382]
                                                 [m_darts_counter_temp3382]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3382)][ny - 2]
                                                 [k_darts_counter_temp3382]
                                                 [m_darts_counter_temp3382]);
                            frct[(i_darts_counter_temp3382)][ny - 2][k_darts_counter_temp3382]
                                [m_darts_counter_temp3382]
                                = frct[(i_darts_counter_temp3382)][ny - 2][k_darts_counter_temp3382]
                                      [m_darts_counter_temp3382]
                                - (*(*dsspm))
                                    * (rsd[(i_darts_counter_temp3382)][ny - 4]
                                          [k_darts_counter_temp3382][m_darts_counter_temp3382]
                                        - 4.
                                            * rsd[(i_darts_counter_temp3382)][ny - 3]
                                                 [k_darts_counter_temp3382]
                                                 [m_darts_counter_temp3382]
                                        + 5.
                                            * rsd[(i_darts_counter_temp3382)][ny - 2]
                                                 [k_darts_counter_temp3382]
                                                 [m_darts_counter_temp3382]);
                        }
                        (*m) = m_darts_counter_temp3382;
                    }
                }
                (*k) = k_darts_counter_temp3382;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets3382[0].decDep();
}
TP3382::TP3382(int in_numThreads, int in_mainCodeletID, TP2241* in_TPParent, int in_initIteration,
    int in_lastIteration, TP3382** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L2_darts3382(new int*[this->numThreads])
    , dsspm_darts3382(new double*[this->numThreads])
    , i_darts3382(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts3382(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jend1_darts3382(new int*[this->numThreads])
    , jst1_darts3382(new int*[this->numThreads])
    , k_darts3382(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts3382(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts3382(new double*[this->numThreads])
    , u21j_darts3382(new double*[this->numThreads])
    , u21jm1_darts3382(new double*[this->numThreads])
    , u31j_darts3382(new double*[this->numThreads])
    , u31jm1_darts3382(new double*[this->numThreads])
    , u41j_darts3382(new double*[this->numThreads])
    , u41jm1_darts3382(new double*[this->numThreads])
    , u51j_darts3382(new double*[this->numThreads])
    , u51jm1_darts3382(new double*[this->numThreads])
    , initIteration3382(in_initIteration)
    , lastIteration3382(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets3382(new _barrierCodelets3382[1])
    , checkInCodelets3383(new _checkInCodelets3383[this->numThreads])
{
    /*Initialize the loop parameters*/
    range3382 = abs(lastIteration3382 - initIteration3382) / 1;
    rangePerCodelet3382 = range3382 / numThreads;
    minIteration3382 = min<int>(lastIteration3382, initIteration3382);
    remainderRange3382 = range3382 % numThreads;
    /*Initialize inputs and vars.*/
    this->L2_darts3382 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->dsspm_darts3382
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jend1_darts3382 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jst1_darts3382 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts3382
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21j_darts3382
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21jm1_darts3382
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31j_darts3382
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31jm1_darts3382
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41j_darts3382
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41jm1_darts3382
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51j_darts3382
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51jm1_darts3382
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets3382[0] = _barrierCodelets3382(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets3383* checkInCodelets3383Ptr = (this->checkInCodelets3383);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets3383);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets3383Ptr) = _checkInCodelets3383(2, 1, this, codeletCounter);
#else
        (*checkInCodelets3383Ptr) = _checkInCodelets3383(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets3383Ptr).decDep();
        checkInCodelets3383Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP3382::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets3383[localID].setID(codeletID);
    this->checkInCodelets3383[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets3383[localID + this->baseNumThreads * i]
            = _checkInCodelets3383(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets3383[localID + this->baseNumThreads * i]
            = _checkInCodelets3383(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets3383[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets3383[localID + this->baseNumThreads * i].decDep();
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
TP3382::~TP3382()
{
    delete[] L2_darts3382;
    delete[] dsspm_darts3382;
    delete[] jend1_darts3382;
    delete[] jst1_darts3382;
    delete[] tmp_darts3382;
    delete[] u21j_darts3382;
    delete[] u21jm1_darts3382;
    delete[] u31j_darts3382;
    delete[] u31jm1_darts3382;
    delete[] u41j_darts3382;
    delete[] u41jm1_darts3382;
    delete[] u51j_darts3382;
    delete[] u51jm1_darts3382;
    delete[] barrierCodelets3382;
    delete[] checkInCodelets3383;
}
/*TP4017: OMPForDirective*/
void TP4017::_barrierCodelets4017::fire(void)
{
    TP4017* myTP = static_cast<TP4017*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets4017[0].decDep();
}
bool TP4017::requestNewRangeIterations4017(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet4017 * codeletID;
        int tempEndRange = rangePerCodelet4017 * (codeletID + 1);
        if (remainderRange4017 != 0) {
            if (codeletID < (uint32_t)remainderRange4017) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange4017;
                tempEndRange += remainderRange4017;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration4017;
        tempEndRange = tempEndRange * 1 + minIteration4017;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration4017 < lastIteration4017) {
            (this->inputsTPParent->i_darts4017[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts4017[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration4017;
        }
    }
    return isThereNewIteration;
}
void TP4017::_checkInCodelets4018::fire(void)
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
    this->inputsTPParent->dsspm_darts4017[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->dsspm_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->q_darts4017[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts4017[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21k_darts4017[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21k_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21km1_darts4017[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21km1_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31k_darts4017[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31k_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31km1_darts4017[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31km1_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41_darts4017[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41k_darts4017[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41k_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41km1_darts4017[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41km1_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51k_darts4017[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51k_darts2241[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51km1_darts4017[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51km1_darts2241[this->getID()]);

    /*printing node 4018: ForStmt*/
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
    double** dsspm = &(this->inputsTPParent->dsspm_darts4017[this->getLocalID()]);
    (void)dsspm /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts4017[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts4017[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts4017[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts4017[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts4017[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts4017[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21k = &(this->inputsTPParent->u21k_darts4017[this->getLocalID()]);
    (void)u21k /*OMP_SHARED_PRIVATE*/;
    double** u21km1 = &(this->inputsTPParent->u21km1_darts4017[this->getLocalID()]);
    (void)u21km1 /*OMP_SHARED_PRIVATE*/;
    double** u31k = &(this->inputsTPParent->u31k_darts4017[this->getLocalID()]);
    (void)u31k /*OMP_SHARED_PRIVATE*/;
    double** u31km1 = &(this->inputsTPParent->u31km1_darts4017[this->getLocalID()]);
    (void)u31km1 /*OMP_SHARED_PRIVATE*/;
    double** u41 = &(this->inputsTPParent->u41_darts4017[this->getLocalID()]);
    (void)u41 /*OMP_SHARED_PRIVATE*/;
    double** u41k = &(this->inputsTPParent->u41k_darts4017[this->getLocalID()]);
    (void)u41k /*OMP_SHARED_PRIVATE*/;
    double** u41km1 = &(this->inputsTPParent->u41km1_darts4017[this->getLocalID()]);
    (void)u41km1 /*OMP_SHARED_PRIVATE*/;
    double** u51k = &(this->inputsTPParent->u51k_darts4017[this->getLocalID()]);
    (void)u51k /*OMP_SHARED_PRIVATE*/;
    double** u51km1 = &(this->inputsTPParent->u51km1_darts4017[this->getLocalID()]);
    (void)u51km1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations4017(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets4017[0].decDep();
        return;
    }
    for (int i_darts_counter_temp4017 = (*i); i_darts_counter_temp4017 <= endRange
         && i_darts_counter_temp4017 <= this->inputsTPParent->lastIteration4017;
         i_darts_counter_temp4017++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp4017 = (*j);
                for (; j_darts_counter_temp4017 <= jend; j_darts_counter_temp4017++) {
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp4017 = (*k);
                        for (; k_darts_counter_temp4017 <= nz - 1; k_darts_counter_temp4017++) {
                            flux[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                [k_darts_counter_temp4017][0]
                                = rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                     [k_darts_counter_temp4017][3];
                            (*(*u41)) = rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                           [k_darts_counter_temp4017][3]
                                / rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                     [k_darts_counter_temp4017][0];
                            (*(*q)) = 0.5
                                * (rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                      [k_darts_counter_temp4017][1]
                                        * rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                             [k_darts_counter_temp4017][1]
                                    + rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                         [k_darts_counter_temp4017][2]
                                        * rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                             [k_darts_counter_temp4017][2]
                                    + rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                         [k_darts_counter_temp4017][3]
                                        * rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                             [k_darts_counter_temp4017][3])
                                / rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                     [k_darts_counter_temp4017][0];
                            flux[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                [k_darts_counter_temp4017][1]
                                = rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                     [k_darts_counter_temp4017][1]
                                * (*(*u41));
                            flux[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                [k_darts_counter_temp4017][2]
                                = rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                     [k_darts_counter_temp4017][2]
                                * (*(*u41));
                            flux[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                [k_darts_counter_temp4017][3]
                                = rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                     [k_darts_counter_temp4017][3]
                                    * (*(*u41))
                                + 0.40000000000000002
                                    * (rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                          [k_darts_counter_temp4017][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                [k_darts_counter_temp4017][4]
                                = (1.3999999999999999
                                          * rsd[(i_darts_counter_temp4017)]
                                               [j_darts_counter_temp4017][k_darts_counter_temp4017]
                                               [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u41));
                        }
                        (*k) = k_darts_counter_temp4017;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp4017 = (*k);
                        for (; k_darts_counter_temp4017 <= nz - 2; k_darts_counter_temp4017++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp4017 = (*m);
                                for (; m_darts_counter_temp4017 < 5; m_darts_counter_temp4017++) {
                                    frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                        [k_darts_counter_temp4017][m_darts_counter_temp4017]
                                        = frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                              [k_darts_counter_temp4017][m_darts_counter_temp4017]
                                        - tz2
                                            * (flux[(i_darts_counter_temp4017)]
                                                   [j_darts_counter_temp4017]
                                                   [k_darts_counter_temp4017 + 1]
                                                   [m_darts_counter_temp4017]
                                                - flux[(i_darts_counter_temp4017)]
                                                      [j_darts_counter_temp4017]
                                                      [k_darts_counter_temp4017 - 1]
                                                      [m_darts_counter_temp4017]);
                                }
                                (*m) = m_darts_counter_temp4017;
                            }
                        }
                        (*k) = k_darts_counter_temp4017;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp4017 = (*k);
                        for (; k_darts_counter_temp4017 <= nz - 1; k_darts_counter_temp4017++) {
                            (*(*tmp)) = 1.
                                / rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                     [k_darts_counter_temp4017][0];
                            (*(*u21k)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                     [k_darts_counter_temp4017][1];
                            (*(*u31k)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                     [k_darts_counter_temp4017][2];
                            (*(*u41k)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                     [k_darts_counter_temp4017][3];
                            (*(*u51k)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                     [k_darts_counter_temp4017][4];
                            (*(*tmp)) = 1.
                                / rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                     [k_darts_counter_temp4017 - 1][0];
                            (*(*u21km1)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                     [k_darts_counter_temp4017 - 1][1];
                            (*(*u31km1)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                     [k_darts_counter_temp4017 - 1][2];
                            (*(*u41km1)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                     [k_darts_counter_temp4017 - 1][3];
                            (*(*u51km1)) = (*(*tmp))
                                * rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                     [k_darts_counter_temp4017 - 1][4];
                            flux[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                [k_darts_counter_temp4017][1]
                                = tz3 * ((*(*u21k)) - (*(*u21km1)));
                            flux[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                [k_darts_counter_temp4017][2]
                                = tz3 * ((*(*u31k)) - (*(*u31km1)));
                            flux[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                [k_darts_counter_temp4017][3]
                                = (4. / 3.) * tz3 * ((*(*u41k)) - (*(*u41km1)));
                            flux[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                [k_darts_counter_temp4017][4]
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
                        (*k) = k_darts_counter_temp4017;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp4017 = (*k);
                        for (; k_darts_counter_temp4017 <= nz - 2; k_darts_counter_temp4017++) {
                            frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                [k_darts_counter_temp4017][0]
                                = frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                      [k_darts_counter_temp4017][0]
                                + dz1 * tz1
                                    * (rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                          [k_darts_counter_temp4017 + 1][0]
                                        - 2.
                                            * rsd[(i_darts_counter_temp4017)]
                                                 [j_darts_counter_temp4017]
                                                 [k_darts_counter_temp4017][0]
                                        + rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                             [k_darts_counter_temp4017 - 1][0]);
                            frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                [k_darts_counter_temp4017][1]
                                = frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                      [k_darts_counter_temp4017][1]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                           [k_darts_counter_temp4017 + 1][1]
                                        - flux[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                              [k_darts_counter_temp4017][1])
                                + dz2 * tz1
                                    * (rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                          [k_darts_counter_temp4017 + 1][1]
                                        - 2.
                                            * rsd[(i_darts_counter_temp4017)]
                                                 [j_darts_counter_temp4017]
                                                 [k_darts_counter_temp4017][1]
                                        + rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                             [k_darts_counter_temp4017 - 1][1]);
                            frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                [k_darts_counter_temp4017][2]
                                = frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                      [k_darts_counter_temp4017][2]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                           [k_darts_counter_temp4017 + 1][2]
                                        - flux[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                              [k_darts_counter_temp4017][2])
                                + dz3 * tz1
                                    * (rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                          [k_darts_counter_temp4017 + 1][2]
                                        - 2.
                                            * rsd[(i_darts_counter_temp4017)]
                                                 [j_darts_counter_temp4017]
                                                 [k_darts_counter_temp4017][2]
                                        + rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                             [k_darts_counter_temp4017 - 1][2]);
                            frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                [k_darts_counter_temp4017][3]
                                = frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                      [k_darts_counter_temp4017][3]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                           [k_darts_counter_temp4017 + 1][3]
                                        - flux[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                              [k_darts_counter_temp4017][3])
                                + dz4 * tz1
                                    * (rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                          [k_darts_counter_temp4017 + 1][3]
                                        - 2.
                                            * rsd[(i_darts_counter_temp4017)]
                                                 [j_darts_counter_temp4017]
                                                 [k_darts_counter_temp4017][3]
                                        + rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                             [k_darts_counter_temp4017 - 1][3]);
                            frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                [k_darts_counter_temp4017][4]
                                = frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                      [k_darts_counter_temp4017][4]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                           [k_darts_counter_temp4017 + 1][4]
                                        - flux[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                              [k_darts_counter_temp4017][4])
                                + dz5 * tz1
                                    * (rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                          [k_darts_counter_temp4017 + 1][4]
                                        - 2.
                                            * rsd[(i_darts_counter_temp4017)]
                                                 [j_darts_counter_temp4017]
                                                 [k_darts_counter_temp4017][4]
                                        + rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                             [k_darts_counter_temp4017 - 1][4]);
                        }
                        (*k) = k_darts_counter_temp4017;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp4017 = (*m);
                        for (; m_darts_counter_temp4017 < 5; m_darts_counter_temp4017++) {
                            frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017][1]
                                [m_darts_counter_temp4017]
                                = frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017][1]
                                      [m_darts_counter_temp4017]
                                - (*(*dsspm))
                                    * (+5.
                                            * rsd[(i_darts_counter_temp4017)]
                                                 [j_darts_counter_temp4017][1]
                                                 [m_darts_counter_temp4017]
                                        - 4.
                                            * rsd[(i_darts_counter_temp4017)]
                                                 [j_darts_counter_temp4017][2]
                                                 [m_darts_counter_temp4017]
                                        + rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                             [3][m_darts_counter_temp4017]);
                            frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017][2]
                                [m_darts_counter_temp4017]
                                = frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017][2]
                                      [m_darts_counter_temp4017]
                                - (*(*dsspm))
                                    * (-4.
                                            * rsd[(i_darts_counter_temp4017)]
                                                 [j_darts_counter_temp4017][1]
                                                 [m_darts_counter_temp4017]
                                        + 6.
                                            * rsd[(i_darts_counter_temp4017)]
                                                 [j_darts_counter_temp4017][2]
                                                 [m_darts_counter_temp4017]
                                        - 4.
                                            * rsd[(i_darts_counter_temp4017)]
                                                 [j_darts_counter_temp4017][3]
                                                 [m_darts_counter_temp4017]
                                        + rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                             [4][m_darts_counter_temp4017]);
                        }
                        (*m) = m_darts_counter_temp4017;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 3;
                        int k_darts_counter_temp4017 = (*k);
                        for (; k_darts_counter_temp4017 <= nz - 4; k_darts_counter_temp4017++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp4017 = (*m);
                                for (; m_darts_counter_temp4017 < 5; m_darts_counter_temp4017++) {
                                    frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                        [k_darts_counter_temp4017][m_darts_counter_temp4017]
                                        = frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                              [k_darts_counter_temp4017][m_darts_counter_temp4017]
                                        - (*(*dsspm))
                                            * (rsd[(i_darts_counter_temp4017)]
                                                  [j_darts_counter_temp4017]
                                                  [k_darts_counter_temp4017 - 2]
                                                  [m_darts_counter_temp4017]
                                                - 4.
                                                    * rsd[(i_darts_counter_temp4017)]
                                                         [j_darts_counter_temp4017]
                                                         [k_darts_counter_temp4017 - 1]
                                                         [m_darts_counter_temp4017]
                                                + 6.
                                                    * rsd[(i_darts_counter_temp4017)]
                                                         [j_darts_counter_temp4017]
                                                         [k_darts_counter_temp4017]
                                                         [m_darts_counter_temp4017]
                                                - 4.
                                                    * rsd[(i_darts_counter_temp4017)]
                                                         [j_darts_counter_temp4017]
                                                         [k_darts_counter_temp4017 + 1]
                                                         [m_darts_counter_temp4017]
                                                + rsd[(i_darts_counter_temp4017)]
                                                     [j_darts_counter_temp4017]
                                                     [k_darts_counter_temp4017 + 2]
                                                     [m_darts_counter_temp4017]);
                                }
                                (*m) = m_darts_counter_temp4017;
                            }
                        }
                        (*k) = k_darts_counter_temp4017;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp4017 = (*m);
                        for (; m_darts_counter_temp4017 < 5; m_darts_counter_temp4017++) {
                            frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017][nz - 3]
                                [m_darts_counter_temp4017]
                                = frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017][nz - 3]
                                      [m_darts_counter_temp4017]
                                - (*(*dsspm))
                                    * (rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                          [nz - 5][m_darts_counter_temp4017]
                                        - 4.
                                            * rsd[(i_darts_counter_temp4017)]
                                                 [j_darts_counter_temp4017][nz - 4]
                                                 [m_darts_counter_temp4017]
                                        + 6.
                                            * rsd[(i_darts_counter_temp4017)]
                                                 [j_darts_counter_temp4017][nz - 3]
                                                 [m_darts_counter_temp4017]
                                        - 4.
                                            * rsd[(i_darts_counter_temp4017)]
                                                 [j_darts_counter_temp4017][nz - 2]
                                                 [m_darts_counter_temp4017]);
                            frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017][nz - 2]
                                [m_darts_counter_temp4017]
                                = frct[(i_darts_counter_temp4017)][j_darts_counter_temp4017][nz - 2]
                                      [m_darts_counter_temp4017]
                                - (*(*dsspm))
                                    * (rsd[(i_darts_counter_temp4017)][j_darts_counter_temp4017]
                                          [nz - 4][m_darts_counter_temp4017]
                                        - 4.
                                            * rsd[(i_darts_counter_temp4017)]
                                                 [j_darts_counter_temp4017][nz - 3]
                                                 [m_darts_counter_temp4017]
                                        + 5.
                                            * rsd[(i_darts_counter_temp4017)]
                                                 [j_darts_counter_temp4017][nz - 2]
                                                 [m_darts_counter_temp4017]);
                        }
                        (*m) = m_darts_counter_temp4017;
                    }
                }
                (*j) = j_darts_counter_temp4017;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets4017[0].decDep();
}
TP4017::TP4017(int in_numThreads, int in_mainCodeletID, TP2241* in_TPParent, int in_initIteration,
    int in_lastIteration, TP4017** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , dsspm_darts4017(new double*[this->numThreads])
    , i_darts4017(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts4017(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts4017(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts4017(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts4017(new double*[this->numThreads])
    , tmp_darts4017(new double*[this->numThreads])
    , u21k_darts4017(new double*[this->numThreads])
    , u21km1_darts4017(new double*[this->numThreads])
    , u31k_darts4017(new double*[this->numThreads])
    , u31km1_darts4017(new double*[this->numThreads])
    , u41_darts4017(new double*[this->numThreads])
    , u41k_darts4017(new double*[this->numThreads])
    , u41km1_darts4017(new double*[this->numThreads])
    , u51k_darts4017(new double*[this->numThreads])
    , u51km1_darts4017(new double*[this->numThreads])
    , initIteration4017(in_initIteration)
    , lastIteration4017(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets4017(new _barrierCodelets4017[1])
    , checkInCodelets4018(new _checkInCodelets4018[this->numThreads])
{
    /*Initialize the loop parameters*/
    range4017 = abs(lastIteration4017 - initIteration4017) / 1;
    rangePerCodelet4017 = range4017 / numThreads;
    minIteration4017 = min<int>(lastIteration4017, initIteration4017);
    remainderRange4017 = range4017 % numThreads;
    /*Initialize inputs and vars.*/
    this->dsspm_darts4017
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->q_darts4017 = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts4017
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21k_darts4017
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21km1_darts4017
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31k_darts4017
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31km1_darts4017
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41_darts4017
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41k_darts4017
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41km1_darts4017
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51k_darts4017
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51km1_darts4017
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets4017[0] = _barrierCodelets4017(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets4018* checkInCodelets4018Ptr = (this->checkInCodelets4018);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets4018);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets4018Ptr) = _checkInCodelets4018(2, 1, this, codeletCounter);
#else
        (*checkInCodelets4018Ptr) = _checkInCodelets4018(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets4018Ptr).decDep();
        checkInCodelets4018Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP4017::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets4018[localID].setID(codeletID);
    this->checkInCodelets4018[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets4018[localID + this->baseNumThreads * i]
            = _checkInCodelets4018(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets4018[localID + this->baseNumThreads * i]
            = _checkInCodelets4018(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets4018[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets4018[localID + this->baseNumThreads * i].decDep();
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
TP4017::~TP4017()
{
    delete[] dsspm_darts4017;
    delete[] q_darts4017;
    delete[] tmp_darts4017;
    delete[] u21k_darts4017;
    delete[] u21km1_darts4017;
    delete[] u31k_darts4017;
    delete[] u31km1_darts4017;
    delete[] u41_darts4017;
    delete[] u41k_darts4017;
    delete[] u41km1_darts4017;
    delete[] u51k_darts4017;
    delete[] u51km1_darts4017;
    delete[] barrierCodelets4017;
    delete[] checkInCodelets4018;
}
/*TP7: TP_jacld*/
void TP7::_checkInCodelets4922::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 4922: DeclStmt*/

    /*printing node 4923: DeclStmt*/

    /*printing node 4924: DeclStmt*/

    /*printing node 4925: DeclStmt*/

    /*printing node 4926: DeclStmt*/

    /*printing node 4927: BinaryOperator*/
    (this->inputsTPParent->r43_darts7[this->getID()]) = (4. / 3.);

    /*printing node 4931: BinaryOperator*/
    (this->inputsTPParent->c1345_darts7[this->getID()])
        = 1.3999999999999999 * 0.10000000000000001 * 1. * 1.3999999999999999;

    /*printing node 4939: BinaryOperator*/
    (this->inputsTPParent->c34_darts7[this->getID()]) = 0.10000000000000001 * 1.;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 4922 nextRegion: 4943 */
    myTP->controlTPParent->checkInCodelets4943[this->getID()].decDep();
}
void TP7::_checkInCodelets4943::fire(void)
{
    /*region 4943 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP4943;
    if (idx < myTP->TPsToUse4943) {
        if (!__sync_val_compare_and_swap(&(myTP->TP4943_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse4943;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse4943;
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
            if (idx == myTP->TPsToUse4943 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse4943 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP4943>(myTP, myTP->codeletsPerTP4943 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP4943Ptr[idx]));
#else
            place<TP4943>(idx, myTP, myTP->codeletsPerTP4943 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP4943Ptr[idx]));
#endif
        } else {
            if (myTP->TP4943Ptr[idx] != nullptr) {
                myTP->TP4943Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP7::_barrierCodelets4943::fire(void)
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
    , TP4943Ptr(new TP4943*[NUMTPS4943])
    , TP4943_alreadyLaunched(new size_t[NUMTPS4943])
    , numTPsSet4943(0)
    , numTPsReady4943(0)
    , TPsToUse4943(NUMTPS4943)
    , codeletsPerTP4943(this->numThreads / NUMTPS4943)
    , totalCodelets4943(this->TPsToUse4943 * this->codeletsPerTP4943)
    , checkInCodelets4922(new _checkInCodelets4922[this->numThreads])
    , checkInCodelets4943(new _checkInCodelets4943[this->numThreads])
    , barrierCodelets4943(new _barrierCodelets4943[1])
{
    barrierCodelets4943[0] = _barrierCodelets4943(NUMTPS4943, NUMTPS4943, this, 0);
    _checkInCodelets4943* checkInCodelets4943Ptr = (this->checkInCodelets4943);
    for (int i = 0; i < NUMTPS4943; i++) {
        TP4943Ptr[i] = nullptr;
        TP4943_alreadyLaunched[i] = 0;
    }
    _checkInCodelets4922* checkInCodelets4922Ptr = (this->checkInCodelets4922);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets4922);
#endif
    for (size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++) {
        (*checkInCodelets4943Ptr) = _checkInCodelets4943(1, 1, this, codeletCounter);
        checkInCodelets4943Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets4922Ptr) = _checkInCodelets4922(2, 1, this, codeletCounter);
#else
        (*checkInCodelets4922Ptr) = _checkInCodelets4922(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets4922Ptr).decDep();
        checkInCodelets4922Ptr++;
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
    delete[] barrierCodelets4943;
    delete[] checkInCodelets4943;
    delete[] checkInCodelets4922;
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
/*TP4943: OMPForDirective*/
void TP4943::_barrierCodelets4943::fire(void)
{
    TP4943* myTP = static_cast<TP4943*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets4943[0].decDep();
}
bool TP4943::requestNewRangeIterations4943(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 1*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet4943 * codeletID;
        int tempEndRange = rangePerCodelet4943 * (codeletID + 1);
        if (remainderRange4943 != 0) {
            if (codeletID < (uint32_t)remainderRange4943) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange4943;
                tempEndRange += remainderRange4943;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration4943;
        tempEndRange = tempEndRange * 1 + minIteration4943;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration4943 < lastIteration4943) {
            (this->inputsTPParent->i_darts4943[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts4943[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration4943;
        }
    }
    return isThereNewIteration;
}
void TP4943::_checkInCodelets4944::fire(void)
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
    this->inputsTPParent->c1345_darts4943[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->c1345_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->c34_darts4943[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->c34_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->k_darts4943[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->k_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->r43_darts4943[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->r43_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp1_darts4943[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp1_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp2_darts4943[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp2_darts7[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp3_darts4943[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp3_darts7[this->getID()]);

    /*printing node 4944: ForStmt*/
    /*var: c1345*/
    /*var: c34*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: r43*/
    /*var: tmp1*/
    /*var: tmp2*/
    /*var: tmp3*/
    double** c1345 = &(this->inputsTPParent->c1345_darts4943[this->getLocalID()]);
    (void)c1345 /*OMP_SHARED_PRIVATE*/;
    double** c34 = &(this->inputsTPParent->c34_darts4943[this->getLocalID()]);
    (void)c34 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts4943[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts4943[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** k = &(this->inputsTPParent->k_darts4943[this->getLocalID()]);
    (void)k /*OMP_SHARED_PRIVATE*/;
    double** r43 = &(this->inputsTPParent->r43_darts4943[this->getLocalID()]);
    (void)r43 /*OMP_SHARED_PRIVATE*/;
    double** tmp1 = &(this->inputsTPParent->tmp1_darts4943[this->getLocalID()]);
    (void)tmp1 /*OMP_SHARED_PRIVATE*/;
    double** tmp2 = &(this->inputsTPParent->tmp2_darts4943[this->getLocalID()]);
    (void)tmp2 /*OMP_SHARED_PRIVATE*/;
    double** tmp3 = &(this->inputsTPParent->tmp3_darts4943[this->getLocalID()]);
    (void)tmp3 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations4943(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets4943[0].decDep();
        return;
    }
    for (int i_darts_counter_temp4943 = (*i); i_darts_counter_temp4943 <= endRange
         && i_darts_counter_temp4943 <= this->inputsTPParent->lastIteration4943;
         i_darts_counter_temp4943++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp4943 = (*j);
                for (; j_darts_counter_temp4943 <= jend; j_darts_counter_temp4943++) {
                    (*(*tmp1))
                        = 1. / u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][0]
                        = 1. + dt * 2. * (tx1 * dx1 + ty1 * dy1 + tz1 * dz1);
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][1] = 0.;
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][2] = 0.;
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][3] = 0.;
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][4] = 0.;
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][0] = dt * 2.
                        * (tx1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                       [(*(*k))][1])
                            + ty1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                       [(*(*k))][1])
                            + tz1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                       [(*(*k))][1]));
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][1] = 1.
                        + dt * 2.
                            * (tx1 * (*(*r43)) * (*(*c34)) * (*(*tmp1))
                                + ty1 * (*(*c34)) * (*(*tmp1)) + tz1 * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx2 + ty1 * dy2 + tz1 * dz2);
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][2] = 0.;
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][3] = 0.;
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][4] = 0.;
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][0] = dt * 2.
                        * (tx1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                       [(*(*k))][2])
                            + ty1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                       [(*(*k))][2])
                            + tz1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                       [(*(*k))][2]));
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][1] = 0.;
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][2] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c34)) * (*(*tmp1))
                                + ty1 * (*(*r43)) * (*(*c34)) * (*(*tmp1))
                                + tz1 * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx3 + ty1 * dy3 + tz1 * dz3);
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][3] = 0.;
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][4] = 0.;
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][0] = dt * 2.
                        * (tx1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                       [(*(*k))][3])
                            + ty1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                       [(*(*k))][3])
                            + tz1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                       [(*(*k))][3]));
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][1] = 0.;
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][2] = 0.;
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][3] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c34)) * (*(*tmp1)) + ty1 * (*(*c34)) * (*(*tmp1))
                                + tz1 * (*(*r43)) * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx4 + ty1 * dy4 + tz1 * dz4);
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][4] = 0.;
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][0] = dt * 2.
                        * (tx1
                                * (-((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp4943)]
                                                [j_darts_counter_temp4943][(*(*k))][1])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp4943)]
                                                [j_darts_counter_temp4943][(*(*k))][2])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp4943)]
                                                [j_darts_counter_temp4943][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                           [(*(*k))][4])
                            + ty1
                                * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp4943)]
                                                [j_darts_counter_temp4943][(*(*k))][1])))
                                    - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp4943)]
                                                [j_darts_counter_temp4943][(*(*k))][2])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp4943)]
                                                [j_darts_counter_temp4943][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                           [(*(*k))][4])
                            + tz1
                                * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp4943)]
                                                [j_darts_counter_temp4943][(*(*k))][1])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp4943)]
                                                [j_darts_counter_temp4943][(*(*k))][2])))
                                    - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp4943)]
                                                [j_darts_counter_temp4943][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                           [(*(*k))][4]));
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][1] = dt * 2.
                        * (tx1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k))]
                                   [1]
                            + ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k))]
                                   [1]
                            + tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k))]
                                   [1]);
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][2] = dt * 2.
                        * (tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k))]
                                   [2]
                            + ty1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k))]
                                   [2]
                            + tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k))]
                                   [2]);
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][3] = dt * 2.
                        * (tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k))]
                                   [3]
                            + ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k))]
                                   [3]
                            + tz1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k))]
                                   [3]);
                    d[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][4] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c1345)) * (*(*tmp1)) + ty1 * (*(*c1345)) * (*(*tmp1))
                                + tz1 * (*(*c1345)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx5 + ty1 * dy5 + tz1 * dz5);
                    (*(*tmp1)) = 1.
                        / u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k)) - 1][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][0] = -dt * tz1 * dz1;
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][1] = 0.;
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][2] = 0.;
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][3] = -dt * tz2;
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][4] = 0.;
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][0] = -dt * tz2
                            * (-(u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                  [(*(*k)) - 1][1]
                                   * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                      [(*(*k)) - 1][3])
                                * (*(*tmp2)))
                        - dt * tz1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                   [(*(*k)) - 1][1]);
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][1] = -dt * tz2
                            * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k)) - 1]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * (*(*c34)) * (*(*tmp1)) - dt * tz1 * dz2;
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][2] = 0.;
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][3] = -dt * tz2
                        * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k)) - 1][1]
                            * (*(*tmp1)));
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][4] = 0.;
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][0] = -dt * tz2
                            * (-(u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                  [(*(*k)) - 1][2]
                                   * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                      [(*(*k)) - 1][3])
                                * (*(*tmp2)))
                        - dt * tz1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                   [(*(*k)) - 1][2]);
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][1] = 0.;
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][2] = -dt * tz2
                            * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k)) - 1]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * ((*(*c34)) * (*(*tmp1))) - dt * tz1 * dz3;
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][3] = -dt * tz2
                        * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k)) - 1][2]
                            * (*(*tmp1)));
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][4] = 0.;
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][0] = -dt * tz2
                            * (-(u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                  [(*(*k)) - 1][3]
                                   * (*(*tmp1)))
                                    * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                        [(*(*k)) - 1][3]
                                        * (*(*tmp1)))
                                + 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                         [(*(*k)) - 1][1]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943][(*(*k)) - 1][1]
                                           + u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                              [(*(*k)) - 1][2]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943][(*(*k)) - 1][2]
                                           + u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                              [(*(*k)) - 1][3]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943][(*(*k)) - 1][3])
                                        * (*(*tmp2))))
                        - dt * tz1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                   [(*(*k)) - 1][3]);
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][1] = -dt * tz2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k)) - 1]
                                [1]
                                * (*(*tmp1))));
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][2] = -dt * tz2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k)) - 1]
                                [2]
                                * (*(*tmp1))));
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][3] = -dt * tz2
                            * (2. - 0.40000000000000002)
                            * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k)) - 1]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * tz1 * dz4;
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][4]
                        = -dt * tz2 * 0.40000000000000002;
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][0] = -dt * tz2
                            * ((0.40000000000000002
                                       * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                           [(*(*k)) - 1][1]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943][(*(*k)) - 1][1]
                                           + u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                              [(*(*k)) - 1][2]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943][(*(*k)) - 1][2]
                                           + u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                              [(*(*k)) - 1][3]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943][(*(*k)) - 1][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                           [(*(*k)) - 1][4]
                                           * (*(*tmp1))))
                                * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                    [(*(*k)) - 1][3]
                                    * (*(*tmp1))))
                        - dt * tz1
                            * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                        [(*(*k)) - 1][1]
                                        * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                           [(*(*k)) - 1][1])
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                        [(*(*k)) - 1][2]
                                        * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                           [(*(*k)) - 1][2])
                                - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                        [(*(*k)) - 1][3]
                                        * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                           [(*(*k)) - 1][3])
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                       [(*(*k)) - 1][4]);
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][1] = -dt * tz2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                    [(*(*k)) - 1][1]
                                    * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                       [(*(*k)) - 1][3])
                                * (*(*tmp2)))
                        - dt * tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k)) - 1]
                               [1];
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][2] = -dt * tz2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                    [(*(*k)) - 1][2]
                                    * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                       [(*(*k)) - 1][3])
                                * (*(*tmp2)))
                        - dt * tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k)) - 1]
                               [2];
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][3] = -dt * tz2
                            * (1.3999999999999999
                                    * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                        [(*(*k)) - 1][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                         [(*(*k)) - 1][1]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943][(*(*k)) - 1][1]
                                           + u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                              [(*(*k)) - 1][2]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943][(*(*k)) - 1][2]
                                           + 3.
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943][(*(*k)) - 1][3]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943][(*(*k)) - 1][3])
                                        * (*(*tmp2))))
                        - dt * tz1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943][(*(*k)) - 1]
                               [3];
                    a[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][4] = -dt * tz2
                            * (1.3999999999999999
                                * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943]
                                    [(*(*k)) - 1][3]
                                    * (*(*tmp1))))
                        - dt * tz1 * (*(*c1345)) * (*(*tmp1)) - dt * tz1 * dz5;
                    (*(*tmp1)) = 1.
                        / u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][0] = -dt * ty1 * dy1;
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][1] = 0.;
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][2] = -dt * ty2;
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][3] = 0.;
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][4] = 0.;
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][0] = -dt * ty2
                            * (-(u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                  [(*(*k))][1]
                                   * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                      [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * ty1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                   [(*(*k))][1]);
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][1] = -dt * ty2
                            * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1][(*(*k))]
                                [2]
                                * (*(*tmp1)))
                        - dt * ty1 * ((*(*c34)) * (*(*tmp1))) - dt * ty1 * dy2;
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][2] = -dt * ty2
                        * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1][(*(*k))][1]
                            * (*(*tmp1)));
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][3] = 0.;
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][4] = 0.;
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][0] = -dt * ty2
                            * (-(u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                  [(*(*k))][2]
                                   * (*(*tmp1)))
                                    * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                        [(*(*k))][2]
                                        * (*(*tmp1)))
                                + 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                         [(*(*k))][1]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943 - 1][(*(*k))][1]
                                           + u[(i_darts_counter_temp4943)]
                                              [j_darts_counter_temp4943 - 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943 - 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp4943)]
                                              [j_darts_counter_temp4943 - 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943 - 1][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * ty1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                   [(*(*k))][2]);
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][1] = -dt * ty2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1][(*(*k))]
                                [1]
                                * (*(*tmp1))));
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][2] = -dt * ty2
                            * ((2. - 0.40000000000000002)
                                * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * ty1 * dy3;
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][3] = -dt * ty2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1][(*(*k))]
                                [3]
                                * (*(*tmp1))));
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][4]
                        = -dt * ty2 * 0.40000000000000002;
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][0] = -dt * ty2
                            * (-(u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                  [(*(*k))][2]
                                   * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                      [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * ty1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                   [(*(*k))][3]);
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][1] = 0.;
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][2] = -dt * ty2
                        * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1][(*(*k))][3]
                            * (*(*tmp1)));
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][3] = -dt * ty2
                            * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1][(*(*k))]
                                [2]
                                * (*(*tmp1)))
                        - dt * ty1 * ((*(*c34)) * (*(*tmp1))) - dt * ty1 * dy4;
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][4] = 0.;
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][0] = -dt * ty2
                            * ((0.40000000000000002
                                       * (u[(i_darts_counter_temp4943)]
                                           [j_darts_counter_temp4943 - 1][(*(*k))][1]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943 - 1][(*(*k))][1]
                                           + u[(i_darts_counter_temp4943)]
                                              [j_darts_counter_temp4943 - 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943 - 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp4943)]
                                              [j_darts_counter_temp4943 - 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943 - 1][(*(*k))][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(i_darts_counter_temp4943)]
                                           [j_darts_counter_temp4943 - 1][(*(*k))][4]
                                           * (*(*tmp1))))
                                * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1
                            * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                          [(*(*k))][1])
                                        * (u[(i_darts_counter_temp4943)]
                                            [j_darts_counter_temp4943 - 1][(*(*k))][1])))
                                - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                          [(*(*k))][2])
                                        * (u[(i_darts_counter_temp4943)]
                                            [j_darts_counter_temp4943 - 1][(*(*k))][2])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                          [(*(*k))][3])
                                        * (u[(i_darts_counter_temp4943)]
                                            [j_darts_counter_temp4943 - 1][(*(*k))][3])))
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                       [(*(*k))][4]);
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][1] = -dt * ty2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                    [(*(*k))][1]
                                    * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                       [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1][(*(*k))]
                               [1];
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][2] = -dt * ty2
                            * (1.3999999999999999
                                    * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                        [(*(*k))][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                         [(*(*k))][1]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943 - 1][(*(*k))][1]
                                           + 3.
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943 - 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943 - 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp4943)]
                                              [j_darts_counter_temp4943 - 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp4943)]
                                                  [j_darts_counter_temp4943 - 1][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * ty1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1][(*(*k))]
                               [2];
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][3] = -dt * ty2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                    [(*(*k))][2]
                                    * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                       [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1][(*(*k))]
                               [3];
                    b[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][4] = -dt * ty2
                            * (1.3999999999999999
                                * (u[(i_darts_counter_temp4943)][j_darts_counter_temp4943 - 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1 * (*(*c1345)) * (*(*tmp1)) - dt * ty1 * dy5;
                    (*(*tmp1)) = 1.
                        / u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][0] = -dt * tx1 * dx1;
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][1] = -dt * tx2;
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][2] = 0.;
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][3] = 0.;
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][0][4] = 0.;
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][0] = -dt * tx2
                            * (-(u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943][(*(*k))]
                                  [1]
                                   * (*(*tmp1)))
                                    * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                        [(*(*k))][1]
                                        * (*(*tmp1)))
                                + 0.40000000000000002 * 0.5
                                    * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                        [(*(*k))][1]
                                            * u[(i_darts_counter_temp4943)-1]
                                               [j_darts_counter_temp4943][(*(*k))][1]
                                        + u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                           [(*(*k))][2]
                                            * u[(i_darts_counter_temp4943)-1]
                                               [j_darts_counter_temp4943][(*(*k))][2]
                                        + u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                           [(*(*k))][3]
                                            * u[(i_darts_counter_temp4943)-1]
                                               [j_darts_counter_temp4943][(*(*k))][3])
                                    * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943][(*(*k))]
                                   [1]);
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][1] = -dt * tx2
                            * ((2. - 0.40000000000000002)
                                * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * tx1 * dx2;
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][2] = -dt * tx2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943][(*(*k))][2]
                                * (*(*tmp1))));
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][3] = -dt * tx2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943][(*(*k))][3]
                                * (*(*tmp1))));
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][1][4]
                        = -dt * tx2 * 0.40000000000000002;
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][0] = -dt * tx2
                            * (-(u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943][(*(*k))]
                                  [1]
                                   * u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                      [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943][(*(*k))]
                                   [2]);
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][1] = -dt * tx2
                        * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943][(*(*k))][2]
                            * (*(*tmp1)));
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][2] = -dt * tx2
                            * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943][(*(*k))][1]
                                * (*(*tmp1)))
                        - dt * tx1 * ((*(*c34)) * (*(*tmp1))) - dt * tx1 * dx3;
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][3] = 0.;
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][2][4] = 0.;
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][0] = -dt * tx2
                            * (-(u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943][(*(*k))]
                                  [1]
                                   * u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                      [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943][(*(*k))]
                                   [3]);
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][1] = -dt * tx2
                        * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943][(*(*k))][3]
                            * (*(*tmp1)));
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][2] = 0.;
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][3] = -dt * tx2
                            * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943][(*(*k))][1]
                                * (*(*tmp1)))
                        - dt * tx1 * ((*(*c34)) * (*(*tmp1))) - dt * tx1 * dx4;
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][3][4] = 0.;
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][0] = -dt * tx2
                            * ((0.40000000000000002
                                       * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                           [(*(*k))][1]
                                               * u[(i_darts_counter_temp4943)-1]
                                                  [j_darts_counter_temp4943][(*(*k))][1]
                                           + u[(i_darts_counter_temp4943)-1]
                                              [j_darts_counter_temp4943][(*(*k))][2]
                                               * u[(i_darts_counter_temp4943)-1]
                                                  [j_darts_counter_temp4943][(*(*k))][2]
                                           + u[(i_darts_counter_temp4943)-1]
                                              [j_darts_counter_temp4943][(*(*k))][3]
                                               * u[(i_darts_counter_temp4943)-1]
                                                  [j_darts_counter_temp4943][(*(*k))][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                           [(*(*k))][4]
                                           * (*(*tmp1))))
                                * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1
                            * (-((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                          [(*(*k))][1])
                                        * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                            [(*(*k))][1])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                          [(*(*k))][2])
                                        * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                            [(*(*k))][2])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                          [(*(*k))][3])
                                        * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                            [(*(*k))][3])))
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                       [(*(*k))][4]);
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][1] = -dt * tx2
                            * (1.3999999999999999
                                    * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                        [(*(*k))][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((3.
                                               * u[(i_darts_counter_temp4943)-1]
                                                  [j_darts_counter_temp4943][(*(*k))][1]
                                               * u[(i_darts_counter_temp4943)-1]
                                                  [j_darts_counter_temp4943][(*(*k))][1]
                                           + u[(i_darts_counter_temp4943)-1]
                                              [j_darts_counter_temp4943][(*(*k))][2]
                                               * u[(i_darts_counter_temp4943)-1]
                                                  [j_darts_counter_temp4943][(*(*k))][2]
                                           + u[(i_darts_counter_temp4943)-1]
                                              [j_darts_counter_temp4943][(*(*k))][3]
                                               * u[(i_darts_counter_temp4943)-1]
                                                  [j_darts_counter_temp4943][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * tx1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943][(*(*k))][1];
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][2] = -dt * tx2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                    [(*(*k))][2]
                                    * u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                       [(*(*k))][1])
                                * (*(*tmp2)))
                        - dt * tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943][(*(*k))][2];
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][3] = -dt * tx2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                    [(*(*k))][3]
                                    * u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                       [(*(*k))][1])
                                * (*(*tmp2)))
                        - dt * tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943][(*(*k))][3];
                    c[(i_darts_counter_temp4943)][j_darts_counter_temp4943][4][4] = -dt * tx2
                            * (1.3999999999999999
                                * (u[(i_darts_counter_temp4943)-1][j_darts_counter_temp4943]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1 * (*(*c1345)) * (*(*tmp1)) - dt * tx1 * dx5;
                }
                (*j) = j_darts_counter_temp4943;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets4943[0].decDep();
}
TP4943::TP4943(int in_numThreads, int in_mainCodeletID, TP7* in_TPParent, int in_initIteration,
    int in_lastIteration, TP4943** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , c1345_darts4943(new double*[this->numThreads])
    , c34_darts4943(new double*[this->numThreads])
    , i_darts4943(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts4943(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts4943(new int*[this->numThreads])
    , r43_darts4943(new double*[this->numThreads])
    , tmp1_darts4943(new double*[this->numThreads])
    , tmp2_darts4943(new double*[this->numThreads])
    , tmp3_darts4943(new double*[this->numThreads])
    , initIteration4943(in_initIteration)
    , lastIteration4943(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets4943(new _barrierCodelets4943[1])
    , checkInCodelets4944(new _checkInCodelets4944[this->numThreads])
{
    /*Initialize the loop parameters*/
    range4943 = abs(lastIteration4943 - initIteration4943) / 1;
    rangePerCodelet4943 = range4943 / numThreads;
    minIteration4943 = min<int>(lastIteration4943, initIteration4943);
    remainderRange4943 = range4943 % numThreads;
    /*Initialize inputs and vars.*/
    this->c1345_darts4943
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->c34_darts4943
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->k_darts4943 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->r43_darts4943
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp1_darts4943
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp2_darts4943
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp3_darts4943
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets4943[0] = _barrierCodelets4943(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets4944* checkInCodelets4944Ptr = (this->checkInCodelets4944);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets4944);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets4944Ptr) = _checkInCodelets4944(2, 1, this, codeletCounter);
#else
        (*checkInCodelets4944Ptr) = _checkInCodelets4944(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets4944Ptr).decDep();
        checkInCodelets4944Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP4943::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets4944[localID].setID(codeletID);
    this->checkInCodelets4944[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets4944[localID + this->baseNumThreads * i]
            = _checkInCodelets4944(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets4944[localID + this->baseNumThreads * i]
            = _checkInCodelets4944(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets4944[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets4944[localID + this->baseNumThreads * i].decDep();
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
TP4943::~TP4943()
{
    delete[] c1345_darts4943;
    delete[] c34_darts4943;
    delete[] k_darts4943;
    delete[] r43_darts4943;
    delete[] tmp1_darts4943;
    delete[] tmp2_darts4943;
    delete[] tmp3_darts4943;
    delete[] barrierCodelets4943;
    delete[] checkInCodelets4944;
}
/*TP8: TP_jacu*/
void TP8::_checkInCodelets7441::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 7441: DeclStmt*/

    /*printing node 7442: DeclStmt*/

    /*printing node 7443: DeclStmt*/

    /*printing node 7444: DeclStmt*/

    /*printing node 7445: DeclStmt*/

    /*printing node 7446: BinaryOperator*/
    (this->inputsTPParent->r43_darts8[this->getID()]) = (4. / 3.);

    /*printing node 7450: BinaryOperator*/
    (this->inputsTPParent->c1345_darts8[this->getID()])
        = 1.3999999999999999 * 0.10000000000000001 * 1. * 1.3999999999999999;

    /*printing node 7458: BinaryOperator*/
    (this->inputsTPParent->c34_darts8[this->getID()]) = 0.10000000000000001 * 1.;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 7441 nextRegion: 7462 */
    myTP->controlTPParent->checkInCodelets7462[this->getID()].decDep();
}
void TP8::_checkInCodelets7462::fire(void)
{
    /*region 7462 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP7462;
    if (idx < myTP->TPsToUse7462) {
        if (!__sync_val_compare_and_swap(&(myTP->TP7462_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(ist - iend) / 1;
            int rangePerCodelet = range / myTP->TPsToUse7462;
            int minIteration = min<int>(ist, iend);
            int remainderRange = range % myTP->TPsToUse7462;
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
            if (idx == myTP->TPsToUse7462 - 1) {
                lastIteration = ist;
            }
#if USEINVOKE == 1
            invoke<TP7462>(myTP, myTP->codeletsPerTP7462 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(myTP->TP7462Ptr[idx]));
#else
            place<TP7462>(idx, myTP, myTP->codeletsPerTP7462 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP7462Ptr[idx]));
#endif
        } else {
            if (myTP->TP7462Ptr[idx] != nullptr) {
                myTP->TP7462Ptr[idx]->dispatchCodelet(this->getID());
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
    , TP7462Ptr(new TP7462*[NUMTPS7462])
    , TP7462_alreadyLaunched(new size_t[NUMTPS7462])
    , numTPsSet7462(0)
    , numTPsReady7462(0)
    , TPsToUse7462(NUMTPS7462)
    , codeletsPerTP7462(this->numThreads / NUMTPS7462)
    , totalCodelets7462(this->TPsToUse7462 * this->codeletsPerTP7462)
    , checkInCodelets7441(new _checkInCodelets7441[this->numThreads])
    , checkInCodelets7462(new _checkInCodelets7462[this->numThreads])
{
    _checkInCodelets7462* checkInCodelets7462Ptr = (this->checkInCodelets7462);
    for (int i = 0; i < NUMTPS7462; i++) {
        TP7462Ptr[i] = nullptr;
        TP7462_alreadyLaunched[i] = 0;
    }
    _checkInCodelets7441* checkInCodelets7441Ptr = (this->checkInCodelets7441);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets7441);
#endif
    for (size_t codeletCounter = 0; codeletCounter < this->numThreads; codeletCounter++) {
        (*checkInCodelets7462Ptr) = _checkInCodelets7462(1, 1, this, codeletCounter);
        checkInCodelets7462Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets7441Ptr) = _checkInCodelets7441(2, 1, this, codeletCounter);
#else
        (*checkInCodelets7441Ptr) = _checkInCodelets7441(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets7441Ptr).decDep();
        checkInCodelets7441Ptr++;
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
    delete[] checkInCodelets7462;
    delete[] checkInCodelets7441;
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
/*TP7462: OMPForDirective*/
bool TP7462::requestNewRangeIterations7462(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 1*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet7462 * codeletID;
        int tempEndRange = rangePerCodelet7462 * (codeletID + 1);
        if (remainderRange7462 != 0) {
            if (codeletID < (uint32_t)remainderRange7462) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange7462;
                tempEndRange += remainderRange7462;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration7462;
        tempEndRange = tempEndRange * 1 + minIteration7462;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration7462 < lastIteration7462) {
            (this->inputsTPParent->i_darts7462[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts7462[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == 0) {
            *endRange = *endRange - 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration7462;
        }
    }
    return isThereNewIteration;
}
void TP7462::_checkInCodelets7463::fire(void)
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
    this->inputsTPParent->c1345_darts7462[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->c1345_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->c34_darts7462[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->c34_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->k_darts7462[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->k_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->r43_darts7462[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->r43_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp1_darts7462[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp1_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp2_darts7462[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp2_darts8[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp3_darts7462[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp3_darts8[this->getID()]);

    /*printing node 7463: ForStmt*/
    /*var: c1345*/
    /*var: c34*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: r43*/
    /*var: tmp1*/
    /*var: tmp2*/
    /*var: tmp3*/
    double** c1345 = &(this->inputsTPParent->c1345_darts7462[this->getLocalID()]);
    (void)c1345 /*OMP_SHARED_PRIVATE*/;
    double** c34 = &(this->inputsTPParent->c34_darts7462[this->getLocalID()]);
    (void)c34 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts7462[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts7462[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** k = &(this->inputsTPParent->k_darts7462[this->getLocalID()]);
    (void)k /*OMP_SHARED_PRIVATE*/;
    double** r43 = &(this->inputsTPParent->r43_darts7462[this->getLocalID()]);
    (void)r43 /*OMP_SHARED_PRIVATE*/;
    double** tmp1 = &(this->inputsTPParent->tmp1_darts7462[this->getLocalID()]);
    (void)tmp1 /*OMP_SHARED_PRIVATE*/;
    double** tmp2 = &(this->inputsTPParent->tmp2_darts7462[this->getLocalID()]);
    (void)tmp2 /*OMP_SHARED_PRIVATE*/;
    double** tmp3 = &(this->inputsTPParent->tmp3_darts7462[this->getLocalID()]);
    (void)tmp3 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations7462(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Find and signal the next codelet*/
        return;
    }
    for (int i_darts_counter_temp7462 = (*i); i_darts_counter_temp7462 >= endRange
         && i_darts_counter_temp7462 >= this->inputsTPParent->lastIteration7462;
         i_darts_counter_temp7462--) {
        {
            {
                /*Loop's init*/
                (*j) = jend;
                int j_darts_counter_temp7462 = (*j);
                for (; j_darts_counter_temp7462 >= jst; j_darts_counter_temp7462--) {
                    (*(*tmp1))
                        = 1. / u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][0]
                        = 1. + dt * 2. * (tx1 * dx1 + ty1 * dy1 + tz1 * dz1);
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][1] = 0.;
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][2] = 0.;
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][3] = 0.;
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][4] = 0.;
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][0] = dt * 2.
                        * (tx1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                       [(*(*k))][1])
                            + ty1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                       [(*(*k))][1])
                            + tz1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                       [(*(*k))][1]));
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][1] = 1.
                        + dt * 2.
                            * (tx1 * (*(*r43)) * (*(*c34)) * (*(*tmp1))
                                + ty1 * (*(*c34)) * (*(*tmp1)) + tz1 * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx2 + ty1 * dy2 + tz1 * dz2);
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][2] = 0.;
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][3] = 0.;
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][4] = 0.;
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][0] = dt * 2.
                        * (tx1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                       [(*(*k))][2])
                            + ty1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                       [(*(*k))][2])
                            + tz1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                       [(*(*k))][2]));
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][1] = 0.;
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][2] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c34)) * (*(*tmp1))
                                + ty1 * (*(*r43)) * (*(*c34)) * (*(*tmp1))
                                + tz1 * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx3 + ty1 * dy3 + tz1 * dz3);
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][3] = 0.;
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][4] = 0.;
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][0] = dt * 2.
                        * (tx1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                       [(*(*k))][3])
                            + ty1
                                * (-(*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                       [(*(*k))][3])
                            + tz1
                                * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                       [(*(*k))][3]));
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][1] = 0.;
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][2] = 0.;
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][3] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c34)) * (*(*tmp1)) + ty1 * (*(*c34)) * (*(*tmp1))
                                + tz1 * (*(*r43)) * (*(*c34)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx4 + ty1 * dy4 + tz1 * dz4);
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][4] = 0.;
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][0] = dt * 2.
                        * (tx1
                                * (-((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp7462)]
                                                [j_darts_counter_temp7462][(*(*k))][1])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp7462)]
                                                [j_darts_counter_temp7462][(*(*k))][2])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp7462)]
                                                [j_darts_counter_temp7462][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                           [(*(*k))][4])
                            + ty1
                                * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp7462)]
                                                [j_darts_counter_temp7462][(*(*k))][1])))
                                    - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp7462)]
                                                [j_darts_counter_temp7462][(*(*k))][2])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp7462)]
                                                [j_darts_counter_temp7462][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                           [(*(*k))][4])
                            + tz1
                                * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                              [(*(*k))][1])
                                            * (u[(i_darts_counter_temp7462)]
                                                [j_darts_counter_temp7462][(*(*k))][1])))
                                    - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                              [(*(*k))][2])
                                            * (u[(i_darts_counter_temp7462)]
                                                [j_darts_counter_temp7462][(*(*k))][2])))
                                    - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                        * (((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                              [(*(*k))][3])
                                            * (u[(i_darts_counter_temp7462)]
                                                [j_darts_counter_temp7462][(*(*k))][3])))
                                    - ((*(*c1345))) * (*(*tmp2))
                                        * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                           [(*(*k))][4]));
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][1] = dt * 2.
                        * (tx1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k))]
                                   [1]
                            + ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k))]
                                   [1]
                            + tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k))]
                                   [1]);
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][2] = dt * 2.
                        * (tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k))]
                                   [2]
                            + ty1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k))]
                                   [2]
                            + tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k))]
                                   [2]);
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][3] = dt * 2.
                        * (tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k))]
                                   [3]
                            + ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k))]
                                   [3]
                            + tz1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k))]
                                   [3]);
                    d[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][4] = 1.
                        + dt * 2.
                            * (tx1 * (*(*c1345)) * (*(*tmp1)) + ty1 * (*(*c1345)) * (*(*tmp1))
                                + tz1 * (*(*c1345)) * (*(*tmp1)))
                        + dt * 2. * (tx1 * dx5 + ty1 * dy5 + tz1 * dz5);
                    (*(*tmp1)) = 1.
                        / u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][0] = -dt * tx1 * dx1;
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][1] = dt * tx2;
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][2] = 0.;
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][3] = 0.;
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][4] = 0.;
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][0] = dt * tx2
                            * (-(u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                  [(*(*k))][1]
                                   * (*(*tmp1)))
                                    * (u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                        [(*(*k))][1]
                                        * (*(*tmp1)))
                                + 0.40000000000000002 * 0.5
                                    * (u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                        [(*(*k))][1]
                                            * u[(i_darts_counter_temp7462) + 1]
                                               [j_darts_counter_temp7462][(*(*k))][1]
                                        + u[(i_darts_counter_temp7462) + 1]
                                           [j_darts_counter_temp7462][(*(*k))][2]
                                            * u[(i_darts_counter_temp7462) + 1]
                                               [j_darts_counter_temp7462][(*(*k))][2]
                                        + u[(i_darts_counter_temp7462) + 1]
                                           [j_darts_counter_temp7462][(*(*k))][3]
                                            * u[(i_darts_counter_temp7462) + 1]
                                               [j_darts_counter_temp7462][(*(*k))][3])
                                    * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                   [(*(*k))][1]);
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][1] = dt * tx2
                            * ((2. - 0.40000000000000002)
                                * (u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * tx1 * dx2;
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][2] = dt * tx2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462][(*(*k))]
                                [2]
                                * (*(*tmp1))));
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][3] = dt * tx2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462][(*(*k))]
                                [3]
                                * (*(*tmp1))));
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][4]
                        = dt * tx2 * 0.40000000000000002;
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][0] = dt * tx2
                            * (-(u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                  [(*(*k))][1]
                                   * u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                      [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                   [(*(*k))][2]);
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][1] = dt * tx2
                        * (u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462][(*(*k))][2]
                            * (*(*tmp1)));
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][2] = dt * tx2
                            * (u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462][(*(*k))]
                                [1]
                                * (*(*tmp1)))
                        - dt * tx1 * ((*(*c34)) * (*(*tmp1))) - dt * tx1 * dx3;
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][3] = 0.;
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][4] = 0.;
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][0] = dt * tx2
                            * (-(u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                  [(*(*k))][1]
                                   * u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                      [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * tx1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                   [(*(*k))][3]);
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][1] = dt * tx2
                        * (u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462][(*(*k))][3]
                            * (*(*tmp1)));
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][2] = 0.;
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][3] = dt * tx2
                            * (u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462][(*(*k))]
                                [1]
                                * (*(*tmp1)))
                        - dt * tx1 * ((*(*c34)) * (*(*tmp1))) - dt * tx1 * dx4;
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][4] = 0.;
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][0] = dt * tx2
                            * ((0.40000000000000002
                                       * (u[(i_darts_counter_temp7462) + 1]
                                           [j_darts_counter_temp7462][(*(*k))][1]
                                               * u[(i_darts_counter_temp7462) + 1]
                                                  [j_darts_counter_temp7462][(*(*k))][1]
                                           + u[(i_darts_counter_temp7462) + 1]
                                              [j_darts_counter_temp7462][(*(*k))][2]
                                               * u[(i_darts_counter_temp7462) + 1]
                                                  [j_darts_counter_temp7462][(*(*k))][2]
                                           + u[(i_darts_counter_temp7462) + 1]
                                              [j_darts_counter_temp7462][(*(*k))][3]
                                               * u[(i_darts_counter_temp7462) + 1]
                                                  [j_darts_counter_temp7462][(*(*k))][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(i_darts_counter_temp7462) + 1]
                                           [j_darts_counter_temp7462][(*(*k))][4]
                                           * (*(*tmp1))))
                                * (u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1
                            * (-((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                          [(*(*k))][1])
                                        * (u[(i_darts_counter_temp7462) + 1]
                                            [j_darts_counter_temp7462][(*(*k))][1])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                          [(*(*k))][2])
                                        * (u[(i_darts_counter_temp7462) + 1]
                                            [j_darts_counter_temp7462][(*(*k))][2])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                          [(*(*k))][3])
                                        * (u[(i_darts_counter_temp7462) + 1]
                                            [j_darts_counter_temp7462][(*(*k))][3])))
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                       [(*(*k))][4]);
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][1] = dt * tx2
                            * (1.3999999999999999
                                    * (u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                        [(*(*k))][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((3.
                                               * u[(i_darts_counter_temp7462) + 1]
                                                  [j_darts_counter_temp7462][(*(*k))][1]
                                               * u[(i_darts_counter_temp7462) + 1]
                                                  [j_darts_counter_temp7462][(*(*k))][1]
                                           + u[(i_darts_counter_temp7462) + 1]
                                              [j_darts_counter_temp7462][(*(*k))][2]
                                               * u[(i_darts_counter_temp7462) + 1]
                                                  [j_darts_counter_temp7462][(*(*k))][2]
                                           + u[(i_darts_counter_temp7462) + 1]
                                              [j_darts_counter_temp7462][(*(*k))][3]
                                               * u[(i_darts_counter_temp7462) + 1]
                                                  [j_darts_counter_temp7462][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * tx1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462][(*(*k))]
                               [1];
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][2] = dt * tx2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                    [(*(*k))][2]
                                    * u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                       [(*(*k))][1])
                                * (*(*tmp2)))
                        - dt * tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462][(*(*k))]
                               [2];
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][3] = dt * tx2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                    [(*(*k))][3]
                                    * u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                       [(*(*k))][1])
                                * (*(*tmp2)))
                        - dt * tx1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462][(*(*k))]
                               [3];
                    a[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][4] = dt * tx2
                            * (1.3999999999999999
                                * (u[(i_darts_counter_temp7462) + 1][j_darts_counter_temp7462]
                                    [(*(*k))][1]
                                    * (*(*tmp1))))
                        - dt * tx1 * (*(*c1345)) * (*(*tmp1)) - dt * tx1 * dx5;
                    (*(*tmp1)) = 1.
                        / u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1][(*(*k))][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][0] = -dt * ty1 * dy1;
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][1] = 0.;
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][2] = dt * ty2;
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][3] = 0.;
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][4] = 0.;
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][0] = dt * ty2
                            * (-(u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                  [(*(*k))][1]
                                   * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                      [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * ty1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                   [(*(*k))][1]);
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][1] = dt * ty2
                            * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1][(*(*k))]
                                [2]
                                * (*(*tmp1)))
                        - dt * ty1 * ((*(*c34)) * (*(*tmp1))) - dt * ty1 * dy2;
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][2] = dt * ty2
                        * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1][(*(*k))][1]
                            * (*(*tmp1)));
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][3] = 0.;
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][4] = 0.;
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][0] = dt * ty2
                            * (-(u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                  [(*(*k))][2]
                                   * (*(*tmp1)))
                                    * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                        [(*(*k))][2]
                                        * (*(*tmp1)))
                                + 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                         [(*(*k))][1]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462 + 1][(*(*k))][1]
                                           + u[(i_darts_counter_temp7462)]
                                              [j_darts_counter_temp7462 + 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462 + 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp7462)]
                                              [j_darts_counter_temp7462 + 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462 + 1][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * ty1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                   [(*(*k))][2]);
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][1] = dt * ty2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1][(*(*k))]
                                [1]
                                * (*(*tmp1))));
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][2] = dt * ty2
                            * ((2. - 0.40000000000000002)
                                * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * ty1 * dy3;
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][3] = dt * ty2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1][(*(*k))]
                                [3]
                                * (*(*tmp1))));
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][4]
                        = dt * ty2 * 0.40000000000000002;
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][0] = dt * ty2
                            * (-(u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                  [(*(*k))][2]
                                   * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                      [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * ty1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                   [(*(*k))][3]);
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][1] = 0.;
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][2] = dt * ty2
                        * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1][(*(*k))][3]
                            * (*(*tmp1)));
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][3] = dt * ty2
                            * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1][(*(*k))]
                                [2]
                                * (*(*tmp1)))
                        - dt * ty1 * ((*(*c34)) * (*(*tmp1))) - dt * ty1 * dy4;
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][4] = 0.;
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][0] = dt * ty2
                            * ((0.40000000000000002
                                       * (u[(i_darts_counter_temp7462)]
                                           [j_darts_counter_temp7462 + 1][(*(*k))][1]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462 + 1][(*(*k))][1]
                                           + u[(i_darts_counter_temp7462)]
                                              [j_darts_counter_temp7462 + 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462 + 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp7462)]
                                              [j_darts_counter_temp7462 + 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462 + 1][(*(*k))][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(i_darts_counter_temp7462)]
                                           [j_darts_counter_temp7462 + 1][(*(*k))][4]
                                           * (*(*tmp1))))
                                * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1
                            * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                          [(*(*k))][1])
                                        * (u[(i_darts_counter_temp7462)]
                                            [j_darts_counter_temp7462 + 1][(*(*k))][1])))
                                - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                          [(*(*k))][2])
                                        * (u[(i_darts_counter_temp7462)]
                                            [j_darts_counter_temp7462 + 1][(*(*k))][2])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                          [(*(*k))][3])
                                        * (u[(i_darts_counter_temp7462)]
                                            [j_darts_counter_temp7462 + 1][(*(*k))][3])))
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                       [(*(*k))][4]);
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][1] = dt * ty2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                    [(*(*k))][1]
                                    * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                       [(*(*k))][2])
                                * (*(*tmp2)))
                        - dt * ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1][(*(*k))]
                               [1];
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][2] = dt * ty2
                            * (1.3999999999999999
                                    * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                        [(*(*k))][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                         [(*(*k))][1]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462 + 1][(*(*k))][1]
                                           + 3.
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462 + 1][(*(*k))][2]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462 + 1][(*(*k))][2]
                                           + u[(i_darts_counter_temp7462)]
                                              [j_darts_counter_temp7462 + 1][(*(*k))][3]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462 + 1][(*(*k))][3])
                                        * (*(*tmp2))))
                        - dt * ty1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1][(*(*k))]
                               [2];
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][3] = dt * ty2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                    [(*(*k))][2]
                                    * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                       [(*(*k))][3])
                                * (*(*tmp2)))
                        - dt * ty1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1][(*(*k))]
                               [3];
                    b[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][4] = dt * ty2
                            * (1.3999999999999999
                                * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462 + 1]
                                    [(*(*k))][2]
                                    * (*(*tmp1))))
                        - dt * ty1 * (*(*c1345)) * (*(*tmp1)) - dt * ty1 * dy5;
                    (*(*tmp1)) = 1.
                        / u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k)) + 1][0];
                    (*(*tmp2)) = (*(*tmp1)) * (*(*tmp1));
                    (*(*tmp3)) = (*(*tmp1)) * (*(*tmp2));
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][0] = -dt * tz1 * dz1;
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][1] = 0.;
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][2] = 0.;
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][3] = dt * tz2;
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][0][4] = 0.;
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][0] = dt * tz2
                            * (-(u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                  [(*(*k)) + 1][1]
                                   * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                      [(*(*k)) + 1][3])
                                * (*(*tmp2)))
                        - dt * tz1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                   [(*(*k)) + 1][1]);
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][1] = dt * tz2
                            * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k)) + 1]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * (*(*c34)) * (*(*tmp1)) - dt * tz1 * dz2;
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][2] = 0.;
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][3] = dt * tz2
                        * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k)) + 1][1]
                            * (*(*tmp1)));
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][1][4] = 0.;
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][0] = dt * tz2
                            * (-(u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                  [(*(*k)) + 1][2]
                                   * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                      [(*(*k)) + 1][3])
                                * (*(*tmp2)))
                        - dt * tz1
                            * (-(*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                   [(*(*k)) + 1][2]);
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][1] = 0.;
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][2] = dt * tz2
                            * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k)) + 1]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * ((*(*c34)) * (*(*tmp1))) - dt * tz1 * dz3;
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][3] = dt * tz2
                        * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k)) + 1][2]
                            * (*(*tmp1)));
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][2][4] = 0.;
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][0] = dt * tz2
                            * (-(u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                  [(*(*k)) + 1][3]
                                   * (*(*tmp1)))
                                    * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                        [(*(*k)) + 1][3]
                                        * (*(*tmp1)))
                                + 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                         [(*(*k)) + 1][1]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462][(*(*k)) + 1][1]
                                           + u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                              [(*(*k)) + 1][2]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462][(*(*k)) + 1][2]
                                           + u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                              [(*(*k)) + 1][3]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462][(*(*k)) + 1][3])
                                        * (*(*tmp2))))
                        - dt * tz1
                            * (-(*(*r43)) * (*(*c34)) * (*(*tmp2))
                                * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                   [(*(*k)) + 1][3]);
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][1] = dt * tz2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k)) + 1]
                                [1]
                                * (*(*tmp1))));
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][2] = dt * tz2
                        * (-0.40000000000000002
                            * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k)) + 1]
                                [2]
                                * (*(*tmp1))));
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][3] = dt * tz2
                            * (2. - 0.40000000000000002)
                            * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k)) + 1]
                                [3]
                                * (*(*tmp1)))
                        - dt * tz1 * ((*(*r43)) * (*(*c34)) * (*(*tmp1))) - dt * tz1 * dz4;
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][3][4]
                        = dt * tz2 * 0.40000000000000002;
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][0] = dt * tz2
                            * ((0.40000000000000002
                                       * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                           [(*(*k)) + 1][1]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462][(*(*k)) + 1][1]
                                           + u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                              [(*(*k)) + 1][2]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462][(*(*k)) + 1][2]
                                           + u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                              [(*(*k)) + 1][3]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462][(*(*k)) + 1][3])
                                       * (*(*tmp2))
                                   - 1.3999999999999999
                                       * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                           [(*(*k)) + 1][4]
                                           * (*(*tmp1))))
                                * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                    [(*(*k)) + 1][3]
                                    * (*(*tmp1))))
                        - dt * tz1
                            * (-((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                          [(*(*k)) + 1][1])
                                        * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                            [(*(*k)) + 1][1])))
                                - ((*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                          [(*(*k)) + 1][2])
                                        * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                            [(*(*k)) + 1][2])))
                                - ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp3))
                                    * (((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                          [(*(*k)) + 1][3])
                                        * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                            [(*(*k)) + 1][3])))
                                - (*(*c1345)) * (*(*tmp2))
                                    * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                       [(*(*k)) + 1][4]);
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][1] = dt * tz2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                    [(*(*k)) + 1][1]
                                    * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                       [(*(*k)) + 1][3])
                                * (*(*tmp2)))
                        - dt * tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k)) + 1]
                               [1];
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][2] = dt * tz2
                            * (-0.40000000000000002
                                * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                    [(*(*k)) + 1][2]
                                    * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                       [(*(*k)) + 1][3])
                                * (*(*tmp2)))
                        - dt * tz1 * ((*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k)) + 1]
                               [2];
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][3] = dt * tz2
                            * (1.3999999999999999
                                    * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                        [(*(*k)) + 1][4]
                                        * (*(*tmp1)))
                                - 0.5 * 0.40000000000000002
                                    * ((u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                         [(*(*k)) + 1][1]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462][(*(*k)) + 1][1]
                                           + u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                              [(*(*k)) + 1][2]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462][(*(*k)) + 1][2]
                                           + 3.
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462][(*(*k)) + 1][3]
                                               * u[(i_darts_counter_temp7462)]
                                                  [j_darts_counter_temp7462][(*(*k)) + 1][3])
                                        * (*(*tmp2))))
                        - dt * tz1 * ((*(*r43)) * (*(*c34)) - (*(*c1345))) * (*(*tmp2))
                            * u[(i_darts_counter_temp7462)][j_darts_counter_temp7462][(*(*k)) + 1]
                               [3];
                    c[(i_darts_counter_temp7462)][j_darts_counter_temp7462][4][4] = dt * tz2
                            * (1.3999999999999999
                                * (u[(i_darts_counter_temp7462)][j_darts_counter_temp7462]
                                    [(*(*k)) + 1][3]
                                    * (*(*tmp1))))
                        - dt * tz1 * (*(*c1345)) * (*(*tmp1)) - dt * tz1 * dz5;
                }
                (*j) = j_darts_counter_temp7462;
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
TP7462::TP7462(int in_numThreads, int in_mainCodeletID, TP8* in_TPParent, int in_initIteration,
    int in_lastIteration, TP7462** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , c1345_darts7462(new double*[this->numThreads])
    , c34_darts7462(new double*[this->numThreads])
    , i_darts7462(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts7462(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts7462(new int*[this->numThreads])
    , r43_darts7462(new double*[this->numThreads])
    , tmp1_darts7462(new double*[this->numThreads])
    , tmp2_darts7462(new double*[this->numThreads])
    , tmp3_darts7462(new double*[this->numThreads])
    , initIteration7462(in_initIteration)
    , lastIteration7462(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , signalNextReady(new int[baseNumThreads])
    , checkInCodelets7463(new _checkInCodelets7463[this->numThreads])
{
    /*Initialize the loop parameters*/
    range7462 = abs(lastIteration7462 - initIteration7462) / 1;
    rangePerCodelet7462 = range7462 / numThreads;
    minIteration7462 = min<int>(lastIteration7462, initIteration7462);
    remainderRange7462 = range7462 % numThreads;
    /*Initialize inputs and vars.*/
    this->c1345_darts7462
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->c34_darts7462
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->k_darts7462 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->r43_darts7462
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp1_darts7462
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp2_darts7462
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp3_darts7462
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    _checkInCodelets7463* checkInCodelets7463Ptr = (this->checkInCodelets7463);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets7463);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
        this->signalNextReady[codeletCounter] = 0;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets7463Ptr) = _checkInCodelets7463(2, 1, this, codeletCounter);
#else
        (*checkInCodelets7463Ptr) = _checkInCodelets7463(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets7463Ptr).decDep();
        checkInCodelets7463Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP7462::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets7463[localID].setID(codeletID);
    this->checkInCodelets7463[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets7463[localID + this->baseNumThreads * i]
            = _checkInCodelets7463(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets7463[localID + this->baseNumThreads * i]
            = _checkInCodelets7463(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets7463[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets7463[localID + this->baseNumThreads * i].decDep();
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
TP7462::~TP7462()
{
    delete[] c1345_darts7462;
    delete[] c34_darts7462;
    delete[] k_darts7462;
    delete[] r43_darts7462;
    delete[] tmp1_darts7462;
    delete[] tmp2_darts7462;
    delete[] tmp3_darts7462;
    delete[] checkInCodelets7463;
}
/*TP9913: OMPParallelDirective*/
void TP9913::_barrierCodelets9913::fire(void)
{
    TP9913* myTP = static_cast<TP9913*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP9913::_checkInCodelets9915::fire(void)
{
    /*Init the vars for this region*/

    /*printing node 9915: DeclStmt*/

    /*printing node 9916: DeclStmt*/
    this->inputsTPParent->sum0_darts9913[this->getID()] = 0.;
    this->inputsTPParent->sum1_darts9913[this->getID()] = 0.;
    this->inputsTPParent->sum2_darts9913[this->getID()] = 0.;
    this->inputsTPParent->sum3_darts9913[this->getID()] = 0.;
    this->inputsTPParent->sum4_darts9913[this->getID()] = 0.;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 9915 nextRegion: 9922 */
    myTP->controlTPParent->checkInCodelets9922[this->getID()].decDep();
}
void TP9913::_checkInCodelets9922::fire(void)
{
    /*Select the thread executing OMPSingleDirective 9922*/
    if (!__sync_val_compare_and_swap(&(myTP->TP9922_alreadyLaunched), 0, 1)) {
        /*Init the vars for this region*/
        /*Initialize the vars of the inlined region*/
        this->inputsTPParent->sum_darts9922
            = (this->inputsTPParent->sum_darts9913) /*OMP_SHARED - VAR INLINED*/;

        /*printing node 9923: ForStmt*/
        {
            /*Loop's init*/
            (this->inputsTPParent->m_darts9922) = 0;
            int m_darts_counter_temp9922 = (this->inputsTPParent->m_darts9922);
            for (; m_darts_counter_temp9922 < 5; m_darts_counter_temp9922++) {
                (*(this->inputsTPParent->sum_darts9922))[m_darts_counter_temp9922] = 0.;
            }
            (this->inputsTPParent->m_darts9922) = m_darts_counter_temp9922;
        }
        /*Signaling next codelet from last stmt in the codelet*/
        /*Signaling omp region's barrier*/
        myTP->controlTPParent->barrierCodelets9922[0].decDep();
    } else {
        /*Signaling omp region's barrier*/
        myTP->barrierCodelets9922[0].decDep();
    }
}
void TP9913::_barrierCodelets9922::fire(void)
{
    TP9913* myTP = static_cast<TP9913*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets9931[codeletsCounter].decDep();
        }
    }
}
void TP9913::_checkInCodelets9931::fire(void)
{
    /*region 9931 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP9931;
    if (idx < myTP->TPsToUse9931) {
        if (!__sync_val_compare_and_swap(&(myTP->TP9931_alreadyLaunched[idx]), 0, 1)) {
            int range = abs((*(this->inputsTPParent->iend_darts9913))
                            - (*(this->inputsTPParent->ist_darts9913)))
                / 1;
            int rangePerCodelet = range / myTP->TPsToUse9931;
            int minIteration = min<int>((*(this->inputsTPParent->iend_darts9913)),
                (*(this->inputsTPParent->ist_darts9913)));
            int remainderRange = range % myTP->TPsToUse9931;
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
            if ((*(this->inputsTPParent->ist_darts9913))
                < (*(this->inputsTPParent->iend_darts9913))) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse9931 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse9931 - 1) {
                lastIteration = (*(this->inputsTPParent->iend_darts9913));
            }
#if USEINVOKE == 1
            invoke<TP9931>(myTP, myTP->codeletsPerTP9931 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration, &(*(this->inputsTPParent->iend_darts9913)),
                &(*(this->inputsTPParent->ist_darts9913)),
                &(*(this->inputsTPParent->jend_darts9913)),
                &(*(this->inputsTPParent->jst_darts9913)),
                &(*(this->inputsTPParent->nz0_darts9913)), &(myTP->TP9931Ptr[idx]));
#else
            place<TP9931>(idx, myTP, myTP->codeletsPerTP9931 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(*(this->inputsTPParent->iend_darts9913)),
                &(*(this->inputsTPParent->ist_darts9913)),
                &(*(this->inputsTPParent->jend_darts9913)),
                &(*(this->inputsTPParent->jst_darts9913)),
                &(*(this->inputsTPParent->nz0_darts9913)), &(myTP->TP9931Ptr[idx]));
#endif
        } else {
            if (myTP->TP9931Ptr[idx] != nullptr) {
                myTP->TP9931Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    } else {
        /*Signaling next codelet region: 9931 nextRegion: 10028 */
        myTP->controlTPParent->checkInCodelets10028[this->getID()].decDep();
    }
}
void TP9913::_checkInCodelets10028::fire(void)
{

    /*printing node 10028: CompoundAssignOperator*/
    TP10026mutex.lock();
    (*(this->inputsTPParent->sum_darts9913))[0]
        += (this->inputsTPParent->sum0_darts9913[this->getID()]);

    /*printing node 10030: CompoundAssignOperator*/
    (*(this->inputsTPParent->sum_darts9913))[1]
        += (this->inputsTPParent->sum1_darts9913[this->getID()]);

    /*printing node 10032: CompoundAssignOperator*/
    (*(this->inputsTPParent->sum_darts9913))[2]
        += (this->inputsTPParent->sum2_darts9913[this->getID()]);

    /*printing node 10034: CompoundAssignOperator*/
    (*(this->inputsTPParent->sum_darts9913))[3]
        += (this->inputsTPParent->sum3_darts9913[this->getID()]);

    /*printing node 10036: CompoundAssignOperator*/
    (*(this->inputsTPParent->sum_darts9913))[4]
        += (this->inputsTPParent->sum4_darts9913[this->getID()]);
    TP10026mutex.unlock();
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 10028 nextRegion: 10038 */
    myTP->controlTPParent->barrierCodelets10038[0].decDep();
}
void TP9913::_barrierCodelets10038::fire(void)
{
    TP9913* myTP = static_cast<TP9913*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets10039[codeletsCounter].decDep();
        }
    }
}
void TP9913::_checkInCodelets10039::fire(void)
{
    /*Select the thread executing OMPSingleDirective 10039*/
    if (!__sync_val_compare_and_swap(&(myTP->TP10039_alreadyLaunched), 0, 1)) {
        /*Init the vars for this region*/
        /*Initialize the vars of the inlined region*/
        this->inputsTPParent->nx0_darts10039
            = (this->inputsTPParent->nx0_darts9913) /*OMP_SHARED - VAR INLINED*/;
        this->inputsTPParent->ny0_darts10039
            = (this->inputsTPParent->ny0_darts9913) /*OMP_SHARED - VAR INLINED*/;
        this->inputsTPParent->nz0_darts10039
            = (this->inputsTPParent->nz0_darts9913) /*OMP_SHARED - VAR INLINED*/;
        this->inputsTPParent->sum_darts10039
            = (this->inputsTPParent->sum_darts9913) /*OMP_SHARED - VAR INLINED*/;

        /*printing node 10040: ForStmt*/
        {
            /*Loop's init*/
            (this->inputsTPParent->m_darts10039) = 0;
            int m_darts_counter_temp10039 = (this->inputsTPParent->m_darts10039);
            for (; m_darts_counter_temp10039 < 5; m_darts_counter_temp10039++) {
                (*(this->inputsTPParent->sum_darts10039))[m_darts_counter_temp10039]
                    = sqrt((*(this->inputsTPParent->sum_darts10039))[m_darts_counter_temp10039]
                        / (((*(this->inputsTPParent->nx0_darts10039)) - 2)
                            * ((*(this->inputsTPParent->ny0_darts10039)) - 2)
                            * ((*(this->inputsTPParent->nz0_darts10039)) - 2)));
            }
            (this->inputsTPParent->m_darts10039) = m_darts_counter_temp10039;
        }
        /*Signaling next codelet from last stmt in the codelet*/
        /*Signaling omp region's barrier*/
        myTP->controlTPParent->barrierCodelets10039[0].decDep();
    } else {
        /*Signaling omp region's barrier*/
        myTP->barrierCodelets10039[0].decDep();
    }
}
void TP9913::_barrierCodelets10039::fire(void)
{
    TP9913* myTP = static_cast<TP9913*>(myTP_);
    myTP->TPParent->barrierCodelets9913[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets9913[0]));
}
TP9913::TP9913(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet, int* in_iend,
    int* in_ist, int* in_jend, int* in_jst, int* in_nx0, int* in_ny0, int* in_nz0, double** in_sum)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , iend_darts9913(in_iend) /*OMP_SHARED - INPUT*/
    , ist_darts9913(in_ist) /*OMP_SHARED - INPUT*/
    , jend_darts9913(in_jend) /*OMP_SHARED - INPUT*/
    , jst_darts9913(in_jst) /*OMP_SHARED - INPUT*/
    , nx0_darts9913(in_nx0) /*OMP_SHARED - INPUT*/
    , ny0_darts9913(in_ny0) /*OMP_SHARED - INPUT*/
    , nz0_darts9913(in_nz0) /*OMP_SHARED - INPUT*/
    , sum_darts9913(in_sum) /*OMP_SHARED - INPUT*/
    , i_darts9913(new int[this->numThreads]) /*VARIABLE*/
    , j_darts9913(new int[this->numThreads]) /*VARIABLE*/
    , k_darts9913(new int[this->numThreads]) /*VARIABLE*/
    , m_darts9913(new int[this->numThreads]) /*VARIABLE*/
    , sum0_darts9913(new double[this->numThreads]) /*VARIABLE*/
    , sum1_darts9913(new double[this->numThreads]) /*VARIABLE*/
    , sum2_darts9913(new double[this->numThreads]) /*VARIABLE*/
    , sum3_darts9913(new double[this->numThreads]) /*VARIABLE*/
    , sum4_darts9913(new double[this->numThreads]) /*VARIABLE*/
    , TP9922_alreadyLaunched(0)
    , TP9931Ptr(new TP9931*[NUMTPS9931])
    , TP9931_alreadyLaunched(new size_t[NUMTPS9931])
    , numTPsSet9931(0)
    , numTPsReady9931(0)
    , TPsToUse9931(NUMTPS9931)
    , codeletsPerTP9931(this->numThreads / NUMTPS9931)
    , totalCodelets9931(this->TPsToUse9931 * this->codeletsPerTP9931)
    , TP10039_alreadyLaunched(0)
    , barrierCodelets9913(new _barrierCodelets9913[1])
    , checkInCodelets9915(new _checkInCodelets9915[this->numThreads])
    , checkInCodelets9922(new _checkInCodelets9922[this->numThreads])
    , barrierCodelets9922(new _barrierCodelets9922[1])
    , checkInCodelets9931(new _checkInCodelets9931[this->numThreads])
    , checkInCodelets10028(new _checkInCodelets10028[this->numThreads])
    , barrierCodelets10038(new _barrierCodelets10038[1])
    , checkInCodelets10039(new _checkInCodelets10039[this->numThreads])
    , barrierCodelets10039(new _barrierCodelets10039[1])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets9913[0] = _barrierCodelets9913(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets10039[0] = _barrierCodelets10039(this->numThreads, this->numThreads, this, 0);
    barrierCodelets10038[0] = _barrierCodelets10038(this->numThreads, this->numThreads, this, 0);
    barrierCodelets9922[0] = _barrierCodelets9922(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets10039* checkInCodelets10039Ptr = (this->checkInCodelets10039);
    _checkInCodelets10028* checkInCodelets10028Ptr = (this->checkInCodelets10028);
    _checkInCodelets9931* checkInCodelets9931Ptr = (this->checkInCodelets9931);
    for (int i = 0; i < NUMTPS9931; i++) {
        TP9931Ptr[i] = nullptr;
        TP9931_alreadyLaunched[i] = 0;
    }
    _checkInCodelets9922* checkInCodelets9922Ptr = (this->checkInCodelets9922);
    _checkInCodelets9915* checkInCodelets9915Ptr = (this->checkInCodelets9915);
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets10039Ptr) = _checkInCodelets10039(1, 1, this, codeletCounter);
        checkInCodelets10039Ptr++;
        (*checkInCodelets10028Ptr) = _checkInCodelets10028(1, 1, this, codeletCounter);
        checkInCodelets10028Ptr++;
        (*checkInCodelets9931Ptr) = _checkInCodelets9931(1, 1, this, codeletCounter);
        checkInCodelets9931Ptr++;
        (*checkInCodelets9922Ptr) = _checkInCodelets9922(1, 1, this, codeletCounter);
        checkInCodelets9922Ptr++;
        (*checkInCodelets9915Ptr) = _checkInCodelets9915(1, 1, this, codeletCounter);
        (*checkInCodelets9915Ptr).decDep();
        checkInCodelets9915Ptr++;
    }
}
TP9913::~TP9913()
{
    delete[] i_darts9913;
    delete[] j_darts9913;
    delete[] k_darts9913;
    delete[] m_darts9913;
    delete[] sum0_darts9913;
    delete[] sum1_darts9913;
    delete[] sum2_darts9913;
    delete[] sum3_darts9913;
    delete[] sum4_darts9913;
    delete[] barrierCodelets9913;
    delete[] barrierCodelets10039;
    delete[] checkInCodelets10039;
    delete[] barrierCodelets10038;
    delete[] checkInCodelets10028;
    delete[] checkInCodelets9931;
    delete[] barrierCodelets9922;
    delete[] checkInCodelets9922;
    delete[] checkInCodelets9915;
}
/*TP9931: OMPForDirective*/
bool TP9931::requestNewRangeIterations9931(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet9931 * codeletID;
        int tempEndRange = rangePerCodelet9931 * (codeletID + 1);
        if (remainderRange9931 != 0) {
            if (codeletID < (uint32_t)remainderRange9931) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange9931;
                tempEndRange += remainderRange9931;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration9931;
        tempEndRange = tempEndRange * 1 + minIteration9931;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration9931 < lastIteration9931) {
            (this->inputsTPParent->i_darts9931[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts9931[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration9931;
        }
    }
    return isThereNewIteration;
}
void TP9931::_checkInCodelets9932::fire(void)
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
    this->inputsTPParent->sum0_darts9931[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->sum0_darts9913[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->sum1_darts9931[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->sum1_darts9913[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->sum2_darts9931[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->sum2_darts9913[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->sum3_darts9931[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->sum3_darts9913[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->sum4_darts9931[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->sum4_darts9913[this->getID()]);

    /*printing node 9932: ForStmt*/
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
    int* i = &(this->inputsTPParent->i_darts9931[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts9931[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* jend = (this->inputsTPParent->jend_darts9931);
    (void)jend /*OMP_SHARED*/;
    int* jst = (this->inputsTPParent->jst_darts9931);
    (void)jst /*OMP_SHARED*/;
    int* k = &(this->inputsTPParent->k_darts9931[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* nz0 = (this->inputsTPParent->nz0_darts9931);
    (void)nz0 /*OMP_SHARED*/;
    double** sum0 = &(this->inputsTPParent->sum0_darts9931[this->getLocalID()]);
    (void)sum0 /*OMP_SHARED_PRIVATE*/;
    double** sum1 = &(this->inputsTPParent->sum1_darts9931[this->getLocalID()]);
    (void)sum1 /*OMP_SHARED_PRIVATE*/;
    double** sum2 = &(this->inputsTPParent->sum2_darts9931[this->getLocalID()]);
    (void)sum2 /*OMP_SHARED_PRIVATE*/;
    double** sum3 = &(this->inputsTPParent->sum3_darts9931[this->getLocalID()]);
    (void)sum3 /*OMP_SHARED_PRIVATE*/;
    double** sum4 = &(this->inputsTPParent->sum4_darts9931[this->getLocalID()]);
    (void)sum4 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations9931(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Find and signal the next codelet*/
        myTP->controlTPParent->TPParent->checkInCodelets10028[this->getID()].decDep();
        return;
    }
    for (int i_darts_counter_temp9931 = (*i); i_darts_counter_temp9931 <= endRange
         && i_darts_counter_temp9931 <= this->inputsTPParent->lastIteration9931;
         i_darts_counter_temp9931++) {
        {
            {
                /*Loop's init*/
                (*j) = (*(jst));
                int j_darts_counter_temp9931 = (*j);
                for (; j_darts_counter_temp9931 <= (*(jend)); j_darts_counter_temp9931++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp9931 = (*k);
                        for (; k_darts_counter_temp9931 <= (*(nz0)) - 2;
                             k_darts_counter_temp9931++) {
                            (*(*sum0)) = (*(*sum0))
                                + rsd[(i_darts_counter_temp9931)][j_darts_counter_temp9931]
                                     [k_darts_counter_temp9931][0]
                                    * rsd[(i_darts_counter_temp9931)][j_darts_counter_temp9931]
                                         [k_darts_counter_temp9931][0];
                            (*(*sum1)) = (*(*sum1))
                                + rsd[(i_darts_counter_temp9931)][j_darts_counter_temp9931]
                                     [k_darts_counter_temp9931][1]
                                    * rsd[(i_darts_counter_temp9931)][j_darts_counter_temp9931]
                                         [k_darts_counter_temp9931][1];
                            (*(*sum2)) = (*(*sum2))
                                + rsd[(i_darts_counter_temp9931)][j_darts_counter_temp9931]
                                     [k_darts_counter_temp9931][2]
                                    * rsd[(i_darts_counter_temp9931)][j_darts_counter_temp9931]
                                         [k_darts_counter_temp9931][2];
                            (*(*sum3)) = (*(*sum3))
                                + rsd[(i_darts_counter_temp9931)][j_darts_counter_temp9931]
                                     [k_darts_counter_temp9931][3]
                                    * rsd[(i_darts_counter_temp9931)][j_darts_counter_temp9931]
                                         [k_darts_counter_temp9931][3];
                            (*(*sum4)) = (*(*sum4))
                                + rsd[(i_darts_counter_temp9931)][j_darts_counter_temp9931]
                                     [k_darts_counter_temp9931][4]
                                    * rsd[(i_darts_counter_temp9931)][j_darts_counter_temp9931]
                                         [k_darts_counter_temp9931][4];
                        }
                        (*k) = k_darts_counter_temp9931;
                    }
                }
                (*j) = j_darts_counter_temp9931;
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
    myTP->controlTPParent->TPParent->checkInCodelets10028[this->getID()].decDep();
}
TP9931::TP9931(int in_numThreads, int in_mainCodeletID, TP9913* in_TPParent, int in_initIteration,
    int in_lastIteration, int* in_iend, int* in_ist, int* in_jend, int* in_jst, int* in_nz0,
    TP9931** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts9931(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iend_darts9931(in_iend) /*OMP_SHARED - INPUT*/
    , ist_darts9931(in_ist) /*OMP_SHARED - INPUT*/
    , j_darts9931(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jend_darts9931(in_jend) /*OMP_SHARED - INPUT*/
    , jst_darts9931(in_jst) /*OMP_SHARED - INPUT*/
    , k_darts9931(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , nz0_darts9931(in_nz0) /*OMP_SHARED - INPUT*/
    , sum0_darts9931(new double*[this->numThreads])
    , sum1_darts9931(new double*[this->numThreads])
    , sum2_darts9931(new double*[this->numThreads])
    , sum3_darts9931(new double*[this->numThreads])
    , sum4_darts9931(new double*[this->numThreads])
    , initIteration9931(in_initIteration)
    , lastIteration9931(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , signalNextReady(new int[baseNumThreads])
    , checkInCodelets9932(new _checkInCodelets9932[this->numThreads])
{
    /*Initialize the loop parameters*/
    range9931 = abs(lastIteration9931 - initIteration9931) / 1;
    rangePerCodelet9931 = range9931 / numThreads;
    minIteration9931 = min<int>(lastIteration9931, initIteration9931);
    remainderRange9931 = range9931 % numThreads;
    /*Initialize inputs and vars.*/
    this->sum0_darts9931
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->sum1_darts9931
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->sum2_darts9931
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->sum3_darts9931
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->sum4_darts9931
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    _checkInCodelets9932* checkInCodelets9932Ptr = (this->checkInCodelets9932);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets9932);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
        this->signalNextReady[codeletCounter] = 0;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets9932Ptr) = _checkInCodelets9932(2, 1, this, codeletCounter);
#else
        (*checkInCodelets9932Ptr) = _checkInCodelets9932(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets9932Ptr).decDep();
        checkInCodelets9932Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP9931::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets9932[localID].setID(codeletID);
    this->checkInCodelets9932[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets9932[localID + this->baseNumThreads * i]
            = _checkInCodelets9932(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets9932[localID + this->baseNumThreads * i]
            = _checkInCodelets9932(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets9932[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets9932[localID + this->baseNumThreads * i].decDep();
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
TP9931::~TP9931()
{
    delete[] sum0_darts9931;
    delete[] sum1_darts9931;
    delete[] sum2_darts9931;
    delete[] sum3_darts9931;
    delete[] sum4_darts9931;
    delete[] checkInCodelets9932;
}
/*TP10830: OMPParallelDirective*/
void TP10830::_barrierCodelets10830::fire(void)
{
    TP10830* myTP = static_cast<TP10830*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP10830::_checkInCodelets10845::fire(void)
{
    /*region 10845 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP10845;
    if (idx < myTP->TPsToUse10845) {
        if (!__sync_val_compare_and_swap(&(myTP->TP10845_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(nx - 1 - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse10845;
            int minIteration = min<int>(nx - 1, 0);
            int remainderRange = range % myTP->TPsToUse10845;
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
            if (idx == myTP->TPsToUse10845 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse10845 - 1) {
                lastIteration = nx - 1;
            }
#if USEINVOKE == 1
            invoke<TP10845>(myTP, myTP->codeletsPerTP10845 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP10845Ptr[idx]));
#else
            place<TP10845>(idx, myTP, myTP->codeletsPerTP10845 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP10845Ptr[idx]));
#endif
        } else {
            if (myTP->TP10845Ptr[idx] != nullptr) {
                myTP->TP10845Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10830::_barrierCodelets10845::fire(void)
{
    TP10830* myTP = static_cast<TP10830*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets10902[codeletsCounter].decDep();
        }
    }
}
void TP10830::_checkInCodelets10902::fire(void)
{

    /*printing node 10902: BinaryOperator*/
    (this->inputsTPParent->L1_darts10830[this->getID()]) = 0;

    /*printing node 10903: BinaryOperator*/
    (this->inputsTPParent->L2_darts10830[this->getID()]) = nx - 1;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 10902 nextRegion: 10905 */
    myTP->controlTPParent->checkInCodelets10905[this->getID()].decDep();
}
void TP10830::_checkInCodelets10905::fire(void)
{
    /*region 10905 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP10905;
    if (idx < myTP->TPsToUse10905) {
        if (!__sync_val_compare_and_swap(&(myTP->TP10905_alreadyLaunched[idx]), 0, 1)) {
            int range = abs((this->inputsTPParent->L2_darts10830[this->getID()])
                            - (this->inputsTPParent->L1_darts10830[this->getID()]))
                / 1;
            int rangePerCodelet = range / myTP->TPsToUse10905;
            int minIteration = min<int>((this->inputsTPParent->L2_darts10830[this->getID()]),
                (this->inputsTPParent->L1_darts10830[this->getID()]));
            int remainderRange = range % myTP->TPsToUse10905;
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
            if ((this->inputsTPParent->L1_darts10830[this->getID()])
                < (this->inputsTPParent->L2_darts10830[this->getID()])) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse10905 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse10905 - 1) {
                lastIteration = (this->inputsTPParent->L2_darts10830[this->getID()]);
            }
#if USEINVOKE == 1
            invoke<TP10905>(myTP, myTP->codeletsPerTP10905 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP10905Ptr[idx]));
#else
            place<TP10905>(idx, myTP, myTP->codeletsPerTP10905 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP10905Ptr[idx]));
#endif
        } else {
            if (myTP->TP10905Ptr[idx] != nullptr) {
                myTP->TP10905Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10830::_barrierCodelets10905::fire(void)
{
    TP10830* myTP = static_cast<TP10830*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets11054[codeletsCounter].decDep();
        }
    }
}
void TP10830::_checkInCodelets11054::fire(void)
{
    /*region 11054 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP11054;
    if (idx < myTP->TPsToUse11054) {
        if (!__sync_val_compare_and_swap(&(myTP->TP11054_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(jend - jst) / 1;
            int rangePerCodelet = range / myTP->TPsToUse11054;
            int minIteration = min<int>(jend, jst);
            int remainderRange = range % myTP->TPsToUse11054;
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
            if (idx == myTP->TPsToUse11054 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse11054 - 1) {
                lastIteration = jend;
            }
#if USEINVOKE == 1
            invoke<TP11054>(myTP, myTP->codeletsPerTP11054 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11054Ptr[idx]));
#else
            place<TP11054>(idx, myTP, myTP->codeletsPerTP11054 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11054Ptr[idx]));
#endif
        } else {
            if (myTP->TP11054Ptr[idx] != nullptr) {
                myTP->TP11054Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10830::_barrierCodelets11054::fire(void)
{
    TP10830* myTP = static_cast<TP10830*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets11691[codeletsCounter].decDep();
        }
    }
}
void TP10830::_checkInCodelets11691::fire(void)
{

    /*printing node 11691: BinaryOperator*/
    (this->inputsTPParent->L1_darts10830[this->getID()]) = 0;

    /*printing node 11692: BinaryOperator*/
    (this->inputsTPParent->L2_darts10830[this->getID()]) = ny - 1;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 11691 nextRegion: 11694 */
    myTP->controlTPParent->checkInCodelets11694[this->getID()].decDep();
}
void TP10830::_checkInCodelets11694::fire(void)
{
    /*region 11694 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP11694;
    if (idx < myTP->TPsToUse11694) {
        if (!__sync_val_compare_and_swap(&(myTP->TP11694_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse11694;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse11694;
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
            if (idx == myTP->TPsToUse11694 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse11694 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP11694>(myTP, myTP->codeletsPerTP11694 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11694Ptr[idx]));
#else
            place<TP11694>(idx, myTP, myTP->codeletsPerTP11694 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11694Ptr[idx]));
#endif
        } else {
            if (myTP->TP11694Ptr[idx] != nullptr) {
                myTP->TP11694Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10830::_barrierCodelets11694::fire(void)
{
    TP10830* myTP = static_cast<TP10830*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets11843[codeletsCounter].decDep();
        }
    }
}
void TP10830::_checkInCodelets11843::fire(void)
{
    /*region 11843 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP11843;
    if (idx < myTP->TPsToUse11843) {
        if (!__sync_val_compare_and_swap(&(myTP->TP11843_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse11843;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse11843;
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
            if (idx == myTP->TPsToUse11843 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse11843 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP11843>(myTP, myTP->codeletsPerTP11843 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11843Ptr[idx]));
#else
            place<TP11843>(idx, myTP, myTP->codeletsPerTP11843 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP11843Ptr[idx]));
#endif
        } else {
            if (myTP->TP11843Ptr[idx] != nullptr) {
                myTP->TP11843Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10830::_barrierCodelets11843::fire(void)
{
    TP10830* myTP = static_cast<TP10830*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets12480[codeletsCounter].decDep();
        }
    }
}
void TP10830::_checkInCodelets12480::fire(void)
{
    /*region 12480 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP12480;
    if (idx < myTP->TPsToUse12480) {
        if (!__sync_val_compare_and_swap(&(myTP->TP12480_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse12480;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse12480;
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
            if (idx == myTP->TPsToUse12480 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse12480 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP12480>(myTP, myTP->codeletsPerTP12480 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP12480Ptr[idx]));
#else
            place<TP12480>(idx, myTP, myTP->codeletsPerTP12480 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP12480Ptr[idx]));
#endif
        } else {
            if (myTP->TP12480Ptr[idx] != nullptr) {
                myTP->TP12480Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP10830::_barrierCodelets12480::fire(void)
{
    TP10830* myTP = static_cast<TP10830*>(myTP_);
    myTP->TPParent->barrierCodelets10830[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets10830[0]));
}
TP10830::TP10830(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts10830(new int[this->numThreads]) /*VARIABLE*/
    , L2_darts10830(new int[this->numThreads]) /*VARIABLE*/
    , i_darts10830(new int[this->numThreads]) /*VARIABLE*/
    , iend1_darts10830(new int[this->numThreads]) /*VARIABLE*/
    , ist1_darts10830(new int[this->numThreads]) /*VARIABLE*/
    , j_darts10830(new int[this->numThreads]) /*VARIABLE*/
    , jend1_darts10830(new int[this->numThreads]) /*VARIABLE*/
    , jst1_darts10830(new int[this->numThreads]) /*VARIABLE*/
    , k_darts10830(new int[this->numThreads]) /*VARIABLE*/
    , m_darts10830(new int[this->numThreads]) /*VARIABLE*/
    , q_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , tmp_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u21_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u21i_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u21im1_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u21j_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u21jm1_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u21k_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u21km1_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u31_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u31i_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u31im1_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u31j_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u31jm1_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u31k_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u31km1_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u41_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u41i_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u41im1_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u41j_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u41jm1_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u41k_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u41km1_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u51i_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u51im1_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u51j_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u51jm1_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u51k_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , u51km1_darts10830(new double[this->numThreads]) /*VARIABLE*/
    , TP10845Ptr(new TP10845*[NUMTPS10845])
    , TP10845_alreadyLaunched(new size_t[NUMTPS10845])
    , numTPsSet10845(0)
    , numTPsReady10845(0)
    , TPsToUse10845(NUMTPS10845)
    , codeletsPerTP10845(this->numThreads / NUMTPS10845)
    , totalCodelets10845(this->TPsToUse10845 * this->codeletsPerTP10845)
    , TP10905Ptr(new TP10905*[NUMTPS10905])
    , TP10905_alreadyLaunched(new size_t[NUMTPS10905])
    , numTPsSet10905(0)
    , numTPsReady10905(0)
    , TPsToUse10905(NUMTPS10905)
    , codeletsPerTP10905(this->numThreads / NUMTPS10905)
    , totalCodelets10905(this->TPsToUse10905 * this->codeletsPerTP10905)
    , TP11054Ptr(new TP11054*[NUMTPS11054])
    , TP11054_alreadyLaunched(new size_t[NUMTPS11054])
    , numTPsSet11054(0)
    , numTPsReady11054(0)
    , TPsToUse11054(NUMTPS11054)
    , codeletsPerTP11054(this->numThreads / NUMTPS11054)
    , totalCodelets11054(this->TPsToUse11054 * this->codeletsPerTP11054)
    , TP11694Ptr(new TP11694*[NUMTPS11694])
    , TP11694_alreadyLaunched(new size_t[NUMTPS11694])
    , numTPsSet11694(0)
    , numTPsReady11694(0)
    , TPsToUse11694(NUMTPS11694)
    , codeletsPerTP11694(this->numThreads / NUMTPS11694)
    , totalCodelets11694(this->TPsToUse11694 * this->codeletsPerTP11694)
    , TP11843Ptr(new TP11843*[NUMTPS11843])
    , TP11843_alreadyLaunched(new size_t[NUMTPS11843])
    , numTPsSet11843(0)
    , numTPsReady11843(0)
    , TPsToUse11843(NUMTPS11843)
    , codeletsPerTP11843(this->numThreads / NUMTPS11843)
    , totalCodelets11843(this->TPsToUse11843 * this->codeletsPerTP11843)
    , TP12480Ptr(new TP12480*[NUMTPS12480])
    , TP12480_alreadyLaunched(new size_t[NUMTPS12480])
    , numTPsSet12480(0)
    , numTPsReady12480(0)
    , TPsToUse12480(NUMTPS12480)
    , codeletsPerTP12480(this->numThreads / NUMTPS12480)
    , totalCodelets12480(this->TPsToUse12480 * this->codeletsPerTP12480)
    , barrierCodelets10830(new _barrierCodelets10830[1])
    , checkInCodelets10845(new _checkInCodelets10845[this->numThreads])
    , barrierCodelets10845(new _barrierCodelets10845[1])
    , checkInCodelets10902(new _checkInCodelets10902[this->numThreads])
    , checkInCodelets10905(new _checkInCodelets10905[this->numThreads])
    , barrierCodelets10905(new _barrierCodelets10905[1])
    , checkInCodelets11054(new _checkInCodelets11054[this->numThreads])
    , barrierCodelets11054(new _barrierCodelets11054[1])
    , checkInCodelets11691(new _checkInCodelets11691[this->numThreads])
    , checkInCodelets11694(new _checkInCodelets11694[this->numThreads])
    , barrierCodelets11694(new _barrierCodelets11694[1])
    , checkInCodelets11843(new _checkInCodelets11843[this->numThreads])
    , barrierCodelets11843(new _barrierCodelets11843[1])
    , checkInCodelets12480(new _checkInCodelets12480[this->numThreads])
    , barrierCodelets12480(new _barrierCodelets12480[1])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets10830[0] = _barrierCodelets10830(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets12480[0] = _barrierCodelets12480(NUMTPS12480, NUMTPS12480, this, 0);
    barrierCodelets11843[0] = _barrierCodelets11843(NUMTPS11843, NUMTPS11843, this, 0);
    barrierCodelets11694[0] = _barrierCodelets11694(NUMTPS11694, NUMTPS11694, this, 0);
    barrierCodelets11054[0] = _barrierCodelets11054(NUMTPS11054, NUMTPS11054, this, 0);
    barrierCodelets10905[0] = _barrierCodelets10905(NUMTPS10905, NUMTPS10905, this, 0);
    barrierCodelets10845[0] = _barrierCodelets10845(NUMTPS10845, NUMTPS10845, this, 0);
    _checkInCodelets12480* checkInCodelets12480Ptr = (this->checkInCodelets12480);
    for (int i = 0; i < NUMTPS12480; i++) {
        TP12480Ptr[i] = nullptr;
        TP12480_alreadyLaunched[i] = 0;
    }
    _checkInCodelets11843* checkInCodelets11843Ptr = (this->checkInCodelets11843);
    for (int i = 0; i < NUMTPS11843; i++) {
        TP11843Ptr[i] = nullptr;
        TP11843_alreadyLaunched[i] = 0;
    }
    _checkInCodelets11694* checkInCodelets11694Ptr = (this->checkInCodelets11694);
    for (int i = 0; i < NUMTPS11694; i++) {
        TP11694Ptr[i] = nullptr;
        TP11694_alreadyLaunched[i] = 0;
    }
    _checkInCodelets11691* checkInCodelets11691Ptr = (this->checkInCodelets11691);
    _checkInCodelets11054* checkInCodelets11054Ptr = (this->checkInCodelets11054);
    for (int i = 0; i < NUMTPS11054; i++) {
        TP11054Ptr[i] = nullptr;
        TP11054_alreadyLaunched[i] = 0;
    }
    _checkInCodelets10905* checkInCodelets10905Ptr = (this->checkInCodelets10905);
    for (int i = 0; i < NUMTPS10905; i++) {
        TP10905Ptr[i] = nullptr;
        TP10905_alreadyLaunched[i] = 0;
    }
    _checkInCodelets10902* checkInCodelets10902Ptr = (this->checkInCodelets10902);
    _checkInCodelets10845* checkInCodelets10845Ptr = (this->checkInCodelets10845);
    for (int i = 0; i < NUMTPS10845; i++) {
        TP10845Ptr[i] = nullptr;
        TP10845_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets12480Ptr) = _checkInCodelets12480(1, 1, this, codeletCounter);
        checkInCodelets12480Ptr++;
        (*checkInCodelets11843Ptr) = _checkInCodelets11843(1, 1, this, codeletCounter);
        checkInCodelets11843Ptr++;
        (*checkInCodelets11694Ptr) = _checkInCodelets11694(1, 1, this, codeletCounter);
        checkInCodelets11694Ptr++;
        (*checkInCodelets11691Ptr) = _checkInCodelets11691(1, 1, this, codeletCounter);
        checkInCodelets11691Ptr++;
        (*checkInCodelets11054Ptr) = _checkInCodelets11054(1, 1, this, codeletCounter);
        checkInCodelets11054Ptr++;
        (*checkInCodelets10905Ptr) = _checkInCodelets10905(1, 1, this, codeletCounter);
        checkInCodelets10905Ptr++;
        (*checkInCodelets10902Ptr) = _checkInCodelets10902(1, 1, this, codeletCounter);
        checkInCodelets10902Ptr++;
        (*checkInCodelets10845Ptr) = _checkInCodelets10845(1, 1, this, codeletCounter);
        (*checkInCodelets10845Ptr).decDep();
        checkInCodelets10845Ptr++;
    }
}
TP10830::~TP10830()
{
    delete[] L1_darts10830;
    delete[] L2_darts10830;
    delete[] i_darts10830;
    delete[] iend1_darts10830;
    delete[] ist1_darts10830;
    delete[] j_darts10830;
    delete[] jend1_darts10830;
    delete[] jst1_darts10830;
    delete[] k_darts10830;
    delete[] m_darts10830;
    delete[] q_darts10830;
    delete[] tmp_darts10830;
    delete[] u21_darts10830;
    delete[] u21i_darts10830;
    delete[] u21im1_darts10830;
    delete[] u21j_darts10830;
    delete[] u21jm1_darts10830;
    delete[] u21k_darts10830;
    delete[] u21km1_darts10830;
    delete[] u31_darts10830;
    delete[] u31i_darts10830;
    delete[] u31im1_darts10830;
    delete[] u31j_darts10830;
    delete[] u31jm1_darts10830;
    delete[] u31k_darts10830;
    delete[] u31km1_darts10830;
    delete[] u41_darts10830;
    delete[] u41i_darts10830;
    delete[] u41im1_darts10830;
    delete[] u41j_darts10830;
    delete[] u41jm1_darts10830;
    delete[] u41k_darts10830;
    delete[] u41km1_darts10830;
    delete[] u51i_darts10830;
    delete[] u51im1_darts10830;
    delete[] u51j_darts10830;
    delete[] u51jm1_darts10830;
    delete[] u51k_darts10830;
    delete[] u51km1_darts10830;
    delete[] barrierCodelets10830;
    delete[] barrierCodelets12480;
    delete[] checkInCodelets12480;
    delete[] barrierCodelets11843;
    delete[] checkInCodelets11843;
    delete[] barrierCodelets11694;
    delete[] checkInCodelets11694;
    delete[] checkInCodelets11691;
    delete[] barrierCodelets11054;
    delete[] checkInCodelets11054;
    delete[] barrierCodelets10905;
    delete[] checkInCodelets10905;
    delete[] checkInCodelets10902;
    delete[] barrierCodelets10845;
    delete[] checkInCodelets10845;
}
/*TP10845: OMPForDirective*/
void TP10845::_barrierCodelets10845::fire(void)
{
    TP10845* myTP = static_cast<TP10845*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets10845[0].decDep();
}
bool TP10845::requestNewRangeIterations10845(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet10845 * codeletID;
        int tempEndRange = rangePerCodelet10845 * (codeletID + 1);
        if (remainderRange10845 != 0) {
            if (codeletID < (uint32_t)remainderRange10845) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange10845;
                tempEndRange += remainderRange10845;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration10845;
        tempEndRange = tempEndRange * 1 + minIteration10845;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration10845 < lastIteration10845) {
            (this->inputsTPParent->i_darts10845[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts10845[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration10845;
        }
    }
    return isThereNewIteration;
}
void TP10845::_checkInCodelets10846::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 10846: ForStmt*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    int* i = &(this->inputsTPParent->i_darts10845[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts10845[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts10845[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts10845[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations10845(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets10845[0].decDep();
        return;
    }
    for (int i_darts_counter_temp10845 = (*i); i_darts_counter_temp10845 <= endRange
         && i_darts_counter_temp10845 <= this->inputsTPParent->lastIteration10845;
         i_darts_counter_temp10845++) {
        {
            {
                /*Loop's init*/
                (*j) = 0;
                int j_darts_counter_temp10845 = (*j);
                for (; j_darts_counter_temp10845 <= ny - 1; j_darts_counter_temp10845++) {
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp10845 = (*k);
                        for (; k_darts_counter_temp10845 <= nz - 1; k_darts_counter_temp10845++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp10845 = (*m);
                                for (; m_darts_counter_temp10845 < 5; m_darts_counter_temp10845++) {
                                    rsd[(i_darts_counter_temp10845)][j_darts_counter_temp10845]
                                       [k_darts_counter_temp10845][m_darts_counter_temp10845]
                                        = -frct[(i_darts_counter_temp10845)]
                                               [j_darts_counter_temp10845]
                                               [k_darts_counter_temp10845]
                                               [m_darts_counter_temp10845];
                                }
                                (*m) = m_darts_counter_temp10845;
                            }
                        }
                        (*k) = k_darts_counter_temp10845;
                    }
                }
                (*j) = j_darts_counter_temp10845;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets10845[0].decDep();
}
TP10845::TP10845(int in_numThreads, int in_mainCodeletID, TP10830* in_TPParent,
    int in_initIteration, int in_lastIteration, TP10845** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts10845(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts10845(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts10845(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts10845(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration10845(in_initIteration)
    , lastIteration10845(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets10845(new _barrierCodelets10845[1])
    , checkInCodelets10846(new _checkInCodelets10846[this->numThreads])
{
    /*Initialize the loop parameters*/
    range10845 = abs(lastIteration10845 - initIteration10845) / 1;
    rangePerCodelet10845 = range10845 / numThreads;
    minIteration10845 = min<int>(lastIteration10845, initIteration10845);
    remainderRange10845 = range10845 % numThreads;
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets10845[0] = _barrierCodelets10845(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets10846* checkInCodelets10846Ptr = (this->checkInCodelets10846);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets10846);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets10846Ptr) = _checkInCodelets10846(2, 1, this, codeletCounter);
#else
        (*checkInCodelets10846Ptr) = _checkInCodelets10846(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets10846Ptr).decDep();
        checkInCodelets10846Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP10845::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets10846[localID].setID(codeletID);
    this->checkInCodelets10846[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets10846[localID + this->baseNumThreads * i]
            = _checkInCodelets10846(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets10846[localID + this->baseNumThreads * i]
            = _checkInCodelets10846(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets10846[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets10846[localID + this->baseNumThreads * i].decDep();
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
TP10845::~TP10845()
{
    delete[] barrierCodelets10845;
    delete[] checkInCodelets10846;
}
/*TP10905: OMPForDirective*/
void TP10905::_barrierCodelets10905::fire(void)
{
    TP10905* myTP = static_cast<TP10905*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets10905[0].decDep();
}
bool TP10905::requestNewRangeIterations10905(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet10905 * codeletID;
        int tempEndRange = rangePerCodelet10905 * (codeletID + 1);
        if (remainderRange10905 != 0) {
            if (codeletID < (uint32_t)remainderRange10905) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange10905;
                tempEndRange += remainderRange10905;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration10905;
        tempEndRange = tempEndRange * 1 + minIteration10905;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration10905 < lastIteration10905) {
            (this->inputsTPParent->i_darts10905[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts10905[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration10905;
        }
    }
    return isThereNewIteration;
}
void TP10905::_checkInCodelets10906::fire(void)
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
    this->inputsTPParent->L1_darts10905[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L1_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L2_darts10905[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->q_darts10905[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21_darts10905[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21_darts10830[this->getID()]);

    /*printing node 10906: ForStmt*/
    /*var: L1*/
    /*var: L2*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: q*/
    /*var: u21*/
    int* i = &(this->inputsTPParent->i_darts10905[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts10905[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts10905[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts10905[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** u21 = &(this->inputsTPParent->u21_darts10905[this->getLocalID()]);
    (void)u21 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations10905(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets10905[0].decDep();
        return;
    }
    for (int i_darts_counter_temp10905 = (*i); i_darts_counter_temp10905 <= endRange
         && i_darts_counter_temp10905 <= this->inputsTPParent->lastIteration10905;
         i_darts_counter_temp10905++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp10905 = (*j);
                for (; j_darts_counter_temp10905 <= jend; j_darts_counter_temp10905++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp10905 = (*k);
                        for (; k_darts_counter_temp10905 <= nz - 2; k_darts_counter_temp10905++) {
                            flux[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                [k_darts_counter_temp10905][0]
                                = u[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                   [k_darts_counter_temp10905][1];
                            (*(*u21)) = u[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                         [k_darts_counter_temp10905][1]
                                / u[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                   [k_darts_counter_temp10905][0];
                            (*(*q)) = 0.5
                                * (u[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                    [k_darts_counter_temp10905][1]
                                        * u[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                           [k_darts_counter_temp10905][1]
                                    + u[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                       [k_darts_counter_temp10905][2]
                                        * u[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                           [k_darts_counter_temp10905][2]
                                    + u[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                       [k_darts_counter_temp10905][3]
                                        * u[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                           [k_darts_counter_temp10905][3])
                                / u[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                   [k_darts_counter_temp10905][0];
                            flux[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                [k_darts_counter_temp10905][1]
                                = u[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                   [k_darts_counter_temp10905][1]
                                    * (*(*u21))
                                + 0.40000000000000002
                                    * (u[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                        [k_darts_counter_temp10905][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                [k_darts_counter_temp10905][2]
                                = u[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                   [k_darts_counter_temp10905][2]
                                * (*(*u21));
                            flux[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                [k_darts_counter_temp10905][3]
                                = u[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                   [k_darts_counter_temp10905][3]
                                * (*(*u21));
                            flux[(i_darts_counter_temp10905)][j_darts_counter_temp10905]
                                [k_darts_counter_temp10905][4]
                                = (1.3999999999999999
                                          * u[(i_darts_counter_temp10905)]
                                             [j_darts_counter_temp10905][k_darts_counter_temp10905]
                                             [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u21));
                        }
                        (*k) = k_darts_counter_temp10905;
                    }
                }
                (*j) = j_darts_counter_temp10905;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets10905[0].decDep();
}
TP10905::TP10905(int in_numThreads, int in_mainCodeletID, TP10830* in_TPParent,
    int in_initIteration, int in_lastIteration, TP10905** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts10905(new int*[this->numThreads])
    , L2_darts10905(new int*[this->numThreads])
    , i_darts10905(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts10905(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts10905(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts10905(new double*[this->numThreads])
    , u21_darts10905(new double*[this->numThreads])
    , initIteration10905(in_initIteration)
    , lastIteration10905(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets10905(new _barrierCodelets10905[1])
    , checkInCodelets10906(new _checkInCodelets10906[this->numThreads])
{
    /*Initialize the loop parameters*/
    range10905 = abs(lastIteration10905 - initIteration10905) / 1;
    rangePerCodelet10905 = range10905 / numThreads;
    minIteration10905 = min<int>(lastIteration10905, initIteration10905);
    remainderRange10905 = range10905 % numThreads;
    /*Initialize inputs and vars.*/
    this->L1_darts10905 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->L2_darts10905 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->q_darts10905
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21_darts10905
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets10905[0] = _barrierCodelets10905(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets10906* checkInCodelets10906Ptr = (this->checkInCodelets10906);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets10906);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets10906Ptr) = _checkInCodelets10906(2, 1, this, codeletCounter);
#else
        (*checkInCodelets10906Ptr) = _checkInCodelets10906(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets10906Ptr).decDep();
        checkInCodelets10906Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP10905::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets10906[localID].setID(codeletID);
    this->checkInCodelets10906[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets10906[localID + this->baseNumThreads * i]
            = _checkInCodelets10906(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets10906[localID + this->baseNumThreads * i]
            = _checkInCodelets10906(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets10906[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets10906[localID + this->baseNumThreads * i].decDep();
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
TP10905::~TP10905()
{
    delete[] L1_darts10905;
    delete[] L2_darts10905;
    delete[] q_darts10905;
    delete[] u21_darts10905;
    delete[] barrierCodelets10905;
    delete[] checkInCodelets10906;
}
/*TP11054: OMPForDirective*/
void TP11054::_barrierCodelets11054::fire(void)
{
    TP11054* myTP = static_cast<TP11054*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets11054[0].decDep();
}
bool TP11054::requestNewRangeIterations11054(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet11054 * codeletID;
        int tempEndRange = rangePerCodelet11054 * (codeletID + 1);
        if (remainderRange11054 != 0) {
            if (codeletID < (uint32_t)remainderRange11054) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange11054;
                tempEndRange += remainderRange11054;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration11054;
        tempEndRange = tempEndRange * 1 + minIteration11054;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration11054 < lastIteration11054) {
            (this->inputsTPParent->j_darts11054[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->j_darts11054[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration11054;
        }
    }
    return isThereNewIteration;
}
void TP11054::_checkInCodelets11055::fire(void)
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
    this->inputsTPParent->L2_darts11054[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->iend1_darts11054[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iend1_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->ist1_darts11054[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->ist1_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts11054[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21i_darts11054[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21i_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21im1_darts11054[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21im1_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31i_darts11054[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31i_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31im1_darts11054[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31im1_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41i_darts11054[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41i_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41im1_darts11054[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41im1_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51i_darts11054[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51i_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51im1_darts11054[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51im1_darts10830[this->getID()]);

    /*printing node 11055: ForStmt*/
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
    int** L2 = &(this->inputsTPParent->L2_darts11054[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts11054[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iend1 = &(this->inputsTPParent->iend1_darts11054[this->getLocalID()]);
    (void)iend1 /*OMP_SHARED_PRIVATE*/;
    int** ist1 = &(this->inputsTPParent->ist1_darts11054[this->getLocalID()]);
    (void)ist1 /*OMP_SHARED_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts11054[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts11054[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts11054[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts11054[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21i = &(this->inputsTPParent->u21i_darts11054[this->getLocalID()]);
    (void)u21i /*OMP_SHARED_PRIVATE*/;
    double** u21im1 = &(this->inputsTPParent->u21im1_darts11054[this->getLocalID()]);
    (void)u21im1 /*OMP_SHARED_PRIVATE*/;
    double** u31i = &(this->inputsTPParent->u31i_darts11054[this->getLocalID()]);
    (void)u31i /*OMP_SHARED_PRIVATE*/;
    double** u31im1 = &(this->inputsTPParent->u31im1_darts11054[this->getLocalID()]);
    (void)u31im1 /*OMP_SHARED_PRIVATE*/;
    double** u41i = &(this->inputsTPParent->u41i_darts11054[this->getLocalID()]);
    (void)u41i /*OMP_SHARED_PRIVATE*/;
    double** u41im1 = &(this->inputsTPParent->u41im1_darts11054[this->getLocalID()]);
    (void)u41im1 /*OMP_SHARED_PRIVATE*/;
    double** u51i = &(this->inputsTPParent->u51i_darts11054[this->getLocalID()]);
    (void)u51i /*OMP_SHARED_PRIVATE*/;
    double** u51im1 = &(this->inputsTPParent->u51im1_darts11054[this->getLocalID()]);
    (void)u51im1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations11054(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets11054[0].decDep();
        return;
    }
    for (int j_darts_counter_temp11054 = (*j); j_darts_counter_temp11054 <= endRange
         && j_darts_counter_temp11054 <= this->inputsTPParent->lastIteration11054;
         j_darts_counter_temp11054++) {
        {
            {
                /*Loop's init*/
                (*k) = 1;
                int k_darts_counter_temp11054 = (*k);
                for (; k_darts_counter_temp11054 <= nz - 2; k_darts_counter_temp11054++) {
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp11054 = (*i);
                        for (; i_darts_counter_temp11054 <= iend; i_darts_counter_temp11054++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp11054 = (*m);
                                for (; m_darts_counter_temp11054 < 5; m_darts_counter_temp11054++) {
                                    rsd[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                                       [k_darts_counter_temp11054][m_darts_counter_temp11054]
                                        = rsd[i_darts_counter_temp11054]
                                             [(j_darts_counter_temp11054)]
                                             [k_darts_counter_temp11054][m_darts_counter_temp11054]
                                        - tx2
                                            * (flux[i_darts_counter_temp11054 + 1]
                                                   [(j_darts_counter_temp11054)]
                                                   [k_darts_counter_temp11054]
                                                   [m_darts_counter_temp11054]
                                                - flux[i_darts_counter_temp11054 - 1]
                                                      [(j_darts_counter_temp11054)]
                                                      [k_darts_counter_temp11054]
                                                      [m_darts_counter_temp11054]);
                                }
                                (*m) = m_darts_counter_temp11054;
                            }
                        }
                        (*i) = i_darts_counter_temp11054;
                    }
                    (*(*L2)) = nx - 1;
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp11054 = (*i);
                        for (; i_darts_counter_temp11054 <= (*(*L2)); i_darts_counter_temp11054++) {
                            (*(*tmp)) = 1.
                                / u[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                                   [k_darts_counter_temp11054][0];
                            (*(*u21i)) = (*(*tmp))
                                * u[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                                   [k_darts_counter_temp11054][1];
                            (*(*u31i)) = (*(*tmp))
                                * u[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                                   [k_darts_counter_temp11054][2];
                            (*(*u41i)) = (*(*tmp))
                                * u[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                                   [k_darts_counter_temp11054][3];
                            (*(*u51i)) = (*(*tmp))
                                * u[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                                   [k_darts_counter_temp11054][4];
                            (*(*tmp)) = 1.
                                / u[i_darts_counter_temp11054 - 1][(j_darts_counter_temp11054)]
                                   [k_darts_counter_temp11054][0];
                            (*(*u21im1)) = (*(*tmp))
                                * u[i_darts_counter_temp11054 - 1][(j_darts_counter_temp11054)]
                                   [k_darts_counter_temp11054][1];
                            (*(*u31im1)) = (*(*tmp))
                                * u[i_darts_counter_temp11054 - 1][(j_darts_counter_temp11054)]
                                   [k_darts_counter_temp11054][2];
                            (*(*u41im1)) = (*(*tmp))
                                * u[i_darts_counter_temp11054 - 1][(j_darts_counter_temp11054)]
                                   [k_darts_counter_temp11054][3];
                            (*(*u51im1)) = (*(*tmp))
                                * u[i_darts_counter_temp11054 - 1][(j_darts_counter_temp11054)]
                                   [k_darts_counter_temp11054][4];
                            flux[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                                [k_darts_counter_temp11054][1]
                                = (4. / 3.) * tx3 * ((*(*u21i)) - (*(*u21im1)));
                            flux[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                                [k_darts_counter_temp11054][2]
                                = tx3 * ((*(*u31i)) - (*(*u31im1)));
                            flux[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                                [k_darts_counter_temp11054][3]
                                = tx3 * ((*(*u41i)) - (*(*u41im1)));
                            flux[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                                [k_darts_counter_temp11054][4]
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
                        (*i) = i_darts_counter_temp11054;
                    }
                    {
                        /*Loop's init*/
                        (*i) = ist;
                        int i_darts_counter_temp11054 = (*i);
                        for (; i_darts_counter_temp11054 <= iend; i_darts_counter_temp11054++) {
                            rsd[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                               [k_darts_counter_temp11054][0]
                                = rsd[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                                     [k_darts_counter_temp11054][0]
                                + dx1 * tx1
                                    * (u[i_darts_counter_temp11054 - 1][(j_darts_counter_temp11054)]
                                        [k_darts_counter_temp11054][0]
                                        - 2.
                                            * u[i_darts_counter_temp11054]
                                               [(j_darts_counter_temp11054)]
                                               [k_darts_counter_temp11054][0]
                                        + u[i_darts_counter_temp11054 + 1]
                                           [(j_darts_counter_temp11054)][k_darts_counter_temp11054]
                                           [0]);
                            rsd[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                               [k_darts_counter_temp11054][1]
                                = rsd[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                                     [k_darts_counter_temp11054][1]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp11054 + 1][(
                                           j_darts_counter_temp11054)][k_darts_counter_temp11054][1]
                                        - flux[i_darts_counter_temp11054]
                                              [(j_darts_counter_temp11054)]
                                              [k_darts_counter_temp11054][1])
                                + dx2 * tx1
                                    * (u[i_darts_counter_temp11054 - 1][(j_darts_counter_temp11054)]
                                        [k_darts_counter_temp11054][1]
                                        - 2.
                                            * u[i_darts_counter_temp11054]
                                               [(j_darts_counter_temp11054)]
                                               [k_darts_counter_temp11054][1]
                                        + u[i_darts_counter_temp11054 + 1]
                                           [(j_darts_counter_temp11054)][k_darts_counter_temp11054]
                                           [1]);
                            rsd[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                               [k_darts_counter_temp11054][2]
                                = rsd[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                                     [k_darts_counter_temp11054][2]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp11054 + 1][(
                                           j_darts_counter_temp11054)][k_darts_counter_temp11054][2]
                                        - flux[i_darts_counter_temp11054]
                                              [(j_darts_counter_temp11054)]
                                              [k_darts_counter_temp11054][2])
                                + dx3 * tx1
                                    * (u[i_darts_counter_temp11054 - 1][(j_darts_counter_temp11054)]
                                        [k_darts_counter_temp11054][2]
                                        - 2.
                                            * u[i_darts_counter_temp11054]
                                               [(j_darts_counter_temp11054)]
                                               [k_darts_counter_temp11054][2]
                                        + u[i_darts_counter_temp11054 + 1]
                                           [(j_darts_counter_temp11054)][k_darts_counter_temp11054]
                                           [2]);
                            rsd[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                               [k_darts_counter_temp11054][3]
                                = rsd[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                                     [k_darts_counter_temp11054][3]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp11054 + 1][(
                                           j_darts_counter_temp11054)][k_darts_counter_temp11054][3]
                                        - flux[i_darts_counter_temp11054]
                                              [(j_darts_counter_temp11054)]
                                              [k_darts_counter_temp11054][3])
                                + dx4 * tx1
                                    * (u[i_darts_counter_temp11054 - 1][(j_darts_counter_temp11054)]
                                        [k_darts_counter_temp11054][3]
                                        - 2.
                                            * u[i_darts_counter_temp11054]
                                               [(j_darts_counter_temp11054)]
                                               [k_darts_counter_temp11054][3]
                                        + u[i_darts_counter_temp11054 + 1]
                                           [(j_darts_counter_temp11054)][k_darts_counter_temp11054]
                                           [3]);
                            rsd[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                               [k_darts_counter_temp11054][4]
                                = rsd[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                                     [k_darts_counter_temp11054][4]
                                + tx3 * 0.10000000000000001 * 1.
                                    * (flux[i_darts_counter_temp11054 + 1][(
                                           j_darts_counter_temp11054)][k_darts_counter_temp11054][4]
                                        - flux[i_darts_counter_temp11054]
                                              [(j_darts_counter_temp11054)]
                                              [k_darts_counter_temp11054][4])
                                + dx5 * tx1
                                    * (u[i_darts_counter_temp11054 - 1][(j_darts_counter_temp11054)]
                                        [k_darts_counter_temp11054][4]
                                        - 2.
                                            * u[i_darts_counter_temp11054]
                                               [(j_darts_counter_temp11054)]
                                               [k_darts_counter_temp11054][4]
                                        + u[i_darts_counter_temp11054 + 1]
                                           [(j_darts_counter_temp11054)][k_darts_counter_temp11054]
                                           [4]);
                        }
                        (*i) = i_darts_counter_temp11054;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp11054 = (*m);
                        for (; m_darts_counter_temp11054 < 5; m_darts_counter_temp11054++) {
                            rsd[1][(j_darts_counter_temp11054)][k_darts_counter_temp11054]
                               [m_darts_counter_temp11054]
                                = rsd[1][(j_darts_counter_temp11054)][k_darts_counter_temp11054]
                                     [m_darts_counter_temp11054]
                                - dssp
                                    * (+5.
                                            * u[1][(j_darts_counter_temp11054)]
                                               [k_darts_counter_temp11054]
                                               [m_darts_counter_temp11054]
                                        - 4.
                                            * u[2][(j_darts_counter_temp11054)]
                                               [k_darts_counter_temp11054]
                                               [m_darts_counter_temp11054]
                                        + u[3][(j_darts_counter_temp11054)]
                                           [k_darts_counter_temp11054][m_darts_counter_temp11054]);
                            rsd[2][(j_darts_counter_temp11054)][k_darts_counter_temp11054]
                               [m_darts_counter_temp11054]
                                = rsd[2][(j_darts_counter_temp11054)][k_darts_counter_temp11054]
                                     [m_darts_counter_temp11054]
                                - dssp
                                    * (-4.
                                            * u[1][(j_darts_counter_temp11054)]
                                               [k_darts_counter_temp11054]
                                               [m_darts_counter_temp11054]
                                        + 6.
                                            * u[2][(j_darts_counter_temp11054)]
                                               [k_darts_counter_temp11054]
                                               [m_darts_counter_temp11054]
                                        - 4.
                                            * u[3][(j_darts_counter_temp11054)]
                                               [k_darts_counter_temp11054]
                                               [m_darts_counter_temp11054]
                                        + u[4][(j_darts_counter_temp11054)]
                                           [k_darts_counter_temp11054][m_darts_counter_temp11054]);
                        }
                        (*m) = m_darts_counter_temp11054;
                    }
                    (*(*ist1)) = 3;
                    (*(*iend1)) = nx - 4;
                    {
                        /*Loop's init*/
                        (*i) = (*(*ist1));
                        int i_darts_counter_temp11054 = (*i);
                        for (; i_darts_counter_temp11054 <= (*(*iend1));
                             i_darts_counter_temp11054++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp11054 = (*m);
                                for (; m_darts_counter_temp11054 < 5; m_darts_counter_temp11054++) {
                                    rsd[i_darts_counter_temp11054][(j_darts_counter_temp11054)]
                                       [k_darts_counter_temp11054][m_darts_counter_temp11054]
                                        = rsd[i_darts_counter_temp11054]
                                             [(j_darts_counter_temp11054)]
                                             [k_darts_counter_temp11054][m_darts_counter_temp11054]
                                        - dssp
                                            * (u[i_darts_counter_temp11054 - 2]
                                                [(j_darts_counter_temp11054)]
                                                [k_darts_counter_temp11054]
                                                [m_darts_counter_temp11054]
                                                - 4.
                                                    * u[i_darts_counter_temp11054 - 1]
                                                       [(j_darts_counter_temp11054)]
                                                       [k_darts_counter_temp11054]
                                                       [m_darts_counter_temp11054]
                                                + 6.
                                                    * u[i_darts_counter_temp11054]
                                                       [(j_darts_counter_temp11054)]
                                                       [k_darts_counter_temp11054]
                                                       [m_darts_counter_temp11054]
                                                - 4.
                                                    * u[i_darts_counter_temp11054 + 1]
                                                       [(j_darts_counter_temp11054)]
                                                       [k_darts_counter_temp11054]
                                                       [m_darts_counter_temp11054]
                                                + u[i_darts_counter_temp11054 + 2]
                                                   [(j_darts_counter_temp11054)]
                                                   [k_darts_counter_temp11054]
                                                   [m_darts_counter_temp11054]);
                                }
                                (*m) = m_darts_counter_temp11054;
                            }
                        }
                        (*i) = i_darts_counter_temp11054;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp11054 = (*m);
                        for (; m_darts_counter_temp11054 < 5; m_darts_counter_temp11054++) {
                            rsd[nx - 3][(j_darts_counter_temp11054)][k_darts_counter_temp11054]
                               [m_darts_counter_temp11054]
                                = rsd[nx - 3][(j_darts_counter_temp11054)]
                                     [k_darts_counter_temp11054][m_darts_counter_temp11054]
                                - dssp
                                    * (u[nx - 5][(j_darts_counter_temp11054)]
                                        [k_darts_counter_temp11054][m_darts_counter_temp11054]
                                        - 4.
                                            * u[nx - 4][(j_darts_counter_temp11054)]
                                               [k_darts_counter_temp11054]
                                               [m_darts_counter_temp11054]
                                        + 6.
                                            * u[nx - 3][(j_darts_counter_temp11054)]
                                               [k_darts_counter_temp11054]
                                               [m_darts_counter_temp11054]
                                        - 4.
                                            * u[nx - 2][(j_darts_counter_temp11054)]
                                               [k_darts_counter_temp11054]
                                               [m_darts_counter_temp11054]);
                            rsd[nx - 2][(j_darts_counter_temp11054)][k_darts_counter_temp11054]
                               [m_darts_counter_temp11054]
                                = rsd[nx - 2][(j_darts_counter_temp11054)]
                                     [k_darts_counter_temp11054][m_darts_counter_temp11054]
                                - dssp
                                    * (u[nx - 4][(j_darts_counter_temp11054)]
                                        [k_darts_counter_temp11054][m_darts_counter_temp11054]
                                        - 4.
                                            * u[nx - 3][(j_darts_counter_temp11054)]
                                               [k_darts_counter_temp11054]
                                               [m_darts_counter_temp11054]
                                        + 5.
                                            * u[nx - 2][(j_darts_counter_temp11054)]
                                               [k_darts_counter_temp11054]
                                               [m_darts_counter_temp11054]);
                        }
                        (*m) = m_darts_counter_temp11054;
                    }
                }
                (*k) = k_darts_counter_temp11054;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets11054[0].decDep();
}
TP11054::TP11054(int in_numThreads, int in_mainCodeletID, TP10830* in_TPParent,
    int in_initIteration, int in_lastIteration, TP11054** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L2_darts11054(new int*[this->numThreads])
    , i_darts11054(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iend1_darts11054(new int*[this->numThreads])
    , ist1_darts11054(new int*[this->numThreads])
    , j_darts11054(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts11054(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts11054(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts11054(new double*[this->numThreads])
    , u21i_darts11054(new double*[this->numThreads])
    , u21im1_darts11054(new double*[this->numThreads])
    , u31i_darts11054(new double*[this->numThreads])
    , u31im1_darts11054(new double*[this->numThreads])
    , u41i_darts11054(new double*[this->numThreads])
    , u41im1_darts11054(new double*[this->numThreads])
    , u51i_darts11054(new double*[this->numThreads])
    , u51im1_darts11054(new double*[this->numThreads])
    , initIteration11054(in_initIteration)
    , lastIteration11054(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets11054(new _barrierCodelets11054[1])
    , checkInCodelets11055(new _checkInCodelets11055[this->numThreads])
{
    /*Initialize the loop parameters*/
    range11054 = abs(lastIteration11054 - initIteration11054) / 1;
    rangePerCodelet11054 = range11054 / numThreads;
    minIteration11054 = min<int>(lastIteration11054, initIteration11054);
    remainderRange11054 = range11054 % numThreads;
    /*Initialize inputs and vars.*/
    this->L2_darts11054 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->iend1_darts11054 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->ist1_darts11054 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts11054
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21i_darts11054
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21im1_darts11054
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31i_darts11054
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31im1_darts11054
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41i_darts11054
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41im1_darts11054
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51i_darts11054
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51im1_darts11054
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets11054[0] = _barrierCodelets11054(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets11055* checkInCodelets11055Ptr = (this->checkInCodelets11055);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets11055);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets11055Ptr) = _checkInCodelets11055(2, 1, this, codeletCounter);
#else
        (*checkInCodelets11055Ptr) = _checkInCodelets11055(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets11055Ptr).decDep();
        checkInCodelets11055Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP11054::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets11055[localID].setID(codeletID);
    this->checkInCodelets11055[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets11055[localID + this->baseNumThreads * i]
            = _checkInCodelets11055(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets11055[localID + this->baseNumThreads * i]
            = _checkInCodelets11055(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets11055[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets11055[localID + this->baseNumThreads * i].decDep();
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
TP11054::~TP11054()
{
    delete[] L2_darts11054;
    delete[] iend1_darts11054;
    delete[] ist1_darts11054;
    delete[] tmp_darts11054;
    delete[] u21i_darts11054;
    delete[] u21im1_darts11054;
    delete[] u31i_darts11054;
    delete[] u31im1_darts11054;
    delete[] u41i_darts11054;
    delete[] u41im1_darts11054;
    delete[] u51i_darts11054;
    delete[] u51im1_darts11054;
    delete[] barrierCodelets11054;
    delete[] checkInCodelets11055;
}
/*TP11694: OMPForDirective*/
void TP11694::_barrierCodelets11694::fire(void)
{
    TP11694* myTP = static_cast<TP11694*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets11694[0].decDep();
}
bool TP11694::requestNewRangeIterations11694(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet11694 * codeletID;
        int tempEndRange = rangePerCodelet11694 * (codeletID + 1);
        if (remainderRange11694 != 0) {
            if (codeletID < (uint32_t)remainderRange11694) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange11694;
                tempEndRange += remainderRange11694;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration11694;
        tempEndRange = tempEndRange * 1 + minIteration11694;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration11694 < lastIteration11694) {
            (this->inputsTPParent->i_darts11694[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts11694[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration11694;
        }
    }
    return isThereNewIteration;
}
void TP11694::_checkInCodelets11695::fire(void)
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
    this->inputsTPParent->L1_darts11694[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L1_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->L2_darts11694[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->q_darts11694[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31_darts11694[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31_darts10830[this->getID()]);

    /*printing node 11695: ForStmt*/
    /*var: L1*/
    /*var: L2*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: q*/
    /*var: u31*/
    int** L1 = &(this->inputsTPParent->L1_darts11694[this->getLocalID()]);
    (void)L1 /*OMP_SHARED_PRIVATE*/;
    int** L2 = &(this->inputsTPParent->L2_darts11694[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts11694[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts11694[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts11694[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts11694[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** u31 = &(this->inputsTPParent->u31_darts11694[this->getLocalID()]);
    (void)u31 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations11694(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets11694[0].decDep();
        return;
    }
    for (int i_darts_counter_temp11694 = (*i); i_darts_counter_temp11694 <= endRange
         && i_darts_counter_temp11694 <= this->inputsTPParent->lastIteration11694;
         i_darts_counter_temp11694++) {
        {
            {
                /*Loop's init*/
                (*j) = (*(*L1));
                int j_darts_counter_temp11694 = (*j);
                for (; j_darts_counter_temp11694 <= (*(*L2)); j_darts_counter_temp11694++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp11694 = (*k);
                        for (; k_darts_counter_temp11694 <= nz - 2; k_darts_counter_temp11694++) {
                            flux[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                [k_darts_counter_temp11694][0]
                                = u[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                   [k_darts_counter_temp11694][2];
                            (*(*u31)) = u[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                         [k_darts_counter_temp11694][2]
                                / u[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                   [k_darts_counter_temp11694][0];
                            (*(*q)) = 0.5
                                * (u[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                    [k_darts_counter_temp11694][1]
                                        * u[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                           [k_darts_counter_temp11694][1]
                                    + u[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                       [k_darts_counter_temp11694][2]
                                        * u[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                           [k_darts_counter_temp11694][2]
                                    + u[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                       [k_darts_counter_temp11694][3]
                                        * u[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                           [k_darts_counter_temp11694][3])
                                / u[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                   [k_darts_counter_temp11694][0];
                            flux[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                [k_darts_counter_temp11694][1]
                                = u[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                   [k_darts_counter_temp11694][1]
                                * (*(*u31));
                            flux[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                [k_darts_counter_temp11694][2]
                                = u[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                   [k_darts_counter_temp11694][2]
                                    * (*(*u31))
                                + 0.40000000000000002
                                    * (u[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                        [k_darts_counter_temp11694][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                [k_darts_counter_temp11694][3]
                                = u[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                   [k_darts_counter_temp11694][3]
                                * (*(*u31));
                            flux[(i_darts_counter_temp11694)][j_darts_counter_temp11694]
                                [k_darts_counter_temp11694][4]
                                = (1.3999999999999999
                                          * u[(i_darts_counter_temp11694)]
                                             [j_darts_counter_temp11694][k_darts_counter_temp11694]
                                             [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u31));
                        }
                        (*k) = k_darts_counter_temp11694;
                    }
                }
                (*j) = j_darts_counter_temp11694;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets11694[0].decDep();
}
TP11694::TP11694(int in_numThreads, int in_mainCodeletID, TP10830* in_TPParent,
    int in_initIteration, int in_lastIteration, TP11694** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L1_darts11694(new int*[this->numThreads])
    , L2_darts11694(new int*[this->numThreads])
    , i_darts11694(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts11694(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts11694(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts11694(new double*[this->numThreads])
    , u31_darts11694(new double*[this->numThreads])
    , initIteration11694(in_initIteration)
    , lastIteration11694(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets11694(new _barrierCodelets11694[1])
    , checkInCodelets11695(new _checkInCodelets11695[this->numThreads])
{
    /*Initialize the loop parameters*/
    range11694 = abs(lastIteration11694 - initIteration11694) / 1;
    rangePerCodelet11694 = range11694 / numThreads;
    minIteration11694 = min<int>(lastIteration11694, initIteration11694);
    remainderRange11694 = range11694 % numThreads;
    /*Initialize inputs and vars.*/
    this->L1_darts11694 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->L2_darts11694 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->q_darts11694
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31_darts11694
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets11694[0] = _barrierCodelets11694(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets11695* checkInCodelets11695Ptr = (this->checkInCodelets11695);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets11695);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets11695Ptr) = _checkInCodelets11695(2, 1, this, codeletCounter);
#else
        (*checkInCodelets11695Ptr) = _checkInCodelets11695(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets11695Ptr).decDep();
        checkInCodelets11695Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP11694::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets11695[localID].setID(codeletID);
    this->checkInCodelets11695[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets11695[localID + this->baseNumThreads * i]
            = _checkInCodelets11695(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets11695[localID + this->baseNumThreads * i]
            = _checkInCodelets11695(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets11695[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets11695[localID + this->baseNumThreads * i].decDep();
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
TP11694::~TP11694()
{
    delete[] L1_darts11694;
    delete[] L2_darts11694;
    delete[] q_darts11694;
    delete[] u31_darts11694;
    delete[] barrierCodelets11694;
    delete[] checkInCodelets11695;
}
/*TP11843: OMPForDirective*/
void TP11843::_barrierCodelets11843::fire(void)
{
    TP11843* myTP = static_cast<TP11843*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets11843[0].decDep();
}
bool TP11843::requestNewRangeIterations11843(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet11843 * codeletID;
        int tempEndRange = rangePerCodelet11843 * (codeletID + 1);
        if (remainderRange11843 != 0) {
            if (codeletID < (uint32_t)remainderRange11843) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange11843;
                tempEndRange += remainderRange11843;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration11843;
        tempEndRange = tempEndRange * 1 + minIteration11843;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration11843 < lastIteration11843) {
            (this->inputsTPParent->i_darts11843[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts11843[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration11843;
        }
    }
    return isThereNewIteration;
}
void TP11843::_checkInCodelets11844::fire(void)
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
    this->inputsTPParent->L2_darts11843[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->L2_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jend1_darts11843[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jend1_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jst1_darts11843[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jst1_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts11843[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21j_darts11843[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21j_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21jm1_darts11843[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21jm1_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31j_darts11843[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31j_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31jm1_darts11843[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31jm1_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41j_darts11843[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41j_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41jm1_darts11843[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41jm1_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51j_darts11843[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51j_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51jm1_darts11843[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51jm1_darts10830[this->getID()]);

    /*printing node 11844: ForStmt*/
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
    int** L2 = &(this->inputsTPParent->L2_darts11843[this->getLocalID()]);
    (void)L2 /*OMP_SHARED_PRIVATE*/;
    int* i = &(this->inputsTPParent->i_darts11843[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts11843[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jend1 = &(this->inputsTPParent->jend1_darts11843[this->getLocalID()]);
    (void)jend1 /*OMP_SHARED_PRIVATE*/;
    int** jst1 = &(this->inputsTPParent->jst1_darts11843[this->getLocalID()]);
    (void)jst1 /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts11843[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts11843[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts11843[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21j = &(this->inputsTPParent->u21j_darts11843[this->getLocalID()]);
    (void)u21j /*OMP_SHARED_PRIVATE*/;
    double** u21jm1 = &(this->inputsTPParent->u21jm1_darts11843[this->getLocalID()]);
    (void)u21jm1 /*OMP_SHARED_PRIVATE*/;
    double** u31j = &(this->inputsTPParent->u31j_darts11843[this->getLocalID()]);
    (void)u31j /*OMP_SHARED_PRIVATE*/;
    double** u31jm1 = &(this->inputsTPParent->u31jm1_darts11843[this->getLocalID()]);
    (void)u31jm1 /*OMP_SHARED_PRIVATE*/;
    double** u41j = &(this->inputsTPParent->u41j_darts11843[this->getLocalID()]);
    (void)u41j /*OMP_SHARED_PRIVATE*/;
    double** u41jm1 = &(this->inputsTPParent->u41jm1_darts11843[this->getLocalID()]);
    (void)u41jm1 /*OMP_SHARED_PRIVATE*/;
    double** u51j = &(this->inputsTPParent->u51j_darts11843[this->getLocalID()]);
    (void)u51j /*OMP_SHARED_PRIVATE*/;
    double** u51jm1 = &(this->inputsTPParent->u51jm1_darts11843[this->getLocalID()]);
    (void)u51jm1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations11843(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets11843[0].decDep();
        return;
    }
    for (int i_darts_counter_temp11843 = (*i); i_darts_counter_temp11843 <= endRange
         && i_darts_counter_temp11843 <= this->inputsTPParent->lastIteration11843;
         i_darts_counter_temp11843++) {
        {
            {
                /*Loop's init*/
                (*k) = 1;
                int k_darts_counter_temp11843 = (*k);
                for (; k_darts_counter_temp11843 <= nz - 2; k_darts_counter_temp11843++) {
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp11843 = (*j);
                        for (; j_darts_counter_temp11843 <= jend; j_darts_counter_temp11843++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp11843 = (*m);
                                for (; m_darts_counter_temp11843 < 5; m_darts_counter_temp11843++) {
                                    rsd[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                                       [k_darts_counter_temp11843][m_darts_counter_temp11843]
                                        = rsd[(i_darts_counter_temp11843)]
                                             [j_darts_counter_temp11843][k_darts_counter_temp11843]
                                             [m_darts_counter_temp11843]
                                        - ty2
                                            * (flux[(i_darts_counter_temp11843)]
                                                   [j_darts_counter_temp11843 + 1]
                                                   [k_darts_counter_temp11843]
                                                   [m_darts_counter_temp11843]
                                                - flux[(i_darts_counter_temp11843)]
                                                      [j_darts_counter_temp11843 - 1]
                                                      [k_darts_counter_temp11843]
                                                      [m_darts_counter_temp11843]);
                                }
                                (*m) = m_darts_counter_temp11843;
                            }
                        }
                        (*j) = j_darts_counter_temp11843;
                    }
                    (*(*L2)) = ny - 1;
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp11843 = (*j);
                        for (; j_darts_counter_temp11843 <= (*(*L2)); j_darts_counter_temp11843++) {
                            (*(*tmp)) = 1.
                                / u[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                                   [k_darts_counter_temp11843][0];
                            (*(*u21j)) = (*(*tmp))
                                * u[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                                   [k_darts_counter_temp11843][1];
                            (*(*u31j)) = (*(*tmp))
                                * u[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                                   [k_darts_counter_temp11843][2];
                            (*(*u41j)) = (*(*tmp))
                                * u[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                                   [k_darts_counter_temp11843][3];
                            (*(*u51j)) = (*(*tmp))
                                * u[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                                   [k_darts_counter_temp11843][4];
                            (*(*tmp)) = 1.
                                / u[(i_darts_counter_temp11843)][j_darts_counter_temp11843 - 1]
                                   [k_darts_counter_temp11843][0];
                            (*(*u21jm1)) = (*(*tmp))
                                * u[(i_darts_counter_temp11843)][j_darts_counter_temp11843 - 1]
                                   [k_darts_counter_temp11843][1];
                            (*(*u31jm1)) = (*(*tmp))
                                * u[(i_darts_counter_temp11843)][j_darts_counter_temp11843 - 1]
                                   [k_darts_counter_temp11843][2];
                            (*(*u41jm1)) = (*(*tmp))
                                * u[(i_darts_counter_temp11843)][j_darts_counter_temp11843 - 1]
                                   [k_darts_counter_temp11843][3];
                            (*(*u51jm1)) = (*(*tmp))
                                * u[(i_darts_counter_temp11843)][j_darts_counter_temp11843 - 1]
                                   [k_darts_counter_temp11843][4];
                            flux[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                                [k_darts_counter_temp11843][1]
                                = ty3 * ((*(*u21j)) - (*(*u21jm1)));
                            flux[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                                [k_darts_counter_temp11843][2]
                                = (4. / 3.) * ty3 * ((*(*u31j)) - (*(*u31jm1)));
                            flux[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                                [k_darts_counter_temp11843][3]
                                = ty3 * ((*(*u41j)) - (*(*u41jm1)));
                            flux[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                                [k_darts_counter_temp11843][4]
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
                        (*j) = j_darts_counter_temp11843;
                    }
                    {
                        /*Loop's init*/
                        (*j) = jst;
                        int j_darts_counter_temp11843 = (*j);
                        for (; j_darts_counter_temp11843 <= jend; j_darts_counter_temp11843++) {
                            rsd[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                               [k_darts_counter_temp11843][0]
                                = rsd[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                                     [k_darts_counter_temp11843][0]
                                + dy1 * ty1
                                    * (u[(i_darts_counter_temp11843)][j_darts_counter_temp11843 - 1]
                                        [k_darts_counter_temp11843][0]
                                        - 2.
                                            * u[(i_darts_counter_temp11843)]
                                               [j_darts_counter_temp11843]
                                               [k_darts_counter_temp11843][0]
                                        + u[(i_darts_counter_temp11843)][j_darts_counter_temp11843
                                            + 1][k_darts_counter_temp11843][0]);
                            rsd[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                               [k_darts_counter_temp11843][1]
                                = rsd[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                                     [k_darts_counter_temp11843][1]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp11843)][j_darts_counter_temp11843
                                           + 1][k_darts_counter_temp11843][1]
                                        - flux[(i_darts_counter_temp11843)]
                                              [j_darts_counter_temp11843][k_darts_counter_temp11843]
                                              [1])
                                + dy2 * ty1
                                    * (u[(i_darts_counter_temp11843)][j_darts_counter_temp11843 - 1]
                                        [k_darts_counter_temp11843][1]
                                        - 2.
                                            * u[(i_darts_counter_temp11843)]
                                               [j_darts_counter_temp11843]
                                               [k_darts_counter_temp11843][1]
                                        + u[(i_darts_counter_temp11843)][j_darts_counter_temp11843
                                            + 1][k_darts_counter_temp11843][1]);
                            rsd[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                               [k_darts_counter_temp11843][2]
                                = rsd[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                                     [k_darts_counter_temp11843][2]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp11843)][j_darts_counter_temp11843
                                           + 1][k_darts_counter_temp11843][2]
                                        - flux[(i_darts_counter_temp11843)]
                                              [j_darts_counter_temp11843][k_darts_counter_temp11843]
                                              [2])
                                + dy3 * ty1
                                    * (u[(i_darts_counter_temp11843)][j_darts_counter_temp11843 - 1]
                                        [k_darts_counter_temp11843][2]
                                        - 2.
                                            * u[(i_darts_counter_temp11843)]
                                               [j_darts_counter_temp11843]
                                               [k_darts_counter_temp11843][2]
                                        + u[(i_darts_counter_temp11843)][j_darts_counter_temp11843
                                            + 1][k_darts_counter_temp11843][2]);
                            rsd[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                               [k_darts_counter_temp11843][3]
                                = rsd[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                                     [k_darts_counter_temp11843][3]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp11843)][j_darts_counter_temp11843
                                           + 1][k_darts_counter_temp11843][3]
                                        - flux[(i_darts_counter_temp11843)]
                                              [j_darts_counter_temp11843][k_darts_counter_temp11843]
                                              [3])
                                + dy4 * ty1
                                    * (u[(i_darts_counter_temp11843)][j_darts_counter_temp11843 - 1]
                                        [k_darts_counter_temp11843][3]
                                        - 2.
                                            * u[(i_darts_counter_temp11843)]
                                               [j_darts_counter_temp11843]
                                               [k_darts_counter_temp11843][3]
                                        + u[(i_darts_counter_temp11843)][j_darts_counter_temp11843
                                            + 1][k_darts_counter_temp11843][3]);
                            rsd[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                               [k_darts_counter_temp11843][4]
                                = rsd[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                                     [k_darts_counter_temp11843][4]
                                + ty3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp11843)][j_darts_counter_temp11843
                                           + 1][k_darts_counter_temp11843][4]
                                        - flux[(i_darts_counter_temp11843)]
                                              [j_darts_counter_temp11843][k_darts_counter_temp11843]
                                              [4])
                                + dy5 * ty1
                                    * (u[(i_darts_counter_temp11843)][j_darts_counter_temp11843 - 1]
                                        [k_darts_counter_temp11843][4]
                                        - 2.
                                            * u[(i_darts_counter_temp11843)]
                                               [j_darts_counter_temp11843]
                                               [k_darts_counter_temp11843][4]
                                        + u[(i_darts_counter_temp11843)][j_darts_counter_temp11843
                                            + 1][k_darts_counter_temp11843][4]);
                        }
                        (*j) = j_darts_counter_temp11843;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp11843 = (*m);
                        for (; m_darts_counter_temp11843 < 5; m_darts_counter_temp11843++) {
                            rsd[(i_darts_counter_temp11843)][1][k_darts_counter_temp11843]
                               [m_darts_counter_temp11843]
                                = rsd[(i_darts_counter_temp11843)][1][k_darts_counter_temp11843]
                                     [m_darts_counter_temp11843]
                                - dssp
                                    * (+5.
                                            * u[(i_darts_counter_temp11843)][1]
                                               [k_darts_counter_temp11843]
                                               [m_darts_counter_temp11843]
                                        - 4.
                                            * u[(i_darts_counter_temp11843)][2]
                                               [k_darts_counter_temp11843]
                                               [m_darts_counter_temp11843]
                                        + u[(i_darts_counter_temp11843)][3]
                                           [k_darts_counter_temp11843][m_darts_counter_temp11843]);
                            rsd[(i_darts_counter_temp11843)][2][k_darts_counter_temp11843]
                               [m_darts_counter_temp11843]
                                = rsd[(i_darts_counter_temp11843)][2][k_darts_counter_temp11843]
                                     [m_darts_counter_temp11843]
                                - dssp
                                    * (-4.
                                            * u[(i_darts_counter_temp11843)][1]
                                               [k_darts_counter_temp11843]
                                               [m_darts_counter_temp11843]
                                        + 6.
                                            * u[(i_darts_counter_temp11843)][2]
                                               [k_darts_counter_temp11843]
                                               [m_darts_counter_temp11843]
                                        - 4.
                                            * u[(i_darts_counter_temp11843)][3]
                                               [k_darts_counter_temp11843]
                                               [m_darts_counter_temp11843]
                                        + u[(i_darts_counter_temp11843)][4]
                                           [k_darts_counter_temp11843][m_darts_counter_temp11843]);
                        }
                        (*m) = m_darts_counter_temp11843;
                    }
                    (*(*jst1)) = 3;
                    (*(*jend1)) = ny - 4;
                    {
                        /*Loop's init*/
                        (*j) = (*(*jst1));
                        int j_darts_counter_temp11843 = (*j);
                        for (; j_darts_counter_temp11843 <= (*(*jend1));
                             j_darts_counter_temp11843++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp11843 = (*m);
                                for (; m_darts_counter_temp11843 < 5; m_darts_counter_temp11843++) {
                                    rsd[(i_darts_counter_temp11843)][j_darts_counter_temp11843]
                                       [k_darts_counter_temp11843][m_darts_counter_temp11843]
                                        = rsd[(i_darts_counter_temp11843)]
                                             [j_darts_counter_temp11843][k_darts_counter_temp11843]
                                             [m_darts_counter_temp11843]
                                        - dssp
                                            * (u[(i_darts_counter_temp11843)]
                                                [j_darts_counter_temp11843 - 2]
                                                [k_darts_counter_temp11843]
                                                [m_darts_counter_temp11843]
                                                - 4.
                                                    * u[(i_darts_counter_temp11843)]
                                                       [j_darts_counter_temp11843 - 1]
                                                       [k_darts_counter_temp11843]
                                                       [m_darts_counter_temp11843]
                                                + 6.
                                                    * u[(i_darts_counter_temp11843)]
                                                       [j_darts_counter_temp11843]
                                                       [k_darts_counter_temp11843]
                                                       [m_darts_counter_temp11843]
                                                - 4.
                                                    * u[(i_darts_counter_temp11843)]
                                                       [j_darts_counter_temp11843 + 1]
                                                       [k_darts_counter_temp11843]
                                                       [m_darts_counter_temp11843]
                                                + u[(i_darts_counter_temp11843)]
                                                   [j_darts_counter_temp11843 + 2]
                                                   [k_darts_counter_temp11843]
                                                   [m_darts_counter_temp11843]);
                                }
                                (*m) = m_darts_counter_temp11843;
                            }
                        }
                        (*j) = j_darts_counter_temp11843;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp11843 = (*m);
                        for (; m_darts_counter_temp11843 < 5; m_darts_counter_temp11843++) {
                            rsd[(i_darts_counter_temp11843)][ny - 3][k_darts_counter_temp11843]
                               [m_darts_counter_temp11843]
                                = rsd[(i_darts_counter_temp11843)][ny - 3]
                                     [k_darts_counter_temp11843][m_darts_counter_temp11843]
                                - dssp
                                    * (u[(i_darts_counter_temp11843)][ny - 5]
                                        [k_darts_counter_temp11843][m_darts_counter_temp11843]
                                        - 4.
                                            * u[(i_darts_counter_temp11843)][ny - 4]
                                               [k_darts_counter_temp11843]
                                               [m_darts_counter_temp11843]
                                        + 6.
                                            * u[(i_darts_counter_temp11843)][ny - 3]
                                               [k_darts_counter_temp11843]
                                               [m_darts_counter_temp11843]
                                        - 4.
                                            * u[(i_darts_counter_temp11843)][ny - 2]
                                               [k_darts_counter_temp11843]
                                               [m_darts_counter_temp11843]);
                            rsd[(i_darts_counter_temp11843)][ny - 2][k_darts_counter_temp11843]
                               [m_darts_counter_temp11843]
                                = rsd[(i_darts_counter_temp11843)][ny - 2]
                                     [k_darts_counter_temp11843][m_darts_counter_temp11843]
                                - dssp
                                    * (u[(i_darts_counter_temp11843)][ny - 4]
                                        [k_darts_counter_temp11843][m_darts_counter_temp11843]
                                        - 4.
                                            * u[(i_darts_counter_temp11843)][ny - 3]
                                               [k_darts_counter_temp11843]
                                               [m_darts_counter_temp11843]
                                        + 5.
                                            * u[(i_darts_counter_temp11843)][ny - 2]
                                               [k_darts_counter_temp11843]
                                               [m_darts_counter_temp11843]);
                        }
                        (*m) = m_darts_counter_temp11843;
                    }
                }
                (*k) = k_darts_counter_temp11843;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets11843[0].decDep();
}
TP11843::TP11843(int in_numThreads, int in_mainCodeletID, TP10830* in_TPParent,
    int in_initIteration, int in_lastIteration, TP11843** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , L2_darts11843(new int*[this->numThreads])
    , i_darts11843(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts11843(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jend1_darts11843(new int*[this->numThreads])
    , jst1_darts11843(new int*[this->numThreads])
    , k_darts11843(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts11843(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts11843(new double*[this->numThreads])
    , u21j_darts11843(new double*[this->numThreads])
    , u21jm1_darts11843(new double*[this->numThreads])
    , u31j_darts11843(new double*[this->numThreads])
    , u31jm1_darts11843(new double*[this->numThreads])
    , u41j_darts11843(new double*[this->numThreads])
    , u41jm1_darts11843(new double*[this->numThreads])
    , u51j_darts11843(new double*[this->numThreads])
    , u51jm1_darts11843(new double*[this->numThreads])
    , initIteration11843(in_initIteration)
    , lastIteration11843(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets11843(new _barrierCodelets11843[1])
    , checkInCodelets11844(new _checkInCodelets11844[this->numThreads])
{
    /*Initialize the loop parameters*/
    range11843 = abs(lastIteration11843 - initIteration11843) / 1;
    rangePerCodelet11843 = range11843 / numThreads;
    minIteration11843 = min<int>(lastIteration11843, initIteration11843);
    remainderRange11843 = range11843 % numThreads;
    /*Initialize inputs and vars.*/
    this->L2_darts11843 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jend1_darts11843 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jst1_darts11843 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts11843
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21j_darts11843
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21jm1_darts11843
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31j_darts11843
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31jm1_darts11843
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41j_darts11843
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41jm1_darts11843
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51j_darts11843
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51jm1_darts11843
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets11843[0] = _barrierCodelets11843(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets11844* checkInCodelets11844Ptr = (this->checkInCodelets11844);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets11844);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets11844Ptr) = _checkInCodelets11844(2, 1, this, codeletCounter);
#else
        (*checkInCodelets11844Ptr) = _checkInCodelets11844(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets11844Ptr).decDep();
        checkInCodelets11844Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP11843::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets11844[localID].setID(codeletID);
    this->checkInCodelets11844[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets11844[localID + this->baseNumThreads * i]
            = _checkInCodelets11844(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets11844[localID + this->baseNumThreads * i]
            = _checkInCodelets11844(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets11844[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets11844[localID + this->baseNumThreads * i].decDep();
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
TP11843::~TP11843()
{
    delete[] L2_darts11843;
    delete[] jend1_darts11843;
    delete[] jst1_darts11843;
    delete[] tmp_darts11843;
    delete[] u21j_darts11843;
    delete[] u21jm1_darts11843;
    delete[] u31j_darts11843;
    delete[] u31jm1_darts11843;
    delete[] u41j_darts11843;
    delete[] u41jm1_darts11843;
    delete[] u51j_darts11843;
    delete[] u51jm1_darts11843;
    delete[] barrierCodelets11843;
    delete[] checkInCodelets11844;
}
/*TP12480: OMPForDirective*/
void TP12480::_barrierCodelets12480::fire(void)
{
    TP12480* myTP = static_cast<TP12480*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets12480[0].decDep();
}
bool TP12480::requestNewRangeIterations12480(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet12480 * codeletID;
        int tempEndRange = rangePerCodelet12480 * (codeletID + 1);
        if (remainderRange12480 != 0) {
            if (codeletID < (uint32_t)remainderRange12480) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange12480;
                tempEndRange += remainderRange12480;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration12480;
        tempEndRange = tempEndRange * 1 + minIteration12480;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration12480 < lastIteration12480) {
            (this->inputsTPParent->i_darts12480[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts12480[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration12480;
        }
    }
    return isThereNewIteration;
}
void TP12480::_checkInCodelets12481::fire(void)
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
    this->inputsTPParent->q_darts12480[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->q_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->tmp_darts12480[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->tmp_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21k_darts12480[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21k_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u21km1_darts12480[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u21km1_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31k_darts12480[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31k_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u31km1_darts12480[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u31km1_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41_darts12480[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41k_darts12480[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41k_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u41km1_darts12480[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u41km1_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51k_darts12480[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51k_darts10830[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->u51km1_darts12480[this->getLocalID()]
        = (double*)&(myTP->TPParent->inputsTPParent->u51km1_darts10830[this->getID()]);

    /*printing node 12481: ForStmt*/
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
    int* i = &(this->inputsTPParent->i_darts12480[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts12480[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts12480[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts12480[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double** q = &(this->inputsTPParent->q_darts12480[this->getLocalID()]);
    (void)q /*OMP_SHARED_PRIVATE*/;
    double** tmp = &(this->inputsTPParent->tmp_darts12480[this->getLocalID()]);
    (void)tmp /*OMP_SHARED_PRIVATE*/;
    double** u21k = &(this->inputsTPParent->u21k_darts12480[this->getLocalID()]);
    (void)u21k /*OMP_SHARED_PRIVATE*/;
    double** u21km1 = &(this->inputsTPParent->u21km1_darts12480[this->getLocalID()]);
    (void)u21km1 /*OMP_SHARED_PRIVATE*/;
    double** u31k = &(this->inputsTPParent->u31k_darts12480[this->getLocalID()]);
    (void)u31k /*OMP_SHARED_PRIVATE*/;
    double** u31km1 = &(this->inputsTPParent->u31km1_darts12480[this->getLocalID()]);
    (void)u31km1 /*OMP_SHARED_PRIVATE*/;
    double** u41 = &(this->inputsTPParent->u41_darts12480[this->getLocalID()]);
    (void)u41 /*OMP_SHARED_PRIVATE*/;
    double** u41k = &(this->inputsTPParent->u41k_darts12480[this->getLocalID()]);
    (void)u41k /*OMP_SHARED_PRIVATE*/;
    double** u41km1 = &(this->inputsTPParent->u41km1_darts12480[this->getLocalID()]);
    (void)u41km1 /*OMP_SHARED_PRIVATE*/;
    double** u51k = &(this->inputsTPParent->u51k_darts12480[this->getLocalID()]);
    (void)u51k /*OMP_SHARED_PRIVATE*/;
    double** u51km1 = &(this->inputsTPParent->u51km1_darts12480[this->getLocalID()]);
    (void)u51km1 /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations12480(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets12480[0].decDep();
        return;
    }
    for (int i_darts_counter_temp12480 = (*i); i_darts_counter_temp12480 <= endRange
         && i_darts_counter_temp12480 <= this->inputsTPParent->lastIteration12480;
         i_darts_counter_temp12480++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp12480 = (*j);
                for (; j_darts_counter_temp12480 <= jend; j_darts_counter_temp12480++) {
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp12480 = (*k);
                        for (; k_darts_counter_temp12480 <= nz - 1; k_darts_counter_temp12480++) {
                            flux[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                [k_darts_counter_temp12480][0]
                                = u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                   [k_darts_counter_temp12480][3];
                            (*(*u41)) = u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                         [k_darts_counter_temp12480][3]
                                / u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                   [k_darts_counter_temp12480][0];
                            (*(*q)) = 0.5
                                * (u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                    [k_darts_counter_temp12480][1]
                                        * u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                           [k_darts_counter_temp12480][1]
                                    + u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                       [k_darts_counter_temp12480][2]
                                        * u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                           [k_darts_counter_temp12480][2]
                                    + u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                       [k_darts_counter_temp12480][3]
                                        * u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                           [k_darts_counter_temp12480][3])
                                / u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                   [k_darts_counter_temp12480][0];
                            flux[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                [k_darts_counter_temp12480][1]
                                = u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                   [k_darts_counter_temp12480][1]
                                * (*(*u41));
                            flux[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                [k_darts_counter_temp12480][2]
                                = u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                   [k_darts_counter_temp12480][2]
                                * (*(*u41));
                            flux[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                [k_darts_counter_temp12480][3]
                                = u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                   [k_darts_counter_temp12480][3]
                                    * (*(*u41))
                                + 0.40000000000000002
                                    * (u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                        [k_darts_counter_temp12480][4]
                                        - (*(*q)));
                            flux[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                [k_darts_counter_temp12480][4]
                                = (1.3999999999999999
                                          * u[(i_darts_counter_temp12480)]
                                             [j_darts_counter_temp12480][k_darts_counter_temp12480]
                                             [4]
                                      - 0.40000000000000002 * (*(*q)))
                                * (*(*u41));
                        }
                        (*k) = k_darts_counter_temp12480;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp12480 = (*k);
                        for (; k_darts_counter_temp12480 <= nz - 2; k_darts_counter_temp12480++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp12480 = (*m);
                                for (; m_darts_counter_temp12480 < 5; m_darts_counter_temp12480++) {
                                    rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                       [k_darts_counter_temp12480][m_darts_counter_temp12480]
                                        = rsd[(i_darts_counter_temp12480)]
                                             [j_darts_counter_temp12480][k_darts_counter_temp12480]
                                             [m_darts_counter_temp12480]
                                        - tz2
                                            * (flux[(i_darts_counter_temp12480)]
                                                   [j_darts_counter_temp12480]
                                                   [k_darts_counter_temp12480 + 1]
                                                   [m_darts_counter_temp12480]
                                                - flux[(i_darts_counter_temp12480)]
                                                      [j_darts_counter_temp12480]
                                                      [k_darts_counter_temp12480 - 1]
                                                      [m_darts_counter_temp12480]);
                                }
                                (*m) = m_darts_counter_temp12480;
                            }
                        }
                        (*k) = k_darts_counter_temp12480;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp12480 = (*k);
                        for (; k_darts_counter_temp12480 <= nz - 1; k_darts_counter_temp12480++) {
                            (*(*tmp)) = 1.
                                / u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                   [k_darts_counter_temp12480][0];
                            (*(*u21k)) = (*(*tmp))
                                * u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                   [k_darts_counter_temp12480][1];
                            (*(*u31k)) = (*(*tmp))
                                * u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                   [k_darts_counter_temp12480][2];
                            (*(*u41k)) = (*(*tmp))
                                * u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                   [k_darts_counter_temp12480][3];
                            (*(*u51k)) = (*(*tmp))
                                * u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                   [k_darts_counter_temp12480][4];
                            (*(*tmp)) = 1.
                                / u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                   [k_darts_counter_temp12480 - 1][0];
                            (*(*u21km1)) = (*(*tmp))
                                * u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                   [k_darts_counter_temp12480 - 1][1];
                            (*(*u31km1)) = (*(*tmp))
                                * u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                   [k_darts_counter_temp12480 - 1][2];
                            (*(*u41km1)) = (*(*tmp))
                                * u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                   [k_darts_counter_temp12480 - 1][3];
                            (*(*u51km1)) = (*(*tmp))
                                * u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                   [k_darts_counter_temp12480 - 1][4];
                            flux[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                [k_darts_counter_temp12480][1]
                                = tz3 * ((*(*u21k)) - (*(*u21km1)));
                            flux[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                [k_darts_counter_temp12480][2]
                                = tz3 * ((*(*u31k)) - (*(*u31km1)));
                            flux[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                [k_darts_counter_temp12480][3]
                                = (4. / 3.) * tz3 * ((*(*u41k)) - (*(*u41km1)));
                            flux[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                [k_darts_counter_temp12480][4]
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
                        (*k) = k_darts_counter_temp12480;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp12480 = (*k);
                        for (; k_darts_counter_temp12480 <= nz - 2; k_darts_counter_temp12480++) {
                            rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                               [k_darts_counter_temp12480][0]
                                = rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                     [k_darts_counter_temp12480][0]
                                + dz1 * tz1
                                    * (u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                        [k_darts_counter_temp12480 - 1][0]
                                        - 2.
                                            * u[(i_darts_counter_temp12480)]
                                               [j_darts_counter_temp12480]
                                               [k_darts_counter_temp12480][0]
                                        + u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                           [k_darts_counter_temp12480 + 1][0]);
                            rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                               [k_darts_counter_temp12480][1]
                                = rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                     [k_darts_counter_temp12480][1]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                           [k_darts_counter_temp12480 + 1][1]
                                        - flux[(i_darts_counter_temp12480)]
                                              [j_darts_counter_temp12480][k_darts_counter_temp12480]
                                              [1])
                                + dz2 * tz1
                                    * (u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                        [k_darts_counter_temp12480 - 1][1]
                                        - 2.
                                            * u[(i_darts_counter_temp12480)]
                                               [j_darts_counter_temp12480]
                                               [k_darts_counter_temp12480][1]
                                        + u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                           [k_darts_counter_temp12480 + 1][1]);
                            rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                               [k_darts_counter_temp12480][2]
                                = rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                     [k_darts_counter_temp12480][2]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                           [k_darts_counter_temp12480 + 1][2]
                                        - flux[(i_darts_counter_temp12480)]
                                              [j_darts_counter_temp12480][k_darts_counter_temp12480]
                                              [2])
                                + dz3 * tz1
                                    * (u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                        [k_darts_counter_temp12480 - 1][2]
                                        - 2.
                                            * u[(i_darts_counter_temp12480)]
                                               [j_darts_counter_temp12480]
                                               [k_darts_counter_temp12480][2]
                                        + u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                           [k_darts_counter_temp12480 + 1][2]);
                            rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                               [k_darts_counter_temp12480][3]
                                = rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                     [k_darts_counter_temp12480][3]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                           [k_darts_counter_temp12480 + 1][3]
                                        - flux[(i_darts_counter_temp12480)]
                                              [j_darts_counter_temp12480][k_darts_counter_temp12480]
                                              [3])
                                + dz4 * tz1
                                    * (u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                        [k_darts_counter_temp12480 - 1][3]
                                        - 2.
                                            * u[(i_darts_counter_temp12480)]
                                               [j_darts_counter_temp12480]
                                               [k_darts_counter_temp12480][3]
                                        + u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                           [k_darts_counter_temp12480 + 1][3]);
                            rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                               [k_darts_counter_temp12480][4]
                                = rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                     [k_darts_counter_temp12480][4]
                                + tz3 * 0.10000000000000001 * 1.
                                    * (flux[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                           [k_darts_counter_temp12480 + 1][4]
                                        - flux[(i_darts_counter_temp12480)]
                                              [j_darts_counter_temp12480][k_darts_counter_temp12480]
                                              [4])
                                + dz5 * tz1
                                    * (u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                        [k_darts_counter_temp12480 - 1][4]
                                        - 2.
                                            * u[(i_darts_counter_temp12480)]
                                               [j_darts_counter_temp12480]
                                               [k_darts_counter_temp12480][4]
                                        + u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                           [k_darts_counter_temp12480 + 1][4]);
                        }
                        (*k) = k_darts_counter_temp12480;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp12480 = (*m);
                        for (; m_darts_counter_temp12480 < 5; m_darts_counter_temp12480++) {
                            rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480][1]
                               [m_darts_counter_temp12480]
                                = rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480][1]
                                     [m_darts_counter_temp12480]
                                - dssp
                                    * (+5.
                                            * u[(i_darts_counter_temp12480)]
                                               [j_darts_counter_temp12480][1]
                                               [m_darts_counter_temp12480]
                                        - 4.
                                            * u[(i_darts_counter_temp12480)]
                                               [j_darts_counter_temp12480][2]
                                               [m_darts_counter_temp12480]
                                        + u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                           [3][m_darts_counter_temp12480]);
                            rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480][2]
                               [m_darts_counter_temp12480]
                                = rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480][2]
                                     [m_darts_counter_temp12480]
                                - dssp
                                    * (-4.
                                            * u[(i_darts_counter_temp12480)]
                                               [j_darts_counter_temp12480][1]
                                               [m_darts_counter_temp12480]
                                        + 6.
                                            * u[(i_darts_counter_temp12480)]
                                               [j_darts_counter_temp12480][2]
                                               [m_darts_counter_temp12480]
                                        - 4.
                                            * u[(i_darts_counter_temp12480)]
                                               [j_darts_counter_temp12480][3]
                                               [m_darts_counter_temp12480]
                                        + u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                           [4][m_darts_counter_temp12480]);
                        }
                        (*m) = m_darts_counter_temp12480;
                    }
                    {
                        /*Loop's init*/
                        (*k) = 3;
                        int k_darts_counter_temp12480 = (*k);
                        for (; k_darts_counter_temp12480 <= nz - 4; k_darts_counter_temp12480++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp12480 = (*m);
                                for (; m_darts_counter_temp12480 < 5; m_darts_counter_temp12480++) {
                                    rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                       [k_darts_counter_temp12480][m_darts_counter_temp12480]
                                        = rsd[(i_darts_counter_temp12480)]
                                             [j_darts_counter_temp12480][k_darts_counter_temp12480]
                                             [m_darts_counter_temp12480]
                                        - dssp
                                            * (u[(i_darts_counter_temp12480)]
                                                [j_darts_counter_temp12480]
                                                [k_darts_counter_temp12480 - 2]
                                                [m_darts_counter_temp12480]
                                                - 4.
                                                    * u[(i_darts_counter_temp12480)]
                                                       [j_darts_counter_temp12480]
                                                       [k_darts_counter_temp12480 - 1]
                                                       [m_darts_counter_temp12480]
                                                + 6.
                                                    * u[(i_darts_counter_temp12480)]
                                                       [j_darts_counter_temp12480]
                                                       [k_darts_counter_temp12480]
                                                       [m_darts_counter_temp12480]
                                                - 4.
                                                    * u[(i_darts_counter_temp12480)]
                                                       [j_darts_counter_temp12480]
                                                       [k_darts_counter_temp12480 + 1]
                                                       [m_darts_counter_temp12480]
                                                + u[(i_darts_counter_temp12480)]
                                                   [j_darts_counter_temp12480]
                                                   [k_darts_counter_temp12480 + 2]
                                                   [m_darts_counter_temp12480]);
                                }
                                (*m) = m_darts_counter_temp12480;
                            }
                        }
                        (*k) = k_darts_counter_temp12480;
                    }
                    {
                        /*Loop's init*/
                        (*m) = 0;
                        int m_darts_counter_temp12480 = (*m);
                        for (; m_darts_counter_temp12480 < 5; m_darts_counter_temp12480++) {
                            rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480][nz - 3]
                               [m_darts_counter_temp12480]
                                = rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                     [nz - 3][m_darts_counter_temp12480]
                                - dssp
                                    * (u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                        [nz - 5][m_darts_counter_temp12480]
                                        - 4.
                                            * u[(i_darts_counter_temp12480)]
                                               [j_darts_counter_temp12480][nz - 4]
                                               [m_darts_counter_temp12480]
                                        + 6.
                                            * u[(i_darts_counter_temp12480)]
                                               [j_darts_counter_temp12480][nz - 3]
                                               [m_darts_counter_temp12480]
                                        - 4.
                                            * u[(i_darts_counter_temp12480)]
                                               [j_darts_counter_temp12480][nz - 2]
                                               [m_darts_counter_temp12480]);
                            rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480][nz - 2]
                               [m_darts_counter_temp12480]
                                = rsd[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                     [nz - 2][m_darts_counter_temp12480]
                                - dssp
                                    * (u[(i_darts_counter_temp12480)][j_darts_counter_temp12480]
                                        [nz - 4][m_darts_counter_temp12480]
                                        - 4.
                                            * u[(i_darts_counter_temp12480)]
                                               [j_darts_counter_temp12480][nz - 3]
                                               [m_darts_counter_temp12480]
                                        + 5.
                                            * u[(i_darts_counter_temp12480)]
                                               [j_darts_counter_temp12480][nz - 2]
                                               [m_darts_counter_temp12480]);
                        }
                        (*m) = m_darts_counter_temp12480;
                    }
                }
                (*j) = j_darts_counter_temp12480;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets12480[0].decDep();
}
TP12480::TP12480(int in_numThreads, int in_mainCodeletID, TP10830* in_TPParent,
    int in_initIteration, int in_lastIteration, TP12480** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts12480(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts12480(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts12480(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts12480(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , q_darts12480(new double*[this->numThreads])
    , tmp_darts12480(new double*[this->numThreads])
    , u21k_darts12480(new double*[this->numThreads])
    , u21km1_darts12480(new double*[this->numThreads])
    , u31k_darts12480(new double*[this->numThreads])
    , u31km1_darts12480(new double*[this->numThreads])
    , u41_darts12480(new double*[this->numThreads])
    , u41k_darts12480(new double*[this->numThreads])
    , u41km1_darts12480(new double*[this->numThreads])
    , u51k_darts12480(new double*[this->numThreads])
    , u51km1_darts12480(new double*[this->numThreads])
    , initIteration12480(in_initIteration)
    , lastIteration12480(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets12480(new _barrierCodelets12480[1])
    , checkInCodelets12481(new _checkInCodelets12481[this->numThreads])
{
    /*Initialize the loop parameters*/
    range12480 = abs(lastIteration12480 - initIteration12480) / 1;
    rangePerCodelet12480 = range12480 / numThreads;
    minIteration12480 = min<int>(lastIteration12480, initIteration12480);
    remainderRange12480 = range12480 % numThreads;
    /*Initialize inputs and vars.*/
    this->q_darts12480
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->tmp_darts12480
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21k_darts12480
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u21km1_darts12480
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31k_darts12480
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u31km1_darts12480
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41_darts12480
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41k_darts12480
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u41km1_darts12480
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51k_darts12480
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->u51km1_darts12480
        = (double**)malloc(sizeof(double*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets12480[0] = _barrierCodelets12480(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets12481* checkInCodelets12481Ptr = (this->checkInCodelets12481);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets12481);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets12481Ptr) = _checkInCodelets12481(2, 1, this, codeletCounter);
#else
        (*checkInCodelets12481Ptr) = _checkInCodelets12481(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets12481Ptr).decDep();
        checkInCodelets12481Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP12480::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets12481[localID].setID(codeletID);
    this->checkInCodelets12481[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets12481[localID + this->baseNumThreads * i]
            = _checkInCodelets12481(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets12481[localID + this->baseNumThreads * i]
            = _checkInCodelets12481(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets12481[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets12481[localID + this->baseNumThreads * i].decDep();
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
TP12480::~TP12480()
{
    delete[] q_darts12480;
    delete[] tmp_darts12480;
    delete[] u21k_darts12480;
    delete[] u21km1_darts12480;
    delete[] u31k_darts12480;
    delete[] u31km1_darts12480;
    delete[] u41_darts12480;
    delete[] u41k_darts12480;
    delete[] u41km1_darts12480;
    delete[] u51k_darts12480;
    delete[] u51km1_darts12480;
    delete[] barrierCodelets12480;
    delete[] checkInCodelets12481;
}
/*TP13231: OMPParallelDirective*/
void TP13231::_barrierCodelets13231::fire(void)
{
    TP13231* myTP = static_cast<TP13231*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP13231::_checkInCodelets13235::fire(void)
{
    /*region 13235 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13235;
    if (idx < myTP->TPsToUse13235) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13235_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(nx - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13235;
            int minIteration = min<int>(nx, 0);
            int remainderRange = range % myTP->TPsToUse13235;
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
            if (idx == myTP->TPsToUse13235 - 1) {
                lastIteration = nx;
            }
#if USEINVOKE == 1
            invoke<TP13235>(myTP, myTP->codeletsPerTP13235 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13235Ptr[idx]));
#else
            place<TP13235>(idx, myTP, myTP->codeletsPerTP13235 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13235Ptr[idx]));
#endif
        } else {
            if (myTP->TP13235Ptr[idx] != nullptr) {
                myTP->TP13235Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13231::_barrierCodelets13235::fire(void)
{
    TP13231* myTP = static_cast<TP13231*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets13286[codeletsCounter].decDep();
        }
    }
}
void TP13231::_checkInCodelets13286::fire(void)
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
void TP13231::_barrierCodelets13286::fire(void)
{
    TP13231* myTP = static_cast<TP13231*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets13328[codeletsCounter].decDep();
        }
    }
}
void TP13231::_checkInCodelets13328::fire(void)
{
    /*region 13328 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13328;
    if (idx < myTP->TPsToUse13328) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13328_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(nx - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13328;
            int minIteration = min<int>(nx, 0);
            int remainderRange = range % myTP->TPsToUse13328;
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
            if (idx == myTP->TPsToUse13328 - 1) {
                lastIteration = nx;
            }
#if USEINVOKE == 1
            invoke<TP13328>(myTP, myTP->codeletsPerTP13328 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13328Ptr[idx]));
#else
            place<TP13328>(idx, myTP, myTP->codeletsPerTP13328 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13328Ptr[idx]));
#endif
        } else {
            if (myTP->TP13328Ptr[idx] != nullptr) {
                myTP->TP13328Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13231::_barrierCodelets13328::fire(void)
{
    TP13231* myTP = static_cast<TP13231*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets13372[codeletsCounter].decDep();
        }
    }
}
void TP13231::_checkInCodelets13372::fire(void)
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
void TP13231::_barrierCodelets13372::fire(void)
{
    TP13231* myTP = static_cast<TP13231*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets13414[codeletsCounter].decDep();
        }
    }
}
void TP13231::_checkInCodelets13414::fire(void)
{
    /*region 13414 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13414;
    if (idx < myTP->TPsToUse13414) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13414_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(ny - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13414;
            int minIteration = min<int>(ny, 0);
            int remainderRange = range % myTP->TPsToUse13414;
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
            if (idx == myTP->TPsToUse13414 - 1) {
                lastIteration = ny;
            }
#if USEINVOKE == 1
            invoke<TP13414>(myTP, myTP->codeletsPerTP13414 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13414Ptr[idx]));
#else
            place<TP13414>(idx, myTP, myTP->codeletsPerTP13414 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13414Ptr[idx]));
#endif
        } else {
            if (myTP->TP13414Ptr[idx] != nullptr) {
                myTP->TP13414Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13231::_barrierCodelets13414::fire(void)
{
    TP13231* myTP = static_cast<TP13231*>(myTP_);
    myTP->TPParent->barrierCodelets13231[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets13231[0]));
}
TP13231::TP13231(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13231(new int[this->numThreads]) /*VARIABLE*/
    , iglob_darts13231(new int[this->numThreads]) /*VARIABLE*/
    , j_darts13231(new int[this->numThreads]) /*VARIABLE*/
    , jglob_darts13231(new int[this->numThreads]) /*VARIABLE*/
    , k_darts13231(new int[this->numThreads]) /*VARIABLE*/
    , TP13235Ptr(new TP13235*[NUMTPS13235])
    , TP13235_alreadyLaunched(new size_t[NUMTPS13235])
    , numTPsSet13235(0)
    , numTPsReady13235(0)
    , TPsToUse13235(NUMTPS13235)
    , codeletsPerTP13235(this->numThreads / NUMTPS13235)
    , totalCodelets13235(this->TPsToUse13235 * this->codeletsPerTP13235)
    , TP13286Ptr(new TP13286*[NUMTPS13286])
    , TP13286_alreadyLaunched(new size_t[NUMTPS13286])
    , numTPsSet13286(0)
    , numTPsReady13286(0)
    , TPsToUse13286(NUMTPS13286)
    , codeletsPerTP13286(this->numThreads / NUMTPS13286)
    , totalCodelets13286(this->TPsToUse13286 * this->codeletsPerTP13286)
    , TP13328Ptr(new TP13328*[NUMTPS13328])
    , TP13328_alreadyLaunched(new size_t[NUMTPS13328])
    , numTPsSet13328(0)
    , numTPsReady13328(0)
    , TPsToUse13328(NUMTPS13328)
    , codeletsPerTP13328(this->numThreads / NUMTPS13328)
    , totalCodelets13328(this->TPsToUse13328 * this->codeletsPerTP13328)
    , TP13372Ptr(new TP13372*[NUMTPS13372])
    , TP13372_alreadyLaunched(new size_t[NUMTPS13372])
    , numTPsSet13372(0)
    , numTPsReady13372(0)
    , TPsToUse13372(NUMTPS13372)
    , codeletsPerTP13372(this->numThreads / NUMTPS13372)
    , totalCodelets13372(this->TPsToUse13372 * this->codeletsPerTP13372)
    , TP13414Ptr(new TP13414*[NUMTPS13414])
    , TP13414_alreadyLaunched(new size_t[NUMTPS13414])
    , numTPsSet13414(0)
    , numTPsReady13414(0)
    , TPsToUse13414(NUMTPS13414)
    , codeletsPerTP13414(this->numThreads / NUMTPS13414)
    , totalCodelets13414(this->TPsToUse13414 * this->codeletsPerTP13414)
    , barrierCodelets13231(new _barrierCodelets13231[1])
    , checkInCodelets13235(new _checkInCodelets13235[this->numThreads])
    , barrierCodelets13235(new _barrierCodelets13235[1])
    , checkInCodelets13286(new _checkInCodelets13286[this->numThreads])
    , barrierCodelets13286(new _barrierCodelets13286[1])
    , checkInCodelets13328(new _checkInCodelets13328[this->numThreads])
    , barrierCodelets13328(new _barrierCodelets13328[1])
    , checkInCodelets13372(new _checkInCodelets13372[this->numThreads])
    , barrierCodelets13372(new _barrierCodelets13372[1])
    , checkInCodelets13414(new _checkInCodelets13414[this->numThreads])
    , barrierCodelets13414(new _barrierCodelets13414[1])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets13231[0] = _barrierCodelets13231(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets13414[0] = _barrierCodelets13414(NUMTPS13414, NUMTPS13414, this, 0);
    barrierCodelets13372[0] = _barrierCodelets13372(NUMTPS13372, NUMTPS13372, this, 0);
    barrierCodelets13328[0] = _barrierCodelets13328(NUMTPS13328, NUMTPS13328, this, 0);
    barrierCodelets13286[0] = _barrierCodelets13286(NUMTPS13286, NUMTPS13286, this, 0);
    barrierCodelets13235[0] = _barrierCodelets13235(NUMTPS13235, NUMTPS13235, this, 0);
    _checkInCodelets13414* checkInCodelets13414Ptr = (this->checkInCodelets13414);
    for (int i = 0; i < NUMTPS13414; i++) {
        TP13414Ptr[i] = nullptr;
        TP13414_alreadyLaunched[i] = 0;
    }
    _checkInCodelets13372* checkInCodelets13372Ptr = (this->checkInCodelets13372);
    for (int i = 0; i < NUMTPS13372; i++) {
        TP13372Ptr[i] = nullptr;
        TP13372_alreadyLaunched[i] = 0;
    }
    _checkInCodelets13328* checkInCodelets13328Ptr = (this->checkInCodelets13328);
    for (int i = 0; i < NUMTPS13328; i++) {
        TP13328Ptr[i] = nullptr;
        TP13328_alreadyLaunched[i] = 0;
    }
    _checkInCodelets13286* checkInCodelets13286Ptr = (this->checkInCodelets13286);
    for (int i = 0; i < NUMTPS13286; i++) {
        TP13286Ptr[i] = nullptr;
        TP13286_alreadyLaunched[i] = 0;
    }
    _checkInCodelets13235* checkInCodelets13235Ptr = (this->checkInCodelets13235);
    for (int i = 0; i < NUMTPS13235; i++) {
        TP13235Ptr[i] = nullptr;
        TP13235_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets13414Ptr) = _checkInCodelets13414(1, 1, this, codeletCounter);
        checkInCodelets13414Ptr++;
        (*checkInCodelets13372Ptr) = _checkInCodelets13372(1, 1, this, codeletCounter);
        checkInCodelets13372Ptr++;
        (*checkInCodelets13328Ptr) = _checkInCodelets13328(1, 1, this, codeletCounter);
        checkInCodelets13328Ptr++;
        (*checkInCodelets13286Ptr) = _checkInCodelets13286(1, 1, this, codeletCounter);
        checkInCodelets13286Ptr++;
        (*checkInCodelets13235Ptr) = _checkInCodelets13235(1, 1, this, codeletCounter);
        (*checkInCodelets13235Ptr).decDep();
        checkInCodelets13235Ptr++;
    }
}
TP13231::~TP13231()
{
    delete[] i_darts13231;
    delete[] iglob_darts13231;
    delete[] j_darts13231;
    delete[] jglob_darts13231;
    delete[] k_darts13231;
    delete[] barrierCodelets13231;
    delete[] barrierCodelets13414;
    delete[] checkInCodelets13414;
    delete[] barrierCodelets13372;
    delete[] checkInCodelets13372;
    delete[] barrierCodelets13328;
    delete[] checkInCodelets13328;
    delete[] barrierCodelets13286;
    delete[] checkInCodelets13286;
    delete[] barrierCodelets13235;
    delete[] checkInCodelets13235;
}
/*TP13235: OMPForDirective*/
void TP13235::_barrierCodelets13235::fire(void)
{
    TP13235* myTP = static_cast<TP13235*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13235[0].decDep();
}
bool TP13235::requestNewRangeIterations13235(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13235 * codeletID;
        int tempEndRange = rangePerCodelet13235 * (codeletID + 1);
        if (remainderRange13235 != 0) {
            if (codeletID < (uint32_t)remainderRange13235) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13235;
                tempEndRange += remainderRange13235;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13235;
        tempEndRange = tempEndRange * 1 + minIteration13235;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13235 < lastIteration13235) {
            (this->inputsTPParent->i_darts13235[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts13235[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13235;
        }
    }
    return isThereNewIteration;
}
void TP13235::_checkInCodelets13236::fire(void)
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
    this->inputsTPParent->iglob_darts13235[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts13231[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->jglob_darts13235[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts13231[this->getID()]);

    /*printing node 13236: ForStmt*/
    /*var: i*/
    /*var: iglob*/
    /*var: j*/
    /*var: jglob*/
    int* i = &(this->inputsTPParent->i_darts13235[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iglob = &(this->inputsTPParent->iglob_darts13235[this->getLocalID()]);
    (void)iglob /*OMP_SHARED_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts13235[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jglob = &(this->inputsTPParent->jglob_darts13235[this->getLocalID()]);
    (void)jglob /*OMP_SHARED_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13235(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13235[0].decDep();
        return;
    }
    for (int i_darts_counter_temp13235 = (*i); i_darts_counter_temp13235 < endRange
         && i_darts_counter_temp13235 < this->inputsTPParent->lastIteration13235;
         i_darts_counter_temp13235++) {
        {
            (*(*iglob)) = (i_darts_counter_temp13235);
            {
                /*Loop's init*/
                (*j) = 0;
                int j_darts_counter_temp13235 = (*j);
                for (; j_darts_counter_temp13235 < ny; j_darts_counter_temp13235++) {
                    (*(*jglob)) = j_darts_counter_temp13235;
                    exact((*(*iglob)), (*(*jglob)), 0,
                        &u[(i_darts_counter_temp13235)][j_darts_counter_temp13235][0][0]);
                    exact((*(*iglob)), (*(*jglob)), nz - 1,
                        &u[(i_darts_counter_temp13235)][j_darts_counter_temp13235][nz - 1][0]);
                }
                (*j) = j_darts_counter_temp13235;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13235[0].decDep();
}
TP13235::TP13235(int in_numThreads, int in_mainCodeletID, TP13231* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13235** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13235(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iglob_darts13235(new int*[this->numThreads])
    , j_darts13235(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jglob_darts13235(new int*[this->numThreads])
    , initIteration13235(in_initIteration)
    , lastIteration13235(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13235(new _barrierCodelets13235[1])
    , checkInCodelets13236(new _checkInCodelets13236[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13235 = abs(lastIteration13235 - initIteration13235) / 1;
    rangePerCodelet13235 = range13235 / numThreads;
    minIteration13235 = min<int>(lastIteration13235, initIteration13235);
    remainderRange13235 = range13235 % numThreads;
    /*Initialize inputs and vars.*/
    this->iglob_darts13235 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->jglob_darts13235 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets13235[0] = _barrierCodelets13235(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13236* checkInCodelets13236Ptr = (this->checkInCodelets13236);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13236);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13236Ptr) = _checkInCodelets13236(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13236Ptr) = _checkInCodelets13236(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13236Ptr).decDep();
        checkInCodelets13236Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13235::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13236[localID].setID(codeletID);
    this->checkInCodelets13236[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13236[localID + this->baseNumThreads * i]
            = _checkInCodelets13236(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13236[localID + this->baseNumThreads * i]
            = _checkInCodelets13236(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13236[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13236[localID + this->baseNumThreads * i].decDep();
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
TP13235::~TP13235()
{
    delete[] iglob_darts13235;
    delete[] jglob_darts13235;
    delete[] barrierCodelets13235;
    delete[] checkInCodelets13236;
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
        = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts13231[this->getID()]);

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
                    exact((*(*iglob)), 0, k_darts_counter_temp13286,
                        &u[(i_darts_counter_temp13286)][0][k_darts_counter_temp13286][0]);
                }
                (*k) = k_darts_counter_temp13286;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13286[0].decDep();
}
TP13286::TP13286(int in_numThreads, int in_mainCodeletID, TP13231* in_TPParent,
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
/*TP13328: OMPForDirective*/
void TP13328::_barrierCodelets13328::fire(void)
{
    TP13328* myTP = static_cast<TP13328*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13328[0].decDep();
}
bool TP13328::requestNewRangeIterations13328(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13328 * codeletID;
        int tempEndRange = rangePerCodelet13328 * (codeletID + 1);
        if (remainderRange13328 != 0) {
            if (codeletID < (uint32_t)remainderRange13328) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13328;
                tempEndRange += remainderRange13328;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13328;
        tempEndRange = tempEndRange * 1 + minIteration13328;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13328 < lastIteration13328) {
            (this->inputsTPParent->i_darts13328[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts13328[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13328;
        }
    }
    return isThereNewIteration;
}
void TP13328::_checkInCodelets13329::fire(void)
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
    this->inputsTPParent->iglob_darts13328[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->iglob_darts13231[this->getID()]);

    /*printing node 13329: ForStmt*/
    /*var: i*/
    /*var: iglob*/
    /*var: k*/
    int* i = &(this->inputsTPParent->i_darts13328[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int** iglob = &(this->inputsTPParent->iglob_darts13328[this->getLocalID()]);
    (void)iglob /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts13328[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13328(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13328[0].decDep();
        return;
    }
    for (int i_darts_counter_temp13328 = (*i); i_darts_counter_temp13328 < endRange
         && i_darts_counter_temp13328 < this->inputsTPParent->lastIteration13328;
         i_darts_counter_temp13328++) {
        {
            (*(*iglob)) = (i_darts_counter_temp13328);
            {
                /*Loop's init*/
                (*k) = 0;
                int k_darts_counter_temp13328 = (*k);
                for (; k_darts_counter_temp13328 < nz; k_darts_counter_temp13328++) {
                    exact((*(*iglob)), ny0 - 1, k_darts_counter_temp13328,
                        &u[(i_darts_counter_temp13328)][ny - 1][k_darts_counter_temp13328][0]);
                }
                (*k) = k_darts_counter_temp13328;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13328[0].decDep();
}
TP13328::TP13328(int in_numThreads, int in_mainCodeletID, TP13231* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13328** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13328(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , iglob_darts13328(new int*[this->numThreads])
    , k_darts13328(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration13328(in_initIteration)
    , lastIteration13328(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13328(new _barrierCodelets13328[1])
    , checkInCodelets13329(new _checkInCodelets13329[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13328 = abs(lastIteration13328 - initIteration13328) / 1;
    rangePerCodelet13328 = range13328 / numThreads;
    minIteration13328 = min<int>(lastIteration13328, initIteration13328);
    remainderRange13328 = range13328 % numThreads;
    /*Initialize inputs and vars.*/
    this->iglob_darts13328 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets13328[0] = _barrierCodelets13328(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13329* checkInCodelets13329Ptr = (this->checkInCodelets13329);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13329);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13329Ptr) = _checkInCodelets13329(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13329Ptr) = _checkInCodelets13329(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13329Ptr).decDep();
        checkInCodelets13329Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13328::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13329[localID].setID(codeletID);
    this->checkInCodelets13329[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13329[localID + this->baseNumThreads * i]
            = _checkInCodelets13329(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13329[localID + this->baseNumThreads * i]
            = _checkInCodelets13329(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13329[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13329[localID + this->baseNumThreads * i].decDep();
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
TP13328::~TP13328()
{
    delete[] iglob_darts13328;
    delete[] barrierCodelets13328;
    delete[] checkInCodelets13329;
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
        = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts13231[this->getID()]);

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
                    exact(0, (*(*jglob)), k_darts_counter_temp13372,
                        &u[0][(j_darts_counter_temp13372)][k_darts_counter_temp13372][0]);
                }
                (*k) = k_darts_counter_temp13372;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13372[0].decDep();
}
TP13372::TP13372(int in_numThreads, int in_mainCodeletID, TP13231* in_TPParent,
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
/*TP13414: OMPForDirective*/
void TP13414::_barrierCodelets13414::fire(void)
{
    TP13414* myTP = static_cast<TP13414*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13414[0].decDep();
}
bool TP13414::requestNewRangeIterations13414(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13414 * codeletID;
        int tempEndRange = rangePerCodelet13414 * (codeletID + 1);
        if (remainderRange13414 != 0) {
            if (codeletID < (uint32_t)remainderRange13414) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13414;
                tempEndRange += remainderRange13414;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13414;
        tempEndRange = tempEndRange * 1 + minIteration13414;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13414 < lastIteration13414) {
            (this->inputsTPParent->j_darts13414[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->j_darts13414[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13414;
        }
    }
    return isThereNewIteration;
}
void TP13414::_checkInCodelets13415::fire(void)
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
    this->inputsTPParent->jglob_darts13414[this->getLocalID()]
        = (int*)&(myTP->TPParent->inputsTPParent->jglob_darts13231[this->getID()]);

    /*printing node 13415: ForStmt*/
    /*var: j*/
    /*var: jglob*/
    /*var: k*/
    int* j = &(this->inputsTPParent->j_darts13414[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int** jglob = &(this->inputsTPParent->jglob_darts13414[this->getLocalID()]);
    (void)jglob /*OMP_SHARED_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts13414[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13414(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13414[0].decDep();
        return;
    }
    for (int j_darts_counter_temp13414 = (*j); j_darts_counter_temp13414 < endRange
         && j_darts_counter_temp13414 < this->inputsTPParent->lastIteration13414;
         j_darts_counter_temp13414++) {
        {
            (*(*jglob)) = (j_darts_counter_temp13414);
            {
                /*Loop's init*/
                (*k) = 0;
                int k_darts_counter_temp13414 = (*k);
                for (; k_darts_counter_temp13414 < nz; k_darts_counter_temp13414++) {
                    exact(nx0 - 1, (*(*jglob)), k_darts_counter_temp13414,
                        &u[nx - 1][(j_darts_counter_temp13414)][k_darts_counter_temp13414][0]);
                }
                (*k) = k_darts_counter_temp13414;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13414[0].decDep();
}
TP13414::TP13414(int in_numThreads, int in_mainCodeletID, TP13231* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13414** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , j_darts13414(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , jglob_darts13414(new int*[this->numThreads])
    , k_darts13414(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration13414(in_initIteration)
    , lastIteration13414(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13414(new _barrierCodelets13414[1])
    , checkInCodelets13415(new _checkInCodelets13415[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13414 = abs(lastIteration13414 - initIteration13414) / 1;
    rangePerCodelet13414 = range13414 / numThreads;
    minIteration13414 = min<int>(lastIteration13414, initIteration13414);
    remainderRange13414 = range13414 % numThreads;
    /*Initialize inputs and vars.*/
    this->jglob_darts13414 = (int**)malloc(sizeof(int*) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    barrierCodelets13414[0] = _barrierCodelets13414(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13415* checkInCodelets13415Ptr = (this->checkInCodelets13415);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13415);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13415Ptr) = _checkInCodelets13415(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13415Ptr) = _checkInCodelets13415(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13415Ptr).decDep();
        checkInCodelets13415Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13414::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13415[localID].setID(codeletID);
    this->checkInCodelets13415[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13415[localID + this->baseNumThreads * i]
            = _checkInCodelets13415(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13415[localID + this->baseNumThreads * i]
            = _checkInCodelets13415(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13415[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13415[localID + this->baseNumThreads * i].decDep();
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
TP13414::~TP13414()
{
    delete[] jglob_darts13414;
    delete[] barrierCodelets13414;
    delete[] checkInCodelets13415;
}
/*TP13903: OMPParallelDirective*/
void TP13903::_barrierCodelets13903::fire(void)
{
    TP13903* myTP = static_cast<TP13903*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP13903::_checkInCodelets13905::fire(void)
{
    /*region 13905 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP13905;
    if (idx < myTP->TPsToUse13905) {
        if (!__sync_val_compare_and_swap(&(myTP->TP13905_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(12 - 0) / 1;
            int rangePerCodelet = range / myTP->TPsToUse13905;
            int minIteration = min<int>(12, 0);
            int remainderRange = range % myTP->TPsToUse13905;
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
            if (idx == myTP->TPsToUse13905 - 1) {
                lastIteration = 12;
            }
#if USEINVOKE == 1
            invoke<TP13905>(myTP, myTP->codeletsPerTP13905 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13905Ptr[idx]));
#else
            place<TP13905>(idx, myTP, myTP->codeletsPerTP13905 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP13905Ptr[idx]));
#endif
        } else {
            if (myTP->TP13905Ptr[idx] != nullptr) {
                myTP->TP13905Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP13903::_barrierCodelets13905::fire(void)
{
    TP13903* myTP = static_cast<TP13903*>(myTP_);
    myTP->TPParent->barrierCodelets13903[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets13903[0]));
}
TP13903::TP13903(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13903(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts13903(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts13903(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts13903(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , TP13905Ptr(new TP13905*[NUMTPS13905])
    , TP13905_alreadyLaunched(new size_t[NUMTPS13905])
    , numTPsSet13905(0)
    , numTPsReady13905(0)
    , TPsToUse13905(NUMTPS13905)
    , codeletsPerTP13905(this->numThreads / NUMTPS13905)
    , totalCodelets13905(this->TPsToUse13905 * this->codeletsPerTP13905)
    , barrierCodelets13903(new _barrierCodelets13903[1])
    , checkInCodelets13905(new _checkInCodelets13905[this->numThreads])
    , barrierCodelets13905(new _barrierCodelets13905[1])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets13903[0] = _barrierCodelets13903(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets13905[0] = _barrierCodelets13905(NUMTPS13905, NUMTPS13905, this, 0);
    _checkInCodelets13905* checkInCodelets13905Ptr = (this->checkInCodelets13905);
    for (int i = 0; i < NUMTPS13905; i++) {
        TP13905Ptr[i] = nullptr;
        TP13905_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets13905Ptr) = _checkInCodelets13905(1, 1, this, codeletCounter);
        (*checkInCodelets13905Ptr).decDep();
        checkInCodelets13905Ptr++;
    }
}
TP13903::~TP13903()
{
    delete[] barrierCodelets13903;
    delete[] barrierCodelets13905;
    delete[] checkInCodelets13905;
}
/*TP13905: OMPForDirective*/
void TP13905::_barrierCodelets13905::fire(void)
{
    TP13905* myTP = static_cast<TP13905*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets13905[0].decDep();
}
bool TP13905::requestNewRangeIterations13905(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet13905 * codeletID;
        int tempEndRange = rangePerCodelet13905 * (codeletID + 1);
        if (remainderRange13905 != 0) {
            if (codeletID < (uint32_t)remainderRange13905) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange13905;
                tempEndRange += remainderRange13905;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration13905;
        tempEndRange = tempEndRange * 1 + minIteration13905;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration13905 < lastIteration13905) {
            (this->inputsTPParent->i_darts13905[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts13905[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration13905;
        }
    }
    return isThereNewIteration;
}
void TP13905::_checkInCodelets13906::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 13906: ForStmt*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    int* i = &(this->inputsTPParent->i_darts13905[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts13905[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts13905[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts13905[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations13905(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets13905[0].decDep();
        return;
    }
    for (int i_darts_counter_temp13905 = (*i); i_darts_counter_temp13905 < endRange
         && i_darts_counter_temp13905 < this->inputsTPParent->lastIteration13905;
         i_darts_counter_temp13905++) {
        {
            {
                /*Loop's init*/
                (*j) = 0;
                int j_darts_counter_temp13905 = (*j);
                for (; j_darts_counter_temp13905 < 12; j_darts_counter_temp13905++) {
                    {
                        /*Loop's init*/
                        (*k) = 0;
                        int k_darts_counter_temp13905 = (*k);
                        for (; k_darts_counter_temp13905 < 5; k_darts_counter_temp13905++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp13905 = (*m);
                                for (; m_darts_counter_temp13905 < 5; m_darts_counter_temp13905++) {
                                    a[(i_darts_counter_temp13905)][j_darts_counter_temp13905]
                                     [k_darts_counter_temp13905][m_darts_counter_temp13905]
                                        = 0.;
                                    b[(i_darts_counter_temp13905)][j_darts_counter_temp13905]
                                     [k_darts_counter_temp13905][m_darts_counter_temp13905]
                                        = 0.;
                                    c[(i_darts_counter_temp13905)][j_darts_counter_temp13905]
                                     [k_darts_counter_temp13905][m_darts_counter_temp13905]
                                        = 0.;
                                    d[(i_darts_counter_temp13905)][j_darts_counter_temp13905]
                                     [k_darts_counter_temp13905][m_darts_counter_temp13905]
                                        = 0.;
                                }
                                (*m) = m_darts_counter_temp13905;
                            }
                        }
                        (*k) = k_darts_counter_temp13905;
                    }
                }
                (*j) = j_darts_counter_temp13905;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets13905[0].decDep();
}
TP13905::TP13905(int in_numThreads, int in_mainCodeletID, TP13903* in_TPParent,
    int in_initIteration, int in_lastIteration, TP13905** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts13905(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts13905(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts13905(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts13905(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration13905(in_initIteration)
    , lastIteration13905(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets13905(new _barrierCodelets13905[1])
    , checkInCodelets13906(new _checkInCodelets13906[this->numThreads])
{
    /*Initialize the loop parameters*/
    range13905 = abs(lastIteration13905 - initIteration13905) / 1;
    rangePerCodelet13905 = range13905 / numThreads;
    minIteration13905 = min<int>(lastIteration13905, initIteration13905);
    remainderRange13905 = range13905 % numThreads;
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets13905[0] = _barrierCodelets13905(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets13906* checkInCodelets13906Ptr = (this->checkInCodelets13906);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets13906);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets13906Ptr) = _checkInCodelets13906(2, 1, this, codeletCounter);
#else
        (*checkInCodelets13906Ptr) = _checkInCodelets13906(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets13906Ptr).decDep();
        checkInCodelets13906Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP13905::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets13906[localID].setID(codeletID);
    this->checkInCodelets13906[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets13906[localID + this->baseNumThreads * i]
            = _checkInCodelets13906(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets13906[localID + this->baseNumThreads * i]
            = _checkInCodelets13906(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets13906[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets13906[localID + this->baseNumThreads * i].decDep();
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
TP13905::~TP13905()
{
    delete[] barrierCodelets13905;
    delete[] checkInCodelets13906;
}
/*TP14002: OMPParallelDirective*/
void TP14002::_barrierCodelets14002::fire(void)
{
    TP14002* myTP = static_cast<TP14002*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP14002::_checkInCodelets14004::fire(void)
{
    /*region 14004 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP14004;
    if (idx < myTP->TPsToUse14004) {
        if (!__sync_val_compare_and_swap(&(myTP->TP14004_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse14004;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse14004;
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
            if (idx == myTP->TPsToUse14004 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse14004 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP14004>(myTP, myTP->codeletsPerTP14004 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP14004Ptr[idx]));
#else
            place<TP14004>(idx, myTP, myTP->codeletsPerTP14004 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(myTP->TP14004Ptr[idx]));
#endif
        } else {
            if (myTP->TP14004Ptr[idx] != nullptr) {
                myTP->TP14004Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP14002::_barrierCodelets14004::fire(void)
{
    TP14002* myTP = static_cast<TP14002*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets14060[codeletsCounter].decDep();
        }
    }
}
void TP14002::_checkInCodelets14060::fire(void)
{

    /*printing node 14060: BinaryOperator*/
    (this->inputsTPParent->k_darts14002[this->getID()]) = 1;

    /*printing node 14061: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts14002[this->getID()]) <= nz - 2) {
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets14059[this->getID()].decDep();
        return;
    } else {
        /*Signal the codelet after the loop from the end condional node.*/
        /*Signaling next codelet region: 14063 nextRegion: 14071 */
        myTP->controlTPParent->barrierCodelets14071[0].decDep();
        return;
    }
}
void TP14002::_checkInCodelets14059::fire(void)
{

    /*printing node 14059: ForStmt*/
    bool haveToLaunch = __sync_bool_compare_and_swap(&(myTP->controlTPParent->TP14059_LoopCounter),
        myTP->controlTPParent->TP14059_LoopCounterPerThread[this->getID()],
        myTP->controlTPParent->TP14059_LoopCounterPerThread[this->getID()] + 1);
    unsigned int iterIdx = myTP->controlTPParent->TP14059_LoopCounterPerThread[this->getID()];
    if (haveToLaunch) {
        this->resetCodelet();
        myTP->controlTPParent->TP14059PtrVec.push_back(nullptr);
        myTP->controlTPParent->TP14059_LoopCounterPerThread[this->getID()] += 1;
        invoke<TP14059>(myTP, myTP->numThreads, this->getID(), myTP, myTP->inputsTPParent,
            &(myTP->controlTPParent->TP14059PtrVec.back()));
    } else {
        if (myTP->controlTPParent->TP14059PtrVec.size() == 0) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP14059PtrVec.size() < (iterIdx + 1)) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP14059PtrVec[iterIdx] == nullptr) {
            this->resetCodelet();
            this->decDep();
            return;
        } else {
            this->resetCodelet();
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14059PtrVec[iterIdx]->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14059PtrVec[iterIdx]->availableCodelets[this->getID()] = 1;
#endif
            myTP->controlTPParent->TP14059_LoopCounterPerThread[this->getID()] += 1;
        }
    }
}
void TP14002::_checkInCodelets14063::fire(void)
{

    /*printing node 14063: UnaryOperator*/
    (this->inputsTPParent->k_darts14002[this->getID()])++;

    /*printing node 14522: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts14002[this->getID()]) <= nz - 2) {
        this->resetCodelet();
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets14059[this->getID()].decDep();
        return;
    } else {
        /*Signal the codelet after the loop from the condtional node.*/
        /*Signaling next codelet region: 14063 nextRegion: 14071 */
        myTP->controlTPParent->barrierCodelets14071[0].decDep();
        return;
    }
}
void TP14002::_barrierCodelets14071::fire(void)
{
    TP14002* myTP = static_cast<TP14002*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets14073[codeletsCounter].decDep();
        }
    }
}
void TP14002::_checkInCodelets14073::fire(void)
{

    /*printing node 14073: BinaryOperator*/
    (this->inputsTPParent->k_darts14002[this->getID()]) = nz - 2;

    /*printing node 14075: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts14002[this->getID()]) >= 1) {
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets14072[this->getID()].decDep();
        return;
    } else {
        /*Signal the codelet after the loop from the end condional node.*/
        /*Signaling next codelet region: 14076 nextRegion: 14080 */
        myTP->controlTPParent->barrierCodelets14080[0].decDep();
        return;
    }
}
void TP14002::_checkInCodelets14072::fire(void)
{

    /*printing node 14072: ForStmt*/
    bool haveToLaunch = __sync_bool_compare_and_swap(&(myTP->controlTPParent->TP14072_LoopCounter),
        myTP->controlTPParent->TP14072_LoopCounterPerThread[this->getID()],
        myTP->controlTPParent->TP14072_LoopCounterPerThread[this->getID()] + 1);
    unsigned int iterIdx = myTP->controlTPParent->TP14072_LoopCounterPerThread[this->getID()];
    if (haveToLaunch) {
        this->resetCodelet();
        myTP->controlTPParent->TP14072PtrVec.push_back(nullptr);
        myTP->controlTPParent->TP14072_LoopCounterPerThread[this->getID()] += 1;
        invoke<TP14072>(myTP, myTP->numThreads, this->getID(), myTP, myTP->inputsTPParent,
            &(myTP->controlTPParent->TP14072PtrVec.back()));
    } else {
        if (myTP->controlTPParent->TP14072PtrVec.size() == 0) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP14072PtrVec.size() < (iterIdx + 1)) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP14072PtrVec[iterIdx] == nullptr) {
            this->resetCodelet();
            this->decDep();
            return;
        } else {
            this->resetCodelet();
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14072PtrVec[iterIdx]->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14072PtrVec[iterIdx]->availableCodelets[this->getID()] = 1;
#endif
            myTP->controlTPParent->TP14072_LoopCounterPerThread[this->getID()] += 1;
        }
    }
}
void TP14002::_checkInCodelets14076::fire(void)
{

    /*printing node 14076: UnaryOperator*/
    (this->inputsTPParent->k_darts14002[this->getID()])--;

    /*printing node 14523: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts14002[this->getID()]) >= 1) {
        this->resetCodelet();
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets14072[this->getID()].decDep();
        return;
    } else {
        /*Signal the codelet after the loop from the condtional node.*/
        /*Signaling next codelet region: 14076 nextRegion: 14080 */
        myTP->controlTPParent->barrierCodelets14080[0].decDep();
        return;
    }
}
void TP14002::_barrierCodelets14080::fire(void)
{
    TP14002* myTP = static_cast<TP14002*>(myTP_);
    {
        for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
            myTP->checkInCodelets14081[codeletsCounter].decDep();
        }
    }
}
void TP14002::_checkInCodelets14081::fire(void)
{
    /*region 14081 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP14081;
    if (idx < myTP->TPsToUse14081) {
        if (!__sync_val_compare_and_swap(&(myTP->TP14081_alreadyLaunched[idx]), 0, 1)) {
            int range = abs(iend - ist) / 1;
            int rangePerCodelet = range / myTP->TPsToUse14081;
            int minIteration = min<int>(iend, ist);
            int remainderRange = range % myTP->TPsToUse14081;
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
            if (idx == myTP->TPsToUse14081 - 1) {
                lastIteration = lastIteration + 1;
            }
            if (idx == myTP->TPsToUse14081 - 1) {
                lastIteration = iend;
            }
#if USEINVOKE == 1
            invoke<TP14081>(myTP, myTP->codeletsPerTP14081 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(*(this->inputsTPParent->tmp_darts14002)),
                &(myTP->TP14081Ptr[idx]));
#else
            place<TP14081>(idx, myTP, myTP->codeletsPerTP14081 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration, &(*(this->inputsTPParent->tmp_darts14002)),
                &(myTP->TP14081Ptr[idx]));
#endif
        } else {
            if (myTP->TP14081Ptr[idx] != nullptr) {
                myTP->TP14081Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    }
}
void TP14002::_barrierCodelets14081::fire(void)
{
    TP14002* myTP = static_cast<TP14002*>(myTP_);
    myTP->TPParent->barrierCodelets14002[0].setDep(0);
    myTP->add(&(myTP->TPParent->barrierCodelets14002[0]));
}
TP14002::TP14002(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet, double* in_tmp)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts14002(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , istep_darts14002(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts14002(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts14002(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts14002(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts14002(in_tmp) /*OMP_SHARED - INPUT*/
    , TP14004Ptr(new TP14004*[NUMTPS14004])
    , TP14004_alreadyLaunched(new size_t[NUMTPS14004])
    , numTPsSet14004(0)
    , numTPsReady14004(0)
    , TPsToUse14004(NUMTPS14004)
    , codeletsPerTP14004(this->numThreads / NUMTPS14004)
    , totalCodelets14004(this->TPsToUse14004 * this->codeletsPerTP14004)
    , TP14059_LoopCounter(0)
    , TP14059_LoopCounterPerThread(new unsigned int[this->numThreads])
    , TP14072_LoopCounter(0)
    , TP14072_LoopCounterPerThread(new unsigned int[this->numThreads])
    , TP14081Ptr(new TP14081*[NUMTPS14081])
    , TP14081_alreadyLaunched(new size_t[NUMTPS14081])
    , numTPsSet14081(0)
    , numTPsReady14081(0)
    , TPsToUse14081(NUMTPS14081)
    , codeletsPerTP14081(this->numThreads / NUMTPS14081)
    , totalCodelets14081(this->TPsToUse14081 * this->codeletsPerTP14081)
    , barrierCodelets14002(new _barrierCodelets14002[1])
    , checkInCodelets14004(new _checkInCodelets14004[this->numThreads])
    , barrierCodelets14004(new _barrierCodelets14004[1])
    , checkInCodelets14060(new _checkInCodelets14060[this->numThreads])
    , checkInCodelets14059(new _checkInCodelets14059[this->numThreads])
    , checkInCodelets14063(new _checkInCodelets14063[this->numThreads])
    , barrierCodelets14071(new _barrierCodelets14071[1])
    , checkInCodelets14073(new _checkInCodelets14073[this->numThreads])
    , checkInCodelets14072(new _checkInCodelets14072[this->numThreads])
    , checkInCodelets14076(new _checkInCodelets14076[this->numThreads])
    , barrierCodelets14080(new _barrierCodelets14080[1])
    , checkInCodelets14081(new _checkInCodelets14081[this->numThreads])
    , barrierCodelets14081(new _barrierCodelets14081[1])
{
    memset((void*)TP14059_LoopCounterPerThread, 0, this->numThreads * sizeof(unsigned int));
    memset((void*)TP14072_LoopCounterPerThread, 0, this->numThreads * sizeof(unsigned int));
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets14002[0] = _barrierCodelets14002(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets14081[0] = _barrierCodelets14081(NUMTPS14081, NUMTPS14081, this, 0);
    barrierCodelets14080[0] = _barrierCodelets14080(this->numThreads, this->numThreads, this, 0);
    barrierCodelets14071[0] = _barrierCodelets14071(this->numThreads, this->numThreads, this, 0);
    barrierCodelets14004[0] = _barrierCodelets14004(NUMTPS14004, NUMTPS14004, this, 0);
    _checkInCodelets14081* checkInCodelets14081Ptr = (this->checkInCodelets14081);
    for (int i = 0; i < NUMTPS14081; i++) {
        TP14081Ptr[i] = nullptr;
        TP14081_alreadyLaunched[i] = 0;
    }
    _checkInCodelets14076* checkInCodelets14076Ptr = (this->checkInCodelets14076);
    _checkInCodelets14072* checkInCodelets14072Ptr = (this->checkInCodelets14072);
    _checkInCodelets14073* checkInCodelets14073Ptr = (this->checkInCodelets14073);
    _checkInCodelets14063* checkInCodelets14063Ptr = (this->checkInCodelets14063);
    _checkInCodelets14059* checkInCodelets14059Ptr = (this->checkInCodelets14059);
    _checkInCodelets14060* checkInCodelets14060Ptr = (this->checkInCodelets14060);
    _checkInCodelets14004* checkInCodelets14004Ptr = (this->checkInCodelets14004);
    for (int i = 0; i < NUMTPS14004; i++) {
        TP14004Ptr[i] = nullptr;
        TP14004_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets14081Ptr) = _checkInCodelets14081(1, 1, this, codeletCounter);
        checkInCodelets14081Ptr++;
        (*checkInCodelets14076Ptr) = _checkInCodelets14076(1, 1, this, codeletCounter);
        checkInCodelets14076Ptr++;
        (*checkInCodelets14072Ptr) = _checkInCodelets14072(1, 1, this, codeletCounter);
        checkInCodelets14072Ptr++;
        (*checkInCodelets14073Ptr) = _checkInCodelets14073(1, 1, this, codeletCounter);
        checkInCodelets14073Ptr++;
        (*checkInCodelets14063Ptr) = _checkInCodelets14063(1, 1, this, codeletCounter);
        checkInCodelets14063Ptr++;
        (*checkInCodelets14059Ptr) = _checkInCodelets14059(1, 1, this, codeletCounter);
        checkInCodelets14059Ptr++;
        (*checkInCodelets14060Ptr) = _checkInCodelets14060(1, 1, this, codeletCounter);
        checkInCodelets14060Ptr++;
        (*checkInCodelets14004Ptr) = _checkInCodelets14004(1, 1, this, codeletCounter);
        (*checkInCodelets14004Ptr).decDep();
        checkInCodelets14004Ptr++;
    }
}
TP14002::~TP14002()
{
    delete[] TP14059_LoopCounterPerThread;
    delete[] TP14072_LoopCounterPerThread;
    delete[] barrierCodelets14002;
    delete[] barrierCodelets14081;
    delete[] checkInCodelets14081;
    delete[] barrierCodelets14080;
    delete[] checkInCodelets14076;
    delete[] checkInCodelets14072;
    delete[] checkInCodelets14073;
    delete[] barrierCodelets14071;
    delete[] checkInCodelets14063;
    delete[] checkInCodelets14059;
    delete[] checkInCodelets14060;
    delete[] barrierCodelets14004;
    delete[] checkInCodelets14004;
}
/*TP14004: OMPForDirective*/
void TP14004::_barrierCodelets14004::fire(void)
{
    TP14004* myTP = static_cast<TP14004*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets14004[0].decDep();
}
bool TP14004::requestNewRangeIterations14004(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet14004 * codeletID;
        int tempEndRange = rangePerCodelet14004 * (codeletID + 1);
        if (remainderRange14004 != 0) {
            if (codeletID < (uint32_t)remainderRange14004) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange14004;
                tempEndRange += remainderRange14004;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration14004;
        tempEndRange = tempEndRange * 1 + minIteration14004;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration14004 < lastIteration14004) {
            (this->inputsTPParent->i_darts14004[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts14004[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration14004;
        }
    }
    return isThereNewIteration;
}
void TP14004::_checkInCodelets14005::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 14005: ForStmt*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    int* i = &(this->inputsTPParent->i_darts14004[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts14004[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts14004[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts14004[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations14004(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets14004[0].decDep();
        return;
    }
    for (int i_darts_counter_temp14004 = (*i); i_darts_counter_temp14004 <= endRange
         && i_darts_counter_temp14004 <= this->inputsTPParent->lastIteration14004;
         i_darts_counter_temp14004++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp14004 = (*j);
                for (; j_darts_counter_temp14004 <= jend; j_darts_counter_temp14004++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp14004 = (*k);
                        for (; k_darts_counter_temp14004 <= nz - 2; k_darts_counter_temp14004++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp14004 = (*m);
                                for (; m_darts_counter_temp14004 < 5; m_darts_counter_temp14004++) {
                                    rsd[(i_darts_counter_temp14004)][j_darts_counter_temp14004]
                                       [k_darts_counter_temp14004][m_darts_counter_temp14004]
                                        = dt
                                        * rsd[(i_darts_counter_temp14004)]
                                             [j_darts_counter_temp14004][k_darts_counter_temp14004]
                                             [m_darts_counter_temp14004];
                                }
                                (*m) = m_darts_counter_temp14004;
                            }
                        }
                        (*k) = k_darts_counter_temp14004;
                    }
                }
                (*j) = j_darts_counter_temp14004;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets14004[0].decDep();
}
TP14004::TP14004(int in_numThreads, int in_mainCodeletID, TP14002* in_TPParent,
    int in_initIteration, int in_lastIteration, TP14004** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts14004(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts14004(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts14004(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts14004(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , initIteration14004(in_initIteration)
    , lastIteration14004(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets14004(new _barrierCodelets14004[1])
    , checkInCodelets14005(new _checkInCodelets14005[this->numThreads])
{
    /*Initialize the loop parameters*/
    range14004 = abs(lastIteration14004 - initIteration14004) / 1;
    rangePerCodelet14004 = range14004 / numThreads;
    minIteration14004 = min<int>(lastIteration14004, initIteration14004);
    remainderRange14004 = range14004 % numThreads;
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets14004[0] = _barrierCodelets14004(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets14005* checkInCodelets14005Ptr = (this->checkInCodelets14005);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets14005);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets14005Ptr) = _checkInCodelets14005(2, 1, this, codeletCounter);
#else
        (*checkInCodelets14005Ptr) = _checkInCodelets14005(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets14005Ptr).decDep();
        checkInCodelets14005Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP14004::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets14005[localID].setID(codeletID);
    this->checkInCodelets14005[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets14005[localID + this->baseNumThreads * i]
            = _checkInCodelets14005(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets14005[localID + this->baseNumThreads * i]
            = _checkInCodelets14005(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets14005[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets14005[localID + this->baseNumThreads * i].decDep();
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
TP14004::~TP14004()
{
    delete[] barrierCodelets14004;
    delete[] checkInCodelets14005;
}
/*TP14059: ForStmt*/
void TP14059::_checkInCodelets14065::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif

    /*printing node 14065: CallExpr*/
    if (!__sync_val_compare_and_swap(&(myTP->controlTPParent->TP14065_alreadyLaunched), 0, 1)) {
        /*Make the function call*/
        invoke<TP_jacld>(myTP, myTP->numThreads, this->getID(),
            &(myTP->controlTPParent->checkInCodelets14066[this->getID()]),
            &(myTP->controlTPParent->TPParent->barrierCodelets14071[0]),
            &(myTP->controlTPParent->TP14065Ptr),
            (this->inputsTPParent->k_darts14002[this->getID()]));
    } else {
        if (myTP->controlTPParent->TP14065Ptr == nullptr) {
            myTP->add(this);
            return;
        } else {
            myTP->controlTPParent->TP14065Ptr->setNewInputs(
                (this->inputsTPParent->k_darts14002[this->getID()]), this->getID());
            myTP->controlTPParent->TP14065Ptr->nextCodeletsjacld[this->getID()]
                = &(myTP->controlTPParent->checkInCodelets14066[this->getID()]);
            myTP->controlTPParent->TP14065Ptr->nextSyncCodeletsjacld[this->getID()]
                = &(myTP->controlTPParent->TPParent->barrierCodelets14071[0]);
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14065Ptr->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14065Ptr->availableCodelets[this->getID()] = 1;
#endif
        }
    }
}
void TP14059::_checkInCodelets14066::fire(void)
{
    if (this->getID() == 0) {
        /*Init the vars for this region*/
        /*Initialize the vars of the inlined region*/
        this->inputsTPParent->k_darts14066
            = &(this->inputsTPParent
                    ->k_darts14002[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;

        /*printing node 14067: CallExpr*/
        printf("jacld complete... #k = %d\n", (*(this->inputsTPParent->k_darts14066)));
        /*Signaling next codelet from last stmt in the codelet*/
        /*Signaling next codelet region: 14066 nextRegion: 14068 */
        myTP->controlTPParent->checkInCodelets14068[this->getID()].decDep();
    } else {
        /*Signaling next codelet region: 14066 nextRegion: 14068 */
        myTP->checkInCodelets14068[this->getID()].decDep();
    }
}
void TP14059::_checkInCodelets14068::fire(void)
{

    /*printing node 14068: CallExpr*/
    if (!__sync_val_compare_and_swap(&(myTP->controlTPParent->TP14068_alreadyLaunched), 0, 1)) {
        /*Make the function call*/
        invoke<TP_blts>(myTP, myTP->numThreads, this->getID(),
            &(myTP->controlTPParent->checkInCodelets14069[this->getID()]),
            &(myTP->controlTPParent->TPParent->barrierCodelets14071[0]),
            &(myTP->controlTPParent->TP14068Ptr), nx, ny, nz,
            (this->inputsTPParent->k_darts14002[this->getID()]), omega, ist, iend, jst, jend, nx0,
            ny0);
    } else {
        if (myTP->controlTPParent->TP14068Ptr == nullptr) {
            myTP->add(this);
            return;
        } else {
            myTP->controlTPParent->TP14068Ptr->setNewInputs(nx, ny, nz,
                (this->inputsTPParent->k_darts14002[this->getID()]), omega, ist, iend, jst, jend,
                nx0, ny0, this->getID());
            myTP->controlTPParent->TP14068Ptr->nextCodeletsblts[this->getID()]
                = &(myTP->controlTPParent->checkInCodelets14069[this->getID()]);
            myTP->controlTPParent->TP14068Ptr->nextSyncCodeletsblts[this->getID()]
                = &(myTP->controlTPParent->TPParent->barrierCodelets14071[0]);
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14068Ptr->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14068Ptr->availableCodelets[this->getID()] = 1;
#endif
        }
    }
}
void TP14059::_checkInCodelets14069::fire(void)
{
    if (this->getID() == 0) {
        /*Init the vars for this region*/
        /*Initialize the vars of the inlined region*/
        this->inputsTPParent->k_darts14069
            = &(this->inputsTPParent
                    ->k_darts14002[this->getLocalID()]) /*OMP_SHARED_PRIVATE - VAR INLINED*/;

        /*printing node 14070: CallExpr*/
        printf("blts complete... #k = %d\n", (*(this->inputsTPParent->k_darts14069)));
        /*Signaling next codelet from last stmt in the codelet*/
        /*The node is the last one in a complex loop, so signal the inc node*/
        myTP->controlTPParent->TPParent->checkInCodelets14063[this->getID()].decDep();
    } else {
        /*The node is the last one in a complex loop, so signal the inc node*/
        myTP->TPParent->checkInCodelets14063[this->getID()].decDep();
    }
}
TP14059::TP14059(int in_numThreads, int in_mainCodeletID, TP14002* in_TPParent,
    TP14002* in_inputsTPParent, TP14059** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(in_inputsTPParent)
    , ptrToThisTP(in_ptrToThisTP)
    , TP14065Ptr(nullptr)
    , TP14065_alreadyLaunched(0)
    , TP14066_alreadyLaunched(0)
    , TP14068Ptr(nullptr)
    , TP14068_alreadyLaunched(0)
    , TP14069_alreadyLaunched(0)
    , checkInCodelets14065(new _checkInCodelets14065[this->numThreads])
    , checkInCodelets14066(new _checkInCodelets14066[this->numThreads])
    , checkInCodelets14068(new _checkInCodelets14068[this->numThreads])
    , checkInCodelets14069(new _checkInCodelets14069[this->numThreads])
{
    /*Initialize Codelets*/
    _checkInCodelets14069* checkInCodelets14069Ptr = (this->checkInCodelets14069);
    _checkInCodelets14068* checkInCodelets14068Ptr = (this->checkInCodelets14068);
    _checkInCodelets14066* checkInCodelets14066Ptr = (this->checkInCodelets14066);
    _checkInCodelets14065* checkInCodelets14065Ptr = (this->checkInCodelets14065);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets14065);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets14069Ptr) = _checkInCodelets14069(1, 1, this, codeletCounter);
        checkInCodelets14069Ptr++;
        (*checkInCodelets14068Ptr) = _checkInCodelets14068(1, 1, this, codeletCounter);
        checkInCodelets14068Ptr++;
        (*checkInCodelets14066Ptr) = _checkInCodelets14066(1, 1, this, codeletCounter);
        checkInCodelets14066Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets14065Ptr) = _checkInCodelets14065(2, 1, this, codeletCounter);
#else
        (*checkInCodelets14065Ptr) = _checkInCodelets14065(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets14065Ptr).decDep();
        checkInCodelets14065Ptr++;
    }
    *(this->ptrToThisTP) = this;
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[this->getID()].decDep();
#else
    this->availableCodelets[this->getID()] = 1;
#endif
}
TP14059::~TP14059()
{
    delete[] checkInCodelets14069;
    delete[] checkInCodelets14068;
    delete[] checkInCodelets14066;
    delete[] checkInCodelets14065;
}
/*TP14072: ForStmt*/
void TP14072::_checkInCodelets14078::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif

    /*printing node 14078: CallExpr*/
    if (!__sync_val_compare_and_swap(&(myTP->controlTPParent->TP14078_alreadyLaunched), 0, 1)) {
        /*Make the function call*/
        invoke<TP_jacu>(myTP, myTP->numThreads, this->getID(),
            &(myTP->controlTPParent->checkInCodelets14079[this->getID()]),
            &(myTP->controlTPParent->TPParent->barrierCodelets14080[0]),
            &(myTP->controlTPParent->TP14078Ptr),
            (this->inputsTPParent->k_darts14002[this->getID()]));
    } else {
        if (myTP->controlTPParent->TP14078Ptr == nullptr) {
            myTP->add(this);
            return;
        } else {
            myTP->controlTPParent->TP14078Ptr->setNewInputs(
                (this->inputsTPParent->k_darts14002[this->getID()]), this->getID());
            myTP->controlTPParent->TP14078Ptr->nextCodeletsjacu[this->getID()]
                = &(myTP->controlTPParent->checkInCodelets14079[this->getID()]);
            myTP->controlTPParent->TP14078Ptr->nextSyncCodeletsjacu[this->getID()]
                = &(myTP->controlTPParent->TPParent->barrierCodelets14080[0]);
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14078Ptr->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14078Ptr->availableCodelets[this->getID()] = 1;
#endif
        }
    }
}
void TP14072::_checkInCodelets14079::fire(void)
{

    /*printing node 14079: CallExpr*/
    if (!__sync_val_compare_and_swap(&(myTP->controlTPParent->TP14079_alreadyLaunched), 0, 1)) {
        /*Make the function call*/
        invoke<TP_buts>(myTP, myTP->numThreads, this->getID(),
            &(myTP->controlTPParent->TPParent->checkInCodelets14076[this->getID()]),
            &(myTP->controlTPParent->TPParent->barrierCodelets14080[0]),
            &(myTP->controlTPParent->TP14079Ptr), nx, ny, nz,
            (this->inputsTPParent->k_darts14002[this->getID()]), omega, ist, iend, jst, jend, nx0,
            ny0);
    } else {
        if (myTP->controlTPParent->TP14079Ptr == nullptr) {
            myTP->add(this);
            return;
        } else {
            myTP->controlTPParent->TP14079Ptr->setNewInputs(nx, ny, nz,
                (this->inputsTPParent->k_darts14002[this->getID()]), omega, ist, iend, jst, jend,
                nx0, ny0, this->getID());
            myTP->controlTPParent->TP14079Ptr->nextCodeletsbuts[this->getID()]
                = &(myTP->controlTPParent->TPParent->checkInCodelets14076[this->getID()]);
            myTP->controlTPParent->TP14079Ptr->nextSyncCodeletsbuts[this->getID()]
                = &(myTP->controlTPParent->TPParent->barrierCodelets14080[0]);
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP14079Ptr->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP14079Ptr->availableCodelets[this->getID()] = 1;
#endif
        }
    }
}
TP14072::TP14072(int in_numThreads, int in_mainCodeletID, TP14002* in_TPParent,
    TP14002* in_inputsTPParent, TP14072** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(in_inputsTPParent)
    , ptrToThisTP(in_ptrToThisTP)
    , TP14078Ptr(nullptr)
    , TP14078_alreadyLaunched(0)
    , TP14079Ptr(nullptr)
    , TP14079_alreadyLaunched(0)
    , checkInCodelets14078(new _checkInCodelets14078[this->numThreads])
    , checkInCodelets14079(new _checkInCodelets14079[this->numThreads])
{
    /*Initialize Codelets*/
    _checkInCodelets14079* checkInCodelets14079Ptr = (this->checkInCodelets14079);
    _checkInCodelets14078* checkInCodelets14078Ptr = (this->checkInCodelets14078);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets14078);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets14079Ptr) = _checkInCodelets14079(1, 1, this, codeletCounter);
        checkInCodelets14079Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets14078Ptr) = _checkInCodelets14078(2, 1, this, codeletCounter);
#else
        (*checkInCodelets14078Ptr) = _checkInCodelets14078(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets14078Ptr).decDep();
        checkInCodelets14078Ptr++;
    }
    *(this->ptrToThisTP) = this;
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[this->getID()].decDep();
#else
    this->availableCodelets[this->getID()] = 1;
#endif
}
TP14072::~TP14072()
{
    delete[] checkInCodelets14079;
    delete[] checkInCodelets14078;
}
/*TP14081: OMPForDirective*/
void TP14081::_barrierCodelets14081::fire(void)
{
    TP14081* myTP = static_cast<TP14081*>(myTP_);
    myTP->controlTPParent->TPParent->barrierCodelets14081[0].decDep();
}
bool TP14081::requestNewRangeIterations14081(int* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        int tempStartRange = rangePerCodelet14081 * codeletID;
        int tempEndRange = rangePerCodelet14081 * (codeletID + 1);
        if (remainderRange14081 != 0) {
            if (codeletID < (uint32_t)remainderRange14081) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange14081;
                tempEndRange += remainderRange14081;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration14081;
        tempEndRange = tempEndRange * 1 + minIteration14081;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration14081 < lastIteration14081) {
            (this->inputsTPParent->i_darts14081[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts14081[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = *endRange + 1;
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration14081;
        }
    }
    return isThereNewIteration;
}
void TP14081::_checkInCodelets14082::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/

    /*printing node 14082: ForStmt*/
    /*var: i*/
    /*var: j*/
    /*var: k*/
    /*var: m*/
    /*var: tmp*/
    int* i = &(this->inputsTPParent->i_darts14081[this->getLocalID()]);
    (void)i /*OMP_PRIVATE*/;
    int* j = &(this->inputsTPParent->j_darts14081[this->getLocalID()]);
    (void)j /*OMP_PRIVATE*/;
    int* k = &(this->inputsTPParent->k_darts14081[this->getLocalID()]);
    (void)k /*OMP_PRIVATE*/;
    int* m = &(this->inputsTPParent->m_darts14081[this->getLocalID()]);
    (void)m /*OMP_PRIVATE*/;
    double* tmp = (this->inputsTPParent->tmp_darts14081);
    (void)tmp /*OMP_SHARED*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations14081(
        (int*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Signaling omp for stmt's barrier*/
        myTP->controlTPParent->barrierCodelets14081[0].decDep();
        return;
    }
    for (int i_darts_counter_temp14081 = (*i); i_darts_counter_temp14081 <= endRange
         && i_darts_counter_temp14081 <= this->inputsTPParent->lastIteration14081;
         i_darts_counter_temp14081++) {
        {
            {
                /*Loop's init*/
                (*j) = jst;
                int j_darts_counter_temp14081 = (*j);
                for (; j_darts_counter_temp14081 <= jend; j_darts_counter_temp14081++) {
                    {
                        /*Loop's init*/
                        (*k) = 1;
                        int k_darts_counter_temp14081 = (*k);
                        for (; k_darts_counter_temp14081 <= nz - 2; k_darts_counter_temp14081++) {
                            {
                                /*Loop's init*/
                                (*m) = 0;
                                int m_darts_counter_temp14081 = (*m);
                                for (; m_darts_counter_temp14081 < 5; m_darts_counter_temp14081++) {
                                    u[(i_darts_counter_temp14081)][j_darts_counter_temp14081]
                                     [k_darts_counter_temp14081][m_darts_counter_temp14081]
                                        = u[(i_darts_counter_temp14081)][j_darts_counter_temp14081]
                                           [k_darts_counter_temp14081][m_darts_counter_temp14081]
                                        + (*(tmp))
                                            * rsd[(i_darts_counter_temp14081)]
                                                 [j_darts_counter_temp14081]
                                                 [k_darts_counter_temp14081]
                                                 [m_darts_counter_temp14081];
                                }
                                (*m) = m_darts_counter_temp14081;
                            }
                        }
                        (*k) = k_darts_counter_temp14081;
                    }
                }
                (*j) = j_darts_counter_temp14081;
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling omp for stmt's barrier*/
    myTP->controlTPParent->barrierCodelets14081[0].decDep();
}
TP14081::TP14081(int in_numThreads, int in_mainCodeletID, TP14002* in_TPParent,
    int in_initIteration, int in_lastIteration, double* in_tmp, TP14081** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , i_darts14081(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , j_darts14081(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , k_darts14081(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , m_darts14081(new int[this->numThreads]) /*OMP_PRIVATE - INPUT*/
    , tmp_darts14081(in_tmp) /*OMP_SHARED - INPUT*/
    , initIteration14081(in_initIteration)
    , lastIteration14081(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , barrierCodelets14081(new _barrierCodelets14081[1])
    , checkInCodelets14082(new _checkInCodelets14082[this->numThreads])
{
    /*Initialize the loop parameters*/
    range14081 = abs(lastIteration14081 - initIteration14081) / 1;
    rangePerCodelet14081 = range14081 / numThreads;
    minIteration14081 = min<int>(lastIteration14081, initIteration14081);
    remainderRange14081 = range14081 % numThreads;
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets14081[0] = _barrierCodelets14081(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets14082* checkInCodelets14082Ptr = (this->checkInCodelets14082);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets14082);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets14082Ptr) = _checkInCodelets14082(2, 1, this, codeletCounter);
#else
        (*checkInCodelets14082Ptr) = _checkInCodelets14082(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets14082Ptr).decDep();
        checkInCodelets14082Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP14081::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets14082[localID].setID(codeletID);
    this->checkInCodelets14082[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets14082[localID + this->baseNumThreads * i]
            = _checkInCodelets14082(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets14082[localID + this->baseNumThreads * i]
            = _checkInCodelets14082(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets14082[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets14082[localID + this->baseNumThreads * i].decDep();
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
TP14081::~TP14081()
{
    delete[] barrierCodelets14081;
    delete[] checkInCodelets14082;
}
